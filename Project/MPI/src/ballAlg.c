#include "ballAlg.h"
#include "../lib/genPoints.h"
#include "../lib/msort.h"
#include "pointArith.h"
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int nProcesses;
int myRank;
int nDims;
int nNodes = 0;
int nid = 1;
int nodesCapacity;
Node *nodes;
MPI_Datatype mpiMedianInfo;

int main(int argc, char *argv[]) {

  double execTime = -omp_get_wtime();
  MPI_Init(&argc, &argv);
  defineMedianInfo();

  MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  nid += myRank;

  int nPoints;
  double **points = getPoints(argc, argv, &nDims, &nPoints, myRank, nProcesses);
  double *pointsValues = *points;

  nodesCapacity = nPoints;
  nodes = (Node *) mallocSafe(sizeof(Node) * nodesCapacity);

  int nTeammates = nProcesses;
  int *teammatesRanks = nProcesses > 1 ? calcInitialTeammates(MY_STATE(nPoints), &nTeammates) : NULL;

#pragma omp parallel
#pragma omp single
  nTeammates == 1
          ? buildTreeOMP(points, nPoints, omp_get_num_threads())
          : buildTree(&points[nPoints], points, nPoints, teammatesRanks, nTeammates);

  if (myRank == 0) {
    if (nNodes > 0) nodes[nNodes - 1].nid = 0;
    execTime += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", execTime);
    fflush(stderr);
  }

  int nNodesGlobal = nNodes;
  if (nProcesses > 1) MPI_Reduce(&nNodes, &nNodesGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myRank == 0) {
    printf("%d %d\n", nDims, nNodesGlobal);
    fflush(stdout);
  }
  dumpTree();

  if (nProcesses > 1) free(teammatesRanks);
  free(pointsValues);
  free(points);
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}

void defineMedianInfo() {
  MedianInfo dummyMedianInfo;
  int lengths[3] = {1, 2, 2};
  MPI_Aint displacements[3];
  MPI_Aint baseAddress;
  MPI_Get_address(&dummyMedianInfo, &baseAddress);
  MPI_Get_address(&dummyMedianInfo.medX, &displacements[0]);
  MPI_Get_address(dummyMedianInfo.medRanks, &displacements[1]);
  MPI_Get_address(dummyMedianInfo.medIdx, &displacements[2]);
  displacements[0] = MPI_Aint_diff(displacements[0], baseAddress);
  displacements[1] = MPI_Aint_diff(displacements[1], baseAddress);
  displacements[2] = MPI_Aint_diff(displacements[2], baseAddress);
  MPI_Datatype datatypes[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
  MPI_Type_create_struct(3, lengths, displacements, datatypes, &mpiMedianInfo);
  MPI_Type_commit(&mpiMedianInfo);
}

int buildTree(double **initialP, double **points, int nPoints, const int *teammatesRanks, int nTeammates) {

  if (nPoints == 0) return -1;

  double *center = (double *) mallocSafe(sizeof(double) * nDims);
  if (nPoints == 1 && nTeammates == 1) {
    copy(points[0], center);
    return newNode(center, 0, -1, -1);
  }

  double *pA = calcFurthestPoint(points, nPoints, *initialP, teammatesRanks, nTeammates, POINT_A_TAG);
  double *pB = calcFurthestPoint(points, nPoints, pA, teammatesRanks, nTeammates, POINT_B_TAG);

  double *subBA = (double *) mallocSafe(sizeof(double) * nDims);
  double *projectionsXs = (double *) mallocSafe(sizeof(double) * nPoints);
  double *projectionsPoints = (double *) mallocSafe(sizeof(double) * nDims * nPoints);
  double **projections = (double **) mallocSafe(sizeof(double *) * nPoints);
  double **pointsTmp = (double **) mallocSafe(sizeof(double *) * nPoints);
  double **pointsL = (double **) mallocSafe(sizeof(double *) * nPoints);
  double **pointsR = (double **) mallocSafe(sizeof(double *) * nPoints);

  sub(pB, pA, subBA);
  const double squaredSubBA = innerProduct(subBA, subBA);

  for (int i = 0; i < nPoints; i++) {
    projections[i] = projectionsPoints + (i * nDims);
    projection(points[i], pA, subBA, squaredSubBA, projections[i]);
  }

  msort(projections, nPoints, pointsTmp);
  for (int i = 0; i < nPoints; i++) projectionsXs[i] = projections[i][0];

  MedianInfo medInfo;
  if (teammatesRanks[0] == myRank) {
    medInfo = bcastMedianInfo(teammatesRanks, nTeammates, projectionsXs, nPoints);
  } else {
    medInfo = recvMedianInfo(teammatesRanks[0], projections, projectionsXs, &nPoints);
  }

  int nPointsL = 0;
  int nPointsR = 0;
  partitionTree(projectionsPoints, medInfo.medX, points, nPoints, pointsL, &nPointsL, pointsR, &nPointsR);
  pointsL = realloc(pointsL, sizeof(double *) * nPointsL);
  pointsR = realloc(pointsR, sizeof(double *) * nPointsR);

  double radius = -1;
  if (teammatesRanks[0] == myRank) {
    calcCenter(medInfo, projections, center);
    calcRadius(points, nPoints, center, teammatesRanks, nTeammates, &radius);
  } else {
    calcCandidateRadius(teammatesRanks[0], points, nPoints, center);
  }

  free(pA);
  free(pB);
  free(subBA);
  free(projectionsXs);
  free(projectionsPoints);
  free(projections);
  free(pointsTmp);

  int myNid = buildTreeLoop(initialP, center, radius, &pointsL, &nPointsL, &pointsR, &nPointsR, teammatesRanks, nTeammates);

  free(pointsL);
  free(pointsR);
  return myNid;
}

void bcastToMyTeam(void *buf, int bufSize, const int *teammatesRanks, int nTeammates, MPI_Datatype datatype, int TAG) {

  MPI_Request request;
  for (int i = 0; i < nTeammates; i++) {
    if (teammatesRanks[i] != myRank) {
      MPI_Isend(buf, bufSize, datatype, teammatesRanks[i], TAG, MPI_COMM_WORLD, &request);
      MPI_Request_free(&request);
    }
  }
}

double *calcFurthestPoint(double **points, int nPoints, const double *pivot, const int *teammatesRanks, int nTeammates, int TAG) {

  double *P = (double *) mallocSafe(sizeof(double) * nDims);
  double *pCmp = (double *) mallocSafe(sizeof(double) * nDims);

  double maxD;
  int iFurthest = calcFurthestIdx(points, nPoints, pivot, &maxD);
  copy(points[iFurthest], P);
  bcastToMyTeam(P, nDims, teammatesRanks, nTeammates, MPI_DOUBLE, TAG);

  for (int i = 0; i < nTeammates; i++) {
    if (teammatesRanks[i] != myRank) {
      MPI_Recv(pCmp, nDims, MPI_DOUBLE, teammatesRanks[i], TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      double d = squareDistance(pivot, pCmp);
      if (d > maxD) {
        maxD = d;
        copy(pCmp, P);
      }
    }
  }

  free(pCmp);
  return P;
}

int calcFurthestIdx(double **points, int nPoints, const double *pivot, double *maxD) {
  int iFurthest = 0;
  *maxD = -1;

  for (int i = 0; i < nPoints; i++) {
    if (points[i] != pivot) {
      double d = squareDistance(points[i], pivot);
      if ((*maxD) < d) {
        *maxD = d;
        iFurthest = i;
      }
    }
  }
  return iFurthest;
}

MedianInfo bcastMedianInfo(const int *teammatesRanks, int nTeammates, double *projectionsXs, int nProjectionsXs) {

  int nTeammatesXsSum = 0;
  int *nTeammatesXs = (int *) calloc(nTeammates, sizeof(int));
  int *iTeammatesXs = (int *) calloc(nTeammates, sizeof(int));
  double **teammatesXs = (double **) mallocSafe(sizeof(double *) * nTeammates);

  teammatesXs[0] = projectionsXs;
  nTeammatesXs[0] = nProjectionsXs;
  nTeammatesXsSum += nTeammatesXs[0];
  for (int i = 1; i < nTeammates; i++) {
    MPI_Recv(&nTeammatesXs[i], 1, MPI_INT, teammatesRanks[i], PROJECTIONS_LEN_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    teammatesXs[i] = (double *) mallocSafe(sizeof(double) * nTeammatesXs[i]);
    MPI_Recv(teammatesXs[i], nTeammatesXs[i], MPI_DOUBLE, teammatesRanks[i], PROJECTIONS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    nTeammatesXsSum += nTeammatesXs[i];
  }
  for (int i = 0; i < nTeammatesXsSum / 2 - 1; i++) {
    iTeammatesXs[teammateMinX(nTeammates, teammatesXs, iTeammatesXs, nTeammatesXs)]++;
  }

  MedianInfo medInfo;
  medInfo.medX = 0;
  for (int i = 0; i < 2; i++) {
    int teammateId = teammateMinX(nTeammates, teammatesXs, iTeammatesXs, nTeammatesXs);
    if (i == 1 || nTeammatesXsSum % 2 == 0) {
      medInfo.medX += teammatesXs[teammateId][iTeammatesXs[teammateId]];
    }
    medInfo.medRanks[i] = teammatesRanks[teammateId];
    medInfo.medIdx[i] = iTeammatesXs[teammateId]++;
  }
  if (nTeammatesXsSum % 2 == 0) medInfo.medX /= 2;
  else {
    medInfo.medRanks[0] = -1;
    medInfo.medIdx[0] = -1;
  }
  bcastToMyTeam(&medInfo, 1, teammatesRanks, nTeammates, mpiMedianInfo, MEDIAN_REQUEST);

  for (int i = 1; i < nTeammates; i++) free(teammatesXs[i]);
  free(nTeammatesXs);
  free(iTeammatesXs);
  free(teammatesXs);
  return medInfo;
}

int teammateMinX(int nTeammates, double **teammatesXs, const int *iTeammatesXs, const int *nTeammatesXs) {
  int iMinX = -1;

  for (int i = 0; i < nTeammates; i++) {
    if (iTeammatesXs[i] < nTeammatesXs[i] && (iMinX < 0 || teammatesXs[i][iTeammatesXs[i]] < teammatesXs[iMinX][iTeammatesXs[iMinX]])) {
      iMinX = i;
    }
  }
  return iMinX;
}

MedianInfo recvMedianInfo(int leaderRank, double **projections, const double *projectionsXs, int *nProjectionsXs) {

  MPI_Request request;
  MPI_Isend(nProjectionsXs, 1, MPI_INT, leaderRank, PROJECTIONS_LEN_TAG, MPI_COMM_WORLD, &request);
  MPI_Request_free(&request);
  MPI_Isend(projectionsXs, (*nProjectionsXs), MPI_DOUBLE, leaderRank, PROJECTIONS_TAG, MPI_COMM_WORLD, &request);
  MPI_Request_free(&request);

  MedianInfo medInfo;
  MPI_Recv(&medInfo, 1, mpiMedianInfo, leaderRank, MEDIAN_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  int i = 0;
  double *medPoints = (double *) mallocSafe(sizeof(double) * nDims * 2);

  for (int j = 0; j < 2; j++) {
    if (medInfo.medRanks[j] == myRank) copy(projections[medInfo.medIdx[j]], &medPoints[nDims * i++]);
  }
  if (i > 0) {
    MPI_Isend(medPoints, nDims * i, MPI_DOUBLE, leaderRank, MEDIAN_REPLY, MPI_COMM_WORLD, &request);
    MPI_Request_free(&request);
  }
  return medInfo;
}

void calcCenter(MedianInfo medInfo, double **projections, double *center) {

  double *medPoints = (double *) mallocSafe(sizeof(double) * nDims * 2);

  if (medInfo.medRanks[0] < 0) {
    if (medInfo.medRanks[1] == myRank) {
      copy(projections[medInfo.medIdx[1]], center);
    } else {
      MPI_Recv(medPoints, nDims, MPI_DOUBLE, medInfo.medRanks[1], MEDIAN_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      copy(medPoints, center);
    }
  } else if (medInfo.medRanks[0] == medInfo.medRanks[1]) {
    if (medInfo.medRanks[0] == myRank) {
      middle(projections[medInfo.medIdx[0]], projections[medInfo.medIdx[1]], center);
    } else {
      MPI_Recv(medPoints, nDims * 2, MPI_DOUBLE, medInfo.medRanks[0], MEDIAN_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      middle(medPoints, &medPoints[nDims], center);
    }
  } else {
    for (int i = 0; i < 2; i++) {
      if (medInfo.medRanks[i] == myRank) {
        copy(projections[medInfo.medIdx[i]], &medPoints[i * nDims]);
      } else {
        MPI_Recv(&medPoints[i * nDims], nDims, MPI_DOUBLE, medInfo.medRanks[i], MEDIAN_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    middle(medPoints, &medPoints[nDims], center);
  }

  free(medPoints);
}

void calcRadius(double **points, int nPoints, double *center, const int *teammatesRanks, int nTeammates, double *radius) {

  bcastToMyTeam(center, nDims, teammatesRanks, nTeammates, MPI_DOUBLE, RADIUS_TAG);

  double maxD;
  int iFurthest = calcFurthestIdx(points, nPoints, center, &maxD);
  *radius = distance(center, points[iFurthest]);

  double candidateRadius;
  for (int i = 0; i < nTeammates; i++) {
    if (teammatesRanks[i] != myRank) {
      MPI_Recv(&candidateRadius, 1, MPI_DOUBLE, teammatesRanks[i], RADIUS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (candidateRadius > (*radius)) *radius = candidateRadius;
    }
  }
}

void calcCandidateRadius(int leaderRank, double **points, int nPoints, double *center) {

  MPI_Recv(center, nDims, MPI_DOUBLE, leaderRank, RADIUS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  double maxD;
  int iFurthest = calcFurthestIdx(points, nPoints, center, &maxD);
  double candidateRadius = distance(center, points[iFurthest]);

  MPI_Request request;
  MPI_Isend(&candidateRadius, 1, MPI_DOUBLE, leaderRank, RADIUS_TAG, MPI_COMM_WORLD, &request);
  MPI_Request_free(&request);
}

void partitionTree(const double *projectionsPoints, double medX,
                   double **points, int nPoints,
                   double **pointsL, int *nPointsL,
                   double **pointsR, int *nPointsR) {

  for (int i = 0; i < nPoints; i++) {
    // projectionsPoints[i * nDims] == (projectionsPoints + (i * nDims))[0]
    if (projectionsPoints[i * nDims] < medX) {
      pointsL[(*nPointsL)++] = points[i];
    } else {
      pointsR[(*nPointsR)++] = points[i];
    }
  }
}

int buildTreeLoop(double **initialP, double *center, double radius, double ***pointsL, int *nPointsL, double ***pointsR, int *nPointsR, const int *teammatesRanks, int nTeammates) {

  int nidL = -1;
  int nidR = -1;
  if (nTeammates == 1) {
    int nThreads = omp_get_num_threads();
#pragma omp task shared(nidL)
    nidL = buildTreeOMP((*pointsL), (*nPointsL), nThreads / 2);
    nidR = buildTreeOMP((*pointsR), (*nPointsR), nThreads - nThreads / 2);
#pragma omp taskwait
    return newNode(center, radius, nidL, nidR);
  }

  int teammateId = 0;
  for (int i = 0; i < nTeammates; i++) {
    if (teammatesRanks[i] == myRank) {
      teammateId = i;
      break;
    }
  }
  if (teammateId % 2 == 0) {
    if (teammateId == nTeammates - 1) {
      exchangePoints(pointsR, nPointsR, NULL, NULL, teammatesRanks[teammateId - 1], false);
    } else {
      exchangePoints(pointsR, nPointsR, pointsL, nPointsL, teammatesRanks[teammateId + 1], false);
    }
  } else {
    exchangePoints(pointsL, nPointsL, pointsR, nPointsR, teammatesRanks[teammateId - 1], true);
    if (teammateId == nTeammates - 2) {
      exchangePoints(NULL, NULL, pointsR, nPointsR, teammatesRanks[teammateId + 1], false);
    }
  }

  int myState = MY_STATE(teammateId % 2 == 0 ? (*nPointsL) : (*nPointsR));
  bcastToMyTeam(&myState, 1, teammatesRanks, nTeammates, MPI_INT, TEAMMATE_STATE_TAG);
  if (myState == FINISHED && teammatesRanks[0] != myRank) return -1;

  int newNTeammates[2] = {0, 0};
  int *newTeammatesRanks[2];
  for (int i = 0; i < 2; i++) {
    newTeammatesRanks[i] = calcNewTeammates(myState, teammatesRanks, nTeammates, &newNTeammates[i], i);
  }

  if (myState == FINISHED && teammatesRanks[0] == myRank) {
    if (newNTeammates[0] > 0) MPI_Recv(&nidL, 1, MPI_INT, newTeammatesRanks[0][0], BRANCH_ID_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (newNTeammates[1] > 0) MPI_Recv(&nidR, 1, MPI_INT, newTeammatesRanks[1][0], BRANCH_ID_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < 2; i++) free(newTeammatesRanks[i]);
    return newNode(center, radius, nidL, nidR);
  }

  if (newTeammatesRanks[teammateId % 2][0] == myRank) {
    *initialP = teammateId % 2 == 0 ? (*pointsL)[0] : (*pointsR)[0];
    bcastToMyTeam(*initialP, nDims, newTeammatesRanks[teammateId % 2], newNTeammates[teammateId % 2], MPI_DOUBLE, INITIAL_POINT_TAG);
  } else {
    MPI_Recv(*initialP, nDims, MPI_DOUBLE, newTeammatesRanks[teammateId % 2][0], INITIAL_POINT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  int myNid;
  if (teammateId % 2 == 0) {
    nidL = buildTree(initialP, (*pointsL), (*nPointsL), newTeammatesRanks[0], newNTeammates[0]);
    myNid = nidL;
    if (newTeammatesRanks[0][0] == myRank) {
      if (teammatesRanks[0] == myRank) {
        if (newNTeammates[1] > 0) {
          MPI_Recv(&nidR, 1, MPI_INT, newTeammatesRanks[1][0], BRANCH_ID_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        myNid = newNode(center, radius, nidL, nidR);
      } else {
        MPI_Request request;
        MPI_Isend(&nidL, 1, MPI_INT, teammatesRanks[0], BRANCH_ID_LEFT, MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);
      }
    }
  } else {
    nidR = buildTree(initialP, (*pointsR), (*nPointsR), newTeammatesRanks[1], newNTeammates[1]);
    myNid = nidR;
    if (newTeammatesRanks[1][0] == myRank) {
      MPI_Request request;
      MPI_Isend(&nidR, 1, MPI_INT, teammatesRanks[0], BRANCH_ID_RIGHT, MPI_COMM_WORLD, &request);
      MPI_Request_free(&request);
    }
  }

  for (int i = 0; i < 2; i++) free(newTeammatesRanks[i]);
  return myNid;
}

void exchangePoints(double ***pointsToSend, int *nPointsToSend, double ***pointsToRecv, int *nPointsToRecv, int teammateRank, bool toMergeLeft) {

  if (nPointsToSend != NULL) {
    MPI_Request request;
    int nFlattedPointsToSend = (*nPointsToSend) * nDims;
    MPI_Isend(&nFlattedPointsToSend, 1, MPI_INT, teammateRank, POINTS_LEN_TAG, MPI_COMM_WORLD, &request);
    MPI_Request_free(&request);
    double *flattedPointsToSend = (double *) mallocSafe(sizeof(double) * nFlattedPointsToSend);
    flat(pointsToSend, nPointsToSend, flattedPointsToSend);
    MPI_Isend(flattedPointsToSend, nFlattedPointsToSend, MPI_DOUBLE, teammateRank, POINTS_TAG, MPI_COMM_WORLD, &request);
    MPI_Request_free(&request);
  }

  if (nPointsToRecv != NULL) {
    int nFlattedPointsToRecv;
    MPI_Recv(&nFlattedPointsToRecv, 1, MPI_INT, teammateRank, POINTS_LEN_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double *flattedPointsToRecv = (double *) mallocSafe(sizeof(double) * nFlattedPointsToRecv);
    MPI_Recv(flattedPointsToRecv, nFlattedPointsToRecv, MPI_DOUBLE, teammateRank, POINTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    unflat(pointsToRecv, nPointsToRecv, flattedPointsToRecv, nFlattedPointsToRecv, toMergeLeft);
  }
}

void flat(double ***points, int *nPoints, double *flattedPoints) {

  for (int i = 0; i < (*nPoints); i++) copy((*points)[i], &flattedPoints[i * nDims]);
}

void unflat(double ***points, int *nPoints, double *flattedPoints, int nFlattedPoints, bool toMergeLeft) {

  *points = (double **) realloc((*points), ((*nPoints) + nFlattedPoints / nDims) * sizeof(double *));
  if (toMergeLeft) {
    for (int i = (*nPoints) - 1; i >= 0; i--) (*points)[i + nFlattedPoints / nDims] = (*points)[i];
    for (int i = 0; i < nFlattedPoints / nDims; i++) (*points)[i] = &flattedPoints[i * nDims];
  } else {
    for (int i = 0; i < nFlattedPoints / nDims; i++) (*points)[(*nPoints) + i] = &flattedPoints[i * nDims];
  }
  *nPoints += nFlattedPoints / nDims;
}

int *calcInitialTeammates(int myState, int *nTeammates) {

  int *teammatesRanks = mallocSafe(sizeof(int) * (*nTeammates));
  for (int i = 0; i < (*nTeammates); i++) teammatesRanks[i] = i;

  bcastToMyTeam(&myState, 1, teammatesRanks, (*nTeammates), MPI_INT, TEAMMATE_STATE_TAG);
  *nTeammates = calcWorkingTeammates(myState, teammatesRanks, (*nTeammates));
  teammatesRanks = (int *) realloc(teammatesRanks, sizeof(int) * (*nTeammates));
  return teammatesRanks;
}

int *calcNewTeammates(int myState, const int *teammatesRanks, int nTeammates, int *newNTeammates, int iParity) {

  int *newTeammatesRanks = (int *) mallocSafe(sizeof(int) * nTeammates);
  for (int i = 0; i < nTeammates; i++) {
    if (i % 2 == iParity) newTeammatesRanks[(*newNTeammates)++] = teammatesRanks[i];
  }
  *newNTeammates = calcWorkingTeammates(myState, newTeammatesRanks, (*newNTeammates));
  newTeammatesRanks = (int *) realloc(newTeammatesRanks, sizeof(int) * (*newNTeammates));
  return newTeammatesRanks;
}

int calcWorkingTeammates(int myState, int *teammatesRanks, int nTeammates) {

  int teammateState;
  int nTeammatesWorking = 0;
  for (int i = 0; i < nTeammates; i++) {
    if (teammatesRanks[i] != myRank) {
      MPI_Recv(&teammateState, 1, MPI_INT, teammatesRanks[i], TEAMMATE_STATE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      teammateState = myState;
    }
    if (teammateState == WORKING) teammatesRanks[nTeammatesWorking++] = teammatesRanks[i];
  }
  return nTeammatesWorking;
}

int newNode(double *center, double radius, int nidL, int nidR) {
  int myNNodes, myNid;

#pragma omp critical(newNode)
  {
    myNNodes = nNodes++;
    myNid = nid;
    nid += nProcesses;

    if (nNodes > nodesCapacity) {
      nodesCapacity *= 2;
      nodes = (Node *) realloc(nodes, sizeof(Node) * nodesCapacity);
      if (nodes == NULL) {
        fprintf(stderr, "FATAL: [realloc]!\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  Node *new = &nodes[myNNodes];
  new->nid = myNid;
  new->center = center;
  new->radius = radius;
  new->nidL = nidL;
  new->nidR = nidR;
  return myNid;
}

void dumpTree() {

  if (myRank != 0) MPI_Recv(NULL, 0, MPI_INT, myRank - 1, PRINT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for (int i = 0; i < nNodes; i++) {
    printf("%d %d %d %.6lf", nodes[i].nid, nodes[i].nidL, nodes[i].nidR, nodes[i].radius);
    printPoint(nodes[i].center, nDims);
    free(nodes[i].center);
  }
  fflush(stdout);
  free(nodes);

  if (myRank != nProcesses - 1) MPI_Send(NULL, 0, MPI_INT, myRank + 1, PRINT_TAG, MPI_COMM_WORLD);
}

void *mallocSafe(size_t size) {

  void *allocBytes = malloc(size);

  if (allocBytes == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }
  return allocBytes;
}
