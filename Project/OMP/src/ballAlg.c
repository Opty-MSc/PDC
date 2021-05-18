#include "ballAlg.h"
#include "../lib/genPoints.h"
#include "../lib/msort.h"
#include "pointArith.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int nDims;
long nNodes = 0;

int main(int argc, char *argv[]) {
  long nPoints;
  Node head;
  double **points = getPoints(argc, argv, &nDims, &nPoints);
  double *pointsValues = *points;

  double execTime = -omp_get_wtime();
#pragma omp parallel
#pragma omp single
  head = buildTree(points, nPoints, omp_get_num_threads());
  execTime += omp_get_wtime();

  fprintf(stderr, "%.1lf\n", execTime);
  printf("%d %ld\n", nDims, nNodes);
  dumpTree(head);

  free(pointsValues);
  free(points);
  exit(EXIT_SUCCESS);
}

Node buildTree(double **points, long nPoints, int nThreads) {

  if (nPoints == 0) return NULL;

  double *center = (double *) mallocSafe(sizeof(double) * nDims);
  if (nPoints == 1) {
    copy(points[0], center);
    return newNode(center, 0, NULL, NULL);
  }

  const long iA = furthest(points, nPoints, points[0]);
  const long iB = furthest(points, nPoints, points[iA]);

  double *subBA = (double *) mallocSafe(sizeof(double) * nDims);
  double *projectionsPoints = (double *) mallocSafe(sizeof(double) * nDims * nPoints);
  double **projections = (double **) mallocSafe(sizeof(double *) * nPoints);
  double **pointsTmp = (double **) mallocSafe(sizeof(double *) * nPoints);

  sub(points[iB], points[iA], subBA);
  const double squaredSubBA = innerProduct(subBA, subBA);

  for (long i = 0; i < nPoints; i++) {
    projections[i] = projectionsPoints + (i * nDims);
    if (i == iA || i == iB) {
      copy(points[i], projections[i]);
    } else {
      projection(points[i], points[iA], subBA, squaredSubBA, projections[i]);
    }
  }

  msort(projections, nPoints, pointsTmp);
  if (nPoints % 2 == 0) {
    middle(projections[nPoints / 2 - 1], projections[nPoints / 2], center);
  } else {
    copy(projections[nPoints / 2], center);
  }
  const double radius = distance(center, points[furthest(points, nPoints, center)]);

  long nPointsL = 0;
  long nPointsR = 0;
  // Reusing Previous Allocated Arrays
  double **pointsL = projections;
  double **pointsR = pointsTmp;

  partitionTree(projectionsPoints, center, points, nPoints, pointsL, &nPointsL, pointsR, &nPointsR);
  pointsL = realloc(pointsL, sizeof(double *) * nPointsL);
  pointsR = realloc(pointsR, sizeof(double *) * nPointsR);

  free(subBA);
  free(projectionsPoints);

  Node L, R;
#pragma omp task shared(L) if (nThreads > 1)
  L = buildTree(pointsL, nPointsL, nThreads / 2);
  R = buildTree(pointsR, nPointsR, nThreads - nThreads / 2);
#pragma omp taskwait

  free(projections);
  free(pointsTmp);

  return newNode(center, radius, L, R);
}

Node newNode(double *center, double radius, Node L, Node R) {
#pragma omp atomic
  nNodes++;
  Node node = (Node) mallocSafe(sizeof(struct Node));
  node->center = center;
  node->radius = radius;
  node->L = L;
  node->R = R;
  return node;
}

long furthest(double **points, long nPoints, const double *pivot) {
  long iFurthest = 0;
  double maxD = -1;

  for (long i = 0; i < nPoints; i++) {
    if (points[i] != pivot) {
      double d = squareDistance(points[i], pivot);
      if (maxD < d) {
        maxD = d;
        iFurthest = i;
      }
    }
  }
  return iFurthest;
}

void partitionTree(const double *projectionsPoints, const double *center,
                   double **points, long nPoints,
                   double **pointsL, long *nPointsL,
                   double **pointsR, long *nPointsR) {

  for (long i = 0; i < nPoints; i++) {
    // projectionsPoints[i * nDims] == (projectionsPoints + (i * nDims))[0]
    if (projectionsPoints[i * nDims] < center[0]) {
      pointsL[(*nPointsL)++] = points[i];
    } else {
      pointsR[(*nPointsR)++] = points[i];
    }
  }
}

int dumpTree(Node head) {
  static long _id = 0;
  if (head == NULL) return -1;

  const long id = _id++;
  const long L = dumpTree(head->L);
  const long R = dumpTree(head->R);

  printf("%ld %ld %ld %.6lf", id, L, R, head->radius);
  printPoint(head->center, nDims);

  free(head->center);
  free(head);
  return id;
}

void *mallocSafe(size_t size) {

  void *allocBytes = malloc(size);

  if (allocBytes == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }
  return allocBytes;
}
