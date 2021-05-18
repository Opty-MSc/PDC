#ifndef MPI_BALL_ALG_H
#define MPI_BALL_ALG_H

#include <mpi.h>
#include <stdbool.h>

// MPI Tags
#define TEAMMATE_STATE_TAG 0
#define POINT_A_TAG 1
#define POINT_B_TAG 2
#define PROJECTIONS_LEN_TAG 3
#define PROJECTIONS_TAG 4
#define MEDIAN_REQUEST 5
#define MEDIAN_REPLY 6
#define RADIUS_TAG 7
#define POINTS_LEN_TAG 8
#define POINTS_TAG 9
#define INITIAL_POINT_TAG 10
#define BRANCH_ID_LEFT 11
#define BRANCH_ID_RIGHT 12
#define PRINT_TAG 13

#define WORKING 0
#define FINISHED 1
#define MY_STATE(nPoints) ((nPoints) == 0 ? FINISHED : WORKING)

typedef struct {
  double *center;
  double radius;
  int nid;
  int nidL;
  int nidR;
} Node;

typedef struct {
  double medX;
  int medRanks[2];
  int medIdx[2];
} MedianInfo;

void defineMedianInfo();
int buildTree(double **initialP, double **points, int nPoints, const int *teammatesRanks, int nTeammates);
void bcastToMyTeam(void *buf, int bufSize, const int *teammatesRanks, int nTeammates, MPI_Datatype datatype, int TAG);
double *calcFurthestPoint(double **points, int nPoints, const double *pivot, const int *teammatesRanks, int nTeammates, int TAG);
int calcFurthestIdx(double **points, int nPoints, const double *pivot, double *maxD);
MedianInfo bcastMedianInfo(const int *teammatesRanks, int nTeammates, double *projectionsXs, int nProjectionsXs);
int teammateMinX(int nTeammates, double **teammatesXs, const int *iTeammatesXs, const int *nTeammatesXs);
MedianInfo recvMedianInfo(int leaderRank, double **projections, const double *projectionsXs, int *nProjectionsXs);
void calcCenter(MedianInfo medInfo, double **projections, double *center);
void calcRadius(double **points, int nPoints, double *center, const int *teammatesRanks, int nTeammates, double *radius);
void calcCandidateRadius(int leaderRank, double **points, int nPoints, double *center);
void partitionTree(const double *projectionsPoints, double medX,
                   double **points, int nPoints,
                   double **pointsL, int *nPointsL,
                   double **pointsR, int *nPointsR);
int buildTreeLoop(double **initialP, double *center, double radius,
                  double ***pointsL, int *nPointsL,
                  double ***pointsR, int *nPointsR,
                  const int *teammatesRanks, int nTeammates);
void exchangePoints(double ***pointsToSend, int *nPointsToSend,
                    double ***pointsToRecv, int *nPointsToRecv,
                    int teammateRank, bool toMergeLeft);
void flat(double ***points, int *nPoints, double *flattedPoints);
void unflat(double ***points, int *nPoints, double *flattedPoints, int nFlattedPoints, bool toMergeLeft);
int *calcInitialTeammates(int myState, int *nTeammates);
int *calcNewTeammates(int myState, const int *teammatesRanks, int nTeammates, int *newNTeammates, int iParity);
int calcWorkingTeammates(int myState, int *teammatesRanks, int nTeammates);
int newNode(double *center, double radius, int nidL, int nidR);
void dumpTree();
void *mallocSafe(size_t size);
int buildTreeOMP(double **points, long nPoints, int nThreads);

#endif
