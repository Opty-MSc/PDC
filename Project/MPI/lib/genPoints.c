#include "genPoints.h"
#include <stdio.h>
#include <stdlib.h>

double **getPoints(int argc, char *argv[], int *nDims, int *nPoints, int myRank, int nProcesses) {

  if (argc != 4) {
    fprintf(stderr, "[%s] Usage: [Number of Dimensions] [Number of Points] [Seed]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  *nDims = atoi(argv[1]);
  if (*nDims < 2) {
    fprintf(stderr, "ILLEGAL: Number of Dimensions (%d)!\n", *nDims);
    exit(EXIT_FAILURE);
  }

  long nPointsGlobal = atol(argv[2]);
  if (nPointsGlobal < 1) {
    fprintf(stderr, "ILLEGAL: Number of Points (%ld)!\n", nPointsGlobal);
    exit(EXIT_FAILURE);
  }

  srandom(atoi(argv[3]));

  int iDiv = nPointsGlobal / nProcesses;
  int iMod = nPointsGlobal % nProcesses;
  int padding = iDiv * myRank + (myRank < iMod ? myRank : iMod);
  *nPoints = iDiv + (myRank < iMod ? 1 : 0);

  double **points = (double **) newPointsArray(*nDims, *nPoints + 1);

  for (int j = 0; j < *nDims; j++) {
    points[0][j] = points[*nPoints][j] = RANGE * ((double) random()) / RAND_MAX;
  }

  for (int i = 1; i < padding; i++) {
    for (int j = 0; j < *nDims; j++) random();
  }

  for (int i = (myRank == 0 ? 1 : 0); i < *nPoints; i++) {
    for (int j = 0; j < *nDims; j++) {
      points[i][j] = RANGE * ((double) random()) / RAND_MAX;
    }
  }

#ifdef DEBUG
  printf("=== Generated Points ===\n");
  for (int i = 0; i < *nPoints; i++) printPoint(points[i], *nDims);
  printf("========================\n");
#endif

  return points;
}

double **newPointsArray(int nDims, int nPoints) {
  double *pointsValues = (double *) malloc(nDims * nPoints * sizeof(double));
  double **points = (double **) malloc(nPoints * sizeof(double *));

  if (pointsValues == NULL || points == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < nPoints; i++) points[i] = &pointsValues[i * nDims];
  return points;
}

void printPoint(double *P, int nDims) {
  for (int i = 0; i < nDims; i++) { printf(" %.6lf", P[i]); }
  printf("\n");
}
