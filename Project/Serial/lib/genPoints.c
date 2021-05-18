#include "genPoints.h"
#include <stdio.h>
#include <stdlib.h>

double **getPoints(int argc, char *argv[], int *nDims, long *nPoints) {

  if (argc != 4) {
    fprintf(stderr, "[%s] Usage: [Number of Dimensions] [Number of Points] [Seed]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  *nDims = atoi(argv[1]);
  if (*nDims < 2) {
    fprintf(stderr, "ILLEGAL: Number of Dimensions (%d)!\n", *nDims);
    exit(EXIT_FAILURE);
  }

  *nPoints = atol(argv[2]);
  if (*nPoints < 1) {
    fprintf(stderr, "ILLEGAL: Number of Points (%ld)!\n", *nPoints);
    exit(EXIT_FAILURE);
  }

  srandom(atoi(argv[3]));

  double **points = (double **) newPointsArray(*nDims, *nPoints);

  for (long i = 0; i < *nPoints; i++)
    for (int j = 0; j < *nDims; j++) points[i][j] = RANGE * ((double) random()) / RAND_MAX;

#ifdef DEBUG
  printf("=== Generated Points ===\n");
  for (long i = 0; i < *nPoints; i++) printPoint(points[i], *nDims);
  printf("========================\n");
#endif

  return points;
}

double **newPointsArray(int nDims, long nPoints) {
  double *pointsValues = (double *) malloc(nDims * nPoints * sizeof(double));
  double **points = (double **) malloc(nPoints * sizeof(double *));

  if (pointsValues == NULL || points == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }

  for (long i = 0; i < nPoints; i++) points[i] = &pointsValues[i * nDims];
  return points;
}

void printPoint(double *P, int nDims) {
  for (int i = 0; i < nDims; i++) { printf(" %.6lf", P[i]); }
  printf("\n");
}
