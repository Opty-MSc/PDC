#ifndef OMP_GEN_POINTS_H
#define OMP_GEN_POINTS_H

#define RANGE 10

double **getPoints(int argc, char *argv[], int *nDims, long *nPoints);
double **newPointsArray(int nDims, long nPoints);
void printPoint(double *P, int nDims);

#endif
