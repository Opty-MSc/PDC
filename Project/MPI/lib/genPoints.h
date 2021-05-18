#ifndef OMP_GEN_POINTS_H
#define OMP_GEN_POINTS_H

#define RANGE 10

double **getPoints(int argc, char *argv[], int *nDims, int *nPoints, int myRank, int nProcesses);
double **newPointsArray(int nDims, int nPoints);
void printPoint(double *P, int nDims);

#endif
