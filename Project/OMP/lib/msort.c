#include "msort.h"
#include <string.h>

void msort(double **points, long nPoints, double **pointsTmp) {
  if (nPoints > 1) {
    msort(points, (nPoints / 2), pointsTmp);
    msort(points + (nPoints / 2), nPoints - (nPoints / 2), pointsTmp);
    merge(points, nPoints, pointsTmp);
  }
}

void merge(double **points, long nPoints, double **pointsTmp) {
  long i = 0;
  long iL = 0;
  long iR = nPoints / 2;

  while (iL < nPoints / 2 && iR < nPoints) {
    // Sorts By The First Coordinate Of The Points
    pointsTmp[i++] = points[points[iL][0] < points[iR][0] ? iL++ : iR++];
  }
  while (iL < nPoints / 2) pointsTmp[i++] = points[iL++];
  while (iR < nPoints) pointsTmp[i++] = points[iR++];
  memcpy(points, pointsTmp, nPoints * sizeof(double *));
}
