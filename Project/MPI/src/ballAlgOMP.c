#include "../lib/msort.h"
#include "ballAlg.h"
#include "pointArith.h"
#include <stdlib.h>

extern int nDims;

int buildTreeOMP(double **points, long nPoints, int nThreads) {

  if (nPoints == 0) return -1;

  double *center = (double *) mallocSafe(sizeof(double) * nDims);
  if (nPoints == 1) {
    copy(points[0], center);
    return newNode(center, 0, -1, -1);
  }

  double maxD;
  const int iA = calcFurthestIdx(points, nPoints, points[0], &maxD);
  const int iB = calcFurthestIdx(points, nPoints, points[iA], &maxD);

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
  const double radius = distance(center, points[calcFurthestIdx(points, nPoints, center, &maxD)]);

  int nPointsL = 0;
  int nPointsR = 0;
  // Reusing Previous Allocated Arrays
  double **pointsL = projections;
  double **pointsR = pointsTmp;

  partitionTree(projectionsPoints, center[0], points, nPoints, pointsL, &nPointsL, pointsR, &nPointsR);
  pointsL = realloc(pointsL, sizeof(double *) * nPointsL);
  pointsR = realloc(pointsR, sizeof(double *) * nPointsR);

  free(subBA);
  free(projectionsPoints);

  int nidL, nidR;
#pragma omp task shared(nidL) if (nThreads > 1)
  nidL = buildTreeOMP(pointsL, nPointsL, nThreads / 2);
  nidR = buildTreeOMP(pointsR, nPointsR, nThreads - nThreads / 2);
#pragma omp taskwait

  free(projections);
  free(pointsTmp);

  return newNode(center, radius, nidL, nidR);
}
