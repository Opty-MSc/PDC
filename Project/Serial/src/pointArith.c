#include "pointArith.h"
#include <math.h>

extern int nDims;

double squareDistance(const double *P1, const double *P2) {
  return sqrt(distance(P1, P2));
}

double distance(const double *P1, const double *P2) {
  double distance = 0;
  for (int i = 0; i < nDims; i++) { distance += (P2[i] - P1[i]) * (P2[i] - P1[i]); }
  return sqrt(distance);
}

void copy(const double *P, double *copyP) {
  for (int i = 0; i < nDims; i++) { copyP[i] = P[i]; }
}

void middle(const double *P1, const double *P2, double *midP) {
  for (int i = 0; i < nDims; i++) { midP[i] = (P1[i] + P2[i]) / 2; }
}

void sum(const double *P1, const double *P2, double *sumP) {
  for (int i = 0; i < nDims; i++) { sumP[i] = P1[i] + P2[i]; }
}

void sub(const double *P1, const double *P2, double *subP) {
  for (int i = 0; i < nDims; i++) { subP[i] = P1[i] - P2[i]; }
}

double innerProduct(const double *P1, const double *P2) {
  double product = 0;
  for (int i = 0; i < nDims; i++) { product += P1[i] * P2[i]; }
  return product;
}

void scale(double scalar, const double *P, double *scaleP) {
  for (int i = 0; i < nDims; i++) { scaleP[i] = scalar * P[i]; }
}

void projection(const double *P, const double *A, const double *subBA, double squaredSubBA, double *projectionP) {
  sub(P, A, projectionP);
  scale(innerProduct(projectionP, subBA) / squaredSubBA, subBA, projectionP);
  sum(projectionP, A, projectionP);
}
