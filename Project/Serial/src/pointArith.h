#ifndef SERIAL_POINT_ARITH_H
#define SERIAL_POINT_ARITH_H

double squareDistance(const double *P1, const double *P2);
double distance(const double *P1, const double *P2);
void copy(const double *P, double *copyP);
void middle(const double *P1, const double *P2, double *midP);
void sum(const double *P1, const double *P2, double *sumP);
void sub(const double *P1, const double *P2, double *subP);
double innerProduct(const double *P1, const double *P2);
void scale(double scalar, const double *P, double *scaleP);
void projection(const double *P, const double *A, const double *subBA, double squaredSubBA, double *projectionP);

#endif
