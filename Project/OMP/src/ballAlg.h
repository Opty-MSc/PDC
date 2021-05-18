#ifndef OMP_BALL_ALG_H
#define OMP_BALL_ALG_H

#include <stddef.h>

typedef struct Node *Node;
struct Node {
  double *center;
  double radius;
  Node L;
  Node R;
};

Node buildTree(double **points, long nPoints, int nThreads);
Node newNode(double *center, double radius, Node L, Node R);
long furthest(double **points, long nPoints, const double *pivot);
void partitionTree(const double *projectionsPoints, const double *center,
                   double **points, long nPoints,
                   double **pointsL, long *nPointsL,
                   double **pointsR, long *nPointsR);
int dumpTree(Node head);
void *mallocSafe(size_t size);

#endif
