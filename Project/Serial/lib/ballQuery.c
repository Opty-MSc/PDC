#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double radius;
  long id;
  long L;
  long R;
} Node;

typedef struct hash {
  long id;
  long index;
  struct hash *next;
} Hash;

int nDims;
long nNodes;
double *P;

Node *tree;
double **center;
Hash **hash;

long currBest;
double minD = 1000000.0;

void allocateHash() {
  hash = (Hash **) calloc(nNodes, sizeof(Hash *));
  if (hash == NULL) {
    fprintf(stderr, "FATAL: [calloc]!\n");
    exit(EXIT_FAILURE);
  }
}

void hashInsert(long treeId, long arrayIdx) {
  long i = treeId % nNodes;

  Hash *new = (Hash *) malloc(sizeof(Hash));
  new->id = treeId;
  new->index = arrayIdx;
  new->next = hash[i];
  hash[i] = new;
}

long hashGetIndex(long id) {
  long i = id % nNodes;

  Hash *h;
  for (h = hash[i]; h != NULL && h->id != id; h = h->next)
    ;

  if (h == NULL) {
    fprintf(stderr, "NOT FOUND: Id %ld!\n", id);
    exit(EXIT_FAILURE);
  }
  return h->index;
}

void allocateTree() {
  tree = (Node *) malloc(nNodes * sizeof(Node));
  double *centerValues = (double *) malloc(nNodes * nDims * sizeof(double));
  center = (double **) malloc(nNodes * sizeof(double *));

  if (centerValues == NULL || center == NULL || tree == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }
  for (long i = 0; i < nNodes; i++)
    center[i] = &centerValues[i * nDims];
}

double distance(double *P1, double *P2) {
  double distance = 0;
  for (int i = 0; i < nDims; i++)
    distance += (P2[i] - P1[i]) * (P2[i] - P1[i]);
  return sqrt(distance);
}

void search_tree(long idx) {
  if (tree[idx].radius == 0.0) {
    double d = distance(center[idx], P);
    if (d < minD) {
      minD = d;
      currBest = idx;
    }
    return;
  }

  long idxL = hashGetIndex(tree[idx].L);
  if (distance(center[idxL], P) - tree[idx].radius < minD)
    search_tree(idxL);
  long idxR = hashGetIndex(tree[idx].R);
  if (distance(center[idxR], P) - tree[idx].radius < minD)
    search_tree(idxR);
}

int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr, "[%s] Usage: [Ball Tree File] [Point]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  FILE *fp = fopen(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "Usage: [fopen] (%s)!\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  if (fscanf(fp, "%d %ld", &nDims, &nNodes) != 2) {
    fprintf(stderr, "FATAL: [fscanf]!\n");
    exit(EXIT_FAILURE);
  }
  if (nDims < 2) {
    fprintf(stderr, "ILLEGAL: Number of Dimensions (%d)!\n", nDims);
    exit(EXIT_FAILURE);
  }
  if (nNodes < 2) {
    fprintf(stderr, "ILLEGAL: Number of Nodes (%ld)!\n", nNodes);
    exit(EXIT_FAILURE);
  }

  if (argc != nDims + 2) {
    fprintf(stderr, "ILLEGAL: Wrong Number of Coordinates!\n");
    exit(EXIT_FAILURE);
  }

  P = (double *) malloc(nDims * sizeof(double));
  if (P == NULL) {
    fprintf(stderr, "FATAL: [malloc]!\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < nDims; i++)
    P[i] = atof(argv[i + 2]);

  allocateHash();
  allocateTree();

  long nodeIdx;
  for (long i = 0; i < nNodes; i++) {
    if (fscanf(fp, "%ld", &nodeIdx) != 1) {
      fprintf(stderr, "FATAL: [fscanf]!\n");
      exit(EXIT_FAILURE);
    }
    hashInsert(nodeIdx, i);
    Node *node = &(tree[i]);
    if (fscanf(fp, "%ld %ld %lf", &(node->L), &(node->R), &(node->radius)) != 3) {
      fprintf(stderr, "FATAL: [fscanf]!\n");
      exit(EXIT_FAILURE);
    }
    for (int j = 0; j < nDims; j++)
      if (fscanf(fp, "%lf", &(center[i][j])) != 1) {
        fprintf(stderr, "FATAL: [fscanf]!\n");
        exit(EXIT_FAILURE);
      }
  }

  search_tree(hashGetIndex(0));

  for (int i = 0; i < nDims; i++)
    printf("%lf ", center[currBest][i]);
  printf("\n");
}
