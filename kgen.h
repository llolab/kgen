// Copyright (c) 2019. Lloyd T. Elliott

#ifndef __KGEN_H_
#define __KGEN_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/sendfile.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef MKL_ILP64
  #include <mkl_lapack.h>
  #include <mkl_blas.h>
  #define INT MKL_INT
#else
  #define INT int
#endif

#define DOUBLE double

#define SWAP(x, y) do { __typeof__(x) SWAP = x; x = y; y = SWAP; } while (0)

#define error(s) { \
  fprintf(stderr, "error: %s (%s:%d).\n", s, __FILE__, __LINE__); \
  exit(-1); \
}

typedef struct {
  double *X;
  INT N;
  INT M;
  INT N0;
  INT M0;
} t_matrix;

typedef struct {
  int N;
  int *nodes;
  DOUBLE *edges;
} t_graph;

typedef struct {
  int M;
  int *e1;
  int *e2;
} t_edges;

typedef struct {
  int N;
  int *nodes;
} t_nodes;

static inline DOUBLE machine_epsilon(void) {
  DOUBLE max = 1.0, min = 0.0, test;
  int i;
  for (i = 0; i < 100000; i++) {
     DOUBLE one_plus_test;
     test = (max + min) / ((DOUBLE)2.0);
     one_plus_test = ((DOUBLE)1.0) + test;
     if (one_plus_test == ((DOUBLE)1.0)) {
        min = test;
     } else {
        max = test;
     }
  }

  return max;
}

static inline void destroy_edges(t_edges edges) {
  free(edges.e1);
  free(edges.e2);
}

static inline void destroy_graph(t_graph G) {
  free(G.nodes);
  free(G.edges);
}

static inline void destroy_nodes(t_nodes n) {
  free(n.nodes);
}

static inline void destroy_matrix(t_matrix A) {
  free(A.X);
}

#define destroy(X) _Generic((X), \
  t_matrix: destroy_matrix, \
  t_graph: destroy_graph, \
  t_nodes: destroy_nodes, \
  t_edges: destroy_edges) (X)

static inline t_matrix clone(t_matrix A) {
  t_matrix R;
  R.N = A.N;
  R.M = A.M;
  R.N0 = A.N0;
  R.M0 = A.M0;
  if ((R.X = (DOUBLE *)malloc(A.N0 * A.M0 * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }
  memcpy(R.X, A.X, A.N0 * A.M0 * sizeof(DOUBLE));
  return R;
}

static inline t_matrix uniform(INT N, INT M) {
  t_matrix A;
  A.N = N;
  A.M = M;
  A.N0 = N;
  A.M0 = M;
  if ((A.X = (DOUBLE *)malloc(A.N0 * A.M0 * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }

  for (INT i = 0; i < N * M; i++) {
    A.X[i] = ((double)rand())/((double)RAND_MAX);
  }

  return A;
}

static inline t_matrix create(INT N, INT M, INT N0, INT M0) {
  t_matrix A;
  A.N = N;
  A.M = M;
  A.N0 = N0;
  A.M0 = M0;
  DOUBLE *X = NULL;
  if ((X = (DOUBLE *)malloc(N0 * M0 * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }

  bzero(X, N0 * M0 * sizeof(DOUBLE));
  A.X = X;

  return A;
}

static inline t_matrix zeros(INT N, INT M) {
  return create(N, M, N, M);
}

static inline void put(t_matrix A, const double f, const int i, const int j) {
  A.X[i + j * A.N0] = f;
}

static inline DOUBLE get(t_matrix A, const int i, const int j) {
  return A.X[i + j * A.N0];
}

static inline void headers(char **V, int n, FILE *fd) {
  int first = 1;
  for (int j = 0; j < n; j++) {
    if (!first) {
      printf(" ");
    } else {
      first = 0;
    }
    fprintf(fd, "%s", V[j]);
  }
  fprintf(fd, "\n");
}

static inline void save(t_matrix A, FILE *fd) {
  for (int i = 0; i < A.N; i++) {
    int first = 1;
    for (int j = 0; j < A.M; j++) {
      if (!first) {
        fprintf(fd, " ");
      } else {
        first = 0;
      }
      fprintf(fd, "%f", (A.X)[i + j * (A.N0)]);
    }
    fprintf(fd, "\n");
  }
}

static inline void save_bin(t_matrix X, char *fname) {
  FILE *fp = fopen(fname, "wb");
  if (fp == NULL) {
    error("Could write to file");
  }

  for (int j = 0; j < X.M; j++) {
    ssize_t size = fwrite(&((X.X)[j * X.N0]), sizeof(DOUBLE), X.N, fp);
    if (size != X.N) {
      error("Could not write column");
    }
  }

  fclose(fp);
}

static inline void save_bin3(t_matrix X, char *fname) {
  FILE *fp = fopen(fname, "wb");
  if (fp == NULL) {
    error("Could not write to file");
  }

  for (int j = 0; j < X.M; j++) {
    ssize_t size = fwrite(&((X.X)[j * X.N0]), sizeof(DOUBLE), j + 1, fp);
    if (size != j + 1) {
      error("Could not write column");
    }
  }

  fclose(fp);
}

static inline INT find_size(char *fname) {
  FILE *fp = fopen(fname, "rb");
  if (fp == NULL) {
    error("Could not open file");
  }

  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  fclose(fp);
  size = size / sizeof(DOUBLE);
  return size;
}

static inline t_matrix load_bin2(char *fname, INT N, INT M, INT N0, INT M0) {
  t_matrix X = create(N, M, N0, M0);
  FILE *fp = fopen(fname, "rb");
  if (fp == NULL) {
    error("Could not open file");
  }

  for (int j = 0; j < X.M; j++) {
    ssize_t size = fread(&((X.X)[j * X.N0]), sizeof(DOUBLE), X.N, fp);
    if (size != X.N) {
      error("Could not read column");
    }
  }

  fclose(fp);
  return X;
}

static inline t_matrix load_bin3(char *fname, INT N, INT M, INT N0, INT M0) {
  t_matrix X = create(N, M, N0, M0);
  FILE *fp = fopen(fname, "rb");
  if (fp == NULL) {
    error("Could not open file");
  }

  for (int j = 0; j < X.M; j++) {
    ssize_t size = fread(&((X.X)[j * X.N0]), sizeof(DOUBLE), j + 1, fp);
    if (size != j + 1) {
      error("Could not read column");
    }
  }

  fclose(fp);
  return X;
}

static inline t_matrix load_bin(char *fname, INT P) {
  t_matrix X = create(P, P, P, P);
  FILE *fp = fopen(fname, "rb");
  if (fp == NULL) {
    error("Could not open file");
  }

  ssize_t size = fread(X.X, sizeof(DOUBLE), P * P, fp);
  if (size != P * P) {
    error("Could not read file");
  }

  fclose(fp);

  return X;
}

static inline t_matrix load(char *fname, int header) {
  FILE *fd = fopen(fname, "r");
  if (fd == NULL) {
    error("Could not open file for reading");
  }

  typedef struct {
    DOUBLE *xs;
    void *next;
  } t_line;

  t_line *first = NULL, *next;
  if ((first = (t_line *)malloc(sizeof(t_line))) == NULL) {
    error("memory");
  }
  int D = 1;
  int c = -1;
  while ((c = fgetc(fd)) != '\n') {
    if (c == EOF) {
      error("Input/output error");
    }

    if (c == ' ') {
      D += 1;
    }
  }

  if (D < 1) {
    error("Wrong arguments to load");
  }
  if ((first->xs = (DOUBLE *)malloc(D * sizeof(DOUBLE))) == NULL) {
    error("Memory limit");
  }
  first->next = NULL;
  DOUBLE f = NAN;
  next = first;
  int d = 0;
  int N = 0;

  if (header == 0) {
    fclose(fd);
    fd = fopen(fname, "r");
    if (fd == NULL) {
      error("Could not open file for reading");
    }
  }

  while (fscanf(fd, "%lf", &f) == 1) {
    if (d == D) {
      if ((next->next = (t_line *)malloc(sizeof(t_line))) == NULL) {
        error("Memory limit reached");
      }
      next = next->next;
      if ((next->xs = (DOUBLE *)malloc(D * sizeof(DOUBLE))) == NULL) {
        error("Memory limit reached");
      }
      next->next = NULL;
      d = 0;
    }

    next->xs[d] = f;
    d += 1;
    if (d == D) {
      N += 1;
    }
  }

  if (fscanf(fd, "%lf", &f) == 1) {
    error("Input/output error");
  }
  fclose(fd);
  if (d != D) {
    error("Input/output error");
  }
  DOUBLE *X;
  if ((X = (DOUBLE *)malloc(N * D * sizeof(DOUBLE))) == NULL) {
    error("memory");
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < D; j++) {
      X[i + j * N] = first->xs[j];
    }

    next = first->next;
    free(first->xs); 
    first->xs = NULL;
    free(first);
    first = next;
  }

  if (N == 0) {
    error("i/o");
  }

  t_matrix result;
  result.X = X;
  result.N = N;
  result.M = D;
  result.N0 = N;
  result.M0 = D;
  return result;
}

static inline void print_matrix(t_matrix A) {
  for (int i = 0; i < A.N; i++) {
    int first = 1;
    for (int j = 0; j < A.M; j++) {
      if (!first) {
        printf(" ");
      } else {
        first = 0;
      }
      printf("%f", (A.X)[i + j * (A.N0)]);
    }
    printf("\n");
  }
}

static inline void print_graph(t_graph G) {
  int N = G.N;
  for (int i = 1; i < N; i++) {
    for (int j = 0; j < i; j++) {
      if (fabs(G.edges[i + N * j]) > 0+1e-15) {
        printf("%d %d %f\n", i, j, G.edges[i + N * j]);
      }
    }
  }
}

#define print(X) _Generic((X), \
  t_matrix: print_matrix, \
  t_graph: print_graph)(X)

void dch1up_(INT *, DOUBLE *, INT *, DOUBLE *, DOUBLE *);

static inline void udate0(t_matrix L, t_matrix V) {
  DOUBLE *w = NULL;
  if ((w = (DOUBLE *)malloc(L.N * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }

  bzero(w, L.N * sizeof(DOUBLE));
  dch1up_(&L.N, L.X, &L.N0, V.X, w);
  free(w);
}

static inline DOUBLE norm(t_matrix A) {
  DOUBLE z = 0.0;
  for (INT i = 0; i < A.N; i++) {
    for (INT j = 0; j < A.M; j++) {
      DOUBLE x = get(A, i, j);
      x = x * x;
      z += x;
    }
  }
  return sqrt(z);
}

void dch1dn_(INT *, DOUBLE *, INT *, DOUBLE *, DOUBLE *, INT *);

static inline void ddate0(t_matrix L, t_matrix V) {
  DOUBLE *w = NULL;
  if ((w = (DOUBLE *)malloc(L.N * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }

  bzero(w, L.N * sizeof(DOUBLE));
  INT info = 0;
  dch1dn_(&L.N, L.X, &L.N0, V.X, w, &info);
  if (info != 0) {
    error("Numerical issue");
  }
  free(w);
}

void dpotrf_(const char *, const INT *, DOUBLE *, const INT *, INT *);

static inline void chol(t_matrix A) {
  if (A.N != A.M) {
    error("Matrix not square");
  }
  
  if (A.N == 0) {
    return;
  }

  INT info = 0;
  char transA = 'U';
  dpotrf_(&transA, &A.N, A.X, &A.N0, &info);
  for (INT i = 0; i < A.M; i++) {
    bzero(&(A.X[i * A.N0 + i + 1]), (A.N - i - 1) * sizeof(DOUBLE));
  }

  if (info != 0) {
    error("Numerical issue");
  }
}

void dgemv_(char *,
  INT *,
  INT *,
  DOUBLE *,
  DOUBLE *,
  INT *,
  DOUBLE *,
  INT *,
  DOUBLE *, 
  DOUBLE *,
  INT *);

void dgemm_(char *,
  char *,
  INT *,
  INT *,
  INT *,
  DOUBLE *,
  DOUBLE *,
  INT *,
  DOUBLE *,
  INT *,
  DOUBLE *,
  DOUBLE *,
  INT *);

static inline void mult(t_matrix A, t_matrix V, t_matrix *R) {
  if (A.M != V.N) {
    error("Dimension mismatch");
  }
  
  if (R->N0 < A.N) {
    error("Dimension mismatch");
  }

  R->N = A.N;
  
  if (R->M0 < V.M) {
    error("Dimension mismatch");
  }

  R->M = V.M;
  
  bzero(R->X, (R->N0 * R->M0) * sizeof(DOUBLE));
  DOUBLE alpha = 1.0;
  DOUBLE beta = 0.0;
  char transA = 'N';

  if (V.M == 1) { 
    INT i = 1;
    dgemv_(&transA,
      &A.M,
      &A.N,
      &alpha,
      A.X,
      &A.M0,
      V.X,
      &i,
      &beta,
      R->X,
      &i);

  } else {
    dgemm_(&transA,
      &transA,
      &A.N,
      &V.M,
      &A.M,
      &alpha,
      A.X,
      &A.N0,
      V.X,
      &V.N0,
      &beta,
      R->X,
      &R->N0);

  }
}

static inline void square(t_matrix A, int dim, t_matrix *R) {
  DOUBLE alpha = 1.0;
  DOUBLE beta = 0.0;
  char transA;
  char transB;
  INT AM = 0, AN = 0, AK = 0;

  switch (dim) {
    case 0:
      if (R->N0 < A.M) {
        error("Dimension mismatch");
      }

      R->N = A.M;
      if (R->M0 < A.M) {
        error("Dimension mismatch");
      }

      R->M = A.M;
      transA = 'T';
      transB = 'N';
      AM = A.M;
      AN = A.M;
      AK = A.N;
      break;

    case 1:
      if (R->N0 < A.N) {
        error("Dimension mismatch");
      }

      R->N = A.N;
      if (R->M0 < A.N) {
        error("Dimension mismatch");
      }

      R->M = A.N;
      transA = 'N';
      transB = 'T';
      AM = A.N;
      AN = A.N;
      AK = A.M;
      break;

    default:
      error("Unknown dimension");
  }
  
  bzero(R->X, (R->N0 * R->M0) * sizeof(DOUBLE));

  dgemm_(&transA,
    &transB,
    &AM,
    &AN,
    &AK,
    &alpha,
    A.X,
    &A.N0,
    A.X,
    &A.N0,
    &beta,
    R->X,
    &R->N0);

}

static inline void scale(t_matrix A, DOUBLE alpha, t_matrix *C) {
  if (C->N0 < A.N) {
    error("Dimension mismatch");
  }

  if (C->M0 < A.M) {
    error("Dimension mismatch");
  }

  C->N = A.N;
  C->M = A.M;

  for (INT i = 0; i < A.N; i++) {
    for (INT j = 0; j < A.M; j++) {
      put(*C, alpha * get(A, i, j), i, j);
    }                   
  }
}

static inline void sub(t_matrix A, t_matrix B, t_matrix *C) {
  if (C->N0 < A.N) {
    error("Dimension mismatch");
  }

  if (C->M0 < A.M) {
    error("Dimension mismatch");
  }

  if (A.N != B.N) {
    error("Dimension mismatch");
  }

  if (A.M != B.M) {
    error("Dimension mismatch");
  }

  C->N = A.N;
  C->M = A.M;

  for (INT i = 0; i < A.N; i++) {
    for (INT j = 0; j < A.M; j++) {
      put(*C, get(A, i, j) - get(B, i, j), i, j);
    }                   
  }
}

static inline int iszero(DOUBLE x) {
  return fabs(x) < 1e-12;
}

static inline void is_nan(t_matrix A) {
  for (int j = 0; j < A.M; j++) {
    for (int i = 0; i < A.N; i++) {
      int x = isnan(get(A, i, j));
      put(A, x, i, j);
    }
  }
}

static inline void add(t_matrix A, t_matrix B, t_matrix *C) {
  if (C->N0 < A.N) {
    error("Dimension mismatch");
  }

  if (C->M0 < A.M) {
    error("Dimension mismatch");
  }

  if (A.N != B.N) {
    error("Dimension mismatch");
  }

  if (A.M != B.M) {
    error("Dimension mismatch");
  }

  C->N = A.N;
  C->M = A.M;
  for (INT i = 0; i < A.N; i++) {
    for (INT j = 0; j < A.M; j++) {
      put(*C, get(A, i, j) + get(B, i, j), i, j);
    }                   
  }
}

void dtrsv_(const char *,
  const char *,
  const char *,
  const INT *,
  const DOUBLE *,
  const INT *,
  DOUBLE *,
  const INT *);

static inline void solve(t_matrix L, t_matrix V) {
  char uplo = 'U';
  char transA = 'N';
  char diag = 'N'; 
  INT incx = 1;
  dtrsv_(&uplo, &transA, &diag, &L.N, L.X, &L.N0, V.X, &incx);
}

static inline void remove_square(t_matrix *A, INT n, INT *which) {
  INT N = A->N;
  INT M = A->M;
  if (N != M) {
    error("Dimension mismatch");
  }

  INT k = 0;
  INT j = 0;
  if (n == 0 || which == NULL) {
    return;
  }

  // TODO: Optimize this
  for (INT i = 0; i < M; i++) { 
    if (i == which[k]) {
      k++; 
      assert(k <= n);
    } else {
      for (INT l = 0; l < N; l++) {
        put(*A, get(*A, i, l), j, l);
      }
      j++;
    }
  }
  
  k = 0;
  j = 0;
  for (INT i = 0; i < M; i++) { 
    if (i == which[k]) {
      k++; 
      assert(k <= n);
    } else {
      for (INT l = 0; l < N; l++) {
        put(*A, get(*A, l, i), l, j);
      }

      j++;
    }
  }

  A->N = N - n;
  A->M = M - n;
}

static inline void delete(t_matrix *A, INT d) {
  if (A->N != A->M) {
    error("Dimension mismatch");
  }

  if (A->N == 0 || A->M == 0) {
    error("Cannot delete further");
  }

  A->N -= 1;
  A->M -= 1;
  if (d == A->N) {
    return;
  }

  for (int j = d; j < A->M; j++) {
    memcpy(&(A->X[j * A->N0]),
      &(A->X[(j + 1) * A->N0]),
      (A->N + 1) * sizeof(DOUBLE));

  }

  for (int j = 0; j < A->M; j++) {
    memmove(&(A->X[d + j * A->N0]),
      &(A->X[(d + 1) + j * A->N0]),
      (A->N - d) * sizeof(DOUBLE));

  }
}

static inline void expand(t_matrix *A, t_matrix V, INT d) {
  A->N += 1;
  A->M += 1;
  if (V.N > V.N0) {
    error("Cannot expand further");
  }
  if (V.M > V.M0) {
    error("Cannot expand further");
  }
  if (V.N != A->N) {
    error("Dimension mismatch");
  }
  if (V.M != 1) {
    error("Dimension mismatch");
  }
  if (A->N != A->M) {
    error("Dimension mismatch");
  }
  for (int j = A->M - 1; j > d; j--) {
    memcpy(&(A->X[j * A->N0]),
      &(A->X[(j - 1) * A->N0]),
      (A->N - 1) * sizeof(DOUBLE));

  }

  for (int j = 0; j < A->M; j++) {
    memmove(&(A->X[(d+1) + j * A->N0]),
      &(A->X[d + j * A->N0]),
      (V.N - d - 1) * sizeof(DOUBLE));

  }
  
  memcpy(&(A->X[d * A->N0]), V.X, V.N * sizeof(DOUBLE));
  for (int i = 0; i < V.N; i++) {
    put(*A, get(V, i, 0), d, i);
  }
}

static inline void udate(t_matrix *L, INT dim, t_matrix V) {
  INT D = L->N;
  if (V.N != L->N + 1) {
    error("Dimension mismatch");
  }
  if (V.N != L->M + 1) {
    error("Dimension mismatch");
  }
  if (V.M != 1) {
    error("Dimension mismatch");
  }
  t_matrix Z = zeros(D + 1, 1);
  expand(L, Z, dim);
  t_matrix v = clone(V);
  put(v, get(v, dim, 0)/2, dim, 0);
  t_matrix u = zeros(D + 1, 1);
  put(u, 1.0, dim, 0);
  DOUBLE vn = norm(v);
  t_matrix x = clone(u);
  t_matrix y = clone(u);
  scale(v, 1.0/vn, &v);
  add(x, v, &x);
  sub(y, v, &y);
  scale(x, sqrt(vn/2.0), &x);
  scale(y, sqrt(vn/2.0), &y);
  udate0(*L, x);
  ddate0(*L, y);
  destroy(Z);
  destroy(v);
  destroy(u);
  destroy(x);
  destroy(y);
}

static inline t_matrix submat(t_matrix A, INT x1, INT y1, INT x2, INT y2) {
  if (x2 == -1) {
    x2 = A.N;
  }

  if (y2 == -1) {
    y2 = A.M;
  }

  INT N = x2 - x1;
  INT M = y2 - y1;
  t_matrix R = zeros(N, M);
  for (INT j = 0; j < M; j++) {
    memcpy(&(R.X[j * R.N0]),
      &(A.X[x1 + (j + y1) * A.N0]),
      N * sizeof(DOUBLE));

  }
  return R;
}

static inline void trans(t_matrix *A) {
  if (A->N == 0 || A->M == 0) {
    error("Cannot transpose empty matrix");
  }

  if (A->N0 == 1 || A->M0 == 1) {
    INT tmp = A->M;
    A->M = A->N;
    A->N = tmp;
    tmp = A->M0;
    A->M0 = A->N0;
    A->N0 = tmp;
    return;
  }
   
  // TODO: Optimize this
  t_matrix B = clone(*A);
  INT tmp = A->M;
  A->M = A->N;
  A->N = tmp;
  tmp = A->M0;
  A->M0 = A->N0;
  A->N0 = tmp;
  for (INT i = 0; i < B.N; i++) {
    for (INT j = 0; j < B.M; j++) {
      put(*A, get(B, i, j), j, i);
    }
  }
  destroy(B);
}

static inline void copy(t_matrix A,
  INT x1,
  INT y1,
  INT x2,
  INT y2,
  t_matrix B) {

  if (x2 == -1) {
    x2 = B.N;
  }
  if (y2 == -1) {
    y2 = B.M;
  }
  for (INT j = y1; j < y2; j++) {
    memcpy(&(B.X[(x1 + j * B.N0)]),
      &(A.X[(j - y1) * A.N0]),
      (x2 - x1) * sizeof(DOUBLE));

  }
}

static inline void ddate(t_matrix *L, INT dim) {
  INT D = L->N;
  if (dim >= D) {
    error("Dimension too high.");
  }
  if (dim == D - 1) {
    delete(L, dim);
    return;
  }

  t_matrix A = submat(*L, dim + 1, dim + 1, -1, -1); 
  t_matrix v = submat(*L, dim, dim + 1, dim + 1,  -1);
  trans(&v);
  udate0(A, v);
  copy(A, dim + 1, dim + 1, -1, -1, *L);
  delete(L, dim);
  destroy(A);
  destroy(v);
}

static inline t_matrix vectorize(t_matrix A) {
  t_matrix R;
  R.M = 1;
  R.N = A.N * A.M;
  R.M0 = 1;
  R.N0 = A.N * A.M;
  if ((R.X = (DOUBLE *)malloc(A.N * A.M * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }
  INT k = 0;

  // TODO: Optimize this
  for (INT j = 0; j < A.M; j++) {
    for (INT i = 0; i < A.N; i++) {
      put(R, get(A, i, j), k, 0);
      k++;
    }
  }
  return R;
}

static inline DOUBLE diff(t_matrix A, t_matrix B) {
  if (A.N < 1) {
    error("Difference on empty matrix");
  }
  if (A.M < 1) {
    error("Difference on empty matrix");
  }
  if (A.N != B.N || A.M != B.M) {
    error("Difference on mismatched matrices");
  }
  DOUBLE f = fabs(get(A, 0, 0) - get(B, 0, 0));
  for (INT j = 0; j < A.M; j++) {
    for (INT i = 0; i < A.N; i++) {
      DOUBLE x = fabs(get(A, i, j) - get(B, i, j));
      if (x > f) {
        f = x;
      }
    }
  }
  return f;
}

static inline t_graph create_graph_matrix(t_matrix A) {
  INT N = A.N;
  int n= N;
  t_graph G;
  G.N = n;
  if (A.N != A.M) {
    error("Adjacency matrix must be square");
  }
  if ((G.nodes = (int *)malloc(N * sizeof(int))) == NULL) {
    error("Memory exhausted");
  }
  for (int i = 0; i < n; i++) {
    G.nodes[i] = 1;
  }
  if ((G.edges = (DOUBLE *)malloc(n * n * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }
  bzero(G.edges, n * n * sizeof(DOUBLE));
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      DOUBLE x = get(A, i, j);
      G.edges[i + n * j] = x; 
      G.edges[j + n * i] = x; 
    }
  }
  return G;
}

static inline t_graph create_graph_int(int N) {
  t_graph G;
  G.N = N;
  if ((G.nodes = (int *)malloc(N * sizeof(int))) == NULL) {
    error("Memory exhausted");
  }

  bzero(G.nodes, N * sizeof(int));
  if ((G.edges = (DOUBLE *)malloc(N * N * sizeof(DOUBLE))) == NULL) {
    error("Memory exhausted");
  }

  bzero(G.edges, N * N * sizeof(DOUBLE));
  return G;
}

#define create_graph(X) _Generic((X), \
  int: create_graph_int, \
  t_matrix: create_graph_matrix)(X)

static inline void add_node(t_graph G, int i) {
  G.nodes[i] = 1;
}


// TODO: Optimize graph operations
static inline void rm_node(t_graph G, int i) {
  G.nodes[i] = 0;
  for (int j = 0; j < G.N; j++) {
    G.edges[i + j * G.N] = 0.0; 
    G.edges[j + i * G.N] = 0.0; 
  }
}

static inline void add_edge(t_graph G, int i, int j, DOUBLE w) {
  G.nodes[i] = 1;
  G.nodes[j] = 1;
  G.edges[i + j * G.N] = w;
  G.edges[j + i * G.N] = w;
}

static inline void remove_edge(t_graph G, int i, int j) {
  G.edges[i + j * G.N] = 0.0;
  G.edges[j + i * G.N] = 0.0;
}

static inline DOUBLE get_weight(t_graph G, int i, int j) {
  return G.edges[i + j * G.N];
}

static inline t_nodes neighbours(t_graph G, int i) {
  int N = 0;
  for (int j = 0; j < G.N; j++) {
    if (G.edges[i + j * G.N] > 1e-12 && i != j) {
      N += 1;
    }
  }
  t_nodes result;
  result.N = N;
  if ((result.nodes = (int *)malloc(N * sizeof(int))) == NULL) {
    error("Memory exhausted");
  }

  int k = 0;
  for (int j = 0; j < G.N; j++) {
    if (G.edges[i + j * G.N] > 1e-12 && i != j) {
      result.nodes[k] = j;
      k++;
    }
  }

  if (k != N) {
    error("Math error");
  }
  return result;
}


static inline t_edges get_edges(t_graph G) {
  int M = 0;
  for (int i = 1; i < G.N; i++) {
    for (int j = 0; j < i; j++) {
      if (G.edges[i + j * G.N] > 0) {
        M++;
      }
    }
  }

  t_edges result;
  result.M = M;
  result.e1 = NULL;
  if ((result.e1 = (int *)malloc(M * sizeof(int))) == NULL) {
    error("Memory exhausted");
  }

  result.e2 = NULL;
  if ((result.e2 = (int *)malloc(M * sizeof(int))) == NULL) {
    error("Memory exhausted");
  }

  M = 0;
  for (int i = 1; i < G.N; i++) {
    for (int j = 0; j < i; j++) {
      if (G.edges[i + j * G.N] > 0) {
        result.e1[M] = i;
        result.e2[M] = j;
        M++;
      }
    }
  }

  return result;
}

// Kruskal algorithm
static inline void minspan(t_graph G) {
  int N = G.N;
  int *indices = NULL;
  int *n1 = NULL;
  int *n2 = NULL;
  DOUBLE *w = NULL;
  int *vertices = NULL;
  int size = (N*(N-1))/2;
  if ((indices = ((int *)malloc(size * sizeof(int)))) == NULL) {
    error("Memory exhausted");
  }
  if ((n1 = ((int *)malloc(size * sizeof(int)))) == NULL) {
    error("Memory exhausted");
  }
  if ((n2 = ((int *)malloc(size * sizeof(int)))) == NULL) {
    error("Memory exhausted");
  }
  if ((w = ((DOUBLE *)malloc(size * sizeof(DOUBLE)))) == NULL) {
    error("Memory exhausted");
  }
  if ((vertices = ((int *)malloc(N * sizeof(int)))) == NULL) {
    error("Memory exhausted");
  }
  int k = 0;
  for (int i = 1; i < N; i++) {
    for (int j = 0; j < i; j++) {
      indices[k] = k;
      n1[k] = i;
      n2[k] = j;
      w[k] = get_weight(G, i, j); 
      k++;
    }
  }

  if (k != size) {
    error("Math error");
  }
  int compare(const void *A, const void *B) {
    int a = *((int *)A);
    int b = *((int *)B);
    DOUBLE w1 = w[a];
    DOUBLE w2 = w[b];
    if (w1 < w2) {
      return -1;
    } else if (w1 > w2) {
      return 1;
    } else {
      return 0;
    }
  }
  
  bzero(G.edges, N * N * sizeof(DOUBLE));
  qsort(indices, size, sizeof(int), compare);
  for (int i = 0; i < N; i++) {
    vertices[i] = i;
  }
  int i = 0;
  for (int j = 0; j < size; j++) {
    if (vertices[n1[indices[j]]] != vertices[n2[indices[j]]]) {
      add_edge(G, n1[indices[j]], n2[indices[j]], w[indices[j]]);
      int z = vertices[n1[indices[j]]];
      i++;
      for (int k = 0; k < N; k++) {
        if (vertices[k] == z) {
          vertices[k] = vertices[n2[indices[j]]];
        }
      }

      if (i == N - 1) {
        break;
      }
    }
  }

  if (i != N - 1) {
    error("Minimum spanning tree error");
  }

  free(indices);
  free(n1);
  free(n2);
  free(w);
  free(vertices);
}

static inline void cp(char *src, char *dst) {
  FILE *fp = fopen(src, "rb");
  if (fp == NULL) {
    error("Could not open file");
  }
  fseek(fp, 0L, SEEK_END);
  ssize_t size = ftell(fp);
  fclose(fp);
  int fd1 = open(src, O_RDONLY);
  if (fd1 < 0) {
    error("Could not open file for reading");
  }
  mode_t mode = S_IRUSR | S_IWUSR;
  int fd2 = open(dst, O_WRONLY | O_CREAT, mode);
  if (fd2 < 0) {
    error("Could not open file for writing");
  }
  ssize_t retval = sendfile(fd2, fd1, NULL, size);
  if (retval != size) {
    error("Could not copy file");
  }
  close(fd1);
  close(fd2);
}

#endif
