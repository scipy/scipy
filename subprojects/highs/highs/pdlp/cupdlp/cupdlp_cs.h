#ifndef CUPDLP_CS_H
#define CUPDLP_CS_H

#include "pdlp/cupdlp/cupdlp_defs.h"

/* sparse matrix in column-oriented form used in reading mps*/
typedef struct cupdlp_cs_sparse {
  int nzmax;
  int m;     /* number of rows */
  int n;     /* number of columns */
  int *p;    /* column pointers (size n+1) or col indices (size nzmax) */
  int *i;    /* row indices, size nzmax */
  double *x; /* numerical values, size nzmax */
  int nz;    /* # of entries in triplet matrix, -1 for compressed-col */
} cupdlp_dcs;

int cupdlp_dcs_entry(cupdlp_dcs *T, int i, int j, double x);
cupdlp_dcs *cupdlp_dcs_compress(const cupdlp_dcs *T);
double cupdlp_dcs_norm(const cupdlp_dcs *A);
int cupdlp_dcs_print(const cupdlp_dcs *A, int brief);

/* utilities */
void *cupdlp_dcs_calloc(int n, size_t size);
void *cupdlp_dcs_free(void *p);
void *cupdlp_dcs_realloc(void *p, int n, size_t size, int *ok);
cupdlp_dcs *cupdlp_dcs_spalloc(int m, int n, int nzmax, int values, int t);
cupdlp_dcs *cupdlp_dcs_spfree(cupdlp_dcs *A);
int cupdlp_dcs_sprealloc(cupdlp_dcs *A, int nzmax);
void *cupdlp_dcs_malloc(int n, size_t size);

/* utilities */
double cupdlp_dcs_cumsum(int *p, int *c, int n);
cupdlp_dcs *cupdlp_dcs_done(cupdlp_dcs *C, void *w, void *x, int ok);
cupdlp_dcs *cupdlp_dcs_transpose(const cupdlp_dcs *A, int values);

#define IS_CSC(A) (A && (A->nz == -1))
#define IS_TRIPLET(A) (A && (A->nz >= 0))
/*--------------------------------------------------------------------------*/

#endif
