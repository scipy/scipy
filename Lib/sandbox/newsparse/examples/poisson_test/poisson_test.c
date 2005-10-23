/* zuse:    cc -o poisson_test -fast -native poisson_test.c mmio.c pcg.c -lm -xlic_lib=sunperf */
/* ru-lt13: cc -O poisson_test.c mmio.c pcg.c -o poisson_test -lblas -lg2c -lm */
#include <assert.h>
#include <stdio.h>
#include "mmio.h"

/* matrix A */
static int n_s;
static double *va_s, *da_s;
static int *ja_s, *ia_s;

/* CONVERT_COO_SSS - convert sparse matrix from COO to SSS format
 */
void convert_COO_SSS(int n, int nz,
                     int *i_coo, int *j_coo, double *v_coo,
                     int **ia, int **ja, double **va, double **da) {
  int i, k, l, t, nnz;
  int *root;

  root = (int *)malloc(n * sizeof(int));
  assert(root);

  /* allocate SSS matrix structure (1st part) */
  (*da) = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i ++) {
    root[i] = 0;
    (*da)[i] = 0.0;
  }

  /* build n linked lists */
  nnz = 0;
  for (k = 0; k < nz; k ++) {
    if (i_coo[k] == j_coo[k])
      /* diagonal element */
      (*da)[i_coo[k]] = v_coo[k];
    else {
      /* off diagonal element */
      if (i_coo[k] < j_coo[k]) {   /* move to lower triangle */
        t = i_coo[k];
        i_coo[k] = j_coo[k];
        j_coo[k] = t;
      }
      i = i_coo[k];                /* link */
      i_coo[k] = root[i];
      root[i] = k;
      nnz ++;
    }
  }

  /* allocate SSS matrix structure (2nd part) */
  (*ia) = (int *)malloc((n+1) * sizeof(int));
  (*va) = (double *)malloc((nnz * sizeof(double)));
  (*ja) = (int *)malloc(nnz * sizeof(int));

  /* fill SSS matrix structure */
  k = 0;
  for (i = 0; i < n; i ++) {
    (*ia)[i] = k;
    l = root[i];
    while (l != 0) {
      (*ja)[k] = j_coo[l];
      (*va)[k] = v_coo[l];
      k ++;
      l = i_coo[l];
    }
  }
  (*ia)[n] = k;
  assert(k == nnz);
  free(root);
}

/* READ_MTX - read symmetric sparse matrix in MatrixMarket format
 */
void read_MTX_SSS(char *fname, int *n,
                  double **va, double **da, int **ja, int **ia) {
  int m, nz, ret_code, i;
  double *v_coo;
  int *i_coo, *j_coo;
  MM_typecode matcode;
  FILE *f;

  f = fopen(fname, "r");
  assert(f != NULL);
  ret_code = mm_read_banner(f, &matcode);
  assert(ret_code == 0);
  assert(mm_is_real(matcode) && mm_is_matrix(matcode) &&
         mm_is_sparse(matcode) && mm_is_symmetric(matcode));
  ret_code = mm_read_mtx_crd_size(f, &m, n, &nz);
  assert(ret_code == 0);
  assert(m == *n);
  /* read COO format */
  i_coo = (int *)malloc(nz * sizeof(int));
  j_coo = (int *)malloc(nz * sizeof(int));
  v_coo = (double *)malloc(nz * sizeof(double));
  assert(i_coo && j_coo && v_coo);
  for (i = 0; i < nz; i ++) {
    fscanf(f, "%d %d %lg\n", &i_coo[i], &j_coo[i], &v_coo[i]);
    i_coo[i]--;  /* adjust from 1-based to 0-based */
    j_coo[i]--;
  }
  fclose(f);
  /* convert to SSS format */
  convert_COO_SSS(*n, nz, i_coo, j_coo, v_coo, ia, ja, va, da);
  free(i_coo); free(j_coo); free(v_coo);
}

/* MATVEC - matrix vector multiplications
 */
void matvec(double *x, double *y) {
  double s, v, xi;
  int i, j, k;
 
  for (i = 0; i < n_s; i ++) {
    xi = x[i];
    s = 0.0;
    for (k = ia_s[i]; k < ia_s[i+1]; k ++) {
      j = ja_s[k];
      v = va_s[k];
      s += v * x[j];
      y[j] += v * xi;
    }
    y[i] = s + da_s[i]*xi;
  }
}

void main () {
  double *x, *b, *work;
  int i;
  double relres;
  int iter, flag;

  read_MTX_SSS("matrices/poi2d_100.mtx", &n_s, &va_s, &da_s, &ja_s, &ia_s);

  x = (double *) malloc(n_s * sizeof(double));
  b = (double *) malloc(n_s * sizeof(double));
  work = (double *) malloc(4*n_s * sizeof(double));
  assert(x != NULL && b != NULL && work != NULL);

  for (i = 0; i < n_s; i ++) {
    x[i] = 0.0;
    b[i] = 1.0;
  }

  printf("Starting PCG solver...\n");
  pcg(n_s, x, b, 1e-12, 2000, 1, &iter, &relres, &flag, work, matvec, NULL);

}
