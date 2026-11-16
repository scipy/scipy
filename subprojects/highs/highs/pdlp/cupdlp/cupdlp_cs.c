#include "pdlp/cupdlp/cupdlp_cs.h"

/* CSparse routines for reading inputs. Referenced from Tim Davis Suite Sparse
 */
void *cupdlp_dcs_malloc(int n, size_t size) {
  return (malloc(MAX(n, 1) * size));
}

/* wrapper for calloc */
void *cupdlp_dcs_calloc(int n, size_t size) {
  return (calloc(MAX(n, 1), size));
}

/* wrapper for free */
void *cupdlp_dcs_free(void *p) {
  if (p) free(p); /* free p if it is not already NULL */
  return (NULL);  /* return NULL to simplify the use of dcupdlp_dcs_free */
}

/* wrapper for realloc */
void *cupdlp_dcs_realloc(void *p, int n, size_t size, int *ok) {
  void *pnew;
  pnew = realloc(p, MAX(n, 1) * size); /* realloc the block */
  *ok = (pnew != NULL);                /* realloc fails if pnew is NULL */
  // return ((*ok) ? pnew : p);           /* return original p if failure */
  if (*ok)
    p = pnew;
  return p;
}

cupdlp_dcs *cupdlp_dcs_spalloc(int m, int n, int nzmax, int values,
                               int triplet) {
  cupdlp_dcs *A = cupdlp_dcs_calloc(
      1, sizeof(cupdlp_dcs)); /* allocate the cupdlp_dcs struct */
  if (!A) return (NULL);      /* out of memory */
  A->m = m;                   /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax = MAX(nzmax, 1);
  A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
  A->p = cupdlp_dcs_malloc(triplet ? nzmax : n + 1, sizeof(int));
  A->i = cupdlp_dcs_malloc(nzmax, sizeof(int));
  A->x = values ? cupdlp_dcs_malloc(nzmax, sizeof(double)) : NULL;
  return ((!A->p || !A->i || (values && !A->x)) ? cupdlp_dcs_spfree(A) : A);
}

/* change the max # of entries sparse matrix */
int cupdlp_dcs_sprealloc(cupdlp_dcs *A, int nzmax) {
  int ok, oki, okj = 1, okx = 1;
  if (!A) return (0);
  if (nzmax <= 0) nzmax = IS_CSC(A) ? (A->p[A->n]) : A->nz;
  nzmax = MAX(nzmax, 1);
  A->i = cupdlp_dcs_realloc(A->i, nzmax, sizeof(int), &oki);
  if (IS_TRIPLET(A)) A->p = cupdlp_dcs_realloc(A->p, nzmax, sizeof(int), &okj);
  if (A->x) A->x = cupdlp_dcs_realloc(A->x, nzmax, sizeof(double), &okx);
  ok = (oki && okj && okx);
  if (ok) A->nzmax = nzmax;
  return (ok);
}

/* free a sparse matrix */
cupdlp_dcs *cupdlp_dcs_spfree(cupdlp_dcs *A) {
  if (!A) return (NULL); /* do nothing if A already NULL */
  cupdlp_dcs_free(A->p);
  cupdlp_dcs_free(A->i);
  cupdlp_dcs_free(A->x);
  return ((cupdlp_dcs *)cupdlp_dcs_free(
      A)); /* free the cupdlp_dcs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
cupdlp_dcs *cupdlp_dcs_done(cupdlp_dcs *C, void *w, void *x, int ok) {
  cupdlp_dcs_free(w); /* free workspace */
  cupdlp_dcs_free(x);
  return (ok ? C
             : cupdlp_dcs_spfree(C)); /* return result if OK, else free it */
}

double cupdlp_dcs_cumsum(int *p, int *c, int n) {
  int i, nz = 0;
  double nz2 = 0;
  if (!p || !c) return (-1); /* check inputs */
  for (i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i]; /* also in double to avoid int overflow */
    c[i] = p[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p[n] = nz;
  return (nz2); /* return sum (c [0..n-1]) */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
int cupdlp_dcs_entry(cupdlp_dcs *T, int i, int j, double x) {
  if (!IS_TRIPLET(T) || i < 0 || j < 0) return (0); /* check inputs */
  if (T->nz >= T->nzmax && !cupdlp_dcs_sprealloc(T, 2 * (T->nzmax))) return (0);
  if (T->x) T->x[T->nz] = x;
  T->i[T->nz] = i;
  T->p[T->nz++] = j;
  T->m = MAX(T->m, i + 1);
  T->n = MAX(T->n, j + 1);
  return (1);
}

cupdlp_dcs *cupdlp_dcs_compress(const cupdlp_dcs *T) {
  int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  double *Cx, *Tx;
  cupdlp_dcs *C;
  if (!IS_TRIPLET(T)) return (NULL); /* check inputs */
  m = T->m;
  n = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C = cupdlp_dcs_spalloc(m, n, nz, Tx != NULL, 0);       /* allocate result */
  w = cupdlp_dcs_calloc(n, sizeof(int));                 /* get workspace */
  if (!C || !w) return (cupdlp_dcs_done(C, w, NULL, 0)); /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (k = 0; k < nz; k++) w[Tj[k]]++; /* column counts */
  cupdlp_dcs_cumsum(Cp, w, n);         /* column pointers */
  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    if (Cx) Cx[p] = Tx[k];
  }
  return (cupdlp_dcs_done(C, w, NULL, 1)); /* success; free w and return C */
}

double cupdlp_dcs_norm(const cupdlp_dcs *A) {
  int p, j, n, *Ap;
  double *Ax;
  double nrm = 0, s;
  if (!IS_CSC(A) || !A->x) return (-1); /* check inputs */
  n = A->n;
  Ap = A->p;
  Ax = A->x;
  for (j = 0; j < n; j++) {
    for (s = 0, p = Ap[j]; p < Ap[j + 1]; p++) s += fabs(Ax[p]);
    nrm = MAX(nrm, s);
  }
  return (nrm);
}

int cupdlp_dcs_print(const cupdlp_dcs *A, int brief) {
  int p, j, m, n, nzmax, nz, *Ap, *Ai;
  double *Ax;
  if (!A) {
    printf("(null)\n");
    return (0);
  }
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  nzmax = A->nzmax;
  nz = A->nz;
  if (nz < 0) {
    printf("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double)m, (double)n,
           (double)nzmax, (double)(Ap[n]), cupdlp_dcs_norm(A));
    for (j = 0; j < n; j++) {
      printf("    col %g : locations %g to %g\n", (double)j, (double)(Ap[j]),
             (double)(Ap[j + 1] - 1));
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        printf("      %g : ", (double)(Ai[p]));
        printf("%50.50e \n", Ax ? Ax[p] : 1);
        if (brief && p > 20) {
          printf("  ...\n");
          return (1);
        }
      }
    }
  } else {
    printf("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double)m, (double)n,
           (double)nzmax, (double)nz);
    for (p = 0; p < nz; p++) {
      printf("    %g %g : ", (double)(Ai[p]), (double)(Ap[p]));
      printf("%g\n", Ax ? Ax[p] : 1);
      if (brief && p > 20) {
        printf("  ...\n");
        return (1);
      }
    }
  }
  return (1);
}

cupdlp_dcs *cupdlp_dcs_transpose(const cupdlp_dcs *A, int values) {
  cupdlp_int p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w;
  double *Cx, *Ax;
  cupdlp_dcs *C;
  if (!IS_CSC(A)) return (NULL); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = cupdlp_dcs_spalloc(n, m, Ap[n], values && Ax, 0);  /* allocate result */
  w = cupdlp_calloc(m, sizeof(cupdlp_int));              /* get workspace */
  if (!C || !w) return (cupdlp_dcs_done(C, w, NULL, 0)); /* out of memory */
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  for (p = 0; p < Ap[n]; p++) w[Ai[p]]++; /* row counts */
  cupdlp_dcs_cumsum(Cp, w, m);            /* row pointers */
  for (j = 0; j < n; j++) {
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      Ci[q = w[Ai[p]]++] = j; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx[q] = Ax[p];
    }
  }
  return (cupdlp_dcs_done(C, w, NULL, 1)); /* success; free w and return C */
}
