cdef extern from "cluster_blas.h":

    void f_sgemm(char *transA, char *transB, int *m, int *n, int *k,
                 float *alpha, float *A, int *lda, float *B, int *ldb,
                 float *beta, float *C, int *ldc) nogil

    void f_dgemm(char *transA, char *transB, int *m, int *n, int *k,
                 double *alpha, double *A, int *lda, double *B, int *ldb,
                 double *beta, double *C, int *ldc) nogil

