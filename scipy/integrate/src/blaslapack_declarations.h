#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

void dgbtrf_(int* m, int* n, int* kl, int* ku, double* ab, int* ldab, int* ipiv, int* info);
void dgbtrs_(char* trans, int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);


#endif
