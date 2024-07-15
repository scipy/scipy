#include <stdio.h>
#include "lapacke.h"

#define n 4

int main()
{
    int lda=n, ldb=n, ldvr=n, ldvl=n, info;
    char jobvl='V', jobvr='V';
    double alphar[n], alphai[n], beta[n];

    double vl[n*n], vr[n*n];
    // int lwork = 156;
    // double work[156];    /* cheat: 156 is the optimal lwork from the actual lapack call*/

    double a[n*n] = {12.0, 28.0, 76.0, 220.0,
                     16.0, 32.0, 80.0, 224.0,
                     24.0, 40.0, 88.0, 232.0,
                     40.0, 56.0, 104.0, 248.0};

    double b[n*n] = {2.0, 4.0, 10.0, 28.0,
                     3.0, 5.0, 11.0, 29.0,
                     5.0, 7.0, 13.0, 31.0,
                     9.0, 11.0, 17.0, 35.0};

    info = LAPACKE_dggev(LAPACK_ROW_MAJOR, jobvl, jobvr, n, a, lda, b,
                         ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr); //, work, lwork, info);

    printf("info = %d\n", info);

    printf("Re(eigv) = ");
    for(int i=0; i < n; i++){
        printf("%f , ", alphar[i] / beta[i] );
    }
    printf("\nIm(eigv = ");
    for(int i=0; i < n; i++){
        printf("%f , ", alphai[i] / beta[i] );
    }
    printf("\n");
}
