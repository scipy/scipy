#ifndef RELAXATION_H
#define RELAXATION_H

#include <assert.h>
#include <iostream>

template<class I, class T>
void gauss_seidel(const I Ap[], 
                  const I Aj[], 
                  const T Ax[],
                        T  x[],
                  const T  b[],
                  const I row_start,
                  const I row_stop,
                  const I row_step)
{
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];
        T rsum = 0;
        T diag = 0;
        
        for(I jj = start; jj < end; jj++){
            I j = Aj[jj];
            if (i == j)
                diag  = Ax[jj];
            else
                rsum += Ax[jj]*x[j];
        }
        
        //TODO raise error? inform user?
        if (diag != 0){
            x[i] = (b[i] - rsum)/diag;
        }
    }
}


template<class I, class T>
void block_gauss_seidel(const I Ap[], 
                        const I Aj[], 
                        const T Ax[],
                              T  x[],
                        const T  b[],
                        const I row_start,
                        const I row_stop,
                        const I row_step,
                        const I blocksize)
{
    const I B2 = blocksize * blocksize;
    
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];

        for(I bi = 0; bi < blocksize; bi++){
            T rsum = 0;
            T diag = 0;

            for(I jj = start; jj < end; jj++){
                I j = Aj[jj];
                const T * block_row = Ax + B2*jj + blocksize*bi;
                const T * block_x   = x + blocksize * j;

                if (i == j){
                    //diagonal block
                    diag = block_row[bi];
                    for(I bj = 0; bj < bi; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                    for(I bj = bi+1; bj < blocksize; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                } else {
                    for(I bj = 0; bj < blocksize; bj++){
                        rsum += block_row[bj] * block_x[bj];
                    }
                }
            }

            //TODO raise error? inform user?
            if (diag != 0){
                x[blocksize*i + bi] = (b[blocksize*i + bi] - rsum)/diag;
            }
        }
    }
}

template<class I, class T>
void jacobi(const I Ap[], 
            const I Aj[], 
            const T Ax[],
                  T  x[],
            const T  b[],
                  T temp[],
            const I row_start,
            const I row_stop,
            const I row_step,
            const T omega)
{
    for(I i = row_start; i != row_stop; i += row_step) {
        temp[i] = x[i];
    }
    
    for(I i = row_start; i != row_stop; i += row_step) {
        I start = Ap[i];
        I end   = Ap[i+1];
        T rsum = 0;
        T diag = 0;
        
        for(I jj = start; jj < end; jj++){
            I j = Aj[jj];
            if (i == j)
                diag  = Ax[jj];
            else
                rsum += Ax[jj]*temp[j];
        }
        
        //TODO raise error? inform user?
        if (diag != 0){ 
            x[i] = (1 - omega) * temp[i] + omega * ((b[i] - rsum)/diag);
        }
    }
}


#endif

