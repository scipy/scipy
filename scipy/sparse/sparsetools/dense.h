#ifndef __DENSE_H__
#define __DENSE_H__

// Simplified BLAS routines and other dense linear algebra functions

/*
 * Level 1
 */

// y += a*x
template <class I, class T>
void axpy(const I n, const T a, const T * x, T * y){
    for(I i = 0; i < n; i++){
        y[i] += a * x[i];
    }
}

// scale a vector in-place
template <class I, class T>
void scal(const I n, const T a, T * x){
    for(I i = 0; i < n; i++){
        x[i] *= a;
    }
}


// dot product
template <class I, class T>
void dot(const I n, const T * x, const T * y){
    T dp = 0;
    for(I i = 0; i < n; i++){
        dp += x[i] * y[i];
    }
    return dp;
}


// vectorize a binary operation
template<class I, class T, class binary_operator>
void vector_binop(const I n, const T * x, const T * y, T * z, 
                  const binary_operator& op)
{
    for(I i = 0; i < n; i++){
        z[i] = op(x[i],y[i]);
    }
}

//template<class I, class T>
//void vector_multiply(const I n, const T * x, const T * y, T * z){
//{
//    vector_binop(n,x,y,z, std::multiplies<T>() );
//}



// Level 2
template <class I, class T>
void gemv(const I m, const I n, const T * A, const T * x, T * y){
    for(I i = 0; i < m; i++){
        T dot = y[i];
        for(I j = 0; j < n; j++){
            dot += A[(npy_intp)n * i + j] * x[j];
        }
        y[i] = dot;
    }
}

// Level 3
template <class I, class T>
void gemm(const I m, const I n, const I k, const T * A, const T * B, T * C){
    for(I i = 0; i < m; i++){
        for(I j = 0; j < n; j++){
            T dot = C[(npy_intp)n * i + j];
            for(I _d = 0; _d < k; _d++){
                dot += A[(npy_intp)k * i + _d] * B[(npy_intp)n * _d + j];
            }
            C[(npy_intp)n * i + j] = dot;
        }
    }
}


#endif
