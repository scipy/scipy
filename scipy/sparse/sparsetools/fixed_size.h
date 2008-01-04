#ifndef FIXED_SIZE_H
#define FIXED_SIZE_H

/*
 * templates for fixed size array and vector arithmetic
 * 
 */



/*
 *  Dot Product
 * 
 */
template<int N, class T>
class _dot
{
    public:
        inline T operator()(const T * V1, const T * V2)
        {
            _dot<N-1,T> d;
            return (*V1) * (*V2) + d(V1 + 1, V2 + 1);
        }
};
template<class T>
class _dot<1,T>
{
    public:
        inline T operator()(const T * V1, const T * V2)
        {
            return (*V1) * (*V2);
        }
};

template<int N, class T>
inline T dot(const T * V1, const T * V2)
{
    _dot<N,T> d;
    return d(V1, V2);
}



/*
 *  Matrix Vector Product
 * 
 */
template<int M, int N, class T>
class _matvec
{
    public:
        inline void operator()(const T * A, const T * X, T * Y)
        {
            *Y += dot<N,T>(A,X);
            _matvec<M-1,N,T> d;
            d(A + N, X, Y + 1);
        }
};
template<int N, class T>
class _matvec<1,N,T>
{
    public:
        inline void operator()(const T * A, const T * X, T * Y)
        {
            *Y += dot<N,T>(A,X);
        }
};

template<int M, int N, class T>
inline void matvec(const T * A, const T * X, T * Y)
{
    _matvec<M,N,T> d;
    d(A,X,Y);
}



template<int N, class T, class bin_op>
class _vec_binop_vec
{
    public:
        inline void operator()(const T * V1, const T * V2, T * V3, const bin_op& op)
        {
            *V3 = op( *V1, *V2 );
            _vec_binop_vec<N-1,T,bin_op> d;
            d(V1 + 1, V2 + 1, V3 + 1, op);
        }
};
template<class T, class bin_op>
class _vec_binop_vec<1,T,bin_op>
{
    public:
        inline void operator()(const T * V1, const T * V2, T * V3, const bin_op& op)
        {
            *V3 = op( *V1, *V2 );
        }
};

template<int N, class T, class bin_op>
inline void vec_binop_vec(const T * V1, const T * V2, T * V3, const bin_op& op)
{
    _vec_binop_vec<N,T,bin_op> d;
    d(V1,V2,V3,op);
}




#endif
