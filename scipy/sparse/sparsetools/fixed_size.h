#ifndef FIXED_SIZE_H
#define FIXED_SIZE_H

/*
 * templates for fixed size vector and matrix arithmetic
 * 
 */



/*
 *  Dot Product
 * 
 */
template<int N, int SX, int SY, class T>
class _dot
{
    public:
        inline T operator()(const T * X, const T * Y)
        {
            _dot<N-1,SX,SY,T> d;
            return (*X) * (*Y) + d(X + SX, Y + SY);
        }
};
template<int SX, int SY, class T>
class _dot<1,SX,SY,T>
{
    public:
        inline T operator()(const T * X, const T * Y)
        {
            return (*X) * (*Y);
        }
};

template<int N, int SX, int SY, class T>
inline T dot(const T * X, const T * Y)
{
    _dot<N,SX,SY,T> d;
    return d(X, Y);
}



/*
 *  Matrix Vector Product Y = A*X
 * 
 */
template<int M, int N, int SX, int SY, class T>
class _matvec
{
    public:
        inline void operator()(const T * A, const T * X, T * Y)
        {
            *Y += dot<N,1,SX,T>(A,X);
            _matvec<M-1,N,SX,SY,T> d;
            d(A + N, X, Y + SY);
        }
};

template<int N, int SX, int SY, class T>
class _matvec<1,N,SX,SY,T>
{
    public:
        inline void operator()(const T * A, const T * X, T * Y)
        {
            *Y += dot<N,1,SX,T>(A,X);
        }
};

template<int M, int N, int SX, int SY, class T>
inline void matvec(const T * A, const T * X, T * Y)
{
    _matvec<M,N,SX,SY,T> d;
    d(A,X,Y);
}


/*
 *  Matrix Matrix Product C = A*B
 *
 *  C is L*N
 *  A is L*M
 *  B is M*N
 * 
 */
template<int L, int M, int N, int U, class T>
class _matmat
{
    public:
        inline void operator()(const T * A, const T * B, T * C)
        {
            matvec<L,M,N,N>(A,B,C);
            
            _matmat<L,M,N,U-1,T> d;
            d(A, B + 1, C + 1);
        }
};
template<int L, int M, int N, class T>
class _matmat<L,M,N,0,T>
{
    public:
        inline void operator()(const T * A, const T * B, T * C)
        {
            matvec<L,M,N,N>(A,B,C);
        }
};

template<int L, int M, int N, class T>
inline void matmat(const T * A, const T * B, T * C)
{
    _matmat<L,M,N,N-1,T> d;
    d(A,B,C);
}



/*
 * Binary vector operation Z = op(X,Y) 
 *
 */

template<int N, class T, class bin_op>
class _vec_binop_vec
{
    public:
        inline void operator()(const T * X, const T * Y, T * Z, const bin_op& op)
        {
            *Z = op( *X, *Y );
            _vec_binop_vec<N-1,T,bin_op> d;
            d(X + 1, Y + 1, Z + 1, op);
        }
};
template<class T, class bin_op>
class _vec_binop_vec<1,T,bin_op>
{
    public:
        inline void operator()(const T * X, const T * Y, T * Z, const bin_op& op)
        {
            *Z = op( *X, *Y );
        }
};

template<int N, class T, class bin_op>
inline void vec_binop_vec(const T * X, const T * Y, T * Z, const bin_op& op)
{
    _vec_binop_vec<N,T,bin_op> d;
    d(X,Y,Z,op);
}




#endif
