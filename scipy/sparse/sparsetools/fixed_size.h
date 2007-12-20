#ifndef FIXED_SIZE_H
#define FIXED_SIZE_H

/*
 * templates for fixed size array and vector arithmetic
 * 
 */

template<int N, class T>
class Dot
{
    public:
        inline T operator()(const T * lhs, const T * rhs)
        {
            Dot<N-1,T> d;
            return (*lhs * *rhs) + d(++lhs, ++rhs);
        }
};

template<class T>
class Dot<0,T>
{
    public:
        inline T operator()(const T * lhs, const T * rhs)
        {
            return 0;
        }
};

    template<int N, class T>
inline T dot(const T * lhs, const T * rhs)
{
    Dot<N,T> d;
    return d(lhs, rhs);
}



#endif
