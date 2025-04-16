/*
 * Helpers for dealing with npy_cdouble/npy_cfloat types.
 */
#pragma once
#include "numpy/npy_math.h"

namespace _numpymath{

/*
 * std::numeric_limits and useful constants.
 */

template<typename T> struct numeric_limits {};

template<>
struct numeric_limits<float>{
    static constexpr double one = 1.0f;
    static constexpr float nan = std::numeric_limits<float>::quiet_NaN();
};

template<>
struct numeric_limits<double>{
    static constexpr double one = 1.0;
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
};


template<>
struct numeric_limits<npy_cfloat>{
    static constexpr npy_cfloat one = {1.0f, 0.0f};
    static constexpr npy_cfloat nan = {std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN()};
};

template<>
struct numeric_limits<npy_cdouble>{
    static constexpr npy_cdouble one = {1.0, 0.0};
    static constexpr npy_cdouble nan = {std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
};



/* 
 * Grab a real part of a possibly complex array.
 * This is for the work queries.
 * There must be a better way, I'm sure.
 */
inline float real_part(float value){ return value; }
inline double real_part(double value){ return value; }
inline float real_part(npy_cfloat value){ return npy_crealf(value); }
inline double real_part(npy_cdouble value){return npy_creal(value); }



/* 
 * Helpers for filling/rearranging matrices
 */


/* identity square matrix generation */
template<typename typ>
static inline void
identity_matrix(typ *matrix, ptrdiff_t n)
{
    ptrdiff_t i;
    /* in IEEE floating point, zeroes are represented as bitwise 0 */
    memset((void *)matrix, 0, n*n*sizeof(typ));

    for (i = 0; i < n; ++i)
    {
        *matrix = numeric_limits<typ>::one;
        matrix += n+1;
    }
}


/* m-by-n matrix full of nans */
template<typename typ>
static inline void
nan_matrix(typ *matrix, ptrdiff_t m, ptrdiff_t n)
{
    for (ptrdiff_t i=0; i < m*n; i++) {
        *(matrix + i) = numeric_limits<typ>::nan;
    }
}


}  // namespace _numpymath
