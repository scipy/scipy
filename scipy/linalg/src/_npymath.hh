/*
 * Helpers for dealing with npy_cdouble/npy_cfloat types.
 */
#pragma once
#include "Python.h"
#include <iostream>
#include <complex>
#include "numpy/npy_math.h"


/*
 * std::numeric_limits and useful constants.
 */

template<typename T> struct numeric_limits {};

template<>
struct numeric_limits<float>{
    static constexpr double zero = 0.0f;
    static constexpr double one = 1.0f;
    static constexpr float nan = std::numeric_limits<float>::quiet_NaN();
    static constexpr float eps = std::numeric_limits<float>::epsilon();
};

template<>
struct numeric_limits<double>{
    static constexpr double zero = 0.0;
    static constexpr double one = 1.0;
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    static constexpr double eps = std::numeric_limits<double>::epsilon();
};


template<>
struct numeric_limits<npy_cfloat>{
    static constexpr npy_cfloat zero = {0.0f, 0.0f};
    static constexpr npy_cfloat one = {1.0f, 0.0f};
    static constexpr npy_cfloat nan = {std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN()};
};

template<>
struct numeric_limits<npy_cdouble>{
    static constexpr npy_cdouble zero = {0.0, 0.0};
    static constexpr npy_cdouble one = {1.0, 0.0};
    static constexpr npy_cdouble nan = {std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
};

/*
 * XXX merge with numeric_limits ?
 */
template<typename T> struct type_traits {};
template<> struct type_traits<float> {
    using real_type = float;
    using value_type = float;
    static constexpr int typenum = NPY_FLOAT;
    static constexpr bool is_complex = false;

};
template<> struct type_traits<double> {
    using real_type = double;
    using value_type = double;
    static constexpr int typenum = NPY_DOUBLE;
    static constexpr bool is_complex = false;
};
template<> struct type_traits<npy_cfloat> {
    using real_type = float;
    using value_type = std::complex<float>; 
    static constexpr int typenum = NPY_CFLOAT;
    static constexpr bool is_complex = true;
};
template<> struct type_traits<npy_cdouble> {
    using real_type = double;
    using value_type = std::complex<double>;
    static constexpr int typenum = NPY_CDOUBLE;
    static constexpr bool is_complex = true;
};


/* 
 * Grab a real part of a possibly complex array.
 * This is for the work queries.
 */
inline float real_part(float value){ return value; }
inline double real_part(double value){ return value; }
inline float real_part(npy_cfloat value){ return npy_crealf(value); }
inline double real_part(npy_cdouble value){return npy_creal(value); }

// XXX: std::conj(double) -> complex(double)
inline float conj(float value){ return value; }
inline double conj(double value){ return value; }
inline npy_cfloat conj(npy_cfloat value){ return npy_cpackf(npy_crealf(value), -npy_cimagf(value)); }
inline npy_cdouble conj(npy_cdouble value){return npy_cpack(npy_creal(value), -npy_cimag(value)); }


/*
 *  abs : XXX no longer needed
 */
inline float abs_(float x) {return fabsf(x);}
inline double abs_(double x) {return fabs(x);}
inline float abs_(npy_cfloat x) {return sqrtf(npy_crealf(x)*npy_crealf(x) + npy_cimagf(x)*npy_cimagf(x));}
inline double abs_(npy_cdouble x) {return sqrt(npy_creal(x)*npy_creal(x) + npy_cimag(x)*npy_cimag(x));}



/*
 * Debug helper: print out an npy_{cfloat,cdouble} value
 */
std::ostream& operator<<(std::ostream& os, npy_cfloat x) {
    os << "(" << npy_crealf(x) << ", " << npy_cimagf(x) << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, npy_cdouble x) {
    os << "(" << npy_creal(x) << ", " << npy_cimag(x) << ")";
    return os;
}
