/*
 * Helpers for dealing with npy_complex128/npy_complex64 types.
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
struct numeric_limits<npy_complex64>{
    static constexpr npy_complex64 zero = {0.0f, 0.0f};
    static constexpr npy_complex64 one = {1.0f, 0.0f};
    static constexpr npy_complex64 nan = {std::numeric_limits<float>::quiet_NaN(),
                                       std::numeric_limits<float>::quiet_NaN()};
};

template<>
struct numeric_limits<npy_complex128>{
    static constexpr npy_complex128 zero = {0.0, 0.0};
    static constexpr npy_complex128 one = {1.0, 0.0};
    static constexpr npy_complex128 nan = {std::numeric_limits<double>::quiet_NaN(),
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
template<> struct type_traits<npy_complex64> {
    using real_type = float;
    using value_type = std::complex<float>;
    static constexpr int typenum = NPY_COMPLEX64;
    static constexpr bool is_complex = true;
};
template<> struct type_traits<npy_complex128> {
    using real_type = double;
    using value_type = std::complex<double>;
    static constexpr int typenum = NPY_COMPLEX128;
    static constexpr bool is_complex = true;
};


/*
 * Grab a real part of a possibly complex array.
 * This is for the work queries.
 */
inline float real_part(float value){ return value; }
inline double real_part(double value){ return value; }
inline float real_part(npy_complex64 value){ return npy_crealf(value); }
inline double real_part(npy_complex128 value){return npy_creal(value); }

// XXX: std::conj(double) -> complex(double)
inline float conj(float value){ return value; }
inline double conj(double value){ return value; }
inline npy_complex64 conj(npy_complex64 value){ return npy_cpackf(npy_crealf(value), -npy_cimagf(value)); }
inline npy_complex128 conj(npy_complex128 value){return npy_cpack(npy_creal(value), -npy_cimag(value)); }


/*
 * Debug helper: print out an npy_{cfloat,cdouble} value
 */
std::ostream& operator<<(std::ostream& os, npy_complex64 x) {
    os << "(" << npy_crealf(x) << ", " << npy_cimagf(x) << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, npy_complex128 x) {
    os << "(" << npy_creal(x) << ", " << npy_cimag(x) << ")";
    return os;
}
