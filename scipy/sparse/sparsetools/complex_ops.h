#ifndef COMPLEX_OPS_H
#define COMPLEX_OPS_H

/*
 *  Functions to handle arithmetic operations on NumPy complex values
 */

#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include "npy_2_complexcompat.h"

template <class c_type, class npy_type>
class complex_wrapper {
    private:
        npy_type complex;
        c_type real() const { return c_type(0); }
        c_type imag() const { return c_type(0); }
        void set_real(const c_type r) { }
        void set_imag(const c_type i) { }

    public:
        /* Constructor */
        complex_wrapper( const c_type r = c_type(0), const c_type i = c_type(0) ){
            set_real(r);
            set_imag(i);
        }
        /* Conversion */
        operator bool() const {
            if (real() == 0 && imag() == 0) {
                return false;
            } else {
                return true;
            }
        }
        /* Operators */
        complex_wrapper operator-() const {
          return complex_wrapper(-real(),-imag());
        }
        complex_wrapper operator+(const complex_wrapper& B) const {
          return complex_wrapper(real() + B.real(), imag() + B.imag());
        }
        complex_wrapper operator-(const complex_wrapper& B) const {
          return complex_wrapper(real() - B.real(), imag() - B.imag());
        }
        complex_wrapper operator*(const complex_wrapper& B) const {
          return complex_wrapper(real() * B.real() - imag() * B.imag(),
                                 real() * B.imag() + imag() * B.real());
        }
        complex_wrapper operator/(const complex_wrapper& B) const {
            complex_wrapper result;
            c_type denom = 1.0 / (B.real() * B.real() + B.imag() * B.imag());
            result.set_real((real() * B.real() + imag() * B.imag()) * denom);
            result.set_imag((imag() * B.real() - real() * B.imag()) * denom);
            return result;
        }
        /* in-place operators */
        complex_wrapper& operator+=(const complex_wrapper & B){
          set_real(real() + B.real());
          set_imag(imag() + B.imag());
          return (*this);
        }
        complex_wrapper& operator-=(const complex_wrapper & B){
          
          set_real(real() - B.real());
          set_imag(imag() - B.imag());
          return (*this);
        }
        complex_wrapper& operator*=(const complex_wrapper & B){
          c_type temp    = real() * B.real() - imag() * B.imag();
          set_imag(real() * B.imag() + imag() * B.real());
          set_real(temp);
          return (*this);
        }
        complex_wrapper& operator/=(const complex_wrapper & B){
          c_type denom   = 1.0 / (B.real() * B.real() + B.imag() * B.imag());
          c_type temp    = (real() * B.real() + imag() * B.imag()) * denom; 
          set_imag((imag() * B.real() - real() * B.imag()) * denom);
          set_real(temp);
          return (*this);
        }
        /* Boolean operations */
        bool operator==(const complex_wrapper& B) const{
          return real() == B.real() && imag() == B.imag();
        }
        bool operator!=(const complex_wrapper& B) const{
          return real() != B.real() || imag() != B.imag();
        }
        bool operator<(const complex_wrapper& B) const{
            if (real() == B.real()){
                return imag() < B.imag();
            } else {
                return real() < B.real();
            }
        }
        bool operator>(const complex_wrapper& B) const{
            if (real() == B.real()){
                return imag() > B.imag();
            } else {
                return real() > B.real();
            }
        }
        bool operator<=(const complex_wrapper& B) const{
            if (real() == B.real()){
                return imag() <= B.imag();
            } else {
                return real() <= B.real();
            }
        }
        bool operator>=(const complex_wrapper& B) const{
            if (real() == B.real()){
                return imag() >= B.imag();
            } else {
                return real() >= B.real();
            }
        }
        template <class T>
        bool operator==(const T& B) const{
          return real() == B && imag() == T(0);
        }
        template <class T>
        bool operator!=(const T& B) const{
          return real() != B || imag() != T(0);
        }
        template <class T>
        bool operator<(const T& B) const{
            if (real() == B) {
                return imag() < T(0);
            } else {
                return real() < B;
            }
        }
        template <class T>
        bool operator>(const T& B) const{
            if (real() == B) {
                return imag() > T(0);
            } else {
                return real() > B;
            }
        }
        template <class T>
        bool operator<=(const T& B) const{
            if (real() == B) {
                return imag() <= T(0);
            } else {
                return real() <= B;
            }
        }
        template <class T>
        bool operator>=(const T& B) const{
            if (real() == B) {
                return imag() >= T(0);
            } else {
                return real() >= B;
            }
        }
        complex_wrapper& operator=(const complex_wrapper& B){
          set_real(B.real());
          set_imag(B.imag());
          return (*this);
        }
        complex_wrapper& operator=(const c_type& B){
          set_real(B);
          set_imag(c_type(0));
          return (*this);
        }
};

template <>
inline float complex_wrapper<float, npy_cfloat>::real() const {
    return npy_crealf(this->complex);
}

template <>
inline void complex_wrapper<float, npy_cfloat>::set_real(const float r) {
    NPY_CSETREALF(&this->complex, r);
}

template <>
inline double complex_wrapper<double, npy_cdouble>::real() const {
    return npy_creal(this->complex);
}

template <>
inline void complex_wrapper<double, npy_cdouble>::set_real(const double r) {
    NPY_CSETREAL(&this->complex, r);
}

template <>
inline long double complex_wrapper<long double, npy_clongdouble>::real() const {
    return npy_creall(this->complex);
}

template <>
inline void complex_wrapper<long double, npy_clongdouble>::set_real(const long double r) {
    NPY_CSETREALL(&this->complex, r);
}

template <>
inline float complex_wrapper<float, npy_cfloat>::imag() const {
    return npy_cimagf(this->complex);
}

template <>
inline void complex_wrapper<float, npy_cfloat>::set_imag(const float i) {
    NPY_CSETIMAGF(&this->complex, i);
}

template <>
inline double complex_wrapper<double, npy_cdouble>::imag() const {
    return npy_cimag(this->complex);
}

template <>
inline void complex_wrapper<double, npy_cdouble>::set_imag(const double i) {
    NPY_CSETIMAG(&this->complex, i);
}

template <>
inline long double complex_wrapper<long double, npy_clongdouble>::imag() const {
    return npy_cimagl(this->complex);
}

template <>
inline void complex_wrapper<long double, npy_clongdouble>::set_imag(const long double i) {
    NPY_CSETIMAGL(&this->complex, i);
}

typedef complex_wrapper<float,npy_cfloat> npy_cfloat_wrapper;
typedef complex_wrapper<double,npy_cdouble> npy_cdouble_wrapper;
typedef complex_wrapper<long double,npy_clongdouble> npy_clongdouble_wrapper;

#endif
