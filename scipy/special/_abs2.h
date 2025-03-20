#ifndef _ABS2_H_
#define _ABS2_H_

#include <numpy/npy_math.h>

template <typename T>
inline auto _abs2(T x) {
  return x*x;
}

template <typename T>
inline auto _cabs2(T z) {
  return z.real*z.real + z.imag*z.imag;
}

npy_float abs2f(npy_float x, npy_float) {return _abs2(x);}
npy_double abs2(npy_double x, npy_double) {return _abs2(x);}
npy_longdouble abs2g(npy_longdouble x, npy_longdouble) {return _abs2(x);}

npy_float cabs2f(npy_cfloat z, npy_float) {return _cabs2(z);}
npy_double cabs2(npy_cdouble z, npy_double) {return _cabs2(z);}
npy_longdouble cabs2g(npy_clongdouble z, npy_longdouble) {return _cabs2(z);}

#endif
