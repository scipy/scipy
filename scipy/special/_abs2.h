#ifndef _ABS2_H_
#define _ABS2_H_

#include <numpy/npy_math.h>

double abs2(double x) {
  return x*x;
}

double cabs2(npy_cdouble z) {
  return z.real*z.real + z.imag*z.imag;
}

#endif