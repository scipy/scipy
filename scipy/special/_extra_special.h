#pragma once

#include <complex>
#include <numpy/npy_math.h>

#include "_lambertw.h"


inline npy_cdouble lambertw_scalar(npy_cdouble zp, long k, double tol) {
    std::complex<double>* z_ptr = reinterpret_cast<std::complex<double>*>(&zp);
    std::complex<double> w = lambertw(*z_ptr, k, tol);
    return *(reinterpret_cast<npy_cdouble*>(&w));
}
