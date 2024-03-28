#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/npy_math.h>

    npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);

#ifdef __cplusplus
}
#endif
