#pragma once

#include <numpy/npy_math.h>

#ifdef __cplusplus
extern "C" {
#endif

    npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);

#ifdef __cplusplus
}
#endif
