#ifndef UFUNCS_PROTO_H
#define UFUNCS_PROTO_H 1
#include "_faddeeva.h"
npy_double faddeeva_dawsn(npy_double);
npy_cdouble faddeeva_dawsn_complex(npy_cdouble);
npy_cdouble faddeeva_erf(npy_cdouble);
npy_cdouble faddeeva_erfc(npy_cdouble);
npy_double faddeeva_erfcx(npy_double);
npy_cdouble faddeeva_erfcx_complex(npy_cdouble);
npy_double faddeeva_erfi(npy_double);
npy_cdouble faddeeva_erfi_complex(npy_cdouble);
npy_cdouble faddeeva_log_ndtr(npy_cdouble);
npy_cdouble faddeeva_ndtr(npy_cdouble);
npy_cdouble faddeeva_w(npy_cdouble);
#include "_wright.h"
npy_cdouble wrightomega(npy_cdouble);
#endif
