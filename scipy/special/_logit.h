#ifndef _LOGIT_H_
#define _LOGIT_H_

#include "numpy/npy_math.h"

npy_float logitf(npy_float x);
npy_double logit(npy_double x);
npy_longdouble logitl(npy_longdouble x);

npy_float expitf(npy_float x);
npy_double expit(npy_double x);
npy_longdouble expitl(npy_longdouble x);

#endif
