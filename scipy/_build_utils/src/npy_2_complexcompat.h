#ifndef NUMPY_CORE_INCLUDE_NUMPY_NPY_2_COMPLEXCOMPAT_H_
#define NUMPY_CORE_INCLUDE_NUMPY_NPY_2_COMPLEXCOMPAT_H_

#include <numpy/npy_math.h>

#ifndef NPY_CSETREALF
#define NPY_CSETREALF(c, r) (c)->real = (r)
#endif
#ifndef NPY_CSETIMAGF
#define NPY_CSETIMAGF(c, i) (c)->imag = (i)
#endif
#ifndef NPY_CSETREAL
#define NPY_CSETREAL(c, r)  (c)->real = (r)
#endif
#ifndef NPY_CSETIMAG
#define NPY_CSETIMAG(c, i)  (c)->imag = (i)
#endif
#ifndef NPY_CSETREALL
#define NPY_CSETREALL(c, r) (c)->real = (r)
#endif
#ifndef NPY_CSETIMAGL
#define NPY_CSETIMAGL(c, i) (c)->imag = (i)
#endif

#endif
