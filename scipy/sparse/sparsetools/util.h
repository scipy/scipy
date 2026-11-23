#ifndef __SPTOOLS_UTIL_H__
#define __SPTOOLS_UTIL_H__

#include <algorithm>
#include <numpy/ndarrayobject.h>

#include "bool_ops.h"
#include "complex_ops.h"

/*
 * Same as std::divides, except return x/0 == 0 for integer types, without
 * raising a SIGFPE.
 */
template <class T>
struct safe_divides {
    T operator() (const T& x, const T& y) const {
        if (y == 0) {
            return 0;
        }
        else {
            return x/y;
        }
    }

    typedef T first_argument_type;
    typedef T second_argument_type;
    typedef T result_type;
};

#define OVERRIDE_safe_divides(typ) \
    template<> inline typ safe_divides<typ>::operator()(const typ& x, const typ& y) const { return x/y; }

OVERRIDE_safe_divides(float)
OVERRIDE_safe_divides(double)
OVERRIDE_safe_divides(long double)
OVERRIDE_safe_divides(npy_cfloat_wrapper)
OVERRIDE_safe_divides(npy_cdouble_wrapper)
OVERRIDE_safe_divides(npy_clongdouble_wrapper)

#undef OVERRIDE_safe_divides

template <class T>
struct maximum {
    T operator() (const T& x, const T& y) const {
        return std::max(x, y);
    }
};

template <class T>
struct minimum {
    T operator() (const T& x, const T& y) const {
        return std::min(x, y);
    }
};

#define SPTOOLS_FOR_EACH_INDEX_TYPE(X) \
  X(npy_int32)                         \
  X(npy_int64)

#define SPTOOLS_FOR_EACH_DATA_TYPE_CODE(X)      \
  X(NPY_BOOL, npy_bool_wrapper)                 \
  X(NPY_BYTE, npy_byte)                         \
  X(NPY_UBYTE, npy_ubyte)                       \
  X(NPY_SHORT, npy_short)                       \
  X(NPY_USHORT, npy_ushort)                     \
  X(NPY_INT, npy_int)                           \
  X(NPY_UINT, npy_uint)                         \
  X(NPY_LONG, npy_long)                         \
  X(NPY_ULONG, npy_ulong)                       \
  X(NPY_LONGLONG, npy_longlong)                 \
  X(NPY_ULONGLONG, npy_ulonglong)               \
  X(NPY_FLOAT, npy_float)                       \
  X(NPY_DOUBLE, npy_double)                     \
  X(NPY_LONGDOUBLE, npy_longdouble)             \
  X(NPY_CFLOAT, npy_cfloat_wrapper)             \
  X(NPY_CDOUBLE, npy_cdouble_wrapper)           \
  X(NPY_CLONGDOUBLE, npy_clongdouble_wrapper)

// If npy_longdouble is distinct from npy_double
#if NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE

#define SPTOOLS_IF_DISTINCT_LONGDOUBLE(x) x

#else  // vvv npy_longdouble is npy_double vvv

#define SPTOOLS_IF_DISTINCT_LONGDOUBLE(x)

#endif

#define SPTOOLS_INVOKE_1(X, T) X(T)
#define SPTOOLS_INVOKE_INT32(X, T) X(npy_int32, T)
#define SPTOOLS_INVOKE_INT64(X, T) X(npy_int64, T)

#define SPTOOLS_FOR_EACH_DATA_TYPE_IMPL(F, X)          \
  F(X, npy_bool_wrapper)                               \
  F(X, npy_byte)                                       \
  F(X, npy_ubyte)                                      \
  F(X, npy_short)                                      \
  F(X, npy_ushort)                                     \
  F(X, npy_int)                                        \
  F(X, npy_uint)                                       \
  F(X, npy_long)                                       \
  F(X, npy_ulong)                                      \
  F(X, npy_longlong)                                   \
  F(X, npy_ulonglong)                                  \
  F(X, npy_float)                                      \
  F(X, npy_double)                                     \
  SPTOOLS_IF_DISTINCT_LONGDOUBLE(F(X, npy_longdouble)) \
  F(X, npy_cfloat_wrapper)                             \
  F(X, npy_cdouble_wrapper)                            \
  F(X, npy_clongdouble_wrapper)

#define SPTOOLS_FOR_EACH_DATA_TYPE(X) \
  SPTOOLS_FOR_EACH_DATA_TYPE_IMPL(SPTOOLS_INVOKE_1, X)

#define SPTOOLS_FOR_EACH_INDEX_DATA_TYPE_COMBINATION(X)    \
  SPTOOLS_FOR_EACH_DATA_TYPE_IMPL(SPTOOLS_INVOKE_INT32, X) \
  SPTOOLS_FOR_EACH_DATA_TYPE_IMPL(SPTOOLS_INVOKE_INT64, X)

#endif
