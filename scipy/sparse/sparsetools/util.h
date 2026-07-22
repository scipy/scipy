#ifndef __SPTOOLS_UTIL_H__
#define __SPTOOLS_UTIL_H__

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

#define SPTOOLS_FOR_EACH_DATA_TYPE(X)           \
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


#define SPTOOLS_FOR_EACH_INDEX_DATA_TYPE_COMBINATION(X) \
  X(npy_int32, npy_bool_wrapper)                        \
  X(npy_int32, npy_byte)                                \
  X(npy_int32, npy_ubyte)                               \
  X(npy_int32, npy_short)                               \
  X(npy_int32, npy_ushort)                              \
  X(npy_int32, npy_int)                                 \
  X(npy_int32, npy_uint)                                \
  X(npy_int32, npy_long)                                \
  X(npy_int32, npy_ulong)                               \
  X(npy_int32, npy_longlong)                            \
  X(npy_int32, npy_ulonglong)                           \
  X(npy_int32, npy_float)                               \
  X(npy_int32, npy_double)                              \
  X(npy_int32, npy_longdouble)                          \
  X(npy_int32, npy_cfloat_wrapper)                      \
  X(npy_int32, npy_cdouble_wrapper)                     \
  X(npy_int32, npy_clongdouble_wrapper)                 \
  X(npy_int64, npy_bool_wrapper)                        \
  X(npy_int64, npy_byte)                                \
  X(npy_int64, npy_ubyte)                               \
  X(npy_int64, npy_short)                               \
  X(npy_int64, npy_ushort)                              \
  X(npy_int64, npy_int)                                 \
  X(npy_int64, npy_uint)                                \
  X(npy_int64, npy_long)                                \
  X(npy_int64, npy_ulong)                               \
  X(npy_int64, npy_longlong)                            \
  X(npy_int64, npy_ulonglong)                           \
  X(npy_int64, npy_float)                               \
  X(npy_int64, npy_double)                              \
  X(npy_int64, npy_longdouble)                          \
  X(npy_int64, npy_cfloat_wrapper)                      \
  X(npy_int64, npy_cdouble_wrapper)                     \
  X(npy_int64, npy_clongdouble_wrapper)

#else  // vvv npy_longdouble is npy_double vvv

#define SPTOOLS_FOR_EACH_DATA_TYPE(X)           \
  X(npy_bool_wrapper)                           \
  X(npy_byte)                                   \
  X(npy_ubyte)                                  \
  X(npy_short)                                  \
  X(npy_ushort)                                 \
  X(npy_int)                                    \
  X(npy_uint)                                   \
  X(npy_long)                                   \
  X(npy_ulong)                                  \
  X(npy_longlong)                               \
  X(npy_ulonglong)                              \
  X(npy_float)                                  \
  X(npy_double)                                 \
  X(npy_cfloat_wrapper)                         \
  X(npy_cdouble_wrapper)                        \
  X(npy_clongdouble_wrapper)


#define SPTOOLS_FOR_EACH_INDEX_DATA_TYPE_COMBINATION(X) \
  X(npy_int32, npy_bool_wrapper)                        \
  X(npy_int32, npy_byte)                                \
  X(npy_int32, npy_ubyte)                               \
  X(npy_int32, npy_short)                               \
  X(npy_int32, npy_ushort)                              \
  X(npy_int32, npy_int)                                 \
  X(npy_int32, npy_uint)                                \
  X(npy_int32, npy_long)                                \
  X(npy_int32, npy_ulong)                               \
  X(npy_int32, npy_longlong)                            \
  X(npy_int32, npy_ulonglong)                           \
  X(npy_int32, npy_float)                               \
  X(npy_int32, npy_double)                              \
  X(npy_int32, npy_cfloat_wrapper)                      \
  X(npy_int32, npy_cdouble_wrapper)                     \
  X(npy_int32, npy_clongdouble_wrapper)                 \
  X(npy_int64, npy_bool_wrapper)                        \
  X(npy_int64, npy_byte)                                \
  X(npy_int64, npy_ubyte)                               \
  X(npy_int64, npy_short)                               \
  X(npy_int64, npy_ushort)                              \
  X(npy_int64, npy_int)                                 \
  X(npy_int64, npy_uint)                                \
  X(npy_int64, npy_long)                                \
  X(npy_int64, npy_ulong)                               \
  X(npy_int64, npy_longlong)                            \
  X(npy_int64, npy_ulonglong)                           \
  X(npy_int64, npy_float)                               \
  X(npy_int64, npy_double)                              \
  X(npy_int64, npy_cfloat_wrapper)                      \
  X(npy_int64, npy_cdouble_wrapper)                     \
  X(npy_int64, npy_clongdouble_wrapper)

#endif

#endif
