#ifndef __SPTOOLS_UTIL_H__
#define __SPTOOLS_UTIL_H__

#include <stdexcept>

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

/*
 * Operator that implements multiple binops to be chosen at runtime.
 * Allows building one copy of the large binop methods to be reused for several
 * operands. Significantly speeds up compilation time and reduces library size.
 * Binops are memory-bound operations, so the performance impact is negligible.
 */
template <class T>
struct multi_op {
    enum op_choice {MULTIPLY, DIVIDE, MAXIMUM, MINIMUM};
    explicit multi_op(op_choice which): which(which) {}

    T operator()(const T& a, const T& b) const {
        switch (which) {
        case MULTIPLY: return (a * b);
        case DIVIDE: return safe_divides<T>()(a, b);
        case MAXIMUM: return std::max<T>(a, b);
        case MINIMUM: return std::min<T>(a, b);
        default: throw std::invalid_argument("Invalid multi_op argument");
        }
    }

protected:
    const op_choice which;
};

/*
 * Comparator version of multi_op.
 */
template <class T>
struct multi_op_bool {
    enum op_choice {NOT_EQUAL_TO, LESS, GREATER, LESS_EQ, GREATER_EQ};
    explicit multi_op_bool(op_choice which): which(which) {}

    bool operator()(const T& a, const T& b) const {
        switch (which) {
        case NOT_EQUAL_TO: return (a != b);
        case GREATER: return (a > b);
        case LESS: return (a < b);
        case GREATER_EQ: return (a >= b);
        case LESS_EQ: return (a <= b);
        default: throw std::invalid_argument("Invalid multi_op_bool argument");
        }
    }

protected:
    const op_choice which;
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
