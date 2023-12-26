#ifndef __SPTOOLS_UTIL_H__
#define __SPTOOLS_UTIL_H__

#include <numeric>

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

/**
 * Zips a key and value iterator.
 * Enables sorting both index and value arrays in-place with std::sort.
 */
template<class RandItA, class RandItB>
class kv_pair_iters {
public:
    using value_type = std::pair<typename std::iterator_traits<RandItA>::value_type,
                                 typename std::iterator_traits<RandItB>::value_type>;
    using difference_type = typename std::iterator_traits<RandItA>::difference_type;
    using pointer = value_type *;
    using iterator_category = std::random_access_iterator_tag;

    struct reference {
        reference(RandItA a, RandItB b) : a(a), b(b) {}

        reference &operator=(reference &&rhs) { *a = std::move(*rhs.a); *b = std::move(*rhs.b); return *this; }
        reference &operator=(const value_type &&rhs) { *a = std::move(rhs.first); *b = std::move(rhs.second); return *this; }
//        operator value_type() && { return std::make_pair(std::move(*a), std::move(*b)); }
        operator value_type() { return std::make_pair(*a, *b); }

        RandItA a;
        RandItB b;

        friend void swap(reference lhs, reference rhs) {
            std::iter_swap(lhs.a, rhs.a);
            std::iter_swap(lhs.b, rhs.b);
        }
    };

    kv_pair_iters() : a{}, b{} {}
    kv_pair_iters(RandItA a, RandItB b) : a(a), b(b) {}
    kv_pair_iters(const kv_pair_iters<RandItA, RandItB> &rhs) : a(rhs.a), b(rhs.b) {}

    kv_pair_iters<RandItA, RandItB> &operator=(const kv_pair_iters &rhs) { a = rhs.a; b = rhs.b; return *this; }

    reference operator*() const { return reference(a, b); }
    reference operator[](difference_type rhs) const { return ref_proxy(a + rhs, b + rhs); }
    // operator-> has no meaning in this application

    bool operator==(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a == rhs.a; }
    bool operator!=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a != rhs.a; }
    bool operator<(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a < rhs.a; }
    bool operator>(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a > rhs.a; }
    bool operator<=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a <= rhs.a; }
    bool operator>=(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a >= rhs.a; }

    kv_pair_iters<RandItA, RandItB> &operator+=(difference_type rhs) { a += rhs; b += rhs; return *this; }
    kv_pair_iters<RandItA, RandItB> &operator-=(difference_type rhs) { a -= rhs; b -= rhs; return *this; }

    kv_pair_iters<RandItA, RandItB> &operator++() { ++a; ++b; return *this; }
    kv_pair_iters<RandItA, RandItB> &operator--() { --a; --b; return *this; }
    kv_pair_iters<RandItA, RandItB> operator++(int) { kv_pair_iters<RandItA, RandItB> ret(*this); ++a; ++b; return ret; }
    kv_pair_iters<RandItA, RandItB> operator--(int) { kv_pair_iters<RandItA, RandItB> ret(*this); --a; --b; return ret; }

    difference_type operator-(const kv_pair_iters<RandItA, RandItB> &rhs) const { return a - rhs.a; }
    kv_pair_iters<RandItA, RandItB> operator-(difference_type rhs) const { return kv_pair_iters(a - rhs, b - rhs); }
    kv_pair_iters<RandItA, RandItB> operator+(difference_type rhs) const { return kv_pair_iters(a + rhs, b + rhs); }

    friend inline kv_pair_iters<RandItA, RandItB> operator+(difference_type lhs, const kv_pair_iters<RandItA, RandItB> &rhs) {
        return kv_pair_iters(lhs + rhs.a, lhs + rhs.b);
    }

protected:
    RandItA a;
    RandItB b;
};

namespace std {
    /**
     * Specialize std::iterator_traits for kv_pair_iters.
     */
    template<class A, class B>
    struct iterator_traits<kv_pair_iters<A, B>> {
        using value_type =        typename kv_pair_iters<A, B>::value_type;
        using difference_type =   typename kv_pair_iters<A, B>::difference_type;
        using pointer =           typename kv_pair_iters<A, B>::pointer;
        using reference =         typename kv_pair_iters<A, B>::reference;
        using iterator_category = typename kv_pair_iters<A, B>::iterator_category;
    };
}

// GCC 8 does not have std::exclusive_scan. Reimplement if necessary.
template<typename InputIt, typename OutputIt, typename T>
OutputIt exclusive_scan(InputIt first, InputIt last, OutputIt result, T init)
{
#if (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 9)
    return std::exclusive_scan(first, last, result, init);
#else
    for (; first != last; ++first) {
        auto v = init;
        init = init + *first;
        *result++ = v;
    }
    return result;
#endif
}

// GCC 8 does not have std::inclusive_scan. Reimplement if necessary.
template<typename InputIt, typename OutputIt>
OutputIt inclusive_scan(InputIt first, InputIt last, OutputIt result)
{
#if (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 9)
    return std::inclusive_scan(first, last, result);
#else
    typename std::iterator_traits<OutputIt>::value_type init = 0;

    for (; first != last; ++first) {
        init = init + *first;
        *result++ = init;
    }
    return result;
#endif
}

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
