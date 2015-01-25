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

#endif
