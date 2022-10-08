#ifndef _LOGIT_H_
#define _LOGIT_H_

#include <cmath>

template <typename T>
inline T _logit(T x) {
    x /= 1 - x;
    return std::log(x);
};


template <typename T>
inline T _expit(T x) {
    return 1 / (1 + std::exp(-x));
};


//
// The logistic sigmoid function 'expit' is
//
//     S(x) = 1/(1 + exp(-x))     = exp(x)/(exp(x) + 1)
//
// so
//
// log S(x) = -log(1 + exp(-x))   = x - log(exp(x) + 1)
//          = -log1p(exp(-x))     = x - log1p(exp(x))
//
// By using -log1p(exp(-x)) for x >= 0 and x - log1p(exp(x))
// for x < 0, we extend the range of x values for which we
// obtain accurate results (compared to the naive implementation
// log(expit(x))).
//
template <typename T>
inline T _log_expit(T x) {
    if (x < 0.0) {
        return x - std::log1p(std::exp(x));
    }
    else {
        return -std::log1p(std::exp(-x));
    }
};


npy_float logitf(npy_float x)  {return _logit(x);};
npy_double logit(npy_double x) {return _logit(x);};
npy_longdouble logitl(npy_longdouble x) {return _logit(x);};

npy_float expitf(npy_float x) {return _expit(x);};
npy_double expit(npy_double x) {return _expit(x);};
npy_longdouble expitl(npy_longdouble x) {return _expit(x);};

npy_float log_expitf(npy_float x) {return _log_expit(x);};
npy_double log_expit(npy_double x) {return _log_expit(x);};
npy_longdouble log_expitl(npy_longdouble x) {return _log_expit(x);};

#endif
