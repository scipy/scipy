#ifndef _LOGIT_H_
#define _LOGIT_H_

#include<cmath>

template <typename T>
inline T _logit(T x) {
    x /= 1 - x;
    return std::log(x);
};


template<typename T>
inline T _expit(T x) {
    return 1 / (1 + std::exp(-x));
};


npy_float logitf(npy_float x)  {return _logit(x);};
npy_double logit(npy_double x) {return _logit(x);};
npy_longdouble logitl(npy_longdouble x) {return _logit(x);};

npy_float expitf(npy_float x) {return _expit(x);};
npy_double expit(npy_double x) {return _expit(x);};
npy_longdouble expitl(npy_longdouble x) {return _expit(x);};


#endif
