#pragma once

#include "cephes/ndtr.h"
#include "../Faddeeva.hh"
#include "config.h"

namespace xsf {

inline double erf(double x){
    return cephes::erf(x);
}

inline float erf(float x){
    return cephes::erf(static_cast<double>(x));
}

inline std::complex<double> erf(std::complex<double> z){
    return Faddeeva::erf(z);
}

inline std::complex<float> erf(std::complex<float> x){
    return static_cast<std::complex<float>>(erf(static_cast<std::complex<double>>(x)));
}

} // namespace xsf
