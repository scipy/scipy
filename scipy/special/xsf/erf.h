#pragma once

#include "../Faddeeva.hh"
#include "cephes/ndtr.h"
#include "config.h"

namespace xsf {

inline double erf(double x) { return cephes::erf(x); }

inline float erf(float x) { return erf(static_cast<double>(x)); }

inline std::complex<double> erf(std::complex<double> z) { return Faddeeva::erf(z); }

inline std::complex<float> erf(std::complex<float> x) {
    return static_cast<std::complex<float>>(erf(static_cast<std::complex<double>>(x)));
}

inline double erfc(double x) { return cephes::erfc(x); }

inline float erfc(float x) { return erfc(static_cast<double>(x)); }

inline std::complex<double> erfc(std::complex<double> z) { return Faddeeva::erfc(z); }

inline std::complex<float> erfc(std::complex<float> x) {
    return static_cast<std::complex<float>>(erfc(static_cast<std::complex<double>>(x)));
}

inline double erfcx(double x) { return Faddeeva::erfcx(x); }

inline float erfcx(float x) { return erfcx(static_cast<double>(x)); }

inline std::complex<double> erfcx(std::complex<double> z) { return Faddeeva::erfcx(z); }

inline std::complex<float> erfcx(std::complex<float> x) {
    return static_cast<std::complex<float>>(erfcx(static_cast<std::complex<double>>(x)));
}

} // namespace xsf
