#pragma once

#include "faddeeva.h"
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

inline double erfi(double x) { return Faddeeva::erfi(x); }

inline float erfi(float x) { return erfi(static_cast<double>(x)); }

inline std::complex<double> erfi(std::complex<double> z) { return Faddeeva::erfi(z); }

inline std::complex<float> erfi(std::complex<float> z) {
    return static_cast<std::complex<float>>(erfi(static_cast<std::complex<double>>(z)));
}

inline double voigt_profile(double x, double sigma, double gamma) {
    const double INV_SQRT_2 = 0.707106781186547524401;
    const double SQRT_2PI = 2.5066282746310002416123552393401042;

    if (sigma == 0) {
        if (gamma == 0) {
            if (std::isnan(x))
                return x;
            if (x == 0)
                return INFINITY;
            return 0;
        }
        return gamma / M_PI / (x * x + gamma * gamma);
    }
    if (gamma == 0) {
        return 1 / SQRT_2PI / sigma * exp(-(x / sigma) * (x / sigma) / 2);
    }

    double zreal = x / sigma * INV_SQRT_2;
    double zimag = gamma / sigma * INV_SQRT_2;
    std::complex<double> z(zreal, zimag);
    std::complex<double> w = Faddeeva::w(z);
    return w.real() / sigma / SQRT_2PI;
}

inline float voigt_profile(float x, float sigma, float gamma) {
    return voigt_profile(static_cast<double>(x), static_cast<double>(sigma), static_cast<double>(gamma));
}

inline std::complex<double> wofz(std::complex<double> z) { return Faddeeva::w(z); }

inline std::complex<float> wofz(std::complex<float> x) {
    return static_cast<std::complex<float>>(wofz(static_cast<std::complex<double>>(x)));
}

inline double dawsn(double x) { return Faddeeva::Dawson(x); }

inline float dawsn(float x) { return dawsn(static_cast<double>(x)); }

inline std::complex<double> dawsn(std::complex<double> z) { return Faddeeva::Dawson(z); }

inline std::complex<float> dawsn(std::complex<float> x) {
    return static_cast<std::complex<float>>(dawsn(static_cast<std::complex<double>>(x)));
}

} // namespace xsf
