#pragma once

#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/chdtr.h"
#include "xsf/cephes/fdtr.h"
#include "xsf/cephes/gdtr.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/incbi.h"
#include "xsf/cephes/kolmogorov.h"
#include "xsf/cephes/nbdtr.h"
#include "xsf/cephes/ndtr.h"
#include "xsf/cephes/ndtri.h"
#include "xsf/cephes/owens_t.h"
#include "xsf/cephes/pdtr.h"
#include "xsf/cephes/tukey.h"
#include "xsf/erf.h"

namespace xsf {

inline double bdtr(double k, int n, double p) { return cephes::bdtr(k, n, p); }

inline double bdtri(double k, int n, double y) { return cephes::bdtri(k, n, y); }

inline double bdtrc(double k, int n, double p) { return cephes::bdtrc(k, n, p); }

inline double chdtr(double df, double x) { return cephes::chdtr(df, x); }

inline double chdtrc(double df, double x) { return cephes::chdtrc(df, x); }

inline double chdtri(double df, double y) { return cephes::chdtri(df, y); }

inline double fdtr(double a, double b, double x) { return cephes::fdtr(a, b, x); }

inline double fdtrc(double a, double b, double x) { return cephes::fdtrc(a, b, x); }

inline double fdtri(double a, double b, double y) { return cephes::fdtri(a, b, y); }

inline double gdtr(double a, double b, double x) { return cephes::gdtr(a, b, x); }

inline double gdtrc(double a, double b, double x) { return cephes::gdtrc(a, b, x); }

inline double kolmogorov(double x) { return cephes::kolmogorov(x); }

inline double kolmogc(double x) { return cephes::kolmogc(x); }

inline double kolmogi(double x) { return cephes::kolmogi(x); }

inline double kolmogci(double x) { return cephes::kolmogci(x); }

inline double kolmogp(double x) { return cephes::kolmogp(x); }

inline double ndtr(double x) { return cephes::ndtr(x); }

inline float ndtr(float x) { return ndtr(static_cast<double>(x)); }

inline std::complex<double> ndtr(std::complex<double> z) { return 0.5 * erfc(-z * M_SQRT1_2); }

inline std::complex<float> ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(ndtr(static_cast<std::complex<double>>(z)));
}

/*
 * Log of the CDF of the normal distribution for double x.
 *
 * Let F(x) be the CDF of the standard normal distribution.
 * This implementation of log(F(x)) is based on the identities
 *
 *   F(x) = erfc(-x/√2)/2
 *        = 1 - erfc(x/√2)/2
 *
 * We use the first formula for x < -1, with erfc(z) replaced
 * by erfcx(z)*exp(-z**2) to ensure high precision for large
 * negative values when we take the logarithm:
 *
 *   log F(x) = log(erfc(-x/√2)/2)
 *            = log(erfcx(-x/√2)/2)*exp(-x**2/2))
 *            = log(erfcx(-x/√2)/2) - x**2/2
 *
 * For x >= -1, we use the second formula for F(x):
 *
 *   log F(x) = log(1 - erfc(x/√2)/2)
 *            = log1p(-erfc(x/√2)/2)
 */
inline double log_ndtr(double x) {
    double t = x * M_SQRT1_2;
    if (x < -1.0) {
        return log(erfcx(-t) / 2) - t * t;
    } else {
        return log1p(-erfc(t) / 2);
    }
}

inline float log_ndtr(float x) { return log_ndtr(static_cast<double>(x)); }

/*
 * Log of the normal CDF for complex arguments.
 *
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses $\erfc(z) = \exp(-z^2) w(iz)$ taking special care to select
 * the principal branch of the log function log( exp(-z^2) w(i z) )
 */
inline std::complex<double> log_ndtr(std::complex<double> z) {
    if (z.real() > 6) {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        std::complex<double> w = -0.5 * erfc(z * M_SQRT1_2);
        if (std::abs(w) < 1e-8) {
            return w;
        }
    }

    z *= -M_SQRT1_2;
    double x = std::real(z);
    double y = std::imag(z);

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$.
     */
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2 * x * y;        // Im(-z^2)

    double im = fmod(mIm_z2, 2.0 * M_PI);
    if (im > M_PI) {
        im -= 2.0 * M_PI;
    }

    std::complex<double> val1 = std::complex<double>(mRe_z2, im);

    std::complex<double> val2 = log(xsf::wofz(complex<double>(-y, x)));
    std::complex<double> result = val1 + val2 - NPY_LOGE2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = imag(result);
    if (im >= M_PI) {
        im -= 2 * M_PI;
    }
    if (im < -M_PI) {
        im += 2 * M_PI;
    }

    return {result.real(), im};
}

inline std::complex<float> log_ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(log_ndtr(static_cast<std::complex<double>>(z)));
}

inline double nbdtr(int k, int n, double p) { return cephes::nbdtr(k, n, p); }

inline double nbdtrc(int k, int n, double p) { return cephes::nbdtrc(k, n, p); }

inline double nbdtri(int k, int n, double p) { return cephes::nbdtri(k, n, p); }

inline double ndtri(double x) { return cephes::ndtri(x); }

inline double owens_t(double h, double a) { return cephes::owens_t(h, a); }

inline double pdtr(double k, double m) { return cephes::pdtr(k, m); }

inline double pdtrc(double k, double m) { return cephes::pdtrc(k, m); }

inline double pdtri(int k, double y) { return cephes::pdtri(k, y); }

inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

inline double tukeylambdacdf(double x, double lmbda) { return cephes::tukeylambdacdf(x, lmbda); }

inline float tukeylambdacdf(float x, double lmbda) {
    return tukeylambdacdf(static_cast<double>(x), static_cast<double>(lmbda));
}

} // namespace xsf
