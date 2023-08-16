#include "_faddeeva.h"

#include <complex>
#include <cmath>

using namespace std;

extern "C" {

npy_cdouble faddeeva_w(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = Faddeeva::w(z);
    return npy_cpack(real(w), imag(w));
}

npy_cdouble faddeeva_erf(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = Faddeeva::erf(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfc(double x)
{
    return Faddeeva::erfc(x);
}

npy_cdouble faddeeva_erfc_complex(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = Faddeeva::erfc(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfcx(double x)
{
    return Faddeeva::erfcx(x);
}

npy_cdouble faddeeva_erfcx_complex(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = Faddeeva::erfcx(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfi(double x)
{
    return Faddeeva::erfi(x);
}

npy_cdouble faddeeva_erfi_complex(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = Faddeeva::erfi(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_dawsn(double x)
{
    return Faddeeva::Dawson(x);
}

npy_cdouble faddeeva_dawsn_complex(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = Faddeeva::Dawson(z);
    return npy_cpack(real(w), imag(w));
}

/*
 * A wrapper for a normal CDF for complex argument
 */

npy_cdouble faddeeva_ndtr(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    z *= M_SQRT1_2;
    complex<double> w = 0.5 * Faddeeva::erfc(-z);
    return npy_cpack(real(w), imag(w));
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
double faddeeva_log_ndtr(double x)
{
    double t = x*M_SQRT1_2;
    if (x < -1.0) {
        return log(faddeeva_erfcx(-t)/2) - t*t;
    }
    else {
        return log1p(-faddeeva_erfc(t)/2);
    }
}

/*
 * Log of the normal CDF for complex arguments.
 *
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses the Faddeva computation, $\erfc(z) = \exp(-z^2) w(iz)$,
 * taking special care to select the principal branch of the log function
 *           log( exp(-z^2) w(i z) )
 */
npy_cdouble faddeeva_log_ndtr_complex(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    if (npy_creal(zp) > 6) {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        complex<double> w = -0.5 * Faddeeva::erfc(z*M_SQRT1_2);
        if (abs(w) < 1e-8) {
            return npy_cpack(real(w), imag(w));
        }
    }

    z *= -M_SQRT1_2;
    double x = real(z), y = imag(z);

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$.
     */
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)

    double im = fmod(mIm_z2, 2.0*M_PI);
    if (im > M_PI) {im -= 2.0*M_PI;}

    complex<double> val1 = complex<double>(mRe_z2, im);

    complex<double> val2 = log(Faddeeva::w(complex<double>(-y, x)));
    complex<double> result = val1 + val2 - NPY_LOGE2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = imag(result);
    if (im >= M_PI){ im -= 2*M_PI; }
    if (im < -M_PI){ im += 2*M_PI; }

    return npy_cpack(real(result), im);
}

double faddeeva_voigt_profile(double x, double sigma, double gamma)
{
    const double INV_SQRT_2 = 0.707106781186547524401;
    const double SQRT_2PI = 2.5066282746310002416123552393401042;

    if(sigma == 0){
        if (gamma == 0){
            if (std::isnan(x))
                return x;
            if (x == 0)
                return INFINITY;
            return 0;
        }
        return gamma / M_PI / (x*x + gamma*gamma);
    }
    if (gamma == 0){
        return 1 / SQRT_2PI / sigma * exp(-(x/sigma)*(x/sigma) / 2);
    }

    double zreal = x / sigma * INV_SQRT_2;
    double zimag = gamma / sigma * INV_SQRT_2;
    std::complex<double> z(zreal, zimag);
    std::complex<double> w = Faddeeva::w(z);
    return real(w) / sigma / SQRT_2PI;
}

}  // extern "C"
