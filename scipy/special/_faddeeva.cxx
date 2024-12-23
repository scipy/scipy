#include "_faddeeva.h"

#include <complex>
#include <cmath>

using namespace std;

extern "C" {

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
        return log(Faddeeva::erfcx(-t)/2) - t*t;
    }
    else {
        return log1p(-Faddeeva::erfc(t)/2);
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

}  // extern "C"
