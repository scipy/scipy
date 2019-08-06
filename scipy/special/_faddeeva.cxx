#include <complex>

#include "_faddeeva.h"

using namespace std;

EXTERN_C_START

npy_cdouble faddeeva_w(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    std::complex<double> w = Faddeeva::w(z);
    return npy_cpack(real(w), imag(w));
}

npy_cdouble faddeeva_erf(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva::erf(z);
    return npy_cpack(real(w), imag(w));
}

npy_cdouble faddeeva_erfc(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva::erfc(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfcx(double x)
{
    return Faddeeva::erfcx(x);
}

npy_cdouble faddeeva_erfcx_complex(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva::erfcx(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfi(double x)
{
    return Faddeeva::erfi(x);
}

npy_cdouble faddeeva_erfi_complex(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva::erfi(z);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_dawsn(double x)
{
    return Faddeeva::Dawson(x);
}

npy_cdouble faddeeva_dawsn_complex(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva::Dawson(z);
    return npy_cpack(real(w), imag(w));
}

/*
 * A wrapper for a normal CDF for complex argument
 */

npy_cdouble faddeeva_ndtr(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    z *= NPY_SQRT1_2;
    complex<double> w = 0.5 * Faddeeva::erfc(-z);
    return npy_cpack(real(w), imag(w));
}

/*
 * Log of the normal CDF for complex arguments.
 * 
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses the Faddeva computation, $\erfc(z) = \exp(-z^2) w(iz)$,
 * taking special care to select the principal branch of the log function
 *           log( exp(-z^2) w(i z) )
 */
npy_cdouble faddeeva_log_ndtr(npy_cdouble zp)
{

    complex<double> z(zp.real, zp.imag);
    if (zp.real > 6) {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        complex<double> w = -0.5 * Faddeeva::erfc(z*NPY_SQRT1_2);
        if (abs(w) < 1e-8) {
            return npy_cpack(real(w), imag(w));
        }
    }

    z *= -NPY_SQRT1_2;
    double x = real(z), y = imag(z);

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$. 
     */
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)

    double im = fmod(mIm_z2, 2.0*NPY_PI);
    if (im > NPY_PI) {im -= 2.0*NPY_PI;}

    complex<double> val1 = complex<double>(mRe_z2, im);

    complex<double> val2 = log(Faddeeva::w(complex<double>(-y, x)));
    complex<double> result = val1 + val2 - NPY_LOGE2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = imag(result);
    if (im >= NPY_PI){ im -= 2*NPY_PI; }
    if (im < -NPY_PI){ im += 2*NPY_PI; }

    return npy_cpack(real(result), im);
}

double faddeeva_voigt(double x, double sigma, double gamma, double mu)
{
    const double INV_SQRT_2 = 0.707106781186547524401;
    const double SQRT_2PI = 2.5066282746310002416123552393401042;

    double zreal = (x-mu) / sigma * INV_SQRT_2;
    double zimag = gamma / sigma * INV_SQRT_2;
    std::complex<double> z(zreal, zimag);
    std::complex<double> w = Faddeeva::w(z);
    return real(w) / sigma / SQRT_2PI;
}

EXTERN_C_END
