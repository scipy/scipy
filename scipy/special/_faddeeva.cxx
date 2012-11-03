#include <complex>

#include "_faddeeva.h"

#define RELERR 0    /* machine precision */

using namespace std;

/*
 * erf implementation (by Steven G. Johnson)
 */
complex<double> Faddeeva_erf(complex<double> z, double relerr)
{
    double x = real(z), y = imag(z);
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return (x >= 0 ? 1.0 : -1.0);
    
    complex<double> mz2(mRe_z2, mIm_z2); // -z^2
    if (x >= 0) {
        if (x < 5e-2 && fabs(y) < 5e-2)
            goto taylor;
        else if (x < 1e-3 && fabs(mIm_z2) < 1e-3)
            goto taylor_erfi;
        return 1.0 - (exp(mz2) * Faddeeva_w(complex<double>(-y,x)));
    }
    else { // x < 0
        if (x > -5e-2 && fabs(y) < 5e-2)
            goto taylor;
        else if (x > -1e-3 && fabs(mIm_z2) < 1e-3)
            goto taylor_erfi;
        return (exp(complex<double>(mRe_z2, mIm_z2))
                * Faddeeva_w(complex<double>(y,-x))) - 1.0;
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    //     erf(z) 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/41 + z^8/216...)
taylor:
    return z * (1.1283791670955125738962
                + mz2 * (0.37612638903183752464
                         + mz2 * (0.11283791670955126365338
                                  + mz2 * (0.026866170645131250268060
                                           + mz2 * 0.0052239776254421875521228))));

    /* for small |x| and small |xy|, 
       use Taylor series to avoid cancellation inaccuracy:
       erf(x+iy) = erf(iy)
       + 2*exp(y^2)/sqrt(pi) *
       [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ... 
       - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
       where:
       erf(iy) = exp(y^2) * Im[w(y)]
    */
taylor_erfi:
    if (x == 0) // important special case (also avoid NaN for large y)
        return complex<double>(x, // preserve sign of 0
                               exp(y*y) * ImFaddeeva_w(y));
    else {
        double x2 = x*x, y2 = y*y;
        double expy2 = exp(y2);
        double sexpy2 = 1.1283791670955125739*expy2; // 2*exp(y^2)/sqrt(pi)
        return complex<double>
            (sexpy2 * x * (1 - x2 * (0.33333333333333333333
                                     + 0.66666666666666666667*y2)
                           + x2*x2 * (0.1 + y2*(0.4+0.13333333333333333333*y2))),
             expy2 * ImFaddeeva_w(y)
             - sexpy2*x2*y * (1 - x2*(0.5+0.33333333333333333333*y2)));
    }
}

complex<double> Faddeeva_erfc(complex<double> z, double relerr)
{
    double x = real(z), y = imag(z);
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return (x >= 0 ? 1.0 : -1.0);
    
    complex<double> mz2(mRe_z2, mIm_z2); // -z^2
    if (x >= 0) {
        return (exp(mz2) * Faddeeva_w(complex<double>(-y,x)));
    }
    else { // x < 0
        return 2.0 - (exp(complex<double>(mRe_z2, mIm_z2))
                * Faddeeva_w(complex<double>(y,-x)));
    }
}


EXTERN_C_START

npy_cdouble faddeeva_w(npy_cdouble zp)
{
    npy_cdouble wp;
    std::complex<double> z(zp.real, zp.imag);
    std::complex<double> w = Faddeeva_w(z, RELERR);

    wp.real = real(w);
    wp.imag = imag(w);
    return wp;
}

npy_cdouble faddeeva_erf(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva_erf(z, RELERR);
    return npy_cpack(real(w), imag(w));
}

npy_cdouble faddeeva_erfc(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = Faddeeva_erfc(z, RELERR);
    return npy_cpack(real(w), imag(w));
}

double faddeeva_erfcx(double zp)
{
    return erfcx(zp);
}

npy_cdouble faddeeva_erfcx_complex(npy_cdouble zp)
{
    return faddeeva_w(npy_cpack(-zp.imag, zp.real));
}

double faddeeva_dawsn(double x)
{
    return 0.88622692545275801364908374167057 // sqrt(pi)/2
        * ImFaddeeva_w(x);
}

EXTERN_C_END
