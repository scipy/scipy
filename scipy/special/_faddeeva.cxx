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

EXTERN_C_END
