#include <complex>

#include <Python.h>

extern "C" {
#include <numpy/npy_math.h>
}

extern std::complex<double> Faddeeva_w(std::complex<double> z, double relerr=0);

extern "C" {

#define RELERR 0    /* machine precision */

npy_cdouble wofz(npy_cdouble zp)
{
    npy_cdouble wp;
    std::complex<double> z(zp.real, zp.imag);
    std::complex<double> w = Faddeeva_w(z, RELERR);

    wp.real = real(w);
    wp.imag = imag(w);
    return wp;
}

}
