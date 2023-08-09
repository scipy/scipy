#include <complex>

#include "_wright.h"

using namespace std;

extern "C" {

npy_cdouble wrightomega(npy_cdouble zp)
{
    complex<double> z(npy_creal(zp), npy_cimag(zp));
    complex<double> w = wright::wrightomega(z);
    return npy_cpack(real(w), imag(w));
}

double wrightomega_real(double x)
{
    return wright::wrightomega_real(x);
}

}  // extern "C"
