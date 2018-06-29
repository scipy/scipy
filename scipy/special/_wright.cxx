#include <complex>

#include "_wright.h"

using namespace std;

EXTERN_C_START

npy_cdouble wrightomega(npy_cdouble zp)
{
    complex<double> z(zp.real, zp.imag);
    complex<double> w = wright::wrightomega(z);
    return npy_cpack(real(w), imag(w));
}

EXTERN_C_END
