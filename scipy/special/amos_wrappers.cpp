#include "amos_wrappers.h"
#include "special/amos.h"

double sin_pi(double x) { return special::sin_pi(x); }

double cos_pi(double x) { return special::cos_pi(x); }

int airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return special::airy_wrap(x, ai, aip, bi, bip);
}

int cairy_wrap(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    return special::cairy_wrap({npy_creal(z), npy_cimag(z)}, reinterpret_cast<std::complex<double> *>(ai),
                               reinterpret_cast<std::complex<double> *>(aip),
                               reinterpret_cast<std::complex<double> *>(bi),
                               reinterpret_cast<std::complex<double> *>(bip));
}

int cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    return special::cairy_wrap_e({npy_creal(z), npy_cimag(z)}, reinterpret_cast<std::complex<double> *>(ai),
                                 reinterpret_cast<std::complex<double> *>(aip),
                                 reinterpret_cast<std::complex<double> *>(bi),
                                 reinterpret_cast<std::complex<double> *>(bip));
}

int cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip) {
    return special::cairy_wrap_e_real(z, ai, aip, bi, bip);
}

npy_cdouble cbesi_wrap(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesi_wrap(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cbesi_wrap_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesi_wrap_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesi_wrap_e_real(double v, double z) { return special::cbesi_wrap_e_real(v, z); }

npy_cdouble cbesj_wrap(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesj_wrap(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesj_wrap_real(double v, double x) { return special::cbesj_wrap_real(v, x); }

npy_cdouble cbesj_wrap_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesj_wrap_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesj_wrap_e_real(double v, double z) { return special::cbesj_wrap_e_real(v, z); }

npy_cdouble cbesy_wrap(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesy_wrap(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesy_wrap_real(double v, double x) { return special::cbesy_wrap_real(v, x); }

npy_cdouble cbesy_wrap_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesy_wrap_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesy_wrap_e_real(double v, double z) { return special::cbesy_wrap_e_real(v, z); }

npy_cdouble cbesk_wrap(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesk_wrap(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cbesk_wrap_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesk_wrap_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double cbesk_wrap_real(double v, double z) { return special::cbesk_wrap_real(v, z); }

double cbesk_wrap_real_int(int n, double z) { return special::cbesk_wrap_real_int(n, z); }

double cbesk_wrap_e_real(double v, double z) { return special::cbesk_wrap_e_real(v, z); }

npy_cdouble cbesh_wrap1(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesh_wrap1(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cbesh_wrap1_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesh_wrap1_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cbesh_wrap2(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesh_wrap2(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cbesh_wrap2_e(double v, npy_cdouble z) {
    std::complex<double> res = special::cbesh_wrap2_e(v, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}
