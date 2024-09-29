#include "ellint_carlson_wrap.hh"
#include "sf_error.h"


static constexpr double ellip_rerr = 5e-16;


extern "C" {


double fellint_RC(double x, double y)
{
    sf_error_t status;
    double res;

    status = static_cast<sf_error_t>(ellint_carlson::rc(x, y,
                                                        ellip_rerr, res));
    sf_error("elliprc (real)", status, NULL);
    return res;
}

npy_cdouble cellint_RC(npy_cdouble x, npy_cdouble y)
{
    sf_error_t status;
    std::complex<double> xx{npy_creal(x), npy_cimag(x)};
    std::complex<double> yy{npy_creal(y), npy_cimag(y)};
    std::complex<double> res;

    status = static_cast<sf_error_t>(ellint_carlson::rc(xx, yy,
                                                        ellip_rerr, res));
    sf_error("elliprc (complex)", status, NULL);
    return npy_cpack(res.real(), res.imag());
}


double fellint_RD(double x, double y, double z)
{
    sf_error_t status;
    double res;

    status = static_cast<sf_error_t>(ellint_carlson::rd(x, y, z,
                                                        ellip_rerr, res));
    sf_error("elliprd (real)", status, NULL);
    return res;
}

npy_cdouble cellint_RD(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    sf_error_t status;
    std::complex<double> xx{npy_creal(x), npy_cimag(x)};
    std::complex<double> yy{npy_creal(y), npy_cimag(y)};
    std::complex<double> zz{npy_creal(z), npy_cimag(z)};
    std::complex<double> res;

    status = static_cast<sf_error_t>(ellint_carlson::rd(xx, yy, zz,
                                                        ellip_rerr, res));
    sf_error("elliprd (complex)", status, NULL);
    return npy_cpack(res.real(), res.imag());
}


double fellint_RF(double x, double y, double z)
{
    sf_error_t status;
    double res;

    status = static_cast<sf_error_t>(ellint_carlson::rf(x, y, z,
                                                        ellip_rerr, res));
    sf_error("elliprf (real)", status, NULL);
    return res;
}

npy_cdouble cellint_RF(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    sf_error_t status;
    std::complex<double> xx{npy_creal(x), npy_cimag(x)};
    std::complex<double> yy{npy_creal(y), npy_cimag(y)};
    std::complex<double> zz{npy_creal(z), npy_cimag(z)};
    std::complex<double> res;

    status = static_cast<sf_error_t>(ellint_carlson::rf(xx, yy, zz,
                                                        ellip_rerr, res));
    sf_error("elliprf (complex)", status, NULL);
    return npy_cpack(res.real(), res.imag());
}


double fellint_RG(double x, double y, double z)
{
    sf_error_t status;
    double res;

    status = static_cast<sf_error_t>(ellint_carlson::rg(x, y, z,
                                                        ellip_rerr, res));
    sf_error("elliprg (real)", status, NULL);
    return res;
}

npy_cdouble cellint_RG(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    sf_error_t status;
    std::complex<double> xx{npy_creal(x), npy_cimag(x)};
    std::complex<double> yy{npy_creal(y), npy_cimag(y)};
    std::complex<double> zz{npy_creal(z), npy_cimag(z)};
    std::complex<double> res;

    status = static_cast<sf_error_t>(ellint_carlson::rg(xx, yy, zz,
                                                        ellip_rerr, res));
    sf_error("elliprg (complex)", status, NULL);
    return npy_cpack(res.real(), res.imag());
}


double fellint_RJ(double x, double y, double z, double p)
{
    sf_error_t status;
    double res;

    status = static_cast<sf_error_t>(ellint_carlson::rj(x, y, z, p,
                                                        ellip_rerr, res));
    sf_error("elliprj (real)", status, NULL);
    return res;
}

npy_cdouble cellint_RJ(npy_cdouble x, npy_cdouble y, npy_cdouble z, npy_cdouble p)
{
    sf_error_t status;
    std::complex<double> xx{npy_creal(x), npy_cimag(x)};
    std::complex<double> yy{npy_creal(y), npy_cimag(y)};
    std::complex<double> zz{npy_creal(z), npy_cimag(z)};
    std::complex<double> pp{npy_creal(p), npy_cimag(p)};
    std::complex<double> res;

    status = static_cast<sf_error_t>(ellint_carlson::rj(xx, yy, zz, pp,
                                                        ellip_rerr, res));
    sf_error("elliprj (complex)", status, NULL);
    return npy_cpack(res.real(), res.imag());
}


}  // extern "C"
