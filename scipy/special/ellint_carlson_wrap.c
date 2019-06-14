#include "ellint_carlson_wrap.h"

npy_double fellint_RC_w(npy_double x, npy_double y)
{
    int status;
    npy_double wres;

    status = fellint_RC(x, y, SF_RERR, &wres);
    sf_error("ellint_RC", status, NULL);

    return wres;
}

npy_cdouble cellint_RC_w(npy_cdouble x, npy_cdouble y)
{
    int status;
    double complex res;
    npy_cdouble wres;

    status = cellint_RC(ELLINT_TO_C(x), ELLINT_TO_C(y), SF_RERR, &res);
    sf_error("ellint_RC", status, NULL);

    ELLINT_MAKE_N(wres, res);
    return wres;
}

npy_double fellint_RD_w(npy_double x, npy_double y, npy_double z)
{
    int status;
    npy_double wres;

    status = fellint_RD(x, y, z, SF_RERR, &wres);
    sf_error("ellint_RD", status, NULL);

    return wres;
}

npy_cdouble cellint_RD_w(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    int status;
    double complex res;
    npy_cdouble wres;

    status = cellint_RD(ELLINT_TO_C(x), ELLINT_TO_C(y), ELLINT_TO_C(z),
                        SF_RERR, &res);
    sf_error("ellint_RD", status, NULL);

    ELLINT_MAKE_N(wres, res);
    return wres;
}

npy_double fellint_RF_w(npy_double x, npy_double y, npy_double z)
{
    int status;
    npy_double wres;

    status = fellint_RF(x, y, z, SF_RERR, &wres);
    sf_error("ellint_RF", status, NULL);

    return wres;
}

npy_cdouble cellint_RF_w(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    int status;
    double complex res;
    npy_cdouble wres;

    status = cellint_RF(ELLINT_TO_C(x), ELLINT_TO_C(y), ELLINT_TO_C(z),
                        SF_RERR, &res);
    sf_error("ellint_RF", status, NULL);

    ELLINT_MAKE_N(wres, res);
    return wres;
}

npy_double fellint_RG_w(npy_double x, npy_double y, npy_double z)
{
    int status;
    npy_double wres;

    status = fellint_RG(x, y, z, SF_RERR, &wres);
    sf_error("ellint_RG", status, NULL);

    return wres;
}

npy_cdouble cellint_RG_w(npy_cdouble x, npy_cdouble y, npy_cdouble z)
{
    int status;
    double complex res;
    npy_cdouble wres;

    status = cellint_RG(ELLINT_TO_C(x), ELLINT_TO_C(y), ELLINT_TO_C(z),
                        SF_RERR, &res);
    sf_error("ellint_RG", status, NULL);

    ELLINT_MAKE_N(wres, res);
    return wres;
}

npy_double fellint_RJ_w(npy_double x, npy_double y, npy_double z, npy_double p)
{
    int status;
    npy_double wres;

    status = fellint_RJ(x, y, z, p, SF_RERR, &wres);
    sf_error("ellint_RJ", status, NULL);

    return wres;
}

npy_cdouble cellint_RJ_w(npy_cdouble x, npy_cdouble y, npy_cdouble z,
                         npy_cdouble p)
{
    int status;
    double complex res;
    npy_cdouble wres;

    status = cellint_RJ(ELLINT_TO_C(x), ELLINT_TO_C(y), ELLINT_TO_C(z),
                        ELLINT_TO_C(p), SF_RERR, &res);
    sf_error("ellint_RJ", status, NULL);

    ELLINT_MAKE_N(wres, res);
    return wres;
}
