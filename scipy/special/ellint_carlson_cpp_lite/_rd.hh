#ifndef ELLINT_RD_GENERIC_GUARD
#define ELLINT_RD_GENERIC_GUARD


#include <algorithm>
#include <cmath>
#include <complex>
#include "ellint_typing.hh"
#include "ellint_argcheck.hh"
#include "ellint_common.hh"
#include "ellint_carlson.hh"


namespace ellint_carlson {

template<typename T>
ExitStatus
rd(const T& x, const T& y, const T& z, const double& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;

    T cct1[6];
    T cct2[6];

    ExitStatus status = ExitStatus::success;
    if ( argcheck::invalid_rerr(rerr, 1.0e-4) )
    {
	res = typing::nan<T>();
	return ExitStatus::bad_rerr;
    }

    if ( argcheck::too_small(z) ||
         !( argcheck::ph_good(x) && argcheck::ph_good(y) &&
            argcheck::ph_good(z) ) )
    {
	res = typing::nan<T>();
	return ExitStatus::bad_args;
    }

    if ( argcheck::isinf(x) || argcheck::isinf(y) || argcheck::isinf(z) )
    {
	res = T(0.0);
	return status;
    }

    if ( argcheck::too_small(x) && argcheck::too_small(y) )
    {
	res = typing::huge<T>();
	return ExitStatus::singular;
    }

    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    cct1[3] = z;
    cct1[4] = z;
    T Am = arithmetic::nsum2(cct1, 5) / (RT)5.0;
    T xm = x;
    T ym = y;
    T zm = z;
    T xxm = Am - xm;
    T yym = Am - ym;
    RT fterm = std::max({std::abs(xxm), std::abs(yym), std::abs(Am - z)}) /
               arithmetic::ocrt(rerr / 5.0);
    RT d4m(1.0);
    T adt(0.0), ade(0.0);

    unsigned int m = 0;
    RT aAm;
    T tmp;
    while ( (aAm = std::abs(Am)) <= fterm ||
	    aAm <= std::max({std::abs(xxm), std::abs(yym), std::abs(Am - zm)}) )
    {
	if ( m > config::max_iter )
	{
	    status = ExitStatus::n_iter;
	    break;
	}
	cct1[0] = cct2[2] = std::sqrt(xm);
	cct1[1] = cct2[0] = std::sqrt(ym);
	cct1[2] = cct2[1] = std::sqrt(zm);
	T lam = arithmetic::ndot2(cct1, cct2, 3);

	tmp = d4m / (cct1[2] * (zm + lam));
	arithmetic::sum2_acc(tmp, adt, ade);

        Am = (Am + lam) * (RT)0.25;
        xm = (xm + lam) * (RT)0.25;
        ym = (ym + lam) * (RT)0.25;
        zm = (zm + lam) * (RT)0.25;
        xxm *= (RT)0.25;
        yym *= (RT)0.25;
        fterm *= (RT)0.25;
        d4m *= (RT)0.25;

        ++m;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = zm;
    cct1[4] = zm;
    Am = arithmetic::nsum2(cct1, 5) / (RT)5.0;
    xxm /= Am;
    yym /= Am;
    T zzm = (xxm + yym) / (RT)(-3.0);
    T xy = xxm * yym;
    T zz2 = zzm * zzm;
    T e2 = xy - zz2 * (RT)6.0;
    T e3 = (xy * (RT)3.0 - zz2 * (RT)8.0) * zzm;
    T e4 = (xy - zz2) * zz2 * (RT)3.0;
    T e5 = xy * zz2 * zzm;
    T t = std::sqrt(Am);
    tmp = d4m / (t * t * t);
    cct1[0] = arithmetic::comp_horner(e2, constants::RDJ_C1);
    cct1[1] = arithmetic::comp_horner(e3, constants::RDJ_C2);
    cct1[2] = arithmetic::comp_horner(e2, constants::RDJ_C3);
    cct1[3] = arithmetic::comp_horner(e2, constants::RDJ_C4);
    cct1[4] = arithmetic::comp_horner(e2, constants::RDJ_C5);
    cct1[5] = e3 * (RT)(constants::RDJ_C5[1]);

    cct2[0] = T(1.0);
    cct2[1] = T(1.0);
    cct2[2] = e3;
    cct2[3] = e4;
    cct2[4] = e5;
    cct2[5] = e4;
    t = arithmetic::dot2(cct1, cct2) / (RT)(constants::RDJ_DENOM) + 1.0;
    tmp *= t;
    tmp += adt * (RT)3.0;

    res = tmp;
    return status;
}

}  /* namespace ellint_carlson */
#endif /* ELLINT_RD_GENERIC_GUARD */
