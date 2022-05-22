#ifndef ELLINT_RF_GENERIC_GUARD
#define ELLINT_RF_GENERIC_GUARD


#include <algorithm>
#include <iterator>
#include <cmath>
#include <complex>
#include "ellint_typing.hh"
#include "ellint_argcheck.hh"
#include "ellint_common.hh"
#include "ellint_carlson.hh"


namespace ellint_carlson {

template<typename T>
inline void agm_update(T& x, T& y)
{
    typedef typing::decplx_t<T> RT;

    T xnext = (x + y) * (RT)0.5;
    T ynext = std::sqrt(x * y);

    x = xnext;
    y = ynext;
}

template<typename T>
static ExitStatus
rf0(const T& x, const T& y, const double& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;
    ExitStatus status = ExitStatus::success;
    double rsq = 2.0 * std::sqrt(rerr);

    T xm = std::sqrt(x);
    T ym = std::sqrt(y);
    unsigned int m = 0;
    while ( std::abs(xm - ym) >= rsq * std::fmin(std::abs(xm), std::abs(ym)) )
    {
	if ( m > config::max_iter )
	{
	    status = ExitStatus::n_iter;
	    break;
	}

	agm_update(xm, ym);

	++m;
    }

    res = (RT)(constants::pi) / (xm + ym);
    return status;
}


template<typename T>
ExitStatus
rf(const T& x, const T& y, const T& z, const double& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;

    T cct1[3];
    T cct2[3];

    ExitStatus status = ExitStatus::success;
#ifndef ELLINT_NO_VALIDATE_RELATIVE_ERROR_BOUND
    if ( argcheck::invalid_rerr(rerr, 3.0e-4) )
    {
	res = typing::nan<T>();
	return ExitStatus::bad_rerr;
    }
#endif

    if ( argcheck::ph_good(x) && argcheck::ph_good(y) && argcheck::ph_good(z) )
    {
	if ( argcheck::isinf(x) || argcheck::isinf(y) || argcheck::isinf(z) )
	{
	    res = T(0.0);
	    return ExitStatus::success;
	}
    } else {
	res = typing::nan<T>();
	return ExitStatus::bad_args;
    }

    cct1[0] = T(x);
    cct1[1] = T(y);
    cct1[2] = T(z);
    std::sort(std::begin(cct1), std::end(cct1), util::abscmp<T>);
    T xm = cct1[0];
    T ym = cct1[1];
    T zm = cct1[2];
    if ( argcheck::too_small(xm) )
    {
	if ( argcheck::too_small(ym) )
	{
	    status = ExitStatus::singular;
	    res = typing::huge<T>();
	    return status;
	} else {
	    T tmpres;
	    status = rf0(ym, zm, rerr, tmpres);
	    res = tmpres - (std::sqrt(xm) / std::sqrt(ym * zm));
	    return status;
	}
    }

    T Am = arithmetic::sum2(cct1) / (RT)3.0;
    T xxm = Am - xm;
    T yym = Am - ym;
    RT fterm = std::abs(std::max({xxm, yym, Am - zm}, util::abscmp<T>)) /
               arithmetic::ocrt(3.0 * rerr);
    unsigned int m = 0;
    RT aAm;
    while ( (aAm = std::abs(Am)) <= fterm ||
	    aAm <= std::abs(std::max({xxm, yym, Am - zm}, util::abscmp<T>)) )
    {
	if ( m > config::max_iter )
	{
	    status = ExitStatus::n_iter;
	    break;
	}
	cct1[0] = cct2[2] = std::sqrt(xm);
	cct1[1] = cct2[0] = std::sqrt(ym);
	cct1[2] = cct2[1] = std::sqrt(zm);
	T lam = arithmetic::dot2(cct1, cct2);
        Am = (Am + lam) * (RT)0.25;
        xm = (xm + lam) * (RT)0.25;
        ym = (ym + lam) * (RT)0.25;
        zm = (zm + lam) * (RT)0.25;
        xxm *= (RT)0.25;
        yym *= (RT)0.25;
        fterm *= (RT)0.25;

        ++m;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid. */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    Am = arithmetic::sum2(cct1) / (RT)3.0;
    xxm /= Am;
    yym /= Am;
    T zzm = -(xxm + yym);
    T e2 = xxm * yym - zzm * zzm;
    T e3 = xxm * (yym * zzm);
    T s = arithmetic::comp_horner(e2, constants::RF_C1);
    s += e3 * (arithmetic::comp_horner(e2, constants::RF_C2) +
               e3 * (RT)(constants::RF_c33));
    s /= (RT)(constants::RF_DENOM);
    s += (RT)1.0;

    res = s / std::sqrt(Am);
    return status;
}


}  /* namespace ellint_carlson */


#endif /* ELLINT_RF_GENERIC_GUARD */
