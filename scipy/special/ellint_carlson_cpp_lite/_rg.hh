#ifndef ELLINT_RG_GENERIC_GUARD
#define ELLINT_RG_GENERIC_GUARD


#include <algorithm>
#include <iterator>
#include <complex>
#include "ellint_typing.hh"
#include "ellint_argcheck.hh"
#include "ellint_common.hh"
#include "ellint_carlson.hh"


/* Forward declaration */
namespace ellint_carlson {
    template<typename T>
    inline void agm_update(T& x, T& y);
}


namespace ellint_carlson {

namespace arithmetic { namespace aux {

template<typename T>
static inline typing::real_only<T, void>
rg_dot2_acc(const T& fac, const T& term, T& acc, T& cor)
{
    fdot2_acc(fac, term, acc, cor);
}

template<typename CT>
static inline typing::cplx_only<CT, void>
rg_dot2_acc(const typing::decplx_t<CT>& fac, const CT& term, CT& acc, CT& cor)
{
    typedef typing::decplx_t<CT> RT;
    RT ar, cr, ai, ci;
    CT p, q;
    eft_prod(fac, term.real(), ar, cr);
    eft_prod(fac, term.imag(), ai, ci);
    eft_sum(acc, CT{ar, ai}, p, q);
    acc = p;
    cor += q + CT{cr, ci};
}

}}  /* namespace ellint_carlson::arithmetic::aux */


template<typename T>
static ExitStatus
rg0(const T& x, const T& y, const double& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;
    ExitStatus status = ExitStatus::success;
    double rsq = 2.0 * std::sqrt(rerr);
    RT fac(0.25);

    T xm = std::sqrt(x);
    T ym = std::sqrt(y);
    T dm = (xm + ym) * (RT)0.5;
    T sum = -dm * dm;
    T cor(0.0);
    dm = xm - ym;
    unsigned int m = 0;
    while ( std::abs(dm) >= rsq * std::fmin(std::abs(xm), std::abs(ym)) )
    {

	if ( m > config::max_iter )
	{
	    status = ExitStatus::n_iter;
	    break;
	}

	agm_update(xm, ym);
	fac *= (RT)2.0;
	dm = xm - ym;
	arithmetic::aux::rg_dot2_acc(fac, dm * dm, sum, cor);

	++m;
    }

    res = (RT)(constants::pi) / (xm + ym);
    res *= -(RT)0.5 * (sum + cor);
    return status;
}


template<typename T>
ExitStatus
rg(const T& x, const T& y, const T& z, const double& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;

    ExitStatus status = ExitStatus::success;
#ifndef ELLINT_NO_VALIDATE_RELATIVE_ERROR_BOUND
    if ( argcheck::invalid_rerr(rerr, 1.0e-4) )
    {
	res = typing::nan<T>();
	return ExitStatus::bad_rerr;
    }
#endif

    T cct[3] = {x, y, z};
    std::sort(std::begin(cct), std::end(cct), util::abscmp<T>);

    if ( (argcheck::isinf(cct[0]) || argcheck::isinf(cct[1]) ||
          argcheck::isinf(cct[2])) &&
	 argcheck::ph_good(cct[0]) && argcheck::ph_good(cct[1]) &&
	 argcheck::ph_good(cct[2]) )
    {
	res = typing::huge<T>();
	return ExitStatus::singular;
    }

    if ( argcheck::too_small(cct[0]) )
    {
	/* Special case -- also covers the case of z ~ zero. */
	if ( argcheck::too_small(cct[1]) )
	{
	    res = std::sqrt(cct[2]) * (RT)0.5;
	    return status;
	} else {
	    /* Special case -- use the AGM algorithm. */
	    status = rg0(cct[1], cct[2], rerr, res);
	    return status;
	}
    }

    T rfv, rdv;
    ExitStatus status_tmp = rf(cct[0], cct[1], cct[2], rerr * 0.5, rfv);
    if ( is_horrible(status_tmp) )
    {
	res = typing::nan<T>();
	return status_tmp;
    }
    status = rd(cct[0], cct[1], cct[2], rerr * 0.5, rdv);
    if ( status_tmp != ExitStatus::success )
    {
	status = status_tmp;
    }
    if ( is_horrible(status) )
    {
	res = typing::nan<T>();
	return status;
    }

    T tmp = cct[2] * rfv;
    tmp += (cct[1] - cct[2]) * (cct[2] - cct[0]) * rdv / (RT)3.0 +
           std::sqrt(cct[0] * cct[1] / cct[2]);
    tmp *= (RT)0.5;

    res = tmp;
    return status;
}


}  /* namespace ellint_carlson */


#endif  /* ELLINT_RG_GENERIC_GUARD */
