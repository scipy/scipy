#ifndef ELLINT_RG_GENERIC_GUARD
#define ELLINT_RG_GENERIC_GUARD


#include <algorithm>
#include <iterator>
#include <complex>
#include "ellint_typing.hh"
#include "ellint_argcheck.hh"
#include "ellint_common.hh"
#include "ellint_carlson.hh"


namespace ellint_carlson {

template<typename T, typename TR>
ExitStatus
rg(const T& x, const T& y, const T& z, const TR& rerr, T& res)
{
    typedef typing::decplx_t<T> RT;

    ExitStatus status = ExitStatus::success;
    if ( argcheck::invalid_rerr(rerr, 1.0e-4) )
    {
	res = typing::nan<T>();
	return ExitStatus::bad_rerr;
    }

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

    /* Special case -- also covers the case of z ~ zero. */
    if ( argcheck::too_small(cct[0]) && argcheck::too_small(cct[1]) )
    {
	res = std::sqrt(cct[2]) * (RT)0.5;
	return status;
    }

    T rfv, rdv;
    ExitStatus status_tmp = rf(cct[0], cct[1], cct[2], rerr * (TR)0.5, rfv);
    if ( is_horrible(status_tmp) )
    {
	res = typing::nan<T>();
	return status_tmp;
    }
    status = rd(cct[0], cct[1], cct[2], rerr * (TR)0.5, rdv);
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
    tmp += (cct[2] - cct[1]) * (cct[2] - cct[0]) * rdv / (RT)(-3.0) +
           std::sqrt(cct[0] * cct[1] / cct[2]);
    tmp *= (RT)0.5;

    res = tmp;
    return status;
}


}  /* namespace ellint_carlson */


#endif  /* ELLINT_RG_GENERIC_GUARD */
