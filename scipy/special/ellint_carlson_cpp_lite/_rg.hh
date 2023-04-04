#ifndef ELLINT_RG_GENERIC_GUARD
#define ELLINT_RG_GENERIC_GUARD


#include <algorithm>
#include <iterator>
#include <complex>
#include "ellint_typing.hh"
#include "ellint_argcheck.hh"
#include "ellint_common.hh"
#include "ellint_carlson.hh"


#define CHECK_STATUS_OR_FAIL()	\
do {	\
    if ( status_tmp != ExitStatus::success )	\
    {	\
	status = status_tmp;	\
    }	\
    if ( is_horrible(status) )	\
    {	\
	res = typing::nan<T>();	\
	return status;	\
    }	\
} while ( 0 )


/* References
 * [1] B. C. Carlson, "Numerical computation of real or complex elliptic
 *     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
 *     https://arxiv.org/abs/math/9409227
 *     https://doi.org/10.1007/BF02198293
 * [2] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
 *     Functions," NIST, US Dept. of Commerce.
 *     https://dlmf.nist.gov/19.16.E1
 *     https://dlmf.nist.gov/19.20.ii
 */


/* Forward declaration */
/* See the file _rf.hh for the definition of agm_update */
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
	/* Ref[1], Eq. 2.39 or Eq. (46) in the arXiv preprint */
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
	if ( argcheck::too_small(cct[1]) )
	{
	    /* Special case -- also covers the case of z ~ zero. */
	    res = std::sqrt(cct[2]) * (RT)0.5;
	    return status;
	} else {
	    /* Special case -- use the AGM algorithm. */
	    status = rg0(cct[1], cct[2], rerr, res);
	    return status;
	}
    }

    /* Ref[2], Eq. 19.21.11 (second identity) <https://dlmf.nist.gov/19.21.E11>
     * i.e. 6R_G() = sum of [ x * (y + z) * R_D() ] over cyclic permutations of
     * (x, y, z).
     * This cyclic form is manifestly symmetric and is preferred over
     * Eq. 19.21.10 ibid.
     * Here we put the three R_D terms in the buffer cct2, and the
     * x * (y + z) = dot({x, x}, {y, z}) terms in cct1. */
    T cct1[3];
    T cct2[3];

    ExitStatus status_tmp;

    status_tmp = rd(y, z, x, rerr, cct2[0]);
    CHECK_STATUS_OR_FAIL();
    status_tmp = rd(z, x, y, rerr, cct2[1]);
    CHECK_STATUS_OR_FAIL();
    status_tmp = rd(x, y, z, rerr, cct2[2]);
    CHECK_STATUS_OR_FAIL();

    /* Fill the cct1 buffer via the intermediate buffers tm1, tm2
     * that are dotted together (using the compensation algorithm). */
    T tm1[2] = {x, x};
    T tm2[2] = {y, z};
    cct1[0] = arithmetic::dot2(tm1, tm2);
    tm1[0] = tm1[1] = y;
    tm2[0] = x;
    cct1[1] = arithmetic::dot2(tm1, tm2);
    tm1[0] = tm1[1] = z;
    tm2[1] = y;
    cct1[2] = arithmetic::dot2(tm1, tm2);

    res = arithmetic::dot2(cct1, cct2) / (RT)6.0;
    return status;
}


}  /* namespace ellint_carlson */


#endif  /* ELLINT_RG_GENERIC_GUARD */
