#ifndef ELLINT_ARITH_AUX_HH_INCLUDED
#define ELLINT_ARITH_AUX_HH_INCLUDED


#include <cstddef>
#include "ellint_typing.hh"


/* Auxiliary floating-point manipulation utilities.
 * Ref:
 *
 * S. M. Rump, T. Ogita, S. Oishi: Accurate Floating-Point Summation.
 *      2005-11-13, Tech. Rep. 05.1, Fac. Inf. Commun. Sci, Hambg. Univ.
 *      Technol., also known as SIAM J. Sci. Comput. Volume 31, Issue 1, pp.
 *      189-224 (2008),
 *      http://www.ti3.tuhh.de/paper/rump/RuOgOi07I.pdf
 *
 * The main purpose of these templates is to implement Algorithm 6.1 (here as
 * acc_sum) which is required by compensated complex floating-point polynomial
 * evaluation.
 */
namespace ellint_carlson { namespace arithmetic { namespace aux
{
    /* Algorithm 3.6 */
    template<typename RT>
    inline typing::real_only<RT, RT>
    next_power_two(const RT& p)
    {
	RT q = p / typing::half_epsilon<RT>();
	RT L = std::abs((q + p) - q);
	if ( L == 0.0 )
	{
	    L = std::abs(p);
	}
	return L;
    }


    /* Algorithm 3.2 */
    template<typename RT>
    inline typing::real_only<RT, void>
    extract_scalar(const RT& sigma, RT& p, RT& q)
    {
	q = (sigma + p) - sigma;
	p -= q;
    }


    /* Algorithm 3.4 */
    template<typename RT, std::size_t LEN>
    inline typing::real_only<RT, void>
    extract_vector(const RT& sigma, RT(& p)[LEN], bool(& mask)[LEN], RT& tau)
    {
	RT q;
	tau = (RT)0.0;
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    if ( mask[i] )
	    {
		extract_scalar(sigma, p[i], q);
		if ( p[i] == 0.0 )
		{
		    mask[i] = false;
		}
		tau += q;
	    }
	}
    }


    /* The functions make, none, and count implements boolean-array
     * initialization, test for all-false, and counting for number of true
     * values. */
    template<typename RT, std::size_t LEN>
    inline typing::real_only<RT, void>
    make(const RT(& buf)[LEN], bool(& mask)[LEN])
    {
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    mask[i] = ( buf[i] != 0.0 );
	}
    }


    template<std::size_t LEN>
    inline bool
    none(const bool(& mask)[LEN])
    {
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    if ( mask[i] )
	    {
		return false;
	    }
	}
	return true;
    }


    template<std::size_t LEN>
    inline std::size_t
    count(const bool(& mask)[LEN])
    {
	std::size_t r = 0;
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    if ( mask[i] )
	    {
		++r;
	    }
	}
	return r;
    }


    /* Seach for maximal absolute value in masked array. */
    template<typename RT, std::size_t LEN>
    inline typing::real_only<RT, RT>
    vmax(const RT(& p)[LEN], const bool(& mask)[LEN])
    {
	RT v(0.0);
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    if ( mask[i] )
	    {
		v = std::max(v, std::abs(p[i]));
	    }
	}
	return v;
    }


    /* Summation for masked array */
    template<typename RT, std::size_t LEN>
    inline typing::real_only<RT, RT>
    masked_sum(const RT(& p)[LEN], const bool(& mask)[LEN])
    {
	RT v(0.0);
	for ( std::size_t i = 0; i < LEN; ++i )
	{
	    if ( mask[i] )
	    {
		v += p[i];
	    }
	}
	return v;
    }


    /* Algorithm 6.1 */
    template<typename RT, std::size_t LEN>
    inline typing::real_only<RT, RT>
    acc_sum(RT(& p)[LEN], bool(& mask)[LEN])
    {
	if ( (LEN == 0) || none(mask) )
	{
	    return (RT)0.0;
	}

	RT mu = vmax(p, mask);
	if ( mu == 0.0 )
	{
	    return mu;
	}

	RT twopM = next_power_two<RT>(RT(count(mask) + 2));
	RT sigma = twopM * next_power_two<RT>(mu);
	RT phi = typing::half_epsilon<RT>() * twopM;
	RT factor = std::numeric_limits<RT>::epsilon() * twopM * twopM;

	RT t = (RT)0.0;
	while ( true )
	{
	    RT tau, tau1, tau2;
	    extract_vector(sigma, p, mask, tau);
	    tau1 = t + tau;
	    if ( ( std::abs(tau1) >= factor * sigma ) ||
		 ( sigma <= std::numeric_limits<RT>::min() ) )
	    {
		/* Dekker "FastTwoSum" algorithm */
		tau2 = tau - (tau1 - t);
		return tau1 + (tau2 + masked_sum(p, mask));
	    }
	    t = tau1;
	    if ( t == 0.0 )
	    {
		return acc_sum(p, mask);
	    }
	    sigma *= phi;
	}
    }
}}}  /* namespace ellint_carlson::arithmetic::aux */


#endif /* ELLINT_ARITH_AUX_HH_INCLUDED */
