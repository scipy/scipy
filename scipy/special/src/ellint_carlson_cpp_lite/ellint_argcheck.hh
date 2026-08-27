#ifndef ELLINT_ARGCHECK_HH_INCLUDED
#define ELLINT_ARGCHECK_HH_INCLUDED


#include <complex>
#include "ellint_typing.hh"
#include "ellint_common.hh"


namespace ellint_carlson { namespace argcheck
{
    template<typename T>
    inline constexpr typing::real_only<T, bool>
    invalid_rerr(const T& r, const double upper)
    {
	return ( r <= 0.0 ) || ( r > upper );
    }


    template<typename T>
    inline typing::real_only<T, bool>
    too_small(const T& x)
    {
	return ( (x == 0.0 ) || ( std::fpclassify(x) == FP_SUBNORMAL ) );
    }

    template<typename T>
    inline typing::cplx_only<T, bool>
    too_small(const T& z)
    {
	return ( too_small(z.real()) && too_small(z.imag()) );
    }


    /* Argument (phase) angle is NOT +/- pi, and is not indeterminate value,
     * checked without computing the actual argument. */
    template<typename T>
    inline constexpr typing::real_only<T, bool>
    ph_good(const T& x)
    {
	/* also handles x being NaN, which compares to false, hence is the
	 * correct result. */
	return ( x >= 0.0 );
    }

    template<typename T>
    inline typing::cplx_only<T, bool>
    ph_good(const T& z)
    {
	typedef typing::decplx_t<T> RT;
	switch ( std::fpclassify(z.imag()) )
	{
	    case FP_NORMAL :
	    case FP_SUBNORMAL :
	    {
		RT r = z.real();
		return ( std::isfinite(r) ||
			(std::isinf(r) && (r > 0.0)) );
	    }
	    case FP_ZERO : return ph_good<RT>(z.real());
	    case FP_NAN : return false;
	    case FP_INFINITE : return ( std::isfinite(z.real()) );
	    default /* other implementations? */ : return false;
	}
    }


    template<typename T>
    inline typing::real_only<T, bool>
    isinf(const T& x)
    {
	return std::isinf(x);
    }

    template<typename T>
    inline typing::cplx_only<T, bool>
    isinf(const T& z)
    {
	return ( std::isinf(z.real()) || std::isinf(z.imag()) );
    }


    /* Two complex numbers are approximately each other's conjugate */
    template<typename T0, typename T1>
    inline bool
    isconj(const T0& x, const T1& y)
    {
	return too_small(x - std::conj(y));
    }

    /* "let two of the variables x, y, z be nonzero and conjugate complex with
     * phase less in magnitude than pi and the third variable be real and
     * nonnegative." */
    template<typename T0, typename T1, typename T2>
    inline bool
    r1conj2(T0 x, T1 y, T2 z)
    {
	return ( isconj(x, y) && !(too_small(x) || too_small(y)) &&
	         too_small(std::imag(z)) && (std::real(z) >= 0.0) &&
		 ph_good(x) && ph_good(y) );
    }
}}  /* namespace ellint_carlson::argcheck */


#endif /* ELLINT_ARGCHECK_HH_INCLUDED */
