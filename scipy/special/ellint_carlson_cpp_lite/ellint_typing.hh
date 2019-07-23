#ifndef ELLINT_TYPING_HH_INCLUDED
#define ELLINT_TYPING_HH_INCLUDED


/* Type-checking support implemented with templates. */


#include <type_traits>
#include <iterator>
#include <limits>
#include <complex>
#include "ellint_common.hh"


namespace ellint_carlson { namespace typing
{
    /* Constructs for confining template-matching to complex math. */
    template<typename T>
    struct is_complex : std::false_type {};

    template<typename RT>
    struct is_complex< std::complex<RT> > :
	std::enable_if< std::is_floating_point<RT>::value,
			std::true_type >::type
    {
	typedef RT rtype;
	typedef std::complex<RT> ctype;
    };

    template<typename T>
    struct is_real_or_complex_fp :
	std::integral_constant< bool,
				std::is_floating_point<T>::value ||
				is_complex<T>::value > {};

    template<typename T>
    struct decplx
    {
	typedef typename std::remove_cv<T>::type type;
    };

    template<typename T>
    struct decplx< std::complex<T> >
    {
	typedef typename std::remove_cv<T>::type type;
    };

    template<typename T>
    using decplx_t = typename decplx<T>::type;
    /* Arithmetic promotion into the most compatible (largest) floating
     * point type, but retaining the numeric domain (i.e. no over-promotion
     * to complex from the real domain) whenever possible. If all types are
     * real, the resulting type is real. Otherwise (with mixed
     * real-and-complex or complex types), the resulting type is complex.
     */
    namespace promote2impl
    {
	/* Basic building blocks. */
	/* Calculate the real, or the *underlying* real, arithmetic
	 * (floating-point) type of two-type promotion. */
	template<typename T0, typename T1>
	using RFPromote2 =
	    decltype(decplx_t<T0>() + decplx_t<T1>());

	/* Complexify the result of the real-promotion. */
	template<typename T0, typename T1>
	using CFPromote2 =
	    decltype(std::complex< RFPromote2<T0, T1> >());


	/* Compute whether we should possibly do the real -> complex domain
	 * upgrade (which should only happen when one of the operands is
	 * already a complex type). */
	template<typename T0, typename T1>
	using Promote2 =
	    typename std::conditional< is_complex<T0>::value ||
	                               is_complex<T1>::value,
			               CFPromote2<T0, T1>,
				       RFPromote2<T0, T1> >::type;
    }  /* namespace promote2impl */
    template<typename T0, typename T1>
    using Promote = promote2impl::Promote2<T0, T1>;


    /* Buffer holding the correction terms for multiplication.  For real
     * floating-point types, this is just the same type ("just" numeric type).
     * For complex types, there must be two buffers, each with length 4, for
     * holding the EFT terms. */
    namespace fpbuftypes
    {
	/* Real type (also primary template), empty struct except for the
	 * typedef referring to the numeric type. */
	template<typename RT>
	struct corrbuf
	{
	    typedef RT type;
	};

	/* Complex type, the typedef refers to the corrbuf type itself.
	 * This works because complex is more specialised than real float. */
	template<typename RT>
	struct corrbuf< std::complex<RT> >
	{
	    /* If selected, this "corrbuf" in the typedef necessarily refers to
	     * the correct type, i.e. this struct with two buffer members. */
	    typedef corrbuf type;
	    RT value_real[NPOLY];
	    RT value_imag[NPOLY];
	};
    }
    /* In the arithmetic functions, "corrbuf" refers to either the bare (and
     * real) numeric type or the specialised one for complex arithmetic, based
     * on template matching. This is a convenient type alias referring back to
     * the (either towards the bare thing or self-referential) type. */
    template<typename T>
    using corrbuf = typename fpbuftypes::corrbuf<T>::type;


    /* For use in function return type declaration, optional arguments,
     * template specialisation, inheritance, etc., where std::enable_if is
     * appropriate. */
    template<typename arg_t, typename ret_t>
    using real_only =
	typename std::enable_if<std::is_floating_point<arg_t>::value,
                                ret_t>::type;

    template<typename arg_t, typename ret_t>
    using cplx_only =
	typename std::enable_if<is_complex<arg_t>::value, ret_t>::type;

    template<typename arg_t, typename ret_t>
    using real_or_cplx =
	typename std::enable_if<is_real_or_complex_fp<arg_t>::value,
                                ret_t>::type;


    /* Constants for different sizes of floating-point number.  The struct
     * simply retrieves the information from the scope of the specialisation
     * std::numeric_limit<T>. For compatibility with c++11, these are
     * implemented as constexpr functions, because c++11 does not support
     * variable templates. */
    template<typename T>
    constexpr
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
    nan()
    {
	return std::numeric_limits<T>::quiet_NaN();
    }

    template<typename T>
    constexpr
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
    huge()
    {
	return std::numeric_limits<T>::infinity();
    }

    template<typename T>
    constexpr
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
    half_epsilon()
    {
	return std::numeric_limits<T>::epsilon() / (T)2.0;
    }

    /* Complexified specialisation, with the complex NaN and infinity as
     * return values. */
    template<typename T>
    constexpr
    typename std::enable_if<is_complex<T>::value, T>::type
    nan()
    {
	typedef decplx_t<T> RT;
	return T{nan<RT>(), nan<RT>()};
    }

    template<typename T>
    constexpr
    typename std::enable_if<is_complex<T>::value, T>::type
    huge()
    {
	typedef decplx_t<T> RT;
	return T{huge<RT>(), (RT)0.0};
    }
}}  /* namespace ellint_carlson::typing */


/* Helper macros for use elsewhere in the implementation functions. */
/* For use with iterator expressions. */
#define PIERCE_ITER_REF(ITEXPR)	\
typename std::remove_reference<decltype(*ITEXPR)>::type

#define JUST_ELEM(ITEXPR)	\
typename std::remove_cv<PIERCE_ITER_REF(ITEXPR)>::type


#endif /* ELLINT_TYPING_HH_INCLUDED */
