#ifndef ELLINT_ARITHMETIC_HH_INCLUDED
#define ELLINT_ARITHMETIC_HH_INCLUDED


/* Algorithms for accurate floating-point summation, dot-product, and
 * polynomial evaluation.
 *
 * Ref:
 *
 * [1] S. Gaillat, V. Ménissier-Morain: Compensated Horner scheme in complex
 *      floating point arithmetic.
 *      2008-07-07, Proc. 8th Conf. Real Number Comput., pp. 133--146
 *      https://www-pequan.lip6.fr/~graillat/papers/rnc08.pdf
 *
 * [2] S. Gaillat, V. Ménissier-Morain: Accurate summation, dot product and
 *      polynomial evaluation in complex floating point arithmetic.
 *      2012-03-30, Inf. Comput., vol. 216 pp. 57--71
 *      https://doi.org/10.1016/j.ic.2011.09.003
 *      https://web.stanford.edu/group/SOL/software/qdotdd/IC2012.pdf
 *
 * [3] S. Graillat, Ph.Langlois, N.Louve: Compensated Horner Scheme.
 *      2005-07-24, Report No. RR2005-04, Université de Perpignan « Via
 *      Domitia »
 *      https://www-pequan.lip6.fr/~jmc/polycopies/Compensation-horner.pdf
 */


#include <iterator>
#include <complex>
#include <cstddef>
#include <cmath>
#include "ellint_typing.hh"
#include "ellint_arith_aux.hh"
#include "ellint_common.hh"


namespace ellint_carlson { namespace arithmetic
{
    /* The functions for arithmetic operations get their arguments by reference
     * (possibly qualified by const) whenever possible. The arguments not
     * qualified by const are meant to be overwritten in the caller. */

    /* Knuth's TwoSum algorithm, EFT on addition, with the sum written to s
     * and the correction term to p. This algorithm is generic enough to
     * write for both real and complex types (because there's no error in the
     * complex "i").
     * Algorithms 2.1 and 3.1 in Ref. [1] */
    template<typename FPT>
    inline typing::real_or_cplx<FPT, void>
    eft_sum(const FPT& x, const FPT& y, FPT& s, FPT& corr)
    {
	s = x + y;
	FPT z = s - x;
	corr = (x - (s - z)) + (y - z);
    }


    /* TwoSum in accumulator style, with sum accumulated to acc and the
     * correction term to corr
     * Algorithm 4.1 in Ref. [2] */
    template<typename FPT>
    inline typing::real_or_cplx<FPT, void>
    sum2_acc(const FPT& summand, FPT& acc, FPT& corr)
    {
	FPT tmp_s, tmp_e;
	/* The tmp variables get their values from eft_sum called-by-reference.
	 */
	eft_sum(summand, acc, tmp_s, tmp_e);
	acc = tmp_s;
	corr += tmp_e;
    }


    /* Sum at most the first n elements of iterable thing, including fixed-size
     * array passed by reference (both real and complex). */
    template<typename IterableT>
    inline auto
    nsum2(const IterableT& x, std::size_t n) -> JUST_ELEM(std::begin(x))
    {
	auto it = std::begin(x);
	typedef JUST_ELEM(it) FPT;
	static_assert(typing::is_real_or_complex_fp<FPT>::value,
		      "the sum2 function only works with real or complex "
		      "floating-point numbers.");
	FPT p(0.0), s(0.0);

	for ( std::size_t i = 0; ( it != std::end(x) && i < n ); ++it, ++i )
	{
	    sum2_acc(*it, p, s);
	}

	return p + s;
    }

    /* Convenient wrapper for summing fixed-size array. */
    template<typename FPT, std::size_t N>
    inline typing::real_or_cplx<FPT, FPT>
    sum2(const FPT(& arr)[N])
    {
	return nsum2(arr, N);
    }


    /* EFT for multiplication, with real-type arguments.
     * NOTE: This implementation uses the FMA function for correct rounding.
     * (mandated by standard since C++11)
     * Algorithm 2.5 in Ref. [1] */
    template<typename RT>
    inline typing::real_only<RT, void>
    eft_prod(const RT& x, const RT& y, RT& prod, typing::corrbuf<RT>& corr)
    {
	prod = x * y;
	corr = std::fma(x, y, -prod);
    }


    /* EFT for multiplication, with complex-type arguments.
     * Algorithm 3.4 in Ref. [1] */
    template<typename CT>
    inline typing::cplx_only<CT, void>
    eft_prod(const CT& x, const CT& y, CT& prod, typing::corrbuf<CT>& corr)
    {
	typedef typename typing::is_complex<CT>::rtype RT;
	RT z1, z2, z3, z4, z5, z6;

	/*  p                  e    f    g
	 * z5 | value_real : [h1, -h2,  h5, xx]
	 * z6 | value_imag : [h3,  h4,  h6, xx]
	 */
	RT a = x.real();
	RT b = x.imag();
	RT c = y.real();
	RT d = y.imag();

	eft_prod(a, c, z1, corr.value_real[0]);
	eft_prod(b, d, z2, corr.value_real[1]);
	eft_prod(a, d, z3, corr.value_imag[0]);
	eft_prod(b, c, z4, corr.value_imag[1]);
	eft_sum(z1, -z2, z5, corr.value_real[2]);
	eft_sum(z3, z4, z6, corr.value_imag[2]);

	corr.value_real[1] = -corr.value_real[1];
	prod = CT{z5, z6};
    }


    /* Accumulator version of dot-product
     * Algorithm 4.3 in Ref. [2] */
    template<typename RT>
    inline typing::real_only<RT, void>
    fdot2_acc(const RT& x, const RT& y, RT& acc, RT& corr)
    {
	RT h, t, q, r;

	eft_prod(x, y, h, r);
	eft_sum(acc, h, t, q);
	acc = t;
	corr += q + r;
    }


    /* Dot-product of at most first n elements, for real floating type only
     * (based on what we see by peeking "begin"). */
    template<typename Iterable>
    inline auto
    ndot2(const Iterable& x, const Iterable& y, std::size_t n)
    -> typename std::enable_if<
	std::is_floating_point<JUST_ELEM(std::begin(x))>::value,
	JUST_ELEM(std::begin(x)) >::type
    {
	auto itx = std::begin(x);
	auto ity = std::begin(y);
	typedef JUST_ELEM(itx) RT;
	RT p(0.0), s(0.0);

	for ( std::size_t i = 0;
	      ( (itx != std::end(x)) && (ity != std::end(y)) && (i < n) );
	      ++itx, ++ity, ++i )
	{
	    fdot2_acc((RT)(*itx), (RT)(*ity), p, s);
	}

	return p + s;
    }


    /* Dot-product, but for complex types.
     * Algorithm 4.4 in Ref. [2] */
    template<typename Iterable>
    inline auto
    ndot2(const Iterable& x, const Iterable& y, std::size_t n)
    -> typename std::enable_if<
	typing::is_complex<JUST_ELEM(std::begin(y))>::value,
	JUST_ELEM(std::begin(x)) >::type
    {
	auto itx = std::begin(x);
	auto ity = std::begin(y); 
	typedef JUST_ELEM(itx) CT;
	typedef typing::decplx_t<CT> RT;

	RT pr, pi, cr, ci;

	pr = pi = cr = ci = (RT)0.0;
	for ( std::size_t i = 0;
	      ( itx != std::end(x) && ity != std::end(y) && i < n );
	      ++itx, ++ity, ++i )
	{
	    RT a = itx->real();
	    RT b = ity->real();
	    RT c = itx->imag();
	    RT d = ity->imag();

	    fdot2_acc(a, b, pr, cr);
	    fdot2_acc(c, -d, pr, cr);
	    fdot2_acc(a, d, pi, ci);
	    fdot2_acc(b, c, pi, ci);
	}

	return CT{pr + cr, pi + ci};
    }

    /* Convenient wrapper for dot-product */
    template<typename FPT, std::size_t M, std::size_t N>
    inline FPT
    dot2(const FPT(& x)[M], const FPT(& y)[N])
    {
	constexpr std::size_t len = (M <= N ? M : N);
	return ndot2(x, y, len);
    }


    /* Polynomial evaluation with compensated Horner scheme.
     * This instance for the real-only inputs is essentially the Algorithm 9 in
     * Ref. [3] */
    template<typename RT>
    inline typing::real_only<RT, RT>
    dcomp_horner(const RT& x, const RT* poly, std::size_t degree)
    {
	RT s(poly[degree]), r(0.0);

	for ( std::ptrdiff_t i = (std::ptrdiff_t)degree - 1; i >= 0; --i )
	{
	    RT ptmp, pp, ps;

	    eft_prod(s, x, ptmp, pp);
	    eft_sum(ptmp, poly[i], s, ps);
	    r = r * x + (pp + ps);
	}
	return s + r;
    }

    /* For complex-valued input: Algorithm 5.4 in Ref. [2].
     * Notice that this makes use of the aux::acc_sum() algorithm, which is
     * named "Accsum" in the referenced paper. */
    template<typename T0, typename T1>
    inline typename std::enable_if< typing::is_complex<T0>::value ||
                                    typing::is_complex<T1>::value,
	                            typing::Promote<T0, T1> >::type
    dcomp_horner(const T0(& x), const T1* poly, std::size_t degree)
    {
	typedef typing::Promote<T0, T1> CT;
	CT s(poly[degree]), r(0.0);

	for ( std::ptrdiff_t i = (std::ptrdiff_t)degree - 1; i >= 0; --i )
	{
	    CT p, tmp;  /* tmp is for copying into the buffer arrays. */
	    typing::corrbuf<CT> ws;
	    bool mask_r[NPOLY];
	    bool mask_i[NPOLY];

	    eft_prod(s, CT(x), p, ws);
	    eft_sum(p, CT(poly[i]), s, tmp);
	    ws.value_real[NPOLY - 1] = tmp.real();
	    ws.value_imag[NPOLY - 1] = tmp.imag();
	    aux::make(ws.value_real, mask_r);
	    aux::make(ws.value_imag, mask_i);
	    r = r * x + CT(aux::acc_sum(ws.value_real, mask_r),
	                   aux::acc_sum(ws.value_imag, mask_i));
	}

	return s + r;
    }

    /* Convenient version for fixed-size array as polynomial coefficients, with
     * type promotion. */
    template<typename T, typename U, std::size_t N>
    inline typing::Promote<T, U>
    comp_horner(const T& x, const U(& poly)[N])
    {
	return dcomp_horner(x, poly, N - 1);
    }


    /* Small integer roots for real types. */
    template<typename RT>
    inline typing::real_only<RT, RT>
    ocrt(const RT& x)  /* 8th (octic) root */
    {
	return std::sqrt(std::sqrt(std::sqrt(x)));
    }

    template<typename RT>
    inline typing::real_only<RT, RT>
    sxrt(const RT& x)  /* 6th (sextic) root */
    {
	return std::sqrt(std::cbrt(x));
    }
}}  /* namespace ellint_carlson::arithmetic  */


#endif /* ELLINT_ARITHMETIC_HH_INCLUDED */
