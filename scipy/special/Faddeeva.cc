/* Copyright (c) 2012 Massachusetts Institute of Technology
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include "Faddeeva.hh"

/* Available at: http://ab-initio.mit.edu/Faddeeva

   Computes various error functions (erf, erfc, erfi, erfcx), 
   including the Dawson integral, in the complex plane, based
   on algorithms for the computation of the Faddeeva function 
              w(z) = exp(-z^2) * erfc(-i*z).
   Given w(z), the error functions are mostly straightforward
   to compute, except for certain regions where we have to
   switch to Taylor expansions to avoid cancellation errors
   [e.g. near the origin for erf(z)].

   To compute the Faddeeva function, we use a combination of two
   algorithms:

   For sufficiently large |z|, we use a continued-fraction expansion
   for w(z) similar to those described in:

      Walter Gautschi, "Efficient computation of the complex error
      function," SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)

      G. P. M. Poppe and C. M. J. Wijers, "More efficient computation
      of the complex error function," ACM Trans. Math. Soft. 16(1),
      pp. 38-46 (1990).

   Unlike those papers, however, we switch to a completely different
   algorithm for smaller |z|:

      Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
      Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15
      (2011).

   (I initially used this algorithm for all z, but it turned out to be
    significantly slower than the continued-fraction expansion for
    larger |z|.  On the other hand, it is competitive for smaller |z|, 
    and is significantly more accurate than the Poppe & Wijers code
    in some regions, e.g. in the vicinity of z=1+1i.)

   Note that this is an INDEPENDENT RE-IMPLEMENTATION of these algorithms,
   based on the description in the papers ONLY.  In particular, I did
   not refer to the authors' Fortran or Matlab implementations, respectively,
   (which are under restrictive ACM copyright terms and therefore unusable
    in free/open-source software).

   Steven G. Johnson, Massachusetts Institute of Technology
   http://math.mit.edu/~stevenj
   October 2012.

    -- Note that Algorithm 916 assumes that the erfc(x) function, 
       or rather the scaled function erfcx(x) = exp(x*x)*erfc(x),
       is supplied for REAL arguments x.   I originally used an
       erfcx routine derived from DERFC in SLATEC, but I have
       since replaced it with a much faster routine written by
       me which uses a combination of continued-fraction expansions
       and a lookup table of Chebyshev polynomials.  For speed,
       I implemented a similar algorithm for Im[w(x)] of real x,
       since this comes up frequently in the other error functions.

   A small test program is included the end, which checks
   the w(z) etc. results against several known values.  To compile
   the test function, compile with -DTEST_FADDEEVA (that is,
   #define TEST_FADDEEVA).

   REVISION HISTORY:
       4 October 2012: Initial public release (SGJ)
       5 October 2012: Revised (SGJ) to fix spelling error,
                       start summation for large x at round(x/a) (> 1)
		       rather than ceil(x/a) as in the original
		       paper, which should slightly improve performance
     		       (and, apparently, slightly improves accuracy)
      19 October 2012: Revised (SGJ) to fix bugs for large x, large -y,
                       and 15<x<26. Performance improvements. Prototype
		       now supplies default value for relerr.
      24 October 2012: Switch to continued-fraction expansion for
                       sufficiently large z, for performance reasons.
		       Also, avoid spurious overflow for |z| > 1e154.
		       Set relerr argument to min(relerr,0.1).
      27 October 2012: Enhance accuracy in Re[w(z)] taken by itself,
                       by switching to Alg. 916 in a region near
		       the real-z axis where continued fractions
		       have poor relative accuracy in Re[w(z)].  Thanks
		       to M. Zaghloul for the tip.
      29 October 2012: Replace SLATEC-derived erfcx routine with
                       completely rewritten code by me, using a very
		       different algorithm which is much faster.
      30 October 2012: Implemented special-case code for real z
                       (where real part is exp(-x^2) and imag part is
		        Dawson integral), using algorithm similar to erfx.
		       Export ImFaddeeva_w function to make Dawson's
		       integral directly accessible.
      3 November 2012: Provide implementations of erf, erfc, erfcx,
                       and Dawson functions in Faddeeva:: namespace,
		       in addition to Faddeeva::w.  Provide header
		       file Faddeeva.hh.
*/

#include <Python.h>
extern "C" {
#include <numpy/npy_math.h>
}

#include <cfloat>
#include <cmath>

#define complex std::complex

/////////////////////////////////////////////////////////////////////////

// use numpy's versions, since std:: versions only available in C++11
#define my_isinf npy_isinf
#define my_isnan npy_isnan
#define my_copysign npy_copysign
#define Inf NPY_INFINITY
#define NaN NPY_NAN

/////////////////////////////////////////////////////////////////////////
// Auxiliary routines to compute other special functions based on w(z)

// compute erfcx(z) = exp(z^2) erfz(z)
complex<double> Faddeeva::erfcx(complex<double> z, double relerr)
{
  return Faddeeva::w(complex<double>(-imag(z), real(z)));
}

// compute the error function erf(x)
double Faddeeva::erf(double x)
{
  double mx2 = -x*x;
  if (mx2 < -750) // underflow
    return (x >= 0 ? 1.0 : -1.0);

  if (x >= 0) {
    if (x < 5e-3) goto taylor;
    return 1.0 - exp(mx2) * Faddeeva::erfcx(x);
  }
  else { // x < 0
    if (x > -5e-3) goto taylor;
    return exp(mx2) * Faddeeva::erfcx(-x) - 1.0;
  }

  // Use Taylor series for small |x|, to avoid cancellation inaccuracy
  //     erf(x) = 2/sqrt(pi) * x * (1 - x^2/3 + x^4/10 - ...)
 taylor:
  return x * (1.1283791670955125739
	      + mx2 * (0.37612638903183752464
		       + mx2 * 0.11283791670955125739));
}

// compute the error function erf(z)
complex<double> Faddeeva::erf(complex<double> z, double relerr)
{
  double x = real(z), y = imag(z);

  if (x == 0) // handle separately for speed & handling of y = Inf or NaN
    return complex<double>(x, // preserve sign of 0
			   /* handle y -> Inf limit manually, since
			      exp(y^2) -> Inf but Im[w(y)] -> 0, so
			      IEEE will give us a NaN when it should be Inf */
			   y*y > 720 ? (y > 0 ? Inf : -Inf)
			   : exp(y*y) * Faddeeva::w_im(y));

  double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
  double mIm_z2 = -2*x*y; // Im(-z^2)
  if (mRe_z2 < -750) // underflow
    return (x >= 0 ? 1.0 : -1.0);

  /* Handle positive and negative x via different formulas,
     using the mirror symmetries of w, to avoid overflow/underflow
     problems from multiplying exponentially large and small quantities. */
  if (x >= 0) {
    if (x < 5e-3) {
      if (fabs(y) < 5e-3)
	goto taylor;
      else if (fabs(mIm_z2) < 5e-3)
	goto taylor_erfi;
    }
    /* don't use complex exp function, since that will produce spurious NaN
       values when multiplying w in an overflow situation. */
    return 1.0 - exp(mRe_z2) *
      (complex<double>(cos(mIm_z2), sin(mIm_z2))
       * Faddeeva::w(complex<double>(-y,x)));
  }
  else { // x < 0
    if (x > -5e-3) { // duplicate from above to avoid fabs(x) call
      if (fabs(y) < 5e-3)
	goto taylor;
      else if (fabs(mIm_z2) < 5e-3)
	goto taylor_erfi;
    }
    else if (my_isnan(x))
      return complex<double>(NaN, y == 0 ? 0 : NaN);
    /* don't use complex exp function, since that will produce spurious NaN
       values when multiplying w in an overflow situation. */
    return exp(mRe_z2) *
      (complex<double>(cos(mIm_z2), sin(mIm_z2))
       * Faddeeva::w(complex<double>(y,-x))) - 1.0;
  }

  // Use Taylor series for small |z|, to avoid cancellation inaccuracy
  //     erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - ...)
 taylor:
  {
    complex<double> mz2(mRe_z2, mIm_z2); // -z^2
    return z * (1.1283791670955125739
		+ mz2 * (0.37612638903183752464
			 + mz2 * 0.11283791670955125739));
  }

  /* for small |x| and small |xy|, 
     use Taylor series to avoid cancellation inaccuracy:
       erf(x+iy) = erf(iy)
          + 2*exp(y^2)/sqrt(pi) *
	    [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ... 
              - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
     where:
        erf(iy) = exp(y^2) * Im[w(y)]
  */
 taylor_erfi:
  {
    double x2 = x*x, y2 = y*y;
    double expy2 = exp(y2);
    return complex<double>
      (expy2 * x * (1.1283791670955125739
		    - x2 * (0.37612638903183752464
			    + 0.75225277806367504925*y2)
		    + x2*x2 * (0.11283791670955125739
			       + y2 * (0.45135166683820502956
				       + 0.15045055561273500986*y2))),
       expy2 * (Faddeeva::w_im(y)
		- x2*y * (1.1283791670955125739 
			  - x2 * (0.56418958354775628695 
				  + 0.37612638903183752464*y2))));
  }
}

// erfi(z) = -i erf(iz)
complex<double> Faddeeva::erfi(complex<double> z, double relerr)
{
  complex<double> e = Faddeeva::erf(complex<double>(-imag(z),real(z)), relerr);
  return complex<double>(imag(e), -real(e));
}

// erfi(x) = -i erf(ix)
double Faddeeva::erfi(double x)
{
  return x*x > 720 ? (x > 0 ? Inf : -Inf)
    : exp(x*x) * Faddeeva::w_im(x);
}

// erfc(x) = 1 - erf(x)
double Faddeeva::erfc(double x)
{
  if (x*x > 750) // underflow
    return (x >= 0 ? 0.0 : 2.0);
  return x >= 0 ? exp(-x*x) * Faddeeva::erfcx(x) 
    : 2. - exp(-x*x) * Faddeeva::erfcx(-x);
}

// erfc(z) = 1 - erf(z)
complex<double> Faddeeva::erfc(complex<double> z, double relerr)
{
  double x = real(z), y = imag(z);

  if (x == 0.)
    return complex<double>(1,
			   /* handle y -> Inf limit manually, since
			      exp(y^2) -> Inf but Im[w(y)] -> 0, so
			      IEEE will give us a NaN when it should be Inf */
			   y*y > 720 ? (y > 0 ? -Inf : Inf)
			   : -exp(y*y) * Faddeeva::w_im(y));
  if (y == 0.) {
    if (x*x > 750) // underflow
      return (x >= 0 ? 0.0 : 2.0);
    return x >= 0 ? exp(-x*x) * Faddeeva::erfcx(x) 
      : 2. - exp(-x*x) * Faddeeva::erfcx(-x);
  }

  double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
  double mIm_z2 = -2*x*y; // Im(-z^2)
  if (mRe_z2 < -750) // underflow
    return (x >= 0 ? 0.0 : 2.0);

  if (x >= 0)
    return exp(complex<double>(mRe_z2, mIm_z2))
      * Faddeeva::w(complex<double>(-y,x), relerr);
  else
    return 2.0 - exp(complex<double>(mRe_z2, mIm_z2))
      * Faddeeva::w(complex<double>(y,-x), relerr);
}

// compute Dawson(x) = sqrt(pi)/2  *  exp(-x^2) * erfi(x)
double Faddeeva::Dawson(double x)
{
  const double spi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
  return spi2 * Faddeeva::w_im(x);
}

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
complex<double> Faddeeva::Dawson(complex<double> z, double relerr)
{
  const double spi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
  double x = real(z), y = imag(z);

  // handle axes separately for speed & proper handling of x or y = Inf or NaN
  if (y == 0)
    return complex<double>(spi2 * Faddeeva::w_im(x),
			   -y); // preserve sign of 0
  if (x == 0) {
    double y2 = y*y;
    if (y2 < 2.5e-5) { // Taylor expansion
      return complex<double>(x, // preserve sign of 0
	 y * (1.
	      + y2 * (0.6666666666666666666666666666666666666667
		      + y2 * 0.2666666666666666666666666666666666666667)));
    }
    return complex<double>(x, // preserve sign of 0
			   spi2 * (y >= 0 
				   ? exp(y2) - Faddeeva::erfcx(y)
				   : Faddeeva::erfcx(-y) - exp(y2)));
  }

  double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
  double mIm_z2 = -2*x*y; // Im(-z^2)
  complex<double> mz2(mRe_z2, mIm_z2); // -z^2

  /* Handle positive and negative x via different formulas,
     using the mirror symmetries of w, to avoid overflow/underflow
     problems from multiplying exponentially large and small quantities. */
  if (y >= 0) {
    if (y < 5e-3) {
      if (fabs(x) < 5e-3)
	goto taylor;
      else if (fabs(mIm_z2) < 5e-3)
	goto taylor_realaxis;
    }
    complex<double> res = exp(mz2) - Faddeeva::w(z);
    return spi2 * complex<double>(-imag(res), real(res));
  }
  else { // y < 0
    if (y > -5e-3) { // duplicate from above to avoid fabs(x) call
      if (fabs(x) < 5e-3)
	goto taylor;
      else if (fabs(mIm_z2) < 5e-3)
	goto taylor_realaxis;
    }
    else if (my_isnan(y))
      return complex<double>(x == 0 ? 0 : NaN, NaN);
    complex<double> res = Faddeeva::w(-z) - exp(mz2);
    return spi2 * complex<double>(-imag(res), real(res));
  }

  // Use Taylor series for small |z|, to avoid cancellation inaccuracy
  //     dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
 taylor:
  return z * (1.
	      + mz2 * (0.6666666666666666666666666666666666666667
		       + mz2 * 0.2666666666666666666666666666666666666667));

  /* for small |y| and small |xy|, 
     use Taylor series to avoid cancellation inaccuracy:
       dawson(x + iy)
        = D + y^2 (D + x - 2Dx^2)
            + y^4 (D/2 + 5x/6 - 2Dx^2 - x^3/3 + 2Dx^4/3)
        + iy [ (1-2Dx) + 2/3 y^2 (1 - 3Dx - x^2 + 2Dx^3)
              + y^4/15 (4 - 15Dx - 9x^2 + 20Dx^3 + 2x^4 - 4Dx^5) ] + ...
     where D = dawson(x) 

     However, for large |x|, 2Dx -> 1 which gives cancellation problems in
     this series (many of the leading terms cancel).  So, for large |x|,
     we need to substitute a continued-fraction expansion for D.

        dawson(x) = 0.5 / (x-0.5/(x-1/(x-1.5/(x-2/(x-2.5/(x...))))))

     The 6 terms shown here seems to be the minimum needed to be
     accurate as soon as the simpler Taylor expansion above starts
     breaking down.  Using this 6-term expansion, factoring out the
     denominator, and simplifying with Maple, we obtain:

      Re dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / x
        = 33 - 28x^2 + 4x^4 + y^2 (18 - 4x^2) + 4 y^4
      Im dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / y
        = -15 + 24x^2 - 4x^4 + 2/3 y^2 (6x^2 - 15) - 4 y^4

     Finally, for |x| > 5e7, we can use a simpler 1-term continued-fraction
     expansion for the real part, and a 2-term expansion for the imaginary
     part.  (This avoids overflow problems for huge |x|.)  This yields:
     
     Re dawson(x + iy) = [1 + y^2 (1 + y^2/2 - (xy)^2/3)] / (2x)
     Im dawson(x + iy) = y [ -1 - 2/3 y^2 + y^4/15 (2x^2 - 4) ] / (2x^2 - 1)

 */
 taylor_realaxis:
  {
    double x2 = x*x;
    if (x2 > 1600) { // |x| > 40
      double y2 = y*y;
      if (x2 > 25e14) {// |x| > 5e7
	double xy2 = (x*y)*(x*y);
	return complex<double>((0.5 + y2 * (0.5 + 0.25*y2
					    - 0.16666666666666666667*xy2)) / x,
			       y * (-1 + y2 * (-0.66666666666666666667
					       + 0.13333333333333333333*xy2
					       - 0.26666666666666666667*y2))
			       / (2*x2 - 1));
      }
      return (1. / (-15 + x2*(90 + x2*(-60 + 8*x2)))) *
	complex<double>(x * (33 + x2 * (-28 + 4*x2)
			     + y2 * (18 - 4*x2 + 4*y2)),
			y * (-15 + x2 * (24 - 4*x2)
			     + y2 * (4*x2 - 10 - 4*y2)));
    }
    else {
      double D = spi2 * Faddeeva::w_im(x);
      double x2 = x*x, y2 = y*y;
      return complex<double>
	(D + y2 * (D + x - 2*D*x2)
	 + y2*y2 * (D * (0.5 - x2 * (2 - 0.66666666666666666667*x2))
		    + x * (0.83333333333333333333
			   - 0.33333333333333333333 * x2)),
	 y * (1 - 2*D*x
	      + y2 * 0.66666666666666666667 * (1 - x2 - D*x * (3 - 2*x2))
	      + y2*y2 * (0.26666666666666666667 -
			 x2 * (0.6 - 0.13333333333333333333 * x2)
			 - D*x * (1 - x2 * (1.3333333333333333333
					    - 0.26666666666666666667 * x2)))));
    }
  }
}

/////////////////////////////////////////////////////////////////////////

// return sinc(x) = sin(x)/x, given both x and sin(x) 
// [since we only use this in cases where sin(x) has already been computed]
static inline double sinc(double x, double sinx) { 
  return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667)*x*x : sinx / x; 
}

// sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
static inline double sinh_taylor(double x) {
  return x * (1 + (x*x) * (0.1666666666666666666667
			   + 0.00833333333333333333333 * (x*x)));
}

static inline double sqr(double x) { return x*x; }

// precomputed table of expa2n2[n-1] = exp(-a2*n*n)
// for double-precision a2 = 0.26865... in Faddeeva::w, below.
static const double expa2n2[] = {
  7.64405281671221563e-01,
  3.41424527166548425e-01,
  8.91072646929412548e-02,
  1.35887299055460086e-02,
  1.21085455253437481e-03,
  6.30452613933449404e-05,
  1.91805156577114683e-06,
  3.40969447714832381e-08,
  3.54175089099469393e-10,
  2.14965079583260682e-12,
  7.62368911833724354e-15,
  1.57982797110681093e-17,
  1.91294189103582677e-20,
  1.35344656764205340e-23,
  5.59535712428588720e-27,
  1.35164257972401769e-30,
  1.90784582843501167e-34,
  1.57351920291442930e-38,
  7.58312432328032845e-43,
  2.13536275438697082e-47,
  3.51352063787195769e-52,
  3.37800830266396920e-57,
  1.89769439468301000e-62,
  6.22929926072668851e-68,
  1.19481172006938722e-73,
  1.33908181133005953e-79,
  8.76924303483223939e-86,
  3.35555576166254986e-92,
  7.50264110688173024e-99,
  9.80192200745410268e-106,
  7.48265412822268959e-113,
  3.33770122566809425e-120,
  8.69934598159861140e-128,
  1.32486951484088852e-135,
  1.17898144201315253e-143,
  6.13039120236180012e-152,
  1.86258785950822098e-160,
  3.30668408201432783e-169,
  3.43017280887946235e-178,
  2.07915397775808219e-187,
  7.36384545323984966e-197,
  1.52394760394085741e-206,
  1.84281935046532100e-216,
  1.30209553802992923e-226,
  5.37588903521080531e-237,
  1.29689584599763145e-247,
  1.82813078022866562e-258,
  1.50576355348684241e-269,
  7.24692320799294194e-281,
  2.03797051314726829e-292,
  3.34880215927873807e-304,
  0.0 // underflow (also prevents reads past array end, below)
};

/////////////////////////////////////////////////////////////////////////

complex<double> Faddeeva::w(complex<double> z, double relerr)
{
  if (real(z) == 0.0)
    return complex<double>(Faddeeva::erfcx(imag(z)), 
			   real(z)); // give correct sign of 0 in imag(w)
  else if (imag(z) == 0)
    return complex<double>(exp(-sqr(real(z))),
			   Faddeeva::w_im(real(z)));

  double a, a2, c;
  if (relerr <= DBL_EPSILON) {
    relerr = DBL_EPSILON;
    a = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
    c = 0.329973702884629072537; // (2/pi) * a;
    a2 = 0.268657157075235951582; // a^2
  }
  else {
    const double pi = 3.14159265358979323846264338327950288419716939937510582;
    if (relerr > 0.1) relerr = 0.1; // not sensible to compute < 1 digit
    a = pi / sqrt(-log(relerr*0.5));
    c = (2/pi)*a;
    a2 = a*a;
  }
  const double x = fabs(real(z));
  const double y = imag(z), ya = fabs(y);

  complex<double> ret(0.,0.); // return value

  double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

#define USE_CONTINUED_FRACTION 1 // 1 to use continued fraction for large |z|

#if USE_CONTINUED_FRACTION
  if (ya > 7 || (x > 6  // continued fraction is faster
		 /* As pointed out by M. Zaghloul, the continued
		    fraction seems to give a large relative error in
		    Re w(z) for |x| ~ 6 and small |y|, so use
		    algorithm 816 in this region: */
		 && (ya > 0.1 || (x > 8 && ya > 1e-10) || x > 28))) {
    
    /* Poppe & Wijers suggest using a number of terms
           nu = 3 + 1442 / (26*rho + 77)
       where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
       (They only use this expansion for rho >= 1, but rho a little less
        than 1 seems okay too.)
       Instead, I did my own fit to a slightly different function
       that avoids the hypotenuse calculation, using NLopt to minimize
       the sum of the squares of the errors in nu with the constraint
       that the estimated nu be >= minimum nu to attain machine precision.
       I also separate the regions where nu == 2 and nu == 1. */
    const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
    double xs = y < 0 ? -real(z) : real(z); // compute for -z if y < 0
    if (x + ya > 4000) { // nu <= 2
      if (x + ya > 1e7) { // nu == 1, w(z) = i/sqrt(pi) / z
	// scale to avoid overflow
	if (x > ya) {
	  double yax = ya / xs; 
	  double denom = ispi / (xs + yax*ya);
	  ret = complex<double>(denom*yax, denom);
	}
	else if (my_isinf(ya))
	  return ((my_isnan(x) || y < 0) 
		  ? complex<double>(NaN,NaN) : complex<double>(0,0));
	else {
	  double xya = xs / ya;
	  double denom = ispi / (xya*xs + ya);
	  ret = complex<double>(denom, denom*xya);
	}
      }
      else { // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
	double dr = xs*xs - ya*ya - 0.5, di = 2*xs*ya;
	double denom = ispi / (dr*dr + di*di);
	ret = complex<double>(denom * (xs*di-ya*dr), denom * (xs*dr+ya*di));
      }
    }
    else { // compute nu(z) estimate and do general continued fraction
      const double c0=3.9, c1=11.398, c2=0.08254, c3=0.1421, c4=0.2023; // fit
      double nu = floor(c0 + c1 / (c2*x + c3*ya + c4));
      double wr = xs, wi = ya;
      for (nu = 0.5 * (nu - 1); nu > 0.4; nu -= 0.5) {
	// w <- z - nu/w:
	double denom = nu / (wr*wr + wi*wi);
	wr = xs - wr * denom;
	wi = ya + wi * denom;
      }
      { // w(z) = i/sqrt(pi) / w:
	double denom = ispi / (wr*wr + wi*wi);
	ret = complex<double>(denom*wi, denom*wr);
      }
    }
    if (y < 0) {
      // use w(z) = 2.0*exp(-z*z) - w(-z), 
      // but be careful of overflow in exp(-z*z) 
      //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
      return 2.0*exp(complex<double>((ya-xs)*(xs+ya), 2*xs*y)) - ret;
    }
    else
      return ret;
  }
#else // !USE_CONTINUED_FRACTION
  if (x + ya > 1e7) { // w(z) = i/sqrt(pi) / z, to machine precision
    const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
    double xs = y < 0 ? -real(z) : real(z); // compute for -z if y < 0
    // scale to avoid overflow
    if (x > ya) {
      double yax = ya / xs; 
      double denom = ispi / (xs + yax*ya);
      ret = complex<double>(denom*yax, denom);
    }
    else {
      double xya = xs / ya;
      double denom = ispi / (xya*xs + ya);
      ret = complex<double>(denom, denom*xya);
    }
    if (y < 0) {
      // use w(z) = 2.0*exp(-z*z) - w(-z), 
      // but be careful of overflow in exp(-z*z) 
      //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
      return 2.0*exp(complex<double>((ya-xs)*(xs+ya), 2*xs*y)) - ret;
    }
    else
      return ret;
  }
#endif // !USE_CONTINUED_FRACTION 

  /* Note: The test that seems to be suggested in the paper is x <
     sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
     underflows to zero and sum1,sum2,sum4 are zero.  However, long
     before this occurs, the sum1,sum2,sum4 contributions are
     negligible in double precision; I find that this happens for x >
     about 6, for all y.  On the other hand, I find that the case
     where we compute all of the sums is faster (at least with the
     precomputed expa2n2 table) until about x=10.  Furthermore, if we
     try to compute all of the sums for x > 20, I find that we
     sometimes run into numerical problems because underflow/overflow
     problems start to appear in the various coefficients of the sums,
     below.  Therefore, we use x < 10 here. */
  else if (x < 10) {
    double prod2ax = 1, prodm2ax = 1;
    double expx2;

    if (my_isnan(y))
      return complex<double>(y,y);
    
    /* Somewhat ugly copy-and-paste duplication here, but I see significant
       speedups from using the special-case code with the precomputed
       exponential, and the x < 5e-4 special case is needed for accuracy. */

    if (relerr == DBL_EPSILON) { // use precomputed exp(-a2*(n*n)) table
      if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
	const double x2 = x*x;
	expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
        // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
	const double ax2 = 1.036642960860171859744*x; // 2*a*x
	const double exp2ax =
	  1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667*ax2));
	const double expm2ax =
	  1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667*ax2));
	for (int n = 1; 1; ++n) {
	  const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
	  prod2ax *= exp2ax;
	  prodm2ax *= expm2ax;
	  sum1 += coef;
	  sum2 += coef * prodm2ax;
	  sum3 += coef * prod2ax;
	  
	  // really = sum5 - sum4
	  sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
	  
	  // test convergence via sum3
	  if (coef * prod2ax < relerr * sum3) break;
	}
      }
      else { // x > 5e-4, compute sum4 and sum5 separately
	expx2 = exp(-x*x);
	const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
	for (int n = 1; 1; ++n) {
	  const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
	  prod2ax *= exp2ax;
	  prodm2ax *= expm2ax;
	  sum1 += coef;
	  sum2 += coef * prodm2ax;
	  sum4 += (coef * prodm2ax) * (a*n);
	  sum3 += coef * prod2ax;
	  sum5 += (coef * prod2ax) * (a*n);
	  // test convergence via sum5, since this sum has the slowest decay
	  if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
	}
      }
    }
    else { // relerr != DBL_EPSILON, compute exp(-a2*(n*n)) on the fly
      const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
      if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
	const double x2 = x*x;
	expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
	for (int n = 1; 1; ++n) {
	  const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
	  prod2ax *= exp2ax;
	  prodm2ax *= expm2ax;
	  sum1 += coef;
	  sum2 += coef * prodm2ax;
	  sum3 += coef * prod2ax;
	  
	  // really = sum5 - sum4
	  sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
	  
	  // test convergence via sum3
	  if (coef * prod2ax < relerr * sum3) break;
	}
      }
      else { // x > 5e-4, compute sum4 and sum5 separately
	expx2 = exp(-x*x);
	for (int n = 1; 1; ++n) {
	  const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
	  prod2ax *= exp2ax;
	  prodm2ax *= expm2ax;
	  sum1 += coef;
	  sum2 += coef * prodm2ax;
	  sum4 += (coef * prodm2ax) * (a*n);
	  sum3 += coef * prod2ax;
	  sum5 += (coef * prod2ax) * (a*n);
	  // test convergence via sum5, since this sum has the slowest decay
	  if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
	}
      }
    }
    const double expx2erfcxy = // avoid spurious overflow for large negative y
      y > -6 // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
      ? expx2*Faddeeva::erfcx(y) : 2*exp(y*y-x*x);
    if (y > 5) { // imaginary terms cancel
      const double sinxy = sin(x*y);
      ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y)
	+ (c*x*expx2) * sinxy * sinc(x*y, sinxy);
    }
    else {
      double xs = real(z);
      const double sinxy = sin(xs*y);
      const double sin2xy = sin(2*xs*y), cos2xy = cos(2*xs*y);
      const double coef1 = expx2erfcxy - c*y*sum1;
      const double coef2 = c*xs*expx2;
      ret = complex<double>(coef1 * cos2xy + coef2 * sinxy * sinc(xs*y, sinxy),
			    coef2 * sinc(2*xs*y, sin2xy) - coef1 * sin2xy);
    }
  }
  else { // x large: only sum3 & sum5 contribute (see above note)    
    if (my_isnan(x))
      return complex<double>(x,x);
    if (my_isnan(y))
      return complex<double>(y,y);

#if USE_CONTINUED_FRACTION
    ret = exp(-x*x); // |y| < 1e-10, so we only need exp(-x*x) term
#else
    if (y < 0) {
      /* erfcx(y) ~ 2*exp(y*y) + (< 1) if y < 0, so
	 erfcx(y)*exp(-x*x) ~ 2*exp(y*y-x*x) term may not be negligible
	 if y*y - x*x > -36 or so.  So, compute this term just in case.
	 We also need the -exp(-x*x) term to compute Re[w] accurately
	 in the case where y is very small. */
      ret = polar(2*exp(y*y-x*x) - exp(-x*x), -2*real(z)*y);
    }
    else
      ret = exp(-x*x); // not negligible in real part if y very small
#endif
    // (round instead of ceil as in original paper; note that x/a > 1 here)
    double n0 = floor(x/a + 0.5); // sum in both directions, starting at n0
    double dx = a*n0 - x;
    sum3 = exp(-dx*dx) / (a2*(n0*n0) + y*y);
    sum5 = a*n0 * sum3;
    double exp1 = exp(4*a*dx), exp1dn = 1;
    int dn;
    for (dn = 1; n0 - dn > 0; ++dn) { // loop over n0-dn and n0+dn terms
      double np = n0 + dn, nm = n0 - dn;
      double tp = exp(-sqr(a*dn+dx));
      double tm = tp * (exp1dn *= exp1); // trick to get tm from tp
      tp /= (a2*(np*np) + y*y);
      tm /= (a2*(nm*nm) + y*y);
      sum3 += tp + tm;
      sum5 += a * (np * tp + nm * tm);
      if (a * (np * tp + nm * tm) < relerr * sum5) goto finish;
    }
    while (1) { // loop over n0+dn terms only (since n0-dn <= 0)
      double np = n0 + dn++;
      double tp = exp(-sqr(a*dn+dx)) / (a2*(np*np) + y*y);
      sum3 += tp;
      sum5 += a * np * tp;
      if (a * np * tp < relerr * sum5) goto finish;
    }
  }
 finish:
  return ret + complex<double>((0.5*c)*y*(sum2+sum3), 
			       (0.5*c)*my_copysign(sum5-sum4, real(z)));
}

/////////////////////////////////////////////////////////////////////////

/* erfcx(x) = exp(x^2) erfc(x) function, for real x, written by
   Steven G. Johnson, October 2012.

   This function combines a few different ideas.

   First, for x > 50, it uses a continued-fraction expansion (same as
   for the Faddeeva function, but with algebraic simplifications for z=i*x).

   Second, for 0 <= x <= 50, it uses Chebyshev polynomial approximations,
   but with two twists:

      a) It maps x to y = 4 / (4+x) in [0,1].  This simple transformation,
         inspired by a similar transformation in the octave-forge/specfun
	 erfcx by Soren Hauberg, results in much faster Chebyshev convergence
	 than other simple transformations I have examined.

      b) Instead of using a single Chebyshev polynomial for the entire
         [0,1] y interval, we break the interval up into 100 equal
	 subintervals, with a switch/lookup table, and use much lower
	 degree Chebyshev polynomials in each subinterval. This greatly
	 improves performance in my tests.

   For x < 0, we use the relationship erfcx(-x) = 2 exp(x^2) - erfc(x),
   with the usual checks for overflow etcetera.

   Performance-wise, it seems to be substantially faster than either
   the SLATEC DERFC function [or an erfcx function derived therefrom]
   or Cody's CALERF function (from netlib.org/specfun), while
   retaining near machine precision in accuracy.  */

/* Given y100=100*y, where y = 4/(4+x) for x >= 0, compute erfc(x).

   Uses a look-up table of 100 different Chebyshev polynomials
   for y intervals [0,0.01], [0.01,0.02], ...., [0.99,1], generated
   with the help of Maple and a little shell script.   This allows
   the Chebyshev polynomials to be of significantly lower degree (about 1/4)
   compared to fitting the whole [0,1] interval with a single polynomial. */
static double erfcx_y100(double y100)
{
  switch ((int) y100) {
case 0: {
double t = 2*y100 - 1;
return 0.70878032454106438663e-3 + (0.71234091047026302958e-3 + (0.35779077297597742384e-5 + (0.17403143962587937815e-7 + (0.81710660047307788845e-10 + (0.36885022360434957634e-12 + 0.15917038551111111111e-14 * t) * t) * t) * t) * t) * t;
}
case 1: {
double t = 2*y100 - 3;
return 0.21479143208285144230e-2 + (0.72686402367379996033e-3 + (0.36843175430938995552e-5 + (0.18071841272149201685e-7 + (0.85496449296040325555e-10 + (0.38852037518534291510e-12 + 0.16868473576888888889e-14 * t) * t) * t) * t) * t) * t;
}
case 2: {
double t = 2*y100 - 5;
return 0.36165255935630175090e-2 + (0.74182092323555510862e-3 + (0.37948319957528242260e-5 + (0.18771627021793087350e-7 + (0.89484715122415089123e-10 + (0.40935858517772440862e-12 + 0.17872061464888888889e-14 * t) * t) * t) * t) * t) * t;
}
case 3: {
double t = 2*y100 - 7;
return 0.51154983860031979264e-2 + (0.75722840734791660540e-3 + (0.39096425726735703941e-5 + (0.19504168704300468210e-7 + (0.93687503063178993915e-10 + (0.43143925959079664747e-12 + 0.18939926435555555556e-14 * t) * t) * t) * t) * t) * t;
}
case 4: {
double t = 2*y100 - 9;
return 0.66457513172673049824e-2 + (0.77310406054447454920e-3 + (0.40289510589399439385e-5 + (0.20271233238288381092e-7 + (0.98117631321709100264e-10 + (0.45484207406017752971e-12 + 0.20076352213333333333e-14 * t) * t) * t) * t) * t) * t;
}
case 5: {
double t = 2*y100 - 11;
return 0.82082389970241207883e-2 + (0.78946629611881710721e-3 + (0.41529701552622656574e-5 + (0.21074693344544655714e-7 + (0.10278874108587317989e-9 + (0.47965201390613339638e-12 + 0.21285907413333333333e-14 * t) * t) * t) * t) * t) * t;
}
case 6: {
double t = 2*y100 - 13;
return 0.98039537275352193165e-2 + (0.80633440108342840956e-3 + (0.42819241329736982942e-5 + (0.21916534346907168612e-7 + (0.10771535136565470914e-9 + (0.50595972623692822410e-12 + 0.22573462684444444444e-14 * t) * t) * t) * t) * t) * t;
}
case 7: {
double t = 2*y100 - 15;
return 0.11433927298290302370e-1 + (0.82372858383196561209e-3 + (0.44160495311765438816e-5 + (0.22798861426211986056e-7 + (0.11291291745879239736e-9 + (0.53386189365816880454e-12 + 0.23944209546666666667e-14 * t) * t) * t) * t) * t) * t;
}
case 8: {
double t = 2*y100 - 17;
return 0.13099232878814653979e-1 + (0.84167002467906968214e-3 + (0.45555958988457506002e-5 + (0.23723907357214175198e-7 + (0.11839789326602695603e-9 + (0.56346163067550237877e-12 + 0.25403679644444444444e-14 * t) * t) * t) * t) * t) * t;
}
case 9: {
double t = 2*y100 - 19;
return 0.14800987015587535621e-1 + (0.86018092946345943214e-3 + (0.47008265848816866105e-5 + (0.24694040760197315333e-7 + (0.12418779768752299093e-9 + (0.59486890370320261949e-12 + 0.26957764568888888889e-14 * t) * t) * t) * t) * t) * t;
}
case 10: {
double t = 2*y100 - 21;
return 0.16540351739394069380e-1 + (0.87928458641241463952e-3 + (0.48520195793001753903e-5 + (0.25711774900881709176e-7 + (0.13030128534230822419e-9 + (0.62820097586874779402e-12 + 0.28612737351111111111e-14 * t) * t) * t) * t) * t) * t;
}
case 11: {
double t = 2*y100 - 23;
return 0.18318536789842392647e-1 + (0.89900542647891721692e-3 + (0.50094684089553365810e-5 + (0.26779777074218070482e-7 + (0.13675822186304615566e-9 + (0.66358287745352705725e-12 + 0.30375273884444444444e-14 * t) * t) * t) * t) * t) * t;
}
case 12: {
double t = 2*y100 - 25;
return 0.20136801964214276775e-1 + (0.91936908737673676012e-3 + (0.51734830914104276820e-5 + (0.27900878609710432673e-7 + (0.14357976402809042257e-9 + (0.70114790311043728387e-12 + 0.32252476000000000000e-14 * t) * t) * t) * t) * t) * t;
}
case 13: {
double t = 2*y100 - 27;
return 0.21996459598282740954e-1 + (0.94040248155366777784e-3 + (0.53443911508041164739e-5 + (0.29078085538049374673e-7 + (0.15078844500329731137e-9 + (0.74103813647499204269e-12 + 0.34251892320000000000e-14 * t) * t) * t) * t) * t) * t;
}
case 14: {
double t = 2*y100 - 29;
return 0.23898877187226319502e-1 + (0.96213386835900177540e-3 + (0.55225386998049012752e-5 + (0.30314589961047687059e-7 + (0.15840826497296335264e-9 + (0.78340500472414454395e-12 + 0.36381553564444444445e-14 * t) * t) * t) * t) * t) * t;
}
case 15: {
double t = 2*y100 - 31;
return 0.25845480155298518485e-1 + (0.98459293067820123389e-3 + (0.57082915920051843672e-5 + (0.31613782169164830118e-7 + (0.16646478745529630813e-9 + (0.82840985928785407942e-12 + 0.38649975768888888890e-14 * t) * t) * t) * t) * t) * t;
}
case 16: {
double t = 2*y100 - 33;
return 0.27837754783474696598e-1 + (0.10078108563256892757e-2 + (0.59020366493792212221e-5 + (0.32979263553246520417e-7 + (0.17498524159268458073e-9 + (0.87622459124842525110e-12 + 0.41066206488888888890e-14 * t) * t) * t) * t) * t) * t;
}
case 17: {
double t = 2*y100 - 35;
return 0.29877251304899307550e-1 + (0.10318204245057349310e-2 + (0.61041829697162055093e-5 + (0.34414860359542720579e-7 + (0.18399863072934089607e-9 + (0.92703227366365046533e-12 + 0.43639844053333333334e-14 * t) * t) * t) * t) * t) * t;
}
case 18: {
double t = 2*y100 - 37;
return 0.31965587178596443475e-1 + (0.10566560976716574401e-2 + (0.63151633192414586770e-5 + (0.35924638339521924242e-7 + (0.19353584758781174038e-9 + (0.98102783859889264382e-12 + 0.46381060817777777779e-14 * t) * t) * t) * t) * t) * t;
}
case 19: {
double t = 2*y100 - 39;
return 0.34104450552588334840e-1 + (0.10823541191350532574e-2 + (0.65354356159553934436e-5 + (0.37512918348533521149e-7 + (0.20362979635817883229e-9 + (0.10384187833037282363e-11 + 0.49300625262222222221e-14 * t) * t) * t) * t) * t) * t;
}
case 20: {
double t = 2*y100 - 41;
return 0.36295603928292425716e-1 + (0.11089526167995268200e-2 + (0.67654845095518363577e-5 + (0.39184292949913591646e-7 + (0.21431552202133775150e-9 + (0.10994259106646731797e-11 + 0.52409949102222222221e-14 * t) * t) * t) * t) * t) * t;
}
case 21: {
double t = 2*y100 - 43;
return 0.38540888038840509795e-1 + (0.11364917134175420009e-2 + (0.70058230641246312003e-5 + (0.40943644083718586939e-7 + (0.22563034723692881631e-9 + (0.11642841011361992885e-11 + 0.55721092871111111110e-14 * t) * t) * t) * t) * t) * t;
}
case 22: {
double t = 2*y100 - 45;
return 0.40842225954785960651e-1 + (0.11650136437945673891e-2 + (0.72569945502343006619e-5 + (0.42796161861855042273e-7 + (0.23761401711005024162e-9 + (0.12332431172381557035e-11 + 0.59246802364444444445e-14 * t) * t) * t) * t) * t) * t;
}
case 23: {
double t = 2*y100 - 47;
return 0.43201627431540222422e-1 + (0.11945628793917272199e-2 + (0.75195743532849206263e-5 + (0.44747364553960993492e-7 + (0.25030885216472953674e-9 + (0.13065684400300476484e-11 + 0.63000532853333333334e-14 * t) * t) * t) * t) * t) * t;
}
case 24: {
double t = 2*y100 - 49;
return 0.45621193513810471438e-1 + (0.12251862608067529503e-2 + (0.77941720055551920319e-5 + (0.46803119830954460212e-7 + (0.26375990983978426273e-9 + (0.13845421370977119765e-11 + 0.66996477404444444445e-14 * t) * t) * t) * t) * t) * t;
}
case 25: {
double t = 2*y100 - 51;
return 0.48103121413299865517e-1 + (0.12569331386432195113e-2 + (0.80814333496367673980e-5 + (0.48969667335682018324e-7 + (0.27801515481905748484e-9 + (0.14674637611609884208e-11 + 0.71249589351111111110e-14 * t) * t) * t) * t) * t) * t;
}
case 26: {
double t = 2*y100 - 53;
return 0.50649709676983338501e-1 + (0.12898555233099055810e-2 + (0.83820428414568799654e-5 + (0.51253642652551838659e-7 + (0.29312563849675507232e-9 + (0.15556512782814827846e-11 + 0.75775607822222222221e-14 * t) * t) * t) * t) * t) * t;
}
case 27: {
double t = 2*y100 - 55;
return 0.53263363664388864181e-1 + (0.13240082443256975769e-2 + (0.86967260015007658418e-5 + (0.53662102750396795566e-7 + (0.30914568786634796807e-9 + (0.16494420240828493176e-11 + 0.80591079644444444445e-14 * t) * t) * t) * t) * t) * t;
}
case 28: {
double t = 2*y100 - 57;
return 0.55946601353500013794e-1 + (0.13594491197408190706e-2 + (0.90262520233016380987e-5 + (0.56202552975056695376e-7 + (0.32613310410503135996e-9 + (0.17491936862246367398e-11 + 0.85713381688888888890e-14 * t) * t) * t) * t) * t) * t;
}
case 29: {
double t = 2*y100 - 59;
return 0.58702059496154081813e-1 + (0.13962391363223647892e-2 + (0.93714365487312784270e-5 + (0.58882975670265286526e-7 + (0.34414937110591753387e-9 + (0.18552853109751857859e-11 + 0.91160736711111111110e-14 * t) * t) * t) * t) * t) * t;
}
case 30: {
double t = 2*y100 - 61;
return 0.61532500145144778048e-1 + (0.14344426411912015247e-2 + (0.97331446201016809696e-5 + (0.61711860507347175097e-7 + (0.36325987418295300221e-9 + (0.19681183310134518232e-11 + 0.96952238400000000000e-14 * t) * t) * t) * t) * t) * t;
}
case 31: {
double t = 2*y100 - 63;
return 0.64440817576653297993e-1 + (0.14741275456383131151e-2 + (0.10112293819576437838e-4 + (0.64698236605933246196e-7 + (0.38353412915303665586e-9 + (0.20881176114385120186e-11 + 0.10310784480000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 32: {
double t = 2*y100 - 65;
return 0.67430045633130393282e-1 + (0.15153655418916540370e-2 + (0.10509857606888328667e-4 + (0.67851706529363332855e-7 + (0.40504602194811140006e-9 + (0.22157325110542534469e-11 + 0.10964842115555555556e-13 * t) * t) * t) * t) * t) * t;
}
case 33: {
double t = 2*y100 - 67;
return 0.70503365513338850709e-1 + (0.15582323336495709827e-2 + (0.10926868866865231089e-4 + (0.71182482239613507542e-7 + (0.42787405890153386710e-9 + (0.23514379522274416437e-11 + 0.11659571751111111111e-13 * t) * t) * t) * t) * t) * t;
}
case 34: {
double t = 2*y100 - 69;
return 0.73664114037944596353e-1 + (0.16028078812438820413e-2 + (0.11364423678778207991e-4 + (0.74701423097423182009e-7 + (0.45210162777476488324e-9 + (0.24957355004088569134e-11 + 0.12397238257777777778e-13 * t) * t) * t) * t) * t) * t;
}
case 35: {
double t = 2*y100 - 71;
return 0.76915792420819562379e-1 + (0.16491766623447889354e-2 + (0.11823685320041302169e-4 + (0.78420075993781544386e-7 + (0.47781726956916478925e-9 + (0.26491544403815724749e-11 + 0.13180196462222222222e-13 * t) * t) * t) * t) * t) * t;
}
case 36: {
double t = 2*y100 - 73;
return 0.80262075578094612819e-1 + (0.16974279491709504117e-2 + (0.12305888517309891674e-4 + (0.82350717698979042290e-7 + (0.50511496109857113929e-9 + (0.28122528497626897696e-11 + 0.14010889635555555556e-13 * t) * t) * t) * t) * t) * t;
}
case 37: {
double t = 2*y100 - 75;
return 0.83706822008980357446e-1 + (0.17476561032212656962e-2 + (0.12812343958540763368e-4 + (0.86506399515036435592e-7 + (0.53409440823869467453e-9 + (0.29856186620887555043e-11 + 0.14891851591111111111e-13 * t) * t) * t) * t) * t) * t;
}
case 38: {
double t = 2*y100 - 77;
return 0.87254084284461718231e-1 + (0.17999608886001962327e-2 + (0.13344443080089492218e-4 + (0.90900994316429008631e-7 + (0.56486134972616465316e-9 + (0.31698707080033956934e-11 + 0.15825697795555555556e-13 * t) * t) * t) * t) * t) * t;
}
case 39: {
double t = 2*y100 - 79;
return 0.90908120182172748487e-1 + (0.18544478050657699758e-2 + (0.13903663143426120077e-4 + (0.95549246062549906177e-7 + (0.59752787125242054315e-9 + (0.33656597366099099413e-11 + 0.16815130613333333333e-13 * t) * t) * t) * t) * t) * t;
}
case 40: {
double t = 2*y100 - 81;
return 0.94673404508075481121e-1 + (0.19112284419887303347e-2 + (0.14491572616545004930e-4 + (0.10046682186333613697e-6 + (0.63221272959791000515e-9 + (0.35736693975589130818e-11 + 0.17862931591111111111e-13 * t) * t) * t) * t) * t) * t;
}
case 41: {
double t = 2*y100 - 83;
return 0.98554641648004456555e-1 + (0.19704208544725622126e-2 + (0.15109836875625443935e-4 + (0.10567036667675984067e-6 + (0.66904168640019354565e-9 + (0.37946171850824333014e-11 + 0.18971959040000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 42: {
double t = 2*y100 - 85;
return 0.10255677889470089531e0 + (0.20321499629472857418e-2 + (0.15760224242962179564e-4 + (0.11117756071353507391e-6 + (0.70814785110097658502e-9 + (0.40292553276632563925e-11 + 0.20145143075555555556e-13 * t) * t) * t) * t) * t) * t;
}
case 43: {
double t = 2*y100 - 87;
return 0.10668502059865093318e0 + (0.20965479776148731610e-2 + (0.16444612377624983565e-4 + (0.11700717962026152749e-6 + (0.74967203250938418991e-9 + (0.42783716186085922176e-11 + 0.21385479360000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 44: {
double t = 2*y100 - 89;
return 0.11094484319386444474e0 + (0.21637548491908170841e-2 + (0.17164995035719657111e-4 + (0.12317915750735938089e-6 + (0.79376309831499633734e-9 + (0.45427901763106353914e-11 + 0.22696025653333333333e-13 * t) * t) * t) * t) * t) * t;
}
case 45: {
double t = 2*y100 - 91;
return 0.11534201115268804714e0 + (0.22339187474546420375e-2 + (0.17923489217504226813e-4 + (0.12971465288245997681e-6 + (0.84057834180389073587e-9 + (0.48233721206418027227e-11 + 0.24079890062222222222e-13 * t) * t) * t) * t) * t) * t;
}
case 46: {
double t = 2*y100 - 93;
return 0.11988259392684094740e0 + (0.23071965691918689601e-2 + (0.18722342718958935446e-4 + (0.13663611754337957520e-6 + (0.89028385488493287005e-9 + (0.51210161569225846701e-11 + 0.25540227111111111111e-13 * t) * t) * t) * t) * t) * t;
}
case 47: {
double t = 2*y100 - 95;
return 0.12457298393509812907e0 + (0.23837544771809575380e-2 + (0.19563942105711612475e-4 + (0.14396736847739470782e-6 + (0.94305490646459247016e-9 + (0.54366590583134218096e-11 + 0.27080225920000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 48: {
double t = 2*y100 - 97;
return 0.12941991566142438816e0 + (0.24637684719508859484e-2 + (0.20450821127475879816e-4 + (0.15173366280523906622e-6 + (0.99907632506389027739e-9 + (0.57712760311351625221e-11 + 0.28703099555555555556e-13 * t) * t) * t) * t) * t) * t;
}
case 49: {
double t = 2*y100 - 99;
return 0.13443048593088696613e0 + (0.25474249981080823877e-2 + (0.21385669591362915223e-4 + (0.15996177579900443030e-6 + (0.10585428844575134013e-8 + (0.61258809536787882989e-11 + 0.30412080142222222222e-13 * t) * t) * t) * t) * t) * t;
}
case 50: {
double t = 2*y100 - 101;
return 0.13961217543434561353e0 + (0.26349215871051761416e-2 + (0.22371342712572567744e-4 + (0.16868008199296822247e-6 + (0.11216596910444996246e-8 + (0.65015264753090890662e-11 + 0.32210394506666666666e-13 * t) * t) * t) * t) * t) * t;
}
case 51: {
double t = 2*y100 - 103;
return 0.14497287157673800690e0 + (0.27264675383982439814e-2 + (0.23410870961050950197e-4 + (0.17791863939526376477e-6 + (0.11886425714330958106e-8 + (0.68993039665054288034e-11 + 0.34101266222222222221e-13 * t) * t) * t) * t) * t) * t;
}
case 52: {
double t = 2*y100 - 105;
return 0.15052089272774618151e0 + (0.28222846410136238008e-2 + (0.24507470422713397006e-4 + (0.18770927679626136909e-6 + (0.12597184587583370712e-8 + (0.73203433049229821618e-11 + 0.36087889048888888890e-13 * t) * t) * t) * t) * t) * t;
}
case 53: {
double t = 2*y100 - 107;
return 0.15626501395774612325e0 + (0.29226079376196624949e-2 + (0.25664553693768450545e-4 + (0.19808568415654461964e-6 + (0.13351257759815557897e-8 + (0.77658124891046760667e-11 + 0.38173420035555555555e-13 * t) * t) * t) * t) * t) * t;
}
case 54: {
double t = 2*y100 - 109;
return 0.16221449434620737567e0 + (0.30276865332726475672e-2 + (0.26885741326534564336e-4 + (0.20908350604346384143e-6 + (0.14151148144240728728e-8 + (0.82369170665974313027e-11 + 0.40360957457777777779e-13 * t) * t) * t) * t) * t) * t;
}
case 55: {
double t = 2*y100 - 111;
return 0.16837910595412130659e0 + (0.31377844510793082301e-2 + (0.28174873844911175026e-4 + (0.22074043807045782387e-6 + (0.14999481055996090039e-8 + (0.87348993661930809254e-11 + 0.42653528977777777779e-13 * t) * t) * t) * t) * t) * t;
}
case 56: {
double t = 2*y100 - 113;
return 0.17476916455659369953e0 + (0.32531815370903068316e-2 + (0.29536024347344364074e-4 + (0.23309632627767074202e-6 + (0.15899007843582444846e-8 + (0.92610375235427359475e-11 + 0.45054073102222222221e-13 * t) * t) * t) * t) * t) * t;
}
case 57: {
double t = 2*y100 - 115;
return 0.18139556223643701364e0 + (0.33741744168096996041e-2 + (0.30973511714709500836e-4 + (0.24619326937592290996e-6 + (0.16852609412267750744e-8 + (0.98166442942854895573e-11 + 0.47565418097777777779e-13 * t) * t) * t) * t) * t) * t;
}
case 58: {
double t = 2*y100 - 117;
return 0.18826980194443664549e0 + (0.35010775057740317997e-2 + (0.32491914440014267480e-4 + (0.26007572375886319028e-6 + (0.17863299617388376116e-8 + (0.10403065638343878679e-10 + 0.50190265831111111110e-13 * t) * t) * t) * t) * t) * t;
}
case 59: {
double t = 2*y100 - 119;
return 0.19540403413693967350e0 + (0.36342240767211326315e-2 + (0.34096085096200907289e-4 + (0.27479061117017637474e-6 + (0.18934228504790032826e-8 + (0.11021679075323598664e-10 + 0.52931171733333333334e-13 * t) * t) * t) * t) * t) * t;
}
case 60: {
double t = 2*y100 - 121;
return 0.20281109560651886959e0 + (0.37739673859323597060e-2 + (0.35791165457592409054e-4 + (0.29038742889416172404e-6 + (0.20068685374849001770e-8 + (0.11673891799578381999e-10 + 0.55790523093333333334e-13 * t) * t) * t) * t) * t) * t;
}
case 61: {
double t = 2*y100 - 123;
return 0.21050455062669334978e0 + (0.39206818613925652425e-2 + (0.37582602289680101704e-4 + (0.30691836231886877385e-6 + (0.21270101645763677824e-8 + (0.12361138551062899455e-10 + 0.58770520160000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 62: {
double t = 2*y100 - 125;
return 0.21849873453703332479e0 + (0.40747643554689586041e-2 + (0.39476163820986711501e-4 + (0.32443839970139918836e-6 + (0.22542053491518680200e-8 + (0.13084879235290858490e-10 + 0.61873153262222222221e-13 * t) * t) * t) * t) * t) * t;
}
case 63: {
double t = 2*y100 - 127;
return 0.22680879990043229327e0 + (0.42366354648628516935e-2 + (0.41477956909656896779e-4 + (0.34300544894502810002e-6 + (0.23888264229264067658e-8 + (0.13846596292818514601e-10 + 0.65100183751111111110e-13 * t) * t) * t) * t) * t) * t;
}
case 64: {
double t = 2*y100 - 129;
return 0.23545076536988703937e0 + (0.44067409206365170888e-2 + (0.43594444916224700881e-4 + (0.36268045617760415178e-6 + (0.25312606430853202748e-8 + (0.14647791812837903061e-10 + 0.68453122631111111110e-13 * t) * t) * t) * t) * t) * t;
}
case 65: {
double t = 2*y100 - 131;
return 0.24444156740777432838e0 + (0.45855530511605787178e-2 + (0.45832466292683085475e-4 + (0.38352752590033030472e-6 + (0.26819103733055603460e-8 + (0.15489984390884756993e-10 + 0.71933206364444444445e-13 * t) * t) * t) * t) * t) * t;
}
case 66: {
double t = 2*y100 - 133;
return 0.25379911500634264643e0 + (0.47735723208650032167e-2 + (0.48199253896534185372e-4 + (0.40561404245564732314e-6 + (0.28411932320871165585e-8 + (0.16374705736458320149e-10 + 0.75541379822222222221e-13 * t) * t) * t) * t) * t) * t;
}
case 67: {
double t = 2*y100 - 135;
return 0.26354234756393613032e0 + (0.49713289477083781266e-2 + (0.50702455036930367504e-4 + (0.42901079254268185722e-6 + (0.30095422058900481753e-8 + (0.17303497025347342498e-10 + 0.79278273368888888890e-13 * t) * t) * t) * t) * t) * t;
}
case 68: {
double t = 2*y100 - 137;
return 0.27369129607732343398e0 + (0.51793846023052643767e-2 + (0.53350152258326602629e-4 + (0.45379208848865015485e-6 + (0.31874057245814381257e-8 + (0.18277905010245111046e-10 + 0.83144182364444444445e-13 * t) * t) * t) * t) * t) * t;
}
case 69: {
double t = 2*y100 - 139;
return 0.28426714781640316172e0 + (0.53983341916695141966e-2 + (0.56150884865255810638e-4 + (0.48003589196494734238e-6 + (0.33752476967570796349e-8 + (0.19299477888083469086e-10 + 0.87139049137777777779e-13 * t) * t) * t) * t) * t) * t;
}
case 70: {
double t = 2*y100 - 141;
return 0.29529231465348519920e0 + (0.56288077305420795663e-2 + (0.59113671189913307427e-4 + (0.50782393781744840482e-6 + (0.35735475025851713168e-8 + (0.20369760937017070382e-10 + 0.91262442613333333334e-13 * t) * t) * t) * t) * t) * t;
}
case 71: {
double t = 2*y100 - 143;
return 0.30679050522528838613e0 + (0.58714723032745403331e-2 + (0.62248031602197686791e-4 + (0.53724185766200945789e-6 + (0.37827999418960232678e-8 + (0.21490291930444538307e-10 + 0.95513539182222222221e-13 * t) * t) * t) * t) * t) * t;
}
case 72: {
double t = 2*y100 - 145;
return 0.31878680111173319425e0 + (0.61270341192339103514e-2 + (0.65564012259707640976e-4 + (0.56837930287837738996e-6 + (0.40035151353392378882e-8 + (0.22662596341239294792e-10 + 0.99891109760000000000e-13 * t) * t) * t) * t) * t) * t;
}
case 73: {
double t = 2*y100 - 147;
return 0.33130773722152622027e0 + (0.63962406646798080903e-2 + (0.69072209592942396666e-4 + (0.60133006661885941812e-6 + (0.42362183765883466691e-8 + (0.23888182347073698382e-10 + 0.10439349811555555556e-12 * t) * t) * t) * t) * t) * t;
}
case 74: {
double t = 2*y100 - 149;
return 0.34438138658041336523e0 + (0.66798829540414007258e-2 + (0.72783795518603561144e-4 + (0.63619220443228800680e-6 + (0.44814499336514453364e-8 + (0.25168535651285475274e-10 + 0.10901861383111111111e-12 * t) * t) * t) * t) * t) * t;
}
case 75: {
double t = 2*y100 - 151;
return 0.35803744972380175583e0 + (0.69787978834882685031e-2 + (0.76710543371454822497e-4 + (0.67306815308917386747e-6 + (0.47397647975845228205e-8 + (0.26505114141143050509e-10 + 0.11376390933333333333e-12 * t) * t) * t) * t) * t) * t;
}
case 76: {
double t = 2*y100 - 153;
return 0.37230734890119724188e0 + (0.72938706896461381003e-2 + (0.80864854542670714092e-4 + (0.71206484718062688779e-6 + (0.50117323769745883805e-8 + (0.27899342394100074165e-10 + 0.11862637614222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 77: {
double t = 2*y100 - 155;
return 0.38722432730555448223e0 + (0.76260375162549802745e-2 + (0.85259785810004603848e-4 + (0.75329383305171327677e-6 + (0.52979361368388119355e-8 + (0.29352606054164086709e-10 + 0.12360253370666666667e-12 * t) * t) * t) * t) * t) * t;
}
case 78: {
double t = 2*y100 - 157;
return 0.40282355354616940667e0 + (0.79762880915029728079e-2 + (0.89909077342438246452e-4 + (0.79687137961956194579e-6 + (0.55989731807360403195e-8 + (0.30866246101464869050e-10 + 0.12868841946666666667e-12 * t) * t) * t) * t) * t) * t;
}
case 79: {
double t = 2*y100 - 159;
return 0.41914223158913787649e0 + (0.83456685186950463538e-2 + (0.94827181359250161335e-4 + (0.84291858561783141014e-6 + (0.59154537751083485684e-8 + (0.32441553034347469291e-10 + 0.13387957943111111111e-12 * t) * t) * t) * t) * t) * t;
}
case 80: {
double t = 2*y100 - 161;
return 0.43621971639463786896e0 + (0.87352841828289495773e-2 + (0.10002929142066799966e-3 + (0.89156148280219880024e-6 + (0.62480008150788597147e-8 + (0.34079760983458878910e-10 + 0.13917107176888888889e-12 * t) * t) * t) * t) * t) * t;
}
case 81: {
double t = 2*y100 - 163;
return 0.45409763548534330981e0 + (0.91463027755548240654e-2 + (0.10553137232446167258e-3 + (0.94293113464638623798e-6 + (0.65972492312219959885e-8 + (0.35782041795476563662e-10 + 0.14455745872000000000e-12 * t) * t) * t) * t) * t) * t;
}
case 82: {
double t = 2*y100 - 165;
return 0.47282001668512331468e0 + (0.95799574408860463394e-2 + (0.11135019058000067469e-3 + (0.99716373005509038080e-6 + (0.69638453369956970347e-8 + (0.37549499088161345850e-10 + 0.15003280712888888889e-12 * t) * t) * t) * t) * t) * t;
}
case 83: {
double t = 2*y100 - 167;
return 0.49243342227179841649e0 + (0.10037550043909497071e-1 + (0.11750334542845234952e-3 + (0.10544006716188967172e-5 + (0.73484461168242224872e-8 + (0.39383162326435752965e-10 + 0.15559069118222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 84: {
double t = 2*y100 - 169;
return 0.51298708979209258326e0 + (0.10520454564612427224e-1 + (0.12400930037494996655e-3 + (0.11147886579371265246e-5 + (0.77517184550568711454e-8 + (0.41283980931872622611e-10 + 0.16122419680000000000e-12 * t) * t) * t) * t) * t) * t;
}
case 85: {
double t = 2*y100 - 171;
return 0.53453307979101369843e0 + (0.11030120618800726938e-1 + (0.13088741519572269581e-3 + (0.11784797595374515432e-5 + (0.81743383063044825400e-8 + (0.43252818449517081051e-10 + 0.16692592640000000000e-12 * t) * t) * t) * t) * t) * t;
}
case 86: {
double t = 2*y100 - 173;
return 0.55712643071169299478e0 + (0.11568077107929735233e-1 + (0.13815797838036651289e-3 + (0.12456314879260904558e-5 + (0.86169898078969313597e-8 + (0.45290446811539652525e-10 + 0.17268801084444444444e-12 * t) * t) * t) * t) * t) * t;
}
case 87: {
double t = 2*y100 - 175;
return 0.58082532122519320968e0 + (0.12135935999503877077e-1 + (0.14584223996665838559e-3 + (0.13164068573095710742e-5 + (0.90803643355106020163e-8 + (0.47397540713124619155e-10 + 0.17850211608888888889e-12 * t) * t) * t) * t) * t) * t;
}
case 88: {
double t = 2*y100 - 177;
return 0.60569124025293375554e0 + (0.12735396239525550361e-1 + (0.15396244472258863344e-3 + (0.13909744385382818253e-5 + (0.95651595032306228245e-8 + (0.49574672127669041550e-10 + 0.18435945564444444444e-12 * t) * t) * t) * t) * t) * t;
}
case 89: {
double t = 2*y100 - 179;
return 0.63178916494715716894e0 + (0.13368247798287030927e-1 + (0.16254186562762076141e-3 + (0.14695084048334056083e-5 + (0.10072078109604152350e-7 + (0.51822304995680707483e-10 + 0.19025081422222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 90: {
double t = 2*y100 - 181;
return 0.65918774689725319200e0 + (0.14036375850601992063e-1 + (0.17160483760259706354e-3 + (0.15521885688723188371e-5 + (0.10601827031535280590e-7 + (0.54140790105837520499e-10 + 0.19616655146666666667e-12 * t) * t) * t) * t) * t) * t;
}
case 91: {
double t = 2*y100 - 183;
return 0.68795950683174433822e0 + (0.14741765091365869084e-1 + (0.18117679143520433835e-3 + (0.16392004108230585213e-5 + (0.11155116068018043001e-7 + (0.56530360194925690374e-10 + 0.20209663662222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 92: {
double t = 2*y100 - 185;
return 0.71818103808729967036e0 + (0.15486504187117112279e-1 + (0.19128428784550923217e-3 + (0.17307350969359975848e-5 + (0.11732656736113607751e-7 + (0.58991125287563833603e-10 + 0.20803065333333333333e-12 * t) * t) * t) * t) * t) * t;
}
case 93: {
double t = 2*y100 - 187;
return 0.74993321911726254661e0 + (0.16272790364044783382e-1 + (0.20195505163377912645e-3 + (0.18269894883203346953e-5 + (0.12335161021630225535e-7 + (0.61523068312169087227e-10 + 0.21395783431111111111e-12 * t) * t) * t) * t) * t) * t;
}
case 94: {
double t = 2*y100 - 189;
return 0.78330143531283492729e0 + (0.17102934132652429240e-1 + (0.21321800585063327041e-3 + (0.19281661395543913713e-5 + (0.12963340087354341574e-7 + (0.64126040998066348872e-10 + 0.21986708942222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 95: {
double t = 2*y100 - 191;
return 0.81837581041023811832e0 + (0.17979364149044223802e-1 + (0.22510330592753129006e-3 + (0.20344732868018175389e-5 + (0.13617902941839949718e-7 + (0.66799760083972474642e-10 + 0.22574701262222222222e-12 * t) * t) * t) * t) * t) * t;
}
case 96: {
double t = 2*y100 - 193;
return 0.85525144775685126237e0 + (0.18904632212547561026e-1 + (0.23764237370371255638e-3 + (0.21461248251306387979e-5 + (0.14299555071870523786e-7 + (0.69543803864694171934e-10 + 0.23158593688888888889e-12 * t) * t) * t) * t) * t) * t;
}
case 97: {
double t = 2*y100 - 195;
return 0.89402868170849933734e0 + (0.19881418399127202569e-1 + (0.25086793128395995798e-3 + (0.22633402747585233180e-5 + (0.15008997042116532283e-7 + (0.72357609075043941261e-10 + 0.23737194737777777778e-12 * t) * t) * t) * t) * t) * t;
}
case 98: {
double t = 2*y100 - 197;
return 0.93481333942870796363e0 + (0.20912536329780368893e-1 + (0.26481403465998477969e-3 + (0.23863447359754921676e-5 + (0.15746923065472184451e-7 + (0.75240468141720143653e-10 + 0.24309291271111111111e-12 * t) * t) * t) * t) * t) * t;
}
case 99: {
double t = 2*y100 - 199;
return 0.97771701335885035464e0 + (0.22000938572830479551e-1 + (0.27951610702682383001e-3 + (0.25153688325245314530e-5 + (0.16514019547822821453e-7 + (0.78191526829368231251e-10 + 0.24873652355555555556e-12 * t) * t) * t) * t) * t) * t;
}
  }
  // we only get here if y = 1, i.e. |x| < 4*eps, in which case
  // erfcx is within 1e-15 of 1..
  return 1.0;
}

double Faddeeva::erfcx(double x)
{
  if (x >= 0) {
    if (x > 50) { // continued-fraction expansion is faster
      const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
      if (x > 5e7) // 1-term expansion, important to avoid overflow
	return ispi / x;
      /* 5-term expansion (rely on compiler for CSE), simplified from:
	        ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))  */
      return ispi*((x*x) * (x*x+4.5) + 2) / (x * ((x*x) * (x*x+5) + 3.75));
    }
    return erfcx_y100(400/(4+x));
  }
  else
    return x < -26.7 ? HUGE_VAL : (x < -6.1 ? 2*exp(x*x) 
				   : 2*exp(x*x) - erfcx_y100(400/(4-x)));
}

/////////////////////////////////////////////////////////////////////////
/* Compute a scaled Dawson integral 
            Faddeeva::w_im(x) = 2*Dawson(x)/sqrt(pi)
   equivalent to the imaginary part w(x) for real x.

   Uses methods similar to the erfcx calculation above: continued fractions
   for large |x|, a lookup table of Chebyshev polynomials for smaller |x|,
   and finally a Taylor expansion for |x|<0.01.
   
   Steven G. Johnson, October 2012. */

/* Given y100=100*y, where y = 1/(1+x) for x >= 0, compute w_im(x).

   Uses a look-up table of 100 different Chebyshev polynomials
   for y intervals [0,0.01], [0.01,0.02], ...., [0.99,1], generated
   with the help of Maple and a little shell script.   This allows
   the Chebyshev polynomials to be of significantly lower degree (about 1/30)
   compared to fitting the whole [0,1] interval with a single polynomial. */
static double w_im_y100(double y100, double x) {
  switch ((int) y100) {
    case 0: {
      double t = 2*y100 - 1;
      return 0.28351593328822191546e-2 + (0.28494783221378400759e-2 + (0.14427470563276734183e-4 + (0.10939723080231588129e-6 + (0.92474307943275042045e-9 + (0.89128907666450075245e-11 + 0.92974121935111111110e-13 * t) * t) * t) * t) * t) * t;
    }
    case 1: {
      double t = 2*y100 - 3;
      return 0.85927161243940350562e-2 + (0.29085312941641339862e-2 + (0.15106783707725582090e-4 + (0.11716709978531327367e-6 + (0.10197387816021040024e-8 + (0.10122678863073360769e-10 + 0.10917479678400000000e-12 * t) * t) * t) * t) * t) * t;
    }
    case 2: {
      double t = 2*y100 - 5;
      return 0.14471159831187703054e-1 + (0.29703978970263836210e-2 + (0.15835096760173030976e-4 + (0.12574803383199211596e-6 + (0.11278672159518415848e-8 + (0.11547462300333495797e-10 + 0.12894535335111111111e-12 * t) * t) * t) * t) * t) * t;
    }
    case 3: {
      double t = 2*y100 - 7;
      return 0.20476320420324610618e-1 + (0.30352843012898665856e-2 + (0.16617609387003727409e-4 + (0.13525429711163116103e-6 + (0.12515095552507169013e-8 + (0.13235687543603382345e-10 + 0.15326595042666666667e-12 * t) * t) * t) * t) * t) * t;
    }
    case 4: {
      double t = 2*y100 - 9;
      return 0.26614461952489004566e-1 + (0.31034189276234947088e-2 + (0.17460268109986214274e-4 + (0.14582130824485709573e-6 + (0.13935959083809746345e-8 + (0.15249438072998932900e-10 + 0.18344741882133333333e-12 * t) * t) * t) * t) * t) * t;
    }
    case 5: {
      double t = 2*y100 - 11;
      return 0.32892330248093586215e-1 + (0.31750557067975068584e-2 + (0.18369907582308672632e-4 + (0.15761063702089457882e-6 + (0.15577638230480894382e-8 + (0.17663868462699097951e-10 + (0.22126732680711111111e-12 + 0.30273474177737853668e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 6: {
      double t = 2*y100 - 13;
      return 0.39317207681134336024e-1 + (0.32504779701937539333e-2 + (0.19354426046513400534e-4 + (0.17081646971321290539e-6 + (0.17485733959327106250e-8 + (0.20593687304921961410e-10 + (0.26917401949155555556e-12 + 0.38562123837725712270e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 7: {
      double t = 2*y100 - 15;
      return 0.45896976511367738235e-1 + (0.33300031273110976165e-2 + (0.20423005398039037313e-4 + (0.18567412470376467303e-6 + (0.19718038363586588213e-8 + (0.24175006536781219807e-10 + (0.33059982791466666666e-12 + 0.49756574284439426165e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 8: {
      double t = 2*y100 - 17;
      return 0.52640192524848962855e-1 + (0.34139883358846720806e-2 + (0.21586390240603337337e-4 + (0.20247136501568904646e-6 + (0.22348696948197102935e-8 + (0.28597516301950162548e-10 + (0.41045502119111111110e-12 + 0.65151614515238361946e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 9: {
      double t = 2*y100 - 19;
      return 0.59556171228656770456e-1 + (0.35028374386648914444e-2 + (0.22857246150998562824e-4 + (0.22156372146525190679e-6 + (0.25474171590893813583e-8 + (0.34122390890697400584e-10 + (0.51593189879111111110e-12 + 0.86775076853908006938e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 10: {
      double t = 2*y100 - 21;
      return 0.66655089485108212551e-1 + (0.35970095381271285568e-2 + (0.24250626164318672928e-4 + (0.24339561521785040536e-6 + (0.29221990406518411415e-8 + (0.41117013527967776467e-10 + (0.65786450716444444445e-12 + 0.11791885745450623331e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 11: {
      double t = 2*y100 - 23;
      return 0.73948106345519174661e-1 + (0.36970297216569341748e-2 + (0.25784588137312868792e-4 + (0.26853012002366752770e-6 + (0.33763958861206729592e-8 + (0.50111549981376976397e-10 + (0.85313857496888888890e-12 + 0.16417079927706899860e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 12: {
      double t = 2*y100 - 25;
      return 0.81447508065002963203e-1 + (0.38035026606492705117e-2 + (0.27481027572231851896e-4 + (0.29769200731832331364e-6 + (0.39336816287457655076e-8 + (0.61895471132038157624e-10 + (0.11292303213511111111e-11 + 0.23558532213703884304e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 13: {
      double t = 2*y100 - 27;
      return 0.89166884027582716628e-1 + (0.39171301322438946014e-2 + (0.29366827260422311668e-4 + (0.33183204390350724895e-6 + (0.46276006281647330524e-8 + (0.77692631378169813324e-10 + (0.15335153258844444444e-11 + 0.35183103415916026911e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 14: {
      double t = 2*y100 - 29;
      return 0.97121342888032322019e-1 + (0.40387340353207909514e-2 + (0.31475490395950776930e-4 + (0.37222714227125135042e-6 + (0.55074373178613809996e-8 + (0.99509175283990337944e-10 + (0.21552645758222222222e-11 + 0.55728651431872687605e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 15: {
      double t = 2*y100 - 31;
      return 0.10532778218603311137e0 + (0.41692873614065380607e-2 + (0.33849549774889456984e-4 + (0.42064596193692630143e-6 + (0.66494579697622432987e-8 + (0.13094103581931802337e-9 + (0.31896187409777777778e-11 + 0.97271974184476560742e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 16: {
      double t = 2*y100 - 33;
      return 0.11380523107427108222e0 + (0.43099572287871821013e-2 + (0.36544324341565929930e-4 + (0.47965044028581857764e-6 + (0.81819034238463698796e-8 + (0.17934133239549647357e-9 + (0.50956666166186293627e-11 + (0.18850487318190638010e-12 + 0.79697813173519853340e-14 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 17: {
      double t = 2*y100 - 35;
      return 0.12257529703447467345e0 + (0.44621675710026986366e-2 + (0.39634304721292440285e-4 + (0.55321553769873381819e-6 + (0.10343619428848520870e-7 + (0.26033830170470368088e-9 + (0.87743837749108025357e-11 + (0.34427092430230063401e-12 + 0.10205506615709843189e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 18: {
      double t = 2*y100 - 37;
      return 0.13166276955656699478e0 + (0.46276970481783001803e-2 + (0.43225026380496399310e-4 + (0.64799164020016902656e-6 + (0.13580082794704641782e-7 + (0.39839800853954313927e-9 + (0.14431142411840000000e-10 + 0.42193457308830027541e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 19: {
      double t = 2*y100 - 39;
      return 0.14109647869803356475e0 + (0.48088424418545347758e-2 + (0.47474504753352150205e-4 + (0.77509866468724360352e-6 + (0.18536851570794291724e-7 + (0.60146623257887570439e-9 + (0.18533978397305276318e-10 + (0.41033845938901048380e-13 - 0.46160680279304825485e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 20: {
      double t = 2*y100 - 41;
      return 0.15091057940548936603e0 + (0.50086864672004685703e-2 + (0.52622482832192230762e-4 + (0.95034664722040355212e-6 + (0.25614261331144718769e-7 + (0.80183196716888606252e-9 + (0.12282524750534352272e-10 + (-0.10531774117332273617e-11 - 0.86157181395039646412e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 21: {
      double t = 2*y100 - 43;
      return 0.16114648116017010770e0 + (0.52314661581655369795e-2 + (0.59005534545908331315e-4 + (0.11885518333915387760e-5 + (0.33975801443239949256e-7 + (0.82111547144080388610e-9 + (-0.12357674017312854138e-10 + (-0.24355112256914479176e-11 - 0.75155506863572930844e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 22: {
      double t = 2*y100 - 45;
      return 0.17185551279680451144e0 + (0.54829002967599420860e-2 + (0.67013226658738082118e-4 + (0.14897400671425088807e-5 + (0.40690283917126153701e-7 + (0.44060872913473778318e-9 + (-0.52641873433280000000e-10 - 0.30940587864543343124e-11 * t) * t) * t) * t) * t) * t) * t;
    }
    case 23: {
      double t = 2*y100 - 47;
      return 0.18310194559815257381e0 + (0.57701559375966953174e-2 + (0.76948789401735193483e-4 + (0.18227569842290822512e-5 + (0.41092208344387212276e-7 + (-0.44009499965694442143e-9 + (-0.92195414685628803451e-10 + (-0.22657389705721753299e-11 + 0.10004784908106839254e-12 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 24: {
      double t = 2*y100 - 49;
      return 0.19496527191546630345e0 + (0.61010853144364724856e-2 + (0.88812881056342004864e-4 + (0.21180686746360261031e-5 + (0.30652145555130049203e-7 + (-0.16841328574105890409e-8 + (-0.11008129460612823934e-9 + (-0.12180794204544515779e-12 + 0.15703325634590334097e-12 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 25: {
      double t = 2*y100 - 51;
      return 0.20754006813966575720e0 + (0.64825787724922073908e-2 + (0.10209599627522311893e-3 + (0.22785233392557600468e-5 + (0.73495224449907568402e-8 + (-0.29442705974150112783e-8 + (-0.94082603434315016546e-10 + (0.23609990400179321267e-11 + 0.14141908654269023788e-12 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 26: {
      double t = 2*y100 - 53;
      return 0.22093185554845172146e0 + (0.69182878150187964499e-2 + (0.11568723331156335712e-3 + (0.22060577946323627739e-5 + (-0.26929730679360840096e-7 + (-0.38176506152362058013e-8 + (-0.47399503861054459243e-10 + (0.40953700187172127264e-11 + 0.69157730376118511127e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 27: {
      double t = 2*y100 - 55;
      return 0.23524827304057813918e0 + (0.74063350762008734520e-2 + (0.12796333874615790348e-3 + (0.18327267316171054273e-5 + (-0.66742910737957100098e-7 + (-0.40204740975496797870e-8 + (0.14515984139495745330e-10 + (0.44921608954536047975e-11 - 0.18583341338983776219e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 28: {
      double t = 2*y100 - 57;
      return 0.25058626331812744775e0 + (0.79377285151602061328e-2 + (0.13704268650417478346e-3 + (0.11427511739544695861e-5 + (-0.10485442447768377485e-6 + (-0.34850364756499369763e-8 + (0.72656453829502179208e-10 + (0.36195460197779299406e-11 - 0.84882136022200714710e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 29: {
      double t = 2*y100 - 59;
      return 0.26701724900280689785e0 + (0.84959936119625864274e-2 + (0.14112359443938883232e-3 + (0.17800427288596909634e-6 + (-0.13443492107643109071e-6 + (-0.23512456315677680293e-8 + (0.11245846264695936769e-9 + (0.19850501334649565404e-11 - 0.11284666134635050832e-12 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 30: {
      double t = 2*y100 - 61;
      return 0.28457293586253654144e0 + (0.90581563892650431899e-2 + (0.13880520331140646738e-3 + (-0.97262302362522896157e-6 + (-0.15077100040254187366e-6 + (-0.88574317464577116689e-9 + (0.12760311125637474581e-9 + (0.20155151018282695055e-12 - 0.10514169375181734921e-12 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 31: {
      double t = 2*y100 - 63;
      return 0.30323425595617385705e0 + (0.95968346790597422934e-2 + (0.12931067776725883939e-3 + (-0.21938741702795543986e-5 + (-0.15202888584907373963e-6 + (0.61788350541116331411e-9 + (0.11957835742791248256e-9 + (-0.12598179834007710908e-11 - 0.75151817129574614194e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 32: {
      double t = 2*y100 - 65;
      return 0.32292521181517384379e0 + (0.10082957727001199408e-1 + (0.11257589426154962226e-3 + (-0.33670890319327881129e-5 + (-0.13910529040004008158e-6 + (0.19170714373047512945e-8 + (0.94840222377720494290e-10 + (-0.21650018351795353201e-11 - 0.37875211678024922689e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 33: {
      double t = 2*y100 - 67;
      return 0.34351233557911753862e0 + (0.10488575435572745309e-1 + (0.89209444197248726614e-4 + (-0.43893459576483345364e-5 + (-0.11488595830450424419e-6 + (0.28599494117122464806e-8 + (0.61537542799857777779e-10 - 0.24935749227658002212e-11 * t) * t) * t) * t) * t) * t) * t;
    }
    case 34: {
      double t = 2*y100 - 69;
      return 0.36480946642143669093e0 + (0.10789304203431861366e-1 + (0.60357993745283076834e-4 + (-0.51855862174130669389e-5 + (-0.83291664087289801313e-7 + (0.33898011178582671546e-8 + (0.27082948188277716482e-10 + (-0.23603379397408694974e-11 + 0.19328087692252869842e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 35: {
      double t = 2*y100 - 71;
      return 0.38658679935694939199e0 + (0.10966119158288804999e-1 + (0.27521612041849561426e-4 + (-0.57132774537670953638e-5 + (-0.48404772799207914899e-7 + (0.35268354132474570493e-8 + (-0.32383477652514618094e-11 + (-0.19334202915190442501e-11 + 0.32333189861286460270e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 36: {
      double t = 2*y100 - 73;
      return 0.40858275583808707870e0 + (0.11006378016848466550e-1 + (-0.76396376685213286033e-5 + (-0.59609835484245791439e-5 + (-0.13834610033859313213e-7 + (0.33406952974861448790e-8 + (-0.26474915974296612559e-10 + (-0.13750229270354351983e-11 + 0.36169366979417390637e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 37: {
      double t = 2*y100 - 75;
      return 0.43051714914006682977e0 + (0.10904106549500816155e-1 + (-0.43477527256787216909e-4 + (-0.59429739547798343948e-5 + (0.17639200194091885949e-7 + (0.29235991689639918688e-8 + (-0.41718791216277812879e-10 + (-0.81023337739508049606e-12 + 0.33618915934461994428e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 38: {
      double t = 2*y100 - 77;
      return 0.45210428135559607406e0 + (0.10659670756384400554e-1 + (-0.78488639913256978087e-4 + (-0.56919860886214735936e-5 + (0.44181850467477733407e-7 + (0.23694306174312688151e-8 + (-0.49492621596685443247e-10 + (-0.31827275712126287222e-12 + 0.27494438742721623654e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 39: {
      double t = 2*y100 - 79;
      return 0.47306491195005224077e0 + (0.10279006119745977570e-1 + (-0.11140268171830478306e-3 + (-0.52518035247451432069e-5 + (0.64846898158889479518e-7 + (0.17603624837787337662e-8 + (-0.51129481592926104316e-10 + (0.62674584974141049511e-13 + 0.20055478560829935356e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 40: {
      double t = 2*y100 - 81;
      return 0.49313638965719857647e0 + (0.97725799114772017662e-2 + (-0.14122854267291533334e-3 + (-0.46707252568834951907e-5 + (0.79421347979319449524e-7 + (0.11603027184324708643e-8 + (-0.48269605844397175946e-10 + (0.32477251431748571219e-12 + 0.12831052634143527985e-13 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 41: {
      double t = 2*y100 - 83;
      return 0.51208057433416004042e0 + (0.91542422354009224951e-2 + (-0.16726530230228647275e-3 + (-0.39964621752527649409e-5 + (0.88232252903213171454e-7 + (0.61343113364949928501e-9 + (-0.42516755603130443051e-10 + (0.47910437172240209262e-12 + 0.66784341874437478953e-14 * t) * t) * t) * t) * t) * t) * t) * t;
    }
    case 42: {
      double t = 2*y100 - 85;
      return 0.52968945458607484524e0 + (0.84400880445116786088e-2 + (-0.18908729783854258774e-3 + (-0.32725905467782951931e-5 + (0.91956190588652090659e-7 + (0.14593989152420122909e-9 + (-0.35239490687644444445e-10 + 0.54613829888448694898e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 43: {
      double t = 2*y100 - 87;
      return 0.54578857454330070965e0 + (0.76474155195880295311e-2 + (-0.20651230590808213884e-3 + (-0.25364339140543131706e-5 + (0.91455367999510681979e-7 + (-0.23061359005297528898e-9 + (-0.27512928625244444444e-10 + 0.54895806008493285579e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 44: {
      double t = 2*y100 - 89;
      return 0.56023851910298493910e0 + (0.67938321739997196804e-2 + (-0.21956066613331411760e-3 + (-0.18181127670443266395e-5 + (0.87650335075416845987e-7 + (-0.51548062050366615977e-9 + (-0.20068462174044444444e-10 + 0.50912654909758187264e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 45: {
      double t = 2*y100 - 91;
      return 0.57293478057455721150e0 + (0.58965321010394044087e-2 + (-0.22841145229276575597e-3 + (-0.11404605562013443659e-5 + (0.81430290992322326296e-7 + (-0.71512447242755357629e-9 + (-0.13372664928000000000e-10 + 0.44461498336689298148e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 46: {
      double t = 2*y100 - 93;
      return 0.58380635448407827360e0 + (0.49717469530842831182e-2 + (-0.23336001540009645365e-3 + (-0.51952064448608850822e-6 + (0.73596577815411080511e-7 + (-0.84020916763091566035e-9 + (-0.76700972702222222221e-11 + 0.36914462807972467044e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 47: {
      double t = 2*y100 - 95;
      return 0.59281340237769489597e0 + (0.40343592069379730568e-2 + (-0.23477963738658326185e-3 + (0.34615944987790224234e-7 + (0.64832803248395814574e-7 + (-0.90329163587627007971e-9 + (-0.30421940400000000000e-11 + 0.29237386653743536669e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 48: {
      double t = 2*y100 - 97;
      return 0.59994428743114271918e0 + (0.30976579788271744329e-2 + (-0.23308875765700082835e-3 + (0.51681681023846925160e-6 + (0.55694594264948268169e-7 + (-0.91719117313243464652e-9 + (0.53982743680000000000e-12 + 0.22050829296187771142e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 49: {
      double t = 2*y100 - 99;
      return 0.60521224471819875444e0 + (0.21732138012345456060e-2 + (-0.22872428969625997456e-3 + (0.92588959922653404233e-6 + (0.46612665806531930684e-7 + (-0.89393722514414153351e-9 + (0.31718550353777777778e-11 + 0.15705458816080549117e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 50: {
      double t = 2*y100 - 101;
      return 0.60865189969791123620e0 + (0.12708480848877451719e-2 + (-0.22212090111534847166e-3 + (0.12636236031532793467e-5 + (0.37904037100232937574e-7 + (-0.84417089968101223519e-9 + (0.49843180828444444445e-11 + 0.10355439441049048273e-12 * t) * t) * t) * t) * t) * t) * t;
    }
    case 51: {
      double t = 2*y100 - 103;
      return 0.61031580103499200191e0 + (0.39867436055861038223e-3 + (-0.21369573439579869291e-3 + (0.15339402129026183670e-5 + (0.29787479206646594442e-7 + (-0.77687792914228632974e-9 + (0.61192452741333333334e-11 + 0.60216691829459295780e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 52: {
      double t = 2*y100 - 105;
      return 0.61027109047879835868e0 + (-0.43680904508059878254e-3 + (-0.20383783788303894442e-3 + (0.17421743090883439959e-5 + (0.22400425572175715576e-7 + (-0.69934719320045128997e-9 + (0.67152759655111111110e-11 + 0.26419960042578359995e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 53: {
      double t = 2*y100 - 107;
      return 0.60859639489217430521e0 + (-0.12305921390962936873e-2 + (-0.19290150253894682629e-3 + (0.18944904654478310128e-5 + (0.15815530398618149110e-7 + (-0.61726850580964876070e-9 + 0.68987888999111111110e-11 * t) * t) * t) * t) * t) * t;
    }
    case 54: {
      double t = 2*y100 - 109;
      return 0.60537899426486075181e0 + (-0.19790062241395705751e-2 + (-0.18120271393047062253e-3 + (0.19974264162313241405e-5 + (0.10055795094298172492e-7 + (-0.53491997919318263593e-9 + (0.67794550295111111110e-11 - 0.17059208095741511603e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 55: {
      double t = 2*y100 - 111;
      return 0.60071229457904110537e0 + (-0.26795676776166354354e-2 + (-0.16901799553627508781e-3 + (0.20575498324332621581e-5 + (0.51077165074461745053e-8 + (-0.45536079828057221858e-9 + (0.64488005516444444445e-11 - 0.29311677573152766338e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 56: {
      double t = 2*y100 - 113;
      return 0.59469361520112714738e0 + (-0.33308208190600993470e-2 + (-0.15658501295912405679e-3 + (0.20812116912895417272e-5 + (0.93227468760614182021e-9 + (-0.38066673740116080415e-9 + (0.59806790359111111110e-11 - 0.36887077278950440597e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 57: {
      double t = 2*y100 - 115;
      return 0.58742228631775388268e0 + (-0.39321858196059227251e-2 + (-0.14410441141450122535e-3 + (0.20743790018404020716e-5 + (-0.25261903811221913762e-8 + (-0.31212416519526924318e-9 + (0.54328422462222222221e-11 - 0.40864152484979815972e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 58: {
      double t = 2*y100 - 117;
      return 0.57899804200033018447e0 + (-0.44838157005618913447e-2 + (-0.13174245966501437965e-3 + (0.20425306888294362674e-5 + (-0.53330296023875447782e-8 + (-0.25041289435539821014e-9 + (0.48490437205333333334e-11 - 0.42162206939169045177e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 59: {
      double t = 2*y100 - 119;
      return 0.56951968796931245974e0 + (-0.49864649488074868952e-2 + (-0.11963416583477567125e-3 + (0.19906021780991036425e-5 + (-0.75580140299436494248e-8 + (-0.19576060961919820491e-9 + (0.42613011928888888890e-11 - 0.41539443304115604377e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 60: {
      double t = 2*y100 - 121;
      return 0.55908401930063918964e0 + (-0.54413711036826877753e-2 + (-0.10788661102511914628e-3 + (0.19229663322982839331e-5 + (-0.92714731195118129616e-8 + (-0.14807038677197394186e-9 + (0.36920870298666666666e-11 - 0.39603726688419162617e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 61: {
      double t = 2*y100 - 123;
      return 0.54778496152925675315e0 + (-0.58501497933213396670e-2 + (-0.96582314317855227421e-4 + (0.18434405235069270228e-5 + (-0.10541580254317078711e-7 + (-0.10702303407788943498e-9 + (0.31563175582222222222e-11 - 0.36829748079110481422e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 62: {
      double t = 2*y100 - 125;
      return 0.53571290831682823999e0 + (-0.62147030670760791791e-2 + (-0.85782497917111760790e-4 + (0.17553116363443470478e-5 + (-0.11432547349815541084e-7 + (-0.72157091369041330520e-10 + (0.26630811607111111111e-11 - 0.33578660425893164084e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 63: {
      double t = 2*y100 - 127;
      return 0.52295422962048434978e0 + (-0.65371404367776320720e-2 + (-0.75530164941473343780e-4 + (0.16613725797181276790e-5 + (-0.12003521296598910761e-7 + (-0.42929753689181106171e-10 + (0.22170894940444444444e-11 - 0.30117697501065110505e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 64: {
      double t = 2*y100 - 129;
      return 0.50959092577577886140e0 + (-0.68197117603118591766e-2 + (-0.65852936198953623307e-4 + (0.15639654113906716939e-5 + (-0.12308007991056524902e-7 + (-0.18761997536910939570e-10 + (0.18198628922666666667e-11 - 0.26638355362285200932e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 65: {
      double t = 2*y100 - 131;
      return 0.49570040481823167970e0 + (-0.70647509397614398066e-2 + (-0.56765617728962588218e-4 + (0.14650274449141448497e-5 + (-0.12393681471984051132e-7 + (0.92904351801168955424e-12 + (0.14706755960177777778e-11 - 0.23272455351266325318e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 66: {
      double t = 2*y100 - 133;
      return 0.48135536250935238066e0 + (-0.72746293327402359783e-2 + (-0.48272489495730030780e-4 + (0.13661377309113939689e-5 + (-0.12302464447599382189e-7 + (0.16707760028737074907e-10 + (0.11672928324444444444e-11 - 0.20105801424709924499e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 67: {
      double t = 2*y100 - 135;
      return 0.46662374675511439448e0 + (-0.74517177649528487002e-2 + (-0.40369318744279128718e-4 + (0.12685621118898535407e-5 + (-0.12070791463315156250e-7 + (0.29105507892605823871e-10 + (0.90653314645333333334e-12 - 0.17189503312102982646e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 68: {
      double t = 2*y100 - 137;
      return 0.45156879030168268778e0 + (-0.75983560650033817497e-2 + (-0.33045110380705139759e-4 + (0.11732956732035040896e-5 + (-0.11729986947158201869e-7 + (0.38611905704166441308e-10 + (0.68468768305777777779e-12 - 0.14549134330396754575e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 69: {
      double t = 2*y100 - 139;
      return 0.43624909769330896904e0 + (-0.77168291040309554679e-2 + (-0.26283612321339907756e-4 + (0.10811018836893550820e-5 + (-0.11306707563739851552e-7 + (0.45670446788529607380e-10 + (0.49782492549333333334e-12 - 0.12191983967561779442e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 70: {
      double t = 2*y100 - 141;
      return 0.42071877443548481181e0 + (-0.78093484015052730097e-2 + (-0.20064596897224934705e-4 + (0.99254806680671890766e-6 + (-0.10823412088884741451e-7 + (0.50677203326904716247e-10 + (0.34200547594666666666e-12 - 0.10112698698356194618e-13 * t) * t) * t) * t) * t) * t) * t;
    }
    case 71: {
      double t = 2*y100 - 143;
      return 0.40502758809710844280e0 + (-0.78780384460872937555e-2 + (-0.14364940764532853112e-4 + (0.90803709228265217384e-6 + (-0.10298832847014466907e-7 + (0.53981671221969478551e-10 + (0.21342751381333333333e-12 - 0.82975901848387729274e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 72: {
      double t = 2*y100 - 145;
      return 0.38922115269731446690e0 + (-0.79249269708242064120e-2 + (-0.91595258799106970453e-5 + (0.82783535102217576495e-6 + (-0.97484311059617744437e-8 + (0.55889029041660225629e-10 + (0.10851981336888888889e-12 - 0.67278553237853459757e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 73: {
      double t = 2*y100 - 147;
      return 0.37334112915460307335e0 + (-0.79519385109223148791e-2 + (-0.44219833548840469752e-5 + (0.75209719038240314732e-6 + (-0.91848251458553190451e-8 + (0.56663266668051433844e-10 + (0.23995894257777777778e-13 - 0.53819475285389344313e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 74: {
      double t = 2*y100 - 149;
      return 0.35742543583374223085e0 + (-0.79608906571527956177e-2 + (-0.12530071050975781198e-6 + (0.68088605744900552505e-6 + (-0.86181844090844164075e-8 + (0.56530784203816176153e-10 + (-0.43120012248888888890e-13 - 0.42372603392496813810e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 75: {
      double t = 2*y100 - 151;
      return 0.34150846431979618536e0 + (-0.79534924968773806029e-2 + (0.37576885610891515813e-5 + (0.61419263633090524326e-6 + (-0.80565865409945960125e-8 + (0.55684175248749269411e-10 + (-0.95486860764444444445e-13 - 0.32712946432984510595e-14 * t) * t) * t) * t) * t) * t) * t;
    }
    case 76: {
      double t = 2*y100 - 153;
      return 0.32562129649136346824e0 + (-0.79313448067948884309e-2 + (0.72539159933545300034e-5 + (0.55195028297415503083e-6 + (-0.75063365335570475258e-8 + (0.54281686749699595941e-10 - 0.13545424295111111111e-12 * t) * t) * t) * t) * t) * t;
    }
    case 77: {
      double t = 2*y100 - 155;
      return 0.30979191977078391864e0 + (-0.78959416264207333695e-2 + (0.10389774377677210794e-4 + (0.49404804463196316464e-6 + (-0.69722488229411164685e-8 + (0.52469254655951393842e-10 - 0.16507860650666666667e-12 * t) * t) * t) * t) * t) * t;
    }
    case 78: {
      double t = 2*y100 - 157;
      return 0.29404543811214459904e0 + (-0.78486728990364155356e-2 + (0.13190885683106990459e-4 + (0.44034158861387909694e-6 + (-0.64578942561562616481e-8 + (0.50354306498006928984e-10 - 0.18614473550222222222e-12 * t) * t) * t) * t) * t) * t;
    }
    case 79: {
      double t = 2*y100 - 159;
      return 0.27840427686253660515e0 + (-0.77908279176252742013e-2 + (0.15681928798708548349e-4 + (0.39066226205099807573e-6 + (-0.59658144820660420814e-8 + (0.48030086420373141763e-10 - 0.20018995173333333333e-12 * t) * t) * t) * t) * t) * t;
    }
    case 80: {
      double t = 2*y100 - 161;
      return 0.26288838011163800908e0 + (-0.77235993576119469018e-2 + (0.17886516796198660969e-4 + (0.34482457073472497720e-6 + (-0.54977066551955420066e-8 + (0.45572749379147269213e-10 - 0.20852924954666666667e-12 * t) * t) * t) * t) * t) * t;
    }
    case 81: {
      double t = 2*y100 - 163;
      return 0.24751539954181029717e0 + (-0.76480877165290370975e-2 + (0.19827114835033977049e-4 + (0.30263228619976332110e-6 + (-0.50545814570120129947e-8 + (0.43043879374212005966e-10 - 0.21228012028444444444e-12 * t) * t) * t) * t) * t) * t;
    }
    case 82: {
      double t = 2*y100 - 165;
      return 0.23230087411688914593e0 + (-0.75653060136384041587e-2 + (0.21524991113020016415e-4 + (0.26388338542539382413e-6 + (-0.46368974069671446622e-8 + (0.40492715758206515307e-10 - 0.21238627815111111111e-12 * t) * t) * t) * t) * t) * t;
    }
    case 83: {
      double t = 2*y100 - 167;
      return 0.21725840021297341931e0 + (-0.74761846305979730439e-2 + (0.23000194404129495243e-4 + (0.22837400135642906796e-6 + (-0.42446743058417541277e-8 + (0.37958104071765923728e-10 - 0.20963978568888888889e-12 * t) * t) * t) * t) * t) * t;
    }
    case 84: {
      double t = 2*y100 - 169;
      return 0.20239979200788191491e0 + (-0.73815761980493466516e-2 + (0.24271552727631854013e-4 + (0.19590154043390012843e-6 + (-0.38775884642456551753e-8 + (0.35470192372162901168e-10 - 0.20470131678222222222e-12 * t) * t) * t) * t) * t) * t;
    }
    case 85: {
      double t = 2*y100 - 171;
      return 0.18773523211558098962e0 + (-0.72822604530339834448e-2 + (0.25356688567841293697e-4 + (0.16626710297744290016e-6 + (-0.35350521468015310830e-8 + (0.33051896213898864306e-10 - 0.19811844544000000000e-12 * t) * t) * t) * t) * t) * t;
    }
    case 86: {
      double t = 2*y100 - 173;
      return 0.17327341258479649442e0 + (-0.71789490089142761950e-2 + (0.26272046822383820476e-4 + (0.13927732375657362345e-6 + (-0.32162794266956859603e-8 + (0.30720156036105652035e-10 - 0.19034196304000000000e-12 * t) * t) * t) * t) * t) * t;
    }
    case 87: {
      double t = 2*y100 - 175;
      return 0.15902166648328672043e0 + (-0.70722899934245504034e-2 + (0.27032932310132226025e-4 + (0.11474573347816568279e-6 + (-0.29203404091754665063e-8 + (0.28487010262547971859e-10 - 0.18174029063111111111e-12 * t) * t) * t) * t) * t) * t;
    }
    case 88: {
      double t = 2*y100 - 177;
      return 0.14498609036610283865e0 + (-0.69628725220045029273e-2 + (0.27653554229160596221e-4 + (0.92493727167393036470e-7 + (-0.26462055548683583849e-8 + (0.26360506250989943739e-10 - 0.17261211260444444444e-12 * t) * t) * t) * t) * t) * t;
    }
    case 89: {
      double t = 2*y100 - 179;
      return 0.13117165798208050667e0 + (-0.68512309830281084723e-2 + (0.28147075431133863774e-4 + (0.72351212437979583441e-7 + (-0.23927816200314358570e-8 + (0.24345469651209833155e-10 - 0.16319736960000000000e-12 * t) * t) * t) * t) * t) * t;
    }
    case 90: {
      double t = 2*y100 - 181;
      return 0.11758232561160626306e0 + (-0.67378491192463392927e-2 + (0.28525664781722907847e-4 + (0.54156999310046790024e-7 + (-0.21589405340123827823e-8 + (0.22444150951727334619e-10 - 0.15368675584000000000e-12 * t) * t) * t) * t) * t) * t;
    }
    case 91: {
      double t = 2*y100 - 183;
      return 0.10422112945361673560e0 + (-0.66231638959845581564e-2 + (0.28800551216363918088e-4 + (0.37758983397952149613e-7 + (-0.19435423557038933431e-8 + (0.20656766125421362458e-10 - 0.14422990012444444444e-12 * t) * t) * t) * t) * t) * t;
    }
    case 92: {
      double t = 2*y100 - 185;
      return 0.91090275493541084785e-1 + (-0.65075691516115160062e-2 + (0.28982078385527224867e-4 + (0.23014165807643012781e-7 + (-0.17454532910249875958e-8 + (0.18981946442680092373e-10 - 0.13494234691555555556e-12 * t) * t) * t) * t) * t) * t;
    }
    case 93: {
      double t = 2*y100 - 187;
      return 0.78191222288771379358e-1 + (-0.63914190297303976434e-2 + (0.29079759021299682675e-4 + (0.97885458059415717014e-8 + (-0.15635596116134296819e-8 + (0.17417110744051331974e-10 - 0.12591151763555555556e-12 * t) * t) * t) * t) * t) * t;
    }
    case 94: {
      double t = 2*y100 - 189;
      return 0.65524757106147402224e-1 + (-0.62750311956082444159e-2 + (0.29102328354323449795e-4 + (-0.20430838882727954582e-8 + (-0.13967781903855367270e-8 + (0.15958771833747057569e-10 - 0.11720175765333333333e-12 * t) * t) * t) * t) * t) * t;
    }
    case 95: {
      double t = 2*y100 - 191;
      return 0.53091065838453612773e-1 + (-0.61586898417077043662e-2 + (0.29057796072960100710e-4 + (-0.12597414620517987536e-7 + (-0.12440642607426861943e-8 + (0.14602787128447932137e-10 - 0.10885859114666666667e-12 * t) * t) * t) * t) * t) * t;
    }
    case 96: {
      double t = 2*y100 - 193;
      return 0.40889797115352738582e-1 + (-0.60426484889413678200e-2 + (0.28953496450191694606e-4 + (-0.21982952021823718400e-7 + (-0.11044169117553026211e-8 + (0.13344562332430552171e-10 - 0.10091231402844444444e-12 * t) * t) * t) * t) * t) * t;
    }
    case 97: {
      double t = 2*y100 - 195;
      return 0.28920121009594899986e-1 + (-0.59271325915413781788e-2 + (0.28796136372768177423e-4 + (-0.30300382596279568642e-7 + (-0.97688275022802329749e-9 + (0.12179215701512592356e-10 - 0.93380988481777777779e-13 * t) * t) * t) * t) * t) * t;
    }
    case 98: {
      double t = 2*y100 - 197;
      return 0.17180782722617876655e-1 + (-0.58123419543161127769e-2 + (0.28591841095380959666e-4 + (-0.37642963496443667043e-7 + (-0.86055809047367300024e-9 + (0.11101709356762665578e-10 - 0.86272947493333333334e-13 * t) * t) * t) * t) * t) * t;
    }
  case 99: case 100: { // use Taylor expansion for small x (|x| <= 0.010101...)
      //  (2/sqrt(pi)) * (x - 2/3 x^3  + 4/15 x^5  - 8/105 x^7) 
      double x2 = x*x;
      return x * (1.1283791670955125739
		  - x2 * (0.75225277806367504925
			  - x2 * (0.30090111122547001970
				  - x2 * 0.085971746064420005629)));
    }
  }
  /* Since 0 <= y100 < 101, this is only reached if x is NaN,
     in which case we should return NaN. */
  return NaN;
}

double Faddeeva::w_im(double x)
{
  if (x >= 0) {
    if (x > 45) { // continued-fraction expansion is faster
      const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
      if (x > 5e7) // 1-term expansion, important to avoid overflow
	return ispi / x;
      /* 5-term expansion (rely on compiler for CSE), simplified from:
	        ispi / (x-0.5/(x-1/(x-1.5/(x-2/x))))  */
      return ispi*((x*x) * (x*x-4.5) + 2) / (x * ((x*x) * (x*x-5) + 3.75));
    }
    return w_im_y100(100/(1+x), x);
  }
  else { // = -Faddeeva::w_im(-x)
    if (x < -45) { // continued-fraction expansion is faster
      const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
      if (x < -5e7) // 1-term expansion, important to avoid overflow
	return ispi / x;
      /* 5-term expansion (rely on compiler for CSE), simplified from:
	        ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))  */
      return ispi*((x*x) * (x*x-4.5) + 2) / (x * ((x*x) * (x*x-5) + 3.75));
    }
    return -w_im_y100(100/(1-x), -x);
  }
}

/////////////////////////////////////////////////////////////////////////

// Compile with -DTEST_FADDEEVA to compile a little test program
#ifdef TEST_FADDEEVA

#include <cstdio>

// compute relative error |b-a|/|a|, handling case of NaN and Inf,
static double relerr(double a, double b) {
  if (my_isnan(a) || my_isnan(b) || my_isinf(a) || my_isinf(b)) {
    if ((my_isnan(a) && !my_isnan(b)) || (!my_isnan(a) && my_isnan(b)) ||
	(my_isinf(a) && !my_isinf(b)) || (!my_isinf(a) && my_isinf(b)) ||
	(my_isinf(a) && my_isinf(b) && a*b < 0))
      return Inf; // "infinite" error
    return 0; // matching infinity/nan results counted as zero error
  }
  if (a == 0)
    return b == 0 ? 0 : Inf;
  else
    return fabs((b-a) / a);
}

int main(void) {
  double errmax_all = 0;
  {
    printf("############# w(z) tests #############\n");
    const int NTST = 57;
    complex<double> z[NTST] = {
      complex<double>(624.2,-0.26123),
      complex<double>(-0.4,3.),
      complex<double>(0.6,2.),
      complex<double>(-1.,1.),
      complex<double>(-1.,-9.),
      complex<double>(-1.,9.),
      complex<double>(-0.0000000234545,1.1234),
      complex<double>(-3.,5.1),
      complex<double>(-53,30.1),
      complex<double>(0.0,0.12345),
      complex<double>(11,1),
      complex<double>(-22,-2),
      complex<double>(9,-28),
      complex<double>(21,-33),
      complex<double>(1e5,1e5),
      complex<double>(1e14,1e14),
      complex<double>(-3001,-1000),
      complex<double>(1e160,-1e159),
      complex<double>(-6.01,0.01),
      complex<double>(-0.7,-0.7),
      complex<double>(2.611780000000000e+01, 4.540909610972489e+03),
      complex<double>(0.8e7,0.3e7),
      complex<double>(-20,-19.8081),
      complex<double>(1e-16,-1.1e-16),
      complex<double>(2.3e-8,1.3e-8),
      complex<double>(6.3,-1e-13),
      complex<double>(6.3,1e-20),
      complex<double>(1e-20,6.3),
      complex<double>(1e-20,16.3),
      complex<double>(9,1e-300),
      complex<double>(6.01,0.11),
      complex<double>(8.01,1.01e-10),
      complex<double>(28.01,1e-300),
      complex<double>(10.01,1e-200),
      complex<double>(10.01,-1e-200),
      complex<double>(10.01,0.99e-10),
      complex<double>(10.01,-0.99e-10),
      complex<double>(1e-20,7.01),
      complex<double>(-1,7.01),
      complex<double>(5.99,7.01),
      complex<double>(1,0),
      complex<double>(55,0),
      complex<double>(-0.1,0),
      complex<double>(1e-20,0),
      complex<double>(0,5e-14),
      complex<double>(0,51),
      complex<double>(Inf,0),
      complex<double>(-Inf,0),
      complex<double>(0,Inf),
      complex<double>(0,-Inf),
      complex<double>(Inf,Inf),
      complex<double>(Inf,-Inf),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(0,NaN),
      complex<double>(NaN,Inf),
      complex<double>(Inf,NaN)
    };
    complex<double> w[NTST] = { /* w(z), computed with WolframAlpha
				   ... note that WolframAlpha is problematic
				   some of the above inputs, so I had to
				   use the continued-fraction expansion
				   in WolframAlpha in some cases, or switch
				   to Maple */
      complex<double>(-3.78270245518980507452677445620103199303131110e-7,
		      0.000903861276433172057331093754199933411710053155),
      complex<double>(0.1764906227004816847297495349730234591778719532788,
		      -0.02146550539468457616788719893991501311573031095617),
      complex<double>(0.2410250715772692146133539023007113781272362309451,
		      0.06087579663428089745895459735240964093522265589350),
      complex<double>(0.30474420525691259245713884106959496013413834051768,
		      -0.20821893820283162728743734725471561394145872072738),
      complex<double>(7.317131068972378096865595229600561710140617977e34,
		      8.321873499714402777186848353320412813066170427e34),
      complex<double>(0.0615698507236323685519612934241429530190806818395,
		      -0.00676005783716575013073036218018565206070072304635),
      complex<double>(0.3960793007699874918961319170187598400134746631,
		      -5.593152259116644920546186222529802777409274656e-9),
      complex<double>(0.08217199226739447943295069917990417630675021771804,
		      -0.04701291087643609891018366143118110965272615832184),
      complex<double>(0.00457246000350281640952328010227885008541748668738,
		      -0.00804900791411691821818731763401840373998654987934),
      complex<double>(0.8746342859608052666092782112565360755791467973338452,
		      0.),
      complex<double>(0.00468190164965444174367477874864366058339647648741,
		      0.0510735563901306197993676329845149741675029197050),
      complex<double>(-0.0023193175200187620902125853834909543869428763219,
		      -0.025460054739731556004902057663500272721780776336),
      complex<double>(9.11463368405637174660562096516414499772662584e304,
		      3.97101807145263333769664875189354358563218932e305),
      complex<double>(-4.4927207857715598976165541011143706155432296e281,
		      -2.8019591213423077494444700357168707775769028e281),
      complex<double>(2.820947917809305132678577516325951485807107151e-6,
		      2.820947917668257736791638444590253942253354058e-6),
      complex<double>(2.82094791773878143474039725787438662716372268e-15,
		      2.82094791773878143474039725773333923127678361e-15),
      complex<double>(-0.0000563851289696244350147899376081488003110150498,
		      -0.000169211755126812174631861529808288295454992688),
      complex<double>(-5.586035480670854326218608431294778077663867e-162,
		      5.586035480670854326218608431294778077663867e-161),
      complex<double>(0.00016318325137140451888255634399123461580248456,
		      -0.095232456573009287370728788146686162555021209999),
      complex<double>(0.69504753678406939989115375989939096800793577783885,
		      -1.8916411171103639136680830887017670616339912024317),
      complex<double>(0.0001242418269653279656612334210746733213167234822,
		      7.145975826320186888508563111992099992116786763e-7),
      complex<double>(2.318587329648353318615800865959225429377529825e-8,
		      6.182899545728857485721417893323317843200933380e-8),
      complex<double>(-0.0133426877243506022053521927604277115767311800303,
		      -0.0148087097143220769493341484176979826888871576145),
      complex<double>(1.00000000000000012412170838050638522857747934,
		      1.12837916709551279389615890312156495593616433e-16),
      complex<double>(0.9999999853310704677583504063775310832036830015,
		      2.595272024519678881897196435157270184030360773e-8),
      complex<double>(-1.4731421795638279504242963027196663601154624e-15,
		      0.090727659684127365236479098488823462473074709),
      complex<double>(5.79246077884410284575834156425396800754409308e-18,
		      0.0907276596841273652364790985059772809093822374),
      complex<double>(0.0884658993528521953466533278764830881245144368,
		      1.37088352495749125283269718778582613192166760e-22),
      complex<double>(0.0345480845419190424370085249304184266813447878,
		      2.11161102895179044968099038990446187626075258e-23),
      complex<double>(6.63967719958073440070225527042829242391918213e-36,
		      0.0630820900592582863713653132559743161572639353),
      complex<double>(0.00179435233208702644891092397579091030658500743634,
		      0.0951983814805270647939647438459699953990788064762),
      complex<double>(9.09760377102097999924241322094863528771095448e-13,
		      0.0709979210725138550986782242355007611074966717),
      complex<double>(7.2049510279742166460047102593255688682910274423e-304,
		      0.0201552956479526953866611812593266285000876784321),
      complex<double>(3.04543604652250734193622967873276113872279682e-44,
		    0.0566481651760675042930042117726713294607499165),
      complex<double>(3.04543604652250734193622967873276113872279682e-44,
		      0.0566481651760675042930042117726713294607499165),
      complex<double>(0.5659928732065273429286988428080855057102069081e-12,
		      0.056648165176067504292998527162143030538756683302),
      complex<double>(-0.56599287320652734292869884280802459698927645e-12,
		      0.0566481651760675042929985271621430305387566833029),
      complex<double>(0.0796884251721652215687859778119964009569455462,
		      1.11474461817561675017794941973556302717225126e-22),
      complex<double>(0.07817195821247357458545539935996687005781943386550,
		      -0.01093913670103576690766705513142246633056714279654),
      complex<double>(0.04670032980990449912809326141164730850466208439937,
		      0.03944038961933534137558064191650437353429669886545),
      complex<double>(0.36787944117144232159552377016146086744581113103176,
		      0.60715770584139372911503823580074492116122092866515),
      complex<double>(0,
		      0.010259688805536830986089913987516716056946786526145),
      complex<double>(0.99004983374916805357390597718003655777207908125383,
		      -0.11208866436449538036721343053869621153527769495574),
      complex<double>(0.99999999999999999999999999999999999999990000,
		      1.12837916709551257389615890312154517168802603e-20),
      complex<double>(0.999999999999943581041645226871305192054749891144158,
		      0),
      complex<double>(0.0110604154853277201542582159216317923453996211744250,
		      0),
      complex<double>(0,0),
      complex<double>(0,0),
      complex<double>(0,0),
      complex<double>(Inf,0),
      complex<double>(0,0),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN)
    };
    double errmax = 0;
    for (int i = 0; i < NTST; ++i) {
      complex<double> fw = Faddeeva::w(z[i],0.);
      double re_err = relerr(real(w[i]), real(fw));
      double im_err = relerr(imag(w[i]), imag(fw));
      printf("w(%g%+gi) = %g%+gi (vs. %g%+gi), re/im rel. err. = %0.2g/%0.2g)\n",
	     real(z[i]),imag(z[i]), real(fw),imag(fw), real(w[i]),imag(w[i]),
	     re_err, im_err);
      if (re_err > errmax) errmax = re_err;
      if (im_err > errmax) errmax = im_err;
    }
    if (errmax > 1e-13) {
      printf("FAILURE -- relative error %g too large!\n", errmax);
      return 1;
    }
    printf("SUCCESS (max relative error = %g)\n", errmax);
    if (errmax > errmax_all) errmax_all = errmax;
  }
  {
    const int NTST = 33;
    complex<double> z[NTST] = {
      complex<double>(1,2),
      complex<double>(-1,2),
      complex<double>(1,-2),
      complex<double>(-1,-2),
      complex<double>(9,-28),
      complex<double>(21,-33),
      complex<double>(1e3,1e3),
      complex<double>(-3001,-1000),
      complex<double>(1e160,-1e159),
      complex<double>(5.1e-3, 1e-8),
      complex<double>(-4.9e-3, 4.95e-3),
      complex<double>(4.9e-3, 0.5),
      complex<double>(4.9e-4, -0.5e1),
      complex<double>(-4.9e-5, -0.5e2),
      complex<double>(5.1e-3, 0.5),
      complex<double>(5.1e-4, -0.5e1),
      complex<double>(-5.1e-5, -0.5e2),
      complex<double>(1e-6,2e-6),
      complex<double>(0,2e-6),
      complex<double>(0,2),
      complex<double>(0,20),
      complex<double>(0,200),
      complex<double>(Inf,0),
      complex<double>(-Inf,0),
      complex<double>(0,Inf),
      complex<double>(0,-Inf),
      complex<double>(Inf,Inf),
      complex<double>(Inf,-Inf),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(0,NaN),
      complex<double>(NaN,Inf),
      complex<double>(Inf,NaN)
    };
    complex<double> w[NTST] = { // erf(z[i]), evaluated with Maple
      complex<double>(-0.5366435657785650339917955593141927494421,
		      -5.049143703447034669543036958614140565553),
      complex<double>(0.5366435657785650339917955593141927494421,
		      -5.049143703447034669543036958614140565553),
      complex<double>(-0.5366435657785650339917955593141927494421,
		      5.049143703447034669543036958614140565553),
      complex<double>(0.5366435657785650339917955593141927494421,
		      5.049143703447034669543036958614140565553),
      complex<double>(0.3359473673830576996788000505817956637777e304,
		      -0.1999896139679880888755589794455069208455e304),
      complex<double>(0.3584459971462946066523939204836760283645e278,
		      0.3818954885257184373734213077678011282505e280),
      complex<double>(0.9996020422657148639102150147542224526887,
		      0.00002801044116908227889681753993542916894856),
      complex<double>(-1, 0),
      complex<double>(1, 0),
      complex<double>(0.005754683859034800134412990541076554934877,
		      0.1128349818335058741511924929801267822634e-7),
      complex<double>(-0.005529149142341821193633460286828381876955,
		      0.005585388387864706679609092447916333443570),
      complex<double>(0.007099365669981359632319829148438283865814,
		      0.6149347012854211635026981277569074001219),
      complex<double>(0.3981176338702323417718189922039863062440e8,
		      -0.8298176341665249121085423917575122140650e10),
      complex<double>(-Inf,
		      -Inf),
      complex<double>(0.007389128308257135427153919483147229573895,
		      0.6149332524601658796226417164791221815139),
      complex<double>(0.4143671923267934479245651547534414976991e8,
		      -0.8298168216818314211557046346850921446950e10),
      complex<double>(-Inf,
		      -Inf),
      complex<double>(0.1128379167099649964175513742247082845155e-5,
		      0.2256758334191777400570377193451519478895e-5),
      complex<double>(0,
		      0.2256758334194034158904576117253481476197e-5),
      complex<double>(0,
		      18.56480241457555259870429191324101719886),
      complex<double>(0,
		      0.1474797539628786202447733153131835124599e173),
      complex<double>(0,
		      Inf),
      complex<double>(1,0),
      complex<double>(-1,0),
      complex<double>(0,Inf),
      complex<double>(0,-Inf),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(0,NaN),
      complex<double>(NaN,NaN),
      complex<double>(NaN,NaN)
    };
#define TST(f)								\
    printf("############# " #f "(z) tests #############\n");		\
    double errmax = 0;							\
    for (int i = 0; i < NTST; ++i) {					\
      complex<double> fw = Faddeeva::f(z[i],0.);			\
      double re_err = relerr(real(w[i]), real(fw));			\
      double im_err = relerr(imag(w[i]), imag(fw));			\
      printf(#f "(%g%+gi) = %g%+gi (vs. %g%+gi), re/im rel. err. = %0.2g/%0.2g)\n", \
	     real(z[i]),imag(z[i]), real(fw),imag(fw), real(w[i]),imag(w[i]), \
	     re_err, im_err);						\
      if (re_err > errmax) errmax = re_err;				\
      if (im_err > errmax) errmax = im_err;				\
    }									\
    if (errmax > 1e-13) {						\
      printf("FAILURE -- relative error %g too large!\n", errmax);	\
      return 1;								\
    }									\
    printf("Checking " #f "(x) special case...\n");			\
    for (int i = 0; i < 10000; ++i) {					\
      double x = pow(10., -300. + i * 600. / (10000 - 1));		\
      double re_err = relerr(Faddeeva::f(x),				\
			     real(Faddeeva::f(complex<double>(x,0.))));	\
      if (re_err > errmax) errmax = re_err;				\
      re_err = relerr(Faddeeva::f(-x),					\
		      real(Faddeeva::f(complex<double>(-x,0.))));	\
      if (re_err > errmax) errmax = re_err;				\
    }									\
    {									\
      double re_err = relerr(Faddeeva::f(Inf),				\
			     real(Faddeeva::f(complex<double>(Inf,0.)))); \
      if (re_err > errmax) errmax = re_err;				\
      re_err = relerr(Faddeeva::f(-Inf),				\
		      real(Faddeeva::f(complex<double>(-Inf,0.))));	\
      if (re_err > errmax) errmax = re_err;				\
      re_err = relerr(Faddeeva::f(NaN),					\
		      real(Faddeeva::f(complex<double>(NaN,0.))));	\
      if (re_err > errmax) errmax = re_err;				\
    }									\
    if (errmax > 1e-13) {						\
      printf("FAILURE -- relative error %g too large!\n", errmax);	\
      return 1;								\
    }									\
    printf("SUCCESS (max relative error = %g)\n", errmax);		\
    if (errmax > errmax_all) errmax_all = errmax

    TST(erf);
  }
  {
    // since erfi just calls through to erf, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
    const int NTST = 1;
    complex<double> z[NTST] = { complex<double>(1.234,0.5678) };
    complex<double> w[NTST] = { // erfi(z[i]), computed with Maple
      complex<double>(1.081032284405373149432716643834106923212,
		      1.926775520840916645838949402886591180834)
    };
    TST(erfi);
  }
  {
    // since erfcx just calls through to w, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
    const int NTST = 1;
    complex<double> z[NTST] = { complex<double>(1.234,0.5678) };
    complex<double> w[NTST] = { // erfcx(z[i]), computed with Maple
      complex<double>(0.3382187479799972294747793561190487832579,
		      -0.1116077470811648467464927471872945833154)
    };
    TST(erfcx);
  }
  {
    const int NTST = 30;
    complex<double> z[NTST] = {
      complex<double>(1,2),
      complex<double>(-1,2),
      complex<double>(1,-2),
      complex<double>(-1,-2),
      complex<double>(9,-28),
      complex<double>(21,-33),
      complex<double>(1e3,1e3),
      complex<double>(-3001,-1000),
      complex<double>(1e160,-1e159),
      complex<double>(5.1e-3, 1e-8),
      complex<double>(0,2e-6),
      complex<double>(0,2),
      complex<double>(0,20),
      complex<double>(0,200),
      complex<double>(2e-6,0),
      complex<double>(2,0),
      complex<double>(20,0),
      complex<double>(200,0),
      complex<double>(Inf,0),
      complex<double>(-Inf,0),
      complex<double>(0,Inf),
      complex<double>(0,-Inf),
      complex<double>(Inf,Inf),
      complex<double>(Inf,-Inf),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(0,NaN),
      complex<double>(NaN,Inf),
      complex<double>(Inf,NaN),
      complex<double>(88,0)
    };
    complex<double> w[NTST] = { // erfc(z[i]), evaluated with Maple
      complex<double>(1.536643565778565033991795559314192749442,
		      5.049143703447034669543036958614140565553),
      complex<double>(0.4633564342214349660082044406858072505579,
		      5.049143703447034669543036958614140565553),
      complex<double>(1.536643565778565033991795559314192749442,
		      -5.049143703447034669543036958614140565553),
      complex<double>(0.4633564342214349660082044406858072505579,
		      -5.049143703447034669543036958614140565553),
      complex<double>(-0.3359473673830576996788000505817956637777e304,
		      0.1999896139679880888755589794455069208455e304),
      complex<double>(-0.3584459971462946066523939204836760283645e278,
		      -0.3818954885257184373734213077678011282505e280),
      complex<double>(0.0003979577342851360897849852457775473112748,
		      -0.00002801044116908227889681753993542916894856),
      complex<double>(2, 0),
      complex<double>(0, 0),
      complex<double>(0.9942453161409651998655870094589234450651,
		      -0.1128349818335058741511924929801267822634e-7),
      complex<double>(1,
		      -0.2256758334194034158904576117253481476197e-5),
      complex<double>(1,
		      -18.56480241457555259870429191324101719886),
      complex<double>(1,
		      -0.1474797539628786202447733153131835124599e173),
      complex<double>(1, -Inf),
      complex<double>(0.9999977432416658119838633199332831406314,
		      0),
      complex<double>(0.004677734981047265837930743632747071389108,
		      0),
      complex<double>(0.5395865611607900928934999167905345604088e-175,
		      0),
      complex<double>(0, 0),
      complex<double>(0, 0),
      complex<double>(2, 0),
      complex<double>(1, -Inf),
      complex<double>(1, Inf),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, 0),
      complex<double>(1, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(0,0)
    };
    TST(erfc);
  }
  {
    const int NTST = 47;
    complex<double> z[NTST] = {
      complex<double>(2,1),
      complex<double>(-2,1),
      complex<double>(2,-1),
      complex<double>(-2,-1),
      complex<double>(-28,9),
      complex<double>(33,-21),
      complex<double>(1e3,1e3),
      complex<double>(-1000,-3001),
      complex<double>(1e-8, 5.1e-3),
      complex<double>(4.95e-3, -4.9e-3),
      complex<double>(0.5, 4.9e-3),
      complex<double>(-0.5e1, 4.9e-4),
      complex<double>(-0.5e2, -4.9e-5),
      complex<double>(0.5e3, 4.9e-6),
      complex<double>(0.5, 5.1e-3),
      complex<double>(-0.5e1, 5.1e-4),
      complex<double>(-0.5e2, -5.1e-5),
      complex<double>(1e-6,2e-6),
      complex<double>(2e-6,0),
      complex<double>(2,0),
      complex<double>(20,0),
      complex<double>(200,0),
      complex<double>(0,4.9e-3),
      complex<double>(0,-5.1e-3),
      complex<double>(0,2e-6),
      complex<double>(0,-2),
      complex<double>(0,20),
      complex<double>(0,-200),
      complex<double>(Inf,0),
      complex<double>(-Inf,0),
      complex<double>(0,Inf),
      complex<double>(0,-Inf),
      complex<double>(Inf,Inf),
      complex<double>(Inf,-Inf),
      complex<double>(NaN,NaN),
      complex<double>(NaN,0),
      complex<double>(0,NaN),
      complex<double>(NaN,Inf),
      complex<double>(Inf,NaN),
      complex<double>(39, 6.4e-5),
      complex<double>(41, 6.09e-5),
      complex<double>(4.9e7, 5e-11),
      complex<double>(5.1e7, 4.8e-11),
      complex<double>(1e9, 2.4e-12),
      complex<double>(1e11, 2.4e-14),
      complex<double>(1e13, 2.4e-16),
      complex<double>(1e300, 2.4e-303)
    };
    complex<double> w[NTST] = { // dawson(z[i]), evaluated with Maple
      complex<double>(0.1635394094345355614904345232875688576839,
		      -0.1531245755371229803585918112683241066853),
      complex<double>(-0.1635394094345355614904345232875688576839,
		      -0.1531245755371229803585918112683241066853),
      complex<double>(0.1635394094345355614904345232875688576839,
		      0.1531245755371229803585918112683241066853),
      complex<double>(-0.1635394094345355614904345232875688576839,
		      0.1531245755371229803585918112683241066853),
      complex<double>(-0.01619082256681596362895875232699626384420,
		      -0.005210224203359059109181555401330902819419),
      complex<double>(0.01078377080978103125464543240346760257008,
		      0.006866888783433775382193630944275682670599),
      complex<double>(-0.5808616819196736225612296471081337245459,
		      0.6688593905505562263387760667171706325749),
      complex<double>(Inf,
		      -Inf),
      complex<double>(0.1000052020902036118082966385855563526705e-7,
		      0.005100088434920073153418834680320146441685),
      complex<double>(0.004950156837581592745389973960217444687524,
		      -0.004899838305155226382584756154100963570500),
      complex<double>(0.4244534840871830045021143490355372016428,
		      0.002820278933186814021399602648373095266538),
      complex<double>(-0.1021340733271046543881236523269967674156,
		      -0.00001045696456072005761498961861088944159916),
      complex<double>(-0.01000200120119206748855061636187197886859,
		      0.9805885888237419500266621041508714123763e-8),
      complex<double>(0.001000002000012000023960527532953151819595,
		      -0.9800058800588007290937355024646722133204e-11),
      complex<double>(0.4244549085628511778373438768121222815752,
		      0.002935393851311701428647152230552122898291),
      complex<double>(-0.1021340732357117208743299813648493928105,
		      -0.00001088377943049851799938998805451564893540),
      complex<double>(-0.01000200120119126652710792390331206563616,
		      0.1020612612857282306892368985525393707486e-7),
      complex<double>(0.1000000000007333333333344266666666664457e-5,
		      0.2000000000001333333333323199999999978819e-5),
      complex<double>(0.1999999999994666666666675199999999990248e-5,
		      0),
      complex<double>(0.3013403889237919660346644392864226952119,
		      0),
      complex<double>(0.02503136792640367194699495234782353186858,
		      0),
      complex<double>(0.002500031251171948248596912483183760683918,
		      0),
      complex<double>(0,0.004900078433419939164774792850907128053308),
      complex<double>(0,-0.005100088434920074173454208832365950009419),
      complex<double>(0,0.2000000000005333333333341866666666676419e-5),
      complex<double>(0,-48.16001211429122974789822893525016528191),
      complex<double>(0,0.4627407029504443513654142715903005954668e174),
      complex<double>(0,-Inf),
      complex<double>(0,0),
      complex<double>(-0,0),
      complex<double>(0, Inf),
      complex<double>(0, -Inf),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, 0),
      complex<double>(0, NaN),
      complex<double>(NaN, NaN),
      complex<double>(NaN, NaN),
      complex<double>(0.01282473148489433743567240624939698290584,
		      -0.2105957276516618621447832572909153498104e-7),
      complex<double>(0.01219875253423634378984109995893708152885,
		      -0.1813040560401824664088425926165834355953e-7),
      complex<double>(0.1020408163265306334945473399689037886997e-7,
		      -0.1041232819658476285651490827866174985330e-25),
      complex<double>(0.9803921568627452865036825956835185367356e-8,
		      -0.9227220299884665067601095648451913375754e-26),
      complex<double>(0.5000000000000000002500000000000000003750e-9,
		      -0.1200000000000000001800000188712838420241e-29),
      complex<double>(5.00000000000000000000025000000000000000000003e-12,
		      -1.20000000000000000000018000000000000000000004e-36),
      complex<double>(5.00000000000000000000000002500000000000000000e-14,
		      -1.20000000000000000000000001800000000000000000e-42),
      complex<double>(5e-301, 0)
    };
    TST(Dawson);
  }
  printf("#####################################\n");
  printf("SUCCESS (max relative error = %g)\n", errmax_all);
}

#endif
