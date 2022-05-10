/*************************** wnchyppr.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2013-12-20
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Calculation of univariate and multivariate Wallenius noncentral 
* hypergeometric probability distribution.
*
* This file contains source code for the class CWalleniusNCHypergeometric 
* and CMultiWalleniusNCHypergeometricMoments defined in stocc.h.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* The file nchyp.pdf, available from www.agner.org/random/theory 
* describes the theory of the calculation methods.
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2002-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <string.h>                    // memcpy function
#include "stocc.h"                     // class definition
#include "erfres.cpp"                  // table of error function residues (Don't precompile this as a header file)

/***********************************************************************
constants
***********************************************************************/
static const double LN2 = 0.693147180559945309417; // log(2)


/***********************************************************************
Log and Exp functions with special care for small x
***********************************************************************/
// These are functions that involve expressions of the types log(1+x)
// and exp(x)-1. These functions need special care when x is small to
// avoid loss of precision. There are three versions of these functions:
// (1) Assembly version in library randomaXX.lib
// (2) Use library functions log1p and expm1 if available
// (3) Use Taylor expansion if none of the above are available

#if defined(__GNUC__) || defined (__INTEL_COMPILER) || defined(HAVE_EXPM1)
// Functions log1p(x) = log(1+x) and expm1(x) = exp(x)-1 are available
// in the math libraries of Gnu and Intel compilers
// and in R.DLL (www.r-project.org).

double pow2_1(double q, double * y0 = 0) {
   // calculate 2^q and (1-2^q) without loss of precision.
   // return value is (1-2^q). 2^q is returned in *y0
   double y, y1;
   q *= LN2;
   if (fabs(q) > 0.1) {
      y = exp(q);                      // 2^q
      y1 = 1. - y;                     // 1-2^q
   }
   else { // Use expm1
      y1 = expm1(q);                   // 2^q-1
      y = y1 + 1;                      // 2^q
      y1 = -y1;                        // 1-2^q
   }
   if (y0) *y0 = y;                    // Return y if not void pointer
   return y1;                          // Return y1
}

double log1mx(double x, double x1) {
   // Calculate log(1-x) without loss of precision when x is small.
   // Parameter x1 must be = 1-x.
   if (fabs(x) > 0.03) {
      return log(x1);
   }
   else { // use log1p(x) = log(1+x)
      return log1p(-x);
   }
}

double log1pow(double q, double x) {
   // calculate log((1-e^q)^x) without loss of precision.
   // Combines the methods of the above two functions.
   double y, y1;

   if (fabs(q) > 0.1) {
      y = exp(q);                      // e^q
      y1 = 1. - y;                     // 1-e^q
   }
   else { // Use expm1
      y1 = expm1(q);                   // e^q-1
      y = y1 + 1;                      // e^q
      y1 = -y1;                        // 1-e^q
   }

   if (y > 0.1) { // (1-y)^x calculated without problem
      return x * log(y1);
   }
   else { // Use log1p
      return x * log1p(-y);
   }
}

#else
// (3)
// Functions log1p and expm1 are not available in MS and Borland compiler
// libraries. Use explicit Taylor expansion when needed.

double pow2_1(double q, double * y0 = 0) {
   // calculate 2^q and (1-2^q) without loss of precision.
   // return value is (1-2^q). 2^q is returned in *y0
   double y, y1, y2, qn, i, ifac;
   q *= LN2;
   if (fabs(q) > 0.1) {
      y = exp(q);
      y1 = 1. - y;
   }
   else { // expand 1-e^q = -summa(q^n/n!) to avoid loss of precision
      y1 = 0;  qn = i = ifac = 1;
      do {
         y2 = y1;
         qn *= q;  ifac *= i++;
         y1 -= qn / ifac;
      }
      while (y1 != y2);
      y = 1.-y1;
   }
   if (y0) *y0 = y;
   return y1;
}

double log1mx(double x, double x1) {
   // Calculate log(1-x) without loss of precision when x is small.
   // Parameter x1 must be = 1-x.
   if (fabs(x) > 0.03) {
      return log(x1);
   }
   else { // expand ln(1-x) = -summa(x^n/n)
      double y, z1, z2, i;
      y = i = 1.;  z1 = 0;
      do {
         z2 = z1;
         y *= x;
         z1 -= y / i++;
      }
      while (z1 != z2);
      return z1;
   }
}

double log1pow(double q, double x) {
   // calculate log((1-e^q)^x) without loss of precision
   // Uses various Taylor expansions to avoid loss of precision
   double y, y1, y2, z1, z2, qn, i, ifac;

   if (fabs(q) > 0.1) {
      y = exp(q);  y1 = 1. - y;
   }
   else { // expand 1-e^q = -summa(q^n/n!) to avoid loss of precision
      y1 = 0;  qn = i = ifac = 1;
      do {
         y2 = y1;
         qn *= q;  ifac *= i++;
         y1 -= qn / ifac;
      }
      while (y1 != y2);
      y = 1. - y1;
   }
   if (y > 0.1) { // (1-y)^x calculated without problem
      return x * log(y1);
   }
   else { // expand ln(1-y) = -summa(y^n/n)
      y1 = i = 1.;  z1 = 0;
      do {
         z2 = z1;
         y1 *= y;
         z1 -= y1 / i++;
      }
      while (z1 != z2);
      return x * z1;
   }
}

#endif

/***********************************************************************
Other shared functions
***********************************************************************/

double LnFacr(double x) {
   // log factorial of non-integer x
   int32_t ix = (int32_t)(x);
   if (x == ix) return LnFac(ix);      // x is integer
   double r, r2, D = 1., f;
   static const double             
      C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
      C1 =  1./12.,
      C3 = -1./360.,
      C5 =  1./1260.,
      C7 = -1./1680.;
   if (x < 6.) {
      if (x == 0 || x == 1) return 0;
      while (x < 6) D *= ++x;
   }
   r  = 1. / x;  r2 = r*r;
   f = (x + 0.5)*log(x) - x + C0 + r*(C1 + r2*(C3 + r2*(C5 + r2*C7)));
   if (D != 1.) f -= log(D);
   return f;
}


double FallingFactorial(double a, double b) {
   // calculates ln(a*(a-1)*(a-2)* ... * (a-b+1))

   if (b < 30 && int(b) == b && a < 1E10) {
      // direct calculation
      double f = 1.;
      for (int i = 0; i < b; i++) f *= a--;
      return log(f);
   }

   if (a > 100.*b && b > 1.) {
      // combine Stirling formulas for a and (a-b) to avoid loss of precision
      double ar = 1./a;
      double cr = 1./(a-b);
      // calculate -log(1-b/a) by Taylor expansion
      double s = 0., lasts, n = 1., ba = b*ar, f = ba;
      do {
         lasts = s;
         s += f/n;
         f *= ba;
         n++;
      }
      while (s != lasts);
      return (a+0.5)*s + b*log(a-b) - b + (1./12.)*(ar-cr)    
         //- (1./360.)*(ar*ar*ar-cr*cr*cr)
         ;
   }
   // use LnFacr function
   return LnFacr(a)-LnFacr(a-b);
}

double Erf (double x) {
   // Calculates the error function erf(x) as a series expansion or
   // continued fraction expansion.
   // This function may be available in math libraries as erf(x)
   static const double rsqrtpi  = 0.564189583547756286948; // 1/sqrt(pi)
   static const double rsqrtpi2 = 1.12837916709551257390;  // 2/sqrt(pi)
   if (x < 0.) return -Erf(-x);
   if (x > 6.) return 1.;
   if (x < 2.4) {
      // use series expansion
      double term;                     // term in summation
      double j21;                      // 2j+1
      double sum = 0;                  // summation
      double xx2 = x*x*2.;             // 2x^2
      int j;  
      term = x;  j21 = 1.;
      for (j=0; j<=50; j++) {          // summation loop
         sum += term;
         if (term <= 1.E-13) break;
         j21 += 2.;
         term *= xx2 / j21;
      }
      return exp(-x*x) * sum * rsqrtpi2;
   }
   else {
      // use continued fraction expansion
      double a, f;
      int n = int(2.25f*x*x - 23.4f*x + 60.84f); // predict expansion degree
      if (n < 1) n = 1;
      a = 0.5 * n;  f = x;
      for (; n > 0; n--) {             // continued fraction loop
         f = x + a / f;
         a -= 0.5;
      }
      return 1. - exp(-x*x) * rsqrtpi / f;
   }
}


int32_t FloorLog2(float x) {
   // This function calculates floor(log2(x)) for positive x.
   // The return value is <= -127 for x <= 0.

   union UfloatInt {  // Union for extracting bits from a float
      float   f;
      int32_t i;
      UfloatInt(float ff) {f = ff;}  // constructor
   };

#if defined(_M_IX86) || defined(__INTEL__) || defined(_M_X64) || defined(__IA64__) || defined(__POWERPC__)
   // Running on a platform known to use IEEE-754 floating point format
   //int32_t n = *(int32_t*)&x;
   int32_t n = UfloatInt(x).i;
   return (n >> 23) - 0x7F;
#else
   // Check if floating point format is IEEE-754
   static const UfloatInt check(1.0f);
   if (check.i == 0x3F800000) {
      // We have the standard IEEE floating point format
      int32_t n = UfloatInt(x).i;
      return (n >> 23) - 0x7F;
   }
   else {
      // Unknown floating point format
      if (x <= 0.f) return -127;
      return (int32_t)floor(log(x)*(1./LN2));
   }
#endif
}


int NumSD (double accuracy) {
   // Gives the length of the integration interval necessary to achieve
   // the desired accuracy when integrating/summating a probability 
   // function, relative to the standard deviation
   // Returns an integer approximation to 2*NormalDistrFractile(accuracy/2)
   static const double fract[] = {
      2.699796e-03, 4.652582e-04, 6.334248e-05, 6.795346e-06, 5.733031e-07,
      3.797912e-08, 1.973175e-09, 8.032001e-11, 2.559625e-12, 6.381783e-14
   };
   int i;
   for (i = 0; i < (int)(sizeof(fract)/sizeof(*fract)); i++) {
      if (accuracy >= fract[i]) break;
   }
   return i + 6;
}


/***********************************************************************
Methods for class CWalleniusNCHypergeometric
***********************************************************************/

CWalleniusNCHypergeometric::CWalleniusNCHypergeometric(int32_t n_, int32_t m_, int32_t N_, double odds_, double accuracy_) {
   // constructor
   accuracy = accuracy_;
   SetParameters(n_, m_, N_, odds_);
}


void CWalleniusNCHypergeometric::SetParameters(int32_t n_, int32_t m_, int32_t N_, double odds) {
   // change parameters
   if (n_ < 0 || n_ > N_ || m_ < 0 || m_ > N_ || odds < 0) FatalError("Parameter out of range in CWalleniusNCHypergeometric");
   n = n_; m = m_; N = N_; omega = odds;          // set parameters
   xmin = m + n - N;  if (xmin < 0) xmin = 0;     // calculate xmin
   xmax = n;  if (xmax > m) xmax = m;             // calculate xmax
   xLastBico = xLastFindpars = -99;               // indicate last x is invalid
   r = 1.;                                        // initialize
}


double CWalleniusNCHypergeometric::mean(void) {
   // find approximate mean
   int iter;                            // number of iterations
   double a, b;                         // temporaries in calculation of first guess
   double mean, mean1;                  // iteration value of mean
   double m1r, m2r;                     // 1/m, 1/m2
   double e1, e2;                       // temporaries
   double g;                            // function to find root of
   double gd;                           // derivative of g
   double omegar;                       // 1/omega

   if (omega == 1.) { // simple hypergeometric
      return (double)(m)*n/N;
   }

   if (omega == 0.) {
      if (n > N-m) FatalError("Not enough items with nonzero weight in CWalleniusNCHypergeometric::mean");
      return 0.;
   }

   if (xmin == xmax) return xmin;

   // calculate Cornfield mean of Fisher noncentral hypergeometric distribution as first guess
   a = (m+n)*omega + (N-m-n); 
   b = a*a - 4.*omega*(omega-1.)*m*n;
   b = b > 0. ? sqrt(b) : 0.;
   mean = (a-b)/(2.*(omega-1.));
   if (mean < xmin) mean = xmin;
   if (mean > xmax) mean = xmax;

   m1r = 1./m;  m2r = 1./(N-m);
   iter = 0;

   if (omega > 1.) {
      do { // Newton Raphson iteration
         mean1 = mean;
         e1 = 1.-(n-mean)*m2r;
         if (e1 < 1E-14) {
            e2 = 0.;     // avoid underflow
         }
         else {
            e2 = pow(e1,omega-1.);
         }
         g = e2*e1 + (mean-m)*m1r;
         gd = e2*omega*m2r + m1r;
         mean -= g / gd;
         if (mean < xmin) mean = xmin;
         if (mean > xmax) mean = xmax;
         if (++iter > 40) {
            FatalError("Search for mean failed in function CWalleniusNCHypergeometric::mean");
         }
      }
      while (fabs(mean1 - mean) > 2E-6);
   }
   else { // omega < 1
      omegar = 1./omega;
      do { // Newton Raphson iteration
         mean1 = mean;
         e1 = 1.-mean*m1r;
         if (e1 < 1E-14) {
            e2 = 0.;   // avoid underflow
         }
         else {
            e2 = pow(e1,omegar-1.);
         }
         g = 1.-(n-mean)*m2r-e2*e1;
         gd = e2*omegar*m1r + m2r;
         mean -= g / gd;
         if (mean < xmin) mean = xmin;
         if (mean > xmax) mean = xmax;
         if (++iter > 40) {
            FatalError("Search for mean failed in function CWalleniusNCHypergeometric::mean");
         }
      }
      while (fabs(mean1 - mean) > 2E-6);
   }
   return mean;
}


double CWalleniusNCHypergeometric::variance(void) {
   // find approximate variance (poor approximation)    
   double my = mean(); // approximate mean
   // find approximate variance from Fisher's noncentral hypergeometric approximation
   double r1 = my * (m-my); double r2 = (n-my)*(my+N-n-m);
   if (r1 <= 0. || r2 <= 0.) return 0.;
   double var = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));
   if (var < 0.) var = 0.;
   return var;
}


double CWalleniusNCHypergeometric::moments(double * mean_, double * var_) {
   // calculate exact mean and variance
   // return value = sum of f(x), expected = 1.
   double y, sy=0, sxy=0, sxxy=0, me1;
   int32_t x, xm, x1;
   const double accuracy = 1E-10f;  // accuracy of calculation
   xm = (int32_t)mean();  // approximation to mean
   for (x = xm; x <= xmax; x++) {
      y = probability(x);
      x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
      sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
      if (y < accuracy && x != xm) break;
   }
   for (x = xm-1; x >= xmin; x--) {
      y = probability(x);
      x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
      sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
      if (y < accuracy) break;
   }

   me1 = sxy / sy;
   *mean_ = me1 + xm;
   y = sxxy / sy - me1 * me1;
   if (y < 0) y=0;
   *var_ = y;
   return sy;
}


int32_t CWalleniusNCHypergeometric::mode(void) {
   // find mode
   int32_t Mode;                       // mode

   if (omega == 1.) { 
      // simple hypergeometric
      int32_t L  = m + n - N;
      int32_t m1 = m + 1, n1 = n + 1;
      Mode = int32_t((double)m1*n1*omega/((m1+n1)*omega-L));
   }
   else {
      // find mode
      double f, f2 = 0.; // f2 = -1.; 
      int32_t xi, x2;
      int32_t xmin = m + n - N;  if (xmin < 0) xmin = 0;  // calculate xmin
      int32_t xmax = n;  if (xmax > m) xmax = m;          // calculate xmax

      Mode = (int32_t)mean();                             // floor(mean)
      if (omega < 1.) {
        if (Mode < xmax) Mode++;                        // ceil(mean)
        x2 = xmin;                                      // lower limit
        if (omega > 0.294 && N <= 10000000) {
          x2 = Mode - 1;}                    // search for mode can be limited
        for (xi = Mode; xi >= x2; xi--) {
          f = probability(xi);
          if (f <= f2) break;
          Mode = xi;  f2 = f;
        }
      }
      else {
        if (Mode < xmin) Mode++; 
        x2 = xmax;                           // upper limit
        if (omega < 3.4 && N <= 10000000) {
          x2 = Mode + 1;}                    // search for mode can be limited
        for (xi = Mode; xi <= x2; xi++) {
          f = probability(xi);
          if (f <= f2) break;
          Mode = xi; f2 = f;
        }
      }
   }
   return Mode;
}


double CWalleniusNCHypergeometric::lnbico() {
   // natural log of binomial coefficients.
   // returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!)
   int32_t x2 = n-x, m2 = N-m;
   if (xLastBico < 0) { // m, n, N have changed
      mFac = LnFac(m) + LnFac(m2);
   }
   if (m < FAK_LEN && m2 < FAK_LEN)  goto DEFLT;
   switch (x - xLastBico) {
  case 0: // x unchanged
     break;
  case 1: // x incremented. calculate from previous value
     xFac += log (double(x) * (m2-x2) / (double(x2+1)*(m-x+1)));
     break;
  case -1: // x decremented. calculate from previous value
     xFac += log (double(x2) * (m-x) / (double(x+1)*(m2-x2+1)));
     break;
  default: DEFLT: // calculate all
     xFac = LnFac(x) + LnFac(x2) + LnFac(m-x) + LnFac(m2-x2);
   }
   xLastBico = x;
   return bico = mFac - xFac;
}


void CWalleniusNCHypergeometric::findpars() {
   // calculate d, E, r, w
   if (x == xLastFindpars) {
      return;    // all values are unchanged since last call
   }

   // find r to center peak of integrand at 0.5
   double dd, d1, z, zd, rr, lastr, rrc, rt, r2, r21, a, b, dummy;
   double oo[2];
   double xx[2] = {static_cast<double>(x), static_cast<double>(n-x)};
   int i, j = 0;
   if (omega > 1.) { // make both omegas <= 1 to avoid overflow
      oo[0] = 1.;  oo[1] = 1./omega;
   }
   else {
      oo[0] = omega;  oo[1] = 1.;
   }
   dd = oo[0]*(m-x) + oo[1]*(N-m-xx[1]);
   d1 = 1./dd;
   E = (oo[0]*m + oo[1]*(N-m)) * d1;
   rr = r;
   if (rr <= d1) rr = 1.2*d1;           // initial guess
   // Newton-Raphson iteration to find r
   do {
      lastr = rr;
      rrc = 1. / rr;
      z = dd - rrc;
      zd = rrc * rrc;
      for (i=0; i<2; i++) {
         rt = rr * oo[i];
         if (rt < 100.) {                  // avoid overflow if rt big
            r21 = pow2_1(rt, &r2);         // r2=2^r, r21=1.-2^r
            a = oo[i] / r21;               // omegai/(1.-2^r)
            b = xx[i] * a;                 // x*omegai/(1.-2^r)
            z  += b;
            zd += b * a * LN2 * r2;
         }
      }
      if (zd == 0) FatalError("can't find r in function CWalleniusNCHypergeometric::findpars");
      rr -= z / zd;
      if (rr <= d1) rr = lastr * 0.125 + d1*0.875;
      if (++j == 70) FatalError("convergence problem searching for r in function CWalleniusNCHypergeometric::findpars");
   }
   while (fabs(rr-lastr) > rr * 1.E-6);
   if (omega > 1) {
      dd *= omega;  rr *= oo[1];
   }
   r = rr;  rd = rr * dd;

   // find peak width
   double ro, k1, k2;
   ro = r * omega;
   if (ro < 300) {                      // avoid overflow
      k1 = pow2_1(ro, &dummy);
      k1 = -1. / k1;
      k1 = omega*omega*(k1+k1*k1);
   }
   else k1 = 0.;
   if (r < 300) {                       // avoid overflow
      k2 = pow2_1(r, &dummy);
      k2 = -1. / k2;
      k2 = (k2+k2*k2);
   }
   else k2 = 0.;
   phi2d = -4.*r*r*(x*k1 + (n-x)*k2);
   if (phi2d >= 0.) {
      FatalError("peak width undefined in function CWalleniusNCHypergeometric::findpars");
      /* wr = r = 0.; */ 
   }
   else {
      wr = sqrt(-phi2d); w = 1./wr;
   }
   xLastFindpars = x;
}


int CWalleniusNCHypergeometric::BernouilliH(int32_t x_, double h, double rh, StochasticLib1 *sto) {
   // This function generates a Bernouilli variate with probability proportional
   // to the univariate Wallenius' noncentral hypergeometric distribution.
   // The return value will be 1 with probability f(x_)/h and 0 with probability
   // 1-f(x_)/h.
   // This is equivalent to calling sto->Bernouilli(probability(x_)/h),
   // but this method is faster. The method used here avoids calculating the
   // Wallenius probability by sampling in the t-domain.
   // rh is a uniform random number in the interval 0 <= rh < h. The function
   // uses additional random numbers generated from sto.
   // This function is intended for use in rejection methods for sampling from
   // the Wallenius distribution. It is called from 
   // StochasticLib3::WalleniusNCHypRatioOfUnifoms in the file stoc3.cpp
   double f0;                      // Lambda*Phi(.5)
   double phideri0;                // phi(.5)/rd
   double qi;                      // 2^(-r*omega[i])
   double qi1;                     // 1-qi
   double omegai[2] = {omega,1.};  // weights for each color
   double romegi;                  // r*omega[i]
   double xi[2] =                  // number of each color sampled
       {static_cast<double>(x_), static_cast<double>(n-x_)};
   double k;                       // adjusted width for majorizing function Ypsilon(t)
   double erfk;                    // erf correction
   double rdm1;                    // rd - 1
   double G_integral;              // integral of majorizing function Ypsilon(t)
   double ts;                      // t sampled from Ypsilon(t) distribution
   double logts;                   // log(ts)
   double rlogts;                  // r*log(ts)
   double fts;                     // Phi(ts)/rd
   double rgts;                    // 1/(Ypsilon(ts)/rd)
   double t2;                      // temporary in calculation of Ypsilon(ts)
   int i, j;                       // loop counters
   static const double rsqrt8 = 0.3535533905932737622; // 1/sqrt(8)
   static const double sqrt2pi = 2.506628274631000454; // sqrt(2*pi)

   x = x_;                         // save x in class object
   lnbico();                       // calculate bico = log(Lambda)
   findpars();                     // calculate r, d, rd, w, E
   if (E > 0.) {
      k = log(E);                  // correction for majorizing function
      k = 1. + 0.0271 * (k * sqrt(k));
   }
   else k = 1.;
   k *= w;                         // w * k   
   rdm1 = rd - 1.;

   // calculate phi(.5)/rd
   phideri0 = -LN2 * rdm1;
   for (i=0; i<2; i++) {
      romegi = r * omegai[i];
      if (romegi > 40.) {
         qi=0.;  qi1 = 1.;           // avoid underflow
      }
      else {
         qi1 = pow2_1(-romegi, &qi);
      }
      phideri0 += xi[i] * log1mx(qi, qi1);
   }

   erfk = Erf(rsqrt8 / k);
   f0 = rd * exp(phideri0 + bico);
   G_integral = f0 * sqrt2pi * k * erfk;

   if (G_integral <= h) {          // G fits under h-hat
      do {
         ts = sto->Normal(0,k);      // sample ts from normal distribution
      }
      while (fabs(ts) >= 0.5);      // reject values outside interval, and avoid ts = 0
      ts += 0.5;                    // ts = normal distributed in interval (0,1)

      for (fts=0., j=0; j<2; j++) { // calculate (Phi(ts)+Phi(1-ts))/2
         logts = log(ts);  rlogts = r * logts; // (ts = 0 avoided above)
         fts += exp(log1pow(rlogts*omega,xi[0]) + log1pow(rlogts,xi[1]) + rdm1*logts + bico);
         ts = 1. - ts;
      }
      fts *= 0.5;

      t2 = (ts-0.5) / k;            // calculate 1/Ypsilon(ts)
      rgts = exp(-(phideri0 + bico - 0.5 * t2*t2));
      return rh < G_integral * fts * rgts;   // Bernouilli variate
   }

   else { // G > h: can't use sampling in t-domain
      return rh < probability(x);
   }
}


/***********************************************************************
methods for calculating probability in class CWalleniusNCHypergeometric
***********************************************************************/

double CWalleniusNCHypergeometric::recursive() {
   // recursive calculation
   // Wallenius noncentral hypergeometric distribution by recursion formula
   // Approximate by ignoring probabilities < accuracy and minimize storage requirement
   const int BUFSIZE = 512;            // buffer size
   double p[BUFSIZE+2];                // probabilities
   double * p1, * p2;                  // offset into p
   double mxo;                         // (m-x)*omega
   double Nmnx;                        // N-m-nu+x
   double y, y1;                       // save old p[x] before it is overwritten
   double d1, d2, dcom;                // divisors in probability formula
   double accuracya;                   // absolute accuracy
   int32_t xi, nu;                     // xi, nu = recursion values of x, n
   int32_t x1, x2;                     // xi_min, xi_max

   accuracya = 0.005f * accuracy;      // absolute accuracy
   p1 = p2 = p + 1;                    // make space for p1[-1]
   p1[-1] = 0.;  p1[0]  = 1.;          // initialize for recursion
   x1 = x2 = 0;
   for (nu=1; nu<=n; nu++) {
      if (n - nu < x - x1 || p1[x1] < accuracya) {
         x1++;                        // increase lower limit when breakpoint passed or probability negligible
         p2--;                        // compensate buffer offset in order to reduce storage space
      }
      if (x2 < x && p1[x2] >= accuracya) {
         x2++;  y1 = 0.;               // increase upper limit until x has been reached
      }
      else {
         y1 = p1[x2];
      }
      if (x1 > x2) return 0.;
      if (p2+x2-p > BUFSIZE) FatalError("buffer overrun in function CWalleniusNCHypergeometric::recursive");

      mxo = (m-x2)*omega;
      Nmnx = N-m-nu+x2+1;
      for (xi = x2; xi >= x1; xi--) {  // backwards loop
         d2 = mxo + Nmnx;
         mxo += omega; Nmnx--;
         d1 = mxo + Nmnx;
         dcom = 1. / (d1 * d2);        // save a division by making common divisor
         y  = p1[xi-1]*mxo*d2*dcom + y1*(Nmnx+1)*d1*dcom;
         y1 = p1[xi-1];                // (warning: pointer alias, can't swap instruction order)
         p2[xi] = y;
      }
      p1 = p2;
   }

   if (x < x1 || x > x2) return 0.;
   return p1[x];
}


double CWalleniusNCHypergeometric::binoexpand() {
   // calculate by binomial expansion of integrand
   // only for x < 2 or n-x < 2 (not implemented for higher x because of loss of precision)
   int32_t x1, m1, m2;
   double o;
   if (x > n/2) { // invert
      x1 = n-x; m1 = N-m; m2 = m; o = 1./omega;
   }
   else {
      x1 = x; m1 = m; m2 = N-m; o = omega;
   }
   if (x1 == 0) {
      return exp(FallingFactorial(m2,n) - FallingFactorial(m2+o*m1,n));
   }    
   if (x1 == 1) {
      double d, e, q, q0, q1;
      q = FallingFactorial(m2,n-1);
      e = o*m1+m2;
      q1 = q - FallingFactorial(e,n);
      e -= o;
      q0 = q - FallingFactorial(e,n);
      d = e - (n-1);
      return m1*d*(exp(q0) - exp(q1));
   }

   FatalError("x > 1 not supported by function CWalleniusNCHypergeometric::binoexpand");
   return 0;
}


double CWalleniusNCHypergeometric::laplace() {
   // Laplace's method with narrow integration interval, 
   // using error function residues table, defined in erfres.cpp
   // Note that this function can only be used when the integrand peak is narrow.
   // findpars() must be called before this function.

   const int COLORS = 2;         // number of colors
   const int MAXDEG = 40;        // arraysize, maximum expansion degree
   int degree;                   // max expansion degree
   double accur;                 // stop expansion when terms below this threshold
   double omegai[COLORS] = {omega, 1.}; // weights for each color
   double xi[COLORS] =           // number of each color sampled
       {static_cast<double>(x), static_cast<double>(n-x)};
   double f0;                    // factor outside integral
   double rho[COLORS];           // r*omegai
   double qi;                    // 2^(-rho)
   double qi1;                   // 1-qi
   double qq[COLORS];            // qi / qi1
   double eta[COLORS+1][MAXDEG+1]; // eta coefficients
   double phideri[MAXDEG+1];     // derivatives of phi
   double PSIderi[MAXDEG+1];     // derivatives of PSI
   double * erfresp;             // pointer to table of error function residues

   // variables in asymptotic summation
   static const double sqrt8  = 2.828427124746190098; // sqrt(8)
   double qqpow;                 // qq^j
   double pow2k;                 // 2^k
   double bino;                  // binomial coefficient  
   double vr;                    // 1/v, v = integration interval
   double v2m2;                  // (2*v)^(-2)
   double v2mk1;                 // (2*v)^(-k-1)
   double s;                     // summation term
   double sum;                   // Taylor sum

   int i;                        // loop counter for color
   int j;                        // loop counter for derivative
   int k;                        // loop counter for expansion degree
   int ll;                       // k/2
   int converg = 0;              // number of consequtive terms below accuracy
   int PrecisionIndex;           // index into ErfRes table according to desired precision

   // initialize
   for (k = 0; k <= 2; k++)  phideri[k] = PSIderi[k] = 0;

   // find rho[i], qq[i], first eta coefficients, and zero'th derivative of phi
   for (i = 0; i < COLORS; i++) {
      rho[i] = r * omegai[i];
      if (rho[i] > 40.) {
         qi=0.;  qi1 = 1.;}                 // avoid underflow
      else {
         qi1 = pow2_1(-rho[i], &qi);}       // qi=2^(-rho), qi1=1.-2^(-rho)
      qq[i] = qi / qi1;                     // 2^(-r*omegai)/(1.-2^(-r*omegai))
      // peak = zero'th derivative
      phideri[0] += xi[i] * log1mx(qi, qi1);
      // eta coefficients
      eta[i][0] = 0.;
      eta[i][1] = eta[i][2] = rho[i]*rho[i];
   }

   // r, rd, and w must be calculated by findpars()
   // zero'th derivative
   phideri[0] -= (rd - 1.) * LN2;
   // scaled factor outside integral
   f0 = rd * exp(phideri[0] + lnbico());

   vr = sqrt8 * w;
   phideri[2] = phi2d;

   // get table according to desired precision
   PrecisionIndex = (-FloorLog2((float)accuracy) - ERFRES_B + ERFRES_S - 1) / ERFRES_S;
   if (PrecisionIndex < 0) PrecisionIndex = 0;
   if (PrecisionIndex > ERFRES_N-1) PrecisionIndex = ERFRES_N-1;
   while (w * NumSDev[PrecisionIndex] > 0.3) { 
      // check if integration interval is too wide
      if (PrecisionIndex == 0) {
         FatalError("Laplace method failed. Peak width too high in function CWalleniusNCHypergeometric::laplace");
         break;}
      PrecisionIndex--;                     // reduce precision to keep integration interval narrow
   }
   erfresp = ErfRes[PrecisionIndex];        // choose desired table

   degree = MAXDEG;                         // max expansion degree
   if (degree >= ERFRES_L*2) degree = ERFRES_L*2-2;

   // set up for starting loop at k=3
   v2m2 = 0.25 * vr * vr;                   // (2*v)^(-2)
   PSIderi[0] = 1.;
   pow2k = 8.;
   sum = 0.5 * vr * erfresp[0];
   v2mk1 = 0.5 * vr * v2m2 * v2m2;
   accur = accuracy * sum;

   // summation loop
   for (k = 3; k <= degree; k++) {
      phideri[k] = 0.;

      // loop for all (2) colors
      for (i = 0; i < COLORS; i++) {
         eta[i][k] = 0.;
         // backward loop for all powers
         for (j = k; j > 0; j--) {
            // find coefficients recursively from previous coefficients
            eta[i][j]  =  eta[i][j]*(j*rho[i]-(k-2)) +  eta[i][j-1]*rho[i]*(j-1);
         }
         qqpow = 1.;
         // forward loop for all powers
         for (j=1; j<=k; j++) {
            qqpow *= qq[i];                 // qq^j
            // contribution to derivative
            phideri[k] += xi[i] * eta[i][j] * qqpow;
         }
      }

      // finish calculation of derivatives
      phideri[k] = -pow2k*phideri[k] + 2*(1-k)*phideri[k-1];

      pow2k *= 2.;    // 2^k

      // loop to calculate derivatives of PSI from derivatives of psi.
      // terms # 0, 1, 2, k-2, and k-1 are zero and not included in loop.
      // The j'th derivatives of psi are identical to the derivatives of phi for j>2, and
      // zero for j=1,2. Hence we are using phideri[j] for j>2 here.
      PSIderi[k] = phideri[k];              // this is term # k
      bino = 0.5 * (k-1) * (k-2);           // binomial coefficient for term # 3
      for (j = 3; j < k-2; j++) {           // loop for remaining nonzero terms (if k>5)
         PSIderi[k] += PSIderi[k-j] * phideri[j] * bino;
         bino *= double(k-j)/double(j);
      }
      if ((k & 1) == 0) {                   // only for even k
         ll = k/2;
         s = PSIderi[k] * v2mk1 * erfresp[ll];
         sum += s;

         // check for convergence of Taylor expansion
         if (fabs(s) < accur) converg++; else converg = 0;
         if (converg > 1) break;

         // update recursive expressions
         v2mk1 *= v2m2;
      }
   }
   // multiply by terms outside integral  
   return f0 * sum;
}


double CWalleniusNCHypergeometric::integrate() {
   // Wallenius non-central hypergeometric distribution function
   // calculation by numerical integration with variable-length steps
   // NOTE: findpars() must be called before this function.
   double s;                           // result of integration step
   double sum;                         // integral
   double ta, tb;                      // subinterval for integration step

   lnbico();                           // compute log of binomial coefficients

   // choose method:
   if (w < 0.02 || (w < 0.1 && (x==m || n-x==N-m) && accuracy > 1E-6)) {
      // normal method. Step length determined by peak width w
      double delta, s1;
      s1 = accuracy < 1E-9 ? 0.5 : 1.;
      delta = s1 * w;                       // integration steplength
      ta = 0.5 + 0.5 * delta;
      sum = integrate_step(1.-ta, ta);      // first integration step around center peak
      do {
         tb = ta + delta;
         if (tb > 1.) tb = 1.;
         s  = integrate_step(ta, tb);       // integration step to the right of peak
         s += integrate_step(1.-tb,1.-ta);  // integration step to the left of peak
         sum += s;
         if (s < accuracy * sum) break;     // stop before interval finished if accuracy reached
         ta = tb;
         if (tb > 0.5 + w) delta *= 2.;     // increase step length far from peak
      }
      while (tb < 1.);
   }
   else {
      // difficult situation. Step length determined by inflection points
      double t1, t2, tinf, delta, delta1;
      sum = 0.;
      // do left and right half of integration interval separately:
      for (t1=0., t2=0.5; t1 < 1.; t1+=0.5, t2+=0.5) { 
         // integrate from 0 to 0.5 or from 0.5 to 1
         tinf = search_inflect(t1, t2);     // find inflection point
         delta = tinf - t1; if (delta > t2 - tinf) delta = t2 - tinf; // distance to nearest endpoint
         delta *= 1./7.;                    // 1/7 will give 3 steps to nearest endpoint
         if (delta < 1E-4) delta = 1E-4;
         delta1 = delta;
         // integrate from tinf forwards to t2
         ta = tinf;
         do {
            tb = ta + delta1;
            if (tb > t2 - 0.25*delta1) tb = t2; // last step of this subinterval
            s = integrate_step(ta, tb);         // integration step
            sum += s;
            delta1 *= 2;                        // double steplength
            if (s < sum * 1E-4) delta1 *= 8.;   // large step when s small
            ta = tb;
         }
         while (tb < t2);
         if (tinf) {
            // integrate from tinf backwards to t1
            tb = tinf;
            do {
               ta = tb - delta;
               if (ta < t1 + 0.25*delta) ta = t1; // last step of this subinterval
               s = integrate_step(ta, tb);        // integration step
               sum += s;
               delta *= 2;                        // double steplength
               if (s < sum * 1E-4) delta *= 8.;   // large step when s small
               tb = ta;}
            while (ta > t1);
         }
      }
   }
   return sum * rd;
}


double CWalleniusNCHypergeometric::integrate_step(double ta, double tb) {
   // integration subprocedure used by integrate()
   // makes one integration step from ta to tb using Gauss-Legendre method.
   // result is scaled by multiplication with exp(bico)
   double ab, delta, tau, ltau, y, sum, taur, rdm1;
   int i;

   // define constants for Gauss-Legendre integration with IPOINTS points
#define IPOINTS  8  // number of points in each integration step

#if   IPOINTS == 3
   static const double xval[3]    = {-.774596669241,0,0.774596668241};
   static const double weights[3] = {.5555555555555555,.88888888888888888,.55555555555555};
#elif IPOINTS == 4
   static const double xval[4]    = {-0.861136311594,-0.339981043585,0.339981043585,0.861136311594},
      static const double weights[4] = {0.347854845137,0.652145154863,0.652145154863,0.347854845137};
#elif IPOINTS == 5
   static const double xval[5]    = {-0.906179845939,-0.538469310106,0,0.538469310106,0.906179845939};
   static const double weights[5] = {0.236926885056,0.478628670499,0.568888888889,0.478628670499,0.236926885056};
#elif IPOINTS == 6
   static const double xval[6]    = {-0.932469514203,-0.661209386466,-0.238619186083,0.238619186083,0.661209386466,0.932469514203};
   static const double weights[6] = {0.171324492379,0.360761573048,0.467913934573,0.467913934573,0.360761573048,0.171324492379};
#elif IPOINTS == 8
   static const double xval[8]    = {-0.960289856498,-0.796666477414,-0.525532409916,-0.183434642496,0.183434642496,0.525532409916,0.796666477414,0.960289856498};
   static const double weights[8] = {0.10122853629,0.222381034453,0.313706645878,0.362683783378,0.362683783378,0.313706645878,0.222381034453,0.10122853629};
#elif IPOINTS == 12
   static const double xval[12]   = {-0.981560634247,-0.90411725637,-0.769902674194,-0.587317954287,-0.367831498998,-0.125233408511,0.125233408511,0.367831498998,0.587317954287,0.769902674194,0.90411725637,0.981560634247};
   static const double weights[12]= {0.0471753363866,0.106939325995,0.160078328543,0.203167426723,0.233492536538,0.249147045813,0.249147045813,0.233492536538,0.203167426723,0.160078328543,0.106939325995,0.0471753363866};
#elif IPOINTS == 16
   static const double xval[16]   = {-0.989400934992,-0.944575023073,-0.865631202388,-0.755404408355,-0.617876244403,-0.458016777657,-0.281603550779,-0.0950125098376,0.0950125098376,0.281603550779,0.458016777657,0.617876244403,0.755404408355,0.865631202388,0.944575023073,0.989400934992};
   static const double weights[16]= {0.027152459411,0.0622535239372,0.0951585116838,0.124628971256,0.149595988817,0.169156519395,0.182603415045,0.189450610455,0.189450610455,0.182603415045,0.169156519395,0.149595988817,0.124628971256,0.0951585116838,0.0622535239372,0.027152459411};
#else
#error // IPOINTS must be a value for which the tables are defined
#endif

   delta = 0.5 * (tb - ta);
   ab = 0.5 * (ta + tb);
   rdm1 = rd - 1.;
   sum = 0;

   for (i = 0; i < IPOINTS; i++) {
      tau = ab + delta * xval[i];
      ltau = log(tau);
      taur = r * ltau;
      // possible loss of precision due to subtraction here:
      y = log1pow(taur*omega,x) + log1pow(taur,n-x) + rdm1*ltau + bico;
      if (y > -50.) sum += weights[i] * exp(y);
   }
   return delta * sum;
}


double CWalleniusNCHypergeometric::search_inflect(double t_from, double t_to) {
   // search for an inflection point of the integrand PHI(t) in the interval
   // t_from < t < t_to
   const int COLORS = 2;                // number of colors
   double t, t1;                        // independent variable
   double rho[COLORS];                  // r*omega[i]
   double q;                            // t^rho[i] / (1-t^rho[i])
   double q1;                           // 1-t^rho[i]
   double xx[COLORS];                   // x[i]
   double zeta[COLORS][4][4];           // zeta[i,j,k] coefficients
   double phi[4];                       // derivatives of phi(t) = log PHI(t)
   double Z2;                           // PHI''(t)/PHI(t)
   double Zd;                           // derivative in Newton Raphson iteration
   double rdm1;                         // r * d - 1
   double tr;                           // 1/t
   double log2t;                        // log2(t)
   double method;                       // 0 for z2'(t) method, 1 for z3(t) method
   int i;                               // color
   int iter;                            // count iterations

   rdm1 = rd - 1.;
   if (t_from == 0 && rdm1 <= 1.) return 0.;//no inflection point
   rho[0] = r*omega;  rho[1] = r;
   xx[0] = x;  xx[1] = n - x;
   t = 0.5 * (t_from + t_to);
   for (i = 0; i < COLORS; i++) {           // calculate zeta coefficients
      zeta[i][1][1] = rho[i];
      zeta[i][1][2] = rho[i] * (rho[i] - 1.);
      zeta[i][2][2] = rho[i] * rho[i];
      zeta[i][1][3] = zeta[i][1][2] * (rho[i] - 2.);
      zeta[i][2][3] = zeta[i][1][2] * rho[i] * 3.;
      zeta[i][3][3] = zeta[i][2][2] * rho[i] * 2.;
   }
   iter = 0;

   do {
      t1 = t;
      tr = 1. / t;
      log2t = log(t)*(1./LN2);
      phi[1] = phi[2] = phi[3] = 0.;
      for (i=0; i<COLORS; i++) {            // calculate first 3 derivatives of phi(t)
         q1 = pow2_1(rho[i]*log2t,&q);
         q /= q1;
         phi[1] -= xx[i] * zeta[i][1][1] * q;
         phi[2] -= xx[i] * q * (zeta[i][1][2] + q * zeta[i][2][2]);
         phi[3] -= xx[i] * q * (zeta[i][1][3] + q * (zeta[i][2][3] + q * zeta[i][3][3]));
      }
      phi[1] += rdm1;
      phi[2] -= rdm1;
      phi[3] += 2. * rdm1;
      phi[1] *= tr;
      phi[2] *= tr * tr;
      phi[3] *= tr * tr * tr;
      method = (iter & 2) >> 1;        // alternate between the two methods
      Z2 = phi[1]*phi[1] + phi[2];
      Zd = method*phi[1]*phi[1]*phi[1] + (2.+method)*phi[1]*phi[2] + phi[3];

      if (t < 0.5) {
         if (Z2 > 0) {
            t_from = t;
         }
         else {
            t_to = t;
         }
         if (Zd >= 0) { 
            // use binary search if Newton-Raphson iteration makes problems
            t = (t_from ? 0.5 : 0.2) * (t_from + t_to);
         }
         else {
            // Newton-Raphson iteration
            t -= Z2 / Zd;
         }
      }
      else {
         if (Z2 < 0) {
            t_from = t;
         }
         else {
            t_to = t;
         }
         if (Zd <= 0) {
            // use binary search if Newton-Raphson iteration makes problems
            t = 0.5 * (t_from + t_to);
         }
         else {
            // Newton-Raphson iteration
            t -= Z2 / Zd;
         }
      }
      if (t >= t_to) t = (t1 + t_to) * 0.5;
      if (t <= t_from) t = (t1 + t_from) * 0.5;
      if (++iter > 20) FatalError("Search for inflection point failed in function CWalleniusNCHypergeometric::search_inflect");
   }
   while (fabs(t - t1) > 1E-5);
   return t;
}


double CWalleniusNCHypergeometric::probability(int32_t x_) {
   // calculate probability function. choosing best method
   x = x_;
   if (x < xmin || x > xmax) return 0.;
   if (xmin == xmax) return 1.;

   if (omega == 1.) { // hypergeometric
      return exp(lnbico() + LnFac(n) + LnFac(N-n) - LnFac(N));
   }

   if (omega == 0.) {
      if (n > N-m) FatalError("Not enough items with nonzero weight in CWalleniusNCHypergeometric::probability");
      return x == 0;
   }

   int32_t x2 = n - x;
   int32_t x0 = x < x2 ? x : x2;
   int em = (x == m || x2 == N-m);

   if (x0 == 0 && n > 500) {
      return binoexpand();
   }

   if (double(n)*x0 < 1000 || (double(n)*x0 < 10000 && (N > 1000.*n || em))) {
      return recursive();
   }

   if (x0 <= 1 && N-n <= 1) {
      return binoexpand();
   }

   findpars();

   if (w < 0.04 && E < 10 && (!em || w > 0.004)) {
      return laplace();
   }

   return integrate();
}


int32_t CWalleniusNCHypergeometric::MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff) {
   // Makes a table of Wallenius noncentral hypergeometric probabilities 
   // table must point to an array of length MaxLength. 
   // The function returns 1 if table is long enough. Otherwise it fills
   // the table with as many correct values as possible and returns 0.
   // The tails are cut off where the values are < cutoff, so that 
   // *xfirst may be > xmin and *xlast may be < xmax.
   // The value of cutoff will be 0.01 * accuracy if not specified.
   // The first and last x value represented in the table are returned in 
   // *xfirst and *xlast. The resulting probability values are returned in 
   // the first (*xfirst - *xlast + 1) positions of table. Any unused part
   // of table may be overwritten with garbage.
   //
   // The function will return the following information when MaxLength = 0:
   // The return value is the desired length of table.
   // *xfirst is 1 if it will be more efficient to call MakeTable than to call
   // probability repeatedly, even if only some of the table values are needed.
   // *xfirst is 0 if it is more efficient to call probability repeatedly.

   double * p1, * p2;                  // offset into p
   double mxo;                         // (m-x)*omega
   double Nmnx;                        // N-m-nu+x
   double y, y1;                       // probability. Save old p[x] before it is overwritten
   double d1, d2, dcom;                // divisors in probability formula
   double area;                        // estimate of area needed for recursion method
   int32_t xi, nu;                     // xi, nu = recursion values of x, n
   int32_t x1, x2;                     // lowest and highest x or xi
   int32_t i1, i2;                     // index into table
   int32_t UseTable;                   // 1 if table method used
   int32_t LengthNeeded;               // Necessary table length

   // special cases
   if (n == 0 || m == 0) {x1 = 0; goto DETERMINISTIC;}
   if (n == N)           {x1 = m; goto DETERMINISTIC;}
   if (m == N)           {x1 = n; goto DETERMINISTIC;}
   if (omega <= 0.) {
      if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
      x1 = 0;
      DETERMINISTIC:
      if (MaxLength == 0) {
         if (xfirst) *xfirst = 1;
         return 1;
      }
      *xfirst = *xlast = x1;
      *table = 1.;
      return 1;
   }

   if (cutoff <= 0. || cutoff > 0.1) cutoff = 0.01 * accuracy;

   LengthNeeded = N - m;               // m2
   if (m < LengthNeeded) LengthNeeded = m;
   if (n < LengthNeeded) LengthNeeded = n; // LengthNeeded = min(m1,m2,n)
   area = double(n)*LengthNeeded;      // Estimate calculation time for table method
   UseTable = area < 5000. || (area < 10000. && N > 1000. * n);

   if (MaxLength <= 0) {
      // Return UseTable and LengthNeeded
      if (xfirst) *xfirst = UseTable;
      i1 = LengthNeeded + 2;           // Necessary table length
      if (!UseTable && i1 > 200) {
         // Calculate necessary table length from standard deviation
         double sd = sqrt(variance()); // calculate approximate standard deviation
         // estimate number of standard deviations to include from normal distribution
         i2 = (int32_t)(NumSD(accuracy) * sd + 0.5);
         if (i1 > i2) i1 = i2;
      }
      return i1;
   }

   if (UseTable && MaxLength > LengthNeeded) {
      // use recursion table method
      p1 = p2 = table + 1;             // make space for p1[-1]
      p1[-1] = 0.;  p1[0] = 1.;        // initialize for recursion
      x1 = x2 = 0;
      for (nu = 1; nu <= n; nu++) {
         if (n - nu < xmin - x1 || p1[x1] < cutoff) {
            x1++;                      // increase lower limit when breakpoint passed or probability negligible
            p2--;                      // compensate buffer offset in order to reduce storage space
         }
         if (x2 < xmax && p1[x2] >= cutoff) {
            x2++;  y1 = 0.;            // increase upper limit until x has been reached
         }
         else {
            y1 = p1[x2];
         }
         if (p2 - table + x2 >= MaxLength || x1 > x2) {
            goto ONE_BY_ONE;           // Error: table length exceeded. Use other method
         }

         mxo = (m-x2)*omega;
         Nmnx = N-m-nu+x2+1;
         for (xi = x2; xi >= x1; xi--) { // backwards loop
            d2 = mxo + Nmnx;
            mxo += omega; Nmnx--;
            d1 = mxo + Nmnx;
            dcom = 1. / (d1 * d2);     // save a division by making common divisor
            y  = p1[xi-1]*mxo*d2*dcom + y1*(Nmnx+1)*d1*dcom;
            y1 = p1[xi-1];             // (warning: pointer alias, can't swap instruction order)
            p2[xi] = y;
         }
         p1 = p2;
      }

      // return results
      i1 = i2 = x2 - x1 + 1;              // desired table length
      if (i2 > MaxLength) i2 = MaxLength; // limit table length
      *xfirst = x1;  *xlast = x1 + i2 - 1;
      if (i2 > 0) memmove(table, table+1, i2*sizeof(table[0]));// copy to start of table
      return i1 == i2;                    // true if table size not reduced
   }

   else {
      // Recursion method would take too much time
      // Calculate values one by one
      ONE_BY_ONE:

      // Start to fill table from the end and down. start with x = floor(mean)
      x2 = (int32_t)mean();
      x1 = x2 + 1;  i1 = MaxLength;
      while (x1 > xmin) {              // loop for left tail
         x1--;  i1--;
         y = probability(x1);
         table[i1] = y;
         if (y < cutoff) break;
         if (i1 == 0) break;
      }
      *xfirst = x1;
      i2 = x2 - x1 + 1; 
      if (i1 > 0 && i2 > 0) { // move numbers down to beginning of table
         memmove(table, table+i1, i2*sizeof(table[0]));
      }
      // Fill rest of table from mean and up
      i2--;
      while (x2 < xmax) {              // loop for right tail
         if (i2 == MaxLength-1) {
            *xlast = x2; return 0;     // table full
         }
         x2++;  i2++;
         y = probability(x2);
         table[i2] = y;
         if (y < cutoff) break;
      }
      *xlast = x2;
      return 1;
   }
}


/***********************************************************************
calculation methods in class CMultiWalleniusNCHypergeometric
***********************************************************************/

CMultiWalleniusNCHypergeometric::CMultiWalleniusNCHypergeometric(int32_t n_, int32_t * m_, double * odds_, int colors_, double accuracy_) {
   // constructor
   accuracy = accuracy_;
   SetParameters(n_, m_, odds_, colors_);
}


void CMultiWalleniusNCHypergeometric::SetParameters(int32_t n_, int32_t * m_, double * odds_, int colors_) {
   // change parameters
   int32_t N1;
   int i;
   n = n_;  m = m_;  omega = odds_;  colors = colors_;
   r = 1.;
   for (N = N1 = 0, i = 0; i < colors; i++) {
      if (m[i] < 0 || omega[i] < 0) FatalError("Parameter negative in constructor for CMultiWalleniusNCHypergeometric");
      N += m[i];
      if (omega[i]) N1 += m[i];
   }
   if (N < n) FatalError("Not enough items in constructor for CMultiWalleniusNCHypergeometric");
   if (N1< n) FatalError("Not enough items with nonzero weight in constructor for CMultiWalleniusNCHypergeometric");
}


void CMultiWalleniusNCHypergeometric::mean(double * mu) {
   // calculate approximate mean of multivariate Wallenius noncentral hypergeometric 
   // distribution. Result is returned in mu[0..colors-1]
   double omeg[MAXCOLORS];              // scaled weights
   double omr;                          // reciprocal mean weight
   double t, t1;                        // independent variable in iteration
   double To, To1;                      // exp(t*omega[i]), 1-exp(t*omega[i])
   double H;                            // function to find root of
   double HD;                           // derivative of H
   double dummy;                        // unused return
   int i;                               // color index
   int iter;                            // number of iterations

   if (n == 0) {
      // needs special case
      for (i = 0; i < colors; i++) {
         mu[i] = 0.;
      }
      return;
   }

   // calculate mean weight
   for (omr=0., i=0; i < colors; i++) omr += omega[i] * m[i];
   omr = N / omr;
   // scale weights to make mean = 1
   for (i = 0; i < colors; i++) omeg[i] = omega[i] * omr;
   // Newton Raphson iteration
   iter = 0;  t = -1.;                  // first guess
   do {
      t1 = t;
      H = HD = 0.;
      // calculate H and HD
      for (i = 0; i < colors; i++) {
         if (omeg[i] != 0.) {
            To1 = pow2_1(t * (1./LN2) * omeg[i], &To);
            H += m[i] * To1;
            HD -= m[i] * omeg[i] * To;
         }
      }
      t -= (H-n) / HD;
      if (t >= 0) t = 0.5 * t1;
      if (++iter > 20) {
         FatalError("Search for mean failed in function CMultiWalleniusNCHypergeometric::mean");
      }
   }
   while (fabs(H - n) > 1E-3);
   // finished iteration. Get all mu[i]
   for (i = 0; i < colors; i++) {
      if (omeg[i] != 0.) {
         To1 = pow2_1(t * (1./LN2) * omeg[i], &dummy);
         mu[i] = m[i] * To1;
      }
      else {
         mu[i] = 0.;
      }
   }
}
/*
void CMultiWalleniusNCHypergeometric::variance(double * var, double * mean_) {
   // calculates approximate variance and mean of multivariate 
   // Wallenius' noncentral hypergeometric distribution 
   // (accuracy is not too good).
   // Variance is returned in variance[0..colors-1].
   // Mean is returned in mean_[0..colors-1] if not NULL.
   // The calculation is reasonably fast.
   double r1, r2;
   double mu[MAXCOLORS];
   int i;

   // Store mean in array mu if mean_ is NULL
   if (mean_ == 0) mean_ = mu;

   // Calculate mean
   mean(mean_);

   // Calculate variance
   for (i = 0; i < colors; i++) {
      r1 = mean_[i] * (m[i]-mean_[i]);
      r2 = (n-mean_[i])*(mean_[i]+N-n-m[i]);
      if (r1 <= 0. || r2 <= 0.) {
         var[i] = 0.;
      }
      else {
         var[i] = N*r1*r2/((N-1)*(m[i]*r2+(N-m[i])*r1));
      }
   }
}
*/

// implementations of different calculation methods
double CMultiWalleniusNCHypergeometric::binoexpand(void) {
   // binomial expansion of integrand
   // only implemented for x[i] = 0 for all but one i
   int i, j, k;
   double W = 0.;                       // total weight
   for (i=j=k=0; i<colors; i++) {
      W += omega[i] * m[i];
      if (x[i]) {
         j=i; k++;                      // find the nonzero x[i]
      }
   }
   if (k > 1) FatalError("More than one x[i] nonzero in CMultiWalleniusNCHypergeometric::binoexpand");
   return exp(FallingFactorial(m[j],n) - FallingFactorial(W/omega[j],n));
}


double CMultiWalleniusNCHypergeometric::lnbico(void) {
   // natural log of binomial coefficients
   bico = 0.;
   int i;
   for (i=0; i<colors; i++) {
      if (x[i] < m[i] && omega[i]) {
         bico += LnFac(m[i]) - LnFac(x[i]) - LnFac(m[i]-x[i]);
      }
   }
   return bico;
}


void CMultiWalleniusNCHypergeometric::findpars(void) {
   // calculate r, w, E
   // calculate d, E, r, w

   // find r to center peak of integrand at 0.5
   double dd;                           // scaled d
   double dr;                           // 1/d

   double z, zd, rr, lastr, rrc, rt, r2, r21, a, b, ro, k1, dummy;
   double omax;                         // highest omega
   double omaxr;                        // 1/omax
   double omeg[MAXCOLORS];              // scaled weights
   int i, j = 0;

   // find highest omega
   for (omax=0., i=0; i < colors; i++) {
      if (omega[i] > omax) omax = omega[i];
   }
   omaxr = 1. / omax;
   dd = E = 0.;
   for (i = 0; i < colors; i++) {
      // scale weights to make max = 1
      omeg[i] = omega[i] * omaxr;
      // calculate d and E
      dd += omeg[i] * (m[i]-x[i]);
      E  += omeg[i] * m[i];
   }
   dr = 1. / dd;
   E *= dr;
   rr = r * omax;
   if (rr <= dr) rr = 1.2 * dr;  // initial guess
   // Newton-Raphson iteration to find r
   do {
      lastr = rr;
      rrc = 1. / rr;
      z = dd - rrc;                    // z(r)
      zd = rrc * rrc;                  // z'(r)
      for (i=0; i<colors; i++) {
         rt = rr * omeg[i];
         if (rt < 100. && rt > 0.) {   // avoid overflow and division by 0
            r21 = pow2_1(rt, &r2);     // r2=2^r, r21=1.-2^r
            a = omeg[i] / r21;         // omegai/(1.-2^r)
            b = x[i] * a;              // x*omegai/(1.-2^r)
            z  += b;
            zd += b * a * r2 * LN2;
         }
      }
      if (zd == 0) FatalError("can't find r in function CMultiWalleniusNCHypergeometric::findpars");
      rr -= z / zd;                    // next r
      if (rr <= dr) rr = lastr * 0.125 + dr * 0.875;
      if (++j == 70) FatalError("convergence problem searching for r in function CMultiWalleniusNCHypergeometric::findpars");
   }
   while (fabs(rr-lastr) > rr * 1.E-5);
   rd = rr * dd;
   r = rr * omaxr;

   // find peak width
   phi2d = 0.;
   for (i=0; i<colors; i++) {
      ro = rr * omeg[i];
      if (ro < 300 && ro > 0.) {       // avoid overflow and division by 0
         k1 = pow2_1(ro, &dummy);
         k1 = -1. / k1;
         k1 = omeg[i] * omeg[i] * (k1 + k1*k1);
      }
      else k1 = 0.;
      phi2d += x[i] * k1;
   }
   phi2d *= -4. * rr * rr;
   if (phi2d > 0.) FatalError("peak width undefined in function CMultiWalleniusNCHypergeometric::findpars");
   wr = sqrt(-phi2d);  w = 1. / wr;
}


double CMultiWalleniusNCHypergeometric::laplace(void) {
   // Laplace's method with narrow integration interval, 
   // using error function residues table, defined in erfres.cpp
   // Note that this function can only be used when the integrand peak is narrow.
   // findpars() must be called before this function.

   const int MAXDEG = 40;              // arraysize
   int degree;                         // max expansion degree
   double accur;                       // stop expansion when terms below this threshold
   double f0;                          // factor outside integral
   double rho[MAXCOLORS];              // r*omegai
   double qi;                          // 2^(-rho)
   double qi1;                         // 1-qi
   double qq[MAXCOLORS];               // qi / qi1
   double eta[MAXCOLORS+1][MAXDEG+1];  // eta coefficients
   double phideri[MAXDEG+1];           // derivatives of phi
   double PSIderi[MAXDEG+1];           // derivatives of PSI
   double * erfresp;                   // pointer to table of error function residues

   // variables in asymptotic summation
   static const double sqrt8  = 2.828427124746190098; // sqrt(8)
   double qqpow;                       // qq^j
   double pow2k;                       // 2^k
   double bino;                        // binomial coefficient  
   double vr;                          // 1/v, v = integration interval
   double v2m2;                        // (2*v)^(-2)
   double v2mk1;                       // (2*v)^(-k-1)
   double s;                           // summation term
   double sum;                         // Taylor sum

   int i;                              // loop counter for color
   int j;                              // loop counter for derivative
   int k;                              // loop counter for expansion degree
   int ll;                             // k/2
   int converg = 0;                    // number of consequtive terms below accuracy
   int PrecisionIndex;                 // index into ErfRes table according to desired precision

   // initialize
   for (k = 0; k <= 2; k++)  phideri[k] = PSIderi[k] = 0;

   // find rho[i], qq[i], first eta coefficients, and zero'th derivative of phi
   for (i = 0; i < colors; i++) {
      rho[i] = r * omega[i];
      if (rho[i] == 0.) continue;
      if (rho[i] > 40.) {
         qi=0.;  qi1 = 1.;             // avoid underflow
      }
      else {
         qi1 = pow2_1(-rho[i], &qi);   // qi=2^(-rho), qi1=1.-2^(-rho)
      }
      qq[i] = qi / qi1;                // 2^(-r*omegai)/(1.-2^(-r*omegai))
      // peak = zero'th derivative
      phideri[0] += x[i] * log1mx(qi, qi1);
      // eta coefficients
      eta[i][0] = 0.;
      eta[i][1] = eta[i][2] = rho[i]*rho[i];
   }

   // d, r, and w must be calculated by findpars()
   // zero'th derivative
   phideri[0] -= (rd - 1.) * LN2;
   // scaled factor outside integral
   f0 = rd * exp(phideri[0] + lnbico());
   // calculate narrowed integration interval
   vr = sqrt8 * w;
   phideri[2] = phi2d;

   // get table according to desired precision
   PrecisionIndex = (-FloorLog2((float)accuracy) - ERFRES_B + ERFRES_S - 1) / ERFRES_S;
   if (PrecisionIndex < 0) PrecisionIndex = 0;
   if (PrecisionIndex > ERFRES_N-1) PrecisionIndex = ERFRES_N-1;
   while (w * NumSDev[PrecisionIndex] > 0.3) { 
      // check if integration interval is too wide
      if (PrecisionIndex == 0) {
         FatalError("Laplace method failed. Peak width too high in function CWalleniusNCHypergeometric::laplace");
         break;
      }
      PrecisionIndex--;                // reduce precision to keep integration interval narrow
   }
   erfresp = ErfRes[PrecisionIndex];   // choose desired table

   degree = MAXDEG;                    // max expansion degree
   if (degree >= ERFRES_L*2) degree = ERFRES_L*2-2;

   // set up for starting loop at k=3
   v2m2 = 0.25 * vr * vr;              // (2*v)^(-2)
   PSIderi[0] = 1.;
   pow2k = 8.;
   sum = 0.5 * vr * erfresp[0];
   v2mk1 = 0.5 * vr * v2m2 * v2m2;
   accur = accuracy * sum;

   // summation loop
   for (k = 3; k <= degree; k++) {
      phideri[k] = 0.;

      // loop for all colors
      for (i = 0; i < colors; i++) {
         if (rho[i] == 0.) continue;
         eta[i][k] = 0.;
         // backward loop for all powers
         for (j = k; j > 0; j--) {
            // find coefficients recursively from previous coefficients
            eta[i][j]  =  eta[i][j]*(j*rho[i]-(k-2)) +  eta[i][j-1]*rho[i]*(j-1);
         }
         qqpow = 1.;
         // forward loop for all powers
         for (j = 1; j <= k; j++) {
            qqpow *= qq[i];   // qq^j
            // contribution to derivative
            phideri[k] += x[i] * eta[i][j] * qqpow;
         }
      }

      // finish calculation of derivatives
      phideri[k] = -pow2k * phideri[k] + 2*(1-k)*phideri[k-1];

      pow2k *= 2.;                     // 2^k

      // loop to calculate derivatives of PSI from derivatives of psi.
      // terms # 0, 1, 2, k-2, and k-1 are zero and not included in loop.
      // The j'th derivatives of psi are identical to the derivatives of phi for j>2, and
      // zero for j=1,2. Hence we are using phideri[j] for j>2 here.
      PSIderi[k] = phideri[k];         // this is term # k
      bino = 0.5 * (k-1) * (k-2);      // binomial coefficient for term # 3
      for (j=3; j < k-2; j++) { // loop for remaining nonzero terms (if k>5)
         PSIderi[k] += PSIderi[k-j] * phideri[j] * bino;
         bino *= double(k-j)/double(j);
      }

      if ((k & 1) == 0) { // only for even k
         ll = k/2;
         s = PSIderi[k] * v2mk1 * erfresp[ll];
         sum += s;

         // check for convergence of Taylor expansion
         if (fabs(s) < accur) converg++; else converg = 0;
         if (converg > 1) break;

         // update recursive expressions
         v2mk1 *= v2m2;
      }
   }

   // multiply by terms outside integral  
   return f0 * sum;
}


double CMultiWalleniusNCHypergeometric::integrate(void) {
   // Wallenius non-central hypergeometric distribution function
   // calculation by numerical integration with variable-length steps
   // NOTE: findpars() must be called before this function.
   double s;                           // result of integration step
   double sum;                         // integral
   double ta, tb;                      // subinterval for integration step

   lnbico();                           // compute log of binomial coefficients

   // choose method:
   if (w < 0.02) {
      // normal method. Step length determined by peak width w
      double delta, s1;
      s1 = accuracy < 1E-9 ? 0.5 : 1.;
      delta = s1 * w;                            // integration steplength
      ta = 0.5 + 0.5 * delta;
      sum = integrate_step(1.-ta, ta);           // first integration step around center peak
      do {
         tb = ta + delta;
         if (tb > 1.) tb = 1.;
         s  = integrate_step(ta, tb);            // integration step to the right of peak
         s += integrate_step(1.-tb,1.-ta);       // integration step to the left of peak
         sum += s;
         if (s < accuracy * sum) break;          // stop before interval finished if accuracy reached
         ta = tb;
         if (tb > 0.5 + w) delta *= 2.;          // increase step length far from peak
      }
      while (tb < 1.);
   }

   else {
      // difficult situation. Step length determined by inflection points
      double t1, t2, tinf, delta, delta1;
      sum = 0.;
      // do left and right half of integration interval separately:
      for (t1=0., t2=0.5; t1 < 1.; t1+=0.5, t2+=0.5) { 
         // integrate from 0 to 0.5 or from 0.5 to 1
         tinf = search_inflect(t1, t2);          // find inflection point
         delta = tinf - t1; if (delta > t2 - tinf) delta = t2 - tinf; // distance to nearest endpoint
         delta *= 1./7.;                         // 1/7 will give 3 steps to nearest endpoint
         if (delta < 1E-4) delta = 1E-4;
         delta1 = delta;
         // integrate from tinf forwards to t2
         ta = tinf;
         do {
            tb = ta + delta1;
            if (tb > t2 - 0.25*delta1) tb = t2;  // last step of this subinterval
            s = integrate_step(ta, tb);          // integration step
            sum += s;
            delta1 *= 2;                         // double steplength
            if (s < sum * 1E-4) delta1 *= 8.;    // large step when s small
            ta = tb;
         }
         while (tb < t2);
         if (tinf) {
            // integrate from tinf backwards to t1
            tb = tinf;
            do {
               ta = tb - delta;
               if (ta < t1 + 0.25*delta) ta = t1; // last step of this subinterval
               s = integrate_step(ta, tb);        // integration step
               sum += s;
               delta *= 2;                       // double steplength
               if (s < sum * 1E-4) delta *= 8.;  // large step when s small
               tb = ta;
            }
            while (ta > t1);
         }
      }
   }
   return sum * rd;
}


double CMultiWalleniusNCHypergeometric::integrate_step(double ta, double tb) {
   // integration subprocedure used by integrate()
   // makes one integration step from ta to tb using Gauss-Legendre method.
   // result is scaled by multiplication with exp(bico)
   double ab, delta, tau, ltau, y, sum, taur, rdm1;
   int i, j;

   // define constants for Gauss-Legendre integration with IPOINTS points
#define IPOINTS  8  // number of points in each integration step

#if   IPOINTS == 3
   static const double xval[3]    = {-.774596669241,0,0.774596668241};
   static const double weights[3] = {.5555555555555555,.88888888888888888,.55555555555555};
#elif IPOINTS == 4
   static const double xval[4]    = {-0.861136311594,-0.339981043585,0.339981043585,0.861136311594},
      static const double weights[4] = {0.347854845137,0.652145154863,0.652145154863,0.347854845137};
#elif IPOINTS == 5
   static const double xval[5]    = {-0.906179845939,-0.538469310106,0,0.538469310106,0.906179845939};
   static const double weights[5] = {0.236926885056,0.478628670499,0.568888888889,0.478628670499,0.236926885056};
#elif IPOINTS == 6
   static const double xval[6]    = {-0.932469514203,-0.661209386466,-0.238619186083,0.238619186083,0.661209386466,0.932469514203};
   static const double weights[6] = {0.171324492379,0.360761573048,0.467913934573,0.467913934573,0.360761573048,0.171324492379};
#elif IPOINTS == 8
   static const double xval[8]    = {-0.960289856498,-0.796666477414,-0.525532409916,-0.183434642496,0.183434642496,0.525532409916,0.796666477414,0.960289856498};
   static const double weights[8] = {0.10122853629,0.222381034453,0.313706645878,0.362683783378,0.362683783378,0.313706645878,0.222381034453,0.10122853629};
#elif IPOINTS == 12
   static const double xval[12]   = {-0.981560634247,-0.90411725637,-0.769902674194,-0.587317954287,-0.367831498998,-0.125233408511,0.125233408511,0.367831498998,0.587317954287,0.769902674194,0.90411725637,0.981560634247};
   static const double weights[12]= {0.0471753363866,0.106939325995,0.160078328543,0.203167426723,0.233492536538,0.249147045813,0.249147045813,0.233492536538,0.203167426723,0.160078328543,0.106939325995,0.0471753363866};
#elif IPOINTS == 16
   static const double xval[16]   = {-0.989400934992,-0.944575023073,-0.865631202388,-0.755404408355,-0.617876244403,-0.458016777657,-0.281603550779,-0.0950125098376,0.0950125098376,0.281603550779,0.458016777657,0.617876244403,0.755404408355,0.865631202388,0.944575023073,0.989400934992};
   static const double weights[16]= {0.027152459411,0.0622535239372,0.0951585116838,0.124628971256,0.149595988817,0.169156519395,0.182603415045,0.189450610455,0.189450610455,0.182603415045,0.169156519395,0.149595988817,0.124628971256,0.0951585116838,0.0622535239372,0.027152459411};
#else
#error // IPOINTS must be a value for which the tables are defined
#endif

   delta = 0.5 * (tb - ta);
   ab = 0.5 * (ta + tb);
   rdm1 = rd - 1.;
   sum = 0;

   for (j = 0; j < IPOINTS; j++) {
      tau = ab + delta * xval[j];
      ltau = log(tau);
      taur = r * ltau;
      y = 0.;
      for (i = 0; i < colors; i++) {
         // possible loss of precision due to subtraction here:
         if (omega[i]) {
            y += log1pow(taur*omega[i],x[i]);   // ln((1-e^taur*omegai)^xi)
         }
      }
      y += rdm1*ltau + bico;
      if (y > -50.) sum += weights[j] * exp(y);
   }
   return delta * sum;
}


double CMultiWalleniusNCHypergeometric::search_inflect(double t_from, double t_to) {
   // search for an inflection point of the integrand PHI(t) in the interval
   // t_from < t < t_to
   double t, t1;                       // independent variable
   double rho[MAXCOLORS];              // r*omega[i]
   double q;                           // t^rho[i] / (1-t^rho[i])
   double q1;                          // 1-t^rho[i]
   double zeta[MAXCOLORS][4][4];       // zeta[i,j,k] coefficients
   double phi[4];                      // derivatives of phi(t) = log PHI(t)
   double Z2;                          // PHI''(t)/PHI(t)
   double Zd;                          // derivative in Newton Raphson iteration
   double rdm1;                        // r * d - 1
   double tr;                          // 1/t
   double log2t;                       // log2(t)
   double method;                      // 0 for z2'(t) method, 1 for z3(t) method
   int i;                              // color
   int iter;                           // count iterations

   rdm1 = rd - 1.;
   if (t_from == 0 && rdm1 <= 1.) return 0.;     //no inflection point
   t = 0.5 * (t_from + t_to);
   for (i = 0; i < colors; i++) {                // calculate zeta coefficients
      rho[i] = r * omega[i];
      zeta[i][1][1] = rho[i];
      zeta[i][1][2] = rho[i] * (rho[i] - 1.);
      zeta[i][2][2] = rho[i] * rho[i];
      zeta[i][1][3] = zeta[i][1][2] * (rho[i] - 2.);
      zeta[i][2][3] = zeta[i][1][2] * rho[i] * 3.;
      zeta[i][3][3] = zeta[i][2][2] * rho[i] * 2.;
   }
   iter = 0;

   do {
      t1 = t;
      tr = 1. / t;
      log2t = log(t)*(1./LN2);
      phi[1] = phi[2] = phi[3] = 0.;
      for (i=0; i<colors; i++) {                 // calculate first 3 derivatives of phi(t)
         if (rho[i] == 0.) continue;
         q1 = pow2_1(rho[i]*log2t,&q);
         q /= q1;
         phi[1] -= x[i] * zeta[i][1][1] * q;
         phi[2] -= x[i] * q * (zeta[i][1][2] + q * zeta[i][2][2]);
         phi[3] -= x[i] * q * (zeta[i][1][3] + q * (zeta[i][2][3] + q * zeta[i][3][3]));
      }
      phi[1] += rdm1;
      phi[2] -= rdm1;
      phi[3] += 2. * rdm1;
      phi[1] *= tr;
      phi[2] *= tr * tr;
      phi[3] *= tr * tr * tr;
      method = (iter & 2) >> 1;                  // alternate between the two methods
      Z2 = phi[1]*phi[1] + phi[2];
      Zd = method*phi[1]*phi[1]*phi[1] + (2.+method)*phi[1]*phi[2] + phi[3];

      if (t < 0.5) {
         if (Z2 > 0) {
            t_from = t;
         }
         else {
            t_to = t;
         }
         if (Zd >= 0) { 
            // use binary search if Newton-Raphson iteration makes problems
            t = (t_from ? 0.5 : 0.2) * (t_from + t_to);
         }
         else {
            // Newton-Raphson iteration
            t -= Z2 / Zd;
         }
      }
      else {
         if (Z2 < 0) {
            t_from = t;
         }
         else {
            t_to = t;
         }
         if (Zd <= 0) {
            // use binary search if Newton-Raphson iteration makes problems
            t = 0.5 * (t_from + t_to);
         }
         else {
            // Newton-Raphson iteration
            t -= Z2 / Zd;
         }
      }
      if (t >= t_to) t = (t1 + t_to) * 0.5;
      if (t <= t_from) t = (t1 + t_from) * 0.5;
      if (++iter > 20) FatalError("Search for inflection point failed in function CMultiWalleniusNCHypergeometric::search_inflect");
   }
   while (fabs(t - t1) > 1E-5);
   return t;
}


double CMultiWalleniusNCHypergeometric::probability(int32_t * x_) {
   // calculate probability function. choosing best method
   int i, j, em;
   int central;
   int32_t xsum;
   x = x_;

   for (xsum = i = 0; i < colors; i++)  xsum += x[i];
   if (xsum != n) {
      FatalError("sum of x values not equal to n in function CMultiWalleniusNCHypergeometric::probability");
   }

   if (colors < 3) { 
      if (colors <= 0) return 1.;
      if (colors == 1) return x[0] == m[0];
      // colors = 2
      if (omega[1] == 0.) return x[0] == m[0];
      return CWalleniusNCHypergeometric(n,m[0],N,omega[0]/omega[1],accuracy).probability(x[0]);
   }

   central = 1;
   for (i = j = em = 0; i < colors; i++) {
      if (x[i] > m[i] || x[i] < 0 || x[i] < n - N + m[i]) return 0.;
      if (x[i] > 0) j++;
      if (omega[i] == 0. && x[i]) return 0.;
      if (x[i] == m[i] || omega[i] == 0.) em++;
      if (i > 0 && omega[i] != omega[i-1]) central = 0;
   }

   if (n == 0 || em == colors) return 1.;

   if (central) {
      // All omega's are equal. 
      // This is multivariate central hypergeometric distribution
      int32_t sx = n,  sm = N;
      double p = 1.;
      for (i = 0; i < colors - 1; i++) {
         // Use univariate hypergeometric (usedcolors-1) times
         p *= CWalleniusNCHypergeometric(sx, m[i], sm, 1.).probability(x[i]);
         sx -= x[i];  sm -= m[i];
      }
      return p;
   }


   if (j == 1) { 
      return binoexpand();
   }

   findpars();
   if (w < 0.04 && E < 10 && (!em || w > 0.004)) {
      return laplace();
   }

   return integrate();
}


/***********************************************************************
Methods for CMultiWalleniusNCHypergeometricMoments
***********************************************************************/

double CMultiWalleniusNCHypergeometricMoments::moments(double * mu, double * variance, int32_t * combinations) {
   // calculates mean and variance of multivariate Wallenius noncentral 
   // hypergeometric distribution by calculating all combinations of x-values.
   // Return value = sum of all probabilities. The deviation of this value 
   // from 1 is a measure of the accuracy.
   // Returns the mean to mean[0...colors-1]
   // Returns the variance to variance[0...colors-1]
   double sumf;                        // sum of all f(x) values
   int32_t msum;                       // temporary sum
   int i;                              // loop counter

   // get approximate mean
   mean(sx);
   // round mean to integers
   for (i=0; i < colors; i++) {
      xm[i] = (int32_t)(sx[i]+0.4999999);
   }

   // set up for recursive loops
   for (i=colors-1, msum=0; i >= 0; i--) {
      remaining[i] = msum;  msum += m[i];
   }
   for (i=0; i<colors; i++)  sx[i] = sxx[i] = 0.;
   sn = 0;

   // recursive loops to calculate sums  
   sumf = loop(n, 0);

   // calculate mean and variance
   for (i = 0; i < colors; i++) {
      mu[i] = sx[i]/sumf;
      variance[i] = sxx[i]/sumf - sx[i]*sx[i]/(sumf*sumf);
   }

   // return combinations and sum
   if (combinations) *combinations = sn;
   return sumf;
}


double CMultiWalleniusNCHypergeometricMoments::loop(int32_t n, int c) {
   // recursive function to loop through all combinations of x-values.
   // used by moments()
   int32_t x, x0;                      // x of color c
   int32_t xmin, xmax;                 // min and max of x[c]
   double s1, s2, sum = 0.;            // sum of f(x) values
   int i;                              // loop counter

   if (c < colors-1) {
      // not the last color
      // calculate min and max of x[c] for given x[0]..x[c-1]
      xmin = n - remaining[c];  if (xmin < 0) xmin = 0;
      xmax = m[c];  if (xmax > n) xmax = n;
      x0 = xm[c];  if (x0 < xmin) x0 = xmin;  if (x0 > xmax) x0 = xmax;
      // loop for all x[c] from mean and up
      for (x = x0, s2 = 0.; x <= xmax; x++) {
         xi[c] = x;
         sum += s1 = loop(n-x, c+1);             // recursive loop for remaining colors
         if (s1 < accuracy && s1 < s2) break;    // stop when values become negligible
         s2 = s1;
      }
      // loop for all x[c] from mean and down
      for (x = x0-1; x >= xmin; x--) {
         xi[c] = x;
         sum += s1 = loop(n-x, c+1);             // recursive loop for remaining colors
         if (s1 < accuracy && s1 < s2) break;    // stop when values become negligible
         s2 = s1;
      }
   }
   else {
      // last color
      xi[c] = n;
      s1 = probability(xi);
      for (i=0; i < colors; i++) {
         sx[i]  += s1 * xi[i];
         sxx[i] += s1 * xi[i] * xi[i];
      }
      sn++;
      sum = s1;
   }
   return sum;
}
