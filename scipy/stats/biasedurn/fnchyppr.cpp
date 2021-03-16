/*************************** fnchyppr.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2014-06-14
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Calculation of univariate and multivariate Fisher's noncentral hypergeometric
* probability distribution.
*
* This file contains source code for the class CFishersNCHypergeometric 
* and CMultiFishersNCHypergeometric defined in stocc.h.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2002-2014 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <string.h>                    // memcpy function
#include "stocc.h"                     // class definition


/***********************************************************************
Methods for class CFishersNCHypergeometric
***********************************************************************/

CFishersNCHypergeometric::CFishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy) {
   // constructor
   // set parameters
   this->n = n;  this->m = m;  this->N = N;
   this->odds = odds;  this->accuracy = accuracy;

   // check validity of parameters
   if (n < 0 || m < 0 || N < 0 || odds < 0. || n > N || m > N) {
      FatalError("Parameter out of range in class CFishersNCHypergeometric");
   }
   if (accuracy < 0) accuracy = 0;
   if (accuracy > 1) accuracy = 1;
   // initialize
   logodds = log(odds);  scale = rsum = 0.;
   ParametersChanged = 1;
   // calculate xmin and xmax
   xmin = m + n - N;  if (xmin < 0) xmin = 0;
   xmax = n;  if (xmax > m) xmax = m;
}


int32_t CFishersNCHypergeometric::mode(void) {
   // Find mode (exact)
   // Uses the method of Liao and Rosen, The American Statistician, vol 55,
   // no 4, 2001, p. 366-369.
   // Note that there is an error in Liao and Rosen's formula. 
   // Replace sgn(b) with -1 in Liao and Rosen's formula. 

   double A, B, C, D;                  // coefficients for quadratic equation
   double x;                           // mode
   int32_t L = m + n - N;
   int32_t m1 = m+1, n1 = n+1;

   if (odds == 1.) { 
      // simple hypergeometric
      x = (m + 1.) * (n + 1.) / (N + 2.);
   }
   else {
      // calculate analogously to Cornfield mean
      A = 1. - odds;
      B = (m1+n1)*odds - L; 
      C = -(double)m1*n1*odds;
      D = B*B -4*A*C;
      D = D > 0. ? sqrt(D) : 0.;
      x = (D - B)/(A+A);
   }
   return (int32_t)x;
}


double CFishersNCHypergeometric::mean(void) {
   // Find approximate mean
   // Calculation analogous with mode
   double a, b;                        // temporaries in calculation
   double mean;                        // mean

   if (odds == 1.) {                   // simple hypergeometric
      return double(m)*n/N;
   }
   // calculate Cornfield mean
   a = (m+n)*odds + (N-m-n); 
   b = a*a - 4.*odds*(odds-1.)*m*n;
   b = b > 0. ? sqrt(b) : 0.;
   mean = (a-b)/(2.*(odds-1.));
   return mean;
}


double CFishersNCHypergeometric::variance(void) {
   // find approximate variance (poor approximation)    
   double my = mean(); // approximate mean
   // find approximate variance from Fisher's noncentral hypergeometric approximation
   double r1 = my * (m-my); double r2 = (n-my)*(my+N-n-m);
   if (r1 <= 0. || r2 <= 0.) return 0.;
   double var = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));
   if (var < 0.) var = 0.;
   return var;
}

double CFishersNCHypergeometric::moments(double * mean_, double * var_) {
   // calculate exact mean and variance
   // return value = sum of f(x), expected = 1.
   double y, sy=0, sxy=0, sxxy=0, me1;
   int32_t x, xm, x1;
   const double accur = 0.1 * accuracy;     // accuracy of calculation
   xm = (int32_t)mean();                      // approximation to mean
   for (x=xm; x<=xmax; x++) {
      y = probability(x);
      x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
      sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
      if (y < accur && x != xm) break;
   }
   for (x=xm-1; x>=xmin; x--) {
      y = probability(x);
      x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
      sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
      if (y < accur) break;
   }
   me1 = sxy / sy;
   *mean_ = me1 + xm;
   y = sxxy / sy - me1 * me1;
   if (y < 0) y=0;
   *var_ = y;
   return sy;
}


double CFishersNCHypergeometric::probability(int32_t x) {
   // calculate probability function
   const double accur = accuracy * 0.1;// accuracy of calculation

   if (x < xmin || x > xmax) return 0;

   if (n == 0) return 1.;

   if (odds == 1.) {
      // central hypergeometric
      return exp(
         LnFac(m)   - LnFac(x)   - LnFac(m-x) +
         LnFac(N-m) - LnFac(n-x) - LnFac((N-m)-(n-x)) -
        (LnFac(N)   - LnFac(n)   - LnFac(N-n)));
   }

   if (odds == 0.) {
      if (n > N-m) FatalError("Not enough items with nonzero weight in CFishersNCHypergeometric::probability");
      return x == 0;
   }

   if (!rsum) {
      // first time. calculate rsum = reciprocal of sum of proportional 
      // function over all probable x values
      int32_t x1, x2;                    // x loop
      double y;                        // value of proportional function
      x1 = (int32_t)mean();              // start at mean
      if (x1 < xmin) x1 = xmin;
      x2 = x1 + 1;
      scale = 0.; scale = lng(x1);     // calculate scale to avoid overflow
      rsum = 1.;                       // = exp(lng(x1)) with this scale
      for (x1--; x1 >= xmin; x1--) {
         rsum += y = exp(lng(x1));     // sum from x1 and down 
         if (y < accur) break;         // until value becomes negligible
      }
      for (; x2 <= xmax; x2++) {       // sum from x2 and up
         rsum += y = exp(lng(x2));
         if (y < accur) break;         // until value becomes negligible
      }
      rsum = 1. / rsum;                // save reciprocal sum
   }
   return exp(lng(x)) * rsum;          // function value
}

double CFishersNCHypergeometric::probabilityRatio(int32_t x, int32_t x0) {
   // Calculate probability ratio f(x)/f(x0)
   // This is much faster than calculating a single probability because
   // rsum is not needed
   double a1, a2, a3, a4, f1, f2, f3, f4;
   int32_t y, dx = x - x0;
   int invert = 0;

   if (x < xmin || x > xmax) return 0.;
   if (x0 < xmin || x0 > xmax) {
      FatalError("Infinity in CFishersNCHypergeometric::probabilityRatio");
   }
   if (dx == 0.) return 1.;
   if (dx < 0.) {
      invert = 1;
      dx = -dx;
      y = x;  x = x0;  x0 = y;
   }
   a1 = m - x0;  a2 = n - x0;  a3 = x;  a4 = N - m - n + x;
   if (dx <= 28 && x <= 100000) {      // avoid overflow
      // direct calculation
      f1 = f2 = 1.;
      // compute ratio of binomials
      for (y = 0; y < dx; y++) {
         f1 *= a1-- * a2--;
         f2 *= a3-- * a4--;
      }
      // compute odds^dx
      f3 = 1.;  f4 = odds;  y = dx;
      while (y) {
         if (f4 < 1.E-100) {
            f3 = 0.;  break;           // avoid underflow
         }
         if (y & 1) f3 *= f4;
         f4 *= f4;
         y = (unsigned long)(y) >> 1;
      }
      f1 = f3 * f1 / f2;
      if (invert) f1 = 1. / f1;
   }
   else {
      // use logarithms
      f1 = FallingFactorial(a1,dx) + FallingFactorial(a2,dx) -
           FallingFactorial(a3,dx) - FallingFactorial(a4,dx) +
           dx * log(odds);
      if (invert) f1 = -f1;
      f1 = exp(f1);
   }
   return f1;
}

double CFishersNCHypergeometric::MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff) {
   // Makes a table of Fisher's noncentral hypergeometric probabilities.
   // Results are returned in the array table of size MaxLength.
   // The values are scaled so that the highest value is 1. The return value
   // is the sum, s, of all the values in the table. The normalized
   // probabilities are obtained by multiplying all values in the table by
   // 1/s.
   // The tails are cut off where the values are < cutoff, so that 
   // *xfirst may be > xmin and *xlast may be < xmax.
   // The value of cutoff will be 0.01 * accuracy if not specified.
   // The first and last x value represented in the table are returned in 
   // *xfirst and *xlast. The resulting probability values are returned in the 
   // first (*xlast - *xfirst + 1) positions of table. If this would require
   // more than MaxLength values then the table is filled with as many 
   // correct values as possible.
   //
   // The function will return the desired length of table when MaxLength = 0.

   double f;                           // probability function value
   double sum;                         // sum of table values
   double a1, a2, b1, b2;              // factors in recursive calculation of f(x)
   int32_t x;                          // x value
   int32_t x1, x2;                     // lowest and highest x
   int32_t i, i0, i1, i2;              // table index
   int32_t mode = this->mode();        // mode
   int32_t L = n + m - N;              // parameter
   int32_t DesiredLength;              // desired length of table

   // limits for x
   x1 = (L > 0) ? L : 0;               // xmin
   x2 = (n < m) ? n : m;               // xmax

   // special cases
   if (x1 == x2) goto DETERMINISTIC;
   if (odds <= 0.) {
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

   if (MaxLength <= 0) {
      // Return UseTable and LengthNeeded
      DesiredLength = x2 - x1 + 1;     // max length of table
      if (DesiredLength > 200) {
         double sd = sqrt(variance()); // calculate approximate standard deviation
         // estimate number of standard deviations to include from normal distribution
         i = (int32_t)(NumSD(accuracy) * sd + 0.5);
         if (DesiredLength > i) DesiredLength = i;
      }
      if (xfirst) *xfirst = 1;         // for analogy with CWalleniusNCHypergeometric::MakeTable
      return DesiredLength;
   }

   // place mode in the table
   if (mode - x1 <= MaxLength/2) {
      // There is enough space for left tail
      i0 = mode - x1;
   }
   else if (x2 - mode <= MaxLength/2) {
      // There is enough space for right tail
      i0 = MaxLength - x2 + mode - 1;
      if (i0 < 0) i0 = 0;
   }
   else {
      // There is not enough space for any of the tails. Place mode in middle of table
      i0 = MaxLength/2;
   }
   // Table start index
   i1 = i0 - mode + x1;  if (i1 < 0) i1 = 0;

   // Table end index
   i2 = i0 + x2 - mode;  if (i2 > MaxLength-1) i2 = MaxLength-1;

   // make center
   table[i0] = sum = f = 1.;

   // make left tail
   x = mode;
   a1 = m + 1 - x;  a2 = n + 1 - x;
   b1 = x;  b2 = x - L;
   for (i = i0 - 1; i >= i1; i--) {
      f *= b1 * b2 / (a1 * a2 * odds); // recursive formula
      a1++;  a2++;  b1--;  b2--;
      sum += table[i] = f;
      if (f < cutoff) {
         i1 = i;  break;               // cut off tail if < accuracy
      }
   }
   if (i1 > 0) {
      // move table down for cut-off left tail
      memcpy(table, table+i1, (i0-i1+1)*sizeof(*table));
      // adjust indices
      i0 -= i1;  i2 -= i1;  i1 = 0;
   }
   // make right tail
   x = mode + 1;
   a1 = m + 1 - x;  a2 = n + 1 - x;
   b1 = x;  b2 = x - L;
   f = 1.;
   for (i = i0 + 1; i <= i2; i++) {
      f *= a1 * a2 * odds / (b1 * b2); // recursive formula
      a1--;  a2--;  b1++;  b2++;
      sum += table[i] = f;
      if (f < cutoff) {
         i2 = i;  break;               // cut off tail if < accuracy
      }
   }
   // x limits
   *xfirst = mode - (i0 - i1);
   *xlast  = mode + (i2 - i0);

   return sum;
}


double CFishersNCHypergeometric::lng(int32_t x) {
   // natural log of proportional function
   // returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!*odds^x)
   int32_t x2 = n - x,  m2 = N - m;
   if (ParametersChanged) {
      mFac = LnFac(m) + LnFac(m2);
      xLast = -99; ParametersChanged = 0;
   }
   if (m < FAK_LEN && m2 < FAK_LEN)  goto DEFLT;
   switch (x - xLast) {
  case 0:   // x unchanged
     break;
  case 1:   // x incremented. calculate from previous value
     xFac += log (double(x) * (m2-x2) / (double(x2+1)*(m-x+1)));
     break;
  case -1:  // x decremented. calculate from previous value
     xFac += log (double(x2) * (m-x) / (double(x+1)*(m2-x2+1)));
     break;
  default: DEFLT: // calculate all
     xFac = LnFac(x) + LnFac(x2) + LnFac(m-x) + LnFac(m2-x2);
   }
   xLast = x;
   return mFac - xFac + x * logodds - scale;
}


/***********************************************************************
calculation methods in class CMultiFishersNCHypergeometric
***********************************************************************/

CMultiFishersNCHypergeometric::CMultiFishersNCHypergeometric(int32_t n_, int32_t * m_, double * odds_, int colors_, double accuracy_) {
   // constructor
   int32_t N1;
   int i;
   // copy parameters
   n = n_;  m = m_;  odds = odds_;  colors = colors_;  accuracy = accuracy_;
   // check if parameters are valid 
   // (Note: there is a more thorough test for validity in the BiasedUrn package)
   for (N = N1 = 0, i = 0; i < colors; i++) {
      if (m[i] < 0 || odds[i] < 0) FatalError("Parameter negative in constructor for CMultiFishersNCHypergeometric");
      N += m[i];
      if (odds[i]) N1 += m[i];
   }
   if (N < n) FatalError("Not enough items in constructor for CMultiFishersNCHypergeometric");
   if (N1< n) FatalError("Not enough items with nonzero weight in constructor for CMultiFishersNCHypergeometric");

   // calculate mFac and logodds
   for (i=0, mFac=0.; i < colors; i++) {
      mFac += LnFac(m[i]);
      logodds[i] = log(odds[i]);
   }
   // initialize
   sn = 0;
}


void CMultiFishersNCHypergeometric::mean(double * mu) {
   // calculates approximate mean of multivariate Fisher's noncentral
   // hypergeometric distribution. Result is returned in mu[0..colors-1].
   // The calculation is reasonably fast.
   // Note: The version in BiasedUrn package deals with unused colors
   double r, r1;                       // iteration variable
   double q;                           // mean of color i
   double W;                           // total weight
   int i;                              // color index
   int iter = 0;                       // iteration counter

   if (colors < 3) {
      // simple cases
      if (colors == 1) mu[0] = n;
      if (colors == 2) {
         mu[0] = CFishersNCHypergeometric(n,m[0],m[0]+m[1],odds[0]/odds[1]).mean();
         mu[1] = n - mu[0];
      }
      return;
   }
   if (n == N) {
      // Taking all balls
      for (i = 0; i < colors; i++) mu[i] = m[i];
      return;
   }

   // initial guess for r
   for (i=0, W=0.; i < colors; i++) W += m[i] * odds[i];
   r = (double)n * N / ((N-n)*W);

   // iteration loop to find r
   do {
      r1 = r;
      for (i=0, q=0.; i < colors; i++) {
         q += m[i] * r * odds[i] / (r * odds[i] + 1.);
      }
      r *= n * (N-q) / (q * (N-n));
      if (++iter > 100) FatalError("convergence problem in function CMultiFishersNCHypergeometric::mean");
   }
   while (fabs(r-r1) > 1E-5);

   // store result
   for (i=0; i < colors; i++) {
      mu[i] = m[i] * r * odds[i] / (r * odds[i] + 1.);
   }
}


void CMultiFishersNCHypergeometric::variance(double * var) {
   // calculates approximate variance of multivariate Fisher's noncentral
   // hypergeometric distribution (accuracy is not too good).
   // Result is returned in variance[0..colors-1].
   // The calculation is reasonably fast.
   // Note: The version in BiasedUrn package deals with unused colors
   double r1, r2;
   double mu[MAXCOLORS];
   int i;
   mean(mu);
   for (i=0; i<colors; i++) {
      r1 = mu[i] * (m[i]-mu[i]);
      r2 = (n-mu[i])*(mu[i]+N-n-m[i]);
      if (r1 <= 0. || r2 <= 0.) {
         var[i] = 0.;
      }
      else {
         var[i] = N*r1*r2/((N-1)*(m[i]*r2+(N-m[i])*r1));
      }
   }
}


double CMultiFishersNCHypergeometric::probability(int32_t * x) {
   // Calculate probability function.
   // Note: The first-time call takes very long time because it requires
   // a calculation of all possible x combinations with probability >
   // accuracy, which may be extreme.
   // The calculation uses logarithms to avoid overflow. 
   // (Recursive calculation may be faster, but this has not been implemented)
   // Note: The version in BiasedUrn package deals with unused colors
   int32_t xsum;  int i, em;
   for (xsum = i = 0; i < colors; i++)  xsum += x[i];
   if (xsum != n) {
      FatalError("sum of x values not equal to n in function CMultiFishersNCHypergeometric::probability");
   }

   for (i = em = 0; i < colors; i++) {
      if (x[i] > m[i] || x[i] < 0 || x[i] < n - N + m[i]) return 0.;
      if (odds[i] == 0. && x[i]) return 0.;
      if (x[i] == m[i] || odds[i] == 0.) em++;
   }

   if (n == 0 || em == colors) return 1.;

   if (sn == 0) SumOfAll();            // first time initialize
   return exp(lng(x)) * rsum;          // function value
}


double CMultiFishersNCHypergeometric::moments(double * mean, double * variance, int32_t * combinations) {
   // calculates mean and variance of the Fisher's noncentral hypergeometric 
   // distribution by calculating all combinations of x-values with
   // probability > accuracy.
   // Return value = 1.
   // Returns the mean in mean[0...colors-1]
   // Returns the variance in variance[0...colors-1]
   // Note: The version in BiasedUrn package deals with unused colors

   int i;                              // color index
   if (sn == 0) {
      // first time initialization includes calculation of mean and variance
      SumOfAll();
   }
   // just copy results
   for (i=0; i < colors; i++) {
      mean[i] = sx[i];
      variance[i] = sxx[i];
   }
   if (combinations) *combinations = sn;
   return 1.;
}


void CMultiFishersNCHypergeometric::SumOfAll() {
   // this function does the very time consuming job of calculating the sum
   // of the proportional function g(x) over all possible combinations of
   // the x[i] values with probability > accuracy. These combinations are 
   // generated by the recursive function loop().
   // The mean and variance are generated as by-products.

   int i;                              // color index
   int32_t msum;                         // sum of m[i]

   // get approximate mean
   mean(sx);
   // round mean to integers
   for (i=0, msum=0; i < colors; i++) {
      msum += xm[i] = (int32_t)(sx[i]+0.4999999);}
   // adjust truncated x values to make the sum = n
   msum -= n;
   for (i = 0; msum < 0; i++) {
      if (xm[i] < m[i]) {
         xm[i]++; msum++;
      }
   }
   for (i = 0; msum > 0; i++) {
      if (xm[i] > 0) {
         xm[i]--; msum--;
      }
   }

   // adjust scale factor to g(mean) to avoid overflow
   scale = 0.; scale = lng(xm);

   // initialize for recursive loops
   sn = 0;
   for (i = colors-1, msum = 0; i >= 0; i--) {
      remaining[i] = msum;  msum += m[i];
   }
   for (i = 0; i < colors; i++) {
      sx[i] = 0;  sxx[i] = 0;
   }

   // recursive loops to calculate sums of g(x) over all x combinations
   rsum = 1. / loop(n, 0);

   // calculate mean and variance
   for (i = 0; i < colors; i++) {
      sxx[i] = sxx[i]*rsum - sx[i]*sx[i]*rsum*rsum;
      sx[i] = sx[i]*rsum;
   }
}


double CMultiFishersNCHypergeometric::loop(int32_t n, int c) {
   // recursive function to loop through all combinations of x-values.
   // used by SumOfAll
   int32_t x, x0;                        // x of color c
   int32_t xmin, xmax;                   // min and max of x[c]
   double s1, s2, sum = 0.;            // sum of g(x) values
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
         sum += s1 = loop(n-x, c+1); // recursive loop for remaining colors
         if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
         s2 = s1;
      }
      // loop for all x[c] from mean and down
      for (x = x0-1; x >= xmin; x--) {
         xi[c] = x;
         sum += s1 = loop(n-x, c+1);   // recursive loop for remaining colors
         if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
         s2 = s1;
      }
   }
   else {
      // last color
      xi[c] = n;
      // sums and squaresums    
      s1 = exp(lng(xi));               // proportional function g(x)
      for (i = 0; i < colors; i++) {   // update sums
         sx[i]  += s1 * xi[i];
         sxx[i] += s1 * xi[i] * xi[i];
      }
      sn++;
      sum += s1;
   }
   return sum;
}


double CMultiFishersNCHypergeometric::lng(int32_t * x) {
   // natural log of proportional function g(x)
   double y = 0.;
   int i;
   for (i = 0; i < colors; i++) {
      y += x[i]*logodds[i] - LnFac(x[i]) - LnFac(m[i]-x[i]);
   }
   return mFac + y - scale;
}
