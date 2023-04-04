/*****************************   stocc.h   **********************************
* Author:        Agner Fog
* Date created:  2004-01-08
* Last modified: 2013-09-20
* Project:       randomc.h
* Source URL:    www.agner.org/random
*
* Description:
* This file contains function prototypes and class declarations for the C++ 
* library of non-uniform random number generators. Most functions are fast and 
* accurate, even for extreme values of the parameters.
*
*
* functions without classes:
* ==========================
*
* void EndOfProgram(void);
* System-specific exit code. You may modify this to make it fit your
* user interface.
*
* void FatalError(const char * ErrorText);
* Used for outputting error messages from the other functions and classes.
* You may have to modify this function to make it fit your user interface.
*
* double Erf (double x);
* Calculates the error function, which is the integral of the normal distribution.
*
* double LnFac(int32_t n);
* Calculates the natural logarithm of the factorial of n.
*
*
* class StochasticLib1:
* ====================
* This class can be derived from any of the uniform random number generators
* defined in randomc.h. StochasticLib1 provides the following non-uniform random 
* variate generators:
*
* int Bernoulli(double p);
* Bernoulli distribution. Gives 0 or 1 with probability 1-p and p.
*
* double Normal(double m, double s);
* Normal distribution with mean m and standard deviation s.
*
* double NormalTrunc(double m, double s, double limit);
* Truncated normal distribution with tails cut off at m +/- limit
*
* int32_t Poisson (double L);
* Poisson distribution with mean L.
*
* int32_t Binomial (int32_t n, double p);
* Binomial distribution. n trials with probability p.
*
* int32_t Hypergeometric (int32_t n, int32_t m, int32_t N);
* Hypergeometric distribution. Taking n items out N, m of which are colored.
*
* void Multinomial (int32_t * destination, double * source, int32_t n, int colors);
* void Multinomial (int32_t * destination, int32_t * source, int32_t n, int colors);
* Multivariate binomial distribution.
*
* void MultiHypergeometric (int32_t * destination, int32_t * source, int32_t n, int colors);
* Multivariate hypergeometric distribution.
*
* void Shuffle(int * list, int min, int n);
* Shuffle a list of integers.
*
*
* class StochasticLib2:
* =====================
* This class is derived from class StochasticLib1. It redefines the functions
* Poisson, Binomial and HyperGeometric.
* In StochasticLib1, these functions are optimized for being called with 
* parameters that vary. In StochasticLib2, the same functions are optimized
* for being called repeatedly with the same parameters. If your parameters
* seldom vary, then StochasticLib2 is faster. The two classes use different
* calculation methods, both of which are accurate.
*
*
* class StochasticLib3:
* =====================
* This class can be derived from either StochasticLib1 or StochasticLib2, 
* whichever is preferred. It contains functions for generating variates with
* the univariate and multivariate Wallenius' and Fisher's noncentral
* hypergeometric distributions.
*
* int32_t WalleniusNCHyp (int32_t n, int32_t m, int32_t N, double odds);
* Sampling from Wallenius' noncentral hypergeometric distribution, which is 
* what you get when taking n items out N, m of which are colored, without 
* replacement, with bias.
*
* int32_t FishersNCHyp (int32_t n, int32_t m, int32_t N, double odds);
* Sampling from Fisher's noncentral hypergeometric distribution which is the
* conditional distribution of independent binomial variates given their sum n.
*
* void MultiWalleniusNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors);
* Sampling from multivariate Wallenius' noncentral hypergeometric distribution.
*
* void MultiFishersNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors);
* Sampling from multivariate Fisher's noncentral hypergeometric distribution.
*
*
* Uniform random number generators (integer and float) are also available, as
* these are inherited from the random number generator class that is the base
* class of StochasticLib1.
*
*
* class CWalleniusNCHypergeometric
* ================================
* This class implements various methods for calculating the probability 
* function and the mean and variance of the univariate Wallenius' noncentral 
* hypergeometric distribution. It is used by StochasticLib3 and can also be 
* used independently.
*
*
* class CMultiWalleniusNCHypergeometric
* =====================================
* This class implements various methods for calculating the probability func-
* tion and the mean of the multivariate Wallenius' noncentral hypergeometric
* distribution. It is used by StochasticLib3 and can also be used independently.
*
*
* class CMultiWalleniusNCHypergeometricMoments
* ============================================
* This class calculates the exact mean and variance of the multivariate
* Wallenius' noncentral hypergeometric probability distribution.
*
*
* class CFishersNCHypergeometric
* ==============================
* This class calculates the probability function and the mean and variance 
* of Fisher's noncentral hypergeometric distribution.
*
*
* class CMultiFishersNCHypergeometric
* ===================================
* This class calculates the probability function and the mean and variance 
* of the multivariate Fisher's noncentral hypergeometric distribution.
*
*
* source code:
* ============
* The code for EndOfProgram and FatalError is found in the file userintf.cpp.
* The code for the functions in StochasticLib1 is found in the file stoc1.cpp.
* The code for the functions in StochasticLib2 is found in the file stoc2.cpp.
* The code for the functions in StochasticLib3 is found in the file stoc3.cpp.
* The code for the functions in CWalleniusNCHypergeometric, 
* CMultiWalleniusNCHypergeometric and CMultiWalleniusNCHypergeometricMoments
* is found in the file wnchyppr.cpp.
* The code for the functions in CFishersNCHypergeometric and 
* CMultiFishersNCHypergeometric is found in the file fnchyppr.cpp
* LnFac is found in stoc1.cpp.
* Erf is found in wnchyppr.cpp.
*
*
* Examples:
* =========
* The file ex-stoc.cpp contains an example of how to use this class library.
*
* The file ex-cards.cpp contains an example of how to shuffle a list of items.
*
* The file ex-lotto.cpp contains an example of how to generate a sequence of
* random integers where no number can occur more than once.
*
* The file testbino.cpp contains an example of sampling from the binomial distribution.
*
* The file testhype.cpp contains an example of sampling from the hypergeometric distribution.
*
* The file testpois.cpp contains an example of sampling from the poisson distribution.
*
* The file testwnch.cpp contains an example of sampling from Wallenius noncentral hypergeometric distribution.
*
* The file testfnch.cpp contains an example of sampling from Fisher's noncentral hypergeometric distribution.
*
* The file testmwnc.cpp contains an example of sampling from the multivariate Wallenius noncentral hypergeometric distribution.
*
* The file testmfnc.cpp contains an example of sampling from the multivariate Fisher's noncentral hypergeometric distribution.
*
* The file evolc.zip contains examples of how to simulate biological evolution using this class library.
*
*
* Documentation:
* ==============
* The file ran-instructions.pdf contains further documentation and 
* instructions for these random number generators.
*
* The file distrib.pdf contains definitions of the standard statistic distributions:
* Bernoulli, Normal, Poisson, Binomial, Hypergeometric, Multinomial, MultiHypergeometric.
*
* The file sampmet.pdf contains theoretical descriptions of the methods used
* for sampling from these distributions.
*
* The file nchyp.pdf, available from www.agner.org/random/, contains
* definitions of the univariate and multivariate Wallenius and Fisher's 
* noncentral hypergeometric distributions and theoretical explanations of 
* the methods for calculating and sampling from these.
*
* Copyright 2004-2013 by Agner Fog. 
* Released under SciPy's license with permission of Agner Fog; see license.txt
*******************************************************************************/

#ifndef STOCC_H
#define STOCC_H

#include <math.h>
#include "randomc.h"

#ifdef R_BUILD
   #include "stocR.h"           // Include this when building R-language interface
#endif


/***********************************************************************
 Choose which uniform random number generator to base these classes on
***********************************************************************/

// STOC_BASE defines which base class to use for the non-uniform
// random number generator classes StochasticLib1, 2, and 3.
#ifndef STOC_BASE
   #ifdef R_BUILD
      // Inherit from StocRBase when building for R-language interface
      #define STOC_BASE StocRBase
   #else
      #define STOC_BASE CRandomMersenne     // C++ Mersenne Twister
      // Or choose any other random number generator base class, for example:
      //#include "randoma.h"
      //#define STOC_BASE CRandomSFMTA      // Binary library SFMT generator
   #endif
#endif

/***********************************************************************
         Other simple functions
***********************************************************************/

double LnFac(int32_t n);               // log factorial (stoc1.cpp)
double LnFacr(double x);               // log factorial of non-integer (wnchyppr.cpp)
double FallingFactorial(double a, double b); // Falling factorial (wnchyppr.cpp)
double Erf (double x);                 // error function (wnchyppr.cpp)
int32_t FloorLog2(float x);            // floor(log2(x)) for x > 0 (wnchyppr.cpp)
int NumSD (double accuracy);           // used internally for determining summation interval


/***********************************************************************
         Constants and tables
***********************************************************************/

// Maximum number of colors in the multivariate distributions
#ifndef MAXCOLORS
   #define MAXCOLORS 32                // You may change this value
#endif

// constant for LnFac function:
static const int FAK_LEN = 1024;       // length of factorial table

// The following tables are tables of residues of a certain expansion
// of the error function. These tables are used in the Laplace method
// for calculating Wallenius' noncentral hypergeometric distribution.
// There are ERFRES_N tables covering desired precisions from
// 2^(-ERFRES_B) to 2^(-ERFRES_E). Only the table that matches the
// desired precision is used. The tables are defined in erfres.h which
// is included in wnchyppr.cpp.

// constants for ErfRes tables:
static const int ERFRES_B = 16;        // begin: -log2 of lowest precision
static const int ERFRES_E = 40;        // end:   -log2 of highest precision
static const int ERFRES_S =  2;        // step size from begin to end
static const int ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1; // number of tables
static const int ERFRES_L = 48;        // length of each table

// tables of error function residues:
extern "C" double ErfRes [ERFRES_N][ERFRES_L];

// number of std. deviations to include in integral to obtain desired precision:
extern "C" double NumSDev[ERFRES_N];


/***********************************************************************
         Class StochasticLib1
***********************************************************************/

class StochasticLib1 : public STOC_BASE {
   // This class encapsulates the random variate generating functions.
   // May be derived from any of the random number generators.
public:
   StochasticLib1 (int seed);          // Constructor
   int Bernoulli(double p);            // Bernoulli distribution
   double Normal(double m, double s);  // Normal distribution
   double NormalTrunc(double m, double s, double limit); // Truncated normal distribution
   int32_t Poisson (double L);         // Poisson distribution
   int32_t Binomial (int32_t n, double p); // Binomial distribution
   int32_t Hypergeometric (int32_t n, int32_t m, int32_t N); // Hypergeometric distribution
   void Multinomial (int32_t * destination, double * source, int32_t n, int colors); // Multinomial distribution
   void Multinomial (int32_t * destination, int32_t * source, int32_t n, int colors);// Multinomial distribution
   void MultiHypergeometric (int32_t * destination, int32_t * source, int32_t n, int colors); // Multivariate hypergeometric distribution
   void Shuffle(int * list, int min, int n); // Shuffle integers

   // functions used internally
protected:
   static double fc_lnpk(int32_t k, int32_t N_Mn, int32_t M, int32_t n); // used by Hypergeometric

   // subfunctions for each approximation method
   int32_t PoissonInver(double L);                         // poisson by inversion
   int32_t PoissonRatioUniforms(double L);                 // poisson by ratio of uniforms
   int32_t PoissonLow(double L);                           // poisson for extremely low L
   int32_t BinomialInver (int32_t n, double p);            // binomial by inversion
   int32_t BinomialRatioOfUniforms (int32_t n, double p);  // binomial by ratio of uniforms
   int32_t HypInversionMod (int32_t n, int32_t M, int32_t N);  // hypergeometric by inversion searching from mode
   int32_t HypRatioOfUnifoms (int32_t n, int32_t M, int32_t N);// hypergeometric by ratio of uniforms method

   // Variables specific to each distribution:
   // Variables used by Normal distribution
   double normal_x2;  int normal_x2_valid;

   // Variables used by Hypergeometric distribution
   int32_t  hyp_n_last, hyp_m_last, hyp_N_last;            // Last values of parameters
   int32_t  hyp_mode, hyp_mp;                              // Mode, mode+1
   int32_t  hyp_bound;                                     // Safety upper bound
   double hyp_a;                                           // hat center
   double hyp_h;                                           // hat width
   double hyp_fm;                                          // Value at mode

   // Variables used by Poisson distribution
   double pois_L_last;                                     // previous value of L
   double pois_f0;                                         // value at x=0 or at mode
   double pois_a;                                          // hat center
   double pois_h;                                          // hat width
   double pois_g;                                          // ln(L)
   int32_t  pois_bound;                                    // upper bound

   // Variables used by Binomial distribution
   int32_t bino_n_last;                                    // last n
   double bino_p_last;                                     // last p
   int32_t bino_mode;                                      // mode
   int32_t bino_bound;                                     // upper bound
   double bino_a;                                          // hat center
   double bino_h;                                          // hat width
   double bino_g;                                          // value at mode
   double bino_r1;                                         // p/(1-p) or ln(p/(1-p))
};


/***********************************************************************
Class StochasticLib2
***********************************************************************/

class StochasticLib2 : public StochasticLib1 {
   // derived class, redefining some functions
public:
   int32_t Poisson (double L);                             // Poisson distribution
   int32_t Binomial (int32_t n, double p);                 // Binomial distribution
   int32_t Hypergeometric(int32_t n, int32_t M, int32_t N);// Hypergeometric distribution
   StochasticLib2(int seed):StochasticLib1(seed){};        // Constructor  

   // subfunctions for each approximation method:
protected:
   int32_t PoissonModeSearch(double L);                    // poisson by search from mode
   int32_t PoissonPatchwork(double L);                     // poisson by patchwork rejection
   static double PoissonF(int32_t k, double l_nu, double c_pm); // used by PoissonPatchwork
   int32_t BinomialModeSearch(int32_t n, double p);        // binomial by search from mode
   int32_t BinomialPatchwork(int32_t n, double p);         // binomial by patchwork rejection
   double BinomialF(int32_t k, int32_t n, double l_pq, double c_pm); // used by BinomialPatchwork
   int32_t HypPatchwork (int32_t n, int32_t M, int32_t N); // hypergeometric by patchwork rejection

   // Variables used by Binomial distribution
   int32_t  Bino_k1, Bino_k2, Bino_k4, Bino_k5;
   double Bino_dl, Bino_dr, Bino_r1, Bino_r2, Bino_r4, Bino_r5, 
      Bino_ll, Bino_lr, Bino_l_pq, Bino_c_pm,
      Bino_f1, Bino_f2, Bino_f4, Bino_f5, 
      Bino_p1, Bino_p2, Bino_p3, Bino_p4, Bino_p5, Bino_p6;

   // Variables used by Poisson distribution
   int32_t  Pois_k1, Pois_k2, Pois_k4, Pois_k5;
   double Pois_dl, Pois_dr, Pois_r1, Pois_r2, Pois_r4, Pois_r5, 
      Pois_ll, Pois_lr, Pois_l_my, Pois_c_pm,
      Pois_f1, Pois_f2, Pois_f4, Pois_f5, 
      Pois_p1, Pois_p2, Pois_p3, Pois_p4, Pois_p5, Pois_p6;

   // Variables used by Hypergeometric distribution
   int32_t  Hyp_L, Hyp_k1, Hyp_k2, Hyp_k4, Hyp_k5;
   double Hyp_dl, Hyp_dr, 
      Hyp_r1, Hyp_r2, Hyp_r4, Hyp_r5, 
      Hyp_ll, Hyp_lr, Hyp_c_pm, 
      Hyp_f1, Hyp_f2, Hyp_f4, Hyp_f5, 
      Hyp_p1, Hyp_p2, Hyp_p3, Hyp_p4, Hyp_p5, Hyp_p6;
};


/***********************************************************************
Class StochasticLib3
***********************************************************************/

class StochasticLib3 : public StochasticLib1 {
   // This class can be derived from either StochasticLib1 or StochasticLib2.
   // Adds more probability distributions
public:
   StochasticLib3(int seed);           // Constructor
   void SetAccuracy(double accur);     // Define accuracy of calculations
   int32_t WalleniusNCHyp (int32_t n, int32_t m, int32_t N, double odds); // Wallenius noncentral hypergeometric distribution
   int32_t FishersNCHyp (int32_t n, int32_t m, int32_t N, double odds); // Fisher's noncentral hypergeometric distribution
   void MultiWalleniusNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Multivariate Wallenius noncentral hypergeometric distribution
   void MultiComplWalleniusNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Multivariate complementary Wallenius noncentral hypergeometric distribution
   void MultiFishersNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Multivariate Fisher's noncentral hypergeometric distribution
   // subfunctions for each approximation method
protected:
   int32_t WalleniusNCHypUrn (int32_t n, int32_t m, int32_t N, double odds); // WalleniusNCHyp by urn model
   int32_t WalleniusNCHypInversion (int32_t n, int32_t m, int32_t N, double odds); // WalleniusNCHyp by inversion method
   int32_t WalleniusNCHypTable (int32_t n, int32_t m, int32_t N, double odds); // WalleniusNCHyp by table method
   int32_t WalleniusNCHypRatioOfUnifoms (int32_t n, int32_t m, int32_t N, double odds); // WalleniusNCHyp by ratio-of-uniforms
   int32_t FishersNCHypInversion (int32_t n, int32_t m, int32_t N, double odds); // FishersNCHyp by inversion
   int32_t FishersNCHypRatioOfUnifoms (int32_t n, int32_t m, int32_t N, double odds); // FishersNCHyp by ratio-of-uniforms

   // variables
   double accuracy;                                        // desired accuracy of calculations

   // Variables for Fisher
   int32_t fnc_n_last, fnc_m_last, fnc_N_last;             // last values of parameters
   int32_t fnc_bound;                                      // upper bound
   double fnc_o_last;
   double fnc_f0, fnc_scale;
   double fnc_a;                                           // hat center
   double fnc_h;                                           // hat width
   double fnc_lfm;                                         // ln(f(mode))
   double fnc_logb;                                        // ln(odds)

   // variables for Wallenius
   int32_t wnc_n_last, wnc_m_last, wnc_N_last;             // previous parameters
   double wnc_o_last;
   int32_t wnc_bound1, wnc_bound2;                         // lower and upper bound
   int32_t wnc_mode;                                       // mode
   double wnc_a;                                           // hat center
   double wnc_h;                                           // hat width
   double wnc_k;                                           // probability value at mode
   int UseChopDown;                                        // use chop down inversion instead
   #define WALL_TABLELENGTH  512                           // max length of table
   double wall_ytable[WALL_TABLELENGTH];                   // table of probability values
   int32_t wall_tablen;                                    // length of table
   int32_t wall_x1;                                        // lower x limit for table
};


/***********************************************************************
Class CWalleniusNCHypergeometric
***********************************************************************/

class CWalleniusNCHypergeometric {
   // This class contains methods for calculating the univariate
   // Wallenius' noncentral hypergeometric probability function
public:
   CWalleniusNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy=1.E-8); // constructor
   void SetParameters(int32_t n, int32_t m, int32_t N, double odds); // change parameters
   double probability(int32_t x);                          // calculate probability function
   int32_t MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                                      // approximate mean
   double variance(void);                                  // approximate variance (poor approximation)
   int32_t mode(void);                                     // calculate mode
   double moments(double * mean, double * var);            // calculate exact mean and variance
   int BernouilliH(int32_t x, double h, double rh, StochasticLib1 *sto); // used by rejection method

   // implementations of different calculation methods
protected:
   double recursive(void);                                 // recursive calculation
   double binoexpand(void);                                // binomial expansion of integrand
   double laplace(void);                                   // Laplace's method with narrow integration interval
   double integrate(void);                                 // numerical integration

   // other subfunctions
   double lnbico(void);                                    // natural log of binomial coefficients
   void findpars(void);                                    // calculate r, w, E
   double integrate_step(double a, double b);              // used by integrate()
   double search_inflect(double t_from, double t_to);      // used by integrate()

   // parameters
   double omega;                                           // Odds
   int32_t n, m, N, x;                                     // Parameters
   int32_t xmin, xmax;                                     // Minimum and maximum x
   double accuracy;                                        // Desired precision
   // parameters used by lnbico
   int32_t xLastBico;
   double bico, mFac, xFac;
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   int32_t xLastFindpars;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometric
***********************************************************************/

class CMultiWalleniusNCHypergeometric {
   // This class encapsulates the different methods for calculating the
   // multivariate Wallenius noncentral hypergeometric probability function
public:
   CMultiWalleniusNCHypergeometric(int32_t n, int32_t * m, double * odds, int colors, double accuracy=1.E-8); // constructor
   void SetParameters(int32_t n, int32_t * m, double * odds, int colors); // change parameters
   double probability(int32_t * x);                        // calculate probability function
   void mean(double * mu);                                 // calculate approximate mean

      // implementations of different calculation methods
protected:
   double binoexpand(void);                                // binomial expansion of integrand
   double laplace(void);                                   // Laplace's method with narrow integration interval
   double integrate(void);                                 // numerical integration

   // other subfunctions
   double lnbico(void);                                    // natural log of binomial coefficients
   void findpars(void);                                    // calculate r, w, E
   double integrate_step(double a, double b);              // used by integrate()
   double search_inflect(double t_from, double t_to);      // used by integrate()

   // parameters
   double * omega;                                         // odds
   double accuracy;                                        // desired accuracy
   int32_t n;                                              // sample size
   int32_t N;                                              // total items in urn
   int32_t * m;                                            // items of each color in urn
   int32_t * x;                                            // items of each color sampled
   int colors;                                             // number of different colors
   int Dummy_align;                                        // filler
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   // generated by lnbico
   double bico;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometricMoments
***********************************************************************/

class CMultiWalleniusNCHypergeometricMoments: public CMultiWalleniusNCHypergeometric {
   // This class calculates the exact mean and variance of the multivariate
   // Wallenius noncentral hypergeometric distribution by calculating all the 
   // possible x-combinations with probability < accuracy
public:
   CMultiWalleniusNCHypergeometricMoments(int32_t n, int32_t * m, double * odds, int colors, double accuracy=1.E-8) 
      : CMultiWalleniusNCHypergeometric(n, m, odds, colors, accuracy) {};
   double moments(double * mean, double * stddev, int32_t * combinations = 0);

protected:
   // functions used internally
   double loop(int32_t n, int c);                          // recursive loops
   // data
   int32_t xi[MAXCOLORS];                                  // x vector to calculate probability of
   int32_t xm[MAXCOLORS];                                  // rounded approximate mean of x[i]
   int32_t remaining[MAXCOLORS];                           // number of balls of color > c in urn
   double sx[MAXCOLORS];                                   // sum of x*f(x)
   double sxx[MAXCOLORS];                                  // sum of x^2*f(x)
   int32_t sn;                                             // number of combinations
};


/***********************************************************************
Class CFishersNCHypergeometric
***********************************************************************/

class CFishersNCHypergeometric {
   // This class contains methods for calculating the univariate Fisher's
   // noncentral hypergeometric probability function
public:
   CFishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy = 1E-8); // constructor
   double probability(int32_t x);                          // calculate probability function
   double probabilityRatio(int32_t x, int32_t x0);         // calculate probability f(x)/f(x0)
   double MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                                      // calculate approximate mean
   double variance(void);                                  // approximate variance
   int32_t mode(void);                                     // calculate mode (exact)
   double moments(double * mean, double * var);            // calculate exact mean and variance

protected:
   double lng(int32_t x);                                  // natural log of proportional function

   // parameters
   double odds;                                            // odds ratio
   double logodds;                                         // ln odds ratio
   double accuracy;                                        // accuracy
   int32_t n, m, N;                                        // Parameters
   int32_t xmin, xmax;                                     // minimum and maximum of x

   // parameters used by subfunctions
   int32_t xLast;
   double mFac, xFac;                                      // log factorials
   double scale;                                           // scale to apply to lng function
   double rsum;                                            // reciprocal sum of proportional function
   int ParametersChanged;
};


/***********************************************************************
Class CMultiFishersNCHypergeometric
***********************************************************************/

class CMultiFishersNCHypergeometric {
   // This class contains functions for calculating the multivariate
   // Fisher's noncentral hypergeometric probability function and its mean and 
   // variance. Warning: the time consumption for first call to 
   // probability or moments is proportional to the total number of
   // possible x combinations, which may be extreme!
public:
   CMultiFishersNCHypergeometric(int32_t n, int32_t * m, double * odds, int colors, double accuracy = 1E-9); // constructor
   double probability(int32_t * x);                        // calculate probability function
   void mean(double * mu);                                 // calculate approximate mean
   void variance(double * var);                            // calculate approximate variance
   double moments(double * mean, double * stddev, int32_t * combinations = 0); // calculate exact mean and variance

protected:
   double lng(int32_t * x);                                // natural log of proportional function
   void SumOfAll(void);                                    // calculates sum of proportional function for all x combinations
   double loop(int32_t n, int c);                          // recursive loops used by SumOfAll
   int32_t n, N;                                           // copy of parameters
   int32_t * m;
   double * odds;
   int colors;
   double logodds[MAXCOLORS];                              // log odds
   double mFac;                                            // sum of log m[i]!
   double scale;                                           // scale to apply to lng function
   double rsum;                                            // reciprocal sum of proportional function
   double accuracy;                                        // accuracy of calculation

   // data used by used by SumOfAll
   int32_t xi[MAXCOLORS];                                  // x vector to calculate probability of
   int32_t xm[MAXCOLORS];                                  // rounded approximate mean of x[i]
   int32_t remaining[MAXCOLORS];                           // number of balls of color > c in urn
   double sx[MAXCOLORS];                                   // sum of x*f(x) or mean
   double sxx[MAXCOLORS];                                  // sum of x^2*f(x) or variance
   int32_t sn;                                             // number of possible combinations of x
};

#endif
