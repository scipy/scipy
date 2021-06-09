/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran_tests.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for testing routines.      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef UNURAN_TESTS_H_SEEN
#define UNURAN_TESTS_H_SEEN
/*---------------------------------------------------------------------------*/

/*
   =NODE  Testing Testing

   =UP TOP [75]

   =DESCRIPTION
     The following routines can be used to test the performance of the
     implemented generators and can be used to verify the implementions.
     They are declared in @file{unuran_tests.h} which has to be included.

   =END

*/

/*---------------------------------------------------------------------------*/
/* possible tests                                                            */
#define UNUR_TEST_ALL      (~0u)     /* make all possible tests */
#define UNUR_TEST_TIME     0x001u    /* estimate time */
#define UNUR_TEST_N_URNG   0x002u    /* count number of calls to urng */
#define UNUR_TEST_N_PDF    0x004u    /* count PDF calls */
#define UNUR_TEST_CHI2     0x008u    /* run chi^2 test for goodness of fit */
#define UNUR_TEST_SAMPLE   0x010u    /* print a sample file */

/*---------------------------------------------------------------------------*/
/* =ROUTINES */

/*---------------------------------------------------------------------------*/

void unur_run_tests( UNUR_PAR *parameters, unsigned tests, FILE *out );
/* 
   Run a battery of tests.
   The following tests are available (use @code{|} to combine these
   tests):
   @table @code
   @item UNUR_TEST_ALL
   run all possible tests.
   @item UNUR_TEST_TIME
   estimate generation times.
   @item UNUR_TEST_N_URNG
   count number of uniform random numbers
   @item UNUR_TEST_N_PDF
   count number of PDF calls
   @item UNUR_TEST_CHI2
   run chi^2 test for goodness of fit
   @item UNUR_TEST_SAMPLE
   print a small sample.
   @end table
   All these tests can be started individually (see below).
*/

/*---------------------------------------------------------------------------*/
/* particular tests                                                          */

void unur_test_printsample( UNUR_GEN *generator, int n_rows, int n_cols, FILE *out );
/* 
   Print a small sample with @var{n_rows} rows and @var{n_cols} columns.
   @var{out} is the output stream to which all results are written.
*/

UNUR_GEN *unur_test_timing( UNUR_PAR *parameters, int log10_samplesize, 
			    double *time_setup, double *time_sample,
			    int verbosity, FILE *out );
/* 
   Timing. @var{parameters} is an parameter object for which setup
   time and marginal generation times have to be measured. The results
   are written into @var{time_setup} and @var{time_sample},
   respectively. @var{log10_samplesize} is the common logarithm of the
   sample size that is used for timing. 

   If @var{verbosity} is TRUE then a small table is printed to
   output stream @var{out} with setup time, marginal generation time and
   average generation times for generating 10, 100, @dots{} random
   variates. All times are given in micro seconds and relative to 
   the generation times for the underlying uniform random number
   (using the UNIF interface) and an exponential distributed 
   random variate using the inversion method.

   The created generator object is returned.
   If a generator object could not be created successfully, then NULL
   is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.

   Notice: All timing results are subject to heavy changes. Reruning
   timings usually results in different results. Minor changes in 
   the source code can cause changes in such timings up to 25 percent.
*/

double unur_test_timing_R( UNUR_PAR *parameters, const char *distrstr, const char *methodstr,
			   double log10_samplesize, double *time_setup, double *time_marginal );
/* 
   Timing. @var{parameters} is an parameter object for which setup
   time and marginal generation times have to be measured. The results
   are written into @var{time_setup} and @var{time_marginal},
   respectively. @var{log10_samplesize} is the common logarithm of the
   sample size that is used for timing. 
   
   Alternatively, one could provide the "parameter object" using
   strings @var{distrstr} and @var{methodstr} as used for unur_makegen_ssu().

   The results are more accurate than those of function
   unur_test_timing() as the timings are computed using linear
   regression with several timings for sample size 1 and
   10^@var{log10_samplesize}.
   For each sample size total generation time (including setup) is
   measured 10 times. Since the these timings can be influenced by
   external effects (like disc sync or handling of interupts) the 2
   fastest and the 3 slowest timings are discarded.
   Intercept and slope for simple linear regression are than stored
   and @unurmath{R^2} is returned.

   In case of an error @code{-100.} is returned.

   Notice: All timing results are subject to heavy changes. Reruning
   timings usually results in different results. Minor changes in 
   the source code can cause changes in such timings up to 25 percent.
*/

double unur_test_timing_uniform( const UNUR_PAR *parameters, int log10_samplesize );
/* */

double unur_test_timing_exponential( const UNUR_PAR *parameters, int log10_samplesize );
/* 
   Marginal generation times for the underlying uniform random number
   (using the UNIF interface) and an exponential distributed 
   random variate using the inversion method. These times are used in
   unur_test_timing() to compute the relative timings results.
*/

double unur_test_timing_total( const UNUR_PAR *parameters, int samplesize, double avg_duration );
/* 
   Timing. @var{parameters} is an parameter object for which average
   times a sample of size @var{samplesize} (including setup) are
   estimated. Thus sampling is repeated and the median of these timings 
   is returned (in micro seconds). The number of iterations is computed
   automatically such that the total amount of time necessary for the
   test ist approximately @var{avg_duration} (given in seconds).
   However, for very slow generator with expensive setup time the time
   necessary for this test may be (much) larger.

   If an error occurs then @code{-1} is returned.

   Notice: All timing results are subject to heavy changes. Reruning
   timings usually results in different results. Minor changes in 
   the source code can cause changes in such timings up to 25 percent.
*/

int unur_test_count_urn( UNUR_GEN *generator, int samplesize, int verbosity, FILE *out );
/* 
   Count used uniform random numbers. It returns the total number of
   uniform random numbers required for a sample of non-uniform random
   variates of size @var{samplesize}. In case of an error @code{-1}
   is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.

   @emph{Notice:} This test uses global variables to store
   counters. Thus it is not thread save.
*/

int unur_test_count_pdf( UNUR_GEN *generator, int samplesize, int verbosity, FILE *out );
/* 
   Count evaluations of PDF and similar functions. It returns the
   total number of evaluations of all such functions required for a
   sample of non-uniform random variates of size @var{samplesize}. 
   If @var{verbosity} is TRUE then a more detailed report is printed
   to the output stream @var{out}.
   In case of an error @code{-1} is returned.
   This test is run on a copy of the given generator object.

   @emph{Notice:} The printed numbers of evaluation should be interpreted
   with care. For example, methods either use the PDF or the logPDF;
   if only the logPDF is given, but a method needs the PDF then both
   the logPDF and the PDF (a wrapper around the logPDF) are called and
   thus one call to the PDF is counted twice.

   @emph{Notice:} This test uses global variables to store function
   pointers and counters. Thus it is not thread save.
*/

int unur_test_par_count_pdf( UNUR_PAR *parameters, int samplesize, int verbosity, FILE *out );
/* 
   Same as unur_test_count_pdf() except that it is run on a parameter
   object. Thus it also prints the number of function evaluations for
   the setup. The temporary created generator object is destroyed
   before the results are returned.
*/

double unur_test_chi2( UNUR_GEN *generator, int intervals, int samplesize, int classmin,
		       int verbosity, FILE *out );
/* 
   Run a Chi^2 test with the @var{generator}. 
   The resulting p-value is returned.

   It works with discrete und continuous univariate distributions.
   For the latter the CDF of the distribution is required.

   @var{intervals} is the number of intervals that is used for
   continuous univariate distributions. @var{samplesize} is the size
   of the sample that is used for testing. If it is set to @code{0}
   then a sample of size @var{intervals}^2 is used (bounded to some
   upper bound).

   @var{classmin} is the minimum number of expected entries per
   class. If a class has to few entries then some classes are joined.

   @var{verbosity} controls the output of the routine. If it is set
   to @code{1} then the result is written to the output stream
   @var{out}. If it is set to @code{2} additionally the list of
   expected and observed data is printed.
   If it is set to @code{3} then all generated numbers are printed.
   There is no output when it is set to @code{0}.

   @emph{Notice:} For multivariate distributions the generated points
   are transformed by the inverse of the Cholesky factor of the
   covariance matrix and the mean vectors (if given for the underlying
   distribution). The marginal distributions of the transformed
   vectors are then tested against the marginal distribution given by
   a unur_distr_cvec_set_marginals() or
   unur_distr_cvec_set_marginal_array() call.
   (Notice that these marginal distributions are never set by default
   for any of the distributions provided by UNU.RAN.)
   Then the Bonferroni corrected p-value of all these tests is returned. 
   However, the test may not be performed correctly if the domain of the 
   underlying distribution is truncated by a
   unur_distr_cvec_set_domain_rect() call and the components of the
   distribution are correlated (i.e. unur_distr_cvec_set_covar() is 
   called with the non-NULL argument). Then it almost surely will fail.
*/

int unur_test_moments( UNUR_GEN *generator, double *moments, int n_moments, int samplesize,
		       int verbosity, FILE *out );
/* 
   Computes the first @var{n_moments} central moments for a sample of
   size @var{samplesize}. The result is stored into the array
   @var{moments}.
   @var{n_moments} must be an integer between @code{1} and @code{4}.
   For multivariate distributions the moments are stored consecutively 
   for each dimension and the provided @var{moments}-array must have 
   a length of at least (@var{n_moments}+1) * @var{dim}, where @var{dim}
   is the dimension of the multivariate distribution.
   The @var{m}'th moment for the @var{d}'th dimension (0<=@var{d}<@var{dim}) 
   is thus stored in the array element 
   @var{moments}[@var{d}*@var{n_moments}+@var{m}]

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

double unur_test_correlation( UNUR_GEN *generator1, UNUR_GEN *generator2,
			      int samplesize, int verbosity, FILE *out );
/* 
   Compute the correlation coefficient between streams from
   @var{generator1} and @var{generator2} for two samples of size
   @var{samplesize}.
   The resultung correlation is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

int unur_test_quartiles( UNUR_GEN *generator,
			 double *q0, double *q1, double *q2, double *q3, double *q4, 
			 int samplesize, int verbosity, FILE *out );
/* 
   Estimate quartiles of sample of size @var{samplesize}. 
   The resulting quantiles are stored in the variables @var{q}:
   @table @var
   @item q0
   minimum
   @item q1
   25%
   @item q2
   median (50%)
   @item q3
   75%
   @item q4
   maximum
   @end table

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

double unur_test_u_error( const UNUR_GEN *generator, 
			  double *max_error, double *MAE, double threshold,
			  int samplesize, int randomized, int testtails,
			  int verbosity, FILE *out );
/*
   Estimate U-error of an inversion method, i.e. 
   @unurmath{error = | CDF^{-1}(U) - U |}, by means of a simple Monte
   Carlo method.
   Maximum and mean absolute errors are stored in @var{max_error} and
   @var{MAE}, respectively.
   The particular computed U-errors should not exceed the given
   @var{threshold}. However, approximization and round-off errors
   might occasionally trigger such an event.
   Thus the function returns a penalty score. It is @code{0.} when the
   U-error never exceed the @var{threshold} value. It roughly gives the
   portion of particular test points where the U-error is too larger.
   However, each such event is weighted with 
   @unurmath{1 + 10 \times (uerror - threshold) / threshold}.

   If @var{randomized} is TRUE a pseudo-random sequence is used for
   the estimation. 

   If @var{randomized} is FALSE then the U-values are choosen
   equidistributed.
   If in addition @var{randomized} is set to TRUE then the tails of
   the distributions are tested with a more dense set of points.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.

   When the domain of the distribution is truncated then the u-error
   might be larger due to rescaling of floating point numbers. Thus
   the observed u-errors are corrected by the corresponding rescaling
   factor.

   The test also works for discrete distributions albeit with some
   restrictions: 
   It does not work correctly with truncated distributions and the
   @var{testtails} flag is ignored.
   Moreover, the value stored in @var{MAE} is rather useless.
   
   In case of an error a negative value is returned. 
*/

/* =END */

/*---------------------------------------------------------------------------*/

int unur_test_cvec_rankcorr( double *rc, UNUR_GEN *gen, int samplesize, int verbose, FILE *out );

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_TESTS_H_SEEN */
/*---------------------------------------------------------------------------*/
