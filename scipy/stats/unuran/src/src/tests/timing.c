/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      timing.c                                                     *
 *                                                                           *
 *   estimate setup and (marginal) generation time                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
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

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/cstd.h>
#include <methods/unif.h>
#include <distr/distr.h>
#include <distributions/unur_distributions.h>
#include <parser/parser.h>
#include <urng/urng.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/
/* define timer */

#if defined(HAVE_GETTIMEOFDAY) && defined(HAVE_SYS_TIME_H)
/* use gettimeofday() command. Not in ANSI C! */
#include <sys/time.h>
static struct timeval tv;
#define _unur_get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )
#else
/* use clock() command. ANSI C but less accurate */
#include <time.h>
#define _unur_get_time() ( (1.e6 * clock()) / CLOCKS_PER_SEC )
#endif

/*---------------------------------------------------------------------------*/
static char test_name[] = "Timing";
/*---------------------------------------------------------------------------*/

static double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int repeat );
/*---------------------------------------------------------------------------*/
/*  estimate average time (in micro seconds) for sampling 1 random variate   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* compare two doubles (needed for sorting) */
inline static int
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

/*---------------------------------------------------------------------------*/

struct unur_gen*
unur_test_timing( struct unur_par *par, 
		  int log10_samplesize, 
		  double *time_setup,
		  double *time_sample,
		  int verbosity,
		  FILE *out )
     /*----------------------------------------------------------------------*/
     /*  init generator and estimate setup and generation time.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par              ... pointer to parameters for generator object    */
     /*   log10_samplesize ... common log of maximal sample size             */
     /*   time_setup       ... time for setup                                */
     /*   time_sample      ... marginal generation time (i.e. for one r.n.)  */
     /*   verbosity        ... verbosity level, 0 = no output, 1 = output    */
     /*   out              ... output stream                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object.                                       */
     /*   setup time and marginal generation time are stored in              */
     /*   setup_time and marginal_time, respectively.                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  int k;
  double x;
  double *vec = NULL;
  double time_uniform, time_exponential;
  double time_start, *time_gen;
  long samples, samplesize, log10_samples;

  /* check parameter */
  _unur_check_NULL(test_name,par,NULL);
  if (log10_samplesize < 2) log10_samplesize = 2;

  /* need an array to store timings */
  time_gen = _unur_xmalloc((log10_samplesize+1) * sizeof(double));

  /* marginal generation time for one unifrom random number */
  time_uniform = unur_test_timing_uniform( par,log10_samplesize );
  /* marginal generation time for one exponential random variate */
  time_exponential = unur_test_timing_exponential( par,log10_samplesize );

  /* we need an array for the vector */
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );

  /* initialize generator (and estimate setup time) */
  time_start = _unur_get_time();
  gen = _unur_init(par);
  *time_setup = _unur_get_time();

  /* init successful ? */
  if (!gen) {
    free (time_gen);
    if (vec) free(vec);
    return NULL;
  }

  /* evaluate generation time */
  samplesize = 10;
  samples = 0;
  for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ ) {

    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( ; samples < samplesize; samples++ )
	k = unur_sample_discr(gen);
      break;
    case UNUR_METH_CONT:
    case UNUR_METH_CEMP:
      for( ; samples < samplesize; samples++ )
	x = unur_sample_cont(gen);
      break;
    case UNUR_METH_VEC:
      for( ; samples < samplesize; samples++ )
	unur_sample_vec(gen,vec);
      break;
    default: /* unknown ! */
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      free(time_gen);
      if (vec) free(vec);
      return NULL;
    }

    time_gen[log10_samples] = _unur_get_time();
    samplesize *= 10;
  }

  /* compute generation times */

  /* marginal generation time */
  *time_sample = (time_gen[log10_samplesize] - time_gen[log10_samplesize-1]) / (0.09 * samplesize);
  /* mean time per random number including setup */
  samplesize = 1;
  for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ ) {
    samplesize *= 10;
    time_gen[log10_samples] = (time_gen[log10_samples] - time_start) / samplesize;
  }
  /* setup time */
  *time_setup -= time_start;
  
  /* now print times */
  if (verbosity) {
    fprintf(out,"\nTIMING:\t\t    usec \t relative to \t relative to\n");
    fprintf(out,"\t\t\t\t uniform\t exponential\n\n");
    /* setup time */
    fprintf(out,"   setup time:\t    %#g \t %#g \t %#g\n",
	    (*time_setup),
	    (*time_setup)/time_uniform,
	    (*time_setup)/time_exponential);
    /* marginal generation time */
    fprintf(out,"   generation time: %#g \t %#g \t %#g\n",
	    (*time_sample),
	    (*time_sample)/time_uniform,
	    (*time_sample)/time_exponential);
    /* generation times */
    fprintf(out,"\n   average generation time for samplesize:\n");
    for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ )
      fprintf(out,"\t10^%ld:\t    %#g \t %#g \t %#g\n",log10_samples,
	      time_gen[log10_samples],
	      time_gen[log10_samples]/time_uniform,
	      time_gen[log10_samples]/time_exponential);
  }

  /* free memory */
  free(time_gen);
  if (vec) free(vec);

  /* return generator object */
  return gen;

} /* end of unur_test_timing() */

/*---------------------------------------------------------------------------*/

double 
unur_test_timing_R( struct unur_par *par, const char *distrstr, const char *methodstr,
		    double log10_samplesize, double *time_setup, double *time_marginal )
     /*----------------------------------------------------------------------*/
     /*  setup time and marginal generation time via linear regression.      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par              ... pointer to parameters for generator object    */
     /*   distrstr         ... distribution described by a string            */
     /*   methodstr        ... chosen method described by a string           */
     /*   log10_samplesize ... common log of maximal sample size             */
     /*   time_setup       ... time for setup                                */
     /*   time_marginal    ... marginal generation time (i.e. for one r.n.)  */
     /*                                                                      */
     /* return:                                                              */
     /*   coefficient of determination R^2                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return -100.                                                       */
     /*----------------------------------------------------------------------*/
{
  const int n_steps = 2;    /* number of sample sizes                        */
  const int n_reps = 10;    /* number of repetitions for each sample size    */

  struct unur_distr *distr_tmp = NULL; /* temporary distribution object      */
  struct unur_par   *par_tmp   = NULL; /* working copy of parameter object   */
  struct unur_gen   *gen_tmp   = NULL; /* temporary generator object         */
  struct unur_slist *mlist     = NULL; /* auxiliary table for _unur_str2par()*/

  int k;
  double x;
  double *vec = NULL;

  double time_start, *time_gen = NULL;
  int sample, rep;
  long samplesize, n;

  double sx = 0.;         /* sums for linear regression */
  double sy = 0.;
  double sxx = 0.;
  double syy = 0.;
  double sxy = 0.;

  double Rsq = -100.;     /* coefficient of determination */

  /* set initial values */
  *time_setup = -100.;
  *time_marginal = -100.;
  Rsq = -100.;

  /* check parameter */
  if (log10_samplesize < 2) log10_samplesize = 2;

  /* create parameter object (if necessary) */
  if (par == NULL) {
    distr_tmp = unur_str2distr(distrstr);
    if (distr_tmp == NULL) goto error;
    par = _unur_str2par( distr_tmp, methodstr, &mlist );
    if (par == NULL) goto error;
  }

  /* we need an array for storing timing results */
  time_gen = _unur_xmalloc( n_reps * sizeof(double) );

  /* we need an array for a random vector */
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );

  /* measure timings for various sample sizes */
  for (sample=0; sample<n_steps; sample++) {
    
    /* compute sample size */
    samplesize = (long) exp(M_LN10 * (1. + sample * (log10_samplesize - 1.) / (n_steps - 1.)));

    for (rep=0; rep<n_reps; rep++) {

      /* make a working copy of parameter object */
      par_tmp = _unur_par_clone(par);

      /* start timer */
      time_start = _unur_get_time();

      /* make generator object (init) */
      gen_tmp = _unur_init(par_tmp);
      if (!gen_tmp) goto error;  /* init failed */

      /* run generator */
      switch (gen_tmp->method & UNUR_MASK_TYPE) {
      case UNUR_METH_DISCR:
	for( n=0; n<samplesize; n++ )
	  k = unur_sample_discr(gen_tmp);
	break;
      case UNUR_METH_CONT:
	for( n=0; n<samplesize; n++ )
	  x = unur_sample_cont(gen_tmp);
	break;
      case UNUR_METH_VEC:
	for( n=0; n<samplesize; n++ )
	  unur_sample_vec(gen_tmp,vec);
	break;
      default: /* unknown ! */
	_unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      }

      /* stop timer */
      time_gen[rep]= _unur_get_time() - time_start;
      
      /* destroy generator object */
      unur_free(gen_tmp);
    }
    
    /* sort timings */
    qsort( time_gen, (size_t)n_reps, sizeof(double), compare_doubles);

    /* update sums for linear regression */
    for (rep=2; rep<n_reps-3; rep++) {
      sx += samplesize;
      sy += time_gen[rep];
      sxx += ((double)samplesize) * ((double)samplesize);
      syy += time_gen[rep] * time_gen[rep];
      sxy += ((double)samplesize) * time_gen[rep];
    }
  }

  /* compute simple linear regression */
  n = n_steps * (n_reps-5);
  *time_marginal = ( n*sxy - sx*sy ) / ( n*sxx - sx*sx );
  *time_setup = sy/n - *time_marginal * sx/n;
  Rsq = ( n*sxy - sx*sy ) / sqrt( (n*sxx - sx*sx) * (n*syy - sy*sy) );
  
 error:
  /* free memory */
  if (distr_tmp) unur_distr_free(distr_tmp);
  if (par)       _unur_par_free(par);
  if (mlist)     _unur_slist_free(mlist);
  if (time_gen)  free(time_gen);
  if (vec)       free(vec);

  /* return result */
  return Rsq;

} /* end of unur_test_timing_R() */

/*---------------------------------------------------------------------------*/

double 
unur_test_timing_total( const UNUR_PAR *par, int samplesize, double avg_duration )
     /*----------------------------------------------------------------------*/
     /*  estimate average time (in micro seconds) for generating a sample    */
     /*  of size `samplesize' (including setup) are estimated.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameters for generator object        */
     /*   samplesize   ... sample size                                       */
     /*   avg_duration ... average total time (in seconds) for running test  */
     /*                                                                      */
     /* return:                                                              */
     /*   median of several runs for generating sample (in micro seconds)    */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
     /*----------------------------------------------------------------------*/
{
  double time_pilot, time_result;
  int size_pilot, size_result;
  int repeat_pilot, repeat_result;
  double time_2nd, d, k;

  /* check parameter */
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0) return -1.;

  /* rescale maximal duration from seconds to micro seconds */
  avg_duration = (avg_duration < 1.e-3) ? 1000. : 1.e6 * avg_duration;

  /* pilot study */
  repeat_pilot = 11 - (int)(log((double)samplesize)/M_LN2);
  if (repeat_pilot<1) repeat_pilot = 1;

  size_pilot = _unur_min(samplesize,1000);

  time_pilot = unur_test_timing_total_run(par, size_pilot, repeat_pilot);
  if (time_pilot < 0) return -1.;   /* init failed */

  if (samplesize > 1000) {
    /* make second pilot study with double sample size */
    time_2nd = unur_test_timing_total_run(par, 2*size_pilot, repeat_pilot);
    if (time_2nd < 0) return -1.;

    /* estimate time for given sample size */
    d = 2*time_pilot - time_2nd;
    if (d<0.) d=0.;
    k = (time_2nd - time_pilot)/size_pilot;
    if (k<=0.) k = time_pilot/size_pilot;
    time_pilot = d + samplesize * k;
  }
  else {
    /* this is not required, but it prevents an error in case of a programming bug */
    d = 0;
    k = time_pilot / size_pilot;
  }
  
  /* now run timing test */

  repeat_result = (int) (avg_duration / time_pilot);
  if (repeat_result > 1000) repeat_result = 1000;
  /* there is no need for more than 1000 repetitions */

  size_result = samplesize;

  if (repeat_result >= 1) {
    repeat_result = _unur_max(4,repeat_result);
    if (repeat_result <= repeat_pilot && size_result == size_pilot) {
      /* there is no need to run this test again */
      time_result = time_pilot;
    }
    else {
      time_result =  unur_test_timing_total_run(par,size_result,repeat_result);
    }
  }
  else {
    /* do not generate the full sample */
    repeat_result = 4;
    size_result = (int) ((avg_duration - d)/k);
    size_result /= 2;
    time_result =  unur_test_timing_total_run(par,size_result,repeat_result);
    time_2nd =  unur_test_timing_total_run(par,2*size_result,repeat_result);
    /* estimate time from shorter sample sizes */
    d = 2*time_result - time_2nd;
    if (d<0.) d=0.;
    k = (time_2nd - time_result)/size_result;
    if (k<=0.) k = time_result/size_result;
    time_result = d + samplesize * k;
  }      

  /* o.k. */
  return time_result;

} /* end of unur_test_timing_total() */

/*---------------------------------------------------------------------------*/

double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int n_repeat )
     /*----------------------------------------------------------------------*/
     /*  estimate average time (in micro seconds) for sampling               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par        ... pointer to parameters for generator object          */
     /*   samplesize ... sample size                                         */
     /*   n_repeat   ... number of samples (repetitions of sampling)         */
     /*                                                                      */
     /* return:                                                              */
     /*   total time in micro seconds                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par_tmp;    /* temporary working copy of parameter object */
  struct unur_gen *gen_tmp;    /* temporary generator object                 */
  double *time;
  int n, rep;
  double time_total;           /* total time for sampling                    */
  double time_start;
  int i,k;
  double x;
  double *vec = NULL;

  /* check parameter */
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0 || n_repeat < 1) 
    return -1.;

  /* we need an array for storing timing results */
  time = _unur_xmalloc( n_repeat * sizeof(double) );

  /* we need an array for a random vector */
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );

  /* make samples */
  for (rep = 0; rep < n_repeat; rep++) {

    /* make a working copy of parameter object */
    par_tmp = _unur_par_clone(par);

    /* start timer */
    time_start = _unur_get_time();

    /* make generator object (init) */
    gen_tmp = _unur_init(par_tmp);
    if (!gen_tmp) {  /* init failed */
      if (vec) free(vec);
      free(time);
      /* wait a little bit till writing log entry is completed */
      for (x=0,i=0; i<100000; i++) x+=i;  
      return -1.;
    }

    /* run generator */
    switch (gen_tmp->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( n=0; n<samplesize; n++ )
	k = unur_sample_discr(gen_tmp);
      break;
    case UNUR_METH_CONT:
      for( n=0; n<samplesize; n++ )
	x = unur_sample_cont(gen_tmp);
      break;
    case UNUR_METH_VEC:
      for( n=0; n<samplesize; n++ )
	unur_sample_vec(gen_tmp,vec);
      break;
    default: /* unknown ! */
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }

    /* stop timer */
    time[rep]= _unur_get_time() - time_start;

    /* destroy generator object */
    unur_free(gen_tmp);
  }

  /* compute median */
  qsort( time, (size_t)n_repeat, sizeof(double), compare_doubles);
  time_total = time[n_repeat/2];

  /* free memory */
  if (vec) free(vec);
  free(time);

  return time_total;
} /* end of unur_test_timing_total_run() */

/*---------------------------------------------------------------------------*/

double
unur_test_timing_uniform( const struct unur_par *par, int log10_samplesize )
     /*----------------------------------------------------------------------*/
     /*  estimate generation time for URNG using UNURAN wrapper.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to paramters for building generator object      */
     /*   log10_samplesize ... common log of maximal sample size             */
     /*                                                                      */
     /* return:                                                              */
     /*   mean generation time                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define TIMING_REPETITIONS (21)

  struct unur_gen *gen_urng;
  static double uniform_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;

  if (uniform_time <= 0.) {  
    /* unknown. have to estimate time first */
    
    /* sample size */
    int samplesize = 1;
    for( j=0; j<log10_samplesize; j++ )
      samplesize *= 10;

    /* make generator object for uniform generator */
    gen_urng = unur_init( unur_unif_new(NULL) );
    _unur_check_NULL( test_name,gen_urng,-1. );
    unur_chg_urng(gen_urng,par->urng);

    /* evaluate marginal generation times */
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(gen_urng);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }

    /* compute median */
    qsort( time, (size_t)TIMING_REPETITIONS, sizeof(double), compare_doubles);

    /* store marginal generation time for uniforms */
    uniform_time = time[TIMING_REPETITIONS/2];

    /* free generator object for uniform random number generator */
    _unur_free(gen_urng);

  }

  return uniform_time;

#undef TIMING_REPETITIONS

} /* end of unur_test_timing_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_test_timing_exponential( const struct unur_par *par, int log10_samplesize )
     /*----------------------------------------------------------------------*/
     /*  estimate generation time for URNG using UNURAN wrapper.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to paramters for building generator object      */
     /*   log10_samplesize ... common log of maximal sample size             */
     /*                                                                      */
     /* return:                                                              */
     /*   mean generation time                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define TIMING_REPETITIONS (21)

  struct unur_distr *unit_distr;
  struct unur_par   *unit_par;
  struct unur_gen   *unit_gen;
  static double exponential_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;

  if (exponential_time <= 0.) {  
    /* unknown. have to estimate time first */
    
    /* sample size */
    int samplesize = 1;
    for( j=0; j<log10_samplesize; j++ )
      samplesize *= 10;

    /* make generator object for uniform generator */
    unit_distr = unur_distr_exponential(NULL,0);
    unit_par = unur_cstd_new(unit_distr);
    unur_cstd_set_variant(unit_par,UNUR_STDGEN_INVERSION);
    unit_gen = unur_init(unit_par); 
    _unur_check_NULL( test_name,unit_gen,-1. );
    unur_chg_urng(unit_gen,par->urng);

    /* evaluate marginal generation times */
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(unit_gen);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }

    /* compute median */
    qsort( time, (size_t)TIMING_REPETITIONS, sizeof(double), compare_doubles);

    /* store marginal generation time for uniforms */
    exponential_time = time[TIMING_REPETITIONS/2];

    /* free generator object for uniform random number generator */
    unur_distr_free(unit_distr);
    unur_free(unit_gen);

  }

  return exponential_time;

#undef TIMING_REPETITIONS

} /* end of unur_test_timing_exponential() */

/*---------------------------------------------------------------------------*/
