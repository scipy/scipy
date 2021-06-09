/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      countpdf.c                                                   *
 *                                                                           *
 *   Count evaluations of PDF and similar functions                          *
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
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/
static char test_name[] = "CountPDF";

/*---------------------------------------------------------------------------*/
/* common variables                                                          */
/*---------------------------------------------------------------------------*/

/* counter for evaluations of PDF and similar functions */
static int counter_pdf = 0;
static int counter_dpdf = 0;
static int counter_pdpdf = 0;
static int counter_logpdf = 0;
static int counter_dlogpdf = 0;
static int counter_pdlogpdf = 0;
static int counter_cdf = 0;
static int counter_hr = 0;
static int counter_pmf = 0;

/* pointer to original functions */
static double (*cont_pdf_to_use)(double x, const struct unur_distr *distr);
static double (*cont_dpdf_to_use)(double x, const struct unur_distr *distr);
static double (*cont_logpdf_to_use)(double x, const struct unur_distr *distr);
static double (*cont_dlogpdf_to_use)(double x, const struct unur_distr *distr);
static double (*cont_cdf_to_use)(double x, const struct unur_distr *distr);
static double (*cont_hr_to_use)(double x, const struct unur_distr *distr);

static double (*discr_pmf_to_use)(int x, const struct unur_distr *distr);
static double (*discr_cdf_to_use)(int x, const struct unur_distr *distr);

static double (*cvec_pdf_to_use)(const double *x, struct unur_distr *distr);
static int (*cvec_dpdf_to_use)(double *result, const double *x, struct unur_distr *distr);
static double (*cvec_pdpdf_to_use)(const double *x, int coord, struct unur_distr *distr);
static double (*cvec_logpdf_to_use)(const double *x, struct unur_distr *distr);
static int (*cvec_dlogpdf_to_use)(double *result, const double *x, struct unur_distr *distr);
static double (*cvec_pdlogpdf_to_use)(const double *x, int coord, struct unur_distr *distr);

/*---------------------------------------------------------------------------*/
/* wrapper for functions                                                     */
/*---------------------------------------------------------------------------*/

static double cont_pdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_pdf; return cont_pdf_to_use(x,distr); }

static double cont_dpdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_dpdf; return cont_dpdf_to_use(x,distr); }

static double cont_logpdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_logpdf; return cont_logpdf_to_use(x,distr); }

static double cont_dlogpdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_dlogpdf; return cont_dlogpdf_to_use(x,distr); }

static double cont_cdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_cdf; return cont_cdf_to_use(x,distr); }

static double cont_hr_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_hr; return cont_hr_to_use(x,distr); }

/*---------------------------------------------------------------------------*/

static double discr_pmf_with_counter( int x, const struct unur_distr *distr ) {
  ++counter_pmf; return discr_pmf_to_use(x,distr); }

static double discr_cdf_with_counter( int x, const struct unur_distr *distr ) {
  ++counter_cdf; return discr_cdf_to_use(x,distr); }

/*---------------------------------------------------------------------------*/

static double cvec_pdf_with_counter( const double *x, struct unur_distr *distr ) {
  ++counter_pdf; return cvec_pdf_to_use(x,distr); }

static int cvec_dpdf_with_counter( double *result, const double *x, struct unur_distr *distr ) {
  ++counter_dpdf; return cvec_dpdf_to_use(result,x,distr); }

static double cvec_pdpdf_with_counter( const double *x, int coord, struct unur_distr *distr ) {
  ++counter_pdpdf; return cvec_pdpdf_to_use(x,coord,distr); }

static double cvec_logpdf_with_counter( const double *x, struct unur_distr *distr ) {
  ++counter_logpdf; return cvec_logpdf_to_use(x,distr); }

static int cvec_dlogpdf_with_counter( double *result, const double *x, struct unur_distr *distr ) {
  ++counter_dlogpdf; return cvec_dlogpdf_to_use(result,x,distr); }

static double cvec_pdlogpdf_with_counter( const double *x, int coord, struct unur_distr *distr ) {
  ++counter_pdlogpdf; return cvec_pdlogpdf_to_use(x,coord,distr); }

/*****************************************************************************/

int
unur_test_count_pdf( struct unur_gen *gen, int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* count evaluations of PDF and similar functions                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of evaluations of PDF and similar functions           */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  int j;
  int count;
  struct unur_gen *genclone;
  struct unur_distr *distr;

  /* check arguments */
  _unur_check_NULL(test_name,gen,-1);

  /* make a copy (clone) of the generator object */
  genclone = _unur_gen_clone(gen);

  /* get pointer to distribution object */
  if (genclone->distr_is_privatecopy) 
    distr = genclone->distr;
  else {
    /* we need a private copy of the distribution object in genclone */
    distr = genclone->distr = _unur_distr_clone( gen->distr );
    genclone->distr_is_privatecopy = TRUE; 
  }

  /* exchange pointer to PDF etc. with counting wrapper */
  switch (distr->type) {
  case UNUR_DISTR_CONT:
    cont_pdf_to_use = distr->data.cont.pdf;
    distr->data.cont.pdf = cont_pdf_with_counter;
    cont_dpdf_to_use = distr->data.cont.dpdf;
    distr->data.cont.dpdf = cont_dpdf_with_counter;
    cont_cdf_to_use = distr->data.cont.cdf;
    distr->data.cont.cdf = cont_cdf_with_counter;
    cont_hr_to_use = distr->data.cont.hr;
    distr->data.cont.hr = cont_hr_with_counter;
    if (distr->data.cont.logpdf) {
      cont_logpdf_to_use = distr->data.cont.logpdf;
      distr->data.cont.logpdf = cont_logpdf_with_counter;
    }
    if (distr->data.cont.dlogpdf) {
      cont_dlogpdf_to_use = distr->data.cont.dlogpdf;
      distr->data.cont.dlogpdf = cont_dlogpdf_with_counter;
    }
    break;
  case UNUR_DISTR_DISCR:
    discr_pmf_to_use = distr->data.discr.pmf;
    distr->data.discr.pmf = discr_pmf_with_counter;
    discr_cdf_to_use = distr->data.discr.cdf;
    distr->data.discr.cdf = discr_cdf_with_counter;
    break;
  case UNUR_DISTR_CVEC:
    cvec_pdf_to_use = distr->data.cvec.pdf;
    distr->data.cvec.pdf = cvec_pdf_with_counter;
    cvec_dpdf_to_use = distr->data.cvec.dpdf;
    distr->data.cvec.dpdf = cvec_dpdf_with_counter;
    cvec_pdpdf_to_use = distr->data.cvec.pdpdf;
    distr->data.cvec.pdpdf = cvec_pdpdf_with_counter;
    if (distr->data.cvec.logpdf) {
      cvec_logpdf_to_use = distr->data.cvec.logpdf;
      distr->data.cvec.logpdf = cvec_logpdf_with_counter;
    }
    if (distr->data.cvec.dlogpdf) {
      cvec_dlogpdf_to_use = distr->data.cvec.dlogpdf;
      distr->data.cvec.dlogpdf = cvec_dlogpdf_with_counter;
    }
    if (distr->data.cvec.pdlogpdf) {
      cvec_pdlogpdf_to_use = distr->data.cvec.pdlogpdf;
      distr->data.cvec.pdlogpdf = cvec_pdlogpdf_with_counter;
    }
    break;
  default:
    if (verbosity)
      fprintf(out,"\nCOUNT-PDF: cannot count PDF for distribution type)\n");
    _unur_free(genclone);
    return -1;
  }

  /* reset counter */
  counter_pdf = 0;
  counter_dpdf = 0;
  counter_pdpdf = 0;
  counter_logpdf = 0;
  counter_dlogpdf = 0;
  counter_pdlogpdf = 0;
  counter_cdf = 0;
  counter_hr = 0;
  counter_pmf = 0;

  /* run generator */
  switch (genclone->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    for( j=0; j<samplesize; j++ )
      _unur_sample_discr(genclone);
    break;
  case UNUR_METH_CONT:
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(genclone);
    break;
  case UNUR_METH_VEC:
    { /* we need an array for the vector */
      double *vec;
      int dim;
      dim = unur_get_dimension(genclone);
      vec = _unur_xmalloc( dim * sizeof(double) );
      for( j=0; j<samplesize; j++ )
	_unur_sample_vec(genclone,vec);
      free(vec);
    }
    break;
  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"cannot run test for method!");
    _unur_free(genclone);
    return -1;
  }

  /* total number of function calls */
  count = ( counter_pdf + counter_dpdf + counter_pdpdf 
	    + counter_logpdf + counter_dlogpdf + counter_pdlogpdf 
	    + counter_cdf + counter_hr + counter_pmf);

  /* print result */
  if (verbosity) {
    fprintf(out,  "\nCOUNT: Running Generator:\n");
    fprintf(out,  "\tfunction calls  (per generated number)\n");
    fprintf(out,  "\ttotal:   %7d  (%g)\n",count,((double)count)/((double) samplesize));
    switch (distr->type) {
    case UNUR_DISTR_CONT:
      fprintf(out,"\tPDF:     %7d  (%g)\n",counter_pdf,((double)counter_pdf)/((double) samplesize));
      fprintf(out,"\tdPDF:    %7d  (%g)\n",counter_dpdf,((double)counter_dpdf)/((double) samplesize));
      fprintf(out,"\tlogPDF:  %7d  (%g)\n",counter_logpdf,((double)counter_logpdf)/((double) samplesize));
      fprintf(out,"\tdlogPDF: %7d  (%g)\n",counter_dlogpdf,((double)counter_dlogpdf)/((double) samplesize));
      fprintf(out,"\tCDF:     %7d  (%g)\n",counter_cdf,((double)counter_cdf)/((double) samplesize));
      fprintf(out,"\tHR:      %7d  (%g)\n",counter_hr,((double)counter_hr)/((double) samplesize));
      break;
    case UNUR_DISTR_CVEC:
      fprintf(out,"\tPDF:     %7d  (%g)\n",counter_pdf,((double)counter_pdf)/((double) samplesize));
      fprintf(out,"\tdPDF:    %7d  (%g)\n",counter_dpdf,((double)counter_dpdf)/((double) samplesize));
      fprintf(out,"\tpdPDF:   %7d  (%g)\n",counter_pdpdf,((double)counter_pdpdf)/((double) samplesize));
      fprintf(out,"\tlogPDF:  %7d  (%g)\n",counter_logpdf,((double)counter_logpdf)/((double) samplesize));
      fprintf(out,"\tdlogPDF: %7d  (%g)\n",counter_dlogpdf,((double)counter_dlogpdf)/((double) samplesize));
      fprintf(out,"\tpdlogPDF:%7d  (%g)\n",counter_dlogpdf,((double)counter_dlogpdf)/((double) samplesize));
      break;
    case UNUR_DISTR_DISCR:
      fprintf(out,"\tPMF:     %7d  (%g)\n",counter_pmf,((double)counter_pmf)/((double) samplesize));
      fprintf(out,"\tCDF:     %7d  (%g)\n",counter_cdf,((double)counter_cdf)/((double) samplesize));
      break;
    }
  }

  /* remove clone */
  _unur_free(genclone);

  /* return total number of evaluations */
  return count;

} /* end of unur_test_count_pdf() */

/*****************************************************************************/

int
unur_test_par_count_pdf( struct unur_par *par, int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* count evaluations of PDF and similar functions                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par        ... pointer to parameter object                         */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of evaluations of PDF and similar functions           */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  int j;
  int count, count_total;
  struct unur_par *parclone;
  struct unur_gen *gen;
  struct unur_distr *distr;

  /* check arguments */
  _unur_check_NULL(test_name,par,-1);

  /* make a copy (clone) of the parameter object */
  parclone = _unur_par_clone(par);
  parclone->distr_is_privatecopy = TRUE;

  /* make a copy (clone) of the distribution object */
  distr = _unur_distr_clone(par->distr);

  /* set pointer to distribution object in parameter object to cloned one */
  parclone->distr = distr;

  /* exchange pointer to PDF etc. with counting wrapper */
  switch (distr->type) {
  case UNUR_DISTR_CONT:
    cont_pdf_to_use = distr->data.cont.pdf;
    distr->data.cont.pdf = cont_pdf_with_counter;
    cont_dpdf_to_use = distr->data.cont.dpdf;
    distr->data.cont.dpdf = cont_dpdf_with_counter;
    cont_cdf_to_use = distr->data.cont.cdf;
    distr->data.cont.cdf = cont_cdf_with_counter;
    cont_hr_to_use = distr->data.cont.hr;
    distr->data.cont.hr = cont_hr_with_counter;
    if (distr->data.cont.logpdf) {
      cont_logpdf_to_use = distr->data.cont.logpdf;
      distr->data.cont.logpdf = cont_logpdf_with_counter;
    }
    if (distr->data.cont.dlogpdf) {
      cont_dlogpdf_to_use = distr->data.cont.dlogpdf;
      distr->data.cont.dlogpdf = cont_dlogpdf_with_counter;
    }
    break;
  case UNUR_DISTR_DISCR:
    discr_pmf_to_use = distr->data.discr.pmf;
    distr->data.discr.pmf = discr_pmf_with_counter;
    discr_cdf_to_use = distr->data.discr.cdf;
    distr->data.discr.cdf = discr_cdf_with_counter;
    break;
  case UNUR_DISTR_CVEC:
    cvec_pdf_to_use = distr->data.cvec.pdf;
    distr->data.cvec.pdf = cvec_pdf_with_counter;
    cvec_dpdf_to_use = distr->data.cvec.dpdf;
    distr->data.cvec.dpdf = cvec_dpdf_with_counter;
    cvec_pdpdf_to_use = distr->data.cvec.pdpdf;
    distr->data.cvec.pdpdf = cvec_pdpdf_with_counter;
    if (distr->data.cvec.logpdf) {
      cvec_logpdf_to_use = distr->data.cvec.logpdf;
      distr->data.cvec.logpdf = cvec_logpdf_with_counter;
    }
    if (distr->data.cvec.dlogpdf) {
      cvec_dlogpdf_to_use = distr->data.cvec.dlogpdf;
      distr->data.cvec.dlogpdf = cvec_dlogpdf_with_counter;
    }
    if (distr->data.cvec.pdlogpdf) {
      cvec_pdlogpdf_to_use = distr->data.cvec.pdlogpdf;
      distr->data.cvec.pdlogpdf = cvec_pdlogpdf_with_counter;
    }
    break;
  default:
    if (verbosity)
      fprintf(out,"\nCOUNT-PDF: cannot count PDF for distribution type)\n");
    _unur_par_free(parclone);
    _unur_distr_free(distr);
    return -1;
  }

  /* reset counter */
  counter_pdf = 0;
  counter_dpdf = 0;
  counter_pdpdf = 0;
  counter_logpdf = 0;
  counter_dlogpdf = 0;
  counter_pdlogpdf = 0;
  counter_cdf = 0;
  counter_hr = 0;
  counter_pmf = 0;

  /* init generator */
  gen = _unur_init(parclone);

  /* total number of function calls */
  count = ( counter_pdf + counter_dpdf + counter_pdpdf
	    + counter_logpdf + counter_dlogpdf + counter_pdlogpdf
	    + counter_cdf + counter_hr + counter_pmf);
  count_total = count;

  /* print result */
  if (verbosity) {
    fprintf(out,  "\nCOUNT: Initializing Generator:\n");
    fprintf(out,  "\tfunction calls\n");
    fprintf(out,  "\ttotal:   %7d\n",count);
    switch (distr->type) {
    case UNUR_DISTR_CONT:
      fprintf(out,"\tPDF:     %7d\n",counter_pdf);
      fprintf(out,"\tdPDF:    %7d\n",counter_dpdf);
      fprintf(out,"\tlogPDF:  %7d\n",counter_logpdf);
      fprintf(out,"\tdlogPDF: %7d\n",counter_dlogpdf);
      fprintf(out,"\tCDF:     %7d\n",counter_cdf);
      fprintf(out,"\tHR:      %7d\n",counter_hr);
      break;
    case UNUR_DISTR_CVEC:
      fprintf(out,"\tPDF:     %7d\n",counter_pdf);
      fprintf(out,"\tdPDF:    %7d\n",counter_dpdf);
      fprintf(out,"\tpdPDF:   %7d\n",counter_pdpdf);
      fprintf(out,"\tlogPDF:  %7d\n",counter_logpdf);
      fprintf(out,"\tdlogPDF: %7d\n",counter_dlogpdf);
      fprintf(out,"\tpdlogPDF:%7d\n",counter_pdlogpdf);
      break;
    case UNUR_DISTR_DISCR:
      fprintf(out,"\tPMF:     %7d\n",counter_pmf);
      fprintf(out,"\tCDF:     %7d\n",counter_cdf);
      break;
    }
  }

  /* reset counter */
  counter_pdf = 0;
  counter_dpdf = 0;
  counter_pdpdf = 0;
  counter_logpdf = 0;
  counter_dlogpdf = 0;
  counter_pdlogpdf = 0;
  counter_cdf = 0;
  counter_hr = 0;
  counter_pmf = 0;

  /* run generator */
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    for( j=0; j<samplesize; j++ )
      _unur_sample_discr(gen);
    break;
  case UNUR_METH_CONT:
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(gen);
    break;
  case UNUR_METH_VEC:
    { /* we need an array for the vector */
      double *vec;
      int dim;
      dim = unur_get_dimension(gen);
      vec = _unur_xmalloc( dim * sizeof(double) );
      for( j=0; j<samplesize; j++ )
	_unur_sample_vec(gen,vec);
      free(vec);
    }
    break;
  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"cannot run test for method!");
  }

  /* total number of function calls */
  count = ( counter_pdf + counter_dpdf + counter_pdpdf
	    + counter_logpdf + counter_dlogpdf + counter_pdlogpdf
	    + counter_cdf + counter_hr + counter_pmf);
  count_total += count;

  /* print result */
  if (verbosity) {
    fprintf(out,  "\nCOUNT: Running Generator:\n");
    fprintf(out,  "\tfunction calls  (per generated number)\n");
    fprintf(out,  "\ttotal:   %7d  (%g)\n",count,((double)count)/((double) samplesize));
    switch (distr->type) {
    case UNUR_DISTR_CONT:
      fprintf(out,"\tPDF:     %7d  (%g)\n",counter_pdf,((double)counter_pdf)/((double) samplesize));
      fprintf(out,"\tdPDF:    %7d  (%g)\n",counter_dpdf,((double)counter_dpdf)/((double) samplesize));
      fprintf(out,"\tlogPDF:  %7d  (%g)\n",counter_logpdf,((double)counter_logpdf)/((double) samplesize));
      fprintf(out,"\tdlogPDF: %7d  (%g)\n",counter_dlogpdf,((double)counter_dlogpdf)/((double) samplesize));
      fprintf(out,"\tCDF:     %7d  (%g)\n",counter_cdf,((double)counter_cdf)/((double) samplesize));
      fprintf(out,"\tHR:      %7d  (%g)\n",counter_hr,((double)counter_hr)/((double) samplesize));
      break;
    case UNUR_DISTR_CVEC:
      fprintf(out,"\tPDF:     %7d  (%g)\n",counter_pdf,((double)counter_pdf)/((double) samplesize));
      fprintf(out,"\tdPDF:    %7d  (%g)\n",counter_dpdf,((double)counter_dpdf)/((double) samplesize));
      fprintf(out,"\tpdPDF:   %7d  (%g)\n",counter_pdpdf,((double)counter_pdpdf)/((double) samplesize));
      fprintf(out,"\tlogPDF:  %7d  (%g)\n",counter_logpdf,((double)counter_logpdf)/((double) samplesize));
      fprintf(out,"\tdlogPDF: %7d  (%g)\n",counter_dlogpdf,((double)counter_dlogpdf)/((double) samplesize));
      fprintf(out,"\tpdlogPDF:%7d  (%g)\n",counter_pdlogpdf,((double)counter_pdlogpdf)/((double) samplesize));
      break;
    case UNUR_DISTR_DISCR:
      fprintf(out,"\tPMF:     %7d  (%g)\n",counter_pmf,((double)counter_pmf)/((double) samplesize));
      fprintf(out,"\tCDF:     %7d  (%g)\n",counter_cdf,((double)counter_cdf)/((double) samplesize));
      break;
    }
  }

  /* free memory */
  _unur_free(gen);
  _unur_distr_free(distr);

  /* return total number of evaluations */
  return count;

} /* end of unur_test_count_pdf() */

/*---------------------------------------------------------------------------*/

