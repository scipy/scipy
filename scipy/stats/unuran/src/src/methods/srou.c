/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      srou.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    simple universal method (ratio-of-uniforms method)           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and mode of a T_{-1/2}-concave distribution                *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      mode of the density                                                  *
 *      area below PDF                                                       *
 *   OPTIONAL:                                                               *
 *      CDF at mode                                                          *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Leydold J. (2001): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions,                        *
 *       ACM Trans. Math. Software 27(1), pp. 66--82.                        *
 *                                                                           *
 *   [2] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *   [3] Leydold J. (2002): Short universal generators via generalized       *
 *       ratio-of-uniforms method, Preprint.                                 *
 *                                                                           *
 *   [4] Wakefield J. C., Gelfand A. E., and Smith A. F. M. (1991):          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method, Statist. Comput. 1(2), pp. 129--133.                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [2] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (V,U) is uniformly distributed in                                   *
 *       A = {(v,u): 0 < u <= sqrt(g(v/u)), x_0 < v/u < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (V,U) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 * The implemented algorithm uses the fact, that for many distribtions,      *
 * A is convex. Then we easily can construct an enveloping rectangle.        *
 * Define                                                                    *
 *    R = {(v,u):  v_l <= v <= v_r, 0 <= u <= u_m},                          *
 *    Q = {(v,u): -v_m <= v <= v_m, 0 <= u <= u_m},                          *
 * where                                                                     *
 *    u_m = sqrt(f(mode)), v_m = (\int f dx) / u_m                           *
 *    v_l = -F(\mode) v_m, v_r = (1-F(mode)) v_m                             *
 * Then                                                                      *
 *    A subset R subset Q                                                    *
 *                                                                           *
 * Thus we can use R to generate whenever the CDF F(mode) at the mode        *
 * is known, and Q otherwise.                                                *
 * Notice, that the rection constant is 2 in the first case and 4 and the    *
 * latter.                                                                   *
 *                                                                           *
 * If F(mode) it known, it is even possible to get an universal squeeze      *
 * (see [1] for details). However its usage is only recommended when         *
 * the PDF is (very) expensive.                                              *
 *                                                                           *
 * When F(mode) is not known the mirror principle can be used. i.e., make    *
 * an enveloping rectangle for f(x)+f(-x). It reduces the rejection constant *
 * to 2 * sqrt(2) at the expense of more evaluations of the PDF.             *
 * Its usage is only recommended when the generation time for the underlying *
 * uniform prng is extreme large.                                            *
 *                                                                           *
 * Distributions with a convex set A are characterized by the following      *
 * theorem that shows a connection to transformed density rejection TDR.     *
 *                                                                           *
 * THEOREM:                                                                  *
 *    A is convex if and only if g is T-concave with transformation          *
 *    T(x) = -1/sqrt(x), i.e., -1/sqrt(g(x)) is a concave function.          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The theory for r!=1 is described in [3].                                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "srou.h"
#include "srou_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define SROU_VARFLAG_VERIFY   0x002u   /* run verify mode                    */
#define SROU_VARFLAG_SQUEEZE  0x004u   /* use universal squeeze if possible  */
#define SROU_VARFLAG_MIRROR   0x008u   /* use mirror principle               */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define SROU_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define SROU_SET_R            0x001u   /* parameter r for power transform.   */
#define SROU_SET_CDFMODE      0x002u   /* CDF at mode is known               */
#define SROU_SET_PDFMODE      0x004u   /* PDF at mode is set                 */

/*---------------------------------------------------------------------------*/

#define GENTYPE "SROU"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_srou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_srou_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_srou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_srou_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_srou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_srou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_srou_sample( struct unur_gen *gen );
static double _unur_srou_sample_mirror( struct unur_gen *gen );
static double _unur_srou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator (case r=1).                                         */
/*---------------------------------------------------------------------------*/

static double _unur_gsrou_sample( struct unur_gen *gen );
static double _unur_gsrou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator (case r>1).                                         */
/*---------------------------------------------------------------------------*/

static int _unur_srou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute universal bounding rectangle (case r=1).                          */
/*---------------------------------------------------------------------------*/

static int _unur_gsrou_envelope( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute parameters for universal bounding envelope (case r>1).            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_srou_debug_init( const struct unur_gen *gen, int is_reinit );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_srou_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_srou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_srou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CONT *
_unur_srou_getSAMPLE( struct unur_gen *gen )
{
  if (gen->variant & SROU_VARFLAG_VERIFY)
    return (gen->set & SROU_SET_R) ? _unur_gsrou_sample_check : _unur_srou_sample_check;
  else {
    if (gen->set & SROU_SET_R)
      return _unur_gsrou_sample;
    else
      return (gen->variant & SROU_VARFLAG_MIRROR) ? _unur_srou_sample_mirror : _unur_srou_sample;
  }
} /* end of _unur_srou_getSAMPLE() */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2    (M_SQRT2)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_srou_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_srou_par) );
  COOKIE_SET(par,CK_SROU_PAR);

  /* copy input */
  par->distr    = distr;           /* pointer to distribution object         */

  /* set default values */
  PAR->r         = 1.;             /* parameter for power transformation     */
  PAR->Fmode     = -1.;            /* CDF at mode (unknown yet   )           */
  PAR->um        = -1.;            /* (square) root of PDF at mode (unknown) */

  par->method   = UNUR_METH_SROU;  /* method and default variant             */
  par->variant  = 0u;              /* default variant                        */
  par->set      = 0u;              /* inidicate default parameters           */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_srou_init;

  return par;

} /* end of unur_srou_new() */

/*****************************************************************************/

int
unur_srou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* set parameter r for power transformation                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   r    ... parameter r                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* check new parameter for generator */
  if (r < 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r < 1");
    return UNUR_ERR_PAR_SET;
  }
  
  if (_unur_isone(r)) {
    /* simple version, same as R is not set */ 
    PAR->r = r;
    par->set &= ~SROU_SET_R;
  }
  else {
    /* for computational reasons r should be at least 1.01 */
    if (r<1.01) r = 1.01;
    PAR->r = r;
    par->set |= SROU_SET_R;
  }

  /* we have to reset the marker for the PDF at the mode */
  par->set &= ~SROU_SET_PDFMODE;

  return UNUR_SUCCESS;

} /* end of unur_srou_set_r() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_cdfatmode( struct unur_par *par, double Fmode )
     /*----------------------------------------------------------------------*/
     /* set value of cdf at mode                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   Fmode ... cdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->Fmode = Fmode;

  /* changelog */
  par->set |= SROU_SET_CDFMODE;

  return UNUR_SUCCESS;

} /* end of unur_srou_set_cdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_pdfatmode( UNUR_PAR *par, double fmode )
     /*----------------------------------------------------------------------*/
     /* Set pdf at mode. if set the PDF at the mode is never changed.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_isfinite(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return UNUR_ERR_PAR_SET;
  }

  /* store date ((square) root of fmode) */
  PAR->um = (par->set & SROU_SET_R) ? pow(fmode,1./(PAR->r+1.)) : sqrt(fmode);

  /* changelog */
  par->set |= SROU_SET_PDFMODE;

  return UNUR_SUCCESS;

} /* end of unur_srou_set_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | SROU_VARFLAG_VERIFY) : (par->variant & (~SROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_srou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= SROU_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~SROU_VARFLAG_VERIFY;

  SAMPLE = _unur_srou_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_srou_chg_verify() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_usesqueeze( struct unur_par *par, int usesqueeze )
     /*----------------------------------------------------------------------*/
     /* set flag for using universal squeeze (default: off)                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   usesqueeze ... 0 = no squeeze,  !0 = use squeeze                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no squeeze is the default                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* we use a bit in variant */
  par->variant = (usesqueeze) ? (par->variant | SROU_VARFLAG_SQUEEZE) : (par->variant & (~SROU_VARFLAG_SQUEEZE));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_srou_set_usesqueeze() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_usemirror( struct unur_par *par, int usemirror )
     /*----------------------------------------------------------------------*/
     /* set flag for using mirror principle (default: off)                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   usemirror ... 0 = no mirror princ.,  !0 = use mirror principle     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   do not use mirror principle is the default                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );

  /* we use a bit in variant */
  par->variant = (usemirror) ? (par->variant | SROU_VARFLAG_MIRROR) : (par->variant & (~SROU_VARFLAG_MIRROR));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_srou_set_usemirror() */

/*****************************************************************************/

int
unur_srou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
     /*----------------------------------------------------------------------*/
     /* change value of cdf at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   Fmode ... cdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  
  /* copy parameters */
  GEN->Fmode = Fmode;

  /* changelog */
  gen->set |= SROU_SET_CDFMODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_srou_chg_cdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_pdfatmode( struct unur_gen *gen, double fmode )
     /*----------------------------------------------------------------------*/
     /* change value of PDF at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_isfinite(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return UNUR_ERR_PAR_SET;
  }

  /* store date ((square) root of fmode) */
  GEN->um = (gen->set & SROU_SET_R) ? pow(fmode,1./(GEN->r+1.)) : sqrt(fmode);

  /* changelog */
  gen->set |= SROU_SET_PDFMODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_srou_chg_pdfatmode() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_srou_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  int rcode;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_SROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SROU_PAR,NULL);

  /* there are no squeezes and no mirrir principle when r is changed */
  if (par->set & SROU_SET_R) {
    par->variant &= ~SROU_VARFLAG_MIRROR;
    par->variant &= ~SROU_VARFLAG_SQUEEZE;
  }

  if (par->set & SROU_SET_CDFMODE)
    /* cdf at mode known -->
       thus it does not make sense to use the mirror principle */
    par->variant &= ~SROU_VARFLAG_MIRROR;
  else
    /* cdf at mode unknown -->
       thus we cannot use universal squeeze */
    par->variant &= ~SROU_VARFLAG_SQUEEZE;

  /* create a new empty generator object */
  gen = _unur_srou_create(par);

  /* free parameters */
  _unur_par_free(par);

  if (!gen) return NULL;

  /* check parameters */
  if (_unur_srou_check_par(gen) != UNUR_SUCCESS) {
    _unur_srou_free(gen); return NULL;
  }

  /* compute universal bounding envelope */
  if (gen->set & SROU_SET_R)
    rcode = _unur_gsrou_envelope( gen );
  else
    rcode = _unur_srou_rectangle( gen );

  if (rcode!=UNUR_SUCCESS) {
    /* error */
    _unur_srou_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_srou_debug_init(gen, FALSE);
#endif

  return gen;

} /* end of _unur_srou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_srou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int rcode;

  /* check parameters */
  if ( (rcode = _unur_srou_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* compute universal bounding envelope */
  if (gen->set & SROU_SET_R)
    rcode = _unur_gsrou_envelope( gen );
  else
    rcode = _unur_srou_rectangle( gen );

  /* (re)set sampling routine */
  SAMPLE = _unur_srou_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
  if (gen->debug & SROU_DEBUG_REINIT) _unur_srou_debug_init(gen,TRUE);
#endif

  return rcode;
} /* end of _unur_srou_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_srou_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_srou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_SROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_srou_getSAMPLE(gen);
  gen->destroy = _unur_srou_free;
  gen->clone = _unur_srou_clone;
  gen->reinit = _unur_srou_reinit;

  /* copy some parameters into generator object */
  GEN->r     = PAR->r;                /* parameter for power transformation    */
  GEN->Fmode = PAR->Fmode;            /* CDF at mode                           */
  GEN->um    = PAR->um;               /* square root of PDF at mode            */

  /* initialize variables */
  GEN->vl = GEN->vr = 0.;
  GEN->xl = GEN->xr = 0.;
  GEN->p = 0.;
  GEN->a = GEN->b = 0.;
  GEN->log_ab = 0.;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_srou_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_srou_create() */

/*---------------------------------------------------------------------------*/

int
_unur_srou_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check for required data: mode */
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)");
    if (unur_distr_cont_upd_mode(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* check for required data: area */
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA)) {
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or CDF at mode");
    DISTR.mode = _unur_max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = _unur_min(DISTR.mode,DISTR.BD_RIGHT);
  }

  return UNUR_SUCCESS;
} /* end of _unur_srou_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_srou_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE  ((struct unur_srou_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_srou_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_srou_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_SROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_srou_free() */

/*****************************************************************************/

double
_unur_srou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X,x,xx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
    U *= GEN->um;
    V = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);

    /* ratio */
    X = V/U;

    /* compute x */
    x = X + DISTR.mode;

    /* inside domain ? */
    if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
      continue;

    /* evaluate squeeze */
    if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	 (X >= GEN->xl) && 
	 (X <= GEN->xr ) && 
	 (U < GEN->um) ) {
      xx = V / (GEN->um - U);
      if ( (xx >= GEN->xl) && (xx <= GEN->xr ) )
	return x;
    }

    /* accept or reject */
    if (U*U <= PDF(x))
      return x;
  }

} /* end of _unur_srou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_srou_sample_mirror( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator, use mirror principle                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X,x,fx,fnx,uu;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
    U *= GEN->um * SQRT2;
    V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN->vr;
    /* vr = vm when the CDF at the mode is not known */

    /* ratio */
    X = V/U;

    /* evaluate PDF */
    x = X + DISTR.mode;
    fx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x);
    uu = U * U;

    /* accept or reject */
    if (uu <= fx)
      return x;

    /* try mirrored PDF */
    x = -X + DISTR.mode;
    fnx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x);
    if (uu <= fx + fnx)
      return x;
  }

} /* end of _unur_srou_sample_mirror() */

/*---------------------------------------------------------------------------*/

double
_unur_srou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,uu,V,X,x,nx,fx,sfx,fnx,xfx,xfnx,xx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,INFINITY);

  if (gen->variant & SROU_VARFLAG_MIRROR) {
    /* use mirror principle */

    while (1) {
      /* generate point uniformly on rectangle */
      while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
      U *= GEN->um * SQRT2;
      V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN->vr;
      /* vr = vm when the CDF at the mode is not known */	

      /* ratio */
      X = V/U;

      /* x values */
      x = X + DISTR.mode;
      nx = -X + DISTR.mode;

      /* evaluate PDF */
      fx  = (x  < DISTR.BD_LEFT || x  > DISTR.BD_RIGHT) ? 0. : PDF(x);
      fnx = (nx < DISTR.BD_LEFT || nx > DISTR.BD_RIGHT) ? 0. : PDF(nx);
      uu = U * U;

      /* check hat */
      xfx  = (x  - DISTR.mode) * sqrt(fx);
      xfnx = (nx - DISTR.mode) * sqrt(fnx);

      if ( ((2.+4.*DBL_EPSILON) * GEN->um*GEN->um < fx + fnx)    /* avoid roundoff error with FP registers */
	   || (xfx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfx > (1.+UNUR_EPSILON) * GEN->vr)
	   || (xfnx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfnx > (1.+UNUR_EPSILON) * GEN->vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");

      /* accept or reject */
      if (uu <= fx)
	return x;
      
      /* try mirrored PDF */
      if (uu <= fx + fnx)
	return nx;
    }
  }

  else { /* do not use mirror principle */

    while (1) {
      /* generate point uniformly on rectangle */
      while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
      U *= GEN->um;
      V = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);

      /* ratio */
      X = V/U;

      /* compute x */
      x = X + DISTR.mode;

      /* inside domain ? */
      if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	continue;

      /* evaluate PDF */
      fx = PDF(x);

      /* the point on the boundary of the region of acceptance
	 in direction X = V/U has the coordinates
	 ( X * sqrt(fx), sqrt(fx) ). */
      sfx = sqrt(fx);
      xfx = X * sfx;

      /* check hat */
      if ( ( sfx > (1.+DBL_EPSILON) * GEN->um )   /* avoid roundoff error with FP registers */
	   || (xfx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfx > (1.+UNUR_EPSILON) * GEN->vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");

      /* evaluate and check squeeze */
      if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	   (X >= GEN->xl) && 
	   (X <= GEN->xr ) && 
	   (U < GEN->um) ) {

	/* check squeeze */
	xx = xfx / (GEN->um - sfx);
	if ( (xx > (1.-UNUR_EPSILON) * GEN->xl) &&
	     (xx < (1.-UNUR_EPSILON) * GEN->xr) )
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");

	/* squeeze acceptance */
	xx = V / (GEN->um - U);
	if ( (xx >= GEN->xl) && (xx <= GEN->xr ) )
	  return x;
      }
      
      /* accept or reject */
      if (U*U <= PDF(x))
	return x;
    }
  }

} /* end of _unur_srou_sample_check() */

/*****************************************************************************/

double
_unur_gsrou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,Ur,V,W,X,Z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,INFINITY);

  while (1) {
    W = GEN->log_ab *_unur_call_urng(gen->urng);
    Z = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
    U = (exp(-W)-1.) * GEN->a/GEN->b;
    V = -Z/(GEN->a + GEN->b*U);
    U *= GEN->um;
    Ur = pow(U,GEN->r);
    X = V/Ur + DISTR.mode;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* accept or reject */
    if (Ur*U <= PDF(X))
      return X;
  }

} /* end of _unur_gsrou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_gsrou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,Ur,V,W,X,x,Z;
  double fx,uf,vf,vhl,vhr;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,INFINITY);

  while (1) {

    W = GEN->log_ab *_unur_call_urng(gen->urng);
    Z = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
    U = (exp(-W)-1.) * GEN->a/GEN->b;
    V = -Z/(GEN->a + GEN->b*U);
    U *= GEN->um;
    Ur = pow(U,GEN->r);
    X = V/Ur;
    
    /* compute x */
    x = X + DISTR.mode;
    /* inside domain ? */
    if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
      continue;

    /* evaluate density */
    fx = PDF(x);

    /* the point on the boundary of the region of acceptance
       in direction X = V/U^r has the coordinates
       ( (x-mode) * (fx)^(r/(r+1)), sqrt[r+1](fx) ). */
    uf = pow(fx,1./(GEN->r+1));
    vf = X * pow(fx,GEN->r/(GEN->r+1.));

    /* the corresponding point on boundary of the enveloping region */
    /* with same u-coordinate.                                      */
    vhl = - GEN->vl /(GEN->a + GEN->b*(uf/GEN->um));
    vhr = - GEN->vr /(GEN->a + GEN->b*(uf/GEN->um));

    /* check hat */
    if ( ( uf > (1.+DBL_EPSILON) * GEN->um )   /* avoid roundoff error with FP registers */
	 || (vf < (1.+UNUR_EPSILON) * vhl )
	 || (vf > (1.+UNUR_EPSILON) * vhr ) )
      {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      }

    /* accept or reject */
    if (Ur*U <= fx)
      return x;
  }

} /* end of _unur_gsrou_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_srou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double vm, fm;             /* width of rectangle, PDF at mode              */

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_SROU_GEN, UNUR_ERR_COOKIE );

  /* compute PDF at mode (if not given by user) */
  if (!(gen->set & SROU_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    if (!_unur_isfinite(fm)) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
      return UNUR_ERR_PAR_SET;
    }
    GEN->um = sqrt(fm);    /* height of rectangle */
  }

  /* width of rectangle */
  vm = DISTR.area / GEN->um;

  if (gen->set & SROU_SET_CDFMODE) {
    /* cdf at mode known */
    GEN->vl = -GEN->Fmode * vm;
    GEN->vr = vm + GEN->vl;
    GEN->xl = GEN->vl/GEN->um;
    GEN->xr = GEN->vr/GEN->um;
  }
  else {
    /* cdf at mode unknown */
    GEN->vl = -vm;
    GEN->vr = vm;
    GEN->xl = GEN->vl/GEN->um;
    GEN->xr = GEN->vr/GEN->um;
    /* we cannot use universal squeeze */
    gen->variant &= ~SROU_VARFLAG_SQUEEZE;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_srou_rectangle() */

/*---------------------------------------------------------------------------*/

int
_unur_gsrou_envelope( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute parameters for universal bounding envelope.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double fm;             /* PDF at mode               */
  double vm;             /* maximal width of envelope */
  double pr;             /* p^r                       */

  double p;              /* short cuts                */
  double r = GEN->r;

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_SROU_GEN, UNUR_ERR_COOKIE );

  if (!(gen->set & SROU_SET_PDFMODE)) {
    /* compute PDF at mode */
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    if (!_unur_isfinite(fm)) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
      return UNUR_ERR_PAR_SET;
    }
    GEN->um = pow(fm,1./(r+1.));    /* height of envelope */
  }

  /* maximal width of envelope */
  vm = DISTR.area / (GEN->r*GEN->um);

  if (gen->set & SROU_SET_CDFMODE) {
    /* cdf at mode known */
    GEN->vl = -GEN->Fmode * vm;
    GEN->vr = vm + GEN->vl;
  }
  else {
    /* cdf at mode unknown */
    GEN->vl = -vm;
    GEN->vr = vm;
  }

  /* construction point for bounding curve */
  GEN->p = p = 1. - 2.187/pow(r + 5 - 1.28/r, 0.9460 );
  pr = pow(p,r);

  /* parameters for bounding envelope */
  GEN->b = (1. - r * pr/p + (r-1.)*pr) / ((pr-1.)*(pr-1));
  GEN->a = -(p-1.)/(pr-1.) - p * GEN->b;
  GEN->log_ab = log(GEN->a/(GEN->a+GEN->b));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_gsrou_envelope() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_srou_debug_init( const struct unur_gen *gen, int is_reinit )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   is_reinit ... if TRUE the generator has been reinitialized         */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_SROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
    fprintf(LOG,"%s: method  = srou (simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: reinit!\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  if (gen->set & SROU_SET_R) {
    fprintf(LOG,"%s: Generalized version: r = %g\n",gen->genid,GEN->r);
    fprintf(LOG,"%s:\n",gen->genid);

    fprintf(LOG,"%s: sampling routine = _unur_gsrou_sample",gen->genid);
    if (gen->variant & SROU_VARFLAG_VERIFY)
      fprintf(LOG,"_check");
    fprintf(LOG,"()\n%s:\n",gen->genid);
  }
  else {
    fprintf(LOG,"%s: Simple version (r = 1)  [default]\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);

    fprintf(LOG,"%s: sampling routine = _unur_srou_sample",gen->genid);
    if (gen->variant & SROU_VARFLAG_VERIFY)
      fprintf(LOG,"_check");
    else if (gen->variant & SROU_VARFLAG_MIRROR)
      fprintf(LOG,"_mirror");
    fprintf(LOG,"()\n%s:\n",gen->genid);
  }

  if (gen->set & SROU_SET_CDFMODE)
    fprintf(LOG,"%s: F(mode) = %g\n",gen->genid,GEN->Fmode);
  else
    fprintf(LOG,"%s: F(mode) unknown\n",gen->genid);

  if (gen->variant & SROU_VARFLAG_SQUEEZE)
    fprintf(LOG,"%s: use universal squeeze\n",gen->genid);
  else
    fprintf(LOG,"%s: no (universal) squeeze\n",gen->genid);

  if (gen->variant & SROU_VARFLAG_MIRROR)
    fprintf(LOG,"%s: use mirror principle\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);

  if (gen->set & SROU_SET_R) {
    fprintf(LOG,"%s: Enveloping region:\n",gen->genid);
    fprintf(LOG,"%s:    um = %g\n",gen->genid,GEN->um);
    fprintf(LOG,"%s:    vl = %g\n",gen->genid,GEN->vl);
    fprintf(LOG,"%s:    vr = %g\n",gen->genid,GEN->vr);
    fprintf(LOG,"%s:    p  = %g\n",gen->genid,GEN->p);
    fprintf(LOG,"%s:    a  = %g\n",gen->genid,GEN->a);
    fprintf(LOG,"%s:    b  = %g\n",gen->genid,GEN->b);
  }
  else {
    fprintf(LOG,"%s: Rectangle:\n",gen->genid);
    fprintf(LOG,"%s:    left upper point  = (%g,%g)\n",gen->genid,GEN->vl,GEN->um);
    fprintf(LOG,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN->vr,GEN->um);
  }

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_srou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_srou_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  double h_area, rc;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %g   %s\n", DISTR.mode,
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = %g\n", DISTR.area);
  if (gen->set & SROU_SET_CDFMODE)
    _unur_string_append(info,"   F(mode)   = %g\n", GEN->Fmode); 
  else
    _unur_string_append(info,"   F(mode)   = [unknown]\n"); 

  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide the \"mode\"");
  }
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: SROU (Simple Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g  %s\n", GEN->r,
		      (gen->set & SROU_SET_R) ? "[generalized version]" : "");
  if (gen->set & SROU_SET_CDFMODE)
    _unur_string_append(info,"   use CDF at mode\n");
  if (gen->variant & SROU_VARFLAG_SQUEEZE)
    _unur_string_append(info,"   use squeeze\n");
  if (gen->variant & SROU_VARFLAG_MIRROR)
    _unur_string_append(info,"   use mirror principle\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  if (gen->set & SROU_SET_R) {
    rc = unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize);
    _unur_string_append(info,"   enveloping rectangle = (%g,%g) x (%g,%g)\n",
			GEN->vl,GEN->vr, 0.,GEN->um);
    _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n", rc);
  }
  else {
    _unur_string_append(info,"   bounding rectangle = (%g,%g) x (%g,%g)\n",
			GEN->vl,GEN->vr, 0.,GEN->um);
    h_area = (GEN->vr - GEN->vl) * GEN->um;
    _unur_string_append(info,"   area(hat) = %g\n", h_area);
    if (gen->set & SROU_SET_CDFMODE) 
      rc = 2.;
    else 
      rc = (gen->variant & SROU_VARFLAG_MIRROR) ? 2.829 : 4.;

    _unur_string_append(info,"   rejection constant = %g\n", rc);
  }
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"     r = %g  %s\n", GEN->r,
			(gen->set & SROU_SET_R) ? "" : "[default]");
    if (gen->set & SROU_SET_CDFMODE)
      _unur_string_append(info,"   cdfatmode = %g\n", GEN->Fmode); 
    else
      _unur_string_append(info,"   cdfatmode = [not set]\n"); 

    if (gen->variant & SROU_VARFLAG_SQUEEZE)
      _unur_string_append(info,"   usesqueeze\n");

    if (gen->variant & SROU_VARFLAG_MIRROR)
      _unur_string_append(info,"   usemirror\n");

    if (gen->variant & SROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
    */
  }


  /* Hints */
  if (help) {
    if ( !(gen->set & SROU_SET_CDFMODE))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"cdfatmode\" to reduce the rejection constant.");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_srou_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
