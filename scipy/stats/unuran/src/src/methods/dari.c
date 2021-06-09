/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dari.c                                                       *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    discrete automatic rejection inversion                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PMF of a T(-1/2)-concave distribution;                         *
 *      produce a value x consistent with its PMF.                           *
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
 *   [1] Hoermann W. and G.Derflinger (1997):                                *
 *       An automatic generator for a large class of discrete unimodal       *
 *       distributions, in A.R. Kaylan and A. Lehmann, ESM 97, pp 139-144    *
 *                                                                           *
 *   [2] Hoermann W. and G.Derflinger (1996):                                *
 *       Rejection-inversion to generate variates from monotone discrete     *
 *       distributions, ACM TOMACS 6(3), 169-184                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dari.h"
#include "dari_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DARI_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DARI_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DARI_SET_CFACTOR        0x001u
#define DARI_SET_TABLESIZE      0x002u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DARI"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dari_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dari_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dari_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dari_sample( struct unur_gen *gen );
static int _unur_dari_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_dari_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute hat.                                                              */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dari_debug_init( struct unur_gen *gen, const char *status );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dari_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       ((struct unur_dari_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dari_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */     

#define PMF(x)    _unur_discr_PMF((x),(gen->distr))    /* call to PMF        */

/*---------------------------------------------------------------------------*/

#define T(x) (-1./sqrt(x))    /* transformation for PMF                      */
#define F(x) (-1./(x))        /* anti-derivative of inverse transformation   */
#define FM(x) (-1./(x))       /* inverse of anti-derivative                  */
#define N0 (GEN->n[0])        /* first position in table                     */

/*---------------------------------------------------------------------------*/

#define _unur_dari_getSAMPLE(gen) \
   ( ((gen)->variant & DARI_VARFLAG_VERIFY) \
     ? _unur_dari_sample_check : _unur_dari_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dari_new( const struct unur_distr *distr )
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
     /*                                                                      */
     /* comment:                                                             */
     /*   if the sum over the PMF is not close to 1 it is necessary to       */
     /*   set pmf_sum to an approximate value of its sum (+/- 30 % is ok).  */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;
  
  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );
  
  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  
  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }

  if (DISTR_IN.domain[0] < 0) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"domain contains negative numbers");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dari_par) );
  COOKIE_SET(par,CK_DARI_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR->c_factor  = 0.664;
  /* optimal value for the normal distribution, which is good for all        */
  /* bell-shaped densities. The minimax approach for that transformation     */
  /* has c_factor = 2.                                                       */

  PAR->squeeze   = 0; /* no squeezes by default as squeezes slow down the     */
  /* sampling for most distributions if PAR->size is big enough. Squeeze is   */
  /* important for the speed only when small samples are required or when    */
  /* the domain of the distribution is very big. (much bigger than PAR->size) */

  PAR->size      = 100; /*size of table that stores the "rejection point" for */
  /* all integers close to the mode when needed the first time while         */
  /* sampling; can speed up the generation considerably.                     */

  par->method   = UNUR_METH_DARI;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dari_init;

  return par;

} /* end of unur_dari_new() */

/*****************************************************************************/

int
unur_dari_set_cpfactor( struct unur_par *par, double cpfactor )
     /*----------------------------------------------------------------------*/
     /* set factor for position of left and right construction point         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   cfactor ... factor                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, DARI );

  /* check new parameter for generator */
  /** TODO: welche werte fuer c sind zulaessig / sinnvoll ? 
  zulaessig ist jedes c>0, man koennte damit falsche Flaechenangaben kompensieren.
  Wenn sum genau bekannt ist, ist ein c > 2 (2 ist der minimax approach) so weit
  ich weiss nie sinnvoll. Ich denke aber, das sollte man besser nicht prinzipiell
  verbieten, hoechstens eine warnung.**/
  if (cpfactor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor <= 0");
    return UNUR_ERR_PAR_SET;
  }

  if (cpfactor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor > 2 not recommended. skip");

  /* store date */
  PAR->c_factor = cpfactor;

  /* changelog */
  par->set |= DARI_SET_CFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_dari_set_cpfactor() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_squeeze( struct unur_par *par, int squeeze )
     /*----------------------------------------------------------------------*/
     /* turn on/off using squeezes                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   squeeze ... 0 = no squeeze,  !0 = use squeeze                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  
  /* check input */
  _unur_check_par_object( par, DARI );

  /* store data */
  PAR->squeeze = squeeze;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_dari_set_squeeze() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_tablesize( struct unur_par *par, int size )
     /*----------------------------------------------------------------------*/
     /* set size of table                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   size   ... table size                                              */
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

  /* check input */
  _unur_check_par_object( par, DARI );

  /* check parameter */
  if (size < 0) {  
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid table size");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store data */
  PAR->size = size;

  /* changelog */
  par->set |= DARI_SET_TABLESIZE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dari_set_tablesize() */
  
/*---------------------------------------------------------------------------*/

int
unur_dari_set_verify( struct unur_par *par, int verify )
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

  /* check input */
  _unur_check_par_object( par, DARI );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | DARI_VARFLAG_VERIFY) : (par->variant & (~DARI_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_dari_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_discr_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= DARI_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~DARI_VARFLAG_VERIFY;

  SAMPLE = _unur_dari_getSAMPLE(gen); 

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_dari_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_dari_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, NULL );

  /* check input */
  if ( par->method != UNUR_METH_DARI ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DARI_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dari_create(par);

  /* free parameters */
  _unur_par_free(par);

  if (!gen) return NULL;

  /* check parameters */
  if (_unur_dari_check_par(gen) != UNUR_SUCCESS) {
    _unur_dari_free(gen); return NULL;
  }

  /* create hat and squeeze (setup procedure) */
  if ( _unur_dari_hat(gen)!=UNUR_SUCCESS ) {
    /* error */
    _unur_dari_free(gen); return NULL;
  }

  /* hat successfully created */
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_dari_debug_init(gen,"INIT completed");
#endif

  /* o.k. */
  return gen;

} /* end of _unur_dari_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dari_reinit( struct unur_gen *gen )
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
  int result;

  /* check parameters */
  if ( (result = _unur_dari_check_par(gen)) != UNUR_SUCCESS)
    return result;

  /* compute hat  */
  if ( (result = _unur_dari_hat( gen )) != UNUR_SUCCESS)
    return result;

  /* (re)set sampling routine */
  SAMPLE = _unur_dari_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & DARI_DEBUG_REINIT)
    _unur_dari_debug_init(gen,"REINIT completed");
#endif

  /* hat successfully created */
  return UNUR_SUCCESS;
} /* end of _unur_dari_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dari_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DARI_PAR,NULL);
  
  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dari_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DARI_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dari_getSAMPLE(gen);
  gen->destroy = _unur_dari_free;
  gen->clone = _unur_dari_clone;
  gen->reinit = _unur_dari_reinit;

  /* copy some parameters into generator object */
  GEN->squeeze = PAR->squeeze;        /* squeeze yes/no?                       */
  GEN->c_factor = PAR->c_factor;      /* constant for choice of design point   */

  /* size of auxiliary table; 0 for none
     it cannot be larger than the given domain (avoid overflow) */
  if ((unsigned)DISTR.BD_RIGHT - (unsigned)DISTR.BD_LEFT < INT_MAX)
    GEN->size = _unur_min(PAR->size,DISTR.BD_RIGHT-DISTR.BD_LEFT+1);
  else /* length of interval > INT_MAX */
    GEN->size = PAR->size;

  /* allocate */
  GEN->hp = (GEN->size > 0) ? _unur_xmalloc( GEN->size * sizeof(double) ) : NULL;
  GEN->hb = (GEN->size > 0) ? _unur_xmalloc( GEN->size * sizeof(char) )   : NULL;

  /* initialize parameters */
  GEN->vt=0.;            /* total volume below hat                           */
  GEN->vc=0.;            /* volume below center part                         */
  GEN->vcr=0.;           /* volume center and right together                 */

  GEN->xsq[0]=0.;        /* value necessary for the squeeze computation      */
  GEN->xsq[1]=0.;        /* value necessary for the squeeze computation      */
  GEN->y[0]=0.;          /* value of the transformed density in points of contact */
  GEN->y[1]=0.;          /* value of the transformed density in points of contact */
  GEN->ys[0]=0.;         /* the slope of the transformed hat                 */
  GEN->ys[1]=0.;         /* the slope of the transformed hat                 */
  GEN->ac[0]=0.;         /* left and right starting point of the uniform hat
                            in the center                                    */
  GEN->ac[1]=0.;         /* left and right starting point of the uniform hat
                            in the center                                    */

  GEN->pm=0.;            /* mode probability                                 */
  GEN->Hat[0]=0.;        /* point where the hat starts for the left and
                            the right tail                                   */
  GEN->Hat[1]=0.;        /* point where the hat starts for the left and
                            the right tail                                   */

  GEN->m=0;              /* mode                                             */
  GEN->x[0]=0;           /* points of contact left and right of the mode     */
  GEN->x[1]=0;           /* points of contact left and right of the mode     */
  GEN->s[0]=0;           /* first and last integer of the center part        */
  GEN->s[1]=0;           /* first and last integer of the center part        */
  GEN->n[0]=0;           /* contains the first and the last i 
                            for which values are stored in table             */
  GEN->n[1]=0;           /* contains the first and the last i 
                            for which values are stored in table             */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dari_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_dari_create() */

/*---------------------------------------------------------------------------*/

int
_unur_dari_check_par( struct unur_gen *gen )
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
    if (unur_distr_discr_upd_mode( gen->distr )!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* check mode. we assume unimodality 
     since otherwise the PMF would not be T-concave! */
  if (DISTR.BD_LEFT > DISTR.mode)
    DISTR.mode = DISTR.BD_LEFT;
  else if (DISTR.BD_RIGHT < DISTR.mode)
    DISTR.mode = DISTR.BD_RIGHT;

  /* check for required data: sum over PMF */
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS)
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF; use default");

  /* sum must not be zero */
  if (DISTR.sum <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"sum <= 0");
    return UNUR_ERR_GEN_DATA;
  }
 
  return UNUR_SUCCESS;
} /* end of _unur_dari_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dari_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dari_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DARI_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy additional data */
  if (GEN->size > 0) {
    CLONE->hp = _unur_xmalloc( GEN->size * sizeof(double) );
    memcpy( CLONE->hp, GEN->hp, GEN->size * sizeof(double) );
    CLONE->hb = _unur_xmalloc( GEN->size * sizeof(char) );
    memcpy( CLONE->hb, GEN->hb, GEN->size * sizeof(char) );
  }
  
  return clone;

#undef CLONE
} /* end of _unur_dari_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_dari_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DARI ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);
  
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */
  
  /* free two auxiliary tables */
  if (GEN->hp)   free(GEN->hp);
  if (GEN->hb)   free(GEN->hb);
  
  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_dari_free() */

/*****************************************************************************/

int
_unur_dari_sample( struct unur_gen *gen )
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
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{
  double U, h;
  double X = 0.;
  int k,i;
  static const int sign[2] = {-1,1};

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DARI_GEN,INT_MAX);

  /* step 1.0 */
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->vt;

    /* step 1.1 */
    if (U<=GEN->vc) {
      X = U * (GEN->ac[1]-GEN->ac[0]) / GEN->vc + GEN->ac[0]; 
      k = (int)(X+0.5);
      i = (k<GEN->m) ? 0 : 1;
      if (GEN->squeeze && sign[i]*(GEN->ac[i]-GEN->s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = 0.5 - PMF(k)/GEN->pm;
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = 0.5-PMF(k)/GEN->pm;
      }
      if (h <= sign[i]*(k-X))
	return k;
    }

    /*step 1.2*/ 
    else {
      if (U<= GEN->vcr) {
	i = 1;
	U -= GEN->vc;
      } 
      else {
	i = 0;
	U -= GEN->vcr;
      }

      U = GEN->Hat[i] + sign[i]*U; 
      X = GEN->x[i] + (FM(U*GEN->ys[i])-GEN->y[i]) / GEN->ys[i];
      k = (int)(X+0.5);

      if (GEN->squeeze && (sign[i]*k <= sign[i]*GEN->x[i]+1) && (GEN->xsq[i] <= sign[i]*(X-k))) 
	return k;

      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k);
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i]-PMF(k);
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
} /* end of _unur_dari_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_dari_sample_check( struct unur_gen *gen )
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
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{
  double U, h;
  double X = 0.;
  double hkm05;
  int k,i;
  static const int sign[2] = {-1,1};
  
  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DARI_GEN,INT_MAX);
  
  /* step 1.0 */
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->vt;

    /* step 1.1 */
    if (U <= GEN->vc) {
      X = U * (GEN->ac[1]-GEN->ac[0]) / GEN->vc + GEN->ac[0]; 
      k = (int)(X+0.5);
      i = (k<GEN->m) ? 0 : 1;
      if (GEN->squeeze && sign[i]*(GEN->ac[i]-GEN->s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = 0.5 - PMF(k)/GEN->pm;
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
	/* CHECKING HAT */
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_log_printf(gen->genid,__FILE__,__LINE__,
			   "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN->pm ); 
        }
	/* end CHECKING HAT */
      }
      else {
	h = 0.5 - PMF(k)/GEN->pm;
	/* CHECKING HAT */
	/* here UNUR_EPSILON can be too small for distributions that 
	   have two neighbouring hats. */
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_log_printf(gen->genid,__FILE__,__LINE__,
			   "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN->pm ); 
        }
	/* end CHECKING HAT */
      }
      if (h <= sign[i]*(k-X))
	return k;
    }

    /*step 1.2*/ 
    else {
      if (U<= GEN->vcr) {
	i = 1;
	U -= GEN->vc;
      } 
      else {
	i = 0;
	U -= GEN->vcr;
      }

      U = GEN->Hat[i] + sign[i]*U; 
      X = GEN->x[i] + (FM(U*GEN->ys[i])-GEN->y[i]) / GEN->ys[i];
      k = (int)(X+0.5);
      /* this is for a very rare case that for k of the tail closest to 
	 the mode an x value farer away than 0.5 is generated. */
      if(k==GEN->s[i]) 
	k += sign[i];

      if (GEN->squeeze && (sign[i]*k <= sign[i]*GEN->x[i]+1) && (GEN->xsq[i] <= sign[i]*(X-k))) 
	return k;

      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if(!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k); 

	  /* CHECKING HAT: (only necessary if(k!=GEN->s+1) as for the border
                            the hat is by construction correct)
	     tests if Hat too low i.e.: (H(k+0.5)-  p_k < H(k-0.5)) */
          if(k != GEN->s[i]+sign[i]) {
            hkm05 = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k-sign[i]*0.5-GEN->x[i])) / GEN->ys[i];
	    if (GEN->hp[k-N0]+UNUR_EPSILON < hkm05) {
	      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			  "for tailpart hat too low, ie hp[k] < H(k-0.5)");
	      _unur_log_printf(gen->genid,__FILE__,__LINE__,
			       "k %d hp  %.20e H(k-0.5) %.20e ", k,GEN->hp[k-N0],hkm05 ); 
            }
	  }
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k);
	/* CHECKING HAT:(only necessary if(k!=GEN->s+1) as for the border
                            the hat is by construction correct)
	   tests if Hat too low i.e.: (H(k+0.5)-p_k < H(k-1/2)) */
        hkm05 = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k-sign[i]*0.5-GEN->x[i])) / GEN->ys[i];
        if(k != GEN->s[i]+sign[i]) {
	  if (h+UNUR_EPSILON < hkm05) {
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			"PMF(i) > hat(i) for tailpart");
	    _unur_log_printf(gen->genid,__FILE__,__LINE__,
			     "k %d h  %.20e H(k-0.5) %.20e ", k,h,hkm05 ); 
	  }
        }
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
  
} /* end of _unur_dari_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_dari_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute hat                                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int sign[2] = {-1,1};
  int b[2], d, i, j;
  double v[2], at[2];
  double t0 = 1.;
  int setup = 1;
  int rep = 1;
  
  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_DARI_GEN, UNUR_ERR_COOKIE );
  
  /* Step 0: setup */
  GEN->m = DISTR.mode;
  b[0] = DISTR.BD_LEFT;
  b[1] = DISTR.BD_RIGHT;
  GEN->pm = PMF(GEN->m);
  d = _unur_max(2, (int)( GEN->c_factor/(GEN->pm/DISTR.sum)));

  /* check mode */
  if (_unur_iszero(GEN->pm)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PMF(mode)=0");
    return UNUR_ERR_GEN_DATA;
  }

  /* step 0.1 */
  do {
    for(i=0; i<=1; i++) {
      GEN->x[i] = GEN->m + sign[i] * d;
      if (sign[i]*GEN->x[i]+1 > sign[i]*b[i]) {
	v[i] = 0; 
	GEN->s[i] = b[i];
      }
      else {
	GEN->y[i] = T( PMF(GEN->x[i]) );
	GEN->ys[i] = sign[i] * (T( PMF(GEN->x[i]+sign[i])) - GEN->y[i]);
	if (GEN->ys[i]*sign[i] > -DBL_EPSILON) {
	  setup = -setup; /* indicate that the hat is not ok */
	  i = 1; 
	}
        else {
	  GEN->s[i] = (int)(0.5+GEN->x[i]+(T(GEN->pm)-GEN->y[i])/GEN->ys[i]);
	  GEN->Hat[i] = ( F(GEN->y[i]+GEN->ys[i]*(GEN->s[i]+sign[i]*1.5-GEN->x[i])) /
			  GEN->ys[i]-sign[i]*PMF(GEN->s[i]+sign[i]) ); 
	  at[i] = GEN->x[i] + (FM(GEN->ys[i]*GEN->Hat[i])-GEN->y[i]) / GEN->ys[i]; 
          if(GEN->squeeze)
	    GEN->xsq[i] = sign[i]*(at[i]-(GEN->s[i]+sign[i]));
	  v[i] = sign[i]*(F(GEN->y[i]+GEN->ys[i]*(b[i]+sign[i]*0.5-GEN->x[i]))/
			  GEN->ys[i]-F(GEN->y[i]+GEN->ys[i]*(at[i]-GEN->x[i]))/GEN->ys[i]);
	}
      }
      if (setup>0)
	GEN->ac[i] = GEN->s[i] + sign[i]*(PMF(GEN->s[i])/GEN->pm-0.5);
    }

    /* step 0.2 */
    if(setup>0) {
      GEN->vc = GEN->pm*(GEN->ac[1]-GEN->ac[0]); 
      GEN->vt = GEN->vc+v[0]+v[1];
      GEN->vcr = GEN->vc+v[1];

      /* step 0.3 */ 
      GEN->n[0] = _unur_max(b[0],GEN->m - GEN->size/2);
      GEN->n[1] = GEN->n[0] + GEN->size - 1;
      if (GEN->n[1] > b[1]) {
	GEN->n[1] = b[1];
	GEN->n[0] = GEN->n[1]- GEN->size + 1;
      }
      /* initialize table */
      for (j=0; j<GEN->size; j++)
	GEN->hb[j] = 0;
    }

    /* setup == 1 first try, up to now ok,  ==2 second try, up to now ok */
    /* setup == -1 first try, not ok,  == -2 second try, not ok          */

    if (setup == 1 || setup == -1) {
      t0= 2. * DISTR.sum;
      if (setup==1 && GEN->vt<=t0)
	rep=0;
      else { 
	setup = 2;
	d = (int) (t0 / GEN->pm);
      }
    }
    else 
      rep=0; 
  } while(rep);

  if (setup == -2 || GEN->vt > 100.*t0 || !(GEN->vt > 0.)) {
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug)
      _unur_dari_debug_init(gen,"RE/INIT failed try again");
#endif
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"Area below hat too large or zero!! possible reasons: PDF, mode or area below PMF wrong;  or PMF not T-concave");
    return UNUR_ERR_GEN_DATA;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_dari_hat() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_dari_debug_init( struct unur_gen *gen, const char *status )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   status ... status of re/init routine                               */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = dari (discrete automatic rejection inversion)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );

  fprintf(LOG,"%s: sampling routine = _unur_dari_sample",gen->genid);
  if (gen->variant & DARI_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s:area below hat: total %f center: %f, left tail %f, right tail %f\n", gen->genid,
	  GEN->vt, GEN->vc, GEN->vt-GEN->vcr, GEN->vcr-GEN->vc);
  fprintf(LOG,"%s: mode %d and mode probability %f\n",gen->genid, GEN->m, GEN->pm); 
  for(i=0;i<=1;i++) {
    fprintf(LOG,"%s:i=%d: x=%d; Hat=%f; ac=%f; s=%d;\n", gen->genid,
	    i, GEN->x[i], GEN->Hat[i], GEN->ac[i], GEN->s[i]);
    fprintf(LOG,"%s:i=%d: xsq=%f; y=%f; ys=%f; n:=%d (for aux.table)\n", gen->genid,
	    i, GEN->xsq[i], GEN->y[i], GEN->ys[i], GEN->n[i]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: %s ************\n",gen->genid, status );
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_dari_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dari_info( struct unur_gen *gen, int help )
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

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PMF\n");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %d   %s\n", DISTR.mode,
                      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   sum(PMF)  = %g   %s\n", DISTR.sum,
                      (distr->set & UNUR_DISTR_SET_PMFSUM) ? "" : "[unknown]");
  _unur_string_append(info,"\n");

  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"mode\".");
    if (!(distr->set & UNUR_DISTR_SET_PMFSUM))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"pmfsum\".");
    _unur_string_append(info,"\n");
  }

  /* method */
  _unur_string_append(info,"method: DARI (Discrete Automatic Rejection Inversion)\n");
  if (GEN->size == 0) 
    _unur_string_append(info,"   no table\n");
  else
    _unur_string_append(info,"   use table of size %d\n", GEN->size);
  if (GEN->squeeze)
    _unur_string_append(info,"   use squeeze\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   sum(hat) = %g\n",GEN->vt);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PMFSUM)
    _unur_string_append(info,"= %g\n", GEN->vt/DISTR.sum);
  else
    _unur_string_append(info,"= %.2f  [approx.]\n", 
			unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   tablesize = %d  %s\n", GEN->size,
			(gen->set & DARI_SET_TABLESIZE) ? "" : "[default]");
    if (GEN->squeeze)
      _unur_string_append(info,"   squeeze = on\n");

    if (gen->set & DARI_SET_CFACTOR)
      _unur_string_append(info,"   cpfactor = %g\n",   GEN->c_factor);

    if (gen->variant & DARI_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
} /* end of _unur_dari_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
