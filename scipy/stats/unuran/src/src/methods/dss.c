/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dss.c                                                        *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    sequential search                                            *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PV or PMF together with sum over PV or PMF,                    *
 *      or the VDF,                                                          *
 *                                                                           *
 *   REQUIRED:  pointer to probability vector, sum over PV                   *
 *              pointer to PMF, sum over PMF                                 *
 *              pointer to CDF                                               *
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
 *   SUMMARY:                                                                *
 *   [1] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This method a varient of the inversion method, i.e.                     *
 *                                                                           *
 *   (1) Generate a random number U ~ U(0,1).                                *
 *   (2) Find largest integer I such that F(I) = P(X<=I) <= U.               *
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
#include "dss.h"
#include "dss_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DSS_VARIANT_NONE       0x000u     /* invalid data                    */
#define DSS_VARIANT_PV         0x001u     /* use PV                          */
#define DSS_VARIANT_PMF        0x002u     /* use PMF                         */
#define DSS_VARIANT_CDF        0x004u     /* use CDF                         */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DSS_DEBUG_REINIT        0x00000010u   /* print parameters after reinit  */
#define DSS_DEBUG_PRINTVECTOR   0x00000100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "DSS"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dss_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dss_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dss_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dss_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dss_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dss_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dss_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dss_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dss_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr      /* data for distribution object      */

#define PAR       ((struct unur_dss_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dss_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

#define PMF(x)    _unur_discr_PMF((x),(gen->distr))   /* call to PMF         */
#define CDF(x)    _unur_discr_CDF((x),(gen->distr))   /* call to CDF         */

/*---------------------------------------------------------------------------*/

#define _unur_dss_getSAMPLE(gen)   (_unur_dss_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dss_new( const struct unur_distr *distr )
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
  unsigned variant = DSS_VARIANT_NONE;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.pv && (distr->set & UNUR_DISTR_SET_PMFSUM))
    variant = DSS_VARIANT_PV;
  else if (DISTR_IN.pmf && (distr->set & UNUR_DISTR_SET_PMFSUM)) 
    variant = DSS_VARIANT_PMF;
  else if (DISTR_IN.cdf)
    variant = DSS_VARIANT_CDF;

  if (variant == DSS_VARIANT_NONE) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV+sum, PMF+sum, or CDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dss_par) );
  COOKIE_SET(par,CK_DSS_PAR);

  /* copy input */
  par->distr       = distr;          /* pointer to distribution object       */

  /* set default values */
  par->method      = UNUR_METH_DSS;  /* method                               */
  par->variant     = variant;        /* variant of method                    */
  par->set         = 0u;             /* inidicate default parameters         */    
  par->urng        = unur_get_default_urng(); /* use default urng            */
  par->urng_aux    = NULL;                    /* no auxilliary URNG required */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dss_init;

  return par;

} /* end of unur_dss_new() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_dss_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;         /* pointer to generator object */
  
  /* check arguments */
  CHECK_NULL(par,NULL);
  
  /* check input */
  if ( par->method != UNUR_METH_DSS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DSS_PAR,NULL);
  
  /* create a new empty generator object */
  gen = _unur_dss_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_dss_debug_init(gen);
#endif

  return gen;
} /* end of _unur_dss_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dss_reinit( struct unur_gen *gen )
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
  if ( (rcode = _unur_dss_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* (re)set sampling routine */
  SAMPLE = _unur_dss_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
  if (gen->debug & DSS_DEBUG_REINIT) _unur_dss_debug_init(gen);
#endif

  return UNUR_SUCCESS;
} /* end of _unur_dss_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dss_create( struct unur_par *par )
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
  struct unur_gen *gen;       /* pointer to generator object */

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSS_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dss_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DSS_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dss_getSAMPLE(gen);
  gen->destroy = _unur_dss_free;
  gen->clone = _unur_dss_clone;
  gen->reinit = _unur_dss_reinit;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dss_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dss_create() */

/*---------------------------------------------------------------------------*/

int
_unur_dss_check_par( struct unur_gen *gen )
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
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    if (DISTR.pv != NULL) break;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV");
    return UNUR_ERR_DISTR_REQUIRED;

  case DSS_VARIANT_PMF:
    if (DISTR.pmf != NULL) break;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF");
    return UNUR_ERR_DISTR_REQUIRED;

  case DSS_VARIANT_CDF:
    if (DISTR.cdf != NULL) return UNUR_SUCCESS;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF");
    return UNUR_ERR_DISTR_REQUIRED;

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* additionally check for required data: sum over PMF */
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return UNUR_ERR_DISTR_REQUIRED;
    }

  return UNUR_SUCCESS;
} /* end of _unur_dss_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dss_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dss_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSS_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_dss_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_dss_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if (!gen) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_DSS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DSS_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_dss_free() */

/*****************************************************************************/

int
_unur_dss_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (sample from random variate)                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{ 
  int J;
  double U;
  double sum;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSS_GEN,INT_MAX);

  switch(gen->variant) {
  case DSS_VARIANT_PV:
    U = DISTR.sum * _unur_call_urng(gen->urng);
    sum = 0.;
    for (J=0; J<DISTR.n_pv; J++) {
      sum += DISTR.pv[J];
      if (sum >= U) break;
    }
    return (J + DISTR.domain[0]);

  case DSS_VARIANT_PMF:
    U = DISTR.sum * _unur_call_urng(gen->urng);
    sum = 0.;
    for (J=DISTR.domain[0]; J<=DISTR.domain[1]; J++) {
      sum += PMF(J);
      if (sum >= U) break;
    }
    return J;

  case DSS_VARIANT_CDF:
    U = _unur_call_urng(gen->urng);
    for (J=DISTR.domain[0]; J<=DISTR.domain[1]; J++) {
      if (CDF(J) >= U) break;
    }
    return J;

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INT_MAX;
  }

} /* end of _unur_dss_sample() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_dss_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = sequential search\n",gen->genid);

  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DSS_DEBUG_PRINTVECTOR));

  fprintf(LOG,"%s: sampling routine = _unur_dss_sample()\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    fprintf(LOG,"use PV\n");  break;
  case DSS_VARIANT_PMF:
    fprintf(LOG,"use PMF\n"); break;
  case DSS_VARIANT_CDF:
    fprintf(LOG,"use CDF\n"); break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  }
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_dss_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dss_info( struct unur_gen *gen, int help )
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

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    _unur_string_append(info,"   functions = PV  [length=%d]\n",DISTR.domain[1]-DISTR.domain[0]+1);
    break;
  case DSS_VARIANT_PMF:
    _unur_string_append(info,"   functions = PMF\n");
    break;
  case DSS_VARIANT_CDF:
    _unur_string_append(info,"   functions = CDF\n");
    break;
  }
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: DSS (Simple Sequential Search)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics: slow\n");
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_dss_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
