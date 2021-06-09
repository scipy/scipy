/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      auto.c                                                       *
 *                                                                           *
 *   selects a method for a given distribution object AUTOmatically          *
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
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "x_gen.h"
#include "auto.h"
#include "auto_struct.h"

#include "cstd.h"
#include "dari.h"
#include "dgt.h"
#include "dstd.h"
#include "empk.h"
#include "hist.h"
#include "mvstd.h"
#include "tdr.h"
#include "vempk.h"

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define AUTO_SET_LOGSS          0x001u

/*---------------------------------------------------------------------------*/

#define GENTYPE "AUTO"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_auto_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_init_cont( struct unur_par *par );
static struct unur_gen *_unur_init_cvec( struct unur_par *par );
static struct unur_gen *_unur_init_discr( struct unur_par *par );
static struct unur_gen *_unur_init_cemp( struct unur_par *par );
static struct unur_gen *_unur_init_cvemp( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize generator object for detected distribution type.               */
/*---------------------------------------------------------------------------*/

/* static struct unur_gen *_unur_auto_create( struct unur_par *par );        */
/* static double _unur_auto_sample( struct unur_gen *gen );                  */
/* static void _unur_auto_free( struct unur_gen *gen);                       */
/*---------------------------------------------------------------------------*/
/* no such functions!                                                        */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_auto_par*)par->datap) /* data for parameter object */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_auto_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... dummy pointer to distribution object                     */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   for testing.                                                       */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL(GENTYPE,distr,NULL);

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_auto_par) );
  COOKIE_SET(par,CK_AUTO_PAR);

  /* copy input */
  par->distr    = distr;           /* pointer to distribution object         */

  /* set default values */
  par->method   = UNUR_METH_AUTO;  /* method and default variant             */
  par->variant  = 0u;              /* default variant                        */
  par->set      = 0u;              /* inidicate default parameters           */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_auto_init;

  return par;

} /* end of unur_auto_new() */

/*---------------------------------------------------------------------------*/

int 
unur_auto_set_logss( UNUR_PAR *par, int logss )
     /*----------------------------------------------------------------------*/
     /* set common logarithm of sample size                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   logss ... common logarithm of sample size                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AUTO );

  if (logss < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"log < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->logss = logss;

  /* changelog */
  par->set |= AUTO_SET_LOGSS;

  return UNUR_SUCCESS;

} /* end of unur_auto_set_logss() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_auto_init( struct unur_par *par )
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
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_AUTO ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_AUTO_PAR,NULL);

  /* get distribution type */
  switch (par->distr->type) {
  case UNUR_DISTR_CONT:
    gen = _unur_init_cont( par );
    break;
  case UNUR_DISTR_CVEC:
    gen = _unur_init_cvec( par );
    break;
  case UNUR_DISTR_DISCR:
    gen = _unur_init_discr( par );
    break;
  case UNUR_DISTR_CEMP:
    gen = _unur_init_cemp( par );
    break;
  case UNUR_DISTR_CVEMP:
    gen = _unur_init_cvemp( par );
    break;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    gen = NULL;
    break;
  }

  /* copy URNGs and debugging flags */
  if (gen) {
    gen->urng = par->urng;
    gen->urng_aux = par->urng_aux;
    gen->debug = par->debug;
  }

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_auto_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Initialize generator object for detected distribution type             **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_init_cont( struct unur_par *par_auto )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for distribution type CONT                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_auto ... pointer to paramters for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;
  struct unur_gen *gen;

  do {
    /* 1st choice: TDR */
    par = unur_tdr_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;

    /* 2nd choice: CSTD */
    par = unur_cstd_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;

  } while (0);

  return gen;
} /* end of _unur_init_cont() */
  
/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_init_cvec( struct unur_par *par_auto )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for distribution type CVEC                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_auto ... pointer to paramters for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;
  struct unur_gen *gen;
  
  /* Choose a method: MVSTD */
  par = unur_mvstd_new(par_auto->distr);
     
  /* Create generator object */
  gen = unur_init(par);

  return gen;
} /* end of _unur_init_cvec() */
  
/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_init_discr( struct unur_par *par_auto )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for distribution type DISCR                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_auto ... pointer to paramters for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;
  struct unur_gen *gen;
  
  do {
    /* 1st choice: DGT */
    if (par_auto->distr->data.discr.pv != NULL) {
      par = unur_dgt_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
    }

    /* 2nd choice: DARI */
    if (par_auto->distr->data.discr.pmf != NULL) {
      par = unur_dari_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
      /* try again with DGT (in this case we have to compute the PV) */
      par = unur_dgt_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
    }

    /* 3rd choice: DSTD */
    par = unur_dstd_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
    
  } while (0);

  return gen;
} /* end of _unur_init_discr() */
  
/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_init_cemp( struct unur_par *par_auto )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for distribution type CEMP                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_auto ... pointer to paramters for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;
  struct unur_gen *gen;
  
  do {
    /* 1st choice: EMPK  [requires raw data] */
    par = unur_empk_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;

    /* 2nd choice: HIST  [requires histogram] */
    par = unur_hist_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;

  } while(0);
     
  return gen;
} /* end of _unur_init_cemp() */
  
/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_init_cvemp( struct unur_par *par_auto )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for distribution type CVEMP                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_auto ... pointer to paramters for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;
  struct unur_gen *gen;

  /* Choose a method: VEMPK */
  par = unur_vempk_new(par_auto->distr);
     
  /* Create generator object */
  gen = unur_init(par);
     
  return gen;
} /* end of _unur_init_cvemp() */
  
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

/** None **/

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
