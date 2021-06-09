/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_methods.c                                              *
 *                                                                           *
 *   Deprecated routines                                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   THESE ROUTINES SHOULD NOT BE USED ANY MORE!                             *
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
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <distr/discr.h>
#include "unur_methods_source.h"
#include "utdr_struct.h"
#include "x_gen.h"
#include "deprecated_methods.h"

/*---------------------------------------------------------------------------*/
#ifdef USE_DEPRECATED_CODE
/*---------------------------------------------------------------------------*/

#define GENTYPE "deprecated"         /* type of generator                    */

#define BD_LEFT   domain[0]      /*  left boundary of domain of distribution */
#define BD_RIGHT  domain[1]      /* right boundary of domain of distribution */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

int 
unur_cstd_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, CSTD, UNUR_ERR_GEN_INVALID );
  if (n_params>0) CHECK_NULL(params, UNUR_ERR_NULL);

  /* set new parameters in distribution object */
  if (unur_distr_cont_set_pdfparams(gen->distr, params,n_params)!=UNUR_SUCCESS)
    return UNUR_ERR_DISTR_SET;

  /* reinit */
  return gen->reinit(gen);
} /* end of unur_cstd_chg_pdfparams() */

/*****************************************************************************/

#define DISTR gen->distr->data.discr /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

int
unur_dari_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_dari_reinit() */
  
/*---------------------------------------------------------------------------*/

int
unur_dari_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  /* set new parameters in distribution object */
  return unur_distr_discr_set_pmfparams(gen->distr,params,n_params);

} /* end of unur_dari_chg_pmfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_dari_chg_domain( struct unur_gen *gen, int left, int right )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity.           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_dari_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_mode( struct unur_gen *gen, int mode )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );
  
  /* copy parameters */
  DISTR.mode = mode;
  
  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dari_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  return unur_distr_discr_upd_mode( gen->distr );
} /* end of unur_dari_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_pmfsum( struct unur_gen *gen, double sum )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change sum over PMF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   sum   ... sum                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );
  
  /* check new parameter for generator */
  if (sum <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"sum <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  DISTR.sum = sum;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dari_chg_pmfsum() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_pmfsum( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute sum over PMF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );

  return unur_distr_discr_upd_pmfsum( gen->distr );
} /* end of unur_dari_upd_pmfsum() */

/*---------------------------------------------------------------------------*/

#undef DISTR

/*****************************************************************************/

#define DISTR gen->distr->data.discr /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_DATA );

  return gen->reinit(gen);
} /* end of unur_dsrou_reinit() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
  
  /* set new parameters in distribution object */
  return unur_distr_discr_set_pmfparams(gen->distr,params,n_params);

} /* end of unur_dsrou_chg_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_mode( struct unur_gen *gen, int mode )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dsrou_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );

  return unur_distr_discr_upd_mode( gen->distr );
} /* end of unur_dsrou_upd_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_dsrou_chg_domain( struct unur_gen *gen, int left, int right )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  
  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;
  
  /* changelog */
  gen->distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN | UNUR_DISTR_SET_MASK_DERIVED );
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_dsrou_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_pmfsum( struct unur_gen *gen, double sum )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change sum over PMF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   sum   ... sum                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
  
  /* check new parameter for generator */
  if (sum <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"sum <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  DISTR.sum = sum;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dsrou_chg_pmfsum() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_upd_pmfsum( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute sum over PMF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );

  return unur_distr_discr_upd_pmfsum( gen->distr );
} /* end of unur_dsrou_upd_pmfsum() */

/*---------------------------------------------------------------------------*/

#undef DISTR

/*****************************************************************************/

int 
unur_dstd_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSTD, UNUR_ERR_GEN_INVALID );
  if (n_params>0) CHECK_NULL(params, UNUR_ERR_NULL);
  
  /* set new parameters in distribution object */
  if (unur_distr_discr_set_pmfparams(gen->distr,params,n_params)!=UNUR_SUCCESS)
    return UNUR_ERR_GEN_DATA;

  /* reinit */
  return gen->reinit(gen);

} /* end of unur_dstd_chg_pmfparams() */

/*****************************************************************************/

int
unur_ninv_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (n_params>0) CHECK_NULL(params, UNUR_ERR_NULL);
  
  /* set new parameters in distribution object */
  if (unur_distr_cont_set_pdfparams(gen->distr,params,n_params)!=UNUR_SUCCESS)
    return UNUR_ERR_GEN_DATA;

  /* reinit */
  return gen->reinit(gen);
} /* end of unur_ninv_chg_pdfparams() */

/*****************************************************************************/

#define DISTR gen->distr->data.cont /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

int
unur_srou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_srou_reinit() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );
  
  /* set new parameters in distribution object */
  return unur_distr_cont_set_pdfparams(gen->distr,params,n_params);

} /* end of unur_srou_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_mode( struct unur_gen *gen, double mode )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_srou_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_mode( gen->distr );
} /* end of unur_srou_upd_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_chg_domain( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN | UNUR_DISTR_SET_MASK_DERIVED );
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_srou_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change area below PDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   area  ... area                                                     */
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
  if (area <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"area <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_srou_chg_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_srou_upd_pdfarea( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute area below PDF of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_pdfarea( gen->distr );
} /* end of unur_srou_upd_pdfarea() */

/*---------------------------------------------------------------------------*/

#undef DISTR

/*****************************************************************************/

#define DISTR gen->distr->data.cont /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

int
unur_ssr_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_ssr_reinit() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );
  
  /* set new parameters in distribution object */
  return unur_distr_cont_set_pdfparams( gen->distr, params,n_params );

} /* end of unur_ssr_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_mode( struct unur_gen *gen, double mode )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_ssr_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_mode( gen->distr );
} /* end of unur_ssr_upd_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_ssr_chg_domain( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN | UNUR_DISTR_SET_MASK_DERIVED );
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_ssr_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change area below PDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   area  ... area                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );
  
  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"area <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_ssr_chg_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_upd_pdfarea( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute area below PDF of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SSR, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_pdfarea( gen->distr );
} /* end of unur_ssr_upd_pdfarea() */


/*---------------------------------------------------------------------------*/

#undef DISTR

/*****************************************************************************/

int
unur_tdr_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_tdr_reinit() */

/*****************************************************************************/

int
unur_tdrgw_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* Remark: method TDRGW has been renamed into ARS!                      */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_tdrgw_reinit() */

/*****************************************************************************/

#define DISTR gen->distr->data.cont /* data for distribution in generator object */
#define GEN   ((struct unur_utdr_gen*)gen->datap) /* data for generator object */

/*---------------------------------------------------------------------------*/

int
unur_utdr_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
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
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  return gen->reinit(gen);
} /* end of unur_utdr_reinit() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );
  
  /* set new parameters in distribution object */
  return unur_distr_cont_set_pdfparams( gen->distr, params,n_params );

} /* end of unur_utdr_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_utdr_chg_domain( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;
  GEN->il = left;
  GEN->ir = right;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_utdr_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_mode( struct unur_gen *gen, double mode )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_utdr_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_mode( gen->distr );
} /* end of unur_utdr_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* change area below PDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   area  ... area                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );
  
  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"area <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_utdr_chg_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_upd_pdfarea( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Deprecated call!                                                     */
     /*----------------------------------------------------------------------*/
     /* recompute area below PDF of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );

  return unur_distr_cont_upd_pdfarea( gen->distr );
} /* end of unur_utdr_upd_pdfarea() */

/*---------------------------------------------------------------------------*/

#undef DISTR
#undef GEN

/*****************************************************************************/
#endif   /* USE_DEPRECATED_CODE */
/*---------------------------------------------------------------------------*/
