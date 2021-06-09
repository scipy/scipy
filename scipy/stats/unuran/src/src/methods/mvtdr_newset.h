/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr_newset.c                                               *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    multivariate transformed density rejection                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given (logarithm of the) PDF of a log-concave distribution;          *
 *      produce a value x consistent with its density.                       *
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

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mvtdr_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (distr->dim < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"dim < 2"); return NULL; }

  if ( ! ((DISTR_IN.pdf && DISTR_IN.dpdf) || (DISTR_IN.logpdf && DISTR_IN.dlogpdf)) ) { 
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"d/(log)PDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mvtdr_par) );
  COOKIE_SET(par,CK_MVTDR_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_MVTDR ;   /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_mvtdr_init;

  /* set default values */
  /* minimum number of triangulation steps */
  PAR->steps_min = 5;

  /* maximum number of cones (at least 2^(dim+T_STEPS_MIN) */
  PAR->max_cones = 10000;

  /* bound for splitting cones */
  PAR->bound_splitting = 1.5;

  /** TODO !! **/
  /* move mode to boundary if |mode - boundary| / length < MODE_TO_BOUNDARY */
  /*   PAR->mode_to_boundary = 0.01; */

  return par;

} /* end of unur_mvtdr_new() */

/*---------------------------------------------------------------------------*/

int 
unur_mvtdr_set_stepsmin( struct unur_par *par, int stepsmin )
     /*----------------------------------------------------------------------*/
     /* set minimum number of triangulation step for each starting cone      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter object                           */
     /*   stepsmin ... minimum number of triangulation steps                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );

  /* check new parameter for generator */
  if (stepsmin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"stepsmin < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->steps_min = stepsmin;

  /* changelog */
  par->set |= MVTDR_SET_STEPSMIN;

  return UNUR_SUCCESS;

} /* end of unur_mvtdr_set_stepsmin() */

/*---------------------------------------------------------------------------*/

int
unur_mvtdr_set_boundsplitting( UNUR_PAR *par, double boundsplitting )
     /*----------------------------------------------------------------------*/
     /* set bound for splitting cones. all cones are splitted when the       */
     /* volume below the hat is greater than boundsplitting times the        */
     /* average volume over all cones.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par            ... pointer to parameter object                     */
     /*   boundsplitting ... maximum number of cones                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );

  /* check new parameter for generator */

  /* store date */
  PAR->bound_splitting = boundsplitting;

  /* changelog */
  par->set |= MVTDR_SET_BOUNDSPLITTING;

  return UNUR_SUCCESS;

} /* end of unur_mvtdr_set_boundsplitting() */

/*---------------------------------------------------------------------------*/

int 
unur_mvtdr_set_maxcones( struct unur_par *par, int maxcones )
     /*----------------------------------------------------------------------*/
     /* set maximum number of cones                                          */
     /* (this number is always increased to 2^(dim+stepsmin) )               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter object                           */
     /*   maxcones ... maximum number of cones                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );

  /* check new parameter for generator */
  /* none here: it is always increased to 2^(dim+stepsmin) during init */

  /* store date */
  PAR->max_cones = maxcones;

  /* changelog */
  par->set |= MVTDR_SET_MAXCONES;

  return UNUR_SUCCESS;

} /* end of unur_mvtdr_set_maxcones() */

/*---------------------------------------------------------------------------*/

int
unur_mvtdr_get_ncones( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get number of cones                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   number of cones ... on success                                     */
     /*   0               ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, MVTDR, 0 );

  return GEN->n_cone;
} /* end of unur_mvtdr_get_ncones() */

/*---------------------------------------------------------------------------*/

double
unur_mvtdr_get_hatvol( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get volume below hat                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   volume   ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, MVTDR, INFINITY );

  return GEN->Htot;
} /* end of unur_mvtdr_get_hatvol() */

/*---------------------------------------------------------------------------*/

int
unur_mvtdr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, MVTDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | MVTDR_VARFLAG_VERIFY) : (par->variant & (~MVTDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_mvtdr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_mvtdr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, MVTDR, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cvec_error) 
    return UNUR_FAILURE;

  /* we use a bit in variant */
  gen->variant = (verify) 
    ? (gen->variant | MVTDR_VARFLAG_VERIFY) 
    : (gen->variant & (~MVTDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_mvtdr_chg_verify() */

/*---------------------------------------------------------------------------*/
