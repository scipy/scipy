/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_newset.c                                                 *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a T-concave distribution                                *
 *      produce a value x consistent with its density                        *
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

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_tdr_new( const struct unur_distr* distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL; }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of PDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_tdr_par) );
  COOKIE_SET(par,CK_TDR_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR->guide_factor        = 2.;    /* size of guide table / number of intervals */

  PAR->c_T                 = -0.5;  /* parameter for transformation (-1. <= c < 0.) */

  PAR->starting_cpoints    = NULL;  /* pointer to array of starting points   */
  PAR->n_starting_cpoints  = 30;    /* number of starting points             */
  PAR->percentiles         = NULL;  /* pointer to array of percentiles       */
  PAR->n_percentiles       = 2;     /* number of percentiles                 */
  PAR->retry_ncpoints      = 50;    /* number of cpoints for second trial of reinit */
  PAR->max_ivs             = 100;   /* maximum number of intervals           */
  PAR->max_ratio           = 0.99;  /* bound for ratio  Atotal / Asqueeze    */
  PAR->bound_for_adding    = 0.5;   /* do not add a new construction point in an interval,
				       where ambigous region is too small, i.e. if 
				       area / ((A_hat - A_squeeze)/number of segments) < bound_for_adding */
  PAR->darsfactor          = 0.99;  /* factor for derandomized ARS           */ 
  PAR->darsrule            = 1;     /* rule for finding splitting points in DARS */
 
  par->method   = UNUR_METH_TDR;                 /* method                   */
  par->variant  = ( TDR_VARFLAG_USECENTER |      /* default variant          */
		    TDR_VARFLAG_USEMODE   |
                    TDR_VARIANT_PS );

  par->set      = 0u;               /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default URNG               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_tdr_init;

  return par;

} /* end of unur_tdr_new() */

/*****************************************************************************/

int
unur_tdr_set_cpoints( struct unur_par *par, int n_stp, const double *stp )
     /*----------------------------------------------------------------------*/
     /* set construction points for hat function                             */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*   stp    ... pointer to array of starting points                     */
     /*              (NULL for changing only the number of default points)   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }

  /* store date */
  PAR->starting_cpoints = stp;
  PAR->n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TDR_SET_N_STP | ((stp) ? TDR_SET_STP : 0);

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_reinit_percentiles( struct unur_par *par, int n_percentiles, const double *percentiles )
     /*----------------------------------------------------------------------*/
     /* set percentiles for construction points for hat function             */
     /* and/or its number for re-initialization                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par           ... pointer to parameter for building generator      */
     /*   n_percentiles ... number of percentiles                            */
     /*   percentiles   ... pointer to array of percentiles                  */
     /*                     (NULL for using a rule of thumb)                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check given percentiles */
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }

  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
    
  if (percentiles) {
    /* percentiles must be strictly monontonically increasing */
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }

  /* store date */
  PAR->percentiles = percentiles;
  PAR->n_percentiles = n_percentiles;

  /* changelog */
  par->set |= TDR_SET_N_PERCENTILES | ((percentiles) ? TDR_SET_PERCENTILES : 0);

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_reinit_percentiles() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_chg_reinit_percentiles( struct unur_gen *gen, int n_percentiles, const double *percentiles )
     /*----------------------------------------------------------------------*/
     /* change percentiles for construction points for hat function          */
     /* and/or its number for re-initialization                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen           ... pointer to generator object                      */
     /*   n_percentiles ... number of percentiles                            */
     /*   percentiles   ... pointer to array of percentiles                  */
     /*                     (NULL for using a rule of thumb)                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  /* check given percentiles */
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }

  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
    
  if (percentiles) {
    /* percentiles must be strictly monontonically increasing */
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }

  /* store date */
  GEN->n_percentiles = n_percentiles;
  GEN->percentiles = _unur_xrealloc( GEN->percentiles, n_percentiles * sizeof(double) );
  if (percentiles) {
    memcpy( GEN->percentiles, percentiles, n_percentiles * sizeof(double) );
  }
  else {
    if (n_percentiles == 2) {
      GEN->percentiles[0] = 0.25;
      GEN->percentiles[1] = 0.75;
    }
    else {
      for (i=0; i<n_percentiles; i++ )
	GEN->percentiles[i] = (i + 1.) / (n_percentiles + 1.);
    }
  }

  /* changelog */
  gen->set |= TDR_SET_N_PERCENTILES | ((percentiles) ? TDR_SET_PERCENTILES : 0);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_chg_reinit_percentiles() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_reinit_ncpoints( struct unur_par *par, int ncpoints )
     /*----------------------------------------------------------------------*/
     /* set number of construction points for second trial of reinit         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator           */
     /*   ncpoints ... number of construction points                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check number */
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->retry_ncpoints = ncpoints;

  /* changelog */
  par->set |= TDR_SET_RETRY_NCPOINTS; 

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_reinit_ncpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_chg_reinit_ncpoints( struct unur_gen *gen, int ncpoints )
     /*----------------------------------------------------------------------*/
     /* change number of construction points for second trial of reinit      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen           ... pointer to generator object                      */
     /*   ncpoints ... number of construction points                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  /* check number */
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  GEN->retry_ncpoints = ncpoints;

  /* changelog */
  gen->set |= TDR_SET_RETRY_NCPOINTS; 

  return UNUR_SUCCESS;

} /* end of unur_tdr_chg_reinit_ncpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= TDR_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_max_sqhratio( struct unur_par *par, double max_ratio )
     /*----------------------------------------------------------------------*/
     /* set bound for ratio A(squeeze) / A(hat)                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ratio ... upper bound for ratio to add a new construction point*/
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1.+DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ratio = max_ratio;

  /* changelog */
  par->set |= TDR_SET_MAX_SQHRATIO;

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tdr_get_sqhratio( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio    ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDR, INFINITY );

  return (GEN->Asqueeze / GEN->Atotal);

} /* end of unur_tdr_get_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tdr_get_hatarea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below hat                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDR, INFINITY );

  return GEN->Atotal;

} /* end of unur_tdr_get_hatarea() */

/*---------------------------------------------------------------------------*/

double
unur_tdr_get_squeezearea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below squeeze                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDR, INFINITY );

  return GEN->Asqueeze;

} /* end of unur_tdr_get_squeezearea() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_max_intervals( struct unur_par *par, int max_ivs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ivs   ... maximum number of intervals                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= TDR_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_is_ARS_running( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check whether more points will be added by adaptive rejection        */
     /* sampling                                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if ARS is still running.                                 */
     /*   FALSE ... otherwise.                                               */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, FALSE );
  _unur_check_gen_object( gen, TDR, FALSE );

  return (GEN->n_ivs < GEN->max_ivs) ? TRUE : FALSE;
} /* end of _unur_tdr_is_ARS_running() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_usecenter( struct unur_par *par, int usecenter )
     /*----------------------------------------------------------------------*/
     /* set flag for using center as construction point                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usecenter ... 0 = do not use,  !0 = use                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using center as construction point is the default                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (usecenter) ? (par->variant | TDR_VARFLAG_USECENTER) : (par->variant & (~TDR_VARFLAG_USECENTER));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_usecenter() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_usemode( struct unur_par *par, int usemode )
     /*----------------------------------------------------------------------*/
     /* set flag for using (exact) mode as construction point                */
     /* (this overwrites "use_center"!)                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usemode   ... 0 = do not use,  !0 = use                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using mode as construction point is the default                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (usemode) ? (par->variant | TDR_VARFLAG_USEMODE) : (par->variant & (~TDR_VARFLAG_USEMODE));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_usemode() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_variant_gw( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use original variant with squeezes as proposed by Gilks & Wild       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_GW;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_variant_gw() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_variant_ps( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use squeezes proportional to the hat function                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_PS;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_variant_ps() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_variant_ia( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use squeezes proportional to the hat function together with a        */
     /* composition method that required less uniform random numbers.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_IA;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_variant_ia() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_usedars( struct unur_par *par, int usedars )
     /*----------------------------------------------------------------------*/
     /* set flag for using DARS (derandomized adaptive rejection sampling).  */
     /* additionally the rule for splitting intervals can be set.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usedars   ... 0 = do not use,  1-3 = use DARS with given rule      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using not using DARS is the default                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (usedars < 0 || usedars > 3) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid rule for DARS");
    return UNUR_ERR_PAR_SET;
  }
    
  /* set rule for DARS */
  PAR->darsrule = usedars;

  /* we use a bit in variant */
  par->variant = (usedars) ? (par->variant | TDR_VARFLAG_USEDARS) : (par->variant & (~TDR_VARFLAG_USEDARS));

  /* changelog */
  par->set |= TDR_SET_USE_DARS;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_usedars() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_darsfactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for derandomized adaptive rejection sampling              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... parameter for DARS                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
    
  /* store date */
  PAR->darsfactor = factor;

  /* changelog */
  par->set |= TDR_SET_DARS_FACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_darsfactor() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_c( struct unur_par *par, double c )
     /*----------------------------------------------------------------------*/
     /* set parameter c for transformation T_c                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   c    ... parameter c                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* check new parameter for generator */
  if (c > 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c > 0");
    return UNUR_ERR_PAR_SET;
  }
  /** TODO: ... **/
/*    if (c <= -1.) { */
/*      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c <= -1 only if domain is bounded. Use `TABL' method then."); */
/*      return 0; */
/*    } */
  /** TODO: ... **/
  if (c < -0.5) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"c < -0.5 not implemented yet");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_iszero(c) && c > -0.5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
    c = -0.5;
  }
    
  /* store date */
  PAR->c_T = c;

  /* changelog */
  par->set |= TDR_SET_C;

  return UNUR_SUCCESS;

} /* end of unur_tdr_set_c() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TDR_VARFLAG_VERIFY) : (par->variant & (~TDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  /* we use a bit in variant */
  gen->variant = (verify) 
    ? (gen->variant | TDR_VARFLAG_VERIFY) 
    : (gen->variant & (~TDR_VARFLAG_VERIFY));

  /* sampling routines */
  SAMPLE = _unur_tdr_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_chg_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_pedantic( struct unur_par *par, int pedantic )
     /*----------------------------------------------------------------------*/
     /* turn pedantic mode on/off                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   pedantic ... 0 = no pedantic mode, !0 = use pedantic mode          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   pedantic is the default                                            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );

  /* we use a bit in variant */
  par->variant = (pedantic) ? (par->variant | TDR_VARFLAG_PEDANTIC) : (par->variant & (~TDR_VARFLAG_PEDANTIC));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdr_set_pedantic() */

/*---------------------------------------------------------------------------*/

int 
unur_tdr_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
     /*                                                                      */
     /* This call does not work for variant IA (immediate acceptance).       */
     /* In this case it switches to variant PS!!                             */
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
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  double Umin, Umax;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  /* we have to disable adaptive rejection sampling */
  if (GEN->max_ivs > GEN->n_ivs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"adaptive rejection sampling disabled for truncated distribution");
    GEN->max_ivs = GEN->n_ivs;
  }

  /* we cannot use immadate acceptance (IA), switch to variant PS instead */
  if ((gen->variant & TDR_VARMASK_VARIANT) == TDR_VARIANT_IA) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"cannot use IA for truncated distribution, switch to PS");
    /* change variante flag */
    gen->variant = (gen->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_PS;
    /* change sampling routine */
    SAMPLE = (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
  }

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain not subset of domain");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain not subset of domain");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* compute CDF at x (with respect to given domain of distribution) */
  Umin = _unur_tdr_eval_cdfhat(gen,left);
  Umax = (right < DISTR.domain[1]) ? _unur_tdr_eval_cdfhat(gen,right) : 1.;

  /* check result */
  if (Umin > Umax) {
    /* this is a serios error that should not happen */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  if (_unur_FP_equal(Umin,Umax)) {
    /* CDF values very close */
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      /* this is very bad */
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }

  /* set bounds for truncated domain and for U (CDF) */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_tdr_chg_truncated() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_eval_cdfhat( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of hat at x (i.e. \int_{-\infty}^x hat(t) dt)           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point at which hat(x) has to be computed                   */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF of hat(x) or                                                   */
     /*   INFINITY in case of error                                          */
     /*                                                                      */
     /* Important:                                                           */
     /*   If gen is a generator object for variant IA (immediate acceptance) */
     /*   then it is treated like variant PS (proportional squeeze)!         */
     /*   This is necessary since variant IA is not a pure rejection         */
     /*   algorithm, but a composition method.                               */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double Aint;
  double cdf;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);

  /* the easy case: left and right boundary of domain */
  if (x <= DISTR.domain[0]) return 0.;
  if (x >= DISTR.domain[1]) return 1.;

  /* there are differencies between variant GW and variant PS */
  switch (gen->variant & TDR_VARMASK_VARIANT) {

  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */

    /* find interval (sequential search) */
    for (iv = GEN->iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
      /* iv->x is left construction point of interval */
      if (x < iv->next->x) break;
    }

    if (iv->next == NULL)
      /* right boundary of domain */
      return 1.;

    /* now iv->x < x <= iv->next->x */

    /* compute are below hat between construction point of tangent and x. */
    /* we have to cases on either side of the intersection point.         */

    if (x < iv->ip) {
      /* left h.s. of intersection point */
      Aint = _unur_tdr_interval_area( gen, iv, iv->dTfx, x);
      if (!_unur_isfinite(Aint)) { 
	/* this should not happen */
	Aint = 0.;
      }
      /* Notice: iv->prev->Acum == iv->Acum - iv->Ahat; */
      cdf = (iv->prev) ? iv->prev->Acum + Aint : Aint;
    }
    else {
      /* right h.s. of intersection point */
      Aint = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, x);
      if (!_unur_isfinite(Aint)) { 
	/* this should not happen */
	Aint = 0.;
      }
      cdf = iv->Acum - Aint;
      if (cdf < 0.) return 0.;
    }

    /* normalize to one (and mind round-off errors) */
    cdf /= GEN->Atotal;
    return ((cdf > 1.) ? 1. : cdf);

    
  case TDR_VARIANT_IA:    /* immediate acceptance */
    /* See comment above */

  case TDR_VARIANT_PS:    /* proportional squeeze */

    /* find interval (sequential search) */
    for (iv = GEN->iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
      if (x <= iv->next->ip) break;
    }
    if (iv->next == NULL)
      /* right boundary of domain */
      return 1.;

    /* now iv->ip < x <= iv->next->ip */

    /* area below hat between construction point and x */
    Aint = _unur_tdr_interval_area( gen, iv, iv->dTfx, x);
    if (!_unur_isfinite(Aint)) { 
      /* this should not happen */
      Aint = 0.;
    }

    /* compute CDF of hat */
    cdf = ((x>iv->x) ? Aint : -Aint) + iv->Acum - iv->Ahatr;

    /* normalize to one (and mind round-off errors) */
    if (cdf < 0.) return 0.;
    cdf /= GEN->Atotal;
    return ((cdf > 1.) ? 1. : cdf);
    
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }

} /* end of _unur_tdr_eval_cdfhat() */

/*****************************************************************************/

