/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_newset.ch                                               *
 *                                                                           *
 *   Routines for creating and changing parameter objects.                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_ninv_new( const struct unur_distr *distr )
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

  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }

  /* if default variant is Newton's method, then we also need the PDF ! */

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_ninv_par) );
  COOKIE_SET(par,CK_NINV_PAR);

  /* copy input */
  par->distr   = distr;            /* pointer to distribution object         */

  /* set default values */
  PAR->max_iter  = 100;            /* maximal number of iterations           */
  PAR->x_resolution = 1.0e-8;      /* maximal tolerated relative x-error     */
  PAR->u_resolution = -1.;         /* maximal tolerated u-error -- DISABLED  */

  /* starting points for numerical inversion */
  PAR->s[0]      = 0.0;     /* regula falsi: left boundary of starting interval
			      newton: starting point                         */
  PAR->s[1]      = 0.0;     /* regula falsi: right boundary of starting interval
			      newton: not used                               */
  /* If s1 and s2 are equal a defaults are used, see below */

  PAR->table_on  = FALSE;   /* Do not use a table for starting points
			      by default.                                    */
 
  par->method   = UNUR_METH_NINV;          /* method and default variant     */
  par->variant  = NINV_VARFLAG_REGULA;     /* Use regula falsi as default 
					      method                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_ninv_init;

  return par;

} /* end of unur_ninv_new() */

/*****************************************************************************/

int
unur_ninv_set_usenewton( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use Newton's method                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   /* use regula falsi instead  */
    return UNUR_ERR_DISTR_REQUIRED;
 }

  /* store date */
  par->variant = NINV_VARFLAG_NEWTON;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_usenewton() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_useregula( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use regula falsi                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* store date */
  par->variant = NINV_VARFLAG_REGULA;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_useregula() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_usebisect( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use bisection method                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* store date */
  par->variant = NINV_VARFLAG_BISECT;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_usebisect() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_max_iter( struct unur_par *par, int max_iter )
     /*----------------------------------------------------------------------*/
     /* set number of maximal iterations                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximal iterations");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_iter = max_iter;

  /* changelog */
  par->set |= NINV_SET_MAX_ITER;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_max_iter() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_max_iter( struct unur_gen *gen, int max_iter )
     /*----------------------------------------------------------------------*/
     /* change number of maximal iterations                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(gen->genid, UNUR_ERR_PAR_SET, "maximal iterations");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  GEN->max_iter = max_iter;

  /* changelog */
  gen->set |= NINV_SET_MAX_ITER;

  /* o.k.  */
  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_max_iter() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated relative x-error                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   x_resolution ... x-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (x_resolution > 0. && x_resolution < 2.*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }

  /* store date */
  PAR->x_resolution = x_resolution;

  /* changelog */
  par->set |= NINV_SET_X_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_x_resolution( struct unur_gen *gen, double x_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated relative x-error                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   x_resolution ... x-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (x_resolution > 0. && x_resolution < DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }

  /* store date */
  GEN->x_resolution = x_resolution;

  /* changelog */
  gen->set |= NINV_SET_X_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated u-error                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... u-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }

  /* store date */
  PAR->u_resolution = u_resolution;

  /* changelog */
  par->set |= NINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_u_resolution( struct unur_gen *gen, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated u-error                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   u_resolution ... u-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }

  /* store date */
  GEN->u_resolution = u_resolution;

  /* changelog */
  gen->set |= NINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_start( struct unur_par *par, double s1, double s2 )
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   other methods: s1, s2       boundary of starting interval          */
     /* arguments that are not used by method are ignored.                   */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   s1    ... left starting point                                      */
     /*   s2    ... right starting point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* store date */
  if ( s1 <= s2 ) {
     PAR->s[0] = s1;
     PAR->s[1] = s2;
  }
  else {
     PAR->s[0] = s2;
     PAR->s[1] = s1;
  }

  /* changelog */
  par->set |= NINV_SET_START;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_start() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_start( struct unur_gen *gen, double s1, double s2 )
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /* arguments that are used by method are ignored.                       */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   s1    ... left starting point                                      */
     /*   s2    ... right starting point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* store date */
  if ( s1 <= s2 ) {
     GEN->s[0] = s1;
     GEN->s[1] = s2;
  }
  else {
     GEN->s[0] = s2;
     GEN->s[1] = s1;
  }

  /* disable table (we now want to use only two starting points) */
  GEN->table_on = FALSE;

  /* compute these points */
  _unur_ninv_compute_start(gen);
  /* this call should not fail */

  /* changelog */
  gen->set |= NINV_SET_START;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_start() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_table( struct unur_par *par, int tbl_pnts )
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* the function unur_ninv_set_start() is overruled                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   tbl_pnts ... number of table points                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  PAR->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  PAR->table_on = TRUE;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_table() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_table( struct unur_gen *gen, int tbl_pnts )
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* set somewhere else will be ignored                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   tbl_pnts ... number of table points                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /*  free(GEN->table);   not freed, because realloc() is used */ 
  /*  free(GEN->f_table); not freed, because realloc() is used */
  GEN->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;

  result = _unur_ninv_create_table(gen); 

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_CHG) 
    if (result==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
#endif
  
  return result;

} /* end of unur_ninv_chg_table() */

/*---------------------------------------------------------------------------*/

int 
unur_ninv_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
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
  double Umin, Umax;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* set bounds of U -- in respect to given bounds */
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;

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

  /* copy new boundaries into generator object */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_CHG) 
    _unur_ninv_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_ninv_chg_truncated() */

/*---------------------------------------------------------------------------*/
