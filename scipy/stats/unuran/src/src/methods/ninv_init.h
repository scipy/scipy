/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_init.ch                                                 *
 *                                                                           *
 *   Routines for initialization and deletion of generator objects.          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_ninv_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* check variant */
  if (par->variant == NINV_VARFLAG_NEWTON && ! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   /* use regula falsi instead  */
  }

  /* create a new empty generator object */    
  gen = _unur_ninv_create(par);
  _unur_par_free(par);
  if (!gen) { return NULL; }
  
  /* check parameters */
  if (_unur_ninv_check_par(gen) != UNUR_SUCCESS) {
    _unur_ninv_free(gen); return NULL;
  }

  /* compute starting points for numerical inversion */
  if (GEN->table_on) {
    /* use a table */
    if (_unur_ninv_create_table(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
  else {
    /* use given points or percentiles */
    if (_unur_ninv_compute_start(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ninv_debug_init(gen);
#endif

  return gen;

} /* end of _unur_ninv_init() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_reinit( struct unur_gen *gen )
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
  if ( (rcode = _unur_ninv_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* compute normalization constant for standard distribution */
  if (DISTR.upd_area != NULL)
    if ((DISTR.upd_area)(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cannot compute normalization constant");
      return UNUR_ERR_GEN_DATA;
    }

  /* regenerate table */
  if (GEN->table != NULL)
    rcode = _unur_ninv_create_table(gen);

  else /* or compute starting points */
    rcode = unur_ninv_chg_start( gen, 0., 0. );

  /* (re)set sampling routine */
  SAMPLE = _unur_ninv_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_REINIT) {
    _unur_distr_cont_debug( gen->distr, gen->genid );
    if (rcode==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
  }
#endif

  return UNUR_SUCCESS;
} /* end of _unur_ninv_reinit() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_ninv_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NINV_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_ninv_getSAMPLE(gen);
  gen->destroy = _unur_ninv_free;
  gen->clone = _unur_ninv_clone;
  gen->reinit = _unur_ninv_reinit;

  /* copy parameters into generator object */
  GEN->max_iter = PAR->max_iter;      /* maximal number of iterations          */
  GEN->x_resolution = PAR->x_resolution; /* maximal tolerated relative x-error */
  GEN->u_resolution = PAR->u_resolution; /* maximal tolerated u-error          */
  GEN->table_on = PAR->table_on;      /* useage of table for starting points   */
  GEN->table_size = PAR->table_size;  /* number of points for table            */
  GEN->s[0] = PAR->s[0];              /* starting points                       */
  GEN->s[1] = PAR->s[1];

  /* init pointer */
  GEN->table = NULL;
  GEN->f_table = NULL;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_ninv_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_ninv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_check_par( struct unur_gen *gen )
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
  /* checking x-error or checking u-error must be enabled */
  if ( GEN->x_resolution < 0. && GEN->u_resolution < 0. ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"both x-resolution and u-resolution negativ. using defaults.");
    GEN->x_resolution = 1.e-8;
  }

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* set bounds of U -- in respect to given bounds                          */
  GEN->CDFmin = GEN->Umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN->CDFmax = GEN->Umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;

  if (_unur_FP_greater(GEN->CDFmin, GEN->CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return UNUR_ERR_GEN_DATA;
  }

  return UNUR_SUCCESS;
} /* end of _unur_ninv_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ninv_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_ninv_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NINV_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy additional data for generator object */
  if (GEN->table) {
    CLONE->table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->table, GEN->table, GEN->table_size * sizeof(double) );
    CLONE->f_table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->f_table, GEN->f_table, GEN->table_size * sizeof(double) );
  }

  return clone;

#undef CLONE
} /* end of _unur_ninv_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free tables */
  if (GEN->table)   free(GEN->table);
  if (GEN->f_table) free(GEN->f_table);

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_ninv_free() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_ninv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a table for starting points and a table with                  */
     /* the corresponding function values                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;
  double x;
  int table_size = GEN->table_size;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);

  GEN->table    = _unur_xrealloc( GEN->table,   table_size * sizeof(double));
  GEN->f_table  = _unur_xrealloc( GEN->f_table, table_size * sizeof(double));

  /* get arbitrary points */
  GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
  GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
  GEN->CDFs[0]  = CDF(GEN->s[0]);
  GEN->CDFs[1]  = CDF(GEN->s[1]);

  /* table can't be used to calculate itself */
  GEN->table_on = FALSE;
  
  /* calculation of the tables   */

  /* left and right boundary */
  GEN->table[0]              = DISTR.domain[0];
  GEN->f_table[0]            = GEN->CDFmin;    /* CDF(DISTR.domain[0]) */
  GEN->table[table_size-1]   = DISTR.domain[1];
  GEN->f_table[table_size-1] = GEN->CDFmax;    /* CDF(DISTR.domain[1]) */
  
  /* all the other points. we compute these points from boundary to center */
  for (i=1; i<table_size/2; i++){

    /* compute table point and CDF at table point */
    x = GEN->CDFmin + i * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[i]   = _unur_ninv_regula(gen,x);
    GEN->f_table[i] = CDF(GEN->table[i]);

    x = GEN->CDFmin + (table_size-i-1) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size-1-i] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size-1-i] = CDF(GEN->table[table_size-1-i]);

    /* set new starting points for computing table points in next step */
    if (GEN->table[i] > -INFINITY) {
      GEN->s[0] = GEN->table[i];
      GEN->CDFs[0] = GEN->f_table[i];
    }
    if (GEN->table[table_size-1-i] < INFINITY) {
      GEN->s[1] = GEN->table[table_size-1-i];
      GEN->CDFs[1] = GEN->f_table[table_size-1-i];
    }

  }  /* end of for()                                                     */

  /* the median point (only if table_size is odd) */
  if (table_size & 1) { 
    x = GEN->CDFmin + (table_size/2) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size/2] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size/2] = CDF(GEN->table[table_size/2]);
  }  

  /* calculation of tables finished  */

  GEN->table_on = TRUE;

  /* o.k. */
  return UNUR_SUCCESS;

}  /* end of _unur_ninv_create_table() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_compute_start( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting points for numerical inversion.                     */
     /* use 5% percentiles.                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double u;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);

  if( GEN->table_on )
    /* we have a table --> nothing to do */
    return UNUR_SUCCESS;

  if( !_unur_FP_same(GEN->s[0],GEN->s[1]) ) {
    /* use given starting points (indicated by s[0] != s[1]) --> nothing to do */
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    return UNUR_SUCCESS;
  }

  switch (gen->variant) {

  case NINV_VARFLAG_BISECT:
  case NINV_VARFLAG_REGULA:

    /* get arbitrary points */
    GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);    

    /* left percentile */
    u = GEN->CDFmin + 0.5*(1.-INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[0] = _unur_ninv_regula(gen,u);
    GEN->CDFs[0] = CDF(GEN->s[0]);

    /* right percentile */
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    u = GEN->CDFmin + 0.5*(1.+INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[1] = _unur_ninv_regula(gen,u);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    
    break;    /* case REGULA end */

  case NINV_VARFLAG_NEWTON:

    /* get arbitrary points */
    GEN->s[0] = _unur_max( DISTR.domain[0], -9.987655 );
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]); 
    GEN->CDFs[1] = CDF(GEN->s[1]);

    /* median */
    u = 0.5 * (GEN->CDFmin + GEN->CDFmax);
    GEN->s[0] = _unur_ninv_regula(gen, u); /* reglua more stable */
    GEN->CDFs[0] = CDF(GEN->s[0]);

    break;    /* case NEWTON end */

  default:

    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;

  }  /* end of switch  */

  /* o.k. */
  return UNUR_SUCCESS;

}  /* end of _unur_ninv_compute_start() */

/*---------------------------------------------------------------------------*/
