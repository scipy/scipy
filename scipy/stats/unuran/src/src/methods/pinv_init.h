/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_init.c                                                  *
 *                                                                           *
 *   Routines for initialization and deletion of generator objects.          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008-2010 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_pinv_init( struct unur_par *par )
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
  double lfc;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_PINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_PINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_pinv_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_pinv_check_par(gen) != UNUR_SUCCESS) {
    _unur_pinv_free(gen); return NULL;
  }

  /* compute rescaling factor for PDF */
  /* (only used when logPDF is given) */
  if (DISTR.logpdf != NULL && (gen->variant & PINV_VARIANT_PDF) ) {
    lfc = UNUR_INFINITY;

    /* use mode if available */
    if ( (gen->distr->set & UNUR_DISTR_SET_MODE) &&
	 !_unur_FP_less(DISTR.mode,DISTR.domain[0]) &&
	 !_unur_FP_greater(DISTR.mode,DISTR.domain[1]) ) {
      lfc = (DISTR.logpdf)(DISTR.mode,gen->distr);
    }

    /* use center otherwise (or if logPDF(mode)==INFINITY) */ 
    if (!_unur_isfinite(lfc))
      lfc = (DISTR.logpdf)(DISTR.center,gen->distr);

    /* rescaling results in more evaluations of the logPDF, */
    /* when the logPDF is approximately 0.                  */
    /* so we only rescale the logPDF when it is too small.  */
    if (lfc < -3.)
      GEN->logPDFconstant = lfc;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_pinv_debug_init_start(gen);
#endif

  /* 1. Preprocessing:                                     */
  /*   find interval for computing Newton interpolation */
  if (_unur_pinv_preprocessing(gen) != UNUR_SUCCESS) {
    /* preprocessing failed */
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }
  
  /* compute table for Newton interpolation */
  if (_unur_pinv_create_table(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }

  /* we do not need the table with CDF values any more. */
  /* thus we may free the allocated memory.             */
  if (! (gen->variant & PINV_VARIANT_KEEPCDF))
    _unur_lobatto_free(&(GEN->aCDF));

  /* make guide table */
  _unur_pinv_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_pinv_debug_init(gen,TRUE);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_pinv_init() */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_reinit( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* re-initialize (existing) generator.                                  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int rcode; */

/*   /\* check parameters *\/ */
/*   if ( (rcode = _unur_pinv_check_par(gen)) != UNUR_SUCCESS) */
/*     return rcode; */

/*   return UNUR_SUCCESS; */
/* } /\* end of _unur_pinv_reinit() *\/ */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_pinv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_PINV_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_pinv_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_PINV_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_pinv_getSAMPLE(gen);
  gen->destroy = _unur_pinv_free;
  gen->clone = _unur_pinv_clone;
  /* gen->reinit = _unur_pinv_reinit; */

  /* copy parameters into generator object */
  GEN->order = PAR->order;            /* order of polynomial                 */
  GEN->smooth = PAR->smooth;          /* smoothness parameter                */
  GEN->u_resolution = PAR->u_resolution; /* maximal error in u-direction     */
  GEN->bleft_par  = PAR->bleft;          /* border of computational domain   */
  GEN->bright_par = PAR->bright;
  GEN->sleft  = PAR->sleft;              /* whether to search for boundary   */
  GEN->sright = PAR->sright;
  GEN->max_ivs = PAR->max_ivs;           /* maximum number of subintervals   */

  /* initialize variables */
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->dleft = -INFINITY;
  GEN->dright = INFINITY;
  GEN->Umax = 1.;
  GEN->iv = NULL;
  GEN->n_ivs = -1;        /* -1 indicates that there are no intervals at all */
  GEN->guide_size = 0; 
  GEN->guide = NULL;
  GEN->area = DISTR.area; /* we use the value in the distribution object as first guess */
  GEN->logPDFconstant = 0.;   /* rescaling constant for logPDF                  */
  GEN->aCDF = NULL;           /* pointer to approximate CDF */

  /* allocate maximal array of intervals */
  /* [ Maybe we could move this into _unur_pinv_interval() ] */
  GEN->iv = _unur_xmalloc(GEN->max_ivs * sizeof(struct unur_pinv_interval) );

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_pinv_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_pinv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_check_par( struct unur_gen *gen )
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
  /* smoothness parameter */
  switch (GEN->smooth) {
  case 2:
    if (GEN->order < 5) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be >= 5 when smoothness equals 2");
      GEN->order = 5;
      gen->set |= PINV_SET_ORDER_COR;
    }
    if (GEN->order % 3 != 2) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be 2 mod 3 when smoothness equals 2");
      GEN->order = 2 + 3 * (GEN->order / 3);
      gen->set |= PINV_SET_ORDER_COR;
    }

    if (DISTR.pdf == NULL || DISTR.dpdf == NULL) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_REQUIRED,"PDF or dPDF --> try smoothness=1 instead");
      GEN->smooth = 1;
      gen->set |= PINV_SET_SMOOTH_COR;
    }
    else {
      break;
    }

  case 1:
    if (GEN->order % 2 != 1) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be odd when smoothness equals 1");
      GEN->order += 1;
      gen->set |= PINV_SET_ORDER_COR;
    }
    if (DISTR.pdf == NULL) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_REQUIRED,"PDF --> use smoothness=0 instead");
      GEN->smooth = 0;
      gen->set |= PINV_SET_SMOOTH_COR;
    }
    else { 
      break;
    }

  case 0:
    break;

  default:
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"smoothness must be 0, 1, or 2");
    GEN->smooth = 0; /* use default */
  }

  /* points for searching computational domain */
  GEN->bleft = _unur_max(GEN->bleft_par,DISTR.domain[0]);
  GEN->bright = _unur_min(GEN->bright_par,DISTR.domain[1]);

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* domain of distribution (used when x with PDF(x)=0 are found) */
  GEN->dleft =  DISTR.domain[0];
  GEN->dright =  DISTR.domain[1];

  /* center of distribution */
  DISTR.center = unur_distr_cont_get_center(gen->distr);
  if (DISTR.center < GEN->dleft || DISTR.center > GEN->dright) {
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,
		"center moved into domain of distribution");
    DISTR.center = _unur_max(DISTR.center,GEN->dleft);
    DISTR.center = _unur_min(DISTR.center,GEN->dright);
  }

  /* check center of distribution */
  if ( (gen->variant & PINV_VARIANT_PDF ) &&
       /* (_unur_pinv_search_center(gen) != UNUR_SUCCESS) ) { */
       (_unur_distr_cont_find_center(gen->distr) != UNUR_SUCCESS) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"PDF(center) <= 0.");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* there is no table of CDF values when the CDF is given */
  if (! (gen->variant & PINV_VARIANT_PDF))
    gen->variant &= ~PINV_VARIANT_KEEPCDF;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_pinv_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_pinv_gen*)clone->datap)

  struct unur_gen *clone;
  int i;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* we do not need the table with CDF values */
  CLONE->aCDF = NULL;

  /* copy coefficients for Newton polynomial */
  CLONE->iv =  _unur_xmalloc((GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  memcpy( CLONE->iv, GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );

  for(i=0; i<=GEN->n_ivs; i++) {
    CLONE->iv[i].ui = _unur_xmalloc( GEN->order * sizeof(double) );
    CLONE->iv[i].zi = _unur_xmalloc( GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].ui, GEN->iv[i].ui, GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].zi, GEN->iv[i].zi, GEN->order * sizeof(double) );
  }

  /* copy guide table */
  CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_pinv_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  int i;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free guide table */
  if (GEN->guide) free (GEN->guide);

  /* free array for subintervals of adaptive Gauss-Lobatto integration */
  _unur_lobatto_free(&(GEN->aCDF));

  /* free tables of coefficients of interpolating polynomials */
  if (GEN->iv) {
    for(i=0; i<=GEN->n_ivs; i++){
      free(GEN->iv[i].ui);
      free(GEN->iv[i].zi);
    }
    free (GEN->iv);
  }

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_pinv_free() */

/*****************************************************************************/

int
_unur_pinv_make_guide_table (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i,j, imax;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  GEN->guide_size = (int) (GEN->n_ivs * PINV_GUIDE_FACTOR);
  if (GEN->guide_size <= 0) GEN->guide_size = 1;
  GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) );

  /* maximum index for array of data */
  imax = GEN->n_ivs;

  /* create guide table */
  i = 0;
  GEN->guide[0] = 0;
  for( j=1; j<GEN->guide_size ;j++ ) {
    while(GEN->iv[i+1].cdfi/GEN->Umax < j/(double)GEN->guide_size && i < imax)
      i++;
    if (i >= imax) break;
    GEN->guide[j]=i;
  }

  /* check i */
  i = _unur_min(i,imax);

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = i;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_make_guide_table() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_eval_PDF (double x, struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* call to PDF.                                                         */
     /* if PDF(x) == INFINITY (i.e., x is a pole of the PDF)                 */
     /* we compute PDF(x+dx) for a small dx instead.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x   ... argument of PDF                                            */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr = gen->distr;
  double fx, dx;
  int i;

  for (i=1; i<=2; i++) {
    /* we try two times:                               */
    /* one time for given x, the second time for x+dx. */

    /* compute PDF(x) */
    if (DISTR.logpdf != NULL)
      fx = exp((DISTR.logpdf)(x,distr) - GEN->logPDFconstant);
    else
      fx = (DISTR.pdf)(x,distr);

    /* check for pole (i.e., PDF(x)==INFINITY) */
    if (fx >= INFINITY) {
      /* pole at x --> move x slightly */
      dx = 2.*fabs(x)*DBL_EPSILON;
      dx = _unur_max(dx,2.*DBL_MIN);
      x += ((x - GEN->bleft) < (GEN->bright - x)) ? dx : -dx; 
    }

    else /* not a pole */
      break;
  }

  return fx;
} /* end of _unur_pinv_eval_PDF() */

/*---------------------------------------------------------------------------*/
