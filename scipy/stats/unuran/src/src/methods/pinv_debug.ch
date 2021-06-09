/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_debug.c                                                 *
 *                                                                           *
 *   Routines for printing debugging information.                            * 
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG filebefore setup starts.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = PINV (Polynomial interpolation based INVerse CDF)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_pinv_sample\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: order of polynomial = %d",gen->genid,GEN->order);
  _unur_print_if_default(gen,PINV_SET_ORDER);
  if(gen->set & PINV_SET_ORDER_COR) fprintf(LOG,"  [corrected]");
  fprintf(LOG,"\n%s: smoothness = %d",gen->genid,GEN->smooth);
  _unur_print_if_default(gen,PINV_SET_SMOOTH);
  if(gen->set & PINV_SET_SMOOTH_COR) fprintf(LOG,"  [corrected]");
  fprintf(LOG,"\n%s: u-resolution = %g",gen->genid,GEN->u_resolution);
  _unur_print_if_default(gen,PINV_SET_U_RESOLUTION);
  fprintf(LOG,"\n%s: maximum number of subintervals = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,PINV_SET_MAX_IVS);

  fprintf(LOG,"\n%s: variant = ",gen->genid);
  if (gen->variant & PINV_VARIANT_PDF)
    fprintf(LOG,"use PDF + Lobatto integration");
  else
    fprintf(LOG,"use CDF");
  _unur_print_if_default(gen,PINV_SET_VARIANT);
  fprintf(LOG,"\n");

  fprintf(LOG,"%s: use Chebyshev points in %s scale",gen->genid,
	  (gen->variant & PINV_VARIANT_UPOINTS) ? "u" : "x");
  _unur_print_if_default(gen,PINV_SET_UPOINTS);
  fprintf(LOG,"\n");

  fprintf(LOG,"%s: keep table of CDF values = %s",gen->genid,
	  (gen->variant & PINV_VARIANT_KEEPCDF) ? "on" : "off");
  _unur_print_if_default(gen,PINV_SET_KEEPCDF);
  fprintf(LOG,"\n");

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} /* end of _unur_pinv_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init( const struct unur_gen *gen, int ok )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   ok  ... exitcode of init call                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: domain of computation = [%g,%g]\n",gen->genid, GEN->bleft,GEN->bright);
  fprintf(LOG,"%s: Umin = 0 [fixed], Umax = %18.16g",gen->genid, GEN->Umax);
  if (_unur_FP_approx(GEN->Umax,1.))
    fprintf(LOG,",  1-Umax = %g",1.-GEN->Umax);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: # Intervals = %d\n",gen->genid,GEN->n_ivs);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_pinv_debug_intervals(gen);

  fprintf(LOG,"%s: initialization %s\n",gen->genid,((ok)?"successful":"*** FAILED ***")); 
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_pinv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_relevant_support (const struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* print relevant domain                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: search for boundaries: left = %s, right = %s\n",gen->genid,
	  GEN->sleft ? "TRUE" : "FALSE", GEN->sright ? "TRUE" : "FALSE");
  fprintf(LOG,"%s: relevant domain = (%g,%g)   [i.e. where PDF > threshold]\n",gen->genid,
	  GEN->bleft,GEN->bright);
  fprintf(LOG,"%s: possible support of distribution = (%g,%g)\n",gen->genid,
	  GEN->dleft,GEN->dright);

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_pinv_debug_relevant_support() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_pdfarea (const struct unur_gen *gen, int approx)
     /*----------------------------------------------------------------------*/
     /* print estimated area below PDF                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   approx ... whether this is an approximate value or accurat up to   */
     /*              requested tolerance                                     */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: area below PDF %s = %19.16g\n",gen->genid,
	  approx ? "(approx.)" : "(accurate)", GEN->area);

  if (GEN->aCDF)
    _unur_lobatto_debug_table(GEN->aCDF, gen, (gen->debug & PINV_DEBUG_ITABLE));

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} /* end of _unur_pinv_debug_pdfarea() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_computational_domain (const struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* print computational domain                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: computational domain = (%g,%g)\n",gen->genid,
	    GEN->bleft,GEN->bright);

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_pinv_debug_computational_domain() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print table of intervals into LOG file                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int n;
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (gen->debug & PINV_DEBUG_TABLE) {
    for (n=0; n<=GEN->n_ivs; n++) {
      fprintf(LOG,"%s: [%3d] xi = %.14g, cdfi = %.14g\n",gen->genid,
	      n, GEN->iv[n].xi, GEN->iv[n].cdfi);

#ifdef PINV_DEVEL
      { 
	int j;
	/* coefficients for constructing Newton interpolation */
	fprintf(LOG,"%s:\tui = %.14g",gen->genid, GEN->iv[n].ui[0]);
	for (j=1; j<GEN->order; j++)
	  fprintf(LOG,", %.14g",GEN->iv[n].ui[j]);
	
	fprintf(LOG,"\n%s:\tzi = %.14g",gen->genid, GEN->iv[n].zi[0]);
	for (j=1; j<GEN->order; j++)
	  fprintf(LOG,", %.14g",GEN->iv[n].zi[j]);

	fprintf(LOG,"\n");
      }	
#endif
    }
  }

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_pinv_debug_intervals() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_create_table (const struct unur_gen *gen,
			       int iter, int n_incr_h, int n_decr_h, int n_use_linear)
     /*----------------------------------------------------------------------*/
     /* print data that have been collected while creating polynomials.      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iter     ... total number of iterations                            */
     /*   n_incr_h ... number of steps where h is increased                  */
     /*   n_decr_h ... number of steps where h is decreased                  */
     /*   n_use_linear .. number of steps where linear interpolation is used */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: Create interpolating polynomials:\n",gen->genid);
  fprintf(LOG,"%s:\t# iterations   = %d\n",gen->genid,iter);
  fprintf(LOG,"%s:\t# increasing h = %d  (%g%%)\n",gen->genid,n_incr_h,(100.*n_incr_h)/iter);
  fprintf(LOG,"%s:\t# decreasing h = %d  (%g%%)\n",gen->genid,n_decr_h,(100.*n_decr_h)/iter);
  fprintf(LOG,"%s:\t# linear       = %d  (%g%%)\n",gen->genid,n_use_linear,(100.*n_use_linear)/iter);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);
} /* end of _unur_pinv_debug_create_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
