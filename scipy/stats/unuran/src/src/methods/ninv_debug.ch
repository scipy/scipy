/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_debug.ch                                                *
 *                                                                           *
 *   Routines for printing debugging information.                            * 
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
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
_unur_ninv_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ninv (numerical inversion of CDF)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_ninv_sample",gen->genid);
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    fprintf(LOG,"_newton\n");
    break;
  case NINV_VARFLAG_BISECT:
    fprintf(LOG,"_bisect\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    fprintf(LOG,"_regula\n");
    break;
  }
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: u-resolution = ",gen->genid);
  if (GEN->u_resolution < 0.) 
    fprintf(LOG,"[disabled]");
  else
  fprintf(LOG,"%g",GEN->u_resolution);
  _unur_print_if_default(gen,NINV_SET_U_RESOLUTION);

  fprintf(LOG,"\n%s: x-resolution = ",gen->genid);
  if (GEN->x_resolution < 0.) 
    fprintf(LOG,"[disabled]");
  else
  fprintf(LOG,"%g",GEN->x_resolution);
  _unur_print_if_default(gen,NINV_SET_X_RESOLUTION);

  fprintf(LOG,"\n%s:\n",gen->genid);

  _unur_ninv_debug_start(gen);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print starting points or table for algorithms into LOG file          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (GEN->table_on) {
    fprintf(LOG,"%s: use table (size = %d)\n",gen->genid,GEN->table_size);
    if (gen->debug & NINV_DEBUG_TABLE)
      for (i=0; i<GEN->table_size; i++)
	fprintf(LOG,"%s:\tx = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->table[i],GEN->f_table[i]);
  }
  else { /* no table */
    fprintf(LOG,"%s: starting points:\n",gen->genid);
    fprintf(LOG,"%s:\ts[0] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[0],GEN->CDFs[0]);
    if (! (gen->variant & NINV_VARFLAG_NEWTON))
      fprintf(LOG,"%s:\ts[1] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[1],GEN->CDFs[1]);
  }

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_start() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_sample( const struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (regula falsi)                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d iterations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN->max_iter);

} /* end of _unur_ninv_debug_sample_regula() */

/*---------------------------------------------------------------------------*/

void 
_unur_ninv_debug_chg_truncated( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(LOG,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);

} /* end of _unur_ninv_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
