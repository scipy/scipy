/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_info.ch                                                 *
 *                                                                           *
 *   Routines for creating info strings.                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_ninv_info( struct unur_gen *gen, int help )
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
  double n_iter;
  int samplesize = 10000;

  int use_newton = (gen->variant==NINV_VARFLAG_NEWTON) ? TRUE : FALSE;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = CDF");
  if (use_newton) 
    _unur_string_append(info," PDF");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n\n");
      
  /* method */
  _unur_string_append(info,"method: NINV (Numerical INVersion)\n");
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    _unur_string_append(info,"   Newton method\n");
    break;
  case NINV_VARFLAG_BISECT:
    _unur_string_append(info,"   Bisection method\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    _unur_string_append(info,"   Regula falsi\n");
    break;
  }
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  n_iter = unur_test_count_pdf(gen,samplesize,FALSE,NULL)/(2.*samplesize);
  if (!use_newton) n_iter *= 2.;
  _unur_string_append(info,"   average number of iterations = %.2f  [approx.]\n", n_iter);

  if (gen->set & NINV_SET_U_RESOLUTION) {
    if (DISTR.cdf) {
      double max_error=1.; double MAE=1.;
      unur_test_u_error(gen, &max_error, &MAE, 1.e-20, 1000, 
		     FALSE, FALSE, FALSE, NULL);
      _unur_string_append(info,"   u-error         <= %g  (mean = %g)  [rough estimate]\n", max_error, MAE);
    }
    else {
      _unur_string_append(info,"   u-error            NA  [requires CDF]\n");
    }
    _unur_string_append(info,  "     [ u-resolution = %g ]\n",GEN->u_resolution);
  }

  if (GEN->table_on) {
    _unur_string_append(info,"   starting points = table of size %d\n", GEN->table_size);
  }
  else {
    _unur_string_append(info,"   starting points = ");
    if (use_newton) {
      _unur_string_append(info,"%g (CDF = %g)  %s\n", GEN->s[0], GEN->CDFs[0],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
    else {
      _unur_string_append(info,"%g, %g  (CDF = %g, %g)   %s\n",
			  GEN->s[0],GEN->s[1], GEN->CDFs[0],GEN->CDFs[1],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
  }
  _unur_string_append(info,"\n");


  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    switch (gen->variant) {
    case NINV_VARFLAG_NEWTON:
      _unur_string_append(info,"   usenewton\n");
      break;
    case NINV_VARFLAG_BISECT:
      _unur_string_append(info,"   usebisect\n");
      break;
    case NINV_VARFLAG_REGULA: default:
      _unur_string_append(info,"   useregula  [default]\n");
      break;
    }

    _unur_string_append(info,"   u_resolution = %g  %s  %s\n", GEN->u_resolution,
			(GEN->u_resolution > 0.) ? "" : "[disabled]", 
			(gen->set & NINV_SET_U_RESOLUTION) ? "" : "[default]");

    _unur_string_append(info,"   x_resolution = %g  %s  %s\n", GEN->x_resolution,
			(GEN->x_resolution > 0.) ? "" : "[disabled]", 
			(gen->set & NINV_SET_X_RESOLUTION) ? "" : "[default]");

    _unur_string_append(info,"   max_iter = %d  %s\n", GEN->max_iter,
			(gen->set & NINV_SET_MAX_ITER) ? "" : "[default]");

    /* Not displayed:
       int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
       int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
    */
    _unur_string_append(info,"\n");
  }


  /* Hints */
  if (help) {
    if (! (gen->set & NINV_SET_X_RESOLUTION) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase accuracy by decreasing \"x_resolution\".");
    if (! (gen->set & NINV_SET_MAX_ITER) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase \"max_iter\" if you encounter problems with accuracy.");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tdr_info() */


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
