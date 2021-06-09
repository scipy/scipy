/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_init.c                                                   *
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

#include "tdr_gw_init.ch"
#include "tdr_ps_init.ch"

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_TDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tdr_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* free parameters */
  _unur_par_free(par);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_tdr_debug_init_start(gen);
#endif

  /* set-up the generator */
  if (_unur_tdr_make_gen( gen ) != UNUR_SUCCESS) {
    _unur_tdr_free(gen); return NULL;
  }

  /* is there any hat at all ? */
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdr_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_tdr_debug_init_finished(gen);
#endif

  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

  /* o.k. */
  return gen;

} /* end of _unur_tdr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_make_gen( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make generator object                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer (almost) empty generator object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* remark:                                                              */
     /*   on success: gen points to the ready-to-use generator object        */
     /*   on error:   gen should not be used at all                          */
     /*----------------------------------------------------------------------*/
{ 
  int i,k;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* get starting points */
  if (_unur_tdr_starting_cpoints(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;

  /* compute intervals for given starting points */
  if (_unur_tdr_starting_intervals(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;

  /* update maximal number of intervals */
  if (GEN->n_ivs > GEN->max_ivs) GEN->max_ivs = GEN->n_ivs;
  
  if (gen->variant & TDR_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TDR_DEBUG_DARS) {
      /* make initial guide table (only necessary for writing debug info) */
      _unur_tdr_make_guide_table(gen);
      /* write info into LOG file */
      _unur_tdr_debug_dars_start(gen);
    }
#endif

    for (i=0; i<3; i++) {
      /* we make several tries */

      /* run DARS */
      if (_unur_tdr_run_dars(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;
    
      /* make initial guide table */
      _unur_tdr_make_guide_table(gen);

      /* check if DARS was completed */
      if (GEN->n_ivs < GEN->max_ivs) {
	/* ran ARS instead */
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_tdr_debug_dars_finished(gen);
#endif
  }
  
  else { /* do not run DARS */
    /* make initial guide table */
    _unur_tdr_make_guide_table(gen);
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end _unur_tdr_make_gen() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_tdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* which transformation */
  if (_unur_iszero(PAR->c_T))
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_LOG;
  else if (_unur_FP_same(PAR->c_T, -0.5))
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_SQRT;
  else
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_POW;

  /** TODO: remove this **/
  if ((gen->variant & TDR_VARMASK_T) == TDR_VAR_T_POW) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"c != 0. and c != -0.5 not implemented!");
    _unur_generic_free(gen);
    return NULL;
  }

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_tdr_getSAMPLE(gen);
  gen->destroy = _unur_tdr_free;
  gen->clone = _unur_tdr_clone;
  gen->reinit = _unur_tdr_reinit;

  /* set all pointers to NULL */
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */
  GEN->c_T = PAR->c_T;                /* parameter for transformation          */
  GEN->darsfactor = PAR->darsfactor;  /* factor for derandomized ARS           */
  GEN->darsrule = PAR->darsrule;      /* rule for finding splitting points in DARS */

  /* bounds for adding construction points */
  GEN->max_ivs = _unur_max(2*PAR->n_starting_cpoints,PAR->max_ivs);  /* maximum number of intervals */
#ifdef UNUR_ENABLE_INFO
  GEN->max_ivs_info = PAR->max_ivs;   /* ... for info string */
#endif
  GEN->max_ratio = PAR->max_ratio;    /* bound for ratio  Atotal / Asqueeze    */
  GEN->bound_for_adding = PAR->bound_for_adding;

  /* get center */
  if ( (gen->distr->set & UNUR_DISTR_SET_CENTER) ||
       (gen->distr->set & UNUR_DISTR_SET_MODE) ) {
    GEN->center = unur_distr_cont_get_center(gen->distr);
    /* center must be in domain */
    GEN->center = _unur_max(GEN->center,DISTR.BD_LEFT);
    GEN->center = _unur_min(GEN->center,DISTR.BD_RIGHT);
    gen->set |= TDR_SET_CENTER;
  }
  else {
    GEN->center = 0.;
    /* we cannot use the center as construction point */
    gen->variant &= ~TDR_VARFLAG_USECENTER;
  }

  /*   mode known and in given domain ?? */
  if ( !(gen->distr->set & UNUR_DISTR_SET_MODE)
       || (DISTR.mode < DISTR.BD_LEFT)
       || (DISTR.mode > DISTR.BD_RIGHT))
    /* we cannot use the mode as construction point */
    gen->variant = gen->variant & (~TDR_VARFLAG_USEMODE);

  /* copy starting points */
  GEN->n_starting_cpoints = PAR->n_starting_cpoints;
  if (PAR->starting_cpoints) {
    GEN->starting_cpoints = _unur_xmalloc( PAR->n_starting_cpoints * sizeof(double) );
    memcpy( GEN->starting_cpoints, PAR->starting_cpoints, PAR->n_starting_cpoints * sizeof(double) );
  }
  else {
    GEN->starting_cpoints = NULL;
  }

  /* copy percentiles */
  GEN->percentiles = NULL;
  if (gen->set & TDR_SET_N_PERCENTILES)
    unur_tdr_chg_reinit_percentiles( gen, PAR->n_percentiles, PAR->percentiles );

  /* copy all other parameters */
  GEN->retry_ncpoints = PAR->retry_ncpoints;   /* number of cpoints for second trial of reinit */

  /* set (default) boundaries for U */
  GEN->Umin = 0.;
  GEN->Umax = 1.;

  /* set default for DARS */
  if (!(gen->set & TDR_SET_USE_DARS) && !PAR->starting_cpoints)
    /* no starting points given by user
       --> enable derandomized ARS      */
    gen->variant |= TDR_VARFLAG_USEDARS;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_tdr_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_tdr_create() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_reinit( struct unur_gen *gen )
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
  struct unur_tdr_interval *iv,*next;
  double *bak_cpoints;
  int bak_n_cpoints;
  int i;
  int n_trials;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );

  /* first trial */
  n_trials = 1;

  /* which construction points should be used ? */
  if (gen->set & TDR_SET_N_PERCENTILES) {
    if (GEN->starting_cpoints==NULL || (GEN->n_starting_cpoints != GEN->n_percentiles)) {
      GEN->n_starting_cpoints = GEN->n_percentiles;
      GEN->starting_cpoints = _unur_xrealloc( GEN->starting_cpoints, GEN->n_percentiles * sizeof(double));
    }
    for (i=0; i<GEN->n_percentiles; i++) {
      GEN->starting_cpoints[i] = unur_tdr_eval_invcdfhat( gen, GEN->percentiles[i], NULL, NULL, NULL );
      if (!_unur_isfinite(GEN->starting_cpoints[i])) 
	/* we cannot use these starting points --> skip to second trial immediately */
	n_trials = 2;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_REINIT)
    _unur_tdr_debug_reinit_start(gen);
#endif

  /* make backup of cpoints */
  bak_n_cpoints = GEN->n_starting_cpoints;
  bak_cpoints = GEN->starting_cpoints;

  for (;; ++n_trials) {
    /* free linked list of intervals */
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
    GEN->iv = NULL;
    GEN->n_ivs = 0;
    GEN->Atotal = 0.;
    GEN->Asqueeze = 0.;

    if (n_trials > 2) {
      /* we have done our best */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points for reinit");
      GEN->n_starting_cpoints = bak_n_cpoints;
      GEN->starting_cpoints = bak_cpoints;
      return UNUR_FAILURE;
    }

    if (n_trials > 1) {
      /* second trial */
      GEN->n_starting_cpoints = GEN->retry_ncpoints;
      GEN->starting_cpoints = NULL;
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (gen->debug & TDR_DEBUG_REINIT)
	_unur_tdr_debug_reinit_retry(gen);
#endif
    }

    /* set-up the generator */
    if (_unur_tdr_make_gen( gen ) != UNUR_SUCCESS)
      continue;

    /* is there any hat at all ? */
    if (GEN->Atotal <= 0.)
      continue;

    /* reinit successful */
    break;

  }

  /* clean up */
  if (n_trials > 1) {
    GEN->n_starting_cpoints = bak_n_cpoints;
    GEN->starting_cpoints = bak_cpoints;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_REINIT)
    if (gen->debug) _unur_tdr_debug_reinit_finished(gen);
#endif

  /* (re)set sampling routine */
  SAMPLE = _unur_tdr_getSAMPLE(gen);

  return UNUR_SUCCESS;
} /* end of _unur_tdr_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_tdr_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_tdr_interval *iv,*next, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tdr_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
      clone_iv->prev = NULL;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
      clone_iv->prev = clone_prev;
    }
    /* next step */
    next = iv->next;
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;

  /* copy starting points */
  if (GEN->starting_cpoints) {
    CLONE->starting_cpoints = _unur_xmalloc( GEN->n_starting_cpoints * sizeof(double) );
    memcpy( CLONE->starting_cpoints, GEN->starting_cpoints, GEN->n_starting_cpoints * sizeof(double) );
  }

  /* copy percentiles */
  if (GEN->percentiles) {
    CLONE->percentiles = _unur_xmalloc( GEN->n_percentiles * sizeof(double) );
    memcpy( CLONE->percentiles, GEN->percentiles, GEN->n_percentiles * sizeof(double) );
  }

  /* make new guide table */
  CLONE->guide = NULL;
  _unur_tdr_make_guide_table(clone);

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tdr_clone() */

/*****************************************************************************/

void
_unur_tdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_tdr_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_tdr_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free list of starting points */
  if (GEN->starting_cpoints) 
    free (GEN->starting_cpoints);

  /* free list of percentiles */
  if (GEN->percentiles) 
    free (GEN->percentiles);

  /* free table */
  if (GEN->guide)  free(GEN->guide);

  /* free other memory not stored in list */
  _unur_generic_free(gen);

} /* end of _unur_tdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_tdr_starting_cpoints( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting intervals.                  */
     /* if not provided as arguments compute these                           */
     /* by means of the "equiangular rule" from AROU.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, fx, fx_last;
  int use_center, use_mode, is_mode, was_mode;
  int i, is_increasing;
  double extra_cpoint;
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  
  /* use mode as construction point ? */
  use_mode = (gen->variant & TDR_VARFLAG_USEMODE) ? TRUE : FALSE;

  /* use center as construction point ? */
  use_center = (!use_mode && (gen->variant & TDR_VARFLAG_USECENTER)) ? TRUE : FALSE;

  /* add extra construction point        */
  /* (use either mode or center or none) */
  extra_cpoint = use_mode ? DISTR.mode : (use_center ? GEN->center : 0. );

  /* reset counter of intervals */
  GEN->n_ivs = 0;

  /* prepare for computing construction points */
  if (!GEN->starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT  - GEN->center);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT - GEN->center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (GEN->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = DISTR.BD_LEFT;
  if (use_mode && DISTR.mode <= x) {
    /* this is the mode of the distribution */
    is_mode = TRUE;
    use_mode = FALSE;  /* do not use the mode again */
    is_increasing = FALSE;
  }
  else if (use_center && GEN->center <= x) {
    is_mode = FALSE;
    use_center = FALSE;     /* do not use the center again */
    is_increasing = TRUE;   /* the center may be left of (unknown) mode */
  }
  else {
    is_mode = FALSE;
    is_increasing = TRUE;
  }
    
  fx = fx_last = _unur_FP_is_minus_infinity(x) ? 0. : PDF(x);
  iv = GEN->iv = _unur_tdr_interval_new( gen, x, fx, is_mode );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */

  /* terminate beginning of list */
  iv->prev = NULL;

  /* now all the other points */
  for( i=0; i<=GEN->n_starting_cpoints; i++ ) {
    was_mode = is_mode;

    /* construction point */
    if (i < GEN->n_starting_cpoints) {
      if (GEN->starting_cpoints) {   
	/* construction points provided by user */
	x = GEN->starting_cpoints[i];
	/* check starting point */
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle ) + GEN->center;
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store 
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /* insert mode or center ? */
    if ((use_mode || use_center) && x >= extra_cpoint) {
      is_mode = use_mode;              /* the next construction point is the mode */
      use_center = use_mode = FALSE;   /* we use the mode only once (of course) */
      if (x>extra_cpoint) {
	x = extra_cpoint;     /* use the mode now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!GEN->starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x == extra_cpoint --> nothing to do */
    }
    else
      is_mode = FALSE;

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of PDF at starting point */
    fx = _unur_FP_is_infinity(x) ? 0. : PDF(x);

    /* check value of PDF at starting point */
    if (!is_increasing && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (is_mode && fx < fx_last * (1.-DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"mode -> ignore");
      continue;
    }
    if (was_mode && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"mode");
      return UNUR_ERR_GEN_DATA;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such point */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<GEN->n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the PDF is constant 0 on all construction points.
	     then we need both boundary points. */
	  iv->x = x;  /* we only have to change x, everything else remains unchanged */
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with PDF(x) > 0 */
	break;
    }
    
    /* need a new interval */
    iv->next = _unur_tdr_interval_new( gen, x, fx, is_mode );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */
    
    /* link into list and skip pointer to current interval */
    iv->next->prev = iv;
    iv = iv->next;

    /* PDF still increasing ? */
    if (is_increasing && fx < fx_last)
      is_increasing = 0;

    /* store last computed values */
    fx_last = fx;

  }

  /* we have left the loop with the right boundary of the support of PDF
     make shure that we will never use iv for sampling. */
  iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
  iv->Acum = INFINITY;
  iv->ip = iv->x;
  iv->fip = iv->fx;
  iv->next = NULL;         /* terminate list */
  --(GEN->n_ivs);           /* we do not count this interval */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_starting_cpoints() */

/*****************************************************************************/

int
_unur_tdr_starting_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    return _unur_tdr_gw_starting_intervals(gen);
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    return _unur_tdr_ps_starting_intervals(gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
} /* end of _unur_tdr_starting_intervals() */

/*****************************************************************************/

int
_unur_tdr_run_dars( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double Atot, Asqueezetot;    /* total area below hat and squeeze, resp. */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  
  /* there is no need to run DARS when the DARS factor is INFINITY */
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_tdr_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }

  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;

  /* now run DARS for different variants */
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    return _unur_tdr_gw_dars(gen);
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    return _unur_tdr_ps_dars(gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

} /* end of _unur_tdr_run_dars() */

/*****************************************************************************/

struct unur_tdr_interval *
_unur_tdr_interval_new( struct unur_gen *gen, double x, double fx, int is_mode )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x       ... left point of new interval                             */
     /*   fx      ... value of PDF at x                                      */
     /*   is_mode ... if TRUE, x is a mode of the PDF                        */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double dfx;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return NULL;
  }
  if (_unur_FP_is_infinity(fx)) {
    /* over flow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN->n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_TDR_IV);

  /* avoid uninitialized variables */
  iv->Acum = iv->Ahat = iv->Ahatr = iv->Asqueeze = 0.;
  iv->ip = iv->fip = iv->sq = 0.;

  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->fx = fx;            /* value of PDF at x */
  
  if (fx<=0.) {           /* --> -INFINITY */
    iv->Tfx = -INFINITY;  /* transformed density */
    iv->dTfx = INFINITY;  /* derivative of transformed density */
    return iv;
  }
  
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    /* transformed density */
    iv->Tfx = log(fx);
    /* derivative of transformed density */
    if (is_mode) { 
      /* we can set dPDF(x) = 0. for the mode */
      iv->dTfx = 0.; break; 
    }
    if (_unur_cont_have_dlogPDF(gen->distr)) {
      iv->dTfx = dlogPDF(x); break; 
    }
    else {
      dfx = dPDF(x);
      if (_unur_iszero(dfx))
	iv->dTfx = 0.;
      else
	iv->dTfx = (1./fx * dfx);   /* possible overflow ? */
    }
    break;

  case TDR_VAR_T_SQRT:
    /* transformed density */
    iv->Tfx = -1./sqrt(fx);
    /* derivative of transformed density */
    if (is_mode) { 
      /* we can set dPDF(x) = 0. for the mode */
      iv->dTfx = 0.; break; }
    if (_unur_cont_have_dlogPDF(gen->distr)) {
      iv->dTfx = -0.5 * iv->Tfx * dlogPDF(x);
      break;
    }
    else {
      dfx = dPDF(x);
      if (_unur_iszero(dfx))
	iv->dTfx = 0.;
      else
	iv->dTfx = (dfx<0.) ? -exp( -M_LN2 - 1.5*log(fx) + log(-dfx))
	  : exp( -M_LN2 - 1.5*log(fx) + log(dfx));
    }
    break;

  case TDR_VAR_T_POW:
    /** TODO **/
    /*      iv->Tfx = -pow(fx,GEN->c_T); */
    /*      iv->dTfx = 0.; */
    break;
  }

  /* the program requires dTfx > -INFINITY */
  if ( !(iv->dTfx > -INFINITY))
    iv->dTfx = INFINITY;

  return iv;

} /* end of _unur_tdr_interval_new() */

/*****************************************************************************/

int
_unur_tdr_tangent_intersection_point( struct unur_gen *gen, struct unur_tdr_interval *iv, double *ipt )
     /*----------------------------------------------------------------------*/
     /* compute cutting point of interval into left and right part.          */
     /* (1) use intersection point of tangents of transformed hat.           */
     /* (2) use mean point if (1) is unstable due to roundoff errors.        */
     /* (3) use boundary point which is closer to the mode. this is          */
     /*     important when the transformed tagents are extremely steep.      */
     /*     (This might cause a serious roundoff error while sampling.)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   ipt ... pointer to intersection point                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* 
     case: there is no tangent at one of the boundary points of the interval
           (then the slope is INFINITY)
     or
     case: the tangents are too steep  (--> case (3))
  */
  if ( iv->dTfx > 1.e+140 ) {
    *ipt = iv->x;        /* intersection point = left boundary of interval */
    return UNUR_SUCCESS; 
  }
  if ( iv->next->dTfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dTfx)) {
    *ipt = iv->next->x;   /* intersection point = right boundary of interval */
    return UNUR_SUCCESS; 
  }
  /** TODO: 1.e+140 (= sqrt(DBL_MAX) / 1.e15) is arbitrary  **/

  /* test for T-concavity */
  if ( _unur_FP_less( iv->dTfx, iv->next->dTfx ) ) {

    /* it might happen because of round-off errors 
       that iv->next->dTfx is almost zero although it should be large.
       thus we ignore this case. */
    if ( fabs(iv->dTfx) < DBL_EPSILON * fabs(iv->next->dTfx) ) {
      *ipt = iv->x;        /* intersection point = left boundary of interval */
      iv->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else if ( fabs(iv->next->dTfx) < DBL_EPSILON * fabs(iv->dTfx) ) {
      *ipt = iv->next->x;   /* intersection point = right boundary of interval */
      iv->next->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else {
/*        fprintf(stdout,"\ndTfx0 = %g < %g = dTfx1 (x0 = %g, x1 = %g)\n", */
/*  	      iv->dTfx,iv->next->dTfx,iv->x,iv->next->x); */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /** TODO: the following test is too sensitve to roundoff errors **/
  /*    if (iv->next->Tfx > iv->x + iv->dTfx*(iv->next->x - iv->x)) { */
  /*      _unur_warning(gen->genid,UNUR_ERR_INIT,"tangent below PDF  not T-concave!"); */
  /*      return UNUR_ERR_INIT; */
  /*    } */
  
  /* case (2): computing intersection of tangents is unstable */
  if (_unur_FP_approx(iv->dTfx, iv->next->dTfx)) {
    /* use mean point */
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }

  /* case (1): compute intersection point of tangents (regular case) */
  *ipt = ( (iv->next->Tfx - iv->Tfx - iv->next->dTfx * iv->next->x + iv->dTfx * iv->x) / 
	   (iv->dTfx - iv->next->dTfx) );

  /* check position of intersection point */
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    /* intersection point of tangents not in interval.
       This is mostly the case for numerical reasons.
       Thus we is the center of the interval instead.
       if the PDF not T-concave, it will catched at a later
       point when we compare slope of tangents and squeeze. */
    *ipt = 0.5 * (iv->x + iv->next->x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_tangent_intersection_point() */

/*****************************************************************************/

double
_unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute area below piece of hat or squeeze in                             */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   area                                                                    */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /* log(x)                                                                    */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(t-x0)) dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /* -1/sqrt(x)                                                                */
     /*   area = | \int_{x0}^x 1/(Tf(x0) + slope*(t-x0))^2 dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = infinity                                      if T(f(x)) >= 0    */
     /*        = | (x-x0) / (Tf(x0)*(Tf(x0)+slope*(x-x0))) |   otherwise          */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double area = 0.;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* if the construction point is at infinity, we cannot compute an area.
     (in this case we should have x == iv->x == INFINITY). */
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return 0.;

  /* unbounded? */
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:   /* T(x) = log(x) */

    if (!_unur_iszero(slope)) {                         
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = iv->fx / slope;
      else {
	double t = slope * (x - iv->x);
	if (fabs(t) > 1.e-6) {
	  if (t > MAXLOG / 10.) {
	    double xdiff = (x>iv->x) ? x - iv->x : iv->x - x;
	    area = exp( log(iv->fx) + log(xdiff) + t - log(t) );
	  }
	  else {
	    area = iv->fx * (x - iv->x) * ( exp(t) - 1. ) / t;
	  }
	}
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  area = iv->fx * (x - iv->x) * (1. + t/2. + t*t/6.);
	else
	  area = iv->fx * (x - iv->x) * (1. + t/2.);
      }
    }
    else { /* hat/squeeze almost constant */
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      else
	area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_SQRT:   /* T(x) = -1./sqrt(x) */

    if (!_unur_iszero(slope)) {
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = 1. / ( iv->Tfx * slope );
      else {
	/* compute value of transformed hat at integration boundary */
	double hx = iv->Tfx + slope * (x - iv->x);
	/* the transformed hat must always be below the x-axis.
	   otherwise the area below the hat in unbounded. */
	if (hx>=0.)
	  return INFINITY; 
	else
	  area = (x - iv->x) / ( iv->Tfx * hx );
      }
    }
    else { /* hat/squeeze almost constant */
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      else
	area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_POW:   /* T(x) = -x^c */
    /** TODO **/
    break;
  }

  return ( (area<0.) ? -area : area );

} /* end of _unur_tdr_interval_area() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_interval_xxarea( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute the interal of x times hat or squeeze ("expected value") in       */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   "expected value"                                                        */
     /*   (to get the real expected value, it must be divided by the area below   */
     /*   the function (hat or squeeze).)                                         */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /* log(x)                                                                    */
     /*   ev = \int_{x0}^x t * \exp(Tf(x0) + slope*(t-x0)) dt                     */
     /*      = 0.5 * f(x0) * (x^2 - x0^2)                            if slope = 0 */
     /*      = f(x0)/slope^2 * (\exp(slope*(x-x0))*(slope*x-1) - (slope*x0-1))    */
     /*                                                             if slope != 0 */
     /*                                                                           */
     /* -1/sqrt(x)                                                                */
     /*   ev = \int_{x0}^x t / (Tf(x0) + slope*(t-x0))^2 dt                       */
     /*      = 0.5 * f(x0) * (x^2 - x0^2)                            if slope = 0 */
     /*      = infinity                         if T(f(x)) >= 0 or |x| = infinity */
     /*      = x0 / (slope*Tf(x0)) - x / (slope*u) + log(u/Tf(x0)) / slope^2      */
     /*        where u = Tf(x0) + slope*(x-x0)                          otherwise */
     /*                                                                           */
     /* To get the right sign we have to return sign(x-x0) * ev.                  */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double ev = 0.;
  double hx,u;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* if the construction point is at infinity, we cannot compute the integral
     (in this case we should have x == iv->x == INFINITY). */
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return 0.;

  /* unbounded? */
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;


  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:    /* T(x) = log(x) */
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x)) {
      ev = iv->fx / (slope*slope) * (1-slope*iv->x);
    }
    else {
      u = (x-iv->x) * slope;

      if (fabs(u) > 1.e-6) {
	ev = iv->fx / (slope*slope) * (exp(u)*(slope*x-1.) - slope*iv->x + 1.);
      }
      else {
	/* use Taylor series */
	/* constant term */
	ev = 0.5 * (x+iv->x);
	if (fabs(u) > 0) {
	  /* 1st order expansion */
	  ev += 1./6. * (2.*x+iv->x) * u;
	  /* 2nd order expansion */
	  ev += 1./24. * (3.*x+iv->x) * u * u;
	}
	ev *= iv->fx * (x-iv->x);
      }
    }
    break;

  case TDR_VAR_T_SQRT:    /* T(x) = -1./sqrt(x) */
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
      /* the integral becomes INFINITY */
      return INFINITY;

    /* compute value of transformed hat at integration boundary */
    hx = iv->Tfx + slope * (x - iv->x);

    if (hx >= 0.)
      /* the transformed hat must always be below the x-axis.
	 otherwise the area below the hat in unbounded. */
      return INFINITY; 

    u = (x-iv->x) * slope / iv->Tfx;

    if (fabs(u) > 1.e-6) {
      ev = ( iv->x / (slope * iv->Tfx) - x / (slope * hx)
	     + log( hx / iv->Tfx ) / (slope*slope) );
    }
    else {
      /* use Taylor series */
      /* constant term */
      ev = 0.5 * (x+iv->x);
      if (fabs(u) > 0) {
	/* 1st order expansion */
	ev -= 1./3. * (2.*x+iv->x) * u;
	/* 2nd order expansion */
	ev += 1./4. * (3.*x+iv->x) * u * u;
      }
      ev *= iv->fx * (x-iv->x);
    }
    break;

  case TDR_VAR_T_POW:    /* T(x) = -x^c */
    /** TODO **/
    break;
  }

  return ((x>iv->x) ? ev : -ev);

} /* end of _unur_tdr_interval_xxarea() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_eval_intervalhat( struct unur_gen *gen, struct unur_tdr_interval *iv, double x )
     /*----------------------------------------------------------------------*/
     /* evaluate hat at x in interval.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval that stores constr. point of tangent   */
     /*   x   ... point at which hat(x) has to be computed                   */
     /*                                                                      */
     /* return:                                                              */
     /*   hat(x) or                                                          */
     /*   0. if x is not finite or                                           */
     /*   INFINITY the hat cannot be computed                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /* comment:                                                             */
     /*   x0    ... construction point of tangent (= iv->x)                  */
     /*                                                                      */
     /* log(x)                                                               */
     /*   hat(x) = f(x0) * exp( Tf'(x0)(x - x_0) )                           */
     /*                                                                      */
     /* 1/sqrt(x)                                                            */
     /*   hat(x) = 1/(Tf(x0) + Tf'(x0)(x - x_0))^2                           */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* we cannot compute the hat at x if any of the parameters are not finite  */
  if ( _unur_FP_is_minus_infinity(iv->Tfx) || _unur_FP_is_infinity(iv->dTfx) )
    return INFINITY;

  /* at +/- infinity the hat should be 0 (or infinity) */
  if ( _unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x) ||
       _unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x) )
    return 0.;

  /* now evaluate hat at x */
  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:
    /* T(x) = log(x) */
    return (iv->fx * exp( iv->dTfx * (x - iv->x) ));

  case TDR_VAR_T_SQRT:
    /* T(x) = -1./sqrt(x) */
    {
      /* compute value of transformed hat at x */
      double hx = iv->Tfx + iv->dTfx * (x - iv->x);
      /* hx must be less than 0 ! */
      return ((hx<0.) ? 1./(hx*hx) : INFINITY);
    }

  case TDR_VAR_T_POW:
    /* T(x) = -1./x^c */
    /** TODO **/
    return INFINITY;

  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }

} /* end of _unur_tdr_eval_intervalhat() */

/*****************************************************************************/

int
_unur_tdr_make_guide_table( struct unur_gen *gen )
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
  struct unur_tdr_interval *iv;
  double Acum, Asqueezecum, Astep;
  int max_guide_size;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    max_guide_size = (GEN->guide_factor > 0.) ? ((int)(GEN->max_ivs * GEN->guide_factor)) : 1;
    if (max_guide_size <= 0) max_guide_size = 1;   /* protect against overflow */
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tdr_interval*) );
  }

  /* first we need cumulated areas in intervals */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN->guide_size = (int)(GEN->n_ivs * GEN->guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   /* this is the last virtual intervall --> do not use */
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;

  return UNUR_SUCCESS;
} /* end of _unur_tdr_make_guide_table() */

/*****************************************************************************/

