/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_gw_init.c                                                *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Initializing routines for variant of Gilks and Wild.                 *
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

int
_unur_tdr_gw_starting_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              /* construction point, value of PDF at x */
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  
  /* compute paramters for all intervals */
  for( iv=GEN->iv; iv->next != NULL; ) {

    /* compute parameters for interval */
    switch (_unur_tdr_gw_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:   /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN->n_ivs);
      
      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make shure that we will never use this segment */
	iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      else
	/* we need a pointer to the previous entry in the list */
	iv->next->prev = iv;
      continue;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area below hat infinite.
       insert new construction point. */
    x = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */

    /* value of PDF at x */
    fx = PDF(x);

    /* add a new interval, but check if we had to used too many intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */


    /* if fx is 0, then we can cut off the tail of the distribution
       (since it must be T-concave)  */
    if (fx <= 0.) {
      if (iv->fx <= 0.) {
	/* cut off left tail */
	iv_new->next = iv->next;
	free(iv); 
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	iv_new->prev = NULL;
	/* compute the parameters for the new left part */
	iv = iv_new;
      }
      else if (iv->next->fx <= 0.) {
	/* cut off right tail */
	free(iv->next);
	--(GEN->n_ivs);	
	iv->next = iv_new;
	iv_new->prev = iv;
	/* compute the paramters for the new part */
	/* (nothing to do here) */
      }
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	free(iv_new);
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    else {
      /* insert new interval into linked list */
      iv_new->prev = iv;
      iv_new->next = iv->next;
      iv->next->prev = iv_new;
      iv->next = iv_new;
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_starting_intervals() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_gw_dars( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling  (Gilks&Wild)           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv, *iv_next;
  double Alimit;               /* threshhold value for splitting interval */
  double x0, x1;               /* boundary of interval */
  double xsp, fxsp;            /* splitting point in interval */
  double xAhatl, xAhatr, xAsqueeze;
  int rule;                    /* id for splitting rule that has been applied */
  int n_splitted;              /* count splitted intervals */
  int splitted;                /* result of splitting routine */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* now split intervals */
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_ivs < GEN->max_ivs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN->n_ivs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_ivs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 
    
    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (iv = GEN->iv; iv->next != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
      
      /* do not exceed the maximum number of intervals */
      if (GEN->n_ivs >= GEN->max_ivs)
	break;

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  /* goto next interval */

      /* store pointer to next interval */
      iv_next = iv->next;

      /* boundary of interval */
      x0 = iv->x;
      x1 = iv->next->x;

      /* get splitting point */
      for (rule = GEN->darsrule; rule <= 3; rule++) {
	switch (rule) {
	case 1:   /* rule 1: expected value */
	  
	  if ( _unur_FP_is_minus_infinity(x0) ||
	       _unur_FP_is_infinity(x1) ||
	       _unur_FP_approx(x0,x1) )
	    /* we do not use the expected value in case of unbounded intervals */
	    continue;  /* try next rule */
	  
	  /* l.h.s. hat */
	  xAhatl = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);
	  /* r.h.s. hat */
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv->next, iv->next->dTfx, iv->ip);
	  /* squeeze */
	  if (iv->Asqueeze > 0.)
	    /* always integrate from point with greater value of transformed density
	       to the other point */
	    xAsqueeze = (iv->Tfx > iv->next->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, iv->sq, x1)
	      : _unur_tdr_interval_xxarea( gen, iv->next, iv->sq, x0);
	  else  /* there is no squeeze */
	    xAsqueeze = 0.;
	  
	  /* check results */
	  if (! (_unur_isfinite(xAhatl) && _unur_isfinite(xAhatr) && _unur_isfinite(xAsqueeze))
	      || _unur_FP_equal(iv->Ahat,iv->Asqueeze) )
	    continue;  /* try next rule */
	  
	  /* compute expected value */
	  xsp = (xAhatl+xAhatr-xAsqueeze) / (iv->Ahat - iv->Asqueeze);
	  break;
	  
	case 2:   /* rule 2: arcmean */
	  xsp = _unur_arcmean(x0,x1);
	  break;
	  
	case 3:  /* rule 3: mean */
	  if (_unur_FP_is_minus_infinity(x0) || _unur_FP_is_infinity(x1))
	    /* no arithmetic mean for unbounded intervals */
	    continue;  /* try next rule */
	  xsp = 0.5 * (x0 + x1);
	  break;
	  
	default:   /* this should not happen */
	  continue;
	}
	
	/* value of PDF at splitting point */
	fxsp = PDF(xsp);
	
	/* now split interval at given point */
	splitted = _unur_tdr_gw_interval_split(gen, iv, xsp, fxsp);

	if (splitted==UNUR_SUCCESS || splitted==UNUR_ERR_INF) {
	  /* splitting successful */
	  if (splitted==UNUR_SUCCESS) ++n_splitted;
	  /* otherwise the area below the hat was not bounded */

	  /* now depending on the location of xps in the interval iv,
	     iv points to the left of the two new intervals,
	     or to the right of the two new intervals.
	     For the first case we have to move the pointer to the
	     new right interval. Then iv will be moved to the next
	     old interval in the list by the for loop.
	     We can distinguish between these two cases by looking 
	     at the iv->next pointer and compare it the pointer in the
	     old unsplitted interval.                                  */
	  if (iv->next != iv_next)
	    iv = iv->next;
	  /* no more splitting points in this interval */
	  break;
	}
	else if (splitted!=UNUR_ERR_SILENT) {
	  /* some serious error occurred */
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  return UNUR_ERR_GEN_CONDITION;
	}
	/* else: could not split construction points */
      }
    }
  
    if (n_splitted == 0) {
      /* we are not successful in splitting any inteval.
	 abort to avoid endless loop */
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }

  /* ratio between squeeze and hat o.k. ? */
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_ivs >= GEN->max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN->max_ivs = GEN->n_ivs;
  }
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_dars() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_gw_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (Gilks & Wild variant)                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... do not add this construction point             */
     /*   UNUR_ERR_INF    ... area = INFINITY                                */
     /*   others          ... error (PDF not T-concave)                      */
     /*----------------------------------------------------------------------*/
{
  double Ahatl;    /* area below hat at left side of intersection point */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents.
     used to partition interval into left hand part (construction point of tangent
     on the left hand boundary) and right hand part (construction point of tangent
     on the left hand boundary). */
  if ( _unur_tdr_tangent_intersection_point(gen,iv,&(iv->ip))!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;

  /* squeeze and area below squeeze */
  if (iv->Tfx > -INFINITY && iv->next->Tfx > -INFINITY) {

    /* we do not compute the slope when the construction points
       are too close. at least 8 significant digits should remain. */
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   /* construction points too close */

    /* slope of transformed squeeze */
    iv->sq = (iv->next->Tfx - iv->Tfx) / (iv->next->x - iv->x);

    /* check squeeze */
    /* we have to take care about round off error.
       the following accepts PDFs with might be a little bit not T_concave */
    if ( ( (iv->sq > iv->dTfx       && (!_unur_FP_approx(iv->sq,iv->dTfx)) ) || 
	   (iv->sq < iv->next->dTfx && (!_unur_FP_approx(iv->sq,iv->next->dTfx)) ) )
	 && iv->next->dTfx < INFINITY ) {
      /* There are big troubles when the density is extremely small. 
	 Then round-off errors may cancel out all significant figures and
	 0 remains. Thus we simply ignore all violations when the 
	 slope of the squeeze or tangent is 0.  */
      if ( !_unur_iszero(iv->sq) && !_unur_iszero(iv->dTfx) && !_unur_iszero(iv->next->dTfx) ) {
      	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      	return UNUR_ERR_GEN_CONDITION;
      }
    }

    /* volume below squeeze */
    /* always integrate from point with greater value of transformed density
       to the other point */
    iv->Asqueeze = (iv->Tfx > iv->next->Tfx)
      ? _unur_tdr_interval_area( gen, iv, iv->sq, iv->next->x)
      : _unur_tdr_interval_area( gen, iv->next, iv->sq, iv->x);
    
    /* check for fatal numerical errors */
    if (!_unur_isfinite(iv->Asqueeze))
      iv->Asqueeze = 0.;

  }
  else {  /* no squeeze */
    iv->sq = 0.;
    iv->Asqueeze = 0.;
  }

  /* volume below hat */
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, iv->ip);

  /* areas below head unbounded ? */
  if (! (_unur_isfinite(Ahatl) && _unur_isfinite(iv->Ahatr)) )
    return UNUR_ERR_INF;

  /* total area */
  iv->Ahat = iv->Ahatr + Ahatl;

  /* check area */
  /* we cannot be more accurate than in the `check squeeze' section */
  if ( iv->Asqueeze > iv->Ahat && !(_unur_FP_approx(iv->Asqueeze, iv->Ahat)) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"A(squeeze) > A(hat). PDF not T-concave!");
    return UNUR_ERR_GEN_CONDITION; 
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_interval_parameter() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_gw_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv_oldl, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* split interval iv_oldl into two intervals at point x                 */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv_oldl ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   fx      ... value of PDF at x                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv_newr;  /* pointer to new interval */
  struct unur_tdr_interval iv_bak;    /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_TDR_IV,UNUR_ERR_COOKIE);

  /* we cannot split point when x is not finite (NaN or Infinity) */
  if (!_unur_isfinite(x)) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not finite (skipped)");
    return UNUR_ERR_SILENT;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_gw_debug_split_start( gen,iv_oldl,x,fx );
#endif

  /* the splitting point must be inside the interval */
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN->n_ivs * (iv_oldl->Ahat - iv_oldl->Asqueeze) / (GEN->Atotal - GEN->Asqueeze))
       < GEN->bound_for_adding)
    return UNUR_ERR_SILENT;

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }

  /* back up data */
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_tdr_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {
    
    /* one of the two boundary points must be 0, too! */
    if (iv_oldl->fx <= 0.) {
      /* chop off left part (it's out of support) */
      iv_oldl->x = x;
    }
    else if (iv_oldl->next->fx <= 0.) {
      /* chop off right part (it's out of support) */
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* compute parameters for chopped interval */
    success = _unur_tdr_gw_interval_parameter(gen, iv_oldl);

    /* we did not add a new interval */
    iv_newr = NULL;
  }

  else {
    
    /* we need a new interval */
    iv_newr = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_newr == NULL) {
      /* PDF(x) < 0 or overflow !! */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* insert into linked list */
    iv_newr->prev = iv_oldl;
    iv_newr->next = iv_oldl->next;
    iv_oldl->next->prev = iv_newr;
    iv_oldl->next = iv_newr;
    
    /* compute parameters for interval */
    success   = _unur_tdr_gw_interval_parameter(gen, iv_oldl);
    success_r = _unur_tdr_gw_interval_parameter(gen, iv_newr);
    
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  
  /* successfull ? */
  if (success!=UNUR_SUCCESS) {
    /* cannot split interval at given point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");

    /* the case of unbounded hat is treated as round-off error for 
       very steep tangents. so we simply do not add this construction point. */

    /* restore old interval */
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_tdr_interval));
    /* remove from linked list; remaines to restore prev pointer in next interval */
    if (iv_oldl->next)
      iv_oldl->next->prev = iv_oldl;

    /* decrement counter for intervals and free unused interval */
    if (iv_newr) {
      --(GEN->n_ivs); 
      free( iv_newr );
    }

  return success;
  }

  /* successful */

  /* update total area below hat and squeeze */
  GEN->Atotal   = ( GEN->Atotal - iv_bak.Ahat
		   + iv_oldl->Ahat + ((iv_newr) ? iv_newr->Ahat : 0.) );
  GEN->Asqueeze = ( GEN->Asqueeze - iv_bak.Asqueeze
		   + iv_oldl->Asqueeze + ((iv_newr) ? iv_newr->Asqueeze : 0. ) );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_SPLIT)
    _unur_tdr_gw_debug_split_stop( gen,iv_oldl,iv_newr );
#endif

  /* when using inside Gibbs sampler Atotal might be 0 */
  if (GEN->Atotal <= 1.e10 * DBL_MIN) {
    _unur_error(gen->genid,UNUR_ERR_ROUNDOFF,"error below hat (almost) 0");
    return UNUR_ERR_ROUNDOFF;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_interval_split() */

/*---------------------------------------------------------------------------*/

