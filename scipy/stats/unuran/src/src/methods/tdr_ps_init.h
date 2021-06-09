/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_ps_init.c                                                *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Initializing routines for variant with proportional squeezes.        *
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
_unur_tdr_ps_starting_intervals( struct unur_gen *gen )
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
  double lb, flb;           /* left boundary point of domain of PDF  */
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* first interval in list */
  iv = GEN->iv;

  /* the left boundary of the domain:
     iv->x in the first interval is always the left boundary of 
     the domain of the PDF */
  lb = iv->x;
  flb = iv->fx;

  /* there is no use for a point iv->x that is not used as a 
     construction point in variants PS and IA.
     (In variant GW it is used to store the left boundary of the domain
     of the PDF)
     Thus we remove it from the list. At such points the slope of the
     tangents to the transformed density is set to INFINITY. */
  if (_unur_FP_is_infinity(iv->dTfx)) {
    GEN->iv = iv->next;
    GEN->iv->prev = NULL;
    free (iv);
    --(GEN->n_ivs);
    iv = GEN->iv;
  }

  /* set left boundary:
     it is stored in iv->ip in the first interval */
  iv->ip = lb;
  iv->fip = flb;

  /* compute paramters for all intervals */
  while (iv) {
    if (iv->next == NULL) {
      /* the last interval in the list*/

      /* analogously to variant GW we want to have a last virtual 
	 (stopping) interval, that simply stores the right boundary
	 point and guaranties that the loop in indexed search stops
	 on an existing interval.
	 However in the case where iv->x is used as construction point
	 we have to add such an interval. */
      if (!_unur_FP_is_infinity(iv->dTfx)) {
	/* get interval */
	iv->next = iv_new = _unur_tdr_interval_new( gen, iv->x, 0., FALSE );
	if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */
	/* link into list */
	iv_new->prev = iv;
	/* copy right boundary of domain */
	iv_new->ip = iv->x;
	iv_new->fip = iv->fx;
	/* make shure that we will never use this interval for sampling. */
	iv->next->Asqueeze = iv->next->Ahat = iv->next->Ahatr = 0.;
	iv->Acum = INFINITY;
	iv->next-> sq = 0.;
	/* we even have to to some additional work */
      }
      else
	/* nothing to do any more */
	break;
	
    }

    /* compute parameters for interval */
    switch (_unur_tdr_ps_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:     /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:    /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:    /* construction points too close */
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

    /* check if we had to used too many intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* area below hat infinite --> insert new construction point.
       We have to find out on which side of the construction point
       the area of the hat is unbounded. */
    if (iv->Ahatr >= INFINITY) {
      /* right hand side */

      /* iv should never be the last (virtual) interval */
      CHECK_NULL(iv->next,0);

      /* use mean point between the construction point and the right
	 boundary of the interval.
	 (The right boundary point might be a better choice
	 but cannot be used in every case.) */
      x = _unur_arcmean(iv->x,iv->next->ip);
      fx = PDF(x);

      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */

      /* if fx is 0, then we can cut off the tail of the distribution
	 (since it must be T-concave)  */
      if (fx <= 0.) {
	if (iv->next->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  free(iv_new);
	  return UNUR_ERR_GEN_CONDITION;
	}

	/* cut off right tail */
	free(iv->next);
	--(GEN->n_ivs);
	iv->next = iv_new;
	iv_new->prev = iv;

      }

      else {
	/* insert the new interval into the linked list after the old one. */
	iv_new->prev = iv;
	iv_new->next = iv->next;
	iv->next->prev = iv_new;
	iv->next = iv_new;
      }
      /* the old interval has to be recomputed. */

    }

    else {
      /* left hand side */

      x = _unur_arcmean(iv->ip,iv->x);
      fx = PDF(x);

      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */


      /* if fx is 0, then we can cut off the tail of the distribution
	 (since it must be T-concave)  */
      if (fx <= 0.) {
	if (iv->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  free(iv_new);
	  return UNUR_ERR_GEN_CONDITION;
	}

	/* cut off left tail */
	iv_new->next = iv->next;
	iv_new->prev = iv->prev;
	iv_new->ip = iv->ip;
	iv_new->fip = iv->fip;
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	/* continue with this new interval */
	free(iv);
	iv = iv_new;
      }

      else {

	if (iv->prev) {
	  /* insert new interval in just before the old unbounded one */
	  iv_tmp = iv->prev;
	  iv_new->prev = iv->prev;
	  iv_new->next = iv;
	  iv->prev->next = iv_new;
	  iv->prev = iv_new;
	  
	  /* make sure that _unur_arcmean(iv->ip,iv->x) is never out of range */
	  iv_new->ip = iv->ip;
	  
	  /* continue with the interval before the old one
	     (neccessary since it will change too). */
	  iv = iv_tmp;
	}
	else { /* iv->prev == NULL */
	  /* insert new interval as first entry in list */
	  iv_new->ip = iv->ip;
	  iv_new->fip = iv->fip;
	  iv_new->prev = NULL;
	  iv_new->next = iv;
	  iv->prev = iv_new;
	  GEN->iv = iv_new;
	  
	  /* continue with this new interval */
	  iv = iv_new;
	}
      }
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_starting_intervals() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_ps_dars( struct unur_gen *gen )
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
  double Adiff;                /* area between hat and squeeze */
  double Alimit;               /* threshhold value for splitting interval */
  double squeeze_gw;           /* slope of squeeze in variant GW */
  double Asqueeze_gw;          /* area below squeeze in variant GW */
  double Ahat_gw;              /* area below hat in variant GW */
  double x0, x1;               /* boundary of interval */
  double xsp, fxsp;            /* splitting point in interval */
  double xAhatl, xAhatr, xAsqueeze_gw;
  int rule;                    /* id for splitting rule that has been applied */
  int n_splitted;              /* count splitted intervals */
  int splitted;                /* result of splitting routine */
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* split intervals */
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
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);

      /* do not exceed the maximum number of intervals */
      if (GEN->n_ivs >= GEN->max_ivs)
	break;

      /* boundary of interval:
	 in opposition to GW the construction point is in the middle of
	 an interval. The boundary points are the intersection points of
	 the tangents. Notice that it does not make sense to have a
	 splitting point near the construction point in the middle of
	 the interval.
	 Thus we use the same intervals as in GW, i.e., the intervals
	 between the construction points. However, now we have to 
	 compute some data and the first and last interval has to 
	 be treated differently.
      */
      if (iv==GEN->iv) {
	/* the first interval */
	x0 = iv->ip;       /* left boundary of interval */
	x1 = iv->x;        /* right boundary of interval */
	Adiff = iv->Ahat - iv->Ahatr; /* area below hat */
	Adiff *= 1. - iv->sq;         /* area between hat and squeeze */
      }
      else {
	/* all the other intervals:
	   we use the intervals between the preceeding construction point
	   and the construction point in the given interval for
	   finding a splitting point.                                */
	x0 = iv->prev->x;       /* left boundary of interval */
	x1 = iv->x;             /* right boundary of interval */
	/* area between hat and squeeze.
	   (notice: we habe two parts in two intervals.)      */
	Adiff = ( (1. - iv->prev->sq) * iv->prev->Ahatr +
		  (1. - iv->sq) * (iv->Ahat - iv->Ahatr) );
      }

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if (Adiff <= Alimit) 
	continue;  /* goto next interval */

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
	  xAhatl = ( (iv->prev == NULL) ? 0. : 
		     _unur_tdr_interval_xxarea( gen, iv->prev, iv->prev->dTfx, iv->ip) );
	  /* r.h.s. hat */
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);

	  /* area below hat */
	  Ahat_gw = ( ((iv->prev == NULL) ? 0. : iv->prev->Ahatr)
		      + (iv->Ahat - iv->Ahatr) );

	  /* squeeze */
	  if (iv->Asqueeze > 0. && iv->prev != NULL) {
	    /* slope of transformed squeeze for variant GW */
	    squeeze_gw = (iv->Tfx - iv->prev->Tfx) / (iv->x - iv->prev->x);

	    /* always integrate from point with greater value of transformed density
	       to the other point */
	    xAsqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_xxarea( gen, iv->prev, squeeze_gw, iv->x);
	    Asqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_area( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_area( gen, iv->prev, squeeze_gw, iv->x);

	    /* check for fatal numerical errors */
	    if (!_unur_isfinite(Asqueeze_gw))  Asqueeze_gw = 0.;

	  }
	  else { /* there is no squeeze */
	    xAsqueeze_gw = 0.;
	    Asqueeze_gw = 0.;
	  }

	  /* check results */
	  if (! (_unur_isfinite(xAhatl) && _unur_isfinite(xAhatr) && _unur_isfinite(xAsqueeze_gw))
	      || _unur_FP_equal(Ahat_gw,Asqueeze_gw) )
	    continue;  /* try next rule */

	  /* compute expected value */
	  xsp = (xAhatl+xAhatr-xAsqueeze_gw) / (Ahat_gw - Asqueeze_gw);
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

	/* store pointer to next interval */
	iv_next = iv->next;

	/* now split interval at given point.
	   we have to take care that xsp might be either in
	   interval iv or in its preceeding interval, depending
	   whether xsp is less than iv->ip or not. */

	splitted = _unur_tdr_ps_interval_split(gen, ((xsp<iv->ip)?iv->prev:iv), xsp, fxsp);

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
	/* else: could not split construction points: too close (?) */
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

} /* end of _unur_tdr_ps_dars() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_ps_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (proportional squeezes)                      */
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
  double hxl, hxr; /* value of hat at left and right point of interval */
  double sq;       /* ration PDF(x) / hat(x) */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents. it is used as boundaries of intervals.
     it is stored together with the right hand (i.e. next) construction point. */
  if (_unur_tdr_tangent_intersection_point(gen,iv,&(iv->next->ip))!=UNUR_SUCCESS)
    return UNUR_ERR_GEN_CONDITION;
  /* value of PDF at intersection point */
  iv->next->fip = _unur_FP_is_infinity(iv->next->ip) ? 0. : PDF(iv->next->ip);

  /* volume below hat */
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->next->ip);

  /* areas below head unbounded ? */
  if (! (_unur_isfinite(Ahatl) && _unur_isfinite(iv->Ahatr)) )
    return UNUR_ERR_INF;

  /* total area */
  iv->Ahat = iv->Ahatr + Ahatl;

  /* compute squeeze:
     squeeze ration = min_{boundary points} PDF(x) / hat(x) */
  
  /* left boundary point */
  hxl = _unur_tdr_eval_intervalhat(gen,iv,iv->ip);
  if (_unur_FP_greater(iv->fip, hxl) ) {
    /* PDF(x) > hat(x); this should not happen */
    if ( (iv->fip < 1.e-50) || _unur_FP_approx(iv->fip, hxl)) {
      /* hat(x) and PDF(x) are approximatly the same, or extremely small.
	 assume round-off error */
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      /* we really have PDF(x) > hat(x); we do not assume a round-off error */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  iv->sq = (_unur_FP_is_infinity(hxl) || hxl <= 0.) ? 0. : iv->fip / hxl;

  /* right boundary point */
  hxr = _unur_tdr_eval_intervalhat(gen,iv,iv->next->ip);
  if (_unur_FP_greater(iv->next->fip, hxr)) {
    /* PDF(x) > hat(x); this should not happen */
    if ((iv->next->fip < 1.e-50) || _unur_FP_approx(iv->next->fip, hxr)) {
      /* hat(x) and PDF(x) are approximatly the same, or extremely small.
	 assume round-off error */
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      /* we really have PDF(x) > hat(x); we do not assume a round-off error */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  sq = (_unur_FP_is_infinity(hxr) || hxr <= 0.) ? 0. : iv->next->fip / hxr;

  /* squeeze */
  if (iv->sq > sq) iv->sq = sq;

  /* area below squeeze */
  iv->Asqueeze = iv->Ahat * iv->sq;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_interval_parameter() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_ps_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* split interval iv into two intervals at point x                      */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv      ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   fx      ... value of PDF at x                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *oldl, *oldr;  /* pointer to old intervals (left, right) */
  struct unur_tdr_interval *iv_new;       /* pointer to new interval */
  struct unur_tdr_interval oldl_bak, oldr_bak; /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL); COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);  COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);

  /* we cannot split point when x is not finite (NaN or Infinity) */
  if (!_unur_isfinite(x)) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not finite (skipped)");
    return UNUR_ERR_SILENT;
  }

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN->n_ivs * (iv->Ahat - iv->Asqueeze) / (GEN->Atotal - GEN->Asqueeze))
       < GEN->bound_for_adding)
    return UNUR_ERR_SILENT;

  /* the splitting point must be inside the interval */
  if (x < iv->ip || x > iv->next->ip) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }

  /* which side of construction point */
  if (x < iv->x) {
    /* left hand side */
    oldl = iv->prev;
    oldr = iv;
  }
  else {
    /* right hand side */
    oldl = iv;
    oldr = iv->next;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_start( gen,oldl,oldr,x,fx );
#endif

  /* back up data */
  if (oldl) memcpy(&oldl_bak, oldl, sizeof(struct unur_tdr_interval));
  memcpy(&oldr_bak, oldr, sizeof(struct unur_tdr_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {

    /* one of the two boundary points must be 0, too! */
    if (oldr->fip <= 0. && oldl==NULL) {
      /* chop off left part (it's out of support) */
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else if (oldr->fip <= 0. && oldr->next==NULL) {
      /* chop off right part (it's out of support) */
      oldr->x = x;
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* we do not add a new interval */
    iv_new = NULL;
  }

  else {

    /* we need a new interval */
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) {
      /* PDF(x) < 0 or overflow !! */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return -1;
    }

    /* insert into linked list */
    iv_new->prev = oldl;
    iv_new->next = oldr;
    oldr->prev = iv_new;
    if (oldl) oldl->next = iv_new;
    
  }

  /* compute parameters for intervals */
  success = UNUR_SUCCESS;
  /* left hand interval */
  if (oldl) {
   success = _unur_tdr_ps_interval_parameter(gen, oldl);
  }
  if( iv_new ) {
    /* middle (newly created) interval) */
    if (!oldl) {
      /* we have to copy the left intersection point from the 
         right hand interval */
      iv_new->ip = oldr->ip;
      iv_new->fip = oldr->fip;
    }
    success_r = _unur_tdr_ps_interval_parameter(gen, iv_new);
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if ( oldr->next ) {
    /* right hand interval */
    success_r = _unur_tdr_ps_interval_parameter(gen, oldr);
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
    if (oldl) memcpy(oldl, &oldl_bak, sizeof(struct unur_tdr_interval));
    memcpy(oldr, &oldr_bak, sizeof(struct unur_tdr_interval));
    /* remove from linked list; remaines to restore prev pointer in next interval */
    oldr->prev = oldl;
    if (oldl) oldl->next = oldr;

    /* decrement counter for intervals and free unused interval */
    if (iv_new) {
      --(GEN->n_ivs); 
      free( iv_new );
    }

  return success;
  }

  /* successful */

  /* we have update the pointer to the list */
  if (oldl == NULL && iv_new)
    /* new first entry */
    GEN->iv = iv_new;

  /* update total area below hat and squeeze */
  GEN->Atotal   = ( GEN->Atotal + (oldr->Ahat - oldr_bak.Ahat)
		   + ((oldl) ? (oldl->Ahat - oldl_bak.Ahat) : 0.)
		   + ((iv_new) ? iv_new->Ahat : 0.) );
  GEN->Asqueeze = ( GEN->Asqueeze + (oldr->Asqueeze - oldr_bak.Asqueeze)
		   + ((oldl) ? (oldl->Asqueeze - oldl_bak.Asqueeze) : 0.)
		   + ((iv_new) ? iv_new->Asqueeze : 0.) );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_stop( gen,oldl,iv_new,oldr );
#endif

  /* when using inside Gibbs sampler Atotal might be 0 */
  if (GEN->Atotal <= 1.e10 * DBL_MIN) {
    _unur_error(gen->genid,UNUR_ERR_ROUNDOFF,"error below hat (almost) 0");
    return UNUR_ERR_ROUNDOFF;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_interval_split() */

/*---------------------------------------------------------------------------*/
