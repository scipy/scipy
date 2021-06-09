/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_regula.ch                                               *
 *                                                                           *
 *   Routines for regula falsi root finding algorithm. c                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *

 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  For a given U find X with F(X) - U = 0.                                  *
 *                                                                           *
 *  Method: stabilized regula falsi preserving sign change;                  *
 *          bisection if slow convergence is detected.                       *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Regula falsi                                                           **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* maximum number of steps to find sign change in Regula Falsi               */
#define MAX_STEPS (100)

/* STEPFAC* (s[1]-s[0]) is used as first step length for finding sign change */
#define STEPFAC  (0.4)

/* for #steps > I_CHANGE_TO_BISEC Regula Falsi is always replaced by bisection*/
#define I_CHANGE_TO_BISEC (50)

/*---------------------------------------------------------------------------*/

double 
_unur_ninv_regula( const struct unur_gen *gen, double u )
     /*---------------------------------------------------------------------*/
     /*   algorithm: regula falsi with bisection steps                      */
     /*                                                                     */
     /*   parameters:                                                       */
     /*      gen ... pointer to generator object                            */
     /*      u   ... random number (uniform distribution)                   */
     /*                                                                     */
     /*   return:                                                           */
     /*     double (sample from random variate)                             */
     /*                                                                     */
     /*   error:                                                            */
     /*     return INFINITY                                                 */
     /*                                                                     */
     /*   Remark:                                                           */
     /*     The routine computes the root of CDF(x)-u                       */
     /*                                                                     */
     /*   Reference:                                                        */
     /*   [1] Neumaier A. (to be published):                                */
     /*       Introduction to numerical analysis,                           */
     /*       Cambridge University Press                                    */
     /*---------------------------------------------------------------------*/
{ 
  double x1, x2, a, xtmp;/* points for regular falsi                        */
  double f1, f2,fa, ftmp;/* function values at x1, x2, xtmp                 */
  double length;         /* oriented length of the interval with sign change*/
  double lengthabs;      /* absolute length of interval                     */
  double lengthsgn;      /* orientation of the Intervalls                   */
  double dx;             /* RF-stepsize                                     */
  int count_nosc = 0;    /* counter for  "no sign change occured"           */
  int i;                 /* loop variable, index                            */
  double min_step_size;  /* minimal step size for regula falsi              */
  double rel_u_resolution; /* relative u resolution                         */

  /* check arguments */
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);

  /* compute relative u resolution */
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
		       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
		       INFINITY );

  /* -- 1. + 2. find bracket -- */

  if ( _unur_ninv_bracket( gen, u, &x1, &f1, &x2, &f2 ) 
       != UNUR_SUCCESS )
    return x2;

  /* sign changes always within [a, x2] */
  a = x1; fa = f1; 

  /* -- 3. secant step, preserve sign change -- */

  /* secant step, preserve sign change */
  for (i=0; TRUE; i++) {

    /* -- 3a. check sign of f1 and f2 -- */

    if ( f1*f2 < 0.) { 
      /* sign change found     */
      count_nosc = 0;   /* reset counter */
      /* f2 always less (better), otherwise exchange */
      if ( fabs(f1) < fabs(f2) ) {
	xtmp = x1; ftmp = f1;
	x1 = x2;   f1 = f2;
	x2 = xtmp; f2 = ftmp;
      }
      /* sign changes within [a, x2] */
      a = x1; fa = f1;
    }
    else {
      /* increment counter for "no sign change occured" */
      count_nosc++;
      /* must not update 'a' */
    }

    /* length of enclosing bracket */
    length = x2 - a;                       /* oriented length */
    lengthabs = fabs(length);              /* absolute length */
    lengthsgn = (length < 0.) ? -1. : 1.;

    /* -- 3b. check stopping criterion -- */

    /* check accuracy goal */
    if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
    			     x2, f2, a, fa ))
      break;

    if (i >= GEN->max_iter)
      /* maximum number of iterations reached --> abort */
      break;

    /* -- 3c. try secant step -- */

    /* step size: secant or bisection step */
    dx = (_unur_FP_same(f1,f2)) ? length/2. : f2*(x2-x1)/(f2-f1) ;  

    /* minimal step size */
    if (GEN->u_resolution < 0.) 
      /* we only look at the x-error, so we do not need shorter steps */
      min_step_size = fabs(x2) * GEN->x_resolution;
    else
      min_step_size = lengthabs * DBL_EPSILON; 

    if ( fabs(dx) < min_step_size ) {
      dx = lengthsgn * 0.99 * min_step_size;

      while ( x2 == x2 - dx ){ /* dx too small  */
        if ( dx != 2.*dx)    /* near limit of calculations */
          dx = 2.*dx;
        else
          dx = length/2.;    /* bisection step   */
      }
    }

    /* bisection step if:                             */  
    /* no sign change   || step leads out of interval */
    if ( count_nosc > 1 || i > I_CHANGE_TO_BISEC ||
	 (lengthabs-GEN->x_resolution*fabs(x2))/(dx*lengthsgn) <= 1. )
      dx = length/2.; /* bisection step        */
  
    /* -- 3c. update point -- */    
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
    
  }  /* end of for-loop */

  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
 
  /* ensure location within given (truncated) domain */
  x2 = _unur_max( x2, DISTR.trunc[0]);
  x2 = _unur_min( x2, DISTR.trunc[1]);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file (in case error) */
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample( gen,u,x2,f2,i );
#endif

  /* finished */
  return x2;

} /* end of _unur_ninv_regula() */

/*****************************************************************************/
/**  Bisection method                                                       **/
/*****************************************************************************/

double 
_unur_ninv_bisect( const struct unur_gen *gen, double u )
     /*---------------------------------------------------------------------*/
     /*   algorithm: bisection method                                       */
     /*                                                                     */
     /*   parameters:                                                       */
     /*      gen ... pointer to generator object                            */
     /*      u   ... random number (uniform distribution)                   */
     /*                                                                     */
     /*   return:                                                           */
     /*     double (sample from random variate)                             */
     /*                                                                     */
     /*   error:                                                            */
     /*     return INFINITY                                                 */
     /*                                                                     */
     /*   Remark:                                                           */
     /*     The routine computes the root of CDF(x)-u                       */
     /*                                                                     */
     /*   References:                                                       */
     /*   [2] Dahlquist, G. and Bj{\"o}rck, {\AA} (2008):                   */
     /*       Numerical methods in scientific computing, Vol 1.,            */
     /*       SIAM, Philadelphia, PA.                                       */
     /*---------------------------------------------------------------------*/
{ 
  double x1, x2, mid=0.;   /* boundary points and mid point of bracket      */
  double f1, f2, fmid;     /* function values at x1, x2, and mid            */
  int i;                   /* loop variable, number of iterations           */
  double rel_u_resolution; /* relative u resolution                         */

  /* check arguments */
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);

  /* compute relative u resolution */
  rel_u_resolution = ( (GEN->u_resolution > 0.) ?
		       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
		       INFINITY );

  /* -- 1. + 2. find bracket -- */

  if ( _unur_ninv_bracket( gen, u, &x1, &f1, &x2, &f2 )
       != UNUR_SUCCESS )
    return x2;

  /* -- 3. bisection steps  */
  for (i=0; i<GEN->max_iter; i++) {
    /* compute new point */
    mid = x1 + (x2-x1)/2.;
    fmid = CDF(mid) - u;

    /* update point */
    if (f1*fmid <= 0) {
      /* keep left endpoint x1 */
      x2 = mid; f2 = fmid;
      /* check accuracy goal */
      if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
			       mid, fmid, x1, f1 ))
	break;
    }
    else {
      /* keep right endpoint x2 */
      x1 = mid; f1 = fmid;
      /* check accuracy goal */
      if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
			       mid, fmid, x2, f2 ))
	break;
    }
  }
  
  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
 
  /* the text book algorithms usually add one more step:              */
  /*   mid = x1 + (x2-x1)/2.;                                         */
  /* however, this causes problems when the maximal tolerated u-error */
  /* has to be ensured!                                               */

  /* ensure location within given (truncated) domain */
  mid = _unur_max( mid, DISTR.trunc[0]);
  mid = _unur_min( mid, DISTR.trunc[1]);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file (in case error) */
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample( gen,u,mid,CDF(mid)-u,i );
#endif

  /* finished */
  return mid;

} /* end of _unur_ninv_bisect() */


/*****************************************************************************/
/**  Auxiliary routines                                                     **/
/*****************************************************************************/

int
_unur_ninv_bracket( const struct unur_gen *gen, double u, 
		    double *xl, double *fl, double *xu, double *fu )
     /*----------------------------------------------------------------------*/
     /* find a bracket (enclosing interval) for root of CDF(x)-u.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   u      ... random number (uniform distribution)                    */
     /*   xl, xu ... store lower and upper boundary of bracket               */
     /*   fl, fu ... store values of CDF(x)-u at xl and xu, resp.            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   *fu is set to "best" result                                        */
     /*----------------------------------------------------------------------*/
{ 
  int i;                 /* loop variable, index                             */
  double x1, x2, xtmp;   /* points for bracket                               */
  double f1, f2, ftmp;   /* function values at x1, x2, xtmp                  */
  double step;           /* enlarges interval til sign change found          */
  int step_count;        /* counts number of steps finding sign change       */

  /* -- 1. initialize starting interval -- */

  if (GEN->table_on) {
    /* -- 1a. use table -- */

    /* 0 <= i <= table_size-2  */
    if ( _unur_FP_same(GEN->CDFmin, GEN->CDFmax) ) {
      /* CDF values in table too close, so we use median point since */
      /* there is no difference between CDF values.  */
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (u - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }

    /* set starting point for regular falsi */
    if ( ! _unur_FP_is_minus_infinity(GEN->table[i]) ){
      x1 = GEN->table[i];
      f1 = GEN->f_table[i]; 
    }
    else{
      x1 = GEN->table[i+1] + (GEN->table[i+1] - GEN->table[i+2]);
      f1 = CDF(x1);
    }

    if( ! _unur_FP_is_infinity(GEN->table[i+1]) ){
      x2 = GEN->table[i+1];
      f2 = GEN->f_table[i+1];
    }
    else{
      x2 = GEN->table[i] + (GEN->table[i] - GEN->table[i-1]);
      f2 = CDF(x2);
    }
  }

  else { 
    /* 1b. -- no table available -- */
    x1 =  GEN->s[0];      /* left boudary of interval */
    f1 =  GEN->CDFs[0];
    x2 =  GEN->s[1];      /* right boudary of interval*/   
    f2 =  GEN->CDFs[1];
  }

  /*  -- 1c. check for ordering of starting points -- */

  if ( x1 >= x2 ) { 
    xtmp = x1; ftmp = f1;
    x1   = x2; f1   = f2;
    x2 = xtmp + fabs(xtmp)*DBL_EPSILON;
    f2 = CDF(x2); 
  }

  /* -- 1d. check for boundary of truncated domain -- */
 
  /* in case of truncated domain there might be better starting points */
  /* ! no problems with INFINITY !  */
  if ( x1 < DISTR.trunc[0] || x1 >= DISTR.trunc[1] ){
    x1 = DISTR.trunc[0];
    f1 = GEN->Umin;    /* = CDF(x1) */
  }
  if ( x2 > DISTR.trunc[1] || x2 <= DISTR.trunc[0] ){
    x2 = DISTR.trunc[1];
    f2 = GEN->Umax;    /* = CDF(x2) */
  }

  /* -- 1z. compute function values at interval boundaries -- */
  f1 -= u;  f2 -= u;

  /* -- 2. search for enclosing bracket (interval where f changes signs) -- */
 
  step = (GEN->s[1]-GEN->s[0]) * STEPFAC;
  step_count = 0;
  while ( f1*f2 > 0. ) {
    /* interval too small -> make it bigger ( + 2^n * gap ) -- */
    if ( f1 > 0. ) { /* lower boundary too big */    
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1  = CDF(x1) - u;
    }
    else {         /* upper boundary too small */
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2  = CDF(x2) - u;
    }

    /* increase step width */
    if (step_count < MAX_STEPS) {
      ++step_count;
      step *= 2.;
      /* safe guard for the case where (GEN->s[1]-GEN->s[0]) is very small*/
      if( step_count > 20 && step < 1.) step = 1.; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
                  "Regula Falsi cannot find interval with sign change");
      *xu = (f1>0.) ? DISTR.trunc[0] : DISTR.trunc[1];
      return UNUR_ERR_GEN_SAMPLING;
    }
  }

  /* o.k. */
  *xl = x1; *xu = x2;
  *fl = f1; *fu = f2;

  return UNUR_SUCCESS;

} /* end of _unur_ninv_bracket() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_accuracy( const struct unur_gen *gen,
		     double x_resol, double u_resol,
		     double x0, double f0, double x1, double f1 )
     /*----------------------------------------------------------------------*/
     /* check accuracy goal for approximate root.                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x_resol ... maximal tolerated x-error                              */
     /*   u_resol ... maximal tolerated u-error                              */
     /*   x0      ... approximate point to be checked                        */
     /*   f0      ... function value at x0                                   */
     /*   x1      ... second boundary point of bracket                       */
     /*               (i.e., function must change sign in interval [x0,x1])  */
     /*   f1      ... function value at x1                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if accuracy goal is reached                              */
     /*   FALSE ... otherwise                                                */
     /*----------------------------------------------------------------------*/
{ 
  int x_goal, u_goal;     /* whether precision goal is reached               */

  if ( x_resol > 0. ) {
    /* check x-error */
    /* we use a combination of absolute and relative x-error: */
    /*    x-error < x-resolution * fabs(x) + x-resolution^2   */
    if ( _unur_iszero(f0) ||          /* exact hit */
  	 fabs(x1-x0) < x_resol * (fabs(x0) + x_resol) ) {
      x_goal = TRUE;
    }
    else if ( _unur_FP_same(f0,f1) ) {
      /* flat region */
      _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
  		    "flat region: accuracy goal in x cannot be reached");
      x_goal = TRUE;
    }
    else
      x_goal = FALSE;
  }
  else {
    /* no check */
    x_goal = TRUE;
  }
  
  if ( GEN->u_resolution > 0. ) {
    /* check u-error */
    /* (we use a slightly smaller maximal tolerated error than given by user) */
    if (fabs(f0) < 0.9 * u_resol) {
      u_goal = TRUE;
    }
    else if ( _unur_FP_same(x0,x1) ) {
      /* sharp peak or pole */
      _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
  		    "sharp peak or pole: accuracy goal in u cannot be reached");
      u_goal = TRUE;
    }
    else
      u_goal = FALSE;
  }
  else {
    u_goal = TRUE;
  }
  
  /* return result */
  return (x_goal && u_goal);

} /* end of _unur_ninv_accuracy() */

/*---------------------------------------------------------------------------*/

#undef MAX_STEPS
#undef STEPFAC
#undef I_CHANGE_TO_BISEC

/*---------------------------------------------------------------------------*/

