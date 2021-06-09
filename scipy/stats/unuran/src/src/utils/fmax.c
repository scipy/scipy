/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      fmax.c                                                       *
 *                                                                           *
 *   Find maximum of a function using Brent's algorithm.                     *
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

#include <unur_source.h>
#include "fmax_source.h"

/*---------------------------------------------------------------------------*/

double
_unur_util_find_max( struct unur_funct_generic fs, /* function structure */ 
                     double interval_min,   /* lower bound of interval   */
		     double interval_max,   /* upper bound of interval   */
		     double guess_max       /* initial guess for maximum */
		   )
     /*----------------------------------------------------------------------*/
     /* find max of univariate continuous distribution numerically           */
     /*                                                                      */
     /* This is achieved by the following steps:                             */
     /* -- Determine a interval containing the max and a third               */
     /*    point within ( x0< x1 < x2 )                                      */
     /*    ++ Determine a region where to search the max; this will be the   */
     /*       interval [max-100, max+100] if this is no contrdiction to      */
     /*       the given interval;`max' is the best known approx to the max   */
     /*    ++ Find a point in the interval with a positve pdf                */
     /*       This is done by two geometric sequences of MAX_SRCH elements   */
     /*       and find two other points within the given domain              */
     /*    ++ Unbounded domains: refine x0, x1, x2 until:                    */
     /*       f(x0) < f(x1) > f(x2)                                          */
     /* -- invoke a maximization-routine to determine the max                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   fs           ... function structure                                */
     /*   interval_min ... left boundary of function definition interval     */
     /*   interval_max ... right boundary of function definition interval    */
     /*   guess_max    ... initial guess for max position                    */
     /*                                                                      */
     /* return:                                                              */
     /*   x        ... approximate position of max on success                */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
#define MAX_SRCH (100)
 
  int i;

  double x[3];   /* max (and x[2]) should be between x[0] and  x[2] ...*/   
  double fx[3];  /* ... and the respective funtion values */
  double max;   /* (approximative) max of the distribution  */
  double max_l; /* lower bound for max search */
  double max_u; /* upper bound for max search */
  double step;

  int unbound_left; 
  int unbound_right; 

  /* first guess for max */
  max = (_unur_FP_is_infinity(guess_max)) ? 0. : guess_max;

  /* determine where to look for the max */
  
  /* unbounded domain */
  if ( _unur_FP_is_minus_infinity(interval_min) &&
       _unur_FP_is_infinity(interval_max) ){

    unbound_left = 1;
    unbound_right = 1;

    x[1]  = max;
    fx[1] = fs.f(x[1], fs.params);    
    max_l = max - 100.0;
    max_u = max + 100.0;

  }
  /* domain unbounded on the right */
  else if ( ! _unur_FP_is_minus_infinity(interval_min) &&
              _unur_FP_is_infinity(interval_max) ){

    unbound_left = 0;
    unbound_right = 1;

    if ( max >= interval_min ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_min;
      max_u = 2 * max - interval_min;
    }
    else{
      x[1] = interval_min + 100.0;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_min;
      max_u = x[1] + 100.0;
    }

  }
  /* domain unbounded on the left */
  else if ( _unur_FP_is_minus_infinity(interval_min) &&
          ! _unur_FP_is_infinity(interval_max) ){

    unbound_left = 1;
    unbound_right = 0;

    if ( max <= interval_max ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_max - 2 * max;
      max_u = interval_max;
    }
    else{
      x[1] = interval_max - 100.0;
      fx[1] = fs.f(x[1], fs.params);
      max_l = x[1] - 100.0;
      max_u = interval_max;      
    }

  }
  /* domain is bounded */
  else {  

    unbound_left = 0;
    unbound_right = 0;

    if ( max >= interval_min && max <= interval_max ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
    }
    else{
      x[1] = interval_min/2.0 + interval_max/2.0;
      fx[1] = fs.f(x[1], fs.params);
    }
    max_l = interval_min;
    max_u = interval_max;

  }

  max = x[1];  /* not exact mode -- best guess */


  /* find point with pdf > 0.0 -- max MAX_SRCH trials */

  /* search on the left side */
  step = pow(x[1]-max_l, 1.0/MAX_SRCH);
  i = 0;  
  while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
    x[1]  = max - pow(step, (double)i);
    fx[1] = fs.f(x[1], fs.params);
    i++;
  }
	 
  /* search on the right side */
  if( _unur_FP_same(0.0, fx[1]) ){
    step = pow(max_u-x[1], 1.0/MAX_SRCH);
    i = 0;
    while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
      x[1]  = max + pow(step, (double)i);
      fx[1] = fs.f(x[1], fs.params);
      i++;
    }
  }
  
  /* no success -- exit routine  */   
  if( _unur_FP_same(fx[1], 0.0) )
     return INFINITY; /* can't find max in flat region  */

  /* x[1] has f > 0 or routines already terminated */ 



  /* determine 3 points in the given domain --
     at least one with pdf > 0                  */
  if ( unbound_left ){

    x[2] = x[1];       fx[2] = fx[1];
    x[1] = x[2] - 1.0; fx[1] = fs.f(x[1], fs.params);
    x[0] = x[2] - 2.0; fx[0] = fs.f(x[0], fs.params);

  }
  else if ( unbound_right ){

    x[0] = x[1];       fx[0] = fx[1];
    x[1] = x[0] + 1.0; fx[1] = fs.f(x[1], fs.params);
    x[2] = x[0] + 2.0; fx[2] = fs.f(x[2], fs.params);

  }
  else{      /* bounded */

    x[0] = interval_min;  fx[0] = fs.f(x[0], fs.params);
    x[2] = interval_max;  fx[2] = fs.f(x[2], fs.params);

    if ( _unur_FP_same(x[1], interval_min) ||
         _unur_FP_same(x[1], interval_max) ){
      x[1]  = interval_min/2.0 + interval_max/2.0; 
      fx[1] = fs.f(x[1], fs.params);
    }
    
  }
  /* points x[i] with their function values determined */

  /* find interval containing the max */

  step = 1.0;
  if ( unbound_right ){
    while(fx[0] <= fx[1] && fx[1] <= fx[2]){ /* on the left side of the max */

      step *= 2.0;
      x[0]  = x[1]; fx[0] = fx[1];
      x[1]  = x[2]; fx[1] = fx[2];
      x[2] += step; fx[2] = fs.f(x[2], fs.params);   
    }
  }

  step = 1.0;  /* reset step size */
  if ( unbound_left ){
    while(fx[0] >= fx[1] && fx[1] >= fx[2]){ /* on the right side of the max */

      step *= 2.0;
      x[2]  = x[1]; fx[2] = fx[1];
      x[1]  = x[0]; fx[1] = fx[0];
      x[0] -= step; fx[0] = fs.f(x[0], fs.params);

    }
  }

  /* now: the max is between x[0] and x[2]   */

  /*    printf("x0: %f, fx0: %e\n", x[0], fx[0]); */
  /*    printf("x1: %f, fx1: %e\n", x[1], fx[1]); */
  /*    printf("x2: %f, fx2: %e\n", x[2], fx[2]); */

  /** TODO: FLT_MIN must be much larger than DBL_MIN **/

  max = _unur_util_brent( fs, x[0], x[2], x[1], FLT_MIN );
  if (!(_unur_FP_is_infinity( max )) ){
    /* mode successfully computed */

  }
  else {
    /* computing max did not work */
    return INFINITY; 
  }

  /* o.k. */
  return max; 

#undef MAX_SRCH
} /* end of _unur_util_find_max() */

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Brent algorithm for maximum-calculation of a continous function           */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*
 *****************************************************************************
 *	    		    C math library                                   *
 * function FMINBR - one-dimensional search for a function minimum           *
 *			  over the given range                               *
 *                                                                           *
 * Author: Oleg Keselyov.                                                    *
 *                                                                           *
 * modified by Josef Leydold (documentation unchanged)                       *
 * computed maximum instead of minimum.                                      *
 *                                                                           *
 * Input                                                                     *
 *	double fminbr(f, a,b,c,tol)                                          *
 *	double a; 			Minimum will be seeked for over      *
 *	double b;  			a range [a,b], a being < b.          *
 *      double c;                       c within (a,b) is first guess        *
 *	double (*f)(double x);		Name of the function whose minimum   *
 *					will be seeked for                   *
 *	double tol;			Acceptable tolerance for the minimum *
 *					location. It have to be positive     *
 *					(e.g. may be specified as EPSILON)   *
 *                                                                           *
 * Output                                                                    *
 *	Fminbr returns an estimate for the minimum location with accuracy    *
 *	3*SQRT_EPSILON*abs(x) + tol.                                         *
 *	The function always obtains a local minimum which coincides with     *
 *	the global one only if a function under investigation being          *
 *	unimodular.                                                          *
 *	If a function being examined possesses no local minimum within       *
 *	the given range, Fminbr returns 'a' (if f(a) < f(b)), otherwise      *
 *	it returns the right range boundary value b.                         *
 *                                                                           *
 * Algorithm                                                                 *
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical    *
 *	computations. M., Mir, 1980, p.202 of the Russian edition            *
 *                                                                           *
 *	The function makes use of the "gold section" procedure combined with *
 *	the parabolic interpolation.                                         *
 *	At every step program operates three abscissae - x,v, and w.         *
 *	x - the last and the best approximation to the minimum location,     *
 *	    i.e. f(x) <= f(a) or/and f(x) <= f(b)                            *
 *	    (if the function f has a local minimum in (a,b), then the both   *
 *	    conditions are fulfiled after one or two steps).                 *
 *	v,w are previous approximations to the minimum location. They may    *
 *	coincide with a, b, or x (although the algorithm tries to make all   *
 *	u, v, and w distinct). Points x, v, and w are used to construct      *
 *	interpolating parabola whose minimum will be treated as a new        *
 *	approximation to the minimum location if the former falls within     *
 *	[a,b] and reduces the range enveloping minimum more efficient than   *
 *	the gold section procedure.                                          *
 *	When f(x) has a second derivative positive at the minimum location   *
 *	(not coinciding with a or b) the procedure converges superlinearly   *
 *	at a rate order about 1.324                                          *
 *                                                                           *
 *****************************************************************************/

/* in case of any error INFINITY is returned */

double
_unur_util_brent(            /* An estimate to the min or max location */
     struct unur_funct_generic fs, /* Function struct                        */
     double a,                     /* Left border | of the range	     */
     double b,                     /* Right border| the min is seeked        */
     double c,                     /* first guess for the min/max            */
     double tol)                   /* Acceptable tolerance                   */
{
#define SQRT_EPSILON  (1.e-7)           /* tolerance for relative error      */
#define MAXIT         (1000)            /* maximum number of iterations      */

  /* compute maximum (instead of minimum as in original implementation ) */
#define f(x) ( (-1) * ((fs.f)(x, fs.params)) )          


  int i;

  double x,v,w;                         /* Abscissae, descr. see above       */
  double fx;                            /* f(x)                              */
  double fv;                            /* f(v)                              */
  double fw;                            /* f(w)                              */
  const double r = (3.-sqrt(5.0))/2;    /* Gold section ratio                */

  /* check arguments */
  CHECK_NULL(fs.f,INFINITY);
  if ( tol < 0. || b <= a || c <= a || b <= c) {
    _unur_error("CMAX",UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }


  /* Origially the third point was computed by golden section. In the        */
  /* modified version it is given as point `c' by the calling function.      */
  /* v = a + r*(b-a);  fv = f(v);        First step - always gold section    */

  v = c;  fv = f(v);                    /* First step */
  x = v;  w = v;
  fx=fv;  fw=fv;

  for(i=0; i < MAXIT; i++) {            /* Main iteration loop	*/
    double range = b-a;                 /* Range over which the minimum is   */
                                        /* seeked for		             */
    double middle_range = (a+b)/2;
    double tol_act =                    /* Actual tolerance                  */
		SQRT_EPSILON*fabs(x) + tol/3;
    double new_step;                    /* Step at this iteration            */
       
    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;                         /* Acceptable approx. is found       */

                                        /* Obtain the gold section step      */
    new_step = r * ( x<middle_range ? b-x : a-x );


                                 /* Decide if the interpolation can be tried */
    if( fabs(x-w) >= tol_act ) {        /* If x and w are distinct           */
                                        /* interpolatiom may be tried	     */
      register double p;                /* Interpolation step is calculated  */
      register double q;                /* as p/q; division operation        */
                                        /* is delayed until last moment      */
      register double t;

      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);

      if( q>(double)0 )                 /* q was calculated with the         */
	p = -p;                         /* opposite sign; make q positive    */
      else                              /* and assign possible minus to	p    */
	q = -q;

      if( fabs(p) < fabs(new_step*q) && /* If x + p/q falls in [a,b] not too */
	  p > q*(a-x+2*tol_act)      && /* close to a and b, and isn't too   */
	  p < q*(b-x-2*tol_act)  )      /* large, it is accepted.            */
	new_step = p/q;                 /* If p/q is too large then the	     */
                                        /* gold section procedure can reduce */
                                        /* [a,b] range to more extent.       */
    }

    if( fabs(new_step) < tol_act ) {    /* Adjust the step to be not less    */
      if( new_step > (double)0 )        /* than tolerance.                   */
	new_step = tol_act;
      else
	new_step = -tol_act;
    }
                                 /* Obtain the next approximation to min     */
    {                            /* and reduce the enveloping range.         */
      register double t = x + new_step;	/* Tentative point for the min       */
      register double ft = f(t);

      if( ft <= fx ) {                  /* t is a better approximation	     */
	if( t < x )   
	  b = x;                        /* Reduce the range so that          */
	else                            /* t would fall within it            */
	  a = x;
      
	v = w;  w = x;  x = t;          /* Assign the best approx to x	     */
	fv=fw;  fw=fx;  fx=ft;
      }
      else {                            /* x remains the better approx       */
	if( t < x )
	  a = t;                        /* Reduce the range enclosing x      */
	else
	  b = t;
      
        if( ft <= fw || _unur_FP_same(w,x) ) {
	  v = w;  w = t;
	  fv=fw;  fw=ft;
        }
        else if( ft<=fv || _unur_FP_same(v,x) || _unur_FP_same(v,w) ) {
	  v = t;
	  fv=ft;
        }
      }
    }                   /* ----- end-of-block ----- */
  }                /* ===== End of for loop ===== */

  /* maximal number of iterations exceeded */
  return INFINITY;

#undef f
#undef MAXIT
#undef SQRT_EPSILON
} /* end of _unur_util_brent() */

/*---------------------------------------------------------------------------*/
