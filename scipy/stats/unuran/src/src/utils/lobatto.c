/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lobatto.c                                                    *
 *                                                                           *
 *   Routines for Gauss-Lobatto integration with 5 points.                   *
 *   We store integral values for subintervals to speed up repeated          *
 *   computations.                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009 Wolfgang Hoermann and Josef Leydold                  *
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
#include "lobatto_source.h"
#include "lobatto_struct.h"

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */

static double 
_unur_lobatto5_simple (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		       double x, double h, double *fx);
/*---------------------------------------------------------------------------*/
/* numerical integration of 'funct' over the interval (x,x+h)                */
/* using Gauss-Lobatto integration with 5 points. (non-adaptive)             */
/*---------------------------------------------------------------------------*/

static double
_unur_lobatto5_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			 double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			 struct unur_lobatto_table *Itable);
/*---------------------------------------------------------------------------*/
/* numerical integration of 'funct' over the interval (x,x+h) using          */
/* adaptive Gauss-Lobatto integration with 5 points for each recursion.      */
/*---------------------------------------------------------------------------*/

static double 
_unur_lobatto5_recursion (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			  double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			  double int1, double fl, double fr, double fc,
			  int *W_accuracy, struct unur_lobatto_table *Itable);
/*---------------------------------------------------------------------------*/
/* run recursion for adaptive Lobatto integration.                           */
/*---------------------------------------------------------------------------*/

static int
_unur_lobatto_table_append (struct unur_lobatto_table *Itable, double x, double u);
/*---------------------------------------------------------------------------*/
/* append entry to table of integral values.                                 */
/*---------------------------------------------------------------------------*/

static void
_unur_lobatto_table_resize (struct unur_lobatto_table *Itable);
/*---------------------------------------------------------------------------*/
/* resize table of integral values.                                          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define FKT(x)  (funct((x),gen))      /* call to integrand */

/*---------------------------------------------------------------------------*/
/* points for Lobatto integration */

#define W1 (0.17267316464601146)   /* = 0.5-sqrt(3/28) */
#define W2 (1.-W1)

/*---------------------------------------------------------------------------*/

double
_unur_lobatto5_simple (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		       double x, double h, double *fx)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of 'funct' over the interval (x,x+h)           */
     /* using Gauss-Lobatto integration with 5 points. (non-adaptive)        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   funct ... integrand                                                */
     /*   gen   ... pointer to generator object                              */
     /*   x     ... left boundary point of interval                          */
     /*   h     ... length of interval                                       */
     /*   fx    ... funct(x) (ignored if NULL or *fx<0)                      */
     /*             set *fx <- funct(x+h)                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*                                                                      */
     /* store:                                                               */
     /*   funct(x+h) in *fx  if fx!=NULL                                     */
     /*----------------------------------------------------------------------*/
{ 
  double fl, fr;

  if (fx==NULL) {
    fl = FKT(x);
    fr = FKT(x+h);
  }
  else {
    fl = (*fx>=0.) ? *fx : FKT(x);
    fr = *fx = FKT(x+h);
  }

  return (9*(fl+fr)+49.*(FKT(x+h*W1)+FKT(x+h*W2))+64*FKT(x+h/2.))*h/180.;
} /* end of _unur_lobatto5_simple() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			double x, double h, double tol, UNUR_LOBATTO_ERROR uerror)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of the 'funct' over the interval (x,x+h)       */
     /* using adaptive Gauss-Lobatto integration with 5 points.              */
     /*----------------------------------------------------------------------*/
{
  return _unur_lobatto5_adaptive(funct,gen,x,h,tol,uerror,NULL); 
} /* end of _unur_lobatto_adaptive() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto5_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen, 
			 double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			 struct unur_lobatto_table *Itable)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of the 'funct' over the interval (x,x+h)       */
     /* using adaptive Gauss-Lobatto integration with 5 points.              */
     /*                                                                      */
     /* Halfs intervals recursively if error is too large.                   */
     /*                                                                      */
     /* The recursion stops when the ABSOLUTE error is less than the         */
     /* respective given tolerances.                                         */
     /*                                                                      */
     /* As a side effect, it stores boundaries and integrals for all         */
     /* subintervals in each recursion step where no further adaptation is   */
     /* required.                                                            */
     /* The values are stored in Itable (unless it equals NULL) in           */
     /* chronological order. It is important to take care when calling this  */
     /* function for several subintervals and 'Itable' is given.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   funct  ... integrand                                               */
     /*   gen    ... pointer to generator object                             */
     /*   x      ... left boundary point of interval                         */
     /*   h      ... length of interval                                      */
     /*   tol    ... tolerated ABSOLUTE error                                */
     /*   Itable ... table for storing integral values (may be NULL)         */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{
  double fl, fc, fr;  /* values of PDF at x, x+h/2, and x+h */
  double int1, int2;  /* estimated values for integral */
  int W_accuracy = FALSE; /* raise flag for printing warning about accuracy */

  /* check length of interval */
  if (_unur_iszero(h))
    return 0.;

  /* arguments which are not finite (inf or NaN) cause infinite recursions */
  if (!_unur_isfinite(x+h)) {
    _unur_error(gen->genid,UNUR_ERR_INF,"boundaries of integration domain not finite");
    return INFINITY;
  }

  /* compute function values */
  fl = FKT(x);
  fc = FKT(x+h/2.);
  fr = FKT(x+h);
 
  /* first estimate for integral on [x,x+h] */
  int1 = (9*(fl+fr)+49.*(FKT(x+h*W1)+FKT(x+h*W2))+64*fc)*h/180.;

  /* run adaptive steps */
  int2 = _unur_lobatto5_recursion(funct,gen,x,h,tol,uerror,int1,fl,fc,fr,&W_accuracy,Itable);

  /* return result */
  if (W_accuracy)
    _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		  "numeric integration did not reach full accuracy");
  return int2;

} /* end of _unur_lobatto5_adaptive() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto5_recursion (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			  double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			  double int1, double fl, double fc, double fr,
			  int *W_accuracy, struct unur_lobatto_table *Itable)
     /*----------------------------------------------------------------------*/
     /* run recursion for adaptive Lobatto integration.                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   funct    ... integrand                                             */
     /*   gen      ... pointer to generator object                           */
     /*   x        ... left boundary point of interval                       */
     /*   h        ... length of interval                                    */
     /*   tol      ... tolerated ABSOLUTE error                              */
     /*   fl       ... PDF at x                                              */
     /*   fc       ... PDF at x+h/2                                          */
     /*   fr       ... PDF at x+h                                            */
     /*   W_accuracy.. warning about accuracy                                */
     /*   Itable   ... table for storing integral values (may be NULL)       */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{
  double flc, frc;    /* values at PDF at x+h/4 and x+3*h/4 */
  double int2;        /* estimated values of integrals */
  double intl, intr;  /* left and right part of int2 */
  double ierror;      /* integration error */

  /* compute function values */
  flc = FKT(x+h/4);
  frc = FKT(x+3*h/4);
 
  /* compute integral on [x,x+h/2] and on [x+h/2,x+h] */
  intl = (9*(fl+fc)+49.*(FKT(x+h*W1*0.5)+FKT(x+h*W2*0.5))+64*flc)*h/360.;
  intr = (9*(fc+fr)+49.*(FKT(x+h*(0.5+W1*0.5))+FKT(x+h*(0.5+W2*0.5)))+64*frc)*h/360.;
  int2 = intl + intr;

  /* integration error */
  if (uerror!=NULL)
    ierror = uerror(gen, fabs(int1-int2), x+h/2.);
  else 
    ierror = fabs(int1-int2);

  /* check whether accuracy goal is reached */
  if (ierror >= tol) {
    /* error above tolerance */
    /* if (_unur_FP_equal(x+h/2.,x) || _unur_FP_equal(int1,int2) ) { */
    if (_unur_FP_equal(x+h/2.,x)) {
      /* we cannot decrease length of subintervals any more */
      *W_accuracy = TRUE;
      /* Remark: Since we are halving intervals, this comparision */
      /* limits the maximal number of iterations to at most 2048. */
    }
    else {
      /* recompute with shorter intervals */
      /* Remark: it is important that the two calls to
	 _unur_lobatto5_recursion() must be done in exactly this order!
	 Otherwise 'Itable' gets corrupted.
	 So we need to two seperate statements and store the results
	 in the temporary variabel 'int2' to prevent the compiler to
	 revert this order.
       */
      int2  = _unur_lobatto5_recursion(funct,gen,x,h/2,tol/1.,uerror,
				      intl,fl,flc,fc, W_accuracy,Itable);
      int2 += _unur_lobatto5_recursion(funct,gen,x+h/2,h/2,tol/1.,uerror,
				       intr,fc,frc,fr, W_accuracy,Itable);
      return int2;
    }
  }

  /* store integral values */
  if (Itable) {
    /* l.h.s. subinterval */
    _unur_lobatto_table_append(Itable, x+h/2., intl);
    /* r.h.s. subinterval */
    _unur_lobatto_table_append(Itable, x+h, intr);
  }
  /* Remark: we do not throw a warning if the table size is exceeded. */

  /* return estimate for integral */
  return int2;

} /* end of _unur_lobatto5_recursion() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto_eval_diff (struct unur_lobatto_table *Itable, double x, double h, double *fx)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of 'funct' over the interval (x,x+h) using     */
     /* table of integral values together with (adaptive) Gauss-Lobatto      */
     /* integration with 5 points.                                           */
     /*                                                                      */
     /* It is important that the pointer points to an entry in the table of  */
     /* integral values with x-value larger than the given 'x'.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable ... table for storing integral values (may be NULL)         */
     /*   x      ... left boundary point of interval                         */
     /*   h      ... length of interval                                      */
     /*   fx     ... funct(x) (ignored if NULL or *fx<0)                     */
     /*              set *fx <- funct(x+h) if _unur_lobatto5_simple called   */
     /*              set *fx <- -1.        otherwise                         */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*                                                                      */
     /* store:                                                               */
     /*   if (fx!=NULL)                                                      */
     /*       *fx = funct(x+h)      if _unur_lobatto5_simple is called       */
     /*       *fx = -1. (=unknown)  if _unur_lobatto5_adaptive is called     */
     /*----------------------------------------------------------------------*/
{
  int cur;                    /* pointer to current position in table */
  double x1;                  /* left and right boundary of a subinterval */
  double Q;                   /* value of integral */
  struct unur_lobatto_nodes *values;
  int n_values;

  /* clear function values if it does not contain the correct value. */
  /* this is in particular necessary whenever we call                */
  /* _unur_lobatto5_adaptive instead if _unur_lobatto5_simple !!     */
#define clear_fx() if(fx!=NULL){*fx=-1.;}

  /* check for invalid NULL pointer */
  CHECK_NULL(Itable,INFINITY);

  /* read data from table */
  values = Itable->values;
  n_values = Itable->n_values;

  /* arguments which are not finite (inf or NaN) cause infinite recursions */
  if (!_unur_isfinite(x+h)) {
    /* _unur_warning(gen->genid,UNUR_ERR_INF,"boundaries of integration domain not finite"); */
    clear_fx();
    return INFINITY;
  }

  /* check for boundaries of integration table */
  if (x < Itable->bleft || x+h > Itable->bright) {
    clear_fx();
    return _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x, h, 
				   Itable->tol, Itable->uerror, NULL);
  }

  /* move pointer for reading to start position in interval */
  cur = Itable->cur_iv;

  /* first entry in interval */
  while (cur < n_values && values[cur].x < x)
    ++cur;

  /* did we find such an entry ? */
  if (cur >= n_values) {
    /* we must use adaptive Lobatto integration if the table for  */
    /* integral values was too small.                             */
    clear_fx();
    return _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x, h, 
				   Itable->tol, Itable->uerror, NULL);
  }

  /* store x value and goto next entry */
  x1 = values[cur].x;
  ++cur;

  /* are there more than one entry in interval ? */
  if (cur >= n_values ||
      values[cur].x > x+h) {
    /* there is at most one entry in the interval [x,x+h]. */
    /* thus we (can) use simple Lobatto integration.       */
    return _unur_lobatto5_simple(Itable->funct, Itable->gen, x, h, fx);
  }

  /* else:                                                         */
  /* there are more than one entries in the interval.              */
  /* thus there is at least one subinterval from adapative Lobatto */
  /* integration within [x,x+h]. we reuse the stored values.       */
  /* for the remaining two subintervals at the boundary [x,x+h]    */
  /* we use simple Gauss-Lobatto integration.                      */

  Q = _unur_lobatto5_simple(Itable->funct, Itable->gen, x, x1-x, fx);
  do {
    Q += values[cur].u;
    /*      = _unur_lobatto5(gen, funct, x1, x2-x1) */
    x1 = values[cur].x;
    ++cur;
  } while (cur < n_values && values[cur].x <= x+h);
  /* now *fx does not hold PDF(x1). so we have to clear it */
  clear_fx();

  /* We have to distinguish two cases: */
  if (cur >= n_values) {
    /* the table of integral values is too small */
    Q += _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x1, x+h-x1,
				 Itable->tol, Itable->uerror, NULL);
  }
  else {
    /* the table is not too small but x+h is outside the computational domain */
    Q += _unur_lobatto5_simple(Itable->funct, Itable->gen, x1, x+h-x1, fx);
  }
      
  return Q;

#undef clear_fx

} /* end of _unur_lobatto_eval_diff() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto_eval_CDF (struct unur_lobatto_table *Itable, double x)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of 'funct' over the interval (-INFINITY, x)    */
     /* using a table of integral values together with                       */
     /* (adaptive) Gauss-Lobatto integration with 5 points.                  */
     /* It is important, that the integration object 'Itable' already exists.*/
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable ... table for storing integral values (may be NULL)         */
     /*   x      ... left boundary point of interval                         */
     /*   h      ... length of interval                                      */
     /*   fx     ... funct(x) (ignored if NULL or *fx<0)                     */
     /*              set *fx <- funct(x+h) if _unur_lobatto5_simple called   */
     /*              set *fx <- -1.        otherwise                         */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{
  struct unur_lobatto_nodes *values;
  int n_values;          /* table size */
  double area;           /* integral over (-INFINITY, INFINITY) */
  double xr;             /* most right point from table */
  double cdf;            /* cdf at x */
  int cur;               /* current interval (position in table) */

  /* check for invalid NULL pointer */
  CHECK_NULL(Itable,INFINITY);

  /* check boundary */
  if (x <= Itable->bleft)  return 0.;
  if (x >= Itable->bright) return 1.;

  /* read data from table */
  values = Itable->values;
  n_values = Itable->n_values;
  area = Itable->integral;

  /* the area must not equal 0 (this should not happen anyway) */
  if (area <= 0.) {
    _unur_error(Itable->gen->genid,UNUR_ERR_NAN,"area below PDF 0.");
    return INFINITY;
  }

  /* sum values over all intervals that are on the l.h.s. of x */
  cdf = 0;
  xr = Itable->bleft;
  for (cdf=0, cur=0; cur < n_values && x > values[cur].x; cur++) {
    cdf += values[cur].u;
    xr = values[cur].x;
  }

  /* integrate from xr to x */
  if (cur >= n_values) {
    cdf += _unur_lobatto5_adaptive(Itable->funct, Itable->gen, xr, x-xr,
				   Itable->tol, Itable->uerror, NULL);
  }
  else {
    cdf += _unur_lobatto5_simple(Itable->funct, Itable->gen, xr, x-xr, NULL);
  }

  /* compute CDF */
  cdf /= area;
  cdf = _unur_max(0., cdf);
  cdf = _unur_min(1., cdf);

  /* return CDF */
  return cdf;

} /* end of _unur_lobatto_eval_CDF() */

/*---------------------------------------------------------------------------*/

double
_unur_lobatto_integral (struct unur_lobatto_table *Itable)
     /*----------------------------------------------------------------------*/
     /* Get value of integral from Lobatto object.                           */
     /* Get integration from Lobatto object.                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable ... table for storing integral values (may be NULL)         */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{
  CHECK_NULL(Itable,INFINITY);
  return Itable->integral;
} /* end of _unur_lobatto_eval_integral() */


/*****************************************************************************/
/*                                                                           */
/*   Table with integral values                                              */
/*                                                                           */
/*****************************************************************************/

struct unur_lobatto_table *
_unur_lobatto_init (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		    double left, double center, double right,
		    double tol, UNUR_LOBATTO_ERROR uerror, int size)
     /*----------------------------------------------------------------------*/
     /* create and initialize table of integral values.                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   funct ... pointer to integrand                                     */
     /*   gen   ... pointer to generator object                              */
     /*   left  ... left boundary of computational domain                    */
     /*   center... central ("typical") point of computational domain        */
     /*   right ... right boundary of computational domain                   */
     /*   tol   ... tolerated ABSOLUTE integration error                     */
     /*   size  ... (maximal) size of table (MUST BE >=2)                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to table                                                   */
     /*----------------------------------------------------------------------*/
{
  struct unur_lobatto_table *Itable;

  /* check argument */
  if (size<2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"size<2");
    return NULL;
  }

  /* allocate memory */
  Itable = _unur_xmalloc( sizeof(struct unur_lobatto_table) );
  Itable->values = _unur_xmalloc(size * sizeof(struct unur_lobatto_nodes) );

  /* set counter */
  Itable->size = size;
  Itable->n_values = 0;
  Itable->cur_iv = 0;

  /* store integrand */
  Itable->funct = funct;
  Itable->gen = gen;
  Itable->bleft = left;
  Itable->bright = right;

  /* tolerated integration error */
  Itable->tol = tol;
  Itable->uerror = uerror;

  /* store left boundary point in table */
  _unur_lobatto_table_append(Itable,left,0.);

  /* compute integral over whole domain                */
  /* Remark: it is important that the two calls to
     _unur_lobatto5_adaptive() must be done in exactly this order!
     Otherwise 'Itable' gets corrupted.
     So we need two seperate statements to prevent the compiler to
     revert this order.
  */
  Itable->integral = 
    _unur_lobatto5_adaptive(funct, gen, left, center-left, tol, uerror, Itable );
  Itable->integral += 
    _unur_lobatto5_adaptive(funct, gen, center, right-center, tol, uerror, Itable );

  /* possibly shrink memory block for table of integral values */
  _unur_lobatto_table_resize(Itable);

  return Itable;
} /* end of _unur_lobatto_table_create() */

/*---------------------------------------------------------------------------*/

int _unur_lobatto_find_linear (struct unur_lobatto_table *Itable, double x)
     /*----------------------------------------------------------------------*/
     /* Find first subinterval where left boundary is not less than x.       */ 
     /* we start at last position found by this routine and continue with    */
     /* linear search.                                                       */
     /* (We set a bookmark in the table of integral values.)                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itables ... table with integral values                             */
     /*   x       ... value to be found                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (Itable != NULL) {
    /* search for first entry in interval. */
    /* we can continue from the position in the last interval. */
    /* otherwise, to restart from the first entry uncomment this line */
    /*    Itable->cur_iv = 0; */
    
    while (Itable->cur_iv < Itable->n_values &&
	   Itable->values[Itable->cur_iv].x < x) 
      ++(Itable->cur_iv);

    return UNUR_SUCCESS;
  }
  else {
    /* nothing to do */
    return UNUR_ERR_SILENT;
  }

} /* end of _unur_lobatto_find_linear() */

/*---------------------------------------------------------------------------*/

int
_unur_lobatto_table_append (struct unur_lobatto_table *Itable, double x, double u)
     /*----------------------------------------------------------------------*/
     /* append entry to the end of the table of integral values.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itables ... table with integral values                             */
     /*   x       ... right boundary of subinterval                          */
     /*   u       ... integral over subinterval                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (Itable==NULL)
    return UNUR_ERR_NULL;

  if (Itable->n_values >= Itable->size - 1)
    /* we do not write a warning here */
    return UNUR_ERR_GENERIC;

  Itable->values[Itable->n_values].x = x;
  Itable->values[Itable->n_values].u = u;
  ++(Itable->n_values);

  return UNUR_SUCCESS;
} /* end of _unur_lobatto_table_append() */

/*---------------------------------------------------------------------------*/

void
_unur_lobatto_table_resize (struct unur_lobatto_table *Itable)
     /*----------------------------------------------------------------------*/
     /* resize table of integral values.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable ... pointer to pointer to table of integral values          */
     /*----------------------------------------------------------------------*/
{
  if (Itable) {
    Itable->size = Itable->n_values;
    Itable->values = _unur_xrealloc(Itable->values,
				    Itable->size * sizeof(struct unur_lobatto_nodes));
  }
  /* else: nothing to do */
  
} /* end of _unur_lobatto_table_resize() */

/*---------------------------------------------------------------------------*/

void
_unur_lobatto_free (struct unur_lobatto_table **Itable)
     /*----------------------------------------------------------------------*/
     /* destroy table of integral values and set pointer to NULL.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable ... pointer to pointer to table with integral values        */
     /*----------------------------------------------------------------------*/
{
  if (*Itable) {
    free ((*Itable)->values);
    free (*Itable);
    *Itable = NULL;
  }
  /* else: nothing to do */

} /* end of _unur_lobatto_free() */

/*---------------------------------------------------------------------------*/

void
_unur_lobatto_debug_table (struct unur_lobatto_table *Itable, const struct unur_gen *gen,
			   int print_Itable )
     /*----------------------------------------------------------------------*/
     /* print size and entries of table of integral values.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable       ... table with integral values                        */
     /*   gen          ... pointer to generator object                       */
     /*   print_Itable ... whether table is printed                          */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int n;

  /* check arguments */
  CHECK_NULL(Itable,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: subintervals for Lobatto integration: %d\n",gen->genid,
	  Itable->n_values - 1);

  for (n=0; print_Itable && n < Itable->n_values; n++) {
    fprintf(LOG,"%s:  [%3d] x = %g, u = %g\n",gen->genid,
	    n, Itable->values[n].x, Itable->values[n].u );
  }

} /* end of _unur_lobatto_debug_table() */

/*---------------------------------------------------------------------------*/

int _unur_lobatto_size_table (struct unur_lobatto_table *Itable)
     /*----------------------------------------------------------------------*/
     /* size of table of integral values.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   Itable       ... table with integral values                        */
     /*----------------------------------------------------------------------*/
{
  return (Itable->n_values - 1);
} /* end of _unur_lobatto_size_table() */

/*---------------------------------------------------------------------------*/
