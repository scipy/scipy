/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_newton.c                                                *
 *                                                                           *
 *   Routines for Newton interpolating polynomial.                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008-2010 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/** Compute coefficients for Newton polynomial                              **/
/*****************************************************************************/

int
_unur_pinv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Create table for Newton interpolation.                               */
     /* The computational domain must be already computed / given.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double utol;               /* tolerated maximum U-error */
  double maxerror;           /* maximum U-error in a particular interval */ 
  double h;                  /* step size / length of intervals */
  int i;                     /* number of interval at work */
  int iter;                  /* number of iterations */
  int cont;                  /* whether we have to continue or loop */

  int smooth;                /* temporary smoothness parameter */
  int use_linear;            /* whether we use linear interpolation instead of polynomials */
  int right_bd;              /* whether we have reach right boundary */
  int use_upoints;           /* whether we use Chebyshev points in u-scale */

  double chebyshev[3][MAX_ORDER+1]; /* Chebyshev points */
  double xval[MAX_ORDER+1];  /* x-values for construction points for Newton polynomial */

  int n_decr_h = 0;          /* number of steps where h is decreased */
  int n_incr_h = 0;          /* number of steps where h is increased */
  int n_use_linear = 0;      /* number of steps where linear interpolation is used */


  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* tolerated U-error */
  utol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION;

  /* initialize step size for subintervals */
  h = (GEN->bright-GEN->bleft)/128.;

  /* initialize array of intervals: starting interval */
  if (_unur_pinv_interval( gen, 0, GEN->bleft, 0.) != UNUR_SUCCESS) 
    return UNUR_ERR_GEN_CONDITION;

  /* compute construction points for Chebyshev points */
  for (smooth=0; smooth<=GEN->smooth; ++smooth)
    _unur_pinv_chebyshev_points(chebyshev[smooth],GEN->order,smooth);

  /* initialize counter and control variables */
  i = 0;                /* start at interval 0 ;-) */
  cont = TRUE;          /* there is at least one iteration */
  use_linear = FALSE;   /* do not use linear interpolation; use polynomial */

  use_upoints = FALSE;  /* do not use Chebyshev points in u-scale.
			 * we use this variable to toggle between 
			 * Chebyshev points in x-scale and u-scale. */

  /* compute intervals and coefficients of Newton polynomials */
  for (iter=0; cont ; iter++) {

    /* check number of iterations */
    if (iter >= PINV_MAX_ITER_IVS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "maximum number of iterations exceeded");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* right boundary reached ? */
    if(!_unur_FP_less(GEN->iv[i].xi+h,GEN->bright)) {
      h = GEN->bright - GEN->iv[i].xi;
      cont = FALSE;    /* we probably can stop after this iteration */
      right_bd = TRUE; /* we have reached right boundary */
    }
    else {
      right_bd = FALSE;
    }

    /* set temporary smoothness parameter */
    smooth = GEN->smooth;

    /* compute interpolating polynomial */
    use_linear = FALSE;
    switch (smooth) {

    case 2:
      _unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
      if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	break;

      smooth = 1;

    case 1:
      if (GEN->order % 2 == 1) {
	/* order must be odd when smoothness equals 1 */
	_unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
	if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	  break;
      }
      smooth = 0;

    case 0:
    default:
      /* compute Newton interpolation polynomial */
      _unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
      if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	break;

    case -1:
      /* area below PDF = 0 or serious round-off errors */
      /*   --> try linear interpolation                 */
      use_linear = TRUE;
    }
    
    if (use_linear) {
      /* use linear interpolation */
      ++n_use_linear;

      if (_unur_pinv_linear_create(gen,GEN->iv+i,xval) != UNUR_SUCCESS) {

	/* area below PDF == 0. */
	if (i==0) { /* left boundary (first interval) */
	  /* cut left boundary */
	  GEN->bleft = GEN->iv[i].xi + h;
	  GEN->iv[i].xi = GEN->bleft;
	  continue;  /* outer loop */
	}
	else if (right_bd) { /* right boudary */
	  GEN->bright = GEN->iv[i].xi;
	  cont = FALSE;
	  break;  /* outer loop */
	}
	else {
	  /* -- no idea what to do now --> abort */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PDF too close to 0 on relevant part of domain --> abort");
	  return UNUR_ERR_GEN_CONDITION;
	}
      }
    }

    /* estimate error of Newton interpolation */
    maxerror = _unur_pinv_newton_maxerror(gen,&(GEN->iv[i]),xval,use_linear);

    if (!(maxerror <= utol)) {
      /* error too large: reduce step size */
      h *= (maxerror > 4.*utol) ? 0.81 : 0.9;
      cont = TRUE;  /* we need another iteration */
      ++n_decr_h;
      use_upoints = FALSE;
      continue;
    }

    /* error is below threshold */

    if (gen->variant & PINV_VARIANT_UPOINTS) {
      /* compute polynomial again with (approximate) Chebyshev points in u-scale */
      if (!use_upoints && !use_linear ) {
	use_upoints = TRUE;
	cont = TRUE;
	continue;
      }
      use_upoints = FALSE;
    }

    /* error is below threshold --> continue with next interval */
    
    /* create next interval */
    if ( _unur_pinv_interval( gen, i+1, GEN->iv[i].xi+h, 
			      GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order-1])
	 /* cdfi holds CDF value at the left border of the interval,                  */
	 /* ui[order-1] holds area below PDF in interval, i.e. CDF(right) - CDF(left) */
	 != UNUR_SUCCESS )
      return UNUR_ERR_GEN_CONDITION;

    /* increase step size for very small errors */
    if (maxerror < 0.3*utol) {
      h *= (maxerror < 0.1*utol) ? 2. : 1.2;
      ++n_incr_h;
    }
      
    /* continue with next interval */
    i++;
  }

  /* update size of array (finish list) */
  _unur_pinv_lastinterval(gen);

  /* set range for uniform random numbers */
  /* Umin = 0, Umax depends on area below PDF, tail cut-off points and round-off errors */
  GEN->Umax = GEN->iv[GEN->n_ivs].cdfi;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_create_table(gen,iter,n_incr_h,n_decr_h,n_use_linear);
#endif

  /* o.k. */
  return UNUR_SUCCESS;
}  /* end of _unur_pinv_create_table() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_chebyshev_points (double *pt, int order, int smooth)
     /*----------------------------------------------------------------------*/
     /* [2a.] Compute Chebyshev points.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pt     ... pointer to array of size (order+1) for storing points   */
     /*   order  ... order of polynomial                                     */
     /*   smooth ... smoothness parameter                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
     /* Chebyshev points rescaled such that the first and last point is      */
     /* moved into 0 and 1, resp.                                            */
     /* We need (order+1) points for a polynomial of given 'order'.          */
     /* For smooth interpolant we need multiple points:                      */
     /*   smooth == 0: single points                                         */
     /*   smooth == 1: double points                                         */
     /*   smooth == 2: tripple points                                        */
     /*----------------------------------------------------------------------*/
{
  int i,k;

  int n_pt = (order+1)/(smooth+1); /* number of points                       */
  double phi = M_PI*0.5/n_pt;     /* constant for computing Chebyshev points */
  
  pt[0] = 0.;
  for (i=1; i<=order; i++) {
    k = i / (smooth+1);
    if (k<n_pt-1)
      pt[i] = sin(k*phi) * sin((k+1)*phi)/cos(phi);
    else
      pt[i] = 1.;
  }

  return UNUR_SUCCESS;
} /* end of _unur_pinv_chebyshev_points() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_cpoints (double *xval, int order, struct unur_pinv_interval *iv,
			   double h, double *chebyshev, int smooth, int use_upoints)
     /*----------------------------------------------------------------------*/
     /* [2a.] Compute points for construct interpolating polynomial.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   xval        ... for storing x-values for constructing polynomial   */
     /*   order       ... order of polynomial                                */
     /*   iv          ... pointer to current interval                        */
     /*   h           ... step size (width of interval)                      */
     /*   chebyshev   ... Chebyshev points                                   */
     /*   smooth      ... smoothness parameter                               */
     /*   use_upoints ... whether we use Chebyshev points in u-scale         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
     /* when the interpolant must be (twice) differentiable we need          */
     /* double (triple) points.                                              */
     /*----------------------------------------------------------------------*/
{
  int k;

  if (use_upoints) {
    /* FIXME: check monotone and != */

    /* use Chebyshev points in u-scale. */
    double hu = iv->ui[order-1];
    for(k=0; k<=order; k++)
      xval[k] = (k % (smooth+1))
       	? xval[k-1]
 	: iv->xi + _unur_pinv_newton_eval(hu * chebyshev[k], iv->ui, iv->zi, order);
  }
  
  else {
    /* use Chebyshev points in x-scale. */
    for(k=0; k<=order; k++) {
      xval[k] = (k % (smooth+1))
      	? xval[k-1]
      	: iv->xi + h * chebyshev[k];
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_cpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval)
     /*----------------------------------------------------------------------*/
     /* 2a. Compute coefficients for Newton interpolation within a           */
     /* subinterval of the domain of the distribution.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   xval ... x-values for constructing polynomial                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double fx = -1.;       /* PDF at x  (-1. --> PDF unknown) */
  double *ui = iv->ui;   /* u-values for Newton interpolation */
  double *zi = iv->zi;   /* coefficients of Newton interpolation */
  double xi, dxi;        /* boundary and length of i-th subinterval */
  double area;           /* integral of PDF over subinterval */
  int i,k;               /* auxiliary variables */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* compute coefficients of interpolating polynomial */

  /* first step: 
   * compute values u_i of nodes and first differences, i.e.,
   * ( CDF^{-1}(u[i]) - CDF^{-1}(u[i-1]) ) / (u[i] - u[i-1])
   * = (x[i] - x[i-1]) / (u[i] - u[i-1])
   */
  for(i=0; i<GEN->order; i++) {

    /* left boundary and length of subinterval for integration */
    xi = xval[i];

    if (!_unur_FP_same(xval[i],xval[i+1])) {
      /* use differences */
      dxi = xval[i+1]-xval[i];

      /* compute integral of PDF in interval (xi,xi+dxi) */
      area = _unur_pinv_Udiff(gen, xi, dxi, &fx);
      
      if (_unur_iszero(area)) return UNUR_ERR_SILENT;
      
      /* construction points for interpolating polynomial for CDF^{-1} */
      ui[i] = (i>0) ? (ui[i-1]+area) : area;
      /* divided differences of values of CDF^{-1} */ 
      zi[i] = dxi/area;
    }
    else {
      /* use derivative:
       *   [F^{-1}(u)]' = 1/f(x)  where x=F^{-1}(u)
       */
      ui[i] = (i>0) ? ui[i-1] : 0.;
      zi[i] = 1./PDF(xi);
    }
  }
  /* Remark: ui[GEN->order-1] is the probability of the interval */

  /* compute second divided differences. */
  /*   k = 1; */
  for(i=GEN->order-1; i>=1; i--) {
    if (!_unur_FP_same(zi[i],zi[i-1]))
      /* use differences */
      zi[i] = (i>1) 
	? (zi[i]-zi[i-1]) / (ui[i]-ui[i-2])
	: (zi[1]-zi[0]) / ui[1];
    else
      /* use second derivative:
       *   [F^{-1}(u)]" = -f'(x) / f(x)^3  where x=F^{-1}(u)
       */
      zi[i] = -0.5 * dPDF(xval[i]) * pow(zi[i],3);
  }

  /* compute all other divided differences. */
  for(k=2; k<GEN->order; k++) {
    for(i=GEN->order-1; i>k; i--) {
      zi[i] = (zi[i]-zi[i-1]) / (ui[i]-ui[i-(k+1)]);
    }
    zi[k] = (zi[k]-zi[k-1]) / ui[k];
  }

  /* check result */
  for (i=0; i<GEN->order; i++) {
    if (!_unur_isfinite(zi[i])) 
      return UNUR_ERR_SILENT;
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_linear_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval)
     /*----------------------------------------------------------------------*/
     /* [2a.] Compute coefficients for linear interpolation within a         */
     /* subinterval of the domain of the distribution.                       */
     /* The coefficients are stored in the array for Newton interpolation    */
     /* of order GEN->order where the coefficients of the higher order terms */
     /* are set to 0.                                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   xval ... x-values for constructing polynomial as used for          */
     /*            _unur_pinv_newton_create                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double *ui = iv->ui;   /* u-values for Newton interpolation */
  double *zi = iv->zi;   /* coefficients of Newton interpolation */
  double x0, x1;         /* lower and upper boundary point of interval */
  double area;           /* integral of PDF over subinterval */
  int i;                 /* auxiliary variable */

  /* boundary of interval */
  x0 = xval[0];
  x1 = xval[GEN->order];

  /* set all coefficents of the Newton polynomial to 0. */
  for (i=1; i<GEN->order; i++) {
    ui[i] = zi[i] = 0.;
  }

  /* area below PDF */
  area = _unur_pinv_Udiff(gen, x0, x1-x0, NULL);

  /* zi[0] contains the slope of the polynomial */
  ui[0] = area;
  zi[0] = (x1 - x0) / area;

  /* ui[GEN->order-1] stores the probability of the interval */
  ui[GEN->order-1] = area;

  /* check result */
  if (!_unur_isfinite(zi[0]))
    return UNUR_ERR_GEN_CONDITION;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_linear_create() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_eval ( double q, double *ui, double *zi, int order )
     /*----------------------------------------------------------------------*/
     /* 2b. evaluate Newton interpolation polynomial using Horner scheme.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   q     ... argument                                                 */
     /*   ui    ... coefficients of polynomial (increasing order)            */
     /*   zi    ... coefficients of polynomial (increasing order)            */
     /*   order ... order of polynomial                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   value of interpolating polynomial at u=q                           */
     /*----------------------------------------------------------------------*/
{
  int k;
  double chi;

  chi = zi[order-1];
  for (k=order-2; k>=0; k--)
    chi = chi*(q-ui[k])+zi[k];

  return (chi*q);
} /* end of _unur_pinv_newton_eval() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv,
			    double *xval, int is_linear)
     /*----------------------------------------------------------------------*/
     /* 2c. Estimate maximal error of Newton interpolation in subinterval.   */
     /*     In addition it makes a simple check for monotonicity of the      */
     /*     inverse CDF.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   xval ... x-values for constructing polynomial                      */
     /*   is_linear ... whether we have linear interpolation                 */
     /*                                                                      */
     /* return:                                                              */
     /*   estimated maximal u-error, or                                      */
     /*   >1000 whenever the inverse CDF is not monotone                     */
     /*----------------------------------------------------------------------*/
{
  double x0 = iv->xi;    /* left boundary point of interval */
  double *ui = iv->ui;   /* u-values for Newton interpolation  */
  double *zi = iv->zi;   /* coefficient of Newton interpolation  */

  double maxerror = 0.;  /* maximum error */
  double uerror;         /* error for given U value */
  double x;              /* x = CDF^{-1}(U) */
  double u;              /* u = CDF(x) */

  double testu[MAX_ORDER]; /* array of U values for testing */ 
  int i;                 /* aux variable */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* check for monotonicity (linear case) */
  if (is_linear && zi[0] < 0.)
    /* not monotone */
    return 1001.;

  /* get U values for test (points with maximal worst case error) */
  if (is_linear) 
    _unur_pinv_linear_testpoints(testu,ui,GEN->order);
  else
    _unur_pinv_newton_testpoints(testu,ui,GEN->order);

  /* calculate the max u-error at the test points */
  for(i=0; i<GEN->order; i++) {
    
    /* inverse CDF for U test point */
    x = _unur_pinv_newton_eval(testu[i], ui, zi, GEN->order);

    /* check for monotonicity (non-linear case) */
    if (!is_linear) {
      if (! (xval[i] <= x0+x && x0+x <= xval[i+1]) )
	if (! _unur_FP_same(xval[i], xval[i+1]))
	  /* not monotone */
	  return 1002.;
    }
    
    /* estimate CDF for interpolated x value */
    if (i==0 || xval==NULL)
      u = _unur_pinv_Udiff(gen, x0, x, NULL);
    else
      u = ui[i-1] + _unur_pinv_Udiff(gen, xval[i], x+x0-xval[i], NULL);

    /* check u-value */
    if (!_unur_isfinite(u))
      return INFINITY;

    /* compute u-error */
    uerror = fabs(u - testu[i]);

    /* update maximal error */
    if (uerror>maxerror) maxerror = uerror;
  }

  /* we have an analytical test for cubic Hermite interpolation */
  if (GEN->order == 3 && GEN->smooth==1 && xval!=NULL &&
      ! _unur_pinv_cubic_hermite_is_monotone(gen,ui,zi,xval))
    /* not monotone */
    return 1003.;
  
  /* return maximal observed u-error */
  return maxerror;
} /* end of _unur_pinv_newton_maxerror() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_testpoints (double *utest, double *ui, int order)
     /*----------------------------------------------------------------------*/
     /* [2c.] calculates the local maxima of the polynomial.                 */
     /* used as control points for error estimate.                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*    utest ... pointer to array for storing control points             */
     /*    ui    ... u-values of interpolation                               */ 
     /*    order ... order of interpolation polynomial                       */
     /*                                                                      */
     /* return:                                                              */
     /*    u-values of control points in the array utest                     */
     /*----------------------------------------------------------------------*/
{
  int k,j,i;
  double sum, qsum,x;
  
  /* compute approximate maxima of error polynomial */
  for(k=0; k<order; k++) {

    /* check for multiple points */
    if ( (k==0 && _unur_iszero(ui[0])) ||
    	 (k>0  && _unur_FP_same(ui[k-1],ui[k])) ) {
      utest[k] = ui[k]; 
      continue;
    }

    /* first approximation: use mean of consecuting construction points */
    x = (k>0) ? 0.5*(ui[k-1]+ui[k]) : 0.5*ui[k];
    
    /* make two iterations for root finding */
    for(j=1; j<=2; j++) {
      sum = 1./x;
      qsum = sum*sum;
      for(i=0; i<order; i++){
	sum += 1./(x-ui[i]);
	qsum += 1./((x-ui[i])*(x-ui[i]));
      }
      x += sum/qsum;
    }

    /* store result in utest */
    utest[k] = x;
  }
  
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_testpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_linear_testpoints (double *utest, double *ui, int order)
     /*----------------------------------------------------------------------*/
     /* [2c.] create table of test points for linear interpolation.          */
     /*                                                                      */
     /* parameters:                                                          */
     /*    utest ... pointer to array for storing control points             */
     /*    ui    ... u-values of interpolation                               */ 
     /*    order ... order of interpolation polynomial                       */
     /*                                                                      */
     /* return:                                                              */
     /*    u-values of control points in the array utest                     */
     /*----------------------------------------------------------------------*/
{
  int k;
  double dx = ui[order-1] / order;  /* distance between points */

  /* use equidistributed points */
  for(k=0; k<order; k++)
    utest[k] = (k+0.5) * dx;

  return UNUR_SUCCESS;
} /* end of _unur_pinv_linear_testpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_cubic_hermite_is_monotone(struct unur_gen *gen, double *ui, double *zi, double *xval)
     /*----------------------------------------------------------------------*/
     /* [2c.] check monotonicity of cubic Hermite interpolation              */
     /*       (order = 3, smoothness = 1)                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer generator object                                  */
     /*   ui   ... u-values of interpolation                                 */ 
     /*   zi   ... coefficients of polynomial (increasing order)             */
     /*   xval ... x-values for constructing polynomial                      */
     /*                                                                      */
     /* Cubic Hermite interpolation is monotonically increasing if           */
     /*   1/f'(x0) <= 3 * (x1-x0) / (u1 - u0)   and                          */
     /*   1/f'(x1) <= 3 * (x1-x0) / (u1 - u0)                                */
     /* where (u0,u1) is the interval for the polynomial and                 */
     /* xi = F^{-1}(ui) are the values of the inverse CDF at the boundary    */
     /* and f denotes the PDF.                                               */
     /*                                                                      */
     /* Notice that u0=0, u1=ui[2], x0=xval[0], and x1=xval[2].              */
     /* 1/f'(x0) = zi[0]                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if cubic is monotone                                     */
     /*   FALSE ... otherwise                                                */
     /*----------------------------------------------------------------------*/
{
  double f0,f1,dq;

  /* we skip the test if computing the bound has too many round-off errors */
  if (_unur_iszero(ui[2])) return TRUE;

  /* difference quotient */
  dq = (xval[2] - xval[0]) / ui[2];

  /* derivative of inverse CDF at boundary points */
  f0 = zi[0];
  f1 = 1./PDF(xval[2]);

  /* test condition */
  return (f0 <= 3. * dq && f1 <= 3. * dq) ? TRUE : FALSE;

} /* end of _unur_pinv_cubic_hermite_is_monotone() */

/*---------------------------------------------------------------------------*/

int 
_unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx )
     /*----------------------------------------------------------------------*/
     /* make a new interval i with left boundary point x and CDF(x).         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   i    ... index (number) of interval                                */
     /*   x    ... left boundary point of new interval                       */
     /*   cdfx ... CDF at x                                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
     /* We use Newtons interpolating polynomial. The construction points are */
     /* created by using Chebyshev points for the x-coordinates.             */
     /*                                                                      */
     /* Alternative implementation would be:                                 */
     /* use Chebyshev points for u-coordinates.                              */
     /* advantage:                                                           */
     /*  - smaller tables as u-coordinates need not be stored.               */
     /*  - fewer intervals required (especially in the tails).               */
     /*  - routine _unur_pinv_newton_testpoints() not required               */
     /*    (test points are fixed Chebyshev points -- closed form).          */ 
     /* disadvantage:                                                        */
     /*  - numerically less stable and more expensive                        */
     /*    as it requires explicit inversion.                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_pinv_interval *iv;

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);

  /* check for free intervals */
  if (i >= GEN->max_ivs) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"maximum number of intervals exceeded");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* set values */
  iv = GEN->iv+i;     /* pointer to interval */
  iv->xi = x;         /* left boundary of interval */
  iv->cdfi = cdfx;    /* CDF at left boundary */
  COOKIE_SET(iv,CK_PINV_IV);

  /* allocate space for coefficients for Newton interpolation */
  iv->ui = _unur_xmalloc( GEN->order * sizeof(double) );
  iv->zi = _unur_xmalloc( GEN->order * sizeof(double) );

  /* update size of array (number of intervals) */
  GEN->n_ivs = i;

  /* set bookmark in table of integral values */
  _unur_lobatto_find_linear(GEN->aCDF,x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_interval() */

/*---------------------------------------------------------------------------*/

int 
_unur_pinv_lastinterval( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* update size of array and set all uninitialized values to 0.          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double *ui, *zi;
  int i;

  /* pointer to last interval */
  struct unur_pinv_interval *last_iv = GEN->iv + GEN->n_ivs;
  ui = last_iv->ui;
  zi = last_iv->zi;
  
  /* update size of array */
  GEN->iv = _unur_xrealloc( GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );

  /* set all entries ui and zi in very last interval to 0 */
  for (i=0; i<GEN->order; i++) {
    ui[i] = 0.;
    zi[i] = 0.;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_lastinterval() */

/*---------------------------------------------------------------------------*/
