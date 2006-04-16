/* tnc : truncated newton bound contrained minimization
         using gradient information, in C */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C implementation of TNBC, a truncated newton minimization
 * package originally developed by Stephen G. Nash in Fortran.
 * 
 * The original source code can be found at :
 * http://iris.gmu.edu/~snash/nash/software/software.html
 * 
 * Copyright for the original TNBC fortran routines:
 * 
 *   TRUNCATED-NEWTON METHOD:  SUBROUTINES
 *     WRITTEN BY:  STEPHEN G. NASH
 *           SCHOOL OF INFORMATION TECHNOLOGY & ENGINEERING
 *           GEORGE MASON UNIVERSITY
 *           FAIRFAX, VA 22030
 */

/*
 * Conversion into C by Elisabeth Nguyen & Jean-Sebastien Roy
 * Modifications by Jean-Sebastien Roy, 2001-2002
 */

static char const rcsid[] =
  "@(#) $Jeannot: tnc.c,v 1.201 2004/04/02 22:36:25 js Exp $";

static char const copyright[] =
  "(c) 2002-2003, Jean-Sebastien Roy (js@jeannot.org)";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tnc.h"

typedef enum
{
  TNC_FALSE = 0,
  TNC_TRUE
} logical;

/*
 * Return code strings
 */

char *tnc_rc_string[10] =
{
  "Memory allocation failed",
  "Invalid parameters (less than 1 dimension)",
  "Infeasible (low bound > up bound)",
  "Local minima reach (|pg| ~= 0)",
  "Converged (|f_n-f_(n-1)| ~= 0)",
  "Maximum number of function evaluations reached",
  "Linear search failed",
  "All lower bounds are equal to the upper bounds",
  "Unable to progress",
  "User requested end of minimization"
};

/*
 * getptc return codes
 */
typedef enum
{
  GETPTC_OK     = 0, /* Suitable point found */
  GETPTC_EVAL   = 1, /* Function evaluation required */
  GETPTC_EINVAL = 2, /* Bad input values */
  GETPTC_FAIL   = 3  /* No suitable point found */
} getptc_rc;

/*
 * linearSearch return codes
 */
typedef enum
{
  LS_OK        = 0, /* Suitable point found */
  LS_MAXFUN    = 1, /* Max. number of function evaluations reach */
  LS_FAIL      = 2, /* No suitable point found */
  LS_USERABORT = 3, /* User requested end of minimization */
  LS_ENOMEM    = 4  /* Memory allocation failed */
} ls_rc;

/*
 * Prototypes
 */
static tnc_rc tnc_minimize(int n, double x[], double *f, double g[],
  tnc_function *function, void *state, 
  double xscale[], double *fscale,
  double low[], double up[], tnc_message messages,
  int maxCGit, int maxnfeval, int *nfeval,
  double eta, double stepmx, double accuracy,
  double fmin, double ftol, double rescale);

static getptc_rc getptcInit(double *reltol, double *abstol, double tnytol, 
  double eta, double rmu, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol);

static getptc_rc getptcIter(double big, double 
  rtsmll, double *reltol, double *abstol, double tnytol, 
  double fpresn, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol);

static void printCurrentIteration(int n, double f, double g[], int niter,
  int nfeval, int pivot[]);

static double initialStep(double fnew, double fmin, double gtp, double smax);

static ls_rc linearSearch(int n, tnc_function *function, void *state,
  double low[], double up[],
  double xscale[], double fscale, int pivot[],
  double eta, double ftol, double xbnd, 
  double p[], double x[], double *f, 
  double *alpha, double gfull[], int maxnfeval, int *nfeval);

static int tnc_direction(double *zsol, double *diagb,
  double *x, double *g, int n,
  int maxCGit, int maxnfeval, int *nfeval, 
  logical upd1, double yksk, double yrsr,
  double *sk, double *yk, double *sr, double *yr,
  logical lreset, tnc_function *function, void *state, 
  double xscale[], double fscale,
  int *pivot, double accuracy,
  double gnorm, double xnorm, double *low, double *up);

static double stepMax(double step, int n, double x[], double p[], int pivot[], 
  double low[], double up[], double xscale[]);

/* Active set of constraints */
static void setContraints(int n, double x[], int pivot[], double xscale[],
  double low[], double up[]);

static logical addConstraint(int n, double x[], double p[], int pivot[],
  double low[], double up[], double xscale[]);

static logical removeConstraint(double gtpnew, double f, 
  double *fLastConstraint, double g[], int pivot[], int n);

static void project(int n, double x[], int pivot[]);

static int hessianTimesVector(double v[], double gv[], int n, 
  double x[], double g[], tnc_function *function, void *state, 
  double xscale[], double fscale,
  double accuracy, double xnorm, double low[], double up[]);

static int msolve(double g[], double *y, int n, 
  double sk[], double yk[], double diagb[], double sr[], 
  double yr[], logical upd1, double yksk, double yrsr, 
  logical lreset);

static void diagonalScaling(int n, double e[], double v[], double gv[],
  double r[]);

static void ssbfgs(int n, double gamma, double sj[], double *hjv,
  double hjyj[], double yjsj, 
  double yjhyj, double vsj, double vhyj, double hjp1v[]);

static int initPreconditioner(double diagb[], double emat[], int n, 
  logical lreset, double yksk, double yrsr,
  double sk[], double yk[], double sr[], double yr[], 
  logical upd1);

/* Scaling */
static void coercex(int n, double x[], double low[], double up[]);
static void unscalex(int n, double x[], double xscale[]);
static void scaleg(int n, double g[], double xscale[], double fscale);
static void scalex(int n, double x[], double xscale[]);
static void projectConstants(int n, double x[], double xscale[]);

/* Machine precision */
static double mchpr1(void);

/* Special blas for incx=incy=1 */
static double ddot1(int n, double dx[], double dy[]);
static void dxpy1(int n, double dx[], double dy[]);
static void daxpy1(int n, double da, double dx[], double dy[]);
static void dcopy1(int n, double dx[], double dy[]);
static double dnrm21(int n, double dx[]);

/* additionnal blas-like functions */
static void dneg1(int n, double v[]);
static double dnrmi1(int n, double v[]);

/*
 * This routine solves the optimization problem
 * 
 *   minimize   f(x)
 *     x
 *   subject to   low <= x <= up
 * 
 * where x is a vector of n real variables. The method used is
 * a truncated-newton algorithm (see "newton-type minimization via
 * the lanczos algorithm" by s.g. nash (technical report 378, math.
 * the lanczos method" by s.g. nash (siam j. numer. anal. 21 (1984),
 * pp. 770-778).  this algorithm finds a local minimum of f(x). It does
 * not assume that the function f is convex (and so cannot guarantee a
 * global solution), but does assume that the function is bounded below.
 * it can solve problems having any number of variables, but it is
 * especially useful when the number of variables (n) is large.
 * 
 */
extern int tnc(int n, double x[], double *f, double g[], tnc_function *function,
  void *state, double low[], double up[], double scale[], int messages, 
  int maxCGit, int maxnfeval, double eta, double stepmx, 
  double accuracy, double fmin, double ftol, double rescale, int *nfeval)
{
  int rc, frc, i, nc, nfeval_local,
    free_low = TNC_FALSE, free_up = TNC_FALSE,
    free_g = TNC_FALSE;
  double *xscale = NULL, fscale, epsmch, rteps;

  if(nfeval==NULL)
  {
    /* Ignore nfeval */
    nfeval = &nfeval_local;
  }
  *nfeval = 0;
  
  /* Version info */
  if (messages & TNC_MSG_VERS)
  {
    fprintf(stderr, "tnc: Version %s, %s\n",TNC_VERSION,copyright);
    fprintf(stderr, "tnc: RCS ID: %s\n",rcsid);
  }

  /* Check for errors in the input parameters */
  if (n < 1)
  {
    rc = TNC_EINVAL;
    goto cleanup;
  }
  
  /* Check bounds arrays */
  if (low == NULL)
  {
    low = malloc(n*sizeof(*low));
    if (low == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_low = TNC_TRUE;
    for (i = 0 ; i < n ; i++) low[i] = -HUGE_VAL;
  }
  if (up == NULL)
  {
    up = malloc(n*sizeof(*up));
    if (up == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_up = TNC_TRUE;
    for (i = 0 ; i < n ; i++) up[i] = HUGE_VAL;
  }

  /* Coherency check */
  for (i = 0 ; i < n ; i++)
  {
    if (low[i] > up [i])
    {
      rc = TNC_INFEASIBLE;
      goto cleanup;
    }
  }

  /* Coerce x into bounds */
  coercex(n, x, low, up);

  if (maxnfeval < 1)
  {
    rc = TNC_MAXFUN;
    goto cleanup;
  }
  
  /* Allocate g if necessary */
  if(g == NULL)
  {
    g = malloc(n*sizeof(*g));
    if (g == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_g = TNC_TRUE;
  }

  /* Initial function evaluation */  
  frc = function(x, f, g, state);
  (*nfeval) ++;
  if (frc)
  {
    rc = TNC_USERABORT;
    goto cleanup;
  }
    
  /* Constant problem ? */
  for (nc = 0, i = 0 ; i < n ; i++)
    if ((low[i] == up[i]) || (scale != NULL && scale[i] == 0.0))
      nc ++;

  if (nc == n)
  {
    rc = TNC_CONSTANT;
    goto cleanup;
  }

  /* Scaling parameters */  
  xscale = malloc(sizeof(*xscale)*n);
  if (xscale == NULL)
  {
    rc = TNC_ENOMEM;
    goto cleanup;
  }
  fscale = 1.0;

  for (i = 0 ; i < n ; i++)
  {
    if (scale != NULL)
    {
      xscale[i] = fabs(scale[i]);
      if (xscale[i] == 0.0)
        low[i] = up[i] = x[i];
    }
    else if (low[i] != -HUGE_VAL && up[i] != HUGE_VAL)
      xscale[i] = up[i] - low[i];
    else
      xscale[i] = 1.0+fabs(x[i]);
  }

  /* Default values for parameters */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);

  if (stepmx < rteps * 10.0) stepmx = 1.0e1;
  if (eta < 0.0 || eta >= 1.0) eta = 0.25;
  if (rescale < 0) rescale = 1.3;
  if (maxCGit < 0) /* maxCGit == 0 is valid */
  {
    maxCGit = n / 2;
    if (maxCGit < 1) maxCGit = 1;
    else if (maxCGit > 50) maxCGit = 50;
  }
  if (maxCGit > n) maxCGit = n;
  if (ftol < 0.0) ftol = 0.0;
  if (accuracy <= epsmch) accuracy = rteps;

  /* Optimisation */
  rc = tnc_minimize(n, x, f, g, function, state,
    xscale, &fscale, low, up, messages, 
    maxCGit, maxnfeval, nfeval, eta, stepmx, accuracy, fmin, ftol, rescale);

cleanup:
  if (messages & TNC_MSG_EXIT)
    fprintf(stderr, "tnc: %s\n", tnc_rc_string[rc - TNC_MINRC]);

  if (xscale) free(xscale);
  if (free_low) free(low);
  if (free_up) free(up);
  if (free_g) free(g);
  
  return rc;
}

/* Coerce x into bounds */
static void coercex(int n, double x[], double low[], double up[])
{
  int i;
  
  for (i = 0 ; i < n ; i++)
  {
    if (x[i]<low[i]) x[i] = low[i];
    else if (x[i]>up[i]) x[i] = up[i];
  }
}

/* Unscale x */
static void unscalex(int n, double x[], double xscale[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    x[i] *= xscale[i];
}

/* Scale x */
static void scalex(int n, double x[], double xscale[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    if (xscale[i]>0.0)
      x[i] /= xscale[i];
}

/* Scale g */
static void scaleg(int n, double g[], double xscale[], double fscale)
{
  int i;
  for (i = 0 ; i < n ; i++)
    g[i] *= xscale[i]*fscale;
}

/* Caculate the pivot vector */
static void setContraints(int n, double x[], int pivot[], double xscale[],
  double low[], double up[])
{
  int i;
  double epsmch;
  
  epsmch = mchpr1();

  for (i = 0; i < n; i++)
  {
    double tol;

    /* tolerances should be better ajusted */
    if (xscale[i] == 0.0)
    {
      pivot[i] = 2;
    }
    else
    {
       tol = epsmch * 10.0 * (fabs(low[i]) + 1.0);
      if ((x[i]*xscale[i] - low[i] <= tol) && low[i] != - HUGE_VAL)
        pivot[i] = -1;
      else
      {
         tol = epsmch * 10.0 * (fabs(up[i]) + 1.0);
        if ((x[i]*xscale[i] - up[i] >= tol) && up[i] != HUGE_VAL)
          pivot[i] = 1;
        else
          pivot[i] = 0;
      }
    }
  }
}

/*
 * This routine is a bounds-constrained truncated-newton method.
 * the truncated-newton method is preconditioned by a limited-memory
 * quasi-newton method (this preconditioning strategy is developed
 * in this routine) with a further diagonal scaling
 * (see routine diagonalscaling).
 */
static tnc_rc tnc_minimize(int n, double x[], 
  double *f, double gfull[], tnc_function *function, void *state, 
  double xscale[], double *fscale,
  double low[], double up[], tnc_message messages, 
  int maxCGit, int maxnfeval, int *nfeval, double eta, double stepmx, 
  double accuracy, double fmin, double ftol, double rescale)
{
  double fLastReset, difnew, epsmch, epsred, oldgtp,
    difold, oldf, rteps, xnorm, newscale,
    gnorm, ustpmax, fLastConstraint, spe, yrsr, yksk,
    *temp = NULL, *sk = NULL, *yk = NULL, *diagb = NULL, *sr = NULL,
    *yr = NULL, *oldg = NULL, *pk = NULL, *g = NULL;
  double alpha = 0.0; /* Default unused value */
  int i, icycle, niter = 0, oldnfeval, *pivot = NULL, frc;
  logical lreset, newcon, upd1, remcon;
  tnc_rc rc = TNC_ENOMEM; /* Default error */

  /* Allocate temporary vectors */  
  oldg = malloc(sizeof(*oldg)*n);
   if (oldg == NULL) goto cleanup;
  g = malloc(sizeof(*g)*n);
   if (g == NULL) goto cleanup;
  temp = malloc(sizeof(*temp)*n);
   if (temp == NULL) goto cleanup;
  diagb = malloc(sizeof(*diagb)*n);
   if (diagb == NULL) goto cleanup;
  pk = malloc(sizeof(*pk)*n);
   if (pk == NULL) goto cleanup;

  sk = malloc(sizeof(*sk)*n);
   if (sk == NULL) goto cleanup;
  yk = malloc(sizeof(*yk)*n);
   if (yk == NULL) goto cleanup;
  sr = malloc(sizeof(*sr)*n);
   if (sr == NULL) goto cleanup;
  yr = malloc(sizeof(*yr)*n);
   if (yr == NULL) goto cleanup;

  pivot = malloc(sizeof(*pivot)*n);
   if (pivot == NULL) goto cleanup;

  /* Initialize variables */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);

  difnew = 0.0;
  epsred = 0.05;
  upd1 = TNC_TRUE;
  icycle = n - 1;
  newcon = TNC_TRUE;
  
  /* Uneeded initialisations */
  lreset = TNC_FALSE;
  yrsr = 0.0;
  yksk = 0.0;
  
  /* Initial scaling */
  scalex(n, x, xscale);
  (*f) *= *fscale;

  /* initial pivot calculation */
  setContraints(n, x, pivot, xscale, low, up);

  dcopy1(n, gfull, g);
  scaleg(n, g, xscale, *fscale);

  /* Test the lagrange multipliers to see if they are non-negative. */
  for (i = 0; i < n; i++)
    if (-pivot[i] * g[i] < 0.0)
      pivot[i] = 0;

  project(n, g, pivot);

  /* Set initial values to other parameters */
  gnorm = dnrm21(n, g);

  fLastConstraint = *f; /* Value at last constraint */
  fLastReset = *f; /* Value at last reset */

  if (messages & TNC_MSG_ITER) fprintf(stderr,
    "  NIT   NF   F                       GTG\n");
  if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
    niter, *nfeval, pivot);

  /* Set the diagonal of the approximate hessian to unity. */
  for (i = 0; i < n; i++) diagb[i] = 1.0;

  /* Start of main iterative loop */
  while(TNC_TRUE)
  {
    /* Tolerance should be user modifiable */
    if (dnrmi1(n, g) <= 1.0e-2*rteps*fabs(*f))
    {
      /* |PG| == 0.0 => local minimum */
      dcopy1(n, gfull, g);
      project(n, g, pivot);
      if (messages & TNC_MSG_INFO) fprintf(stderr,
        "tnc: |pg| = %g -> local minimum\n",dnrmi1(n, g));
      rc = TNC_LOCALMINIMUM;
      break;
    }

    /* Terminate if more than maxnfeval evaluations have been made */
    if (*nfeval >= maxnfeval)
    {
      rc = TNC_MAXFUN;
      break;
    }

    /* Rescale function if necessary */
    newscale = dnrmi1(n, g);
    if ((newscale > epsmch) && (fabs(log10(newscale)) > rescale))
    {
      newscale = 1.0/newscale;
      
      *f *= newscale;
      *fscale *= newscale;
      gnorm *= newscale;
      fLastConstraint *= newscale;
      fLastReset *= newscale;
      difnew *= newscale;

      for (i = 0; i < n; i++) g[i] *= newscale;
      for (i = 0; i < n; i++) diagb[i] = 1.0;

      upd1 = TNC_TRUE;
      icycle = n - 1;
      newcon = TNC_TRUE;

      if (messages & TNC_MSG_INFO) fprintf(stderr, 
        "tnc: fscale = %g\n", *fscale);
    }

    dcopy1(n, x, temp);
    project(n, temp, pivot);
    xnorm = dnrm21(n, temp);
    oldnfeval = *nfeval;

    /* Compute the new search direction */
    frc = tnc_direction(pk, diagb, x, g, n, maxCGit, maxnfeval, nfeval,
      upd1, yksk, yrsr, sk, yk, sr, yr,
      lreset, function, state, xscale, *fscale,
      pivot, accuracy, gnorm, xnorm, low, up);

    if (frc == -1)
    {
      rc = TNC_ENOMEM;
      break;
    }

    if (frc)
    {
      rc = TNC_USERABORT;
      break;
    }

    if (!newcon)
    {
      if (!lreset)
      {
        /* Compute the accumulated step and its corresponding gradient
          difference. */
        dxpy1(n, sk, sr);
        dxpy1(n, yk, yr);
        icycle++;
      }
      else
      {
        /* Initialize the sum of all the changes */
        dcopy1(n, sk, sr);
        dcopy1(n, yk, yr);
        fLastReset = *f;
        icycle = 1;
      }
    }

    dcopy1(n, g, oldg);
    oldf = *f;
    oldgtp = ddot1(n, pk, g);

    /* Maximum unconstrained step length */
    ustpmax = stepmx / (dnrm21(n, pk) + epsmch);

    /* Maximum constrained step length */
    spe = stepMax(ustpmax, n, x, pk, pivot, low, up, xscale);
    
    if (spe > 0.0)
    {
      ls_rc lsrc;
      /* Set the initial step length */
      alpha = initialStep(*f, fmin / (*fscale), oldgtp, spe);

      /* Perform the linear search */
      lsrc = linearSearch(n, function, state, low, up,
        xscale, *fscale, pivot,
        eta, ftol, spe, pk, x, f, &alpha, gfull, maxnfeval, nfeval);

      if (lsrc == LS_ENOMEM)
      {
        rc = TNC_ENOMEM;
        break;
      }

      if (lsrc == LS_USERABORT)
      {
        rc = TNC_USERABORT;
        break;
      }

      /* If we went up to the maximum unconstrained step, increase it */
      if (alpha >= 0.9 * ustpmax)
      {
        stepmx *= 1e2;
        if (messages & TNC_MSG_INFO) fprintf(stderr,
          "tnc: stepmx = %g\n", stepmx);
      }

      /* If we went up to the maximum constrained step,
         a new constraint was encountered */
      if (alpha - spe >= -epsmch * 10.0)
      {
        newcon = TNC_TRUE;
      }
      else
      {
        /* Break if the linear search has failed to find a lower point */
        if (lsrc != LS_OK)
        {
          if (lsrc == LS_MAXFUN) rc = TNC_MAXFUN;
          else rc = TNC_LSFAIL;
          break;
        }
        newcon = TNC_FALSE;
      }
    }
    else
    {
      /* Maximum constrained step == 0.0 => new constraint */
      newcon = TNC_TRUE;
    }

    if (newcon)
    {
      if(!addConstraint(n, x, pk, pivot, low, up, xscale))
      {
        if(*nfeval == oldnfeval)
        {
          rc = TNC_NOPROGRESS;
          break;
        }
      }
      fLastConstraint = *f;
    }

    niter++;

    /* Set up parameters used in convergence and resetting tests */
    difold = difnew;
    difnew = oldf - *f;

    /* If this is the first iteration of a new cycle, compute the
       percentage reduction factor for the resetting test */
    if (icycle == 1)
    {
      if (difnew > difold * 2.0) epsred += epsred;
      if (difnew < difold * 0.5) epsred *= 0.5;
    }

    dcopy1(n, gfull, g);
    scaleg(n, g, xscale, *fscale);

    dcopy1(n, g, temp);
    project(n, temp, pivot);
    gnorm = dnrm21(n, temp);

    /* Reset pivot */
    remcon = removeConstraint(oldgtp, *f, &fLastConstraint, g, pivot, n);

    if (!remcon && !newcon)
    {
      /* No constraint removed & no new constraint : test for convergence */
      if (fabs(difnew) <= ftol*epsmch*0.5*(fabs(oldf)+fabs(*f)))
      {
        if (messages & TNC_MSG_INFO) fprintf(stderr, 
          "tnc: |fn-fn-1] = %g -> convergence\n",fabs(difnew));
        rc = TNC_CONVERGED;
        break;
      }
    }

    project(n, g, pivot);

    if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
      niter, *nfeval, pivot);

    /* Compute the change in the iterates and the corresponding change in the
      gradients */
    if (!newcon)
    {
      for (i = 0; i < n; i++)
      {
        yk[i] = g[i] - oldg[i];
        sk[i] = alpha * pk[i];
      }

      /* Set up parameters used in updating the preconditioning strategy */
      yksk = ddot1(n, yk, sk);
      
      if (icycle == (n - 1) || difnew < epsred * (fLastReset - *f))
        lreset = TNC_TRUE;
      else
      {
        yrsr = ddot1(n, yr, sr);
        if (yrsr <= 0.0) lreset = TNC_TRUE;
        else lreset = TNC_FALSE;
      }
      upd1 = TNC_FALSE;
    }
  }

  if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
    niter, *nfeval, pivot);

  /* Unscaling */
  unscalex(n, x, xscale);
  coercex(n, x, low, up);
  (*f) /= *fscale;

cleanup: 
  if (oldg) free(oldg);
  if (g) free(g);
  if (temp) free(temp);
  if (diagb) free(diagb);
  if (pk) free(pk);

  if (sk) free(sk);
  if (yk) free(yk);
  if (sr) free(sr);
  if (yr) free(yr);
  
  if (pivot) free(pivot);

  return rc;
}

/* Print the results of the current iteration */
static void printCurrentIteration(int n, double f, double g[], int niter,
  int nfeval, int pivot[])
{
  int i;
  double gtg;

  gtg = 0.0;
  for (i = 0; i < n; i++)
    if (pivot[i] == 0)
      gtg += g[i] * g[i];

  fprintf(stderr, " %4d %4d %22.15E  %15.8E\n", niter, nfeval, f, gtg);
}

/*
 * Set x[i] = 0.0 if direction i is currently constrained 
 */
static void project(int n, double x[], int pivot[])
{
  int i;
  for (i = 0; i < n; i++)
    if (pivot[i] != 0)
      x[i] = 0.0;
}

/*
 * Set x[i] = 0.0 if direction i is constant
 */
static void projectConstants(int n, double x[], double xscale[])
{
  int i;
  for (i = 0; i < n; i++)
    if (xscale[i] == 0.0)
      x[i] = 0.0;
}

/*
 * Compute the maximum allowable step length
 */
static double stepMax(double step, int n, double x[], double dir[],
  int pivot[], double low[], double up[], double xscale[])
{
  int i;
  double t;

  /* Constrained maximum step */
  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (dir[i] != 0.0))
    {
      if (dir[i] < 0.0)
      {
        t = low[i]/xscale[i] - x[i];
        if (t > step * dir[i]) step = t / dir[i];
      }
      else
      {
        t = up[i]/xscale[i] - x[i];
        if (t < step * dir[i]) step = t / dir[i];
      }
    }
  }
  
  return step;
}

/*
 * Update the constraint vector pivot if a new constraint is encountered
 */
static logical addConstraint(int n, double x[], double p[], int pivot[],
  double low[], double up[], double xscale[])
{
  int i, newcon = TNC_FALSE;
  double tol, epsmch;

  epsmch = mchpr1();

  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (p[i] != 0.0))
    {
       if (p[i] < 0.0 && low[i] != - HUGE_VAL)
      {
         tol = epsmch * 10.0 * (fabs(low[i]) + 1.0);
        if (x[i]*xscale[i] - low[i] <= tol)
        {
          pivot[i] = -1;
          x[i] = low[i]/xscale[i];
          newcon = TNC_TRUE;
        }
      }
      else if (up[i] != HUGE_VAL)
      {
        tol = epsmch * 10.0 * (fabs(up[i]) + 1.0);
        if (up[i] - x[i]*xscale[i] <= tol)
        {
          pivot[i] = 1;
          x[i] = up[i]/xscale[i];
          newcon = TNC_TRUE;
        }
      }
    }
  }
  return newcon;
}

/*
 * Check if a constraint is no more active
 */
static logical removeConstraint(double gtpnew, double f, 
  double *fLastConstraint, double g[], int pivot[], int n)
{
  double cmax, t;
  int imax, i;
  logical ltest;

  imax = -1;
  cmax = 0.0;
  ltest = (*fLastConstraint - f) <= (gtpnew * -0.5);
  for (i = 0; i < n; i++)
  {
    if (pivot[i] != 2)
    {
      t = -pivot[i] * g[i];
      if (t < 0.0)
      {
        if ((!ltest) && (cmax > t))
        {
          cmax = t;
          imax = i;
        }
      }
    }
  }

  if (imax != -1)
  {
    pivot[imax] = 0;
    *fLastConstraint = f;
    return TNC_TRUE;
  }
  else
    return TNC_FALSE;

/*
 * For details, see gill, murray, and wright (1981, p. 308) and
 * fletcher (1981, p. 116). The multiplier tests (here, testing
 * the sign of the components of the gradient) may still need to
 * modified to incorporate tolerances for zero.
 */
}

/*
 * This routine performs a preconditioned conjugate-gradient
 * iteration in order to solve the newton equations for a search
 * direction for a truncated-newton algorithm.
 * When the value of the quadratic model is sufficiently reduced,
 * the iteration is terminated.
 */
static int tnc_direction(double *zsol, double *diagb,
  double *x, double g[], int n,
  int maxCGit, int maxnfeval, int *nfeval, 
  logical upd1, double yksk, double yrsr,
  double *sk, double *yk, double *sr, double *yr,
  logical lreset, tnc_function *function, void *state, 
  double xscale[], double fscale,
  int *pivot, double accuracy,
  double gnorm, double xnorm, double low[], double up[])
{
  double alpha, beta, qold, qnew, rhsnrm, tol, vgv, rz, rzold, qtest, pr, gtp;
  int i, k, frc;
  /* Temporary vectors */
  double *r = NULL, *zk = NULL, *v = NULL, *emat = NULL, *gv = NULL;

  /* No CG it. => dir = -grad */
  if (maxCGit == 0)
  {
    dcopy1(n, g, zsol);
    dneg1(n, zsol);
    project(n, zsol, pivot);
    return 0;
  }

  /* General initialization */
  rhsnrm = gnorm;
  tol = 1e-12;
  qold = 0.0;
  rzold = 0.0; /* Uneeded */

  frc = -1; /* ENOMEM here */
  r = malloc(sizeof(*r)*n); /* Residual */
  if (r == NULL) goto cleanup;
  v = malloc(sizeof(*v)*n);
  if (v == NULL) goto cleanup;
  zk = malloc(sizeof(*zk)*n);
  if (zk == NULL) goto cleanup;
  emat = malloc(sizeof(*emat)*n); /* Diagonal preconditoning matrix */
  if (emat == NULL) goto cleanup;
  gv = malloc(sizeof(*gv)*n); /* hessian times v */
  if (gv == NULL) goto cleanup;

  /* Initialization for preconditioned conjugate-gradient algorithm */
  frc = initPreconditioner(diagb, emat, n, lreset, yksk, yrsr, sk, yk, sr, yr,
    upd1);
  if (frc) goto cleanup;

  for (i = 0; i < n; i++)
  {
    r[i] = -g[i];
    v[i] = 0.0;
    zsol[i] = 0.0; /* Computed search direction */
  }

  /* Main iteration */
  for (k = 0; k < maxCGit; k++)
  {
    /* CG iteration to solve system of equations */
    project(n, r, pivot);
    frc = msolve(r, zk, n, sk, yk, diagb, sr, yr, upd1, yksk, yrsr, lreset);
    if (frc) goto cleanup;
    project(n, zk, pivot);
    rz = ddot1(n, r, zk);

    if ((rz / rhsnrm < tol) || ((*nfeval) >= (maxnfeval-1)))
    {
      /* Truncate algorithm in case of an emergency
         or too many function evaluations */
      if (k == 0)
      {
        dcopy1(n, g, zsol);
        dneg1(n, zsol);
        project(n, zsol, pivot);
      }
      break;
    }
    if (k == 0) beta = 0.0;
    else beta = rz / rzold;

    for (i = 0; i < n; i++)
      v[i] = zk[i] + beta * v[i];

    project(n, v, pivot);
    frc = hessianTimesVector(v, gv, n, x, g, function, state,
      xscale, fscale, accuracy, xnorm, low, up);
    ++(*nfeval);
    if (frc) goto cleanup;
    project(n, gv, pivot);

    vgv = ddot1(n, v, gv);
    if (vgv / rhsnrm < tol)
    {
      /* Truncate algorithm in case of an emergency */
      if (k == 0)
      {
        frc = msolve(g, zsol, n, sk, yk, diagb, sr, yr, upd1, yksk, yrsr,
          lreset);
        if (frc) goto cleanup;
        dneg1(n, zsol);
        project(n, zsol, pivot);
      }
      break;
    }
    diagonalScaling(n, emat, v, gv, r);

    /* Compute linear step length */
    alpha = rz / vgv;

    /* Compute current solution and related vectors */
    daxpy1(n, alpha, v, zsol);
    daxpy1(n, -alpha, gv, r);

    /* Test for convergence */
    gtp = ddot1(n, zsol, g);
    pr = ddot1(n, r, zsol);
    qnew = (gtp + pr) * 0.5;
    qtest = (k + 1) * (1.0 - qold / qnew);
    if (qtest <= 0.5) break;

    /* Perform cautionary test */
    if (gtp > 0.0)
    {
      /* Truncate algorithm in case of an emergency */
      daxpy1(n, -alpha, v, zsol);
      break;
    }

    qold = qnew;
    rzold = rz;
  }

  /* Terminate algorithm */
  /* Store (or restore) diagonal preconditioning */
  dcopy1(n, emat, diagb);

cleanup:
  if (r) free(r);
  if (v) free(v);
  if (zk) free(zk);
  if (emat) free(emat);
  if (gv) free(gv);
  return frc;
}

/* 
 * Update the preconditioning matrix based on a diagonal version
 * of the bfgs quasi-newton update.
 */
static void diagonalScaling(int n, double e[], double v[], double gv[],
  double r[])
{
  int i;
  double vr, vgv;

  vr = 1.0/ddot1(n, v, r);
  vgv = 1.0/ddot1(n, v, gv);
  for (i = 0; i < n; i++)
  {
    e[i] += - r[i]*r[i]*vr + gv[i]*gv[i]*vgv;
    if (e[i] <= 1e-6) e[i] = 1.0;
  }
}

/*
 * Returns the length of the initial step to be taken along the
 * vector p in the next linear search.
 */
static double initialStep(double fnew, double fmin, double gtp, double smax)
{
  double d, alpha;

  d = fabs(fnew - fmin);
  alpha = 1.0;
  if (d * 2.0 <= -(gtp) && d >= mchpr1()) alpha = d * -2.0 / gtp;
  if (alpha >= smax) alpha = smax;

  return alpha;
}

/*
 * Hessian vector product through finite differences
 */
static int hessianTimesVector(double v[], double gv[], int n, 
  double x[], double g[], tnc_function *function, void *state, 
  double xscale[], double fscale,
  double accuracy, double xnorm, double low[], double up[])
{
  double dinv, f, delta, *xv;
  int i, frc;
  
  xv = malloc(sizeof(*xv)*n);
  if (xv == NULL) return -1;

  delta = accuracy * (xnorm + 1.0);
  for (i = 0; i < n; i++)
    xv[i] = x[i] + delta * v[i];

  unscalex(n, xv, xscale);
  coercex(n, xv, low, up);
  frc = function(xv, &f, gv, state);
  free(xv);
  if (frc) return 1;
  scaleg(n, gv, xscale, fscale);

  dinv = 1.0 / delta;
  for (i = 0; i < n; i++)
    gv[i] = (gv[i] - g[i]) * dinv;
    
  projectConstants(n, gv, xscale);

  return 0;
}

/*
 * This routine acts as a preconditioning step for the 
 * linear conjugate-gradient routine. It is also the 
 * method of computing the search direction from the 
 * gradient for the non-linear conjugate-gradient code. 
 * It represents a two-step self-scaled bfgs formula. 
 */
static int msolve(double g[], double y[], int n, 
  double sk[], double yk[], double diagb[], double sr[], 
  double yr[], logical upd1, double yksk, double yrsr, 
  logical lreset)
{
  double ghyk, ghyr, yksr, ykhyk, ykhyr, yrhyr, rdiagb, gsr, gsk;
  int i, frc;
  double *hg = NULL, *hyk = NULL, *hyr = NULL;

  if (upd1)
  {
    for (i = 0; i < n; i++) y[i] = g[i] / diagb[i];
    return 0;
  }

  frc = -1;
  gsk = ddot1(n, g, sk);
  hg = malloc(sizeof(*hg)*n);
  if (hg == NULL) goto cleanup;
  hyr = malloc(sizeof(*hyr)*n);
  if (hyr == NULL) goto cleanup;
  hyk = malloc(sizeof(*hyk)*n);
  if (hyk == NULL) goto cleanup;
  frc = 0;

  /* Compute gh and hy where h is the inverse of the diagonals */
  if (lreset)
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
    }
    ykhyk = ddot1(n, yk, hyk);
    ghyk = ddot1(n, g, hyk);
    ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
      hyr[i] = yr[i] * rdiagb;
    }
    gsr = ddot1(n, g, sr);
    ghyr = ddot1(n, g, hyr);
    yrhyr = ddot1(n, yr, hyr);
    ssbfgs(n, 1.0, sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
    yksr = ddot1(n, yk, sr);
    ykhyr = ddot1(n, yk, hyr);
    ssbfgs(n, 1.0, sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
    ykhyk = ddot1(n, hyk, yk);
    ghyk = ddot1(n, hyk, g);
    ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }
  
cleanup:
  if (hg) free(hg);
  if (hyk) free(hyk);
  if (hyr) free(hyr);
  
  return frc;
}

/*
 * Self-scaled BFGS
 */
static void ssbfgs(int n, double gamma, double sj[], double hjv[],
  double hjyj[], double yjsj, 
  double yjhyj, double vsj, double vhyj, double hjp1v[])
{
  double beta, delta;
  int i;

  if (yjsj == 0.0)
  {
    delta = 0.0;
    beta = 0.0;
  }
  else
  {
    delta = (gamma * yjhyj / yjsj + 1.0) * vsj / yjsj - gamma * vhyj / yjsj;
    beta = -gamma * vsj / yjsj;
  }

  for (i = 0; i < n; i++)
    hjp1v[i] = gamma * hjv[i] + delta * sj[i] + beta * hjyj[i];
}

/*
 * Initialize the preconditioner
 */
static int initPreconditioner(double diagb[], double emat[], int n, 
  logical lreset, double yksk, double yrsr,
  double sk[], double yk[], double sr[], double yr[], 
  logical upd1)
{
  double srds, yrsk, td, sds;
  int i;
  double *bsk;

  if (upd1)
  {
    dcopy1(n, diagb, emat);
    return 0;
  }
  
  bsk = malloc(sizeof(*bsk)*n);
  if (bsk == NULL) return -1;
  
  if (lreset) 
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sk[i];
    sds = ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      emat[i] = td - td * td * sk[i] * sk[i] / sds + yk[i] * yk[i] / yksk;
    }
  }
  else
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sr[i];
    sds = ddot1(n, sr, bsk);
    srds = ddot1(n, sk, bsk);
    yrsk = ddot1(n, yr, sk);
    if (yrsr == 0.0) yrsr = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      bsk[i] = td * sk[i] - bsk[i] * srds / sds + yr[i] * yrsk / yrsr;
      emat[i] = td - td * td * sr[i] * sr[i] / sds + yr[i] * yr[i] / yrsr;
    }
    sds = ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
      emat[i] = emat[i] - bsk[i] * bsk[i] / sds + yk[i] * yk[i] / yksk;
  }
  
  free(bsk);
  return 0;
}


/*
 * Line search algorithm of gill and murray
 */
static ls_rc linearSearch(int n, tnc_function *function, void *state,
  double low[], double up[],
  double xscale[], double fscale, int pivot[],
  double eta, double ftol, double xbnd,
  double p[], double x[], double *f,
  double *alpha, double gfull[], int maxnfeval, int *nfeval)
{
  double b1, big, tol, rmu, fpresn, fu, gu, fw, gw, gtest1, gtest2,
    oldf, fmin, gmin, rtsmll, step, a, b, e, u, ualpha, factor, scxbnd, xw,
    epsmch, reltol, abstol, tnytol, pe, xnorm, rteps;
  double *temp = NULL, *tempgfull = NULL, *newgfull = NULL;
  int maxlsit = 64, i, itcnt, frc;
  ls_rc rc;
  getptc_rc itest;
  logical braktd;
  
  rc = LS_ENOMEM;
  temp = malloc(sizeof(*temp)*n);
  if (temp == NULL) goto cleanup;
  tempgfull = malloc(sizeof(*tempgfull)*n);
  if (tempgfull == NULL) goto cleanup;
  newgfull = malloc(sizeof(*newgfull)*n);
  if (newgfull == NULL) goto cleanup;

  dcopy1(n, gfull, temp);
  scaleg(n, temp, xscale, fscale);
  gu = ddot1(n, temp, p);

  dcopy1(n, x, temp);
  project(n, temp, pivot);
  xnorm = dnrm21(n, temp);

  /* Compute the absolute and relative tolerances for the linear search */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);
  pe = dnrm21(n, p) + epsmch;
  reltol = rteps * (xnorm + 1.0) / pe;
  abstol = -epsmch * (1.0 + fabs(*f)) / (gu - epsmch);

  /* Compute the smallest allowable spacing between points in the linear
    search */
  tnytol = epsmch * (xnorm + 1.0) / pe;

  rtsmll = epsmch;
  big = 1.0 / (epsmch * epsmch);
  itcnt = 0;

  /* Set the estimated relative precision in f(x). */
  fpresn = epsmch * ftol;

  u = *alpha;
  fu = *f;
  fmin = *f;
  rmu = 1e-4;

  /* Setup */
  itest = getptcInit(&reltol, &abstol, tnytol, eta, rmu, 
    xbnd, &u, &fu, &gu, alpha, &fmin, &gmin, &xw, &fw, &gw, &a, &b,
    &oldf, &b1, &scxbnd, &e, &step, &factor, &braktd, &gtest1, &gtest2, &tol);

  /* If itest == GETPTC_EVAL, the algorithm requires the function value to be 
    calculated */
  while(itest == GETPTC_EVAL)
  {
    /* Test for too many iterations or too many function evals */
    if ((++itcnt > maxlsit) || ((*nfeval) >= maxnfeval)) break;

    ualpha = *alpha + u;
    for (i = 0; i < n; i++)
      temp[i] = x[i] + ualpha * p[i];

    /* Function evaluation */
    unscalex(n, temp, xscale);
    coercex(n, temp, low, up);

    frc = function(temp, &fu, tempgfull, state);
    ++(*nfeval);
    if (frc)
    {
      rc = LS_USERABORT;
      goto cleanup;
    }

    fu *= fscale;

    dcopy1(n, tempgfull, temp);
    scaleg(n, temp, xscale, fscale);
    gu = ddot1(n, temp, p);

    itest = getptcIter(big, rtsmll, &reltol, &abstol, tnytol, fpresn,
      xbnd, &u, &fu, &gu, alpha, &fmin, &gmin, &xw, &fw, &gw, &a, &b,
      &oldf, &b1, &scxbnd, &e, &step, &factor, &braktd, &gtest1, &gtest2, &tol);
    
    /* New best point ? */
    if (*alpha == ualpha)
      dcopy1(n, tempgfull, newgfull);
  }

  if (itest == GETPTC_OK)
  {
    /* A successful search has been made */
    *f = fmin;
    daxpy1(n, *alpha, p, x);
    dcopy1(n, newgfull, gfull);
    rc = LS_OK;
  }
  /* Too many iterations ? */
  else if (itcnt > maxlsit) rc = LS_FAIL;
  /* If itest=GETPTC_FAIL or GETPTC_EINVAL a lower point could not be found */
  else if (itest != GETPTC_EVAL) rc = LS_FAIL;
  /* Too many function evaluations */
  else rc = LS_MAXFUN;

cleanup:
  if (temp) free(temp);
  if (tempgfull) free(tempgfull);
  if (newgfull) free(newgfull);

  return rc;
}

/*
 * getptc, an algorithm for finding a steplength, called repeatedly by
 * routines which require a step length to be computed using cubic
 * interpolation. The parameters contain information about the interval
 * in which a lower point is to be found and from this getptc computes a
 * point at which the function can be evaluated by the calling program.
 */
static getptc_rc getptcInit(double *reltol, double *abstol, double tnytol, 
  double eta, double rmu, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol)
{
  /* Check input parameters */
  if (*u <= 0.0 || xbnd <= tnytol || *gu > 0.0) return GETPTC_EINVAL;
  if (xbnd < *abstol) *abstol = xbnd;
  *tol = *abstol;

  /* a and b define the interval of uncertainty, x and xw are points */
  /* with lowest and second lowest function values so far obtained. */
  /* initialize a,smin,xw at origin and corresponding values of */
  /* function and projection of the gradient along direction of search */
  /* at values for latest estimate at minimum. */

  *a = 0.0;
  *xw = 0.0;
  *xmin = 0.0;
  *oldf = *fu;
  *fmin = *fu;
  *fw = *fu;
  *gw = *gu;
  *gmin = *gu;
  *step = *u;
  *factor = 5.0;

  /* The minimum has not yet been bracketed. */
  *braktd = TNC_FALSE;

  /* Set up xbnd as a bound on the step to be taken. (xbnd is not computed */
  /* explicitly but scxbnd is its scaled value.) Set the upper bound */
  /* on the interval of uncertainty initially to xbnd + tol(xbnd). */
  *scxbnd = xbnd;
  *b = *scxbnd + *reltol * fabs(*scxbnd) + *abstol;
  *e = *b + *b;
  *b1 = *b;

  /* Compute the constants required for the two convergence criteria. */
  *gtest1 = -rmu * *gu;
  *gtest2 = -eta * *gu;

  /* If the step is too large, replace by the scaled bound (so as to */
  /* compute the new point on the boundary). */
  if (*step >= *scxbnd)
  {  
    *step = *scxbnd;
    /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
    *scxbnd -= (*reltol * fabs(xbnd) + *abstol) / (1.0 + *reltol);
  }
  *u = *step;
  if (fabs(*step) < *tol && *step < 0.0) *u = -(*tol);
  if (fabs(*step) < *tol && *step >= 0.0) *u = *tol;
  return GETPTC_EVAL;
}

static getptc_rc getptcIter(double big, double 
  rtsmll, double *reltol, double *abstol, double tnytol, 
  double fpresn, double xbnd,
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol)
{
  double abgw, absr, p, q, r, s, scale, denom, 
    a1, d1, d2, sumsq, abgmin, chordm, chordu,
    xmidpt, twotol;
  logical convrg;

  /* Update a,b,xw, and xmin */
  if (*fu <= *fmin)
  {
    /* If function value not increased, new point becomes next */
    /* origin and other points are scaled accordingly. */
    chordu = *oldf - (*xmin + *u) * *gtest1;
    if (*fu > chordu)
    {
      /* The new function value does not satisfy the sufficient decrease */
      /* criterion. prepare to move the upper bound to this point and */
      /* force the interpolation scheme to either bisect the interval of */
      /* uncertainty or take the linear interpolation step which estimates */
      /* the root of f(alpha)=chord(alpha). */

      chordm = *oldf - *xmin * *gtest1;
      *gu = -(*gmin);
      denom = chordm - *fmin;
      if (fabs(denom) < 1e-15)
      {
        denom = 1e-15;
        if (chordm - *fmin < 0.0) denom = -denom;
      }
      if (*xmin != 0.0) *gu = *gmin * (chordu - *fu) / denom;
      *fu = 0.5 * *u * (*gmin + *gu) + *fmin;
      if (*fu < *fmin) *fu = *fmin;
    }
    else
    {
      *fw = *fmin;
      *fmin = *fu;
      *gw = *gmin;
      *gmin = *gu;
      *xmin += *u;
      *a -= *u;
      *b -= *u;
      *xw = -(*u);
      *scxbnd -= *u;
      if (*gu <= 0.0)
      {
        *a = 0.0;
      }
      else
      {
        *b = 0.0;
        *braktd = TNC_TRUE;
      }
      *tol = fabs(*xmin) * *reltol + *abstol;
      goto ConvergenceCheck;
    }
  }

  /* If function value increased, origin remains unchanged */
  /* but new point may now qualify as w. */
  if (*u < 0.0)
    *a = *u;
  else
  {
    *b = *u;
    *braktd = TNC_TRUE;
  }
  *xw = *u;
  *fw = *fu;
  *gw = *gu;

ConvergenceCheck:
  twotol = *tol + *tol;
  xmidpt = 0.5 * (*a + *b);

  /* Check termination criteria */
  convrg = (fabs(xmidpt) <= twotol - 0.5 * (*b - *a)) || 
    (fabs(*gmin) <= *gtest2 && *fmin < *oldf && ((fabs(*xmin - xbnd) > *tol) ||
    (! (*braktd))));
  if (convrg)
  {
    if (*xmin != 0.0) return GETPTC_OK;

    /*
     * If the function has not been reduced, check to see that the relative
     * change in f(x) is consistent with the estimate of the delta-
     * unimodality constant, tol. If the change in f(x) is larger than
     * expected, reduce the value of tol.
     */
    if (fabs(*oldf - *fw) <= fpresn * 0.5 * (fabs(*fw) + fabs(*oldf)))
      return GETPTC_FAIL;
    *tol = 0.1 * *tol;
    if (*tol < tnytol) return GETPTC_FAIL;
    *reltol = 0.1 * *reltol;
    *abstol = 0.1 * *abstol;
    twotol = 0.1 * twotol;
  }

  /* Continue with the computation of a trial step length */
  r = 0.0;
  q = 0.0;
  s = 0.0;
  if (fabs(*e) > *tol)
  {
    /* Fit cubic through xmin and xw */
    r = 3.0 * (*fmin - *fw) / *xw + *gmin + *gw;
    absr = fabs(r);
    q = absr;
    if (*gw != 0.0 && *gmin != 0.0)
    {
      /* Compute the square root of (r*r - gmin*gw) in a way
         which avoids underflow and overflow. */
      abgw = fabs(*gw);
      abgmin = fabs(*gmin);
      s = sqrt(abgmin) * sqrt(abgw);
      if (*gw / abgw * *gmin > 0.0) 
      {
        if (r >= s || r <= -s)
        {
          /* Compute the square root of r*r - s*s */
          q = sqrt(fabs(r + s)) * sqrt(fabs(r - s));
        }
        else
        {
          r = 0.0;
          q = 0.0;
          goto MinimumFound;
        }
      }
      else
      {
        /* Compute the square root of r*r + s*s. */
        sumsq = 1.0;
        p = 0.0;
        if (absr >= s)
        {
          /* There is a possibility of underflow. */
          if (absr > rtsmll) p = absr * rtsmll;
          if (s >= p)
          {
            double value = s / absr;
            sumsq = 1.0 + value * value;
          }
          scale = absr;
        }
        else
        {
          /* There is a possibility of overflow. */
          if (s > rtsmll) p = s * rtsmll;
          if (absr >= p)
          {
            double value = absr / s;
            sumsq = 1.0 + value * value;
          }
          scale = s;
        }
        sumsq = sqrt(sumsq);
        q = big;
        if (scale < big / sumsq) q = scale * sumsq;
      }
    }

    /* Compute the minimum of fitted cubic */
    if (*xw < 0.0) q = -q;
    s = *xw * (*gmin - r - q);
    q = *gw - *gmin + q + q;
    if (q > 0.0) s = -s;
    if (q <= 0.0) q = -q;
    r = *e;
    if (*b1 != *step || *braktd) *e = *step;
  }

MinimumFound:
  /* Construct an artificial bound on the estimated steplength */
  a1 = *a;
  *b1 = *b;
  *step = xmidpt;
  if ( (! *braktd) || ((*a == 0.0 && *xw < 0.0) || (*b == 0.0 && *xw > 0.0)) )
  {
    if (*braktd)
    {
      /* If the minimum is not bracketed by 0 and xw the step must lie
         within (a1,b1). */
      d1 = *xw;
      d2 = *a;
      if (*a == 0.0) d2 = *b;
      /* This line might be : */
      /* if (*a == 0.0) d2 = *e */
      *u = -d1 / d2;
      *step = 5.0 * d2 * (0.1 + 1.0 / *u) / 11.0;
      if (*u < 1.0) *step = 0.5 * d2 * sqrt(*u);
    }
    else
    {
      *step = -(*factor) * *xw;
      if (*step > *scxbnd) *step = *scxbnd;
      if (*step != *scxbnd) *factor = 5.0 * *factor;
    }
    /* If the minimum is bracketed by 0 and xw the step must lie within (a,b) */
    if (*step <= 0.0) a1 = *step;
    if (*step > 0.0) *b1 = *step;
  }

/*
 *   Reject the step obtained by interpolation if it lies outside the
 *   required interval or it is greater than half the step obtained
 *   during the last-but-one iteration.
 */
  if (fabs(s) <= fabs(0.5 * q * r) || s <= q * a1 || s >= q * *b1)
    *e = *b - *a;
  else
  {
    /* A cubic interpolation step */
    *step = s / q;

    /* The function must not be evaluated too close to a or b. */
    if (*step - *a < twotol || *b - *step < twotol)
    {
      if (xmidpt <= 0.0)
        *step = -(*tol);
      else
        *step = *tol;
    }
  }

  /* If the step is too large, replace by the scaled bound (so as to */
  /* compute the new point on the boundary). */
  if (*step >= *scxbnd)
  {  
    *step = *scxbnd;
    /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
    *scxbnd -= (*reltol * fabs(xbnd) + *abstol) / (1.0 + *reltol);
  }
  *u = *step;
  if (fabs(*step) < *tol && *step < 0.0) *u = -(*tol);
  if (fabs(*step) < *tol && *step >= 0.0) *u = *tol;
  return GETPTC_EVAL;
}

/*
 * Return epsmch, where epsmch is the smallest possible
 * power of 2 such that 1.0 + epsmch > 1.0
 */
static double mchpr1(void)
{
  static double epsmch = 0.0;
  
  if (epsmch == 0.0)
  {
    double eps = 1.0;
    while((1.0 + (eps*0.5)) > 1.0)
      eps *= 0.5;
    epsmch = eps;
  }

  return epsmch;
}

/* Blas like routines */

/* dy+=dx */
static void dxpy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += dx[i];
}

/* dy+=da*dx */
static void daxpy1(int n, double da, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += da*dx[i];
}

/* Copy dx -> dy */
/* Could use memcpy */
static void dcopy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] = dx[i];
}

/* Negate */
static void dneg1(int n, double v[])
{
  int i;
  for (i = 0; i < n; i++)
    v[i] = -v[i];
}

/* Dot product */
static double ddot1(int n, double dx[], double dy[])
{
  int i;
  double dtemp = 0.0;
  for (i = 0; i < n; i++)
    dtemp += dy[i]*dx[i];
  return dtemp;
}

/* Infinity norm */
static double dnrmi1(int n, double v[])
{
  int i;
  double dtemp, dmax;
  for (dmax = fabs(v[0]), i = 1; i < n; i++)
    if ((dtemp = fabs(v[i])) > dmax) dmax = dtemp;
  return dmax;
}

/* Euclidian norm */
static double dnrm21(int n, double dx[])
{
  int i;
  double dssq = 1.0, dscale = 0.0;

  for (i = 0; i < n; i++)
  {
    if (dx[i] != 0.0)
    {
      double dabsxi = fabs(dx[i]);
      if (dscale<dabsxi)
      {
        /* Normalization to prevent overflow */
        double ratio = dscale/dabsxi;
        dssq = 1.0 + dssq*ratio*ratio;
        dscale = dabsxi;
      }
      else
      {
        double ratio = dabsxi/dscale;
        dssq += ratio*ratio;
      }
    }
  }

  return dscale*sqrt(dssq);
}
