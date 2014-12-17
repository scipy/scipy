/* tnc : truncated newton bound constrained minimization
         using gradient information, in C */

/*
 * Copyright (c) 2002-2005, Jean-Sebastien Roy (js@jeannot.org)
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

/* $Jeannot: tnc.h,v 1.55 2005/01/28 18:27:31 js Exp $ */

#ifndef _TNC_
#define _TNC_

#define TNC_VERSION "1.3"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Verbosity level
 */
typedef enum {
  TNC_MSG_NONE = 0, /* No messages */
  TNC_MSG_ITER = 1, /* One line per iteration */
  TNC_MSG_INFO = 2, /* Informational messages */
  TNC_MSG_VERS = 4, /* Version info */
  TNC_MSG_EXIT = 8, /* Exit reasons */

  TNC_MSG_ALL = TNC_MSG_ITER | TNC_MSG_INFO
    | TNC_MSG_VERS | TNC_MSG_EXIT /* All messages */
} tnc_message;

/*
 * Possible return values for tnc
 */
typedef enum
{
  TNC_MINRC        = -3, /* Constant to add to get the rc_string */
  TNC_ENOMEM       = -3, /* Memory allocation failed */
  TNC_EINVAL       = -2, /* Invalid parameters (n<0) */
  TNC_INFEASIBLE   = -1, /* Infeasible (low bound > up bound) */
  TNC_LOCALMINIMUM =  0, /* Local minima reach (|pg| ~= 0) */
  TNC_FCONVERGED   =  1, /* Converged (|f_n-f_(n-1)| ~= 0) */
  TNC_XCONVERGED   =  2, /* Converged (|x_n-x_(n-1)| ~= 0) */
  TNC_MAXFUN       =  3, /* Max. number of function evaluations reach */
  TNC_LSFAIL       =  4, /* Linear search failed */
  TNC_CONSTANT     =  5, /* All lower bounds are equal to the upper bounds */
  TNC_NOPROGRESS   =  6, /* Unable to progress */
  TNC_USERABORT    =  7  /* User requested end of minization */
} tnc_rc;

/*
 * Return code strings
 * use tnc_rc_string[rc - TNC_MINRC] to get the message associated with
 * return code rc.
 */

extern const char *const tnc_rc_string[11];

/*
 * A function as required by tnc
 * state is a void pointer provided to the function at each call
 *
 * x     : on input, then vector of variables (should not be modified)
 * f     : on output, the value of the function
 * g     : on output, the value of the gradient
 * state : on input, the value of the state variable as provided to tnc
 *
 * must returns 0 if no error occurs or 1 to immediately end the minimization.
 *
 */
typedef int tnc_function(double x[], double *f, double g[], void *state);

/*
 * A callback function accepting x as input parameter along with the state
 * pointer.
 */
typedef void tnc_callback(double x[], void *state);

/*
 * tnc : minimize a function with variables subject to bounds, using
 *       gradient information.
 *
 * n         : number of variables (must be >= 0)
 * x         : on input, initial estimate ; on output, the solution
 * f         : on output, the function value at the solution
 * g         : on output, the gradient value at the solution
 *             g should be an allocated vector of size n or NULL,
 *             in which case the gradient value is not returned.
 * function  : the function to minimize (see tnc_function)
 * state     : used by function (see tnc_function)
 * low, up   : the bounds
 *             set low[i] to -HUGE_VAL to remove the lower bound
 *             set up[i] to HUGE_VAL to remove the upper bound
 *             if low == NULL, the lower bounds are removed.
 *             if up == NULL, the upper bounds are removed.
 * scale     : scaling factors to apply to each variable
 *             if NULL, the factors are up-low for interval bounded variables
 *             and 1+|x] for the others.
 * offset    : constant to substract to each variable
 *             if NULL, the constant are (up+low)/2 for interval bounded
 *             variables and x for the others.
 * messages  : see the tnc_message enum
 * maxCGit   : max. number of hessian*vector evaluation per main iteration
 *             if maxCGit == 0, the direction chosen is -gradient
 *             if maxCGit < 0, maxCGit is set to max(1,min(50,n/2))
 * maxnfeval : max. number of function evaluation
 * eta       : severity of the line search. if < 0 or > 1, set to 0.25
 * stepmx    : maximum step for the line search. may be increased during call
 *             if too small, will be set to 10.0
 * accuracy  : relative precision for finite difference calculations
 *             if <= machine_precision, set to sqrt(machine_precision)
 * fmin      : minimum function value estimate
 * ftol      : precision goal for the value of f in the stoping criterion
 *             if ftol < 0.0, ftol is set to accuracy
 * xtol      : precision goal for the value of x in the stopping criterion
 *             (after applying x scaling factors)
 *             if xtol < 0.0, xtol is set to sqrt(machine_precision)
 * pgtol     : precision goal for the value of the projected gradient in the
 *             stopping criterion (after applying x scaling factors)
 *             if pgtol < 0.0, pgtol is set to 1e-2 * sqrt(accuracy)
 *             setting it to 0.0 is not recommended
 * rescale   : f scaling factor (in log10) used to trigger f value rescaling
 *             if 0, rescale at each iteration
 *             if a big value, never rescale
 *             if < 0, rescale is set to 1.3
 * nfeval    : on output, the number of function evaluations.
 *             ignored if nfeval==NULL.
 *
 * The tnc function returns a code defined in the tnc_rc enum.
 * On output, x, f and g may be very slightly out of sync because of scaling.
 *
 */
extern int tnc(int n, double x[], double *f, double g[],
  tnc_function *function, void *state,
  double low[], double up[], double scale[], double offset[],
  int messages, int maxCGit, int maxnfeval, double eta, double stepmx,
  double accuracy, double fmin, double ftol, double xtol, double pgtol,
  double rescale, int *nfeval, int *niter, tnc_callback *callback);

#ifdef __cplusplus
}
#endif

#endif /* _TNC_ */
