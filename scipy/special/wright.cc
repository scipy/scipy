#include "wright.hh"

/**********************************************************************/
/* wrightomega is the simple routine for evaluating the wright omega  */
/* function.                                                          */
/*                                                                    */
/* Calling:                                                           */
/*    w = wrightomega(z)                                              */
/*                                                                    */
/* Input:                                                             */
/*   z  --  double complex                                            */
/*          Value to evaluate Wrightomega(z) at.                      */
/*                                                                    */
/* Output:                                                            */
/*   w  --  double complex                                            */
/*          Value of Wrightomega(z)                                   */
/*                                                                    */
/**********************************************************************/

/*
  Also published as ACM TOMS 917; relicensed as BSD by the author.

  Copyright (C) Piers Lawrence.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


/**********************************************************************/
/* wrightomega_ext is the extended routine for evaluating the wright  */
/* omega function with the option of extracting the last update step, */
/* the penultimate residual and the condition number estimate.        */
/*                                                                    */
/* Calling:                                                           */
/*   success = wrightomega_ext(z,w,e,r,cond);                         */
/*                                                                    */
/* Input:                                                             */
/*   z  --  double complex                                            */
/*          Value to evaluate Wrightomega(z) at.                      */
/*                                                                    */
/*   w  --  double complex*                                           */
/*          Pointer to return value of Wrightomega(z).                */
/*                                                                    */
/*   cond  --  double complex*                                        */
/*         Pointer to the condition number estimate. If NULL the      */
/*         condition number is not calculated.                        */
/*                                                                    */
/* Output: returns 0 on successful exit.                               */
/**********************************************************************/

#include <Python.h>
extern "C" {
#include <numpy/npy_math.h>
#include "sf_error.h"
#include "_c99compat.h"
#include "_round.h"
}

#include <cmath>
#include <cfloat>

using std::complex;

#define NaN NPY_NAN
#define TWOITERTOL DBL_EPSILON

const complex<double> I(0.0, 1.0);


int
wright::wrightomega_ext(complex<double> z, complex<double> *w,
			complex<double> *cond)
{
  double pi = NPY_PI, s = 1.0;
  double x, y, ympi, yppi, near;
  complex<double> e, r, pz, wp1, t, fac;


  /* extract real and imaginary parts of z */
  x=real(z);
  y=imag(z);

  /* compute if we are near the branch cuts */
  ympi=y-pi;
  yppi=y+pi;
  near=0.1e-1;

  /* Test for floating point exceptions */
  /*****************************/
  /* NaN output for NaN input  */
  /*****************************/
  if(sc_isnan(x) || sc_isnan(y))
    {
      *w = complex<double>(NaN, NaN);
      return 0;
    }
  /*********************************/
  /* Signed zeros between branches */
  /*********************************/
  else if(sc_isinf(x) && (x < 0.0) && (-pi < y) && (y<= pi))
    {
      if (fabs(y) <= pi/2.0)
        {
	  if (y >= 0)
	    {
	      *w = complex<double>(0.0, 0.0);
	    }
	  else
	    {
	      *w = complex<double>(0.0, -0.0);
	    }
	}
      else
        {
	  if (y >= 0)
	    {
	      *w = complex<double>(-0.0, 0.0);
	    }
	  else
	    {
	      *w = complex<double>(-0.0, -0.0);
	    }
        }
      return 0;
    }
  /**************************/
  /* Asymptotic for large z */
  /**************************/
  else if(sc_isinf(x) || sc_isinf(y))
    {
      *w = complex<double>(x, y);
      return 0;
    }

  /******************************************/
  /* Test If exactly on the singular points */
  /******************************************/
  if((x==-1.0) && (fabs(y)==pi))
    {
      *w = complex<double>(-1.0, 0.0);
      return 0;
    }


  /* Choose approximation based on region */
  /**********************************/
  /* Region 1: upper branch point   */
  /* Series about z=-1+Pi*I         */
  /**********************************/
  if ((-2.0<x && x<=1.0 && 1.0<y && y< 2.0*pi))
    {
      pz=conj(sqrt(conj(2.0*(z+1.0-I*pi))));
      *w=-1.0+(I+(1.0/3.0+(-1.0/36.0*I+(1.0/270.0+1.0/4320.0*I*pz)*pz)*pz)*pz)*pz;
    }
  /**********************************/
  /* Region 2: lower branch point   */
  /* Series about z=-1-Pi*I         */
  /**********************************/
  else if ((-2.0<x && x<=1.0 && -2.0*pi<y && y<-1.0))
    {
      pz=conj(sqrt(conj(2.0*(z+1.0+I*pi))));
      *w=-1.0+(-I+(1.0/3.0+(1.0/36.0*I+(1.0/270.0-1.0/4320.0*I*pz)*pz)*pz)*pz)*pz;
    }
  /*********************************/
  /* Region 3: between branch cuts */
  /* Series: About -infinity       */
  /*********************************/
  else if (x <= -2.0 && -pi < y && y <= pi)
    {
      pz=exp(z);
      *w=(1.0+(-1.0+(3.0/2.0+(-8.0/3.0+125.0/24.0*pz)*pz)*pz)*pz)*pz;
      if (*w == 0.0)
	{
	  sf_error("wrightomega", SF_ERROR_UNDERFLOW, "underflow in exponential series");
	  /* Skip the iterative scheme because it computes log(*w) */
	  if (cond != NULL)
	    {
	      *cond = z/(1.0+*w);
	    }
	  return 0;
	}
    }
  /**********************/
  /* Region 4: Mushroom */
  /* Series about z=1   */
  /**********************/
  else if (((-2.0 < x) && (x<=1.0) && (-1.0 <= y) && (y <= 1.0))
           || ((-2.0 < x) && (x-0.1e1)*(x-0.1e1)+y*y<=pi*pi))
    {
      pz=z-1.0;
      *w=1.0/2.0+1.0/2.0*z+(1.0/16.0+(-1.0/192.0+(-1.0/3072.0+13.0/61440.0*pz)*pz)*pz)*pz*pz;
    }
  /*************************/
  /* Region 5: Top wing    */
  /* Negative log series   */
  /*************************/
  else if (x<=-0.105e1 && pi<y && y-pi<=-0.75e0*(x+0.1e1))
    {
      t=z-I*pi;
      pz=log(-t);
      *w = t - pz;
      fac = pz/t;
      *w += fac;
      fac /= t;
      *w += fac*(0.5*pz - 1.0);
      fac /= t;
      *w += fac*(pz*pz/3.0 - 3.0*pz/2.0 + 1.0);
      if (abs(z) > 1e50)
	/* Series is accurate and the iterative scheme could overflow */
	{
	  if (cond != NULL)
	    {
	      *cond = z/(1.0+*w);
	    }
	  return 0;
	}
    }
  /***************************/
  /* Region 6: Bottom wing   */
  /* Negative log series     */
  /***************************/
  else if (x<=-0.105e1 && 0.75e0*(x+0.1e1)< y+pi && y+pi<=0.0e0)
    {
      t=z+I*pi;
      pz=log(-t);
      *w = t - pz;
      fac = pz/t;
      *w += fac;
      fac /= t;
      *w += fac*(0.5*pz - 1.0);
      fac /= t;
      *w += fac*(pz*pz/3.0 - 3.0*pz/2.0 + 1.0);
      if (abs(z) > 1e50)
	/* Series is accurate and the iterative scheme could overflow */
	{
	  if (cond != NULL)
	    {
	      *cond = z/(1.0+*w);
	    }
	  return 0;
	}
    }
  /************************************/
  /* Region 7: Everywhere else        */
  /* Series solution about infinity   */
  /************************************/
  else
    {
      pz=log(z);
      *w = z - pz;
      fac = pz/z;
      *w += fac;
      fac /= z;
      *w += fac*(0.5*pz - 1.0);
      fac /= z;
      *w += fac*(pz*pz/3.0 - 3.0*pz/2.0 + 1.0);
      if (abs(z) > 1e50)
	/* Series is accurate and the iterative scheme could overflow */
	{
	  if (cond != NULL)
	    {
	      *cond = z/(1.0+*w);
	    }
	  return 0;
	}
    }

  /**********************************/
  /* Regularize if near branch cuts */
  /**********************************/
  if (x <= -0.1e1 + near && (fabs(ympi) <= near || fabs(yppi) <= near))
    {
      s = -1.0;
      if (fabs(ympi) <= near)
        {
          /* Recompute ympi with directed rounding */
          ympi = add_round_up(y, -pi);

          if( ympi <= 0.0)
            {
              ympi = add_round_down(y, -pi);
            }

          z = x + I*ympi;
        }
      else
        {
          /* Recompute yppi with directed rounding */
          yppi = add_round_up(y, pi);

          if( yppi <= 0.0)
            {
              yppi = add_round_down(y, pi);
            }

          z = x + I*yppi;
        }
    }

  /*****************/
  /* Iteration one */
  /*****************/
  *w=s**w;
  r=z-s**w-log(*w);
  wp1=s**w+1.0;
  e=r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
  *w=*w*(1.0+e);

  /*****************/
  /* Iteration two */
  /*****************/
  if(abs((2.0**w**w-8.0**w-1.0)*pow(abs(r),4.0)) >= TWOITERTOL*72.0*pow(abs(wp1),6.0) )
    {
      r=z-s**w-log(*w);
      wp1=s**w+1.0;
      e=r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
      *w=*w*(1.0+e);
    }

  /***********************/
  /* Undo regularization */
  /***********************/
  *w=s**w;

  /***************************************************/
  /*  Provide condition number estimate if requested */
  /***************************************************/
  if(cond != NULL)
    {
      *cond = z/(1.0+*w);
    }

  return 0;
}

complex<double>
wright::wrightomega(complex<double> z)
{
  complex<double> w;
  wrightomega_ext(z,&w,NULL);
  return w;
}


double
wright::wrightomega_real(double x)
{
  double w, wp1, e, r;

  /* NaN output for NaN input  */
  if (sc_isnan(x))
    {
      return x;
    }

  /* Positive infinity is asymptotically x */
  /* Negative infinity is zero */
  if (sc_isinf(x))
    {
      if (x > 0.0)
	{
	  return x;
	}
      else
	{
	  return 0.0;
	}
    }

  if (x < -50.0)
    {
      /*
       * Skip the iterative scheme because it exp(x) is already
       * accurate to double precision.
       */
      w = exp(x);
      if (w == 0.0)
        {
          sf_error("wrightomega", SF_ERROR_UNDERFLOW, "underflow in exponential series");
        }
      return w;
    }
  if (x > 1e20)
    {
      /*
       * Skip the iterative scheme because the result is just x to
       * double precision
       */
      return x;
    }

  /* Split int three distinct intervals (-inf,-2), [-2,1), [1,inf) */
  if (x < -2.0)
    {
      /* exponential is approx < 1.3e-1 accurate */
      w = exp(x);
    }
  else if (x < 1)
    {
      /* on [-2,1) approx < 1.5e-1 accurate */
      w = exp(2.0*(x-1.0)/3.0);
    }
  else
    {
      /* infinite series with 2 terms approx <1.7e-1 accurate */
      w = log(x);
      w = x - w + w/x;
    }

  /* Iteration one of Fritsch, Shafer, and Crowley (FSC) iteration */
  r = x - w - log(w);
  wp1 = w + 1.0;
  e = r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
  w = w*(1.0+e);

  /* Iteration two (if needed based on the condition number) */
  if (fabs((2.0*w*w-8.0*w-1.0)*pow(fabs(r),4.0)) >= TWOITERTOL*72.0*pow(fabs(wp1),6.0))
    {
      r = x - w - log(w);
      wp1 = w+1.0;
      e = r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
      w = w*(1.0+e);
    }

  return w;
}
