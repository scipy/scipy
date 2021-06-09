/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_specfunct_source.h                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         prototypes and macros for using special functions like erf(),     *
 *         gamma(), beta(), etc., which are imported from other packages.    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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
#ifndef UNUR_SPECFUNCT_SOURCE_H_SEEN
#define UNUR_SPECFUNCT_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Prototypes for special functions like erf(), gamma(), beta(), etc.      *
 *   which are imported from other packages.                                 *
 *                                                                           *
 *   We use the package CEPHES/DOUBLE for computing these functions          *
 *   (available from NETLIB, http://www.netlib.org/cephes/                   *
 *   Copyright 1984 - 1994 by Stephen L. Moshier                             *
 *                                                                           *
 *   Alternatively, we also can use the functions from the Rmath library     *
 *   from the R project for statistical computing, http://www.R-project.org/ *
 *                                                                           *
 *****************************************************************************/

/* We define macros for special functions.
 *
 * The following macros must be defined:
 *
 *   _unur_SF_incomplete_beta   ... incomplete beta integral
 *   _unur_SF_ln_gamma          ... logarithm of gamma function
 *   _unur_SF_ln_factorial      ... logarithm of factorial
 *   _unur_SF_incomplete_gamma  ... incomplete gamma function
 *   _unur_SF_cdf_normal        ... CDF of normal distribution
 *   _unur_SF_invcdf_normal     ... inverse CDF of normal distribution
 *
 *---------------------------------------------------------------------------*/

#ifdef HAVE_LIBRMATH

/*---------------------------------------------------------------------------*/
/* Routines from the Rmath library (R project).                              */
/*---------------------------------------------------------------------------*/

/* we have to distinguish between two cases: */
#  ifdef R_UNURAN
/*   Rmath for 'Runuran': nothing special to do. */
#  else
/*   Rmath standalone library. */
#    define MATHLIB_STANDALONE
#  endif

/* include Rmath header file */
#  include <Rmath.h>

/* we have to #undef some macros from Rmath.h */
#ifdef trunc
#undef trunc
#endif

#ifdef beta
#undef beta
#endif

/* ......................................................................... */

/* incomplete beta integral */
#define _unur_SF_incomplete_beta(x,a,b)   pbeta((x),(a),(b),TRUE,FALSE)

/* logarithm of gamma function */
#define _unur_SF_ln_gamma(x)              lgammafn(x)

/* logarithm of factorial */
#define _unur_SF_ln_factorial(x)          lgammafn((x)+1.)

/* incomplete gamma function */
#define _unur_SF_incomplete_gamma(x,a)    pgamma(x,a,1.,TRUE,FALSE)

/* modified Bessel function K_nu of second kind (AKA third kind) */
#define _unur_SF_bessel_k(x,nu)           bessel_k((x),(nu),1)

/* Normal distribution */
#define _unur_SF_cdf_normal(x)            pnorm((x),0.,1.,TRUE,FALSE)
#define _unur_SF_invcdf_normal(u)         qnorm((u),0.,1.,TRUE,FALSE)

/* ..........................................................................*/

/* Beta Distribution */
#define _unur_SF_invcdf_beta(u,p,q)       qbeta((u),(p),(q),TRUE,FALSE)

/* F Distribution */
#define _unur_SF_cdf_F(x,nua,nub)         pf((x),(nua),(nub),TRUE,FALSE)
#define _unur_SF_invcdf_F(u,nua,nub)      qf((u),(nua),(nub),TRUE,FALSE)

/* Gamma Distribution */
#define _unur_SF_invcdf_gamma(u,shape,scale)  qgamma((u),(shape),(scale),TRUE,FALSE)

/* Student t Distribution */
#define _unur_SF_cdf_student(x,nu)         pt((x),(nu),TRUE,FALSE)
#define _unur_SF_invcdf_student(u,nu)      qt((u),(nu),TRUE,FALSE)

/* Binomial Distribution */
#define _unur_SF_invcdf_binomial(u,n,p)   qbinom((u),(n),(p),TRUE,FALSE)

/* Hypergeometric Distribution */
#define _unur_SF_cdf_hypergeometric(x,N,M,n)  phyper((x),(M),(N)-(M),(n),TRUE,FALSE)
#define _unur_SF_invcdf_hypergeometric(u,N,M,n)  qhyper((u),(M),(N)-(M),(n),TRUE,FALSE)

/* Negative Binomial Distribution */
#define _unur_SF_cdf_negativebinomial(x,n,p)      pnbinom((x),(n),(p),TRUE,FALSE)
#define _unur_SF_invcdf_negativebinomial(u,n,p)   qnbinom((u),(n),(p),TRUE,FALSE)

/* Poisson Distribution */
#define _unur_SF_invcdf_poisson(u,theta)   qpois((u),(theta),TRUE,FALSE)

/*---------------------------------------------------------------------------*/
/* end: Rmath library (R project)                                            */
/*---------------------------------------------------------------------------*/

#else

/*---------------------------------------------------------------------------*/
/* Routines from the CEPHES library.                                         */
/*---------------------------------------------------------------------------*/

/* incomplete beta integral */
double _unur_cephes_incbet(double a, double b, double x);
#define _unur_SF_incomplete_beta(x,a,b)   _unur_cephes_incbet((a),(b),(x))

/* logarithm of gamma function */
double _unur_cephes_lgam(double x);
#define _unur_SF_ln_gamma(x)              _unur_cephes_lgam(x)

/* logarithm of factorial */
#define _unur_SF_ln_factorial(x)          _unur_cephes_lgam((x)+1.)

/* incomplete gamma function */
double _unur_cephes_igam(double a, double x);
#define _unur_SF_incomplete_gamma(x,a)    _unur_cephes_igam((a),(x))

/* normal distribution function */
double _unur_cephes_ndtr(double x);
#define _unur_SF_cdf_normal(x)            _unur_cephes_ndtr(x)

/* inverse of normal distribution function */
double _unur_cephes_ndtri(double x);
#define _unur_SF_invcdf_normal(x)         _unur_cephes_ndtri(x)

/*---------------------------------------------------------------------------*/
/* end: CEPHES library                                                       */
/*---------------------------------------------------------------------------*/

#endif

/*****************************************************************************
 *                                                                           *
 *   Replacement for missing (system) functions                              *
 *                                                                           *
 *****************************************************************************/

#if !HAVE_DECL_LOG1P
/* log(1+x) */
/* (replacement for missing C99 function log1p) */
double _unur_log1p(double x);
#endif

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_SPECFUNCT_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
