/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_vmt.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VMT                                *
 *         (Vector Matrix Transformation)                                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  THIS METHOD AND THE CORRESPONDING ROUTINES SHOULD NOT BE USED ANY MORE!  *
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

/* 
   =deprecatedMETHOD  VMT   Vector Matrix Transformation

   =UP  Methods_for_CVEC

   =REQUIRED mean vector, covariance matrix, standardized marginal distributions

   =SPEED Set-up: slow,
          Sampling: depends on dimension

   =REINIT supported

   =REF  [HLD04: Sect.11.1.6, Alg.11.3.]

   =DESCRIPTION
      VMT generates random vectors for distributions with given mean
      vector mu and covariance matrix @unurmath{Sigma.} It produces random vectors
      of the form @unurmath{X = L Y + mu,} where @unurmath{L} is the
      Cholesky factor of Sigma, i.e. @unurmath{L L^t = Sigma,} and
      @unurmath{Y} has independent components of the same 
      distribution with mean 0 and standard deviation 1.
      
      The method VMT has been implemented especially to sample from a
      multinormal distribution. Nevertheless, it can also be used (or
      abused) for other distributions. However, notice that the given
      standardized marginal distributions are not checked; i.e.
      if the given distributions do not have mean 0 and variance 1
      then @unurmath{mu} and @unurmath{Sigma} are not the mean vector
      and covariance matrix, respectively, of the resulting
      distribution.

      @strong{Important:} Notice that except for the multinormal
      distribution the given marginal distribution are distorted by
      the transformation using the Cholesky matrix. Thus for other
      (non-multinormal) distributions this method should only be used
      when everything else fails and some approximate results which
      might even be not entirely correct are better than no results.
      A much better method is the NORTA (NORmal To Anything) method, 
      see @ref{NORTA}.

   =HOWTOUSE
      Create a multivariate generator object, set mean vector and
      covariance matrix by means of the unur_distr_cvec_set_mean() and  
      unur_distr_cvec_set_covar() call. Set standard marginal 
      distributions using unur_distr_cvec_set_stdmarginals() , 
      unur_distr_cvec_set_stdmarginal_array() , or 
      unur_distr_cvec_set_stdmarginal_list().
      (Do not use the corresponding calls for the (non-standard)
      marginal distributions).

      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

      The method can also be used with multivariate distribution with
      a truncated domain, i.e., where the domain has been set by a 
      unur_distr_cvec_set_domain_rect() call. However, it then uses a
      simple rejection method that can have extremely poor rejection
      constant especially when dimension is (even moderately) high.

      There are no optional parameters for this method.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vmt_new( const UNUR_DISTR *distribution );
/* 
   Get parameters for generator.
*/

/* =END */
/*---------------------------------------------------------------------------*/


