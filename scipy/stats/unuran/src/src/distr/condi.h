/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: condi.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         id   CONDI  (continuous full conditional distribution)            *
 *         type CONT   (continuous univariate distribution)                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/* 
   =NODEX   CONDI   Continuous univariate full conditional distribution

   =UP Distribution_objects [35]

   =DESCRIPTION
      Full conditional distribution for a given continuous
      multivariate distributiion. The condition is a position vector
      and either a variable that is variated or a vector that
      indicates the direction on which the random vector can variate.

      There is a subtle difference between using direction
      vector and using the @var{k}-th variable.
      When a direction vector is given the PDF of the conditional
      distribution is defined by
      @unurmath{f(t) = PDF(pos + t\cdot dir).}
      When a variable is selected the full conditional distribution
      with all other variables fixed is used.

      This is a special case of a continuous univariate distribution
      and thus they have most of these parameters (with the exception
      that functions cannot be changed). Additionally,

      @itemize @minus
      @item there is a call to extract the underlying multivariate
            distribution,

      @item and a call to handle the variables that are fixed and the
            direction for changing the random vector.

      @end itemize

      This distibution type is primarily used for evaluation the
      conditional distribution and its derivative (as required for,
      e.g., the Gibbs sampler). The density is not normalized (i.e. does
      not integrate to one). Mode and area are not available and it
      does not make sense to use any call to set or change parameters
      except the ones given below.

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous full conditional
   distributions (CONDI).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_condi_new( const UNUR_DISTR *distribution, const double *pos, const double *dir, int k );
/* 
   Create an object for full conditional distribution for the given
   @var{distribution}. The condition is given by a position vector
   @var{pos} and either the @var{k}-th variable that is variated or
   the vector @var{dir} that contains the direction on which the
   random vector can variate.

   @var{distribution} must be a pointer to a multivariate continuous
   distribution. 
   @var{pos} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension of the underlying distribution object.
   @var{dir} must be a pointer to an array if size @code{dim} or NULL.
   @var{k} must be in the range @code{0, @dots{}, dim-1}.
   If the @var{k}-th variable is used, @var{dir} must be set to NULL. 

   @emph{Notice:} There is a subtle difference between using direction
   vector @var{dir} and using the @var{k}-th variable.
   When @var{dir} is given, the current position @var{pos} is mapped into
   0 of the conditional distribution and the derivative is taken from
   the function PDF(@var{pos}+t*@var{dir}) w.r.t. @i{t}.
   On the other hand, when the coordinate @var{k} is used (i.e., when
   @var{dir} is set to NULL), the full conditional distribution of the
   distribution is considered (as used for the Gibbs sampler).
   In particular, the current point is just projected into the
   one-dimensional subspace without mapping it into the point 0.

   @emph{Notice:} If a coordinate @var{k} is used, then the @var{k}-th
   partial derivative is used if it as available. Otherwise the
   gradient is computed and the @var{k}-th component is returned.

   The resulting generator object is of the same type as of a
   unur_distr_cont_new() call.
*/

int unur_distr_condi_set_condition( struct unur_distr *distribution, const double *pos, const double *dir, int k );
/* 
   Set/change condition for conditional @var{distribution}. 
   Change values of fixed variables to @var{pos} and use direction
   @var{dir} or @var{k}-th variable of conditional @var{distribution}.

   @var{pos} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension of the underlying distribution object.
   @var{dir} must be a pointer to an array if size @code{dim} or NULL.
   @var{k} must be in the range @code{0, @dots{}, dim-1}.
   If the @var{k}-th variable is used, @var{dir} must be set to NULL. 

   @emph{Notice:} There is a subtle difference between using direction
   vector @var{dir} and using the @var{k}-th variable.
   When @var{dir} is given, the current position @var{pos} is mapped into
   0 of the conditional distribution and the derivative is taken from
   the function PDF(@var{pos}+t*@var{dir}) w.r.t. @i{t}.
   On the other hand, when the coordinate @var{k} is used (i.e., when
   @var{dir} is set to NULL), the full conditional distribution of the
   distribution is considered (as used for the Gibbs sampler).
   In particular, the current point is just projected into the
   one-dimensional subspace without mapping it into the point 0.
*/

int unur_distr_condi_get_condition( struct unur_distr *distribution, const double **pos, const double **dir, int *k );
/* 
   Get condition for conditional @var{distribution}. 
   The values for the fixed variables are stored in @var{pos}, which
   must be a pointer to an array of size @code{dim}.
   The condition is stored in @var{dir} and @var{k}, respectively.

   @emph{Important:} Do @strong{not} change the entries in @var{pos}
   and @var{dir}!
*/

const UNUR_DISTR *unur_distr_condi_get_distribution( const UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/

/* =END */

/*---------------------------------------------------------------------------*/
