/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_distr.h                                                *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes and macros for deprecated routines            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   THESE ROUTINES SHOULD NOT BE USED ANY MORE!                             *
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
#ifndef UNUR_DEPRECATED_DISTR_H_SEEN
#define UNUR_DEPRECATED_DISTR_H_SEEN
/*---------------------------------------------------------------------------*/
/* Set and get standardized marginal distributions.                          */

int unur_distr_cvec_set_stdmarginals( UNUR_DISTR *distribution, UNUR_DISTR *marginal );
/* 
   Sets standardized marginal distributions
   of the given @var{distribution} to the same @var{marginal}
   distribution object. The @var{marginal} distribution must be an
   instance of a continuous univariate distribution object.
   In conjunction with unur_distr_cvec_set_covar()
   and unur_distr_cvec_set_mean() the standardized marginals must be
   used, i.e., they should have mean 0 and standard deviation 1 
   (if both exist for the given marginal distribution).
   Notice that the marginal distribution is copied into the
   @var{distribution} object.
*/

int unur_distr_cvec_set_stdmarginal_array( UNUR_DISTR *distribution, UNUR_DISTR **marginals );
/* 
   Analogously to the above unur_distr_cvec_set_stdmarginals() calls.
   However, now an array @var{marginals} of the pointers to each of
   the marginal distributions must be given. It @strong{must} be an
   array of size @code{dim}, where @code{dim} is the dimension
   returned by unur_distr_get_dim(). 

   @emph{Notice}: Local copies for each of the entries are stored in
   the @var{distribution} object. If some of these entries are
   identical (i.e. contain the same pointer), then for each of these a
   new copy is made.
*/

int unur_distr_cvec_set_stdmarginal_list( UNUR_DISTR *distribution, ... );
/* 
   Similar to the above unur_distr_cvec_set_stdmarginal_array() call.
   However, now the pointers to the particular marginal distributions
   can be given as parameter and does not require an array of
   pointers. Additionally the given distribution objects are
   immediately destroyed. Thus calls like unur_distr_normal() can be 
   used as arguments. 
   (With unur_distr_cvec_set_marginal_array() the result of such call
   has to be stored in a pointer since it has to be freed afterwarts
   to avoid memory leaks!)

   If one of the given pointer to marginal distributions is the NULL
   pointer then the marginal distributions of @var{distribution} are
   not set (or previous settings are not changed) and an error code is
   returned.

   @strong{Important:} All distribution objects given in the argument
   list are destroyed!
*/

const UNUR_DISTR *unur_distr_cvec_get_stdmarginal( const UNUR_DISTR *distribution, int n );
/* 
   Get pointer to the @var{n}-th (standardized) marginal distribution
   object from the given multivariate @var{distribution}. 
   If this does not exist, NULL is returned. 
   The marginal distributions are enumerated from @code{1}
   to @code{dim}, where @code{dim} is the dimension
   returned by unur_distr_get_dim(). 
*/

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_DEPRECATED_DISTR_H_SEEN */
/*---------------------------------------------------------------------------*/
