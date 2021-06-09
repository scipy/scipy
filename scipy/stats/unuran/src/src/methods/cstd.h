/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cstd.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method CSTD                               *
 *         (wrapper for special generators for                               *
 *         Continuous STanDard distributions)                                *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/* 
   =METHOD  CSTD   Continuous STandarD distributions

   =UP  Methods_for_CONT

   =REQUIRED standard distribution from UNU.RAN library
      (@pxref{Stddist,,Standard distributions}) or continuous
      distribution with inverse CDF.

   =SPEED Set-up: fast, Sampling: depends on distribution and generator

   =REINIT supported

   =DESCRIPTION
      CSTD is a wrapper for special generators for continuous
      univariate standard distributions. It only works for
      distributions in the UNU.RAN library of standard distributions
      (@pxref{Stddist,,Standard distributions})
      or for continuous distributions where the inverse CDF is given.
      If a distribution object is provided that is build from scratch,
      it must provide the inverse CDF. Then CSTD implements the
      inversion method. Otherwise, the NULL pointer is returned.

      For some distributions more than one special generator
      is possible. 

   =HOWTOUSE
      Create a distribution object for a standard distribution
      from the UNU.RAN library 
      (@pxref{Stddist,,Standard distributions}),
      or create a continuous distribution object and set the function
      for the inverse CDF using unur_distr_cont_set_invcdf().
      For some distributions more than one special generator
      (@emph{variants}) is possible. These can be choosen by a
      unur_cstd_set_variant() call. For possible variants 
      @pxref{Stddist,,Standard distributions}.
      However the following are common to all distributions:

      @table @code
      @item UNUR_STDGEN_DEFAULT
      the default generator.                      
      @item UNUR_STDGEN_FAST
      the fastest available special generator.
      @item UNUR_STDGEN_INVERSION
      the inversion method (if available).
      @end table
      
      Notice that the variant @code{UNUR_STDGEN_FAST} for a special
      generator may be slower than one of the universal algorithms!
      Additional variants may exist for particular distributions.
      
      Sampling from truncated distributions (which can be constructed by 
      changing the default domain of a distribution by means of
      unur_distr_cont_set_domain() or unur_cstd_chg_truncated() calls)
      is possible but requires the inversion method. Moreover the CDF
      of the distribution must be implemented.
   
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Common variants for all special generators                                */

#define UNUR_STDGEN_DEFAULT   0        /* default algorithm (don't change 0!)*/
#define UNUR_STDGEN_INVERSION (~0u)    /* inversion method                   */
#define UNUR_STDGEN_FAST      (0)      /* fastest algorithm                  */

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_cstd_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. It requires a distribution object 
   for a continuous univariant distribution from the 
   UNU.RAN library of standard distributions 
   (@pxref{Stddist,,Standard distributions}).

   Using a truncated distribution is allowed only if the inversion method
   is available and selected by the unur_cstd_set_variant() call immediately 
   after creating the parameter object.
   Use a unur_distr_cont_set_domain() call to get a truncated distribution.
   To change the domain of a (truncated) distribution of a generator use the
   unur_cstd_chg_truncated() call.
*/

/*...........................................................................*/

int unur_cstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/*
   Set variant (special generator) for sampling from a given distribution.
   For possible variants 
   @pxref{Stddist,,Standard distributions}.

   Common variants are @code{UNUR_STDGEN_DEFAULT} for the default generator,
   @code{UNUR_STDGEN_FAST} for (one of the) fastest implemented
   special generators, and @code{UNUR_STDGEN_INVERSION} for the
   inversion method (if available). 
   If the selected variant number is not implemented, then an error code is
   returned and the variant is not changed.
*/

/*...........................................................................*/

int unur_cstd_chg_truncated( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the (truncated) distribution.
   This is only possible if the inversion method is used.
   Otherwise this call has no effect and an error code is returned.

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call.

   It is not required to run unur_reinit() after this call has been used.

   @emph{Important:} If the CDF is (almost) the same for @var{left} and 
   @var{right} and (almost) equal to @code{0} or @code{1}, then the truncated 
   domain is not chanced and the call returns an error code.

   @emph{Notice:} If the parameters of the distribution has been changed
   it is recommended to set the truncated domain again, since the
   former call might change the domain of the distribution but not
   update the values for the boundaries of the truncated
   distribution.
*/

/* =END */
/*---------------------------------------------------------------------------*/

/* Yet not documented! */

double unur_cstd_eval_invcdf( const UNUR_GEN *generator, double u );
/*
   Evaluate inverse CDF at @var{u}. However, this requires that
   @var{generator} implements an inversion method.
   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).

   @emph{Notice}: When the domain has been truncated by a  
   unur_cstd_chg_truncated() call then the inverse CDF of the
   truncated distribution is returned.
*/

/*---------------------------------------------------------------------------*/

