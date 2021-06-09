/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mvstd.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MVSTD                              *
 *         (wrapper for special generators for                               *
 *         MultiVariate continuous STandarD distributions)                   *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold             *
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
   =METHOD  MVSTD   MultiVariate continuous STandarD distributions

   =UP  Methods_for_CVEC

   =REQUIRED standard distribution from UNU.RAN library
      (@pxref{Stddist,,Standard distributions}).

   =SPEED depends on distribution and generator

   =REINIT supported

   =DESCRIPTION
      MVSTD is a wrapper for special generators for multivariate
      continuous standard distributions. It only works for
      distributions in the UNU.RAN library of standard distributions
      (@pxref{Stddist,,Standard distributions}).
      If a distribution object is provided that is build from scratch,
      or if no special generator for the given standard distribution is
      provided, the NULL pointer is returned.

   =HOWTOUSE
      Create a distribution object for a standard distribution
      from the UNU.RAN library (@pxref{Stddist,,Standard distributions}).
      
      Sampling from truncated distributions (which can be constructed by 
      changing the default domain of a distribution by means of
      unur_distr_cvec_set_domain_rect() call) is not possible.
   
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mvstd_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. It requires a distribution object 
   for a multivariate continuous distribution from the 
   UNU.RAN library of standard distributions 
   (@pxref{Stddist,,Standard distributions}).
   Using a truncated distribution is not possible.
*/

/* =END */
/*---------------------------------------------------------------------------*/
