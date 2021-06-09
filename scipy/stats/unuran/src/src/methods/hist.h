/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hist.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HIST                               *
 *         (HISTogram of empirical distribution)                             *
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

/* 
   =METHOD  HIST   HISTogramm of empirical distribution

   =UP  Methods_for_CEMP

   =REQUIRED histogram 

   =SPEED Set-up: moderate,
          Sampling: fast

   =REINIT not implemented

   =DESCRIPTION
      Method HIST generates random variates from an empirical distribution
      that is given as histogram. Sampling is done using the inversion
      method. 

      If observed (raw) data are provided we recommend method EMPK
      (@pxref{EMPK,,EMPirical distribution with Kernel smoothing})
      instead of compting a histogram as this reduces information.

   =HOWTOUSE
      Method HIST uses empirical distributions that are given as a
      histgram. There are no optional parameters.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hist_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

/* =END */

/*---------------------------------------------------------------------------*/

