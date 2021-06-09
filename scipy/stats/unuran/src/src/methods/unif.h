/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unif.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method UNIF                               *
 *         (passes UNIForm random numbers through UNU.RAN framework;         *
 *         for testing only)                                                 *
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
   =METHOD  UNIF  wrapper for UNIForm random number generator

   =UP  Methods_for_UNID

   =REINIT supported

   =DESCRIPTION
      UNIF is a simple wrapper that makes it possible to use a uniform
      random number generator as a UNU.RAN generator. There are no
      parameters for this method.

   =HOWTOUSE
      Create a generator object with NULL as argument. The created generator
      object returns raw random numbers from the underlying uniform 
      random number generator.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_unif_new( const UNUR_DISTR *dummy );
/* 
   Get default parameters for generator.                                     
   UNIF does not need a distribution object. @var{dummy} is not used and
   can (should) be set to NULL. It is used to keep the API consistent.
*/

/* =END */

/*---------------------------------------------------------------------------*/












