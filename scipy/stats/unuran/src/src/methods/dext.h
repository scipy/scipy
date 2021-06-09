/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dext.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DEXT                               *
 *         (wrapper for Continuous EXTernal generators)                      *
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
   =METHOD  DEXT   wrapper for Discrete EXTernal generators

   =UP  Methods_for_DISCR

   =REQUIRED routine for sampling discrete random variates

   =SPEED depends on external generator

   =REINIT supported

   =DESCRIPTION
      Method DEXT is a wrapper for external generators for discrete
      univariate distributions. It allows the usage of external
      random variate generators within the UNU.RAN framework.

   =HOWTOUSE
      The following steps are required to use some external generator
      within the UNU.RAN framework (some of these are optional):

      @enumerate
      @item
      Make an empty generator object using a unur_dext_new() call.
      The argument @var{distribution} is optional and can be replaced
      by NULL. However, it is required if you want to pass
      parameters of the generated distribution to the external
      generator or for running some validation tests provided by
      UNU.RAN.

      @item
      Create an initialization routine of type 
      @code{int (*init)(UNUR_GEN *gen)} and plug it into the generator
      object using the unur_dext_set_init() call. Notice that the
      @var{init} routine must return @code{UNUR_SUCCESS} when it has
      been executed successfully and @code{UNUR_FAILURE} otherwise.
      It is possible to get the size of and the pointer to the array
      of parameters of the underlying distribution object by the
      respective calls unur_dext_get_ndistrparams() and
      unur_dext_get_distrparams().
      Parameters for the external generator that are computed in the 
      @var{init} routine can be stored in a single array or structure
      which is available by the unur_dext_get_params() call.

      Using an @var{init} routine is optional and can be omitted.

      @item
      Create a sampling routine of type 
      @code{int (*sample)(UNUR_GEN *gen)} and plug it into the
      generator object using the unur_dext_set_sample() call.

      Uniform random numbers are provided by the unur_sample_urng()
      call. Do not use your own implementation of a uniform random
      number generator directly. If you want to use your own random
      number generator we recommend to use the UNU.RAN interface (see
      @pxref{URNG,,Using uniform random number generators}). 

      The array or structure that contains parameters for the external
      generator that are computed in the @var{init} routine are
      available using the unur_dext_get_params() call.

      Using a @var{sample} routine is of course obligatory.
      @end enumerate

      It is possible to change the parameters and the domain of the
      chosen distribution and run unur_reinit() to reinitialize the
      generator object. The @var{init} routine is then called again.

      Here is a short example that demonstrates the application of
      this method by means of the geometric distribution:

      @smallexample
      @include ref_example_dext.texi
      @end smallexample

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dext_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. 
*/

/*...........................................................................*/

int unur_dext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
/*
   Set initialization routine for external generator. Inside the

   @emph{Important:} The routine @var{init} must return
   @code{UNUR_SUCCESS} when the generator was initialized successfully
   and @code{UNUR_FAILURE} otherwise.

   Parameters that are computed in the @var{init} routine can be
   stored in an array or structure that is avaiable by means of the 
   unur_dext_get_params() call. Parameters of the underlying
   distribution object can be obtained by the
   unur_dext_get_distrparams() call.
*/

int unur_dext_set_sample( UNUR_PAR *parameters, int (*sample)(UNUR_GEN *gen) );
/*
   Set sampling routine for external generator.

   @emph{Important:}
   Use @code{unur_sample_urng(gen)} to get a uniform random number.
   The pointer to the array or structure that contains the parameters
   that are precomputed in the @var{init} routine are available by
   @code{unur_dext_get_params(gen,0)}.
   Additionally one can use the unur_dext_get_distrparams() call.
*/

void *unur_dext_get_params( UNUR_GEN *generator, size_t size );
/*
   Get pointer to memory block for storing parameters of external
   generator. A memory block of size @var{size} is automatically (re-)
   allocated if necessary and the pointer to this block is stored in
   the @var{generator} object. If one only needs the pointer to this
   memory block set @var{size} to @code{0}.

   Notice, that @var{size} is the size of the memory block and not the
   length of an array.

   @emph{Important:} This rountine should only be used in the
   initialization and sampling routine of the external generator. 
*/

double *unur_dext_get_distrparams( UNUR_GEN *generator );
/* */

int unur_dext_get_ndistrparams( UNUR_GEN *generator );
/*
   Get size of and pointer to array of parameters of underlying
   distribution in @var{generator} object.

   @emph{Important:} These rountines should only be used in the
   initialization and sampling routine of the external generator. 
*/

/*...........................................................................*/


/* =END */
/*---------------------------------------------------------------------------*/
