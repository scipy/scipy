/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         generator objects.                                                *
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
   =NODE  Methods_all  Routines for all generator objects

   =UP Methods [0]

   =DESCRIPTION
      Routines for all generator objects.

   =END
*/

/*---------------------------------------------------------------------------*/
/* (Re-) Initialize generators                                               */

/* =ROUTINES */

UNUR_GEN *unur_init( UNUR_PAR *parameters );
/*
  Initialize a generator object. All necessary information must be
  stored in the parameter object.

  @strong{Important:} If an error has occurred a NULL pointer is
  return. This must not be used for the sampling routines (this causes a
  segmentation fault). 

  @strong{Always} check whether the call was successful or not!

  @emph{Important:} This call destroys the @var{parameter} object
  automatically. Thus it is not necessary/allowed to free it.
*/

/*...........................................................................*/

int unur_reinit( UNUR_GEN *generator );
/*
   Update an existing generator object after the underlying
   distribution has been modified (using unur_get_distr() together
   with corresponding set calls.
   It @strong{must} be executed before sampling using this generator
   object is continued as otherwise it produces an invalid sample or
   might even cause a segmentation fault.

   @emph{Important}: Currently not all methods allow reinitialization,
   see the description of the particular method (keyword @i{Reinit}).

   @emph{Important}: Reinitialization of the generator object might fail.
   Thus one @strong{must} check the return code:

   @table @asis 
   @item @code{UNUR_SUCCESS (0x0u)}
   success (no error)
   @item @code{UNUR_ERR_NO_REINIT}
   reinit routine not implemented.
   @item other values
   some error has occured while trying to reinitialize the generator
   object. 
   @end table

   @emph{Important}: When reinitialization fails then sampling routines
   always return INFINITY (for continuous distributions) or @code{0} (for
   discrete distributions), respectively.
   However, it is still possible to change the underlying distribution
   and try to reinitialize again.

   @emph{Important}: When one tries to run unur_reinit(), but reinitialization
   is not implemented, then the generator object cannot be used any more
   and must be destroyed and a new one has to be built from scratch.
   
*/

/*---------------------------------------------------------------------------*/
/* Sample from generator                                                     */

int    unur_sample_discr(UNUR_GEN *generator);
/* */

double unur_sample_cont(UNUR_GEN *generator);
/* */

int    unur_sample_vec(UNUR_GEN *generator, double *vector);
/* */

int    unur_sample_matr(UNUR_GEN *generator, double *matrix);
/*
   Sample from generator object. The three routines depend on the type
   of the generator object (discrete or continuous univariate
   distribution, multivariate distribution, or random matrix).
   
   @emph{Notice:} UNU.RAN uses arrays of @code{double}s to handle
   matrices. There the rows of the matrix are stored consecutively.
   
   @emph{Notice:} The routines unur_sample_vec() and unur_sample_matr()
   return @code{UNUR_SUCCESS} if generation was successful and 
   some error code otherwise.
   
   @strong{Important:} These routines do @strong{not} check whether 
   @var{generator} is an invalid NULL pointer.
*/

/*---------------------------------------------------------------------------*/
/* Estimate quantiles                                                        */

double unur_quantile ( UNUR_GEN *generator, double U );
/* 
   Compute the @var{U} quantile of a continuous distribution using a
   @var{generator} object that implements an (approximate) inversion
   methods.

   The following methods are currently available:
   @itemize
   @item 
   HINV, @pxref{HINV}.
   @item 
   NINV, @pxref{NINV}.
   @item
   PINV, @pxref{PINV}.
   @item
   CSTD, @pxref{CSTD}. @*
   This requires that @var{generator} implements an inversion method.
   @item
   DGT, @pxref{DGT}. @*
   The return value is (of course) type casted to @code{double}.
   @end itemize

   @strong{Important:} This routine does @strong{not} check whether 
   @var{generator} is an invalid NULL pointer.

   In case of an error UNUR_INFINITY or INT_MAX 
   (depending on the type of @var{generator}) is returned.
*/

/*---------------------------------------------------------------------------*/
/* Destroy (free) generator object                                           */

void  unur_free( UNUR_GEN *generator );
/*
  Destroy (free) the given generator object.
*/

/*---------------------------------------------------------------------------*/
/* Get information about generator object                                    */

const char *unur_gen_info( UNUR_GEN *generator, int help );
/* 
   Get a string with informations about the given @var{generator}.
   These informations allow some fine tuning of the generation method. 
   If @var{help} is TRUE, some hints on setting parameters are given.

   This function is intented for using in interactive environments 
   (like @t{R}).

   If an error occurs, then NULL is returned.
*/

/*---------------------------------------------------------------------------*/
/* Get dimension of generator for (multivariate) distribution                */

int unur_get_dimension( const UNUR_GEN *generator );
/*
  Get the number of dimension of a (multivariate) distribution.
  For a univariate distribution @code{1} is return.
*/

/*---------------------------------------------------------------------------*/

const char *unur_get_genid( const UNUR_GEN *generator );
/*
  Get identifier string for generator.
*/

/*---------------------------------------------------------------------------*/

UNUR_DISTR *unur_get_distr( const UNUR_GEN *generator );
/* 
   Get pointer to distribution object from generator object. 
   This function can be used to change the parameters of the distribution
   and reinitialize the generator object.
   Notice that currently @strong{not all} generating methods have a
   reinitialize routine. 
   This function should be used with extreme care. Changing the distribution
   is changed and using the generator object without reinitializing
   might cause wrong samples or segmentation faults.
   Moreover, if the corresponding generator object is freed, the
   pointer must not be used.

   @strong{Important:} The returned distribution object must not 
   be freed. If the distribution object is changed then one @strong{must}
   run unur_reinit()!
*/

/*---------------------------------------------------------------------------*/

int unur_set_use_distr_privatecopy( UNUR_PAR *parameters, int use_privatecopy );
/* 
   Set flag whether the generator object should make a private copy of
   the given distribution object or just stores the pointer to this
   distribution object. Values for @var{use_privatecopy}:
   @table @code
   @item TRUE
   make a private copy (default)
   @item FALSE
   do not make a private copy and store pointer to given (external)
   distribution object.
   @end table

   By default, generator objects keep their own private copy of the
   given distribution object. Thus the generator object can be handled
   independently from other UNU.RAN objects (with uniform random number
   generators as the only exception). When the generator object is
   initialized the given distribution object is cloned and stored.

   However, in some rare situations it can be useful when only the
   pointer to the given distribution object is stored without making a
   private copy. A possible example is when only one random variate has
   to be drawn from the distribution.
   This behavior can be achieved when @var{use_localcopy} is set to
   FALSE.

   @strong{Warning!}
   Using a pointer to the external distribution object instead of a
   private copy must be done with @strong{extreme care}!
   When the distrubtion object is changed or freed then the generator
   object does not work any more, might case a segmentation fault, or
   (even worse) produces garbage. 
   On the other hand, when the generator object is initialized or used
   to draw a random sampling the distribution object may be changed.

   @emph{Notice:}
   The prototypes of all @code{unur_<method>_new} calls use a 
   @code{const} qualifier for the distribution argument.
   However, if @var{use_privatecopy} is set to FALSE this qualifier is
   discarded and the distribution might be changed.

   @strong{Important!}
   If @var{use_localcopy} is set to FALSE and the corresponding
   distribution object is changed then one must run unur_reinit()
   on the generator object.
   (Notice that currently not all generation methods support
   reinitialization.)

   Default: @var{use_privatecopy} is TRUE.

*/

/* =END */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

UNUR_GEN *unur_gen_clone( const UNUR_GEN *gen );
void unur_par_free( UNUR_PAR *par);

/*---------------------------------------------------------------------------*/
