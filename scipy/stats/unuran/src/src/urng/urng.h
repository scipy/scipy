/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros and function prototypes for using uniform         *
 *         random number generators inside UNU.RAN.                          *
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
#ifndef URNG_H_SEEN
#define URNG_H_SEEN
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG  Using uniform random number generators

   =UP TOP [50]

   =DESCRIPTION
      UNU.RAN is designed to work with many sources of (pseudo-) random
      numbers or low discrepancy numbers (so called quasi-random
      numbers) for almost all tasks in discrete event simulation,
      (quasi-) Monte Carlo integration or any other stochastic
      methods. Hence UNU.RAN uses pointers to access uniform (pseudo-)
      random number generators (URNG). 

      Each UNU.RAN (non-uniform random variate) generator object has a
      pointer to a URNG object. Thus each UNU.RAN generator object may
      have its own (independent) URNG or several generator objects can
      share the same URNG.

      If no URNG is provided for a parameter or generator object a default
      generator is used which is the same for all generators. This URNG is
      defined in @file{unuran_config.h} at compile time and can be
      changed at runtime.

      UNU.RAN uses a unified interface for all sources of random numbers.
      Unfortunately, the API for random number generators, like the
      @file{GSL} (GNU Scientific Library), Otmar Lendl's @file{prng}
      (Pseudo random number generators), or a single function
      implemented by the user herself, are quite different. 
      Hence an object of type @code{UNUR_URNG} is introduced to store
      the URNG. Thus it is possible to handle different sources of
      such URNGs with the unified API. It is inspired from similar to
      Pierre L'Ecuyers @file{RngStreams} library: 

      @itemize @minus
      @item seed the random number generator;
      @item get a uniform random number;
      @item reset the URNG;
      @item skip to the begining next substream;
      @item sample antithetic numbers;
      @item delete the URNG object.
      @end itemize

      The routine to create a URNG depends on the chosen random number
      generator (i.e. library). Nevertheless, there exist wrapper
      functions to simplify this task.

      Currently the following sources of uniform random numbers are
      directly supported (i.e., there exist wrapper functions). 
      Of course other random number generation libraries can be used.

      @enumerate

      @item
      @code{FVOID}
      
      URNGs of type @code{double uniform(void *state)}.
      The argument @var{state} can be simply ignored in the
      implementation of @code{uniform} when a global state variable is
      used.
      UNU.RAN contains some build-in URNGs of this type in directory
      @file{src/uniform/}.
      
      @item
      @code{PRNG}
      
      URNGs from Otmar Lendl's @code{prng} library. It provides a very
      flexible way to sample form arbitrary URNGs by means of an object
      oriented programing paradigma. Similarly to the UNU.RAN library
      independent generator objects can be build and used.
      
      This library has been developed by the pLab group at the university
      of Salzburg (Austria, EU) and implemented by Otmar Lendl.
      It is available from
      @uref{http://statmath.wu.ac.at/prng/}
      or from the pLab site at
      @uref{http://random.mat.sbg.ac.at/}.
      
      This interface must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-prng}.

      @item
      @code{RNGSTREAM}

      Pierre L'Ecuyer's @code{RngStream} library for multiple 
      independent streams of pseudo-random numbers. 
      A GNU-style package is available from
      @uref{http://statmath.wu.ac.at/software/RngStreams/}.

      This interface must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-rngstream}.

      @item
      @code{GSL}

      URNG from the GNU Scientific Library (GSL).
      It is available from
      @uref{http://www.gnu.org/software/gsl/}.                               

      This interface must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-gsl}.

      @end enumerate

   =HOWTOUSE
      Each UNU.RAN generator object has a pointer to a uniform
      (pseudo-) random number generator (URNG). It can be set via the
      unur_set_urng() call. It is also possible to read this pointer
      via unur_get_urng() or change the URNG for an existing generator
      object by means of unur_chg_urng().

      It is important to note that these calls only copy the pointer
      to the URNG object into the generator object.

      If no URNG is provided for a parameter or generator object a default
      URNG is used which is the same for all generators. This URNG is
      defined in @file{unuran_config.h} at compile time. A pointer to
      this default URNG can be obtained via unur_get_default_urng().
      Nevertheless, it is also possible to change this default URNG by
      another one at runtime by means of the unur_set_default_urng()
      call. However, this only takes effect for new parameter objects.

      Some generating methods provide the possibility of correlation
      induction. For this feature a second auxiliary URNG is required.
      It can be set and changed by unur_set_urng_aux() and
      unur_chg_urng_aux() calls, respectively. Since the auxiliary
      URNG is by default the same as the main URNG, the
      auxiliary URNG must be set after any unur_set_urng() or
      unur_chg_urng() call! Since in special cases mixing of two URNG
      might cause problems, we supply a default auxiliary generator
      that can be used by a unur_use_urng_aux_default() call (after
      the main URNG has been set). This default auxiliary generator
      can be changed with analogous calls as the (main) default
      uniform generator. 

      Uniform random number generators form different sources have
      different programming interfaces. Thus UNU.RAN stores all
      information about a particular uniform random number generator
      in a structure of type @code{UNUR_URNG}. Before a URNG can be
      used with UNU.RAN an appropriate object has to be created ba a
      unur_urng_new() call. 
      This call takes two arguments: the pointer to the sampling
      routine of the generator and a pointer to a possible argument
      that stores the state of the generator. The function must be of
      type @code{double (*sampleunif)(void *params)}, but functions
      without any argument also work.
      Additionally one can set pointers to functions for reseting or
      jumping the streams generated by the URNG by the corresponding
      @code{set} calls.

      UNU.RAN provides a unified API to all sources of random numbers.
      Notice, however, that not all functions work for all random
      number generators (as the respective library has not implemented
      the corresponding feature).

      There are wrapper functions for some libraries of uniform random
      number generators to simplify the task of creating a UNU.RAN
      object for URNGs.
      These functions must be compiled into UNU.RAN using the
      corresponding configure flags (see description of the respective
      interface below).

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subheading Set and get default uniform RNGs
*/

UNUR_URNG *unur_get_default_urng( void );
/*
  Get the pointer to the default URNG. The default URNG is used by all
  generators where no URNG was set explicitly by a unur_set_urng()
  call.
*/

UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );
/*
  Change the default URNG that is used for new parameter objects.
  It returns the pointer to the old default URNG that has been used.
*/


UNUR_URNG *unur_set_default_urng_aux( UNUR_URNG *urng_new );
/* */

UNUR_URNG *unur_get_default_urng_aux( void );
/*
  Analogous calls for default auxiliary generator.
*/

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subheading Set, change and get uniform RNGs in generator objects
*/

int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
/*
  Use the URNG @code{urng} for the new generator. This overrides the
  default URNG. It also sets the auxiliary URNG to @code{urng}.

  @emph{Important}: For multivariate distributions that use 
  marginal distributions this call does not work properly.
  It is then better first to create the generator object (by
  a unur_init() call) and then change the URNG by means of 
  unur_chg_urng().
*/

UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
/*
  Change the URNG for the given generator. It returns the pointer to
  the old URNG that has been used by the generator.
  It also changes the auxiliary URNG to @code{urng} and thus it
  overrides the last unur_chg_urng_aux() call.
*/

UNUR_URNG *unur_get_urng( UNUR_GEN *generator );
/*
  Get the pointer to the URNG that is used by the @var{generator}.
  This is usefull if two generators should share the same URNG.
*/

int unur_set_urng_aux( UNUR_PAR *parameters, UNUR_URNG *urng_aux );
/*
  Use the auxiliary URNG @code{urng_aux} for the new generator. 
  (Default is the default URNG or the URNG from the last
  unur_set_urng() call. Thus if the auxiliary generator should be
  different to the main URNG, unur_set_urng_aux() must be called after
  unur_set_urng(). 
  The auxiliary URNG is used as second stream of uniform random
  number for correlation induction.
  It is not possible to set an auxiliary URNG for a method that does
  not need one. In this case an error code is returned.
*/

int unur_use_urng_aux_default( UNUR_PAR *parameters );
/* 
   Use the default auxiliary URNG.
   (It must be set after unur_get_urng().)
   It is not possible to set an auxiliary URNG for a method that does
   not use one (i.e. the call returns an error code).
*/

int unur_chgto_urng_aux_default( UNUR_GEN *generator );
/*
   Switch to default auxiliary URNG.
   (It must be set after unur_get_urng().)
   It is not possible to set an auxiliary URNG for a method that does
   not use one (i.e. the call returns an error code).
*/

UNUR_URNG *unur_chg_urng_aux( UNUR_GEN *generator, UNUR_URNG *urng_aux );
/*
  Change the auxiliary URNG for the given @var{generator}. It returns
  the pointer to the old auxiliary URNG that has been used by the
  generator. It has to be called after each unur_chg_urng() when the 
  auxiliary URNG should be different from the main URNG.
  It is not possible to change the auxiliary URNG for a method that
  does not use one (i.e. the call NULL).
*/

UNUR_URNG *unur_get_urng_aux( UNUR_GEN *generator );
/*
  Get the pointer to the auxiliary URNG that is used by the
  @var{generator}. This is usefull if two generators should share the same
  URNG.
*/

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subheading Handle uniform RNGs

   @emph{Notice:} Some of the below function calls do not work for
   every source of random numbers since not every library has
   implemented these features.

*/

double unur_urng_sample (UNUR_URNG *urng);
/*
   Get a uniform random number from @var{urng}.
   If the NULL pointer is given, the default uniform generator is
   used. 
*/

double unur_sample_urng (UNUR_GEN *gen);
/* 
   Get a uniform random number from the underlying uniform
   random number generator of generator @var{gen}.
   If the NULL pointer is given, the default uniform generator is
   used. 
*/

int unur_urng_sample_array (UNUR_URNG *urng, double *X, int dim);
/*
   Set array @var{X} of length @var{dim} with uniform random numbers 
   sampled from generator @var{urng}. If @var{urng} is the NULL
   pointer, the default uniform generator is used.

   @emph{Important:} 
   If @var{urng} is based on a point set generator (this is the case
   for generators of low discrepance point sets as used in quasi-Monte
   Carlo methods) it has a ``natural dimension'' @i{s}. 
   In this case either only the first @i{s} entries of @var{X} are
   filled (if @i{s} < @var{dim}), or the first @var{dim} coordinates
   of the generated point are filled. 

   The called returns the actual number of entries filled. In case of
   an error @code{0} is returned.
*/

int unur_urng_reset (UNUR_URNG *urng);
/*
   Reset @var{urng} object. 
   The routine tries two ways to reset the generator (in this order):

   @enumerate
   @item
     It uses the reset function given by an unur_urng_set_reset()
     call. 

   @item
     It uses the seed given by the last unur_urng_seed() call (which
     requires a seeding function given by a unur_urng_set_seed()
     call). 
   @end enumerate

   If neither of the two methods work resetting of the generator is
   not possible and an error code is returned.

   If the NULL pointer is given, the default uniform generator is
   reset.
*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_URNG_UNURAN
/*---------------------------------------------------------------------------*/

int unur_urng_sync (UNUR_URNG *urng);
/* 
   Jump into defined state ("sync") of the generator. This is useful
   when point generators are used where the coordinates are 
   sampled via unur_urng_sample(). Then this call can be used to 
   jump to the first coordinate of the next generated point.
*/

int unur_urng_seed (UNUR_URNG *urng, unsigned long seed);
/*
   Set @var{seed} for generator @var{urng}.
   It returns an error code if this is not possible for the given
   URNG. If the NULL pointer is given, the default uniform generator is
   seeded (if possible).

   @emph{Notice}: Seeding should be done only once for a particular
   generator (except for resetting it to the initial state).
   Expertise is required when multiple seeds are used to get independent
   streams. Thus we recommend appropriate libraries for this task,
   e.g. Pierre L'Ecuyer's @file{RngStreams} package. For this library
   only a package seed can be set and thus the unur_urng_seed() call
   will not have any effect to generators of this type. Use
   unur_urng_reset() or unur_urng_rngstream_new() instead, depending
   whether one wants to reset the stream or get a new stream that is
   independent from the previous ones.
*/

int unur_urng_anti (UNUR_URNG *urng, int anti);
/*
   Switch to antithetic random numbers in @var{urng}.
   It returns an error code if this is not possible for the given
   URNG.

   If the NULL pointer is given, the antithetic flag of the default
   uniform generator is switched (if possible).
*/

int unur_urng_nextsub (UNUR_URNG *urng);
/*
   Jump to start of the next substream of @var{urng}.
   It returns an error code if this is not possible for the given
   URNG.

   If the NULL pointer is given, the default uniform generator is set
   to the start of the next substream (if possible).
*/

int unur_urng_resetsub (UNUR_URNG *urng);
/*
   Jump to start of the current substream of @var{urng}.
   It returns an error code if this is not possible for the given
   URNG.

   If the NULL pointer is given, the default uniform generator is set
   to the start of the current substream (if possible).
*/

int unur_gen_sync (UNUR_GEN *generator);
/* */

int unur_gen_seed (UNUR_GEN *generator, unsigned long seed);
/* */

int unur_gen_anti (UNUR_GEN *generator, int anti);
/* */

int unur_gen_reset (UNUR_GEN *generator);
/* */

int unur_gen_nextsub (UNUR_GEN *generator);
/* */

int unur_gen_resetsub (UNUR_GEN *generator);
/* 
   Analogous to unur_urng_sync(), unur_urng_seed(), unur_urng_anti(), 
   unur_urng_reset(), unur_urng_nextsub(), and unur_urng_resetsub(),
   but act on the URNG object used by the @var{generator} object.

   @emph{Warning:} These calls should be used with care as it
   influences all generator objects that share the same URNG object!
*/

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subheading API to create a new URNG object

   @emph{Notice:} These functions are provided to built a 
   UNUR_URNG object for a particular external random number
   generator from scratch. For some libraries that contain random
   number generators (like the GSL) there are special calls,
   e.g. unur_urng_gsl_new(), to get such an object. Then there is no
   need to change the UNUR_URNG object as it already contains all
   available features. 

   If you have a particular library for random number generators you
   can either write wrapper function like those in
   @file{src/uniform/urng_gsl.c} or write an email to the authors of
   UNU.RAN to write it for you.
*/

UNUR_URNG *unur_urng_new( double (*sampleunif)(void *state), void *state );
/*
   Get a new URNG object. 
   @var{sampleunif} is a function to the uniform sampling routine,
   @var{state} a pointer to its arguments which usually contains the
   state variables of the generator.

   Functions @var{sampleunif} with a different type for @var{p} or
   without an argument at all also work. A typecast might be necessary
   to avoid compiler warnings or error messages.

   For functions @var{sampleunif} that does not have any argument
   should use NULL for @var{state}.
   
   @emph{Important:} @var{sampleunif} must not be the NULL pointer.

   There are appropriate calls that simplifies the task of creating
   URNG objects for some libraries with uniform random number
   generators, see below.
*/

void unur_urng_free (UNUR_URNG *urng);
/* 
   Destroy @var{urng} object.
   It returns an error code if this is not possible. 

   If the NULL is given, this function does nothing.

   @emph{Warning:} This call must be used with care. The @var{urng}
   object must not be used by any existing generator object!
   It is designed to work in conjunction with the wrapper functions
   to create URNG objects for generators of a particular library.
   Thus an object created by an unur_urng_prng_new() call can be
   simply destroyed by an unur_urng_free() call.
*/

int unur_urng_set_sample_array( UNUR_URNG *urng, unsigned int (*samplearray)(void *state, double *X, int dim) );
/*
   Set function to fill array @var{X} of length @var{dim} with random
   numbers generated by generator @var{urng} (if available).
*/

int unur_urng_set_sync( UNUR_URNG *urng, void (*sync)(void *state) );
/* 
   Set function for jumping into a defined state (``sync'').
*/

int unur_urng_set_seed( UNUR_URNG *urng, void (*setseed)(void *state, unsigned long seed) );
/*
   Set function to seed generator @var{urng} (if available).
*/

int unur_urng_set_anti( UNUR_URNG *urng, void (*setanti)(void *state, int anti) );
/*
   Set function to switch the antithetic flag of generator @var{urng}
   (if available).
*/

int unur_urng_set_reset( UNUR_URNG *urng, void (*reset)(void *state) );
/* 
   Set function for reseting the uniform random number generator
   @var{urng} (if available).
*/

int unur_urng_set_nextsub( UNUR_URNG *urng, void (*nextsub)(void *state) );
/*
   Set function that allows jumping to start of the next substream of
   @var{urng} (if available).
*/

int unur_urng_set_resetsub( UNUR_URNG *urng, void (*resetsub)(void *state) );
/*
   Set function that allows jumping to start of the current substream
   of @var{urng} (if available).
*/

int unur_urng_set_delete( UNUR_URNG *urng, void (*fpdelete)(void *state) );
/*
   Set function for destroying @var{urng} (if available).
*/

/*---------------------------------------------------------------------------*/
#endif   /* end defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* URNG_H_SEEN */
/*---------------------------------------------------------------------------*/
