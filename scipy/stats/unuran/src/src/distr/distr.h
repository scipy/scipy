/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects         *
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
   =NODEX  Distribution_objects  Handling distribution objects

   =UP TOP [30]

   =DESCRIPTION
      Objects of type @code{UNUR_DISTR} are used for handling
      distributions. All data about a distribution are stored in this
      object. UNU.RAN provides functions that return instances of such
      objects for standard distributions 
      (@pxref{Stddist,,Standard distributions}).
      It is then possible to change these distribution objects by
      various set calls. Moreover, it is possible to build a
      distribution object entirely from scratch. For this purpose
      there exists @command{unur_distr_<type>_new} calls that
      return an empty object of this type for each object type
      (eg. univariate contiuous) which can be filled with the
      appropriate set calls.

      UNU.RAN distinguishes between several types of distributions,
      each of which has its own sets of possible parameters (for
      details see the corresponding sections):
      @itemize @minus
      @item continuous univariate distributions
      @item continuous univariate order statistics
      @item continuous empirical univariate distributions
      @item continuous multivariate distributions
      @item continuous empirical multivariate distributions
      @item matrix distributions
      @item discrete univariate distributions
      @end itemize
      
      Notice that there are essential data about a distribution,
      eg. the PDF, a list of (shape, scale, location) parameters for
      the distribution, and the domain of (the possibly truncated)
      distribution. And there exist parameters that are/can be
      derived from these, eg. the mode of the distribution or the area
      below the given PDF (which need not be normalized for many
      methods). UNU.RAN keeps track of parameters which are
      known. Thus if one of the essential parameters is changed all
      derived parameters are marked as unknown and must be set again
      if these are required for the chosen generation method.
      Additionally to set calls there are calls for updating derived
      parameters for objects provided by the UNU.RAN library of standard
      distributions (one for each parameter to avoid computational
      overhead since not all parameters are required for all generator
      methods). 

      All parameters of distribution objects can be read by
      corresponding get calls.

      Every generator object has its own copy of a distribution object
      which is accessible by a unur_get_distr() call. Thus the
      parameter for this distribution can be read. However,
      @strong{never} extract the distribution object out of a
      generator object and run one of the set calls on it to modify
      the distribution.  (How should the poor generator object know
      what has happend?) Instead there exist calls for each of the
      generator methods that change particular parameters of the
      internal copy of the distribution object.

   =HOWTOUSE
      UNU.RAN collects all data required for a particular generation
      method in a @emph{distribution object}. There are two ways to
      get an instance of a distributions object: 
      @enumerate
      @item 
      Build a distribtion from scratch, by means of 
      the corresponding @command{unur_distr_<type>_new} call, 
      where @command{<type>} is the type of the distribution as
      listed in the below subsections.

      @item
      Use the corresponding @command{unur_distr_<name>_new} call
      to get prebuild distribution from the UNU.RAN library of standard
      distributions.
      Here @command{<name>} is the name of the
      standard distribution in @ref{Stddist,,Standard distributions}.
      @end enumerate

      In either cases the corresponding 
      @command{unur_distr_<type>_set_<param>} calls to set the
      necessary parameters @command{<param>} (case 1), or 
      change the values of the standard distribution in case 2 (if
      this makes sense for you). In the latter case @command{<type>}
      is the type to which the standard distribution belongs to.
      These @command{set} calls return @code{UNUR_SUCCESS} when the
      correspondig parameter has been set successfully. Otherwise an
      error code is returned.
      
      The parameters of a distribution are divided into
      @emph{essential} and @emph{derived} parameters.

      Notice, that there are some restrictions in setting parameters
      to avoid possible confusions.
      Changing essential parameters marks derived parameters as
      @code{unknown}. Some of the parameters cannot be changed any
      more when already set; some parameters block each others.
      In such a case a new instance of a distribution object has to be
      build. 

      Additionally @command{unur_distr_<type>_upd_<param>} calls can
      be used for updating derived parameters for objects provided by
      the UNU.RAN library of standard distributions.

      All parameters of a distribution object get be read by means of
      @command{unur_distr_<type>_get_<param>} calls.

      Every distribution object be identified by its @code{name} which
      is a string of arbitrary characters provided by the user. For
      standard distribution it is automatically set to
      @command{<name>} in the corresponding @command{new} call. It can
      be changed to any other string.

   =EON
*/

/*---------------------------------------------------------------------------*/

/*
   =NODEX   AllDistr   Functions for all kinds of distribution objects

   =UP Distribution_objects [05]

   =DESCRIPTION
     The calls in this section can be applied to all distribution
     objects. 

     @itemize @minus
     @item Destroy @command{free} an instance of a generator object.

     @item Ask for the @command{type} of a generator object.

     @item Ask for the @command{dimension} of a generator object.

     @item Deal with the @command{name} (identifier string) of a generator object. 
     @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x010u,     /* univariate continuous distribution       */ 
  UNUR_DISTR_CEMP  = 0x011u,     /* empirical univ. cont. distr. (a sample)  */ 
  UNUR_DISTR_CVEC  = 0x110u,     /* mulitvariate continuous distribution     */ 
  UNUR_DISTR_CVEMP = 0x111u,     /* empirical multiv. cont. distr. (sample)  */ 
  UNUR_DISTR_MATR  = 0x210u,     /* matrix distribution                      */ 
  UNUR_DISTR_DISCR = 0x020u      /* univariate discrete distribution         */ 
};

/*---------------------------------------------------------------------------*/

/*
  Parameters common to all distributions.
*/

/* =ROUTINES */

void unur_distr_free( UNUR_DISTR *distribution );
/* 
   Destroy the @var{distribution} object.
*/


int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
/* */

const char *unur_distr_get_name( const UNUR_DISTR *distribution );
/* 
   Set and get @var{name} of @var{distribution}. The @var{name} can be
   an arbitrary character string. It can be used to identify generator
   objects for the user. It is used by UNU.RAN when printing
   information of the distribution object into a log files.
*/


int unur_distr_get_dim( const UNUR_DISTR *distribution );
/* 
   Get number of components of a random vector (its dimension) the
   @var{distribution}. 

   For univariate distributions it returns dimension @code{1}.

   For matrix distributions it returns the number of components
   (i.e., number of rows times number of columns).
   When the respective numbers of rows and columns are needed use
   unur_distr_matr_get_dim() instead.
*/


unsigned int unur_distr_get_type( const UNUR_DISTR *distribution );
/* 
   Get type of @var{distribution}. 
   Possible types are
   @table @code
   @item UNUR_DISTR_CONT
   univariate continuous distribution
   @item UNUR_DISTR_CEMP
   empirical continuous univariate distribution (i.e. a sample)
   @item UNUR_DISTR_CVEC
   continuous mulitvariate distribution
   @item UNUR_DISTR_CVEMP
   empirical continuous multivariate distribution (i.e. a vector sample)
   @item UNUR_DISTR_DISCR
   discrete univariate distribution
   @item UNUR_DISTR_MATR
   matrix distribution
   @end table

   Alternatively the @command{unur_distr_is_<TYPE>}
   calls can be used.
*/

int unur_distr_is_cont( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a continuous univariate distribution.
*/

int unur_distr_is_cvec( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a continuous multivariate distribution.
*/

int unur_distr_is_cemp( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is an empirical continuous univariate distribution,
   i.e. a sample.
*/

int unur_distr_is_cvemp( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is an empirical continuous multivariate
   distribution.
*/

int unur_distr_is_discr( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a discrete univariate distribution.
*/

int unur_distr_is_matr( const UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a matrix distribution.
*/


int unur_distr_set_extobj( UNUR_DISTR *distribution, const void *extobj );
/* 
   Store a pointer to an external object. This might be usefull if
   the PDF, PMF, CDF or other functions used to implement a particular
   distribution a parameter set that cannot be stored as doubles
   (e.g. pointers to some structure that holds information of the distribution).
   
   @strong{Important:} 
   When UNU.RAN copies this distribution object into the generator object,
   then the address @var{extobj} that this pointer contains is simply copied.
   Thus the generator holds an address of a non-private object!
   Once the generator object has been created any change in the external 
   object might effect the generator object.
   
   @strong{Warning:} 
   External objects must be used with care. Once the generator object has been 
   created or the distribution object has been copied you @emph{must not}
   destroy this external object.
*/

const void *unur_distr_get_extobj( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the external object.

   @emph{Important:} 
   Changing this object must be done with with extreme care.
*/

/* =END */

/*---------------------------------------------------------------------------*/

UNUR_DISTR *unur_distr_clone( const UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/



