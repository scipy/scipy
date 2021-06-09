/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: error.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *      Global variables and prototypes for error handling.                  *
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
   =NODE  Error_reporting  Error reporting

   =UP Error_Debug [30] 

   =DESCRIPTION
      UNU.RAN routines report an error whenever they cannot perform the
      requested task. For example, applying transformed density
      rejection to a distribution that violates the T-concavity
      condition, or trying to set a parameter that is out of range,
      result in an error message.
      It might also happen that the setup fails for transformed density
      rejection for a T-concave distribution with some extreme density
      function simply because of round-off errors that makes the
      generation of a hat function numerically impossible.
      Situations like this may happen when using black box algorithms and
      you should check the return values of all routines.
      
      All @command{..._set_...}, and @command{..._chg_...} calls
      return @code{UNUR_SUCCESS} if they could be executed
      successfully. Otherwise, some error codes are returned if it was
      not possible to set or change the desired parameters,
      e.g. because the given values are out of range, or simply
      because the set call does not work for the chosen method.
      
      All routines that return a pointer to the requested object will
      return a NULL pointer in case of error.
      (Thus you should always check the pointer to avoid possible
      segmentation faults. Sampling routines usually do not check the
      given pointer to the generator object.)
      
      The library distinguishes between two major classes of error:
      
      @table @emph

      @item (fatal) errors:
      The library was not able to construct the
      requested object. 
      
      @item warnings:
      Some problems encounters while constructing a generator
      object. The routine has tried to solve the problem but the resulting
      object might not be what you want. For example, chosing a special
      variant of a method does not work and the initialization routine
      might switch to another variant. Then the generator produces random
      variates of the requested distribution but correlation induction
      is not possible. However, it also might happen that 
      changing the domain of a distribution has failed. Then the generator
      produced random variates with too large/too small range, i.e. their
      distribution is not correct.
      @end table

      It is obvious from the example that this distinction between errors
      and warning is rather crude and sometimes arbitrary. 
      
      UNU.RAN routines use the global variable @var{unur_errno} to
      report errors, completely analogously to @var{errno} in the ANSI
      C standard library.
      (However this approach is not thread-safe. There can 
      be only one instance of a global variable per program. Different
      threads of execution may overwrite @var{unur_errno}
      simultaneously). 
      Thus when an error occurs the caller of the routine can examine the
      error code in @var{unur_errno} to get more details about the
      reason why a routine failed. You get a short
      description of the error by a unur_get_strerror() call.
      All the error code numbers have prefix @code{UNUR_ERR_} and expand
      to non-zero constant unsigned integer values. 
      Error codes are divided into six main groups, 
      see @ref{Errno,,Error codes}.

      Alternatively, the variable @var{unur_errno} can also read by a 
      unur_get_errno() call and can be reset by the unur_reset_errno()
      call (this is in particular required for the Windows version of the
      library).

      Additionally, there exists a error handler
      (@pxref{Error_handlers,,Error handlers}) that is invoked in case
      of an error.


      In addition to reporting errors by setting error codes in
      @var{unur_errno}, the library also has an error handler
      function. This function is called by other library functions
      when they report an error, just before they return to the
      caller (@pxref{Error_handlers,,Error handlers}).
      The default behavior of the error handler is to print a short
      message:

      @smallexample
      AROU.004: [error] arou.c:1500 - (generator) condition for method violated:
      AROU.004: ..>  PDF not unimodal
      @end smallexample

      The purpose of the error handler is to provide a function
      where a breakpoint can be set that will catch library errors when
      running under the debugger. It is not intended for use in production
      programs, which should handle any errors using the return codes.

   =END
*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

extern int unur_errno;
/*
  Global variable for reporting diagnostics of error.
*/

int unur_get_errno ( void );
/* 
   Get current value of global variable @var{unur_errno}.
*/

void unur_reset_errno ( void );
/* 
   Reset global variable @var{unur_errno} to @code{UNUR_SUCCESS} 
   (i.e., no errors occured).
*/

const char *unur_get_strerror ( const int errnocode );
/*
  Get a short description for error code value.
*/

/* =END */

/* =EON */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*
   =NODE  Error_handlers  Error handlers

   =UP Error_Debug [50] 

   =DESCRIPTION
      The default behavior of the UNU.RAN error handler is to print a
      short message onto the output stream, usually a logfile
      (@pxref{Output_streams,,Output streams}), e.g.,

      @smallexample
      AROU.004: [error] arou.c:1500 - (generator) condition for method violated:
      AROU.004: ..>  PDF not unimodal
      @end smallexample

      This error handler can be switched off using the
      unur_set_error_handler_off() call, or replace it by a new one.
      Thus it allows to set a breakpoint that will catch library errors when
      running under the debugger. It also can be used to redirect
      error messages when UNU.RAN is included in general purpose
      libraries or in interactive programming environments.

      @deftp {Data Type} UNUR_ERROR_HANDLER

      This is the type of UNU.RAN error handler functions. An error
      handler will be passed six arguments which specify 
      the identifier of the object where the error occured (a string),
      the name of the source file in which it occurred (also a string), 
      the line number in that file (an integer),
      the type of error (a string: @code{"error"} or @code{"warning"}),
      the error number (an integert), and
      the reason for the error (a string). 
      The source file and line number are set at compile time
      using the @code{__FILE__} and @code{__LINE__} directives in the
      preprocessor. 
      The error number can be translated into a short description
      using a unur_get_strerror() call.
      An error handler function returns type @code{void}.
      
      Error handler functions should be defined like this,
      @example
      void my_handler( 
      @ @ @ @ @ @ @ @ @ @ @ const char *objid, 
      @ @ @ @ @ @ @ @ @ @ @ const char *file,
      @ @ @ @ @ @ @ @ @ @ @ int line,
      @ @ @ @ @ @ @ @ @ @ @ const char *errortype,
      @ @ @ @ @ @ @ @ @ @ @ int unur_errno, 
      @ @ @ @ @ @ @ @ @ @ @ const char *reason )
      @end example
      @end deftp

      To request the use of your own error handler you need the call
      unur_set_error_handler().

   =END

*/

/* =ROUTINES */


UNUR_ERROR_HANDLER *unur_set_error_handler( UNUR_ERROR_HANDLER *new_handler );
/* 
   This function sets a new error handler, @var{new_handler}, for the
   UNU.RAN library routines. The previous handler is returned (so that you
   can restore it later). Note that the pointer to a user defined
   error handler function is stored in a static variable, so there
   can be only one error handler per program. This function should
   be not be used in multi-threaded programs except to set up a
   program-wide error handler from a master thread. 

   To use the default behavior set the error handler to NULL.
*/

UNUR_ERROR_HANDLER *unur_set_error_handler_off( void );
/* 
   This function turns off the error handler by defining an error
   handler which does nothing (except of setting @var{unur_errno}. 
   The previous handler is returned (so that you can restore it later).
*/

/* =END */

/* =EON */

/*---------------------------------------------------------------------------*/

/*

      @example
      void my_handler( const char *objid, const char *file, int line, 
                                 const char *errortype, int unur_errno, const char *reason );
      void handler (const char * reason, 
              const char * file, 
              int line, 
              int gsl_errno)
      @end example
*/
