/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_errno.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines error codes.                                              *
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
#ifndef UNUR_ERRNO_H_SEEN
#define UNUR_ERRNO_H_SEEN
/*---------------------------------------------------------------------------*/

/*
   =NODE  Errno  Error codes

   =UP Error_Debug [40] 

   =DESCRIPTION

      @subsubheading List of error codes

      @itemize @bullet
      @item Procedure executed successfully (no error)
      @ftable @code
      @item UNUR_SUCCESS (0x0u)
      success (no error)
      @end ftable

      @item Errors that occurred while handling distribution objects.
      @ftable @code
      @item UNUR_ERR_DISTR_SET
      set failed (invalid parameter).
      @item UNUR_ERR_DISTR_GET
      get failed (parameter not set).
      @item UNUR_ERR_DISTR_NPARAMS
      invalid number of parameters.
      @item UNUR_ERR_DISTR_DOMAIN
      parameter(s) out of domain.
      @item UNUR_ERR_DISTR_GEN
      invalid variant for special generator.
      @item UNUR_ERR_DISTR_REQUIRED
      incomplete distribution object, entry missing.
      @item UNUR_ERR_DISTR_UNKNOWN
      unknown distribution, cannot handle.
      @item UNUR_ERR_DISTR_INVALID
      invalid distribution object.
      @item UNUR_ERR_DISTR_DATA
      data are missing.
      @item UNUR_ERR_DISTR_PROP
      desired property does not exist
      @end ftable

      @item Errors that occurred while handling parameter objects.
      @ftable @code
      @item UNUR_ERR_PAR_SET 
      set failed (invalid parameter)
      @item UNUR_ERR_PAR_VARIANT
      invalid variant -> using default
      @item UNUR_ERR_PAR_INVALID
      invalid parameter object
      @end ftable

      @item Errors that occurred while handling generator objects.
      @ftable @code
      @item UNUR_ERR_GEN
      error with generator object.
      @item UNUR_ERR_GEN_DATA
      (possibly) invalid data.
      @item UNUR_ERR_GEN_CONDITION
      condition for method violated.
      @item UNUR_ERR_GEN_INVALID
      invalid generator object.
      @item UNUR_ERR_GEN_SAMPLING
      sampling error.
      @item UNUR_ERR_NO_REINIT
      reinit routine not implemented.
      @item UNUR_ERR_NO_QUANTILE
      quantile routine not implemented.
      @end ftable

      @item Errors that occurred while handling URNG objects.
      @ftable @code
      @item UNUR_ERR_URNG
      generic error with URNG object.
      @item UNUR_ERR_URNG_MISS
      missing functionality.
      @end ftable

      @item Errors that occurred while parsing strings.
      @ftable @code
      @item UNUR_ERR_STR
      error in string.
      @item UNUR_ERR_STR_UNKNOWN
      unknown keyword.
      @item UNUR_ERR_STR_SYNTAX
      syntax error.
      @item UNUR_ERR_STR_INVALID
      invalid parameter.
      @item UNUR_ERR_FSTR_SYNTAX
      syntax error in function string.
      @item UNUR_ERR_FSTR_DERIV
      cannot derivate function.
      @end ftable

      @item Other run time errors.
      @ftable @code
      @item UNUR_ERR_DOMAIN
      argument out of domain.
      @item UNUR_ERR_ROUNDOFF
      (serious) round-off error.
      @item UNUR_ERR_MALLOC
      virtual memory exhausted.
      @item UNUR_ERR_NULL
      invalid NULL pointer.
      @item UNUR_ERR_COOKIE
      invalid cookie.
      @item UNUR_ERR_GENERIC
      generic error.
      @item UNUR_ERR_SILENT
      silent error (no error message).
      @item UNUR_ERR_INF
      infinity occured.
      @item UNUR_ERR_NAN
      NaN occured.
      @item UNUR_ERR_COMPILE
      Requested routine requires different compilation switches.
      Recompilation of library necessary.
      @item UNUR_ERR_SHOULD_NOT_HAPPEN
      Internal error, that should not happen.
      Please report this bug!
      @end ftable

      @end itemize

   =END

   =EON
*/

/*---------------------------------------------------------------------------*/

enum { 

  /** procedure executed successfully **/
  UNUR_SUCCESS            = 0x00,     /* exited successfully                 */                          

  /** procedure executed with error (for internal use)  **/
  UNUR_FAILURE            = 0x01,     /* failure                             */

  /** distribution object **/
  /*
    @code{UNUR_ERR_DISTR_...}
    Errors that occurred while handling distribution objects.
  */
  UNUR_ERR_DISTR_SET      = 0x11,     /* set failed (invalid parameter)      */
  UNUR_ERR_DISTR_GET      = 0x12,     /* get failed (parameter not set)      */
  UNUR_ERR_DISTR_NPARAMS  = 0x13,     /* invalid number of parameters        */
  UNUR_ERR_DISTR_DOMAIN   = 0x14,     /* parameter out of domain             */
  UNUR_ERR_DISTR_GEN      = 0x15,     /* invalid variant for special generator */
  UNUR_ERR_DISTR_REQUIRED = 0x16,     /* incomplete distribution object, entry missing */
  UNUR_ERR_DISTR_UNKNOWN  = 0x17,     /* unknown distribution, cannot handle */
  UNUR_ERR_DISTR_INVALID  = 0x18,     /* invalid distribution object         */
  UNUR_ERR_DISTR_DATA     = 0x19,     /* data are missing                    */
  UNUR_ERR_DISTR_PROP     = 0x20,     /* desired property does not exist     */

  /** parameter object **/
  /*
    @code{UNUR_ERR_PAR_...}
    Errors that occurred while handling parameter objects.
  */
  UNUR_ERR_PAR_SET        = 0x21,     /* set failed (invalid parameter)      */
  UNUR_ERR_PAR_VARIANT    = 0x22,     /* invalid variant -> using default    */
  UNUR_ERR_PAR_INVALID    = 0x23,     /* invalid parameter object            */

  /** generator object **/
  /*
    @code{UNUR_ERR_GEN_...}
    Errors that occurred while handling generator objects.
  */
  UNUR_ERR_GEN            = 0x31,     /* bit for generator object            */
  UNUR_ERR_GEN_DATA       = 0x32,     /* (possible) invalid data             */
  UNUR_ERR_GEN_CONDITION  = 0x33,     /* condition for method violated       */
  UNUR_ERR_GEN_INVALID    = 0x34,     /* invalid generator object            */
  UNUR_ERR_GEN_SAMPLING   = 0x35,     /* sampling error                      */
  UNUR_ERR_NO_REINIT      = 0x36,     /* reinit not implemented              */
  UNUR_ERR_NO_QUANTILE    = 0x37,     /* qunantile not implemented           */

  /** uniform random number generator (URNG) object **/
  /*
    @code{UNUR_ERR_URNG_...}
    Errors that occurred while handling URNG objects.
  */

  UNUR_ERR_URNG           = 0x41,     /* generic error with URNG object      */
  UNUR_ERR_URNG_MISS      = 0x42,     /* missing functionality               */

  /** string parser **/
  /*
    @code{UNUR_ERR_STR_...}
    Errors that occurred while parsing strings.
  */
  UNUR_ERR_STR            = 0x51,     /* error in stringparser               */
  UNUR_ERR_STR_UNKNOWN    = 0x52,     /* unknown key word in string          */
  UNUR_ERR_STR_SYNTAX     = 0x53,     /* syntax error in string              */
  UNUR_ERR_STR_INVALID    = 0x54,     /* invalid parameter in argument       */
  UNUR_ERR_FSTR_SYNTAX    = 0x55,     /* syntax error in function parser     */
  UNUR_ERR_FSTR_DERIV     = 0x56,     /* cannot derivate function            */

  /** misc **/
  /*
    @code{UNUR_ERR_...}
    Other errors.
  */
  UNUR_ERR_DOMAIN         = 0x61,     /* argument out of domain              */
  UNUR_ERR_ROUNDOFF       = 0x62,     /* (serious) round-off error           */
  UNUR_ERR_MALLOC         = 0x63,     /* virtual memory exhausted            */
  UNUR_ERR_NULL           = 0x64,     /* invalid NULL pointer                */ 
  UNUR_ERR_COOKIE         = 0x65,     /* invalid cookie                      */
  UNUR_ERR_GENERIC        = 0x66,     /* generic error                       */
  UNUR_ERR_SILENT         = 0x67,     /* silent error (no error message)     */
  UNUR_ERR_INF            = 0x68,     /* infinity occured                    */
  UNUR_ERR_NAN            = 0x69,     /* NaN occured                         */

  /** compilation switches **/
  /*
    @code{UNUR_ERR_COMPILE}
    Requested routine requires different compilation switches.
    Recompilation of library necessary.
  */
  UNUR_ERR_COMPILE        = 0xa0,     /* not available, recompile library    */

  /** this should not happen **/
  /*
    @code{UNUR_ERR_SHOULD_NOT_HAPPEN}
    Internal error. This should not happen. 
    Please make a bug report.
  */
  UNUR_ERR_SHOULD_NOT_HAPPEN = 0xf0  /* error should not happen, report this! */

};

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_ERRNO_H_SEEN */
/*---------------------------------------------------------------------------*/


