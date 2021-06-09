/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_methods.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes and macros for deprecated routines            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   THESE ROUTINES SHOULD NOT BE USED ANY MORE!                             *
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
#ifndef UNUR_DEPRECATED_METHODS_H_SEEN
#define UNUR_DEPRECATED_METHODS_H_SEEN
/*---------------------------------------------------------------------------*/

int unur_cstd_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. If the given parameters are invalid for the distribution,
   no parameters are set.
   Notice that optional parameters are (re-)set to their default values if 
   not given for UNURAN standard distributions.
*/

/*---------------------------------------------------------------------------*/

int unur_dari_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions has been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_dari_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. Notice that this call simply copies the parameters into
   the generator object. Thus if fewer parameters are provided then
   the remaining parameters are left unchanged.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.

   @emph{Important:} The given parameters are not checked against
   domain errors; in opposition to the @command{unur_<distr>_new} calls.
*/

int unur_dari_chg_domain( UNUR_GEN *generator, int left, int right );
/* 
   Change the left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_dari_chg_mode() call is required.
   (There is no domain checking as in the unur_init() call.)
   Use @code{INT_MIN} and @code{INT_MAX} for (minus) infinity.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_chg_mode( UNUR_GEN *generator, int mode );
/* 
   Change mode of distribution.
   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/


int unur_dari_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. This call only works well
   when a distribution object from the UNURAN library of standard
   distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise a (slow) numerical mode finder is called.
   If no mode can be found, then an error code is returnded and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_chg_pmfsum( UNUR_GEN *generator, double sum );
/* 
   Change sum over the PMF of distribution.
   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_upd_pmfsum( UNUR_GEN *generator );
/* 
   Recompute sum over the PMF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise an error code is returned and @code{unur_errno} 
   is set to @code{UNUR_ERR_DISTR_DATA}.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

/*---------------------------------------------------------------------------*/

int unur_dsrou_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distribution have been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_dsrou_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_dsrou_chg_domain( UNUR_GEN *generator, int left, int right );
/* 
   Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_dsrou_chg_mode() is required.
   (There is no checking whether the domain is set or not as in the
   unur_init() call.)
*/

int unur_dsrou_chg_mode( UNUR_GEN *generator, int mode );
/* 
   Change mode of distribution.
   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_chg_pmfsum( UNUR_GEN *generator, double sum );
/* 
   Change sum over PMF of distribution.
   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_upd_pmfsum( UNUR_GEN *generator );
/*
   Recompute the sum over the the PMF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

/*---------------------------------------------------------------------------*/

int unur_dstd_chg_pmfparams( UNUR_GEN *gen, double *params, int n_params );
/*
   Change array of parameters of the distribution in a given generator
   object. If the given parameters are invalid for the distribution,
   no parameters are set.
   Notice that optional parameters are (re-)set to their default values if 
   not given for UNURAN standard distributions.

   @emph{Important:} Integer parameter must be given as doubles.
*/

/*---------------------------------------------------------------------------*/

int unur_ninv_chg_pdfparams(UNUR_GEN *generator, double *params, int n_params);
/*
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/ 

/*---------------------------------------------------------------------------*/

int unur_srou_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions have been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_srou_chg_mode() is required.
   (There is no checking whether the domain is set or not as in the
   unur_init() call.)
*/

int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
/* 
   Change mode of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
/* 
   Change area below PDF of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_upd_pdfarea( UNUR_GEN *generator );
/*
   Recompute the area below the PDF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

/*---------------------------------------------------------------------------*/

int unur_ssr_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distribution has been changed (see below).
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_ssr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_ssr_chg_domain( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the distribution.  
   If the mode changes when the domain of the distribution is 
   changed, then a correspondig unur_ssr_chg_mode() is required.
   (There is no domain checking as in the unur_init() call.)
*/

int unur_ssr_chg_mode( UNUR_GEN *generator, double mode );
/* 
   Change mode of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_chg_pdfarea( UNUR_GEN *generator, double area );
/* 
   Change area below PDF of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_upd_pdfarea( UNUR_GEN *generator );
/*
   Recompute the area below the PDF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

/*---------------------------------------------------------------------------*/

int unur_tdr_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters of the
   distribution has been changed. 
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is
   returned. In case of a failure an error code is returned. Then
   @var{generator} cannot be used before another successful reinit
   (with proper parameters for the underlying distribution).
*/

/*---------------------------------------------------------------------------*/

int unur_tdrgw_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters of the
   distribution has been changed. 
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is
   returned. In case of a failure an error code is returned. Then
   @var{generator} cannot be used before another successful reinit
   (with proper parameters for the underlying distribution).
*/

/*---------------------------------------------------------------------------*/

int unur_utdr_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions has been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.

   @emph{Important:} Do not use the @var{generator} object for
   sampling after a failed reinit, since otherwise it may produce
   garbage.
*/

int unur_utdr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_utdr_chg_domain( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_utdr_chg_mode() is required.
   (There is no domain checking as in the unur_init() call.)
*/

int unur_utdr_chg_mode( UNUR_GEN *generator, double mode );
/* 
   Change mode of distribution.
   unur_utdr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_utdr_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_utdr_chg_pdfarea( UNUR_GEN *generator, double area );
/* 
   Change area below PDF of distribution.
   unur_utdr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_utdr_upd_pdfarea( UNUR_GEN *generator );
/*
   Recompute the area below the PDF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_DEPRECATED_METHODS_H_SEEN */
/*---------------------------------------------------------------------------*/
