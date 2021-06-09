/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdrgw.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TDRGW                              *
 *         (Transformed Density Rejection - Gilks & Wild variant)            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  THIS METHOD HAS BEEN RENAMED TO 'ARS'.                                   *
 *  THE CORRESPONDING ROUTINES SHOULD NOT BE USED ANY MORE!                  *
 *                                                                           *
 *  Please simply replace 'tdrgw' by 'ars' to get the new function names.    *
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
   =deprecatedMETHOD  TDRGW   Transformed Density Rejection - Gilks&Wild variant

   =UP  Methods_for_CONT

   =REQUIRED concave logPDF, derivative of logPDF

   =OPTIONAL mode

   =SPEED Set-up: fast, Sampling: slow

   =REINIT supported

   =REF  [GWa92] [HLD04: Cha.4]

   =DESCRIPTION
      Same as method ARS.
      Please simply replace 'tdrgw' by 'ars' to get the new function names.

   =HOWTOUSE
      Same as method ARS.
      Please simply replace 'tdrgw' by 'ars' to get the new function names.
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

#define unur_tdrgw_new( distr ) \
   (unur_ars_new( distr ))

#define unur_tdrgw_set_max_intervals( par,max_ivs ) \
   (unur_ars_set_max_intervals( (par),(max_ivs) ))

#define unur_tdrgw_set_cpoints( par,n_cpoints,cpoints ) \
   (unur_ars_set_cpoints( (par),(n_cpoints),(cpoints) ))

#define unur_tdrgw_set_reinit_percentiles( par,n_percentiles,percentiles ) \
   (unur_ars_set_reinit_percentiles( (par),(n_percentiles),(percentiles) ))

#define unur_tdrgw_chg_reinit_percentiles( gen,n_percentiles,percentiles ) \
   (unur_ars_chg_reinit_percentiles( (gen),(n_percentiles),(percentiles) ))

#define unur_tdrgw_set_reinit_ncpoints( par,ncpoints ) \
   (unur_ars_set_reinit_ncpoints( (par),(ncpoints) ))

#define unur_tdrgw_chg_reinit_ncpoints( gen,ncpoints ) \
   (unur_ars_chg_reinit_ncpoints( (gen),(ncpoints) ))

#define unur_tdrgw_set_verify( par,verify ) \
   (unur_ars_set_verify( (par),(verify) ))

#define unur_tdrgw_chg_verify( gen,verify ) \
   (unur_ars_chg_verify( (gen),(verify) ))

#define unur_tdrgw_set_pedantic( par,pedantic ) \
   (unur_ars_set_pedantic( (par),(pedantic) ))

#define unur_tdrgw_get_loghatarea( gen ) \
   (unur_ars_get_loghatarea( gen ))

#define unur_tdrgw_eval_invcdfhat( gen,u ) \
   (unur_ars_eval_invcdfhat( (gen),(u) ))

/*---------------------------------------------------------------------------*/
