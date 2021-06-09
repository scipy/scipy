/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_info.ch                                                  *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Routines for creating info strings.                                  *
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
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_tdr_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);

  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }
  
  if (help) {
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
      
  /* method */
  _unur_string_append(info,"method: TDR (Transformed Density Rejection)\n");
  _unur_string_append(info,"   variant   = ");
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    _unur_string_append(info,"GW (original Gilks & Wild)\n"); break;
  case TDR_VARIANT_PS:
    _unur_string_append(info,"PS (proportional squeeze)\n"); break;
  case TDR_VARIANT_IA:
    _unur_string_append(info,"IA (immediate acceptance)\n"); break;
  }
  /* used transformation */
  _unur_string_append(info,"   T_c(x)    = ");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_string_append(info,"log(x)  ... c = 0\n"); break;
  case TDR_VAR_T_SQRT:
    _unur_string_append(info,"-1/sqrt(x)  ... c = -1/2\n"); break;
  case TDR_VAR_T_POW:
    _unur_string_append(info,"-x^(%g)  ... c = %g\n",GEN->c_T,GEN->c_T); break;
  }
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g\n", GEN->Atotal);

  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal/DISTR.area);
  else
    _unur_string_append(info,"<= %g\n", GEN->Atotal/GEN->Asqueeze);

  _unur_string_append(info,"   area ratio squeeze/hat = %g\n",
		      GEN->Asqueeze/GEN->Atotal);

  _unur_string_append(info,"   # intervals = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  
  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    switch (gen->variant & TDR_VARMASK_VARIANT) {
    case TDR_VARIANT_GW:
      _unur_string_append(info,"   variant_gw = on\n"); break;
    case TDR_VARIANT_PS:
      _unur_string_append(info,"   variant_ps = on  [default]\n"); break;
    case TDR_VARIANT_IA:
      _unur_string_append(info,"   variant_ia = on\n"); break;
    }

    _unur_string_append(info,"   c = %g  %s\n", GEN->c_T,
			(gen->set & TDR_SET_C) ? "" : "[default]");

    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & TDR_SET_MAX_SQHRATIO) ? "" : "[default]");

    _unur_string_append(info,"   max_intervals = %d  %s\n", GEN->max_ivs_info,
			(gen->set & TDR_SET_MAX_IVS) ? "" : "[default]");

    if (gen->variant & TDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    if (gen->variant & TDR_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");

    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_tdr_set_usedars( UNUR_PAR *parameters, int usedars );
       int unur_tdr_set_darsfactor( UNUR_PAR *parameters, double factor );
       int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
       int unur_tdr_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
       int unur_tdr_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
       int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
       int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
       int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
    */
  }


  /* Hints */
  if (help) {
    if ( (gen->variant & TDR_VARMASK_VARIANT) != TDR_VARIANT_IA) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may use \"variant_ia\" for faster generation times."); 
    if ( !(gen->set & TDR_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_intervals\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tdr_info() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
