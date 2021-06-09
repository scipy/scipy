/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_info.c                                                  *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
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
_unur_tabl_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   mode      = %g   %s\n", unur_distr_cont_get_mode(distr),
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = ");
  if (gen->distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"%g\n", DISTR.area);
  else 
    _unur_string_append(info,"[not set: use 1.0]\n");
  _unur_string_append(info,"\n");
  
  /*   if (help) { */
  /*   _unur_string_append(info,"\n"); */
  /*   } */
      
  /* method */
  _unur_string_append(info,"method: TABL (Ahrens' TABLe Method)\n");
  _unur_string_append(info,"   variant   = ");
  if (gen->variant & TABL_VARIANT_IA) 
    _unur_string_append(info,"immediate acceptance [ia = on]\n");
  else
    _unur_string_append(info,"acceptance/rejection [ia = off]\n");
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
    if (gen->variant & TABL_VARIANT_IA) 
      _unur_string_append(info,"   variant_ia = on  [default]\n");
    else
      _unur_string_append(info,"   variant_ia = off\n");

    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & TABL_SET_MAX_SQHRATIO) ? "" : "[default]");

    _unur_string_append(info,"   max_intervals = %d  %s\n", GEN->max_ivs_info,
			(gen->set & TABL_SET_MAX_IVS) ? "" : "[default]");

    if (gen->variant & TABL_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    if (gen->variant & TABL_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");

    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_tabl_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
       int unur_tabl_set_nstp( UNUR_PAR *parameters, int n_stp );
       int unur_tabl_set_useear( UNUR_PAR *parameters, int useear );
       int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
       int unur_tabl_set_usedars( UNUR_PAR *parameters, int usedars );
       int unur_tabl_set_darsfactor( UNUR_PAR *parameters, double factor );
       int unur_tabl_set_variant_splitmode( UNUR_PAR *parameters, unsigned splitmode );
       int unur_tabl_set_slopes( UNUR_PAR *parameters, const double *slopes, int n_slopes );
       int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
       int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
    */
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & TABL_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio)
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_intervals\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tabl_info() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
