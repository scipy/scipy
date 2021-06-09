/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr_info.c                                                 *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    multivariate transformed density rejection                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given (logarithm of the) PDF of a log-concave distribution;          *
 *      produce a value x consistent with its density.                       *
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

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

/*****************************************************************************/
/**  Info string                                                            **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_info( struct unur_gen *gen, int help )
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
  int samplesize = 10000;
  double rc;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_distr_cvec_info_domain(gen);

  if ( distr->set & UNUR_DISTR_SET_MODE ) {
    _unur_string_append(info,"   mode      = ");
    _unur_distr_info_vector( gen, DISTR.mode, GEN->dim);
  }
  _unur_string_append(info,"\n");

  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, GEN->center, GEN->dim);
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]");
    else
      _unur_string_append(info,"  [default]");
  }
  _unur_string_append(info,"\n\n");
  
  if (help) {
    if ( !(distr->set & UNUR_DISTR_SET_MODE) ) 
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can set the mode to improve the rejection constant.");
    _unur_string_append(info,"\n");
  } 

  /* method */
  _unur_string_append(info,"method: MVTDR (Multi-Variate Transformed Density Rejection)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   volume(hat) = %g\n", GEN->Htot);

  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFVOLUME) {
    _unur_string_append(info,"= %g\n", GEN->Htot / DISTR.volume);
  }
  else {
    rc = unur_test_count_urn(gen,samplesize,0,NULL)/((1.+GEN->dim)*samplesize);
    _unur_string_append(info,"= %.2f  [approx.]\n", rc);
  }

  _unur_string_append(info,"   # cones = %d\n", GEN->n_cone);
  _unur_string_append(info,"   # vertices = %d\n", GEN->n_vertex);
  if (GEN->steps_min == GEN->n_steps)
    _unur_string_append(info,"   triangulation levels = %d\n", GEN->n_steps);
  else
    _unur_string_append(info,"   triangulation levels = %d-%d\n", GEN->steps_min, GEN->n_steps);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");

    _unur_string_append(info,"   stepsmin = %d  %s\n", GEN->steps_min,
 			(gen->set & MVTDR_SET_STEPSMIN) ? "" : "[default]");

    _unur_string_append(info,"   maxcones = %d  %s\n", GEN->max_cones,
 			(gen->set & MVTDR_SET_MAXCONES) ? "" : "[default]");

    _unur_string_append(info,"   boundsplitting = %g  %s\n", GEN->bound_splitting,
 			(gen->set & MVTDR_SET_BOUNDSPLITTING) ? "" : "[default]");

    if (gen->variant & MVTDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & MVTDR_SET_STEPSMIN) )
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can increase \"stepsmin\" to improve the rejection constant." );
    if (GEN->max_cones <= GEN->n_cone)
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can increase \"maxcones\" to improve the rejection constant." );
    if ( !(gen->set & MVTDR_SET_BOUNDSPLITTING) )
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can change \"boundsplitting\" to change the creating of the hat function." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_mvtdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
