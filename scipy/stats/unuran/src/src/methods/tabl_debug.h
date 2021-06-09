/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_debug.c                                                 *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(LOG,"%s:\n",gen->genid);

/*---------------------------------------------------------------------------*/

void
_unur_tabl_debug_init_start( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into LOG file                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  empty_line();
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  if (gen->variant & TABL_VARIANT_IA) 
    fprintf(LOG,"immediate acceptance ... IA\n");
  else
    fprintf(LOG,"acceptance/rejection ... RH\n");
  empty_line();

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_tabl_",gen->genid);
  if (gen->variant & TABL_VARIANT_IA) 
    fprintf(LOG,"ia"); 
  else
    fprintf(LOG,"rh");
  if (gen->variant & TABL_VARFLAG_VERIFY)
    fprintf(LOG,"_sample_check()\n");
  else
    fprintf(LOG,"_sample()\n");

  fprintf(LOG,"%s: computation interval = (%g, %g)\n",gen->genid,GEN->bleft,GEN->bright);

  fprintf(LOG,"%s: maximum number of intervals = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,TABL_SET_MAX_IVS);
  fprintf(LOG,"\n");

  if (gen->variant & TABL_VARFLAG_USEEAR) {
    fprintf(LOG,"%s: use equal area rule\n",gen->genid);
    fprintf(LOG,"%s:\tarea fraction for equal area rule = %g ",gen->genid,PAR->area_fract);
    _unur_print_if_default(par,TABL_SET_AREAFRACTION);
    fprintf(LOG,"\n");
  }

  if (gen->variant & TABL_VARFLAG_USEDARS) {
    fprintf(LOG,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(gen,TABL_SET_USE_DARS);
    fprintf(LOG,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(gen,TABL_SET_DARS_FACTOR);
  }
  else {
    fprintf(LOG,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(gen,TABL_SET_USE_DARS);
  }
  fprintf(LOG,"\n");

  fprintf(LOG,"%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,GEN->max_ratio*100.);
  _unur_print_if_default(gen,TABL_SET_MAX_SQHRATIO);
  fprintf(LOG,"\n");

  fprintf(LOG,"%s: split intervals at ",gen->genid);
  switch( gen->variant & TABL_VARMASK_SPLIT ) {
  case TABL_VARFLAG_SPLIT_MEAN:
    fprintf(LOG,"mean point");
    break;
  case TABL_VARFLAG_SPLIT_ARC:
    fprintf(LOG,"\"arcmean\" point");
    break;
  case TABL_VARFLAG_SPLIT_POINT:
  default:
    fprintf(LOG,"sample point");
    break;
  }
  fprintf(LOG," when using adaptive sampling.\n");
  empty_line();

  fprintf(LOG,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(LOG,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,TABL_SET_GUIDEFACTOR);
  fprintf(LOG,"\n");
  empty_line();

  if (par->set & TABL_SET_SLOPES) {
    fprintf(LOG,"%s: slopes = %d\n",gen->genid,PAR->n_slopes);
    for (i=0; i<PAR->n_slopes; i++) {
      if ( PAR->slopes[2*i] > PAR->slopes[2*i+1] )
	fprintf(LOG,"%s:   (+)  ",gen->genid);
      else
	fprintf(LOG,"%s:   (-)  ",gen->genid);
      fprintf(LOG,"< %#g, %#g >\n", PAR->slopes[2*i], PAR->slopes[2*i+1] );
    }
  }
  else
    fprintf(LOG,"%s: no slopes given. compute from domain and mode.\n",gen->genid);

  if (PAR->cpoints != NULL) {
    fprintf(LOG,"%s: cpoints for slopes (%d):\n",gen->genid,PAR->n_cpoints);
    fprintf(LOG,"%s:\t %g",gen->genid,(PAR->cpoints)[0]);
    for (i=1; i<PAR->n_cpoints; i++)
      fprintf(LOG,", %g",(PAR->cpoints)[i]);
    fprintf(LOG,"\n");
  }

  fprintf(LOG,"%s: number of starting intervals (approx.) = %d",gen->genid,PAR->n_stp);
  _unur_print_if_default(par,TABL_SET_N_STP);
  fprintf(LOG,"\n");
  empty_line();

} /* end of _unur_tabl_debug_init_start() */

/*****************************************************************************/

void
_unur_tabl_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into LOG file                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  _unur_tabl_debug_intervals(gen,"INIT completed",TRUE);

  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  empty_line();

} /* end of _unur_tabl_debug_init_finished() */

/*****************************************************************************/

void 
_unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print header before runniung DARS into LOG file                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: DARS started **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fflush(LOG);
} /* end of _unur_tabl_debug_dars_start() */

/*****************************************************************************/

void
_unur_tabl_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into LOG file           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  empty_line();
  if (gen->status == UNUR_SUCCESS) {
    fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
    empty_line();
    _unur_tabl_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
  }
  empty_line();

  fflush(LOG);

} /* end of _unur_tabl_debug_free() */

/*****************************************************************************/

void
_unur_tabl_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into LOG file                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen         ... pointer to generator object                        */
     /*   header      ... header for table                                   */
     /*   print_areas ... whether table of areas should be printed           */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  struct unur_tabl_interval *iv;
  double sAsqueeze, Atotal;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (header)
    fprintf(LOG,"%s: %s\n",gen->genid,header);

  fprintf(LOG,"%s: intervals = %d\n",gen->genid,GEN->n_ivs);
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(LOG,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID);
      fprintf(LOG,"%s:[%3d]: (%s)   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      gen->genid, i, (iv->xmin <= iv->xmax)?"+":"-",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }

  if (!print_areas) return;

  if (GEN->Atotal <= 0.) {
    fprintf(LOG,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(LOG,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    Atotal = -1.;   /* to avoid floating point exceptions */
  }
  else {
    Atotal = GEN->Atotal;
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(LOG,"%s:Areas in intervals:\n",gen->genid);
    fprintf(LOG,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID); 
      sAsqueeze += iv->Asqueeze;
      fprintf(LOG,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
	      iv->Ahat, iv->Ahat * 100. / Atotal, 
	      iv->Acum, iv->Acum * 100. / Atotal);
    }
    fprintf(LOG,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(LOG,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / Atotal, Atotal);
    empty_line();
  }
    
  /* summary of areas */
  fprintf(LOG,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(LOG,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(LOG,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);

  empty_line();

} /* end of _unur_tabl_debug_intervals */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
