/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr_debug.c                                                *
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
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = MVTDR (Multi-Variate Transformed Density Rejection)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_mvtdr_sample_cvec()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: PDF(center) = %g\n",gen->genid, GEN->pdfcenter);
  fprintf(LOG,"%s: bound for splitting cones = %g * mean volume\n",gen->genid,GEN->bound_splitting);
  fprintf(LOG,"%s: maximum number of cones = %d\n",gen->genid,GEN->max_cones);

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_mvtdr_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_init_finished( const struct unur_gen *gen, int successful )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   successful ... whether creation of was succesfull                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (!successful) {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
  }

  /* gamma variates */
  fprintf(LOG,"%s: upper bound for gamma variates = %g\n",gen->genid,GEN->max_gamma);

  /* triangulation steps */
  fprintf(LOG,"%s: minimum triangulation level = %d\n",gen->genid,GEN->steps_min);
  fprintf(LOG,"%s: maximum triangulation level = %d\n",gen->genid,GEN->n_steps);

  /* number of vertices and cones */
  fprintf(LOG,"%s: number of cones = %d\n",gen->genid,GEN->n_cone);
  fprintf(LOG,"%s: number of vertices = %d\n",gen->genid,GEN->n_vertex);

  /* volume below hat */
  fprintf(LOG,"%s: volume below hat = %g",gen->genid,GEN->Htot);
  if (gen->distr->set & UNUR_DISTR_SET_PDFVOLUME)
    fprintf(LOG,"\t[ hat/pdf ratio = %g ]",GEN->Htot/DISTR.volume);
  fprintf(LOG,"\n");

  if (gen->debug & MVTDR_DEBUG_VERTEX)
    _unur_mvtdr_debug_vertices(gen);

  if (gen->debug & MVTDR_DEBUG_CONE)
    _unur_mvtdr_debug_cones(gen);

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fflush(LOG);

} /* end of _unur_mvtdr_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_vertices( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print list of vertices                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  VERTEX *vt;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: List of vertices: %d\n",gen->genid,GEN->n_vertex);

  for (vt = GEN->vertex; vt != NULL; vt = vt->next) {
    fprintf(LOG,"%s: [%4d] = ( %g",gen->genid,vt->index,vt->coord[0]);
    for (i=1; i<GEN->dim; i++)
      fprintf(LOG,", %g",vt->coord[i]);
    fprintf(LOG," )\n");
  }

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_mvtdr_debug_vertices() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_cones( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print list of cones.                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  CONE *c;
  int i,n;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: List of cones: %d\n",gen->genid,GEN->n_cone);

  for (c = GEN->cone,n=0; c != NULL; c = c->next, n++) {
    fprintf(LOG,"%s: [%4d|%2d] = { %d", gen->genid, n, c->level, (c->v[0])->index);
    for (i=1; i<GEN->dim; i++)
      fprintf(LOG,", %d",(c->v[i])->index);
    fprintf(LOG," }\n");
    fprintf(LOG,"%s:\tgv     = ( %g", gen->genid, c->gv[0]);
    for (i=1; i<GEN->dim; i++)
      fprintf(LOG,", %g", c->gv[i]);
    fprintf(LOG," )\n");
    fprintf(LOG,"%s:\tHi     = %g\t[ %g%% ]\n", gen->genid, c->Hi, 100.*c->Hi/GEN->Htot);
    fprintf(LOG,"%s:\ttp     = %g\n", gen->genid, c->tp);
    fprintf(LOG,"%s:\tf(tp)  = %g\n", gen->genid, exp(c->Tfp));
    fprintf(LOG,"%s:\theight = %g\n", gen->genid, c->height);
  }

  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);

} /* end of _unur_mvtdr_debug_cones() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
