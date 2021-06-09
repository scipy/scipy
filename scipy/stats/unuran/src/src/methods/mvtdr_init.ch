/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr_init.c                                                 *
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
 *   Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold             *
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
/**  Initialzation: Create Hat                                              **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mvtdr_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_MVTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mvtdr_create(par);

  /* free parameters */
  _unur_par_free(par);

  /* check generator object */ 
  if (!gen) return NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_mvtdr_debug_init_start(gen);
#endif

  /* check data */
  if (!(GEN->pdfcenter > 0.)) {
    _unur_error(gen->genid,UNUR_ERR_DISTR_DOMAIN,"center out of support of PDF");
    _unur_mvtdr_free(gen); return NULL;
  }

  /* make hat function */
  if(_unur_mvtdr_create_hat(gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create hat");
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_mvtdr_debug_init_finished(gen, FALSE);
#endif
    _unur_mvtdr_free(gen); return NULL;
  }

  /* compute upper bound for gamma variates */
  _unur_mvtdr_max_gamma(gen);

  /* we need an auxiliary generator for gamma random variates */
  GEN_GAMMA = _unur_mvtdr_gammagen( gen, (double)(GEN->dim) );
  if ( GEN_GAMMA == NULL ) {
      _unur_mvtdr_free(gen); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_mvtdr_debug_init_finished(gen, TRUE);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_mvtdr_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_mvtdr_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  
  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mvtdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MVTDR_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim; 

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mvtdr_sample_cvec;
  gen->destroy = _unur_mvtdr_free;
  gen->clone = _unur_mvtdr_clone;

  /* initialize counter and check given parameters */
  GEN->n_steps = 0;   /* no triangulation steps yet */
  GEN->steps_min = _unur_max( 0, PAR->steps_min );  /* minimum number of triangulation steps */
  /* check maximal number of cones */
  if ( (1 << (GEN->dim + GEN->steps_min)) > PAR->max_cones) {
    /*     WARNING( "number of cones raised to 2^(dim + T_STEPS_MIN)" ); */
    PAR->max_cones = 1 << (GEN->dim + GEN->steps_min);
  }
  GEN->max_gamma = INFINITY;            /* upper bound for gamma variaties */

  /* initialize  pointers to lists */
  GEN->cone = NULL;
  GEN->last_cone = NULL;
  GEN->n_cone = 0;                      /* number cones */
  GEN->max_cones = PAR->max_cones;      /* maximum number of cones */
  GEN->bound_splitting = PAR->bound_splitting;    /* bound for splitting cones */

  GEN->vertex = NULL;
  GEN->last_vertex = NULL;
  GEN->n_vertex = 0;                    /* maximum number of vertices */

  GEN->etable = NULL;                   /* pointer to edge table */
  GEN->etable_size = 0;                 /* size of edge table */

  GEN->guide = NULL;
  GEN->guide_size = 0;
  
  /* initialize working arrays: */
  /*   point on simples */
  GEN->S         = malloc( GEN->dim * sizeof(double) );
  /*   vector g (direction of sweeping plane) */
  GEN->g         = malloc( GEN->dim * sizeof(double) );
  /*   coordinates of touching point of hat */
  GEN->tp_coord  = malloc( GEN->dim * sizeof(double) );
  /*   coordinates of touching point of hat moved into center */
  GEN->tp_mcoord = malloc( GEN->dim * sizeof(double) );
  /*   gradient of transformed density at tp */
  GEN->tp_Tgrad  = malloc( GEN->dim * sizeof(double) );

  if (GEN->S==NULL || GEN->g==NULL || GEN->tp_coord==NULL || 
      GEN->tp_mcoord==NULL || GEN->tp_Tgrad==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,"");
    _unur_mvtdr_free(gen); return NULL;
  }

  /* get center of the distribution and its PDF */
  GEN->center = unur_distr_cvec_get_center(gen->distr);
  GEN->pdfcenter = PDF(GEN->center);

  /* whether we have set a domain for the distribution */
  GEN->has_domain = (gen->distr->set & UNUR_DISTR_SET_DOMAIN) ? TRUE : FALSE;
 
#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_mvtdr_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_mvtdr_create() */

/*****************************************************************************/

struct unur_gen *
_unur_mvtdr_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE  ((struct unur_mvtdr_gen*)clone->datap)

  struct unur_gen *clone;
  int error = FALSE;
  size_t size;
  VERTEX *vt, **vtindex;
  CONE *c;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MVTDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy data */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  /* working arrays */
  size = GEN->dim * sizeof(double);
  CLONE->S         = malloc(size);
  CLONE->g         = malloc(size);
  CLONE->tp_coord  = malloc(size);
  CLONE->tp_mcoord = malloc(size);
  CLONE->tp_Tgrad  = malloc(size);
  vtindex = malloc(GEN->n_vertex * sizeof (VERTEX *));

  if (CLONE->S==NULL || CLONE->g==NULL || CLONE->tp_coord==NULL || 
      CLONE->tp_mcoord==NULL || CLONE->tp_Tgrad==NULL || vtindex==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,"");
    if (vtindex) free (vtindex);
    _unur_mvtdr_free(clone); return NULL;
  }

  /* copy data */
  if (GEN->S) memcpy( CLONE->S, GEN->S, size );
  if (GEN->g) memcpy( CLONE->g, GEN->g, size );
  if (GEN->tp_coord) memcpy( CLONE->tp_coord, GEN->tp_coord, size );
  if (GEN->tp_mcoord) memcpy( CLONE->tp_mcoord, GEN->tp_mcoord, size );
  if (GEN->tp_Tgrad) memcpy( CLONE->tp_Tgrad, GEN->tp_Tgrad, size );

  /* clear lists in clone */
  CLONE->vertex = NULL;  CLONE->n_vertex = 0;
  CLONE->cone = NULL;    CLONE->n_cone = 0;
  CLONE->guide = NULL;

  /* copy list of vertices */
  for (vt = GEN->vertex; vt != NULL; vt = vt->next) {
    VERTEX *vtc = _unur_mvtdr_vertex_new( clone );
    if (vtc == NULL) {
      error = TRUE; break; }
    memcpy(vtc->coord, vt->coord, size);
    vtc->index = vt->index;
    vtindex[vt->index] = vtc;
  }

  /* copy list of cones */
  for (c = GEN->cone; c != NULL && !error; c = c->next) {
    CONE *cc, *cc_next;
    VERTEX **v;
    double *center, *gv;
    int i;
    cc = _unur_mvtdr_cone_new( clone );
    if (cc == NULL) {
      error = TRUE; break; }
    cc_next = cc->next;
    center = cc->center;
    gv = cc->gv;
    v = cc->v;
    memcpy(cc,c,sizeof(CONE));
    memcpy(center, c->center, size);
    memcpy(gv, c->gv, size);
    for (i=0; i<GEN->dim; i++)
      v[i] = vtindex[(c->v[i])->index];
    cc->next = cc_next;
    cc->center = center;
    cc->gv = gv;
    cc->v = v;
  }

  /* make new guide table */
  if (_unur_mvtdr_make_guide_table(clone) != UNUR_SUCCESS)
    error = TRUE;

  /* clear auxiliary array */
  free (vtindex);

  if (error == TRUE) {
    _unur_mvtdr_free(clone); return NULL;
  }

  return clone;

#undef CLONE
} /* end of _unur_mvtdr_clone() */

/*****************************************************************************/

void
_unur_mvtdr_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  VERTEX *vt, *vt_next;
  CONE *c, *c_next;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_MVTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* clear lists: */

  /* hash table for edges */
  _unur_mvtdr_etable_free(gen);

  /* linked list of vertices */
  for (vt = GEN->vertex; vt != NULL; vt = vt_next) {
    vt_next = vt->next;
    free (vt->coord);    /* coordinates of vertex */
    free (vt);
  }

  /* linked list of cones */
  for (c = GEN->cone; c != NULL; c = c_next) {
      c_next = c->next;
      free (c->v);        /* list of vertices of the cone */
      free (c->center);   /* barycenter of cone */
      free (c->gv);       /* <g,v> for all vertices v */
      free (c);
  }

  /* guide table */
  if (GEN->guide) free (GEN->guide);

  /* working arrays */
  if (GEN->S)         free (GEN->S);
  if (GEN->g)         free (GEN->g);
  if (GEN->tp_coord)  free (GEN->tp_coord);
  if (GEN->tp_mcoord) free (GEN->tp_mcoord);
  if (GEN->tp_Tgrad)  free (GEN->tp_Tgrad);

  _unur_generic_free(gen);

} /* end of _unur_mvtdr_free() */

/*****************************************************************************/

struct unur_gen *
_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha )
     /*----------------------------------------------------------------------*/
     /* create a gamma random variate generator with shape parameter alpha.  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to MVTDR generator object                        */
     /*   alpha ... shape parameter                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *gammadistr;
  struct unur_par   *gammapar;
  struct unur_gen   *gammagen;
  double shape;

  /* make generator object */
  shape = alpha;
  gammadistr = unur_distr_gamma(&shape,1);
  if (_unur_isfinite(GEN->max_gamma)) {
    unur_distr_cont_set_domain(gammadistr,0.,GEN->max_gamma);
  }
  gammapar = unur_tdr_new( gammadistr );
  unur_tdr_set_usedars( gammapar, TRUE );
  unur_tdr_set_max_sqhratio( gammapar, MVTDR_TDR_SQH_RATIO);
  if (! GEN->has_domain) {
    /* we do not need a truncated gamma distribution. */
    /* hence we can use immediate acceptance.          */
    unur_tdr_set_variant_ia( gammapar );
  }
  gammagen = unur_init( gammapar );
  _unur_distr_free( gammadistr );

  /* check result */
  if (gammagen == NULL) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
		"Cannot create aux Gamma generator");
    return NULL;
  }
  
  /* uniform random number generator and debugging flags */
  gammagen->urng = gen->urng;
  gammagen->debug = gen->debug;

  return gammagen;

} /* end of _unur_mvtdr_gammagen() */


/*****************************************************************************/
/*                                                                           */
/*   Hat.                                                                    */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_create_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create hat function.                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int step;            /* triangulation steps */
  double Hi_bound;     /* lower bound on Hi for splitting cone */
  CONE *c;
  int n_splitted;

  /* vertices of initial cones */
  if( _unur_mvtdr_initial_vertices(gen) != UNUR_SUCCESS ) 
    return UNUR_FAILURE;

  /* initial cones */
  if( _unur_mvtdr_initial_cones(gen) != UNUR_SUCCESS ) 
    return UNUR_FAILURE;

  /* execute minimal number of triangulation steps */
  for( step = 1; step <= GEN->steps_min; step++ ) {
    if (_unur_mvtdr_triangulate(gen,step,TRUE) < 0)
      return UNUR_FAILURE;
  }

  /* compute optimal distance of touching points and volume Hi */
  for( c = GEN->cone; c != NULL; c = c->next )
    _unur_mvtdr_tp_find (gen,c);

  /* cones with invalid hats (or too large volumes) must be split */
  while( _unur_mvtdr_triangulate(gen,step,FALSE) > 0 ) {
    if (GEN->n_cone > GEN->max_cones)
      return UNUR_FAILURE;
    step++;
  }

  /* maximum number of triangulations yet */
  GEN->n_steps = step-1;

  /* compute cumulated volumes in all cones */
  GEN->Htot = 0.;                 /* accumulated sum of volumes */
  for( c=GEN->cone; c!=NULL; c=c->next ) {
    /* volume below hat */
    GEN->Htot += c->Hi;           /* volume below hat */
    c->Hsum = GEN->Htot;          /* accumulated sum of volumes */
  }

  /* split until stopping criterion in reached */
  while (1) {

    /* bound for splitting cones */
    /* do until all cones have approx same hat volumes */
    Hi_bound = GEN->bound_splitting * GEN->Htot / GEN->n_cone;

    /* and now check all the cones again */
    GEN->Htot = 0.;
    n_splitted = 0;
    for( c=GEN->cone; c!=NULL; c=c->next ) {   /* all cones */
      while( Hi_bound < c->Hi && GEN->n_cone < GEN->max_cones ) { 
	/* we (must) split the cone again */
	if (_unur_mvtdr_cone_split(gen,c,c->level+1) != UNUR_SUCCESS)
	  return UNUR_FAILURE;
	++n_splitted;
	/* and compute optimal touching point */
	_unur_mvtdr_tp_find (gen,c);
	_unur_mvtdr_tp_find (gen,GEN->last_cone);
      }
      GEN->Htot += c->Hi;           /* volume below hat */
      c->Hsum = GEN->Htot;          /* accumulated sum of volumes */
      if( c == GEN->last_cone ) break;
    }

    /* stop when maximal number of cones is reached or no cones are split */
    if (!n_splitted || GEN->n_cone >= GEN->max_cones) break;
  }

  /* create guide table for finding cones */
  if (_unur_mvtdr_make_guide_table(gen) != UNUR_SUCCESS) 
    return UNUR_FAILURE;

  /* we do not need the hash table generated in triangulate() cone any more */
  if (GEN->dim > 2)
    _unur_mvtdr_etable_free(gen);

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_create_hat() */


/*****************************************************************************/
/*                                                                           */
/*   CONES.                                                                  */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_initial_cones( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get initial cones                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i,j,k;
  CONE *c;
  int max_c;       /* maximum number of cones */
  VERTEX *vt;
  VERTEX **ivtl;   /* list of initial vertices */
  int dim = GEN->dim;

  int error = FALSE;
  int have_negative_index = FALSE;
  int cone_out_of_domain;

  /* make array of initial vertices */
  ivtl = malloc(2 * dim * sizeof(VERTEX*));
  if( ivtl==NULL ) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; 
  }
  for (vt = GEN->vertex, i=0; i < 2*GEN->dim && vt!=NULL; vt = vt->next, i++) {
    if (vt->index<0) have_negative_index = TRUE;
    ivtl[i] = vt;
  }

  /* we have (at most) 2^dim initial cones */
  max_c = 1 << dim;

  /* we need vertices, index, volume for each cone */
  for( k=0; k<max_c; k++ ) {

    if (have_negative_index) {
      /* first check vertices of cone for negative index */
      cone_out_of_domain = FALSE;
      for( i=0; i < dim; i++ ) {
	if ( (!((k>>i)&1) && ivtl[i]->index < 0) ||
	     ( ((k>>i)&1) && ivtl[i+dim]->index < 0) ) {
	  cone_out_of_domain = TRUE;  break;
	}
      }
      if (cone_out_of_domain)  /* cone is not inside domain --> skip */
	continue;
    }

    /* get new (empty) cone object */
    c = _unur_mvtdr_cone_new(gen);
    if (c==NULL) { error = TRUE; break; }

    /* this is level 0 of triangulation */
    c->level = 0;

    /* each cone is incident to 'dim' edges.                           */
    /* The i-th edge of the cone is either GEN->v[i] or GEN->v[dim+i], */
    /* (the latter is equal to (-1)*(GEN->v[i])).                      */
    /* The indices of the vertices must be in ascending order.         */ 
    j = 0;
    for( i=0; i < dim; i++ )
      if (!((k>>i)&1))  (c->v)[j++] = ivtl[i];
    for( i=0; i < dim && j < dim; i++ )
      if ( ((k>>i)&1))  (c->v)[j++] = ivtl[i + dim];

    /* determinant and volume of triangle * ((dim-1)!) */
    c->logdetf = 0.;

    /* touching point not known yet */
    c->tp = -1.;   /* > 0. if and only if tp is computed !! */
  }

  /* free list of initial vertices */
  free (ivtl);

  if (error==TRUE) return UNUR_ERR_MALLOC;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_initial_cones() */

/*---------------------------------------------------------------------------*/

CONE *
_unur_mvtdr_cone_new( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* allocate new cone and append it to linked list of all cones.         */
     /* increment counter for cones.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to newly allocated cone                                    */
     /*----------------------------------------------------------------------*/
{
  CONE *c; 
  
  /* allocate memory */
  c = malloc(sizeof(CONE));
  if (c==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }

  /* insert into list of cones */
  if (GEN->cone == NULL)
    GEN->last_cone = GEN->cone = c;
  else
    GEN->last_cone = GEN->last_cone->next = c;
  c->next = NULL;

  /* list of vertices of the cone */
  c->v      = malloc( GEN->dim * sizeof(VERTEX *));

  /* barycenter of cone */
  c->center = malloc( GEN->dim * sizeof(double));

  /* <g,v> for all vertices v */
  c->gv     = malloc( GEN->dim * sizeof(double));

  if (c->v==NULL || c->center==NULL || c->gv==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }

  /* the cone is unbounded when created */
  c->height = UNUR_INFINITY;

  /* mark as invalid */
  c->tp = -1.;
  c->Hi = INFINITY;

  /* and update counter */
  ++(GEN->n_cone);

  /* return pointer to next vertex */
  return c;

} /* end of _unur_mvtdr_cone_new() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c )
     /*----------------------------------------------------------------------*/
     /* computer center of cone and normalize the corresponding vector       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   c   ... cone for which center has to be computed                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i,k;
  double norm;
  int dim = GEN->dim;
  
  /* compute sum of all vertices and square of its norm */
  norm = 0.;
  for( i=0; i<dim; i++ ) {
    c->center[i] = 0.;
    for( k=0; k<dim; k++ )
      c->center[i] += (c->v[k])->coord[i];       /* dim * barycenter */
    norm += c->center[i] * c->center[i];         /* norm ^2 */
  }
  
  /* norm --> 1 */
  norm = sqrt(norm);
  for( i=0; i<dim; i++ )
    c->center[i] /= norm;

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_cone_center() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c )
     /*----------------------------------------------------------------------*/
     /* compute parameters for hat for a cone                                */
     /* (expect touching point and volume below hat)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   c   ... cone for which parameters have to be computed              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double Tf,f;                  /* (transformed) density */
  double Tderf;                 /* T'(f(x)) */
  int i;                        /* aux variable */
  int dim = GEN->dim;           /* dimension */

  double *g = GEN->g;               /* Vector g (direction of sweeping plane) */
  double *coord = GEN->tp_coord;    /* coordinates of touching point */
  double *mcoord = GEN->tp_mcoord;  /* coordinates of touching point moved into center */
  double *Tgrad = GEN->tp_Tgrad;    /* gradient of transformed density */

  /* we consider values for the PDF and of derivated numbers that are
     too small as equal 0: */
  double tolerance = TOLERANCE * GEN->pdfcenter / dim;

  /* coordinates of touching point */
  for( i=0; i<dim; i++ ) {
    coord[i] = c->tp * c->center[i];
    mcoord[i] = coord[i] + GEN->center[i];
  }

  /* density at construction point */
  if( DISTR.logpdf != NULL ) {
    c->Tfp = Tf = logPDF(mcoord);
    if (! _unur_isfinite(Tf))
      return UNUR_ERR_DISTR_DOMAIN;
  }
  else {
    f = PDF(mcoord);
    /* check density */
    if( f < tolerance )
      /* assume f = 0. */
      return UNUR_ERR_DISTR_DOMAIN;
    /* transformed density */
    c->Tfp = Tf = T(f);
  }

  /* gradient of logPDF */
  if( DISTR.dlogpdf != NULL ) {
    /* gradient of logPDF available */
    dlogPDF(Tgrad,mcoord);
  }
  else {
    /* grad( T(f(x) ) = T'(f(x)) * grad(f(x)) */
    dPDF(Tgrad,mcoord);
    Tderf = T_deriv(exp(Tf));
    for( i=0; i<dim; i++ )
      /* grad( T(f(x) ) = T'(f(x)) * grad(f(x)) */
      Tgrad[i] *= Tderf;
  }

  /* parameters alpha and beta */
  c->alpha = Tf - _unur_vector_scalar_product(dim,Tgrad,coord);
  c->beta = _unur_vector_norm(dim,Tgrad);

  /* |Tgrad| must not be too small */
  if( c->beta < tolerance )
    return UNUR_FAILURE;

  /* vector g = - grad(T(f)) / |grad(T(f))| */
  for( i=0; i<dim; i++ )
    g[i] = - Tgrad[i] / c->beta;

  /* <g,v> for each vertex v of cone and    */
  /* parameter a1 for volume of cone        */
  c->logai = c->logdetf;
  for( i=0; i<dim; i++ ) {
    c->gv[i] = _unur_vector_scalar_product(dim,g,(c->v[i])->coord);   /* <g,v> */
    if( c->gv[i] < tolerance )
      /* too small: would result in severe numerical errors */
      return UNUR_FAILURE;
    else
      c->logai -= log(c->gv[i]);
  }

  /* this is expensive for calculation for every touching point !!! */
  /* Maybe we only perform this computation when at the very end    */
  /* and use height=INFINITY for finding construction points.       */
  if (_unur_mvtdr_cone_height(gen,c) != UNUR_SUCCESS) 
    return UNUR_FAILURE;

  /* return error code */
  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_cone_params() */

/*---------------------------------------------------------------------------*/

double
_unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c )
     /*----------------------------------------------------------------------*/
     /* calculate log of volume below hat for given touching point.          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   c   ... cone for which volume below hat has to be computed         */
     /*                                                                      */
     /* return:                                                              */
     /*   success               ... logarithm of volume below hat            */
     /*   PDF(tp)==0            ... -INFINITY                                */
     /*   error(hat unbounded?) ... +INFINITY                                */
     /*----------------------------------------------------------------------*/
{
  double logH;

  /* compute parameters for cone */
  switch ( _unur_mvtdr_cone_params(gen,c) ) {
  case UNUR_SUCCESS:
    break;
  case UNUR_ERR_DISTR_DOMAIN:
    /* construction point out of support */
    return -INFINITY;
  default:
    /* something is wrong: beta = 0 and/or <g,v> <= 0 */
    return INFINITY;
  }

  /* compute log of volume below hat */
  logH = c->alpha - GEN->dim * log(c->beta) + c->logai;

  if (_unur_isfinite(c->height)) {
    /* there is a change to the paper, eq.(15):          */
    /*   'logai' does not contain '(dim-1)!' .           */
    /* thus we have to modify eq.(31) accordingly.       */
    /* Remark: this is a rather expensive computation.   */
    /* It could be omitted if 'c->beta*c->height' is     */
    /* large (not too small).                            */
    if (c->height < 1.e-50) 
      return -INFINITY;
    else
      logH += log(_unur_SF_incomplete_gamma(c->beta*c->height,(double)GEN->dim));
  }

  /* check for numerical errors (alpha or beta too small) */
  return (_unur_isfinite(logH)) ? logH : INFINITY;

} /* end of _unur_mvtdr_cone_logH() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step )
     /*----------------------------------------------------------------------*/
     /* split a cone along "oldest" edge                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   c    ... cone which has to be split                                */
     /*   step ... triangulation level                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  CONE *newc;       /* new cone */
  VERTEX *newv;      /* new vertex */
  int dim = GEN->dim;  /* dimension */
  int i;

  if (dim == 2)
    /* there is only one edge in this cone */
    newv = _unur_mvtdr_vertex_on_edge(gen,c->v);
  else
    /* find "oldest" edge, read center of edge from table (or computer) */
    newv = _unur_mvtdr_etable_find_or_insert(gen,c->v);
  if (newv==NULL) return UNUR_FAILURE;
  
  /* construct two new cones */

  /* first cone */
  newc = _unur_mvtdr_cone_new(gen);   /* new cone */
  if (newc==NULL) return UNUR_ERR_MALLOC;
  newc->level = step;                 /* triangulation level */
  for (i=0; i<dim-1; i++)
    newc->v[i] = c->v[i+1];           /* copy list of vertices to new cone */
  newc->v[dim-1] = newv;              /* add new vertex */
  newc->logdetf = c->logdetf - log(2.*newv->norm);  /* log of det of spanning vectors */
  newc->tp = c->tp;                   /* distance of touching point remains unchanged */

  /* second cone */
  c->level = step;                    /* triangulation level */
  for (i=0; i<dim-2; i++)
    c->v[i+1] = c->v[i+2];            /* shift list of vertices */
  (c->v)[dim-1] = newv;               /* add new vertex */
  c->logdetf = newc->logdetf;         /* the determinant */

  /* store maximal triangulation level for debugging */
  GEN->n_steps = _unur_max(GEN->n_steps, step); 

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_cone_split() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_triangulate( struct unur_gen *gen, int step, int all )
     /*----------------------------------------------------------------------*/
     /* make one triangulation step                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   step ... level of triangulation steps                              */
     /*   all  ... whether all (TRUE) cones have to be split                 */
     /*            or only those cones where volumes are too big             */ 
     /*                                                                      */
     /* return:                                                              */
     /*   number of new cones                                                */
     /*   -1 in case of error                                                */
     /*----------------------------------------------------------------------*/
{
  int k,nc;
  CONE *c;
  int dim = GEN->dim;  /* dimension */

  if (dim > 2) {
    /* We need a hash table for storing the edges.                                 */
    /* size of table = maximum number of vertices in current triangulation cycle.  */
    /* only required at begining of new tiangulation iteration, i.e. step = dim-1. */
    /* length of cycle dim-1.                                                      */
    /* only necessary if dim > 2                                                   */
    if( step % (dim-1) == 1 )
      if( _unur_mvtdr_etable_new(gen, _unur_mvtdr_number_vertices(gen, (step/(dim-1)+1)*(dim-1) ))
	  != UNUR_SUCCESS )
	return -1;
  }

  /*   number of cones before triangulation */
  nc = GEN->n_cone;

  /*   triangulate every cone */
  for( k=0, c=GEN->cone; k<nc; k++ ) {
    if( all ) {
      if (_unur_mvtdr_cone_split(gen,c,step) != UNUR_SUCCESS)
	return -1;
    }
    else if ( c->tp < 0. ) {
      if (_unur_mvtdr_cone_split(gen,c,step) != UNUR_SUCCESS)
	return -1;
      _unur_mvtdr_tp_find (gen,c);
      _unur_mvtdr_tp_find (gen,GEN->last_cone);
    }
    /* next cone */
    c = c->next;
  }

  /* return number of new cones */
  return (GEN->n_cone - nc);

} /* end of _unur_mvtdr_triangulate() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_height( struct unur_gen *gen, CONE *c )
     /*----------------------------------------------------------------------*/
     /* calculate height of pyramid (cone) using simplex algorithm           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   c    ... cone for which we compute height                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* working array for simplex tableau: double A[N+1][N+1]; */
#define A(i,j)  (AA[(dim+1)*(i)+(j)])
  double *AA;
  int dim = GEN->dim;

  /* coordinates of lower left and upper right vertices of domain */
#define ll(i)   (domain[2*(i)]   - GEN->center[(i)])
#define ur(i)   (domain[2*(i)+1] - GEN->center[(i)])
  double *domain;
 
  int i,j,row,ipc,ipr;
  double pc,pr,ratio;
  double sgn;

  /* bounded domain ? */
  if (! GEN->has_domain)
    /* nothing to do (c->height is set to INIFINITY at creating time) */
    return UNUR_SUCCESS;
  
  /* get rectangular domain */
  if (DISTR.domainrect == NULL) {
    _unur_error(gen->genid,UNUR_ERR_DISTR_DOMAIN,"no domain given");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  domain = DISTR.domainrect;

  /* allocate memory for simplex tableau */
  AA = _unur_xmalloc( (dim+1)*(dim+1)*sizeof(double) );

  /* set initial matrix for simplex algorithm */
  for( i=0, row=0; i<dim; i++ ) {
    /* loop over all coordinate directions */
    
    /* get sign of coordinates (these must be either all nonnegative or nonpositive)  */
    for( j=0, sgn=0.; j<dim; j++ ) {
      if( (c->v[j])->coord[i] > 0. ) {
	sgn = +1.; break;
      }
      if( (c->v[j])->coord[i] < 0. ) {
	sgn = -1.; break;
      }
    }
    /* check sign */
    if (_unur_iszero(sgn)) continue;

    /* coefficients of inequalities */
    for( j=0; j<dim; j++ )
      A(row,j) = sgn * (c->v[j])->coord[i];
    A(row,dim) = (sgn > 0.) ? ur(i) : -(ll(i));
    row++;
  }

  /* objective function */
  for( j=0; j<dim; j++ )
    A(row,j) = -(c->gv[j]);
  A(row,dim) = 0.;

  /* find maximum */
  while( 1 ) {

    /* pivot column */
    for( j=0,pc=0.,ipc=-1; j<dim; j++) {
      if( A(row,j) < pc ) {
        ipc = j;
        pc = A(row,ipc);
      }
    }

    /* stop ? */
    if( ipc == -1 ) {   /* no pivot column found */
      c->height = A(row,dim);
      break;
    }

    /* find pivot row */
    for( i=0,pr=-1.,ipr=-1; i<row; i++ ) {
      if( A(i,ipc) <= 0. ) continue;
      ratio = A(i,dim) / A(i,ipc);
      if( pr < 0 || pr > ratio ) {
        ipr = i;
        pr = ratio;
      }
    }
    
    /* unbounded ? */
    if( ipr == -1 ) {
      c->height = UNUR_INFINITY;
      break;
    }

    /* make pivot step */
    /* the rest */
    for( i=0; i<=row; i++ )
      if( i!= ipr )
        for( j=0; j<dim+1; j++ )
          if( j!= ipc ) 
            A(i,j) -= A(ipr,j) * A(i,ipc) / A(ipr,ipc);
    /* pivot column */
    for( i=0; i<=row; i++ )
      if( i != ipr )
        A(i,ipc) /= -A(ipr,ipc);
    /* pivot row */
    for( j=0; i<dim; i++ )
      if( j != ipc )
        A(ipr,j) /= A(ipr,ipc);
    /* pivot element */
    A(ipr,ipc) = 1./A(ipr,ipc);
  }

  /* free working space (simplex tableau) */
  free (AA);

  /* check result */
  if (_unur_isnan(c->height))  c->height = UNUR_INFINITY;

  /* o.k. */
  return UNUR_SUCCESS;

#undef A
#undef ll
#undef ur
} /* end of _unur_mvtdr_cone_height() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_max_gamma( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute upper bound for gamma variates                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double max, tmp;
  CONE *c;

  if (!GEN->has_domain) {
    GEN->max_gamma = INFINITY;
  }
  else {
    max = 0.;
    for (c = GEN->cone; c != NULL; c = c->next) {
      tmp = c->height * c->beta;
      max = _unur_max(max, tmp);
    }
    GEN->max_gamma = (max > 0.) ? max : INFINITY;
  }

  return UNUR_SUCCESS;
} /* end of _unur_mvtdr_max_gamma() */


/*****************************************************************************/
/*                                                                           */
/*   Optimal distance for touching points                                    */
/*                                                                           */
/*****************************************************************************/

double
_unur_mvtdr_tp_min (double t, void *p )
     /*----------------------------------------------------------------------*/
     /* volume function.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   t ... location of touching point                                   */
     /*   p ... pointer to arguments of volume functions                     */
     /*                                                                      */
     /* return:                                                              */
     /*   logarithm of the volume below the hat                              */
     /*----------------------------------------------------------------------*/
{
  /* unpack arguments */
  TP_ARG *a = p;
  /* set new construction point for hat */
  (a->c)->tp = a->t = t;
  /* compute volume below hat in cone and return result */
  a->logH = _unur_mvtdr_cone_logH (a->gen, a->c);

  /* status of result */
  switch (_unur_isinf(a->logH)) {
  case -1:
    a->logH = INFINITY;
    a->status = MVTDR_CONE_DOMAIN;
    break;
  case 1:
    a->status = MVTDR_CONE_INVALID;
    break;
  case 0:
  default:
    a->status = MVTDR_CONE_OK;
  }
  
  if( a->status != MVTDR_CONE_OK )
    /* we mark this case by setting tp = -1 */
    (a->c)->tp = -1.;

  return a->logH;
} /* end of _unur_mvtdr_tp_min() */

/*---------------------------------------------------------------------------*/

double
_unur_mvtdr_tp_min_aux(double t, void *p)
     /*----------------------------------------------------------------------*/
     /* auxiliary function to be used with _unur_util_brent().               */
     /* _unur_util_brent() maximizes functions, so we need the negative      */
     /* of the volume function.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   t ... location of touching point                                   */
     /*   p ... pointer to arguments of volume functions                     */
     /*                                                                      */
     /* return:                                                              */
     /*   negative of logarithm of the volume below the hat                  */
     /*----------------------------------------------------------------------*/
{
  return (- _unur_mvtdr_tp_min(t, p) );
}

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c )
     /*----------------------------------------------------------------------*/
     /* find optimal touching point for cone using Brent's method            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   c   ... cone for which touching point has to be computed           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_funct_generic tpaux;
  TP_ARG a[3];    /* left, middle and right point of bracket for Brent algorithm */
  int i;

  /* compute center */
  _unur_mvtdr_cone_center(gen,c);

  /* pack arguments for finding minimum */
  for (i=0; i<3; i++) { a[i].c = c; a[i].gen = gen; }

  /* find proper touching point */
  switch (_unur_mvtdr_tp_search(gen,a)) {
  case UNUR_SUCCESS:
    break;
  case UNUR_ERR_DISTR_DOMAIN:
    /* cone not in support of PDF */
    c->tp = 0.;
    c->Hi = 0.;
    c->height = 0.;
    return UNUR_ERR_DISTR_DOMAIN;
  case UNUR_FAILURE:
  default:
    /* no proper point found */
    return UNUR_FAILURE;
  }

  /* find "bracket" for Brent's algorithm */
  switch( _unur_mvtdr_tp_bracket(gen,a) ) {      /* searching for intervall that contains minimum */
  case TP_BRACKET:                 /* bracket found */
    /* make auxiliary function for Brent's algorithms */
    tpaux.f = _unur_mvtdr_tp_min_aux;
    tpaux.params = a+1;
    c->tp = _unur_util_brent( tpaux, a[0].t, a[2].t, a[1].t, FIND_TP_TOL);
    c->Hi = exp(a[1].logH);
    break;                         /* c->tp already set by tp_min() */
  case TP_LEFT:                    /* minimum in left point */
    c->tp = a[0].t;
    c->Hi = exp(a[0].logH);
    break;
  case TP_MIDDLE:                  /* minimum in middle point */
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */
    c->Hi = exp(a[1].logH);
    break;
  case TP_RIGHT:                   /* minimum in right point */
    c->tp = a[2].t;
    c->Hi = exp(a[2].logH);
    break;
  default:                         /* no proper touching point found */
    c->tp = -1.;
    return UNUR_FAILURE;
  }

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_tp_find() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_tp_search( struct unur_gen *gen ATTRIBUTE__UNUSED, TP_ARG *a )
     /*----------------------------------------------------------------------*/
     /* search for proper touching point.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   a   ... arguments for function that should be minimized            */
     /*           a[0], a[1], a[2] ... left, middle and right point of       */
     /*           bracket used for Brent's algorithm,                        */
     /*           see _unur_util_brent().                                    */
     /*           the result is stored in a[1].                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define N_STEPS  (10)  /* max number of steps for finding proper touching point */

  int i;                           /* aux counter */
  int is_unbounded_cone = FALSE;  /* boolean */
  double start = 1.;               /* starting point for search routine */

  /** search from 0 --> 1 (towards infinity) **/

  /* initialize boundary of intervall */
  a[0].t = 0.;     /* x[0] must >= 0. */
  a[1].t = start;  /* starting point for searching proper touching point */
  a[2].t = -1.;    /* not known. marked by setting to -1. */
                  
  for( i=1; i <= N_STEPS; i++ ) {
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */

    if (a[1].status == MVTDR_CONE_OK) {
      return UNUR_SUCCESS;
    }
    else if (a[1].status == MVTDR_CONE_DOMAIN) {
      /* touching point out of domain */
      break;
    }
    else {
      /* not a proper touching point */
      is_unbounded_cone = TRUE;
      a[0].t = a[1].t;
      a[1].t *= 2.;
    }
  }
      
  /** search from 1 --> 0 **/

  /* initialize boundary of intervall */
  a[0].t = 0.;        /* x[0] must >= 0. */
  a[1].t = start/2.;  /* starting point for searching proper touching point */
  a[2].t = 1.;        /* t[2] must >= t[1]. */

  for( i=0;; i++ ) {
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */

    if (a[1].status == MVTDR_CONE_OK) {
      return UNUR_SUCCESS;
    }
    else if (a[1].status == MVTDR_CONE_DOMAIN) {
      /* touching point out of domain */
      if (a[1].t < 1.e-20 ) {
	/* the bound 1.e-20 is rather arbirary */
	return (is_unbounded_cone) ? UNUR_FAILURE : UNUR_ERR_DISTR_DOMAIN;
      }
      a[2].t = a[1].t;
      a[1].t /= 10.;
    }
    else {
      if (i > N_STEPS) {
	/* no proper touching point found               */
	/* --> giving up and split cone before next try */
	return UNUR_FAILURE;
      }
      is_unbounded_cone = TRUE;
      a[2].t = a[1].t;
      a[1].t /= 2.;
    }
  }

  /* no proper touching point found */
  /* return UNUR_FAILURE; */

#undef N_STEPS
} /* end of _unur_mvtdr_tp_search() */

/*-----------------------------------------------------------------*/

int 
_unur_mvtdr_tp_bracket( struct unur_gen *gen ATTRIBUTE__UNUSED, TP_ARG *a )
     /*----------------------------------------------------------------------*/
     /* search for proper bracket for Brent's algorithm                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   a   ... arguments for function that should be minimized            */
     /*           a[0], a[1], a[2] ... left, middle and right point of       */
     /*           bracket used for Brent's algorithm,                        */
     /*           see _unur_util_brent().                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   TP_LEFT    ... minimum in left point                               */
     /*   TP_MIDDLE  ... use middle point                                    */
     /*   TP_RIGHT   ... minimum in right point                              */
     /*   TP_BRACKET ... bracket found                                       */
     /*----------------------------------------------------------------------*/
{
#define N_STEPS  (10)  /* max number of steps for finding bracket */

  int i;                 /* aux variable */
  double tleft, tright;  /* left boundary of searching region */

  /** left point of intervall **/

  /* initialize boundary of searching region and set starting point */
  tleft = a[0].t;
  a[0].t = a[1].t / 2.;

  /* search */
  for( i=1; i <= _unur_max(1,N_STEPS); i++ ) {
    _unur_mvtdr_tp_min(a[0].t, a);  /* volume function */

    if( a[0].status != MVTDR_CONE_OK ) {
      /* a[0] not a proper touching point */
      tleft = a[0].t;                     /* change boundary of searching region */
      a[0].t += (a[1].t - a[0].t) / 2.;   /* try another one */
    }

    else if( a[0].logH <= a[1].logH ) {
      /* a[0] is proper touching point, but ... */
      a[2].t = a[1].t; a[2].logH = a[1].logH; a[2].status = MVTDR_CONE_OK;
      a[1].t = a[0].t; a[1].logH = a[0].logH; a[1].status = MVTDR_CONE_OK;
      a[0].t = tleft + (a[0].t - tleft)*0.5;
    }
    else  /* all right: a[0].logH > a[1].logH */
      break;
  }

  /* search successful ? */
  if( a[0].status != MVTDR_CONE_OK )
    /* no proper touching point on left side --> use middle point */
    return TP_MIDDLE;
  if( a[0].logH <= a[1].logH )
    /* vol(left) <= vol(middle) */
    return TP_LEFT;
  
  /** right point of intervall **/

  /* initialize a[2] if necessary */
  if( a[2].t < 0. )
    a[2].t = 1.1 * a[1].t;
  tright = -1.;    /* no right boundary known yet */
  tleft = a[1].t;

  /* search */
  for( i=1; i <= _unur_max(1,N_STEPS); i++ ) {
    _unur_mvtdr_tp_min(a[2].t, a+2);  /* volume function */
    if( a[2].status != MVTDR_CONE_OK ) {
      /* a[2] not a proper touching point */
      tright = a[2].t;
      a[2].t = (tleft + a[2].t) * 0.5;   /* try another one */
    }
    else if( a[2].logH <= a[1].logH ) {
      /* move right */
      tleft = a[2].t;
      a[2].t = (tright < 0.) ? a[2].t * 2. : (tright + a[2].t) * 0.5;
    }
    else    /* all right: vol(right) > vol(middle) */
      break;
  }

  /* search successful ? */
  if( a[2].status != MVTDR_CONE_OK )
    /* no proper touching point on right side --> use middle point */
    return TP_MIDDLE; 
  if( a[2].logH <= a[1].logH )
    /* f(right) <= f(middle) */
    return TP_RIGHT;

  /* we have found a bracket */
  return TP_BRACKET;

#undef N_STEPS
} /* end of _unur_mvtdr_tp_bracket() */

/*****************************************************************************/
/*                                                                           */
/*   VERTICES.                                                               */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_initial_vertices( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get vertices of initial cones.                                       */
     /* these are the vertices (0,...,0,+/- 1, 0, ..., 0).                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  VERTEX *vt;
  int i,k;
  double d;
  double *domain;

  /* get rectangular domain */
  domain = ( (GEN->has_domain && DISTR.domainrect != NULL) 
	     ? DISTR.domainrect : NULL );

  /* unit vectors e_k  and -e_k */
  for ( d=1.; d > -2.; d -= 2.) {
    /* '+'-sign and '-'-sign */
    for( k=0; k<GEN->dim; k++ ) {
      vt = _unur_mvtdr_vertex_new(gen);
      if (vt==NULL) return UNUR_FAILURE;

      for( i=0; i<GEN->dim; i++ ) {
	/* coordinates */
	(vt->coord)[i] = (i==k) ? d : 0.;
      }
      /* all vectors have norm 1. */
      vt->norm = 1.;

      /* vertices outside domain of PDF get a negative index */
      if ( domain ) {
	if ( (d<0 && _unur_FP_equal(GEN->center[k],domain[2*k]))   ||
	     (d>0 && _unur_FP_equal(GEN->center[k],domain[2*k+1])) )
	  vt->index = -vt->index -1;
      }
    }
  }

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_initial_vertices() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_vertex_new( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* allocate new vertex and append it to linked list of all vertices.    */
     /* increment counter for vertices.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to newly allocated vertex                                  */
     /*----------------------------------------------------------------------*/
{
  VERTEX *v;

  /* allocate memory */
  v = malloc(sizeof(VERTEX));
  if (v==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }

  /* insert into list of cones */
  if (GEN->vertex == NULL) {
    GEN->last_vertex = GEN->vertex = v;
  }
  else {
    GEN->last_vertex = GEN->last_vertex->next = v;
  }
  v->next = NULL;
  
  /* coordinates of vertex */
  v->coord = malloc(GEN->dim * sizeof(double));
  if (v->coord==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }

  /* index of vertex */
  v->index = GEN->n_vertex;
  /* and update counter */
  ++(GEN->n_vertex);

  /* return pointer to next vertex */
  return GEN->last_vertex;

} /* end of _unur_mvtdr_vertex_new() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl )
     /*----------------------------------------------------------------------*/
     /* compute new vertex on edge (i.e., its barycenter)                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vl  ... arraqy of the two end vertices of the edge.                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new vertex                                              */
     /*----------------------------------------------------------------------*/
{
  int i;
  VERTEX *newv;          /* pointer to new vertex */

  /* get an empty vertex */
  newv = _unur_mvtdr_vertex_new(gen);
  if (newv==NULL) return NULL;

  /* barycenter of edge */
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] =
      0.5 * ( ((vl[0])->coord)[i] + ((vl[1])->coord)[i] );

  /* norm */
  newv->norm = _unur_vector_norm(GEN->dim, newv->coord);

  /* norm --> 1 */
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] /= newv->norm;

  return newv;

} /* end of _unur_mvtdr_vertex_on_edge() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_number_vertices( struct unur_gen *gen, int level )
     /*----------------------------------------------------------------------*/
     /* calculate (approximate) number of vertices in a given triangulation  */
     /* level (when all initial cones are splitted).                         */
     /* These numbers are used for the size of the hash table for edges.     */
     /* (WARNING! The numbers are found by computer experiments.)            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   level ... number of triangulation steps                            */
     /*                                                                      */
     /* return:                                                              */
     /*   number of vertices                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  if (level < 0 || GEN->dim < 2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return -1;  /** errorcode **/
  }

  switch (GEN->dim) {

  case 2: {
    return (1 << (level+2));
  }
  case 3: {
    static int nv[]={ 6, 13, 18, 40, 66,142,258,538,1026,2098,4098,8290,16386,32962,65538,131458,262146};
    return nv[_unur_min(level,16)];
  }
  case 4: {
    static int nv[]={ 8, 19, 25, 32, 80,128,192,456, 824,1408,3120,5968,11008,23264,45600, 87552};
    return nv[_unur_min(level,15)];
  }
  case 5: {
    static int nv[]={10, 26, 33, 41, 50,140,220,321, 450,1186,2158,3636, 5890,13970,27130};
    return nv[_unur_min(level,14)];
  }
  case 6: {
    static int nv[]={12, 34, 42, 51, 61, 72,224,348, 501, 681, 912,2660, 4896, 8254};
    return nv[_unur_min(level,13)];
  }
  case 7: {
    static int nv[]={14, 43, 52, 62, 73, 85, 98,336, 518, 743, 985,1289, 1666};
    return nv[_unur_min(level,12)];
  }
  case 8: {
    static int nv[]={16, 53, 63, 74, 86, 99,113,128, 480, 736,1059};
    return nv[_unur_min(level,10)];
  }
  case 9: {
    static int nv[]={18, 64, 75, 87,100,114,129,145, 162, 660};
    return nv[_unur_min(level,9)];
  }
  case 10: {
    static int nv[]={20, 76, 88,101,115,130,146,163, 181, 200};
    return nv[_unur_min(level,9)];
  }
  case 11: {
    static int nv[]={22, 89,102,116,131,147,164,182, 201, 221, 242};
    return nv[_unur_min(level,10)];
  }
  default: { /* dim >= 12 */
    static int nv[]={24,103,117,132,148,165,183,202, 222, 243, 265, 288};
    return nv[_unur_min(level,11)];
  }
  }

} /* end of _unur_mvtdr_number_vertices() */


/*****************************************************************************/
/*                                                                           */
/*   hash table for storing EDGES.                                           */
/*                                                                           */
/*****************************************************************************/

/* hash function for edge table */
#define _unur_mvtdr_etable_hash(x,y)  ( (3*((x)+(y))/2) % GEN->etable_size )

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_etable_new( struct unur_gen *gen, int size )
     /*----------------------------------------------------------------------*/
     /* make a new hash table.                                               */
     /* destroy old table if it exists.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   size ... size of hash table                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int n;

  /* first clear hash table (if necessary) */
  _unur_mvtdr_etable_free(gen);

  /* set size of edge table */
  GEN->etable_size = size;

  /* make root */
  GEN->etable = malloc( size * sizeof(E_TABLE*) );
  if (GEN->etable==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; }

  /* initialize table */
  for (n = 0; n< size; n++) 
    GEN->etable[n] = NULL;
  
  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_etable_new() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_etable_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* free hash table.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int i;
  E_TABLE *et, *et_next;

  if( GEN->etable == NULL )
    /* nothing to do */
    return;

  /* clear all branches */
  for( i=0; i<GEN->etable_size; i++ ) {
    /* free branch in table (linked list) */
    for (et = GEN->etable[i]; et != NULL; et = et_next) {
      et_next = et->next;
      free (et);
    }
  }

  /* clear root */
  free( GEN->etable );
  GEN->etable = NULL;
  GEN->etable_size = 0;

} /* end if _unur_mvtdr_etable_free() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx )
     /*----------------------------------------------------------------------*/
     /* search for an edge in the hash table.                                */
     /* the edge is given by its end vertices vidx[0] and vidx[1].           */
     /* (if vidx is an ordered array of the vertices of a cone then this     */
     /* edge is the "oldest" edge of that cone (i.e., the end vertices have  */
     /* the smallest indices of all vertices of the cone.)                   */
     /*                                                                      */
     /* if the edge is found, a pointer to the vertex that corresponds to    */
     /* the barycenter of the edge.                                          */
     /*                                                                      */
     /* if the edge is not found in the table, it is created and inserted    */
     /* into the table, the barycenter is computed.                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   vidx ... array of pointers to vertices                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  E_TABLE *pet, *pet_last;  /* table entry */
  int idx[2];               /* store indices */
  int hidx;                 /* hash number */

  /* check pointer */
  CHECK_NULL(GEN->etable,NULL);

  /* hash number */
  idx[0] = vidx[0]->index;
  idx[1] = vidx[1]->index;
  hidx = _unur_mvtdr_etable_hash(idx[0],idx[1]);

  /* get branch of hash table */
  pet = pet_last = *(GEN->etable + hidx);

  /* now find entry in branch */
  while( pet != NULL ) {
    if( pet->index[0] == idx[0] && pet->index[1] == idx[1] )
      break;   /* found ! */
    pet_last = pet;
    pet =  pet->next;
  }

  if( pet == NULL ) {
    /* we have not found the index */

    /* new entry in hash table */
    pet = malloc( sizeof(E_TABLE) );
    if (pet==NULL) {
      _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }

    pet->next = NULL;
    if (pet_last == NULL)
      *(GEN->etable + hidx) = pet;
    else
      pet_last->next = pet;

    /* insert data of new edge */
    /* indices of incident vertices */
    pet->index[0] = idx[0];
    pet->index[1] = idx[1];

    /* compute new vertex */
    pet->vertex = _unur_mvtdr_vertex_on_edge(gen,vidx);
  }

  /* return pointer to (new) edge */
  return pet->vertex;

} /* end of _unur_mvtdr_etable_find_or_insert() */

/*****************************************************************************/

int
_unur_mvtdr_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create guide table.                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int j;
  CONE *c;

  /* memory for the guide table */
  GEN->guide_size = GEN->n_cone * GUIDE_TABLE_SIZE;
  GEN->guide = malloc (GEN->guide_size * sizeof(CONE*));
  if (GEN->guide==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; }
  /* initialize table */
  for( j = 0; j < GEN->guide_size ; j++ )
    GEN->guide[j] = NULL;

  /* make table */
  for( c=GEN->cone, j=0; c!=NULL && j<GEN->guide_size; j++ ) {
    while( c->Hsum / GEN->Htot < (double) j / GEN->guide_size )
      c=c->next;
    (GEN->guide)[j] = c;
    if( c == GEN->last_cone ) break;
  }

  /* is there an error ? */
  if( j<GEN->guide_size )
    /* this should not happen */
    for( ; j<GEN->guide_size; j++ )
      (GEN->guide)[j] = GEN->last_cone;

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_make_guide_table() */

/*---------------------------------------------------------------------------*/
