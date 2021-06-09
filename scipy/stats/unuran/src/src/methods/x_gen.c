/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen.c                                                           *
 *                                                                           *
 *   miscelleanous routines for manipulation generator objects               *
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

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/matr.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/dgt.h>
#include <methods/dstd.h>
#include <methods/dstd_struct.h>
#include <methods/hinv.h>
#include <methods/mixt.h>
#include <methods/mixt_struct.h>
#include <methods/ninv.h>
#include <methods/pinv.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Call Init, Sampling, and Free functions                                **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_gen *unur_init( struct unur_par *par )
{                
  _unur_check_NULL(NULL,par,NULL);
  return (par->init(par));
} /* end of unur_init() */

/*---------------------------------------------------------------------------*/

int unur_reinit( struct unur_gen *gen )
{
  int status = UNUR_SUCCESS;
  _unur_check_NULL(NULL,gen,UNUR_ERR_NULL);

  if (gen->reinit) {
    status = gen->reinit(gen);
    if (status == UNUR_SUCCESS) return status;
  }
  else {
    _unur_error(gen->genid,UNUR_ERR_NO_REINIT,"");
    status = UNUR_ERR_NO_REINIT;
  }

  /* error: change sampling routine */
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    gen->sample.discr = _unur_sample_discr_error;
    break;
  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
    gen->sample.cont = _unur_sample_cont_error;
    break;
  case UNUR_METH_VEC:
  case UNUR_METH_CVEMP:
    gen->sample.cvec = _unur_sample_cvec_error;
    break;
  case UNUR_METH_MAT:
    gen->sample.matr = _unur_sample_matr_error;
    break;
  default:
    _unur_error("reinit",UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  }

  return status;
} /* end of unur_reinit() */

/*---------------------------------------------------------------------------*/

int
unur_sample_discr( struct unur_gen *gen )
{
  CHECK_NULL(gen,0);
  return (gen->sample.discr(gen));
} /* end of unur_sample_discr() */

double
unur_sample_cont( struct unur_gen *gen )
{
  CHECK_NULL(gen,INFINITY);
  return (gen->sample.cont(gen));
} /* end of unur_sample_cont() */

int
unur_sample_vec( struct unur_gen *gen, double *vector )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  return (gen->sample.cvec(gen,vector));
} /* end of unur_sample_vec() */

int
unur_sample_matr( struct unur_gen *gen, double *matrix )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  return (gen->sample.matr(gen,matrix));
} /* end of unur_sample_matr() */

/*---------------------------------------------------------------------------*/
/* Estimate quantiles                                                        */

double
unur_quantile ( struct unur_gen *gen, double U )
{
  /* check arguments */
  CHECK_NULL(gen,FALSE);

  /* Remark:
   * We DO NOT check the argument U here 
   * (i.e. whether 0<=U<=1 holds)
   */
  switch (gen->method) {
  case UNUR_METH_HINV:
    return unur_hinv_eval_approxinvcdf(gen,U);

  case UNUR_METH_NINV:
    return unur_ninv_eval_approxinvcdf(gen,U);

  case UNUR_METH_PINV:
    return unur_pinv_eval_approxinvcdf(gen,U);

  case UNUR_METH_CSTD:
    if (((struct unur_cstd_gen*)gen->datap)->is_inversion)
      return unur_cstd_eval_invcdf(gen,U);
    break;

  case UNUR_METH_MIXT:
    if (((struct unur_mixt_gen*)gen->datap)->is_inversion)
      return unur_mixt_eval_invcdf(gen,U);
    break;

  case UNUR_METH_DGT:
    return ((double) unur_dgt_eval_invcdf(gen,U));

  case UNUR_METH_DSTD:
    if (((struct unur_dstd_gen*)gen->datap)->is_inversion)
      return unur_dstd_eval_invcdf(gen,U);
    break;
  }
  
  /* default: */
  _unur_error(gen->genid,UNUR_ERR_NO_QUANTILE,"");
  return UNUR_INFINITY;

} /* end of unur_quantile() */

/*---------------------------------------------------------------------------*/

int
_unur_gen_is_inversion ( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check for type of generator object                                   */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,FALSE);

  switch (gen->method) {
  case UNUR_METH_HINV:
  case UNUR_METH_NINV:
  case UNUR_METH_PINV:
  case UNUR_METH_DGT:
    return TRUE;

  case UNUR_METH_CSTD:
    return (((struct unur_cstd_gen*)gen->datap)->is_inversion);

  case UNUR_METH_MIXT:
    return (((struct unur_mixt_gen*)gen->datap)->is_inversion);

  default:
    return FALSE;
  }
} /* end of _unur_gen_is_inversion() */

/*---------------------------------------------------------------------------*/
/* aux routines when no sampling routine is available                         */

int
_unur_sample_discr_error( struct unur_gen *gen ATTRIBUTE__UNUSED )
{
  unur_errno = UNUR_ERR_GEN_CONDITION;
  return 0;
} /* end of _unur_sample_discr_error() */

double
_unur_sample_cont_error( struct unur_gen *gen ATTRIBUTE__UNUSED )
{
  unur_errno = UNUR_ERR_GEN_CONDITION;
  return INFINITY;
} /* end of _unur_sample_cont_error() */

int
_unur_sample_cvec_error( struct unur_gen *gen, double *vec )
{ 
  int d;
  unur_errno = UNUR_ERR_GEN_CONDITION;
  for (d=0; d<(gen->distr->dim); d++) vec[d] = INFINITY;
  return UNUR_FAILURE;
} /* end of _unur_sample_cvec_error() */

int
_unur_sample_matr_error( struct unur_gen *gen, double *mat )
{ 
  int n_rows, n_cols, dim, j;
  
  unur_errno = UNUR_ERR_GEN_CONDITION;
  unur_distr_matr_get_dim(gen->distr, &n_rows, &n_cols );
  dim = n_rows * n_cols;
  for (j=0; j<dim; j++)
    mat[j] = INFINITY;
  return UNUR_FAILURE;
} /* end of _unur_sample_matr_error() */

/*---------------------------------------------------------------------------*/

void
unur_free( struct unur_gen *gen )
{                
  if (gen) gen->destroy(gen);
} /* end of unur_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Get data about generator object                                        **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
const char *
unur_gen_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* return pointer to character string that contains information about   */
     /* the given generator object.                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   help ... whether to print additional comments                      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to character string                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#ifdef UNUR_ENABLE_INFO

  /* check arguments */
  _unur_check_NULL("",gen,NULL);

  if (gen->info) {
    /* prepare generator object for creating info string */
    if (gen->infostr == NULL) 
      /* either allocate memory block */
      gen->infostr = _unur_string_new();
    else 
      /* or clear string object (i.e. reset pointer) */
      _unur_string_clear(gen->infostr);

    /* create info string */
    gen->info((struct unur_gen*) gen, help);

    /* return info string */
    return gen->infostr->text;
  }
  else {
    return NULL;
  }
#else

  return "INFO string not enable";

#endif
} /* end of unur_gen_info() */
/*---------------------------------------------------------------------------*/

int
unur_get_dimension( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get dimension of generator for multivariate distribution             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   dimension of distribution                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  return (gen->distr->dim);
} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/

const char *
unur_get_genid( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get generator id                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator id                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->genid;
} /* end of unur_get_genid() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_get_distr( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object from generator object             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->distr;
} /* end of unur_get_distr() */

/*---------------------------------------------------------------------------*/

int 
unur_set_use_distr_privatecopy( struct unur_par *par, int use_privatecopy )
     /*----------------------------------------------------------------------*/
     /* Set flag whether the generator object should make a private copy or  */
     /* just stores the pointer to the distribution object                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par             ... pointer to parameter object                    */
     /*   use_privatecopy ... TRUE = use private copy                        */
     /*                       FALSE = store pointer to given distr. object   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*                                                                      */
     /* WARNING!                                                             */
     /* Using a pointer to the external distribution object instead of a     */
     /* private copy must be applied with EXTREME CARE!                      */
     /* When the distrubtion object is changed or freed then the generator   */
     /* object does not work any more, might case a segmentation fault, or   */
     /* (even worse) produces garbage.                                       */
     /* On the other hand, when the generator object is initialized or used  */
     /* to draw a random sampling the distribution object may be changed.    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL("",par,UNUR_ERR_NULL);

  par->distr_is_privatecopy = use_privatecopy;
  return UNUR_SUCCESS;
} /* end of unur_set_use_distr_privatecopy() */


/*****************************************************************************/
/**                                                                         **/
/**  Copy (clone) generator object                                          **/
/**                                                                         **/
/*****************************************************************************/

struct unur_gen *
unur_gen_clone( const struct unur_gen *gen )
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
  /* check arguments */
  _unur_check_NULL( "Clone", gen, NULL );
  _unur_check_NULL( "Clone", gen->clone, NULL );

  return (gen->clone(gen));
} /* end of unur_gen_clone() */


/*****************************************************************************/
/**                                                                         **/
/**  Create and free parameter objects                                      **/
/**                                                                         **/
/*****************************************************************************/

struct unur_par *
_unur_par_new( size_t s)
     /*----------------------------------------------------------------------*/
     /* create new parameter object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   s ... size of data structure                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par = _unur_xmalloc( sizeof(struct unur_par) );
  par->datap = _unur_xmalloc(s);
  par->s_datap = s;

  /* set defaults for distribution object */
  par->distr_is_privatecopy = TRUE;   /* use private copy of distribution object */

  return par;
} /* end of _unur_par_new() */

/*---------------------------------------------------------------------------*/

struct unur_par *
_unur_par_clone( const struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* copy parameter object                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter object                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *clone;

  _unur_check_NULL("clone", par, NULL);

  clone = _unur_xmalloc( sizeof(struct unur_par) );
  memcpy (clone, par, sizeof(struct unur_par));

  clone->datap = _unur_xmalloc(par->s_datap);
  memcpy (clone->datap, par->datap, par->s_datap);

  return clone;
} /* end of unur_par_free() */

/*---------------------------------------------------------------------------*/

void 
unur_par_free( struct unur_par *par)
     /*----------------------------------------------------------------------*/
     /* free parameter object                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter object                                */
     /*----------------------------------------------------------------------*/
{
  _unur_check_NULL("free", par, RETURN_VOID );
  _unur_par_free(par);
} /* end of unur_par_free() */

/*****************************************************************************/
/**                                                                         **/
/**  Create, copy (clone) and free generator objects                        **/
/**                                                                         **/
/*****************************************************************************/

struct unur_gen *
_unur_generic_create( struct unur_par *par, size_t s )
     /*----------------------------------------------------------------------*/
     /* create new generic generator object                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   s   ... size of data structure                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* allocate memory for generator object */
  gen = _unur_xmalloc( sizeof(struct unur_gen) );
  gen->datap = _unur_xmalloc(s);
  gen->s_datap = s;

  /* copy distribution object into generator object */
  gen->distr_is_privatecopy = par->distr_is_privatecopy;
  if (gen->distr_is_privatecopy) 
    gen->distr = (par->distr) ? _unur_distr_clone(par->distr) : NULL;
  else
    gen->distr = (struct unur_distr *) par->distr;

  /* initialize function pointers */
  gen->destroy = NULL;              /* destructor      */ 
  gen->clone = NULL;                /* clone generator */
  gen->reinit = NULL;               /* reinit routine  */ 

  /* copy some parameters into generator object */
  gen->method = par->method;        /* indicates method and variant          */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */
  gen->urng_aux = par->urng_aux;    /* pointer to auxilliary URNG            */

  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_list = NULL;         /* no auxilliary generator objects       */
  gen->n_gen_aux_list = 0;

  /* status of generator object */
  gen->status = UNUR_FAILURE;       /* not successfully created yet          */

#ifdef UNUR_ENABLE_INFO
  gen->infostr = NULL;              /* pointer to info string                */
  gen->info = NULL;                 /* routine that return info string       */
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_generic_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_generic_clone( const struct unur_gen *gen, const char *type )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generic generator object                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   type ... type of generator (string)                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *clone;

  /* allocate memory for generator object and copy main part */
  clone = _unur_xmalloc( sizeof(struct unur_gen) );
  memcpy( clone, gen, sizeof(struct unur_gen) );
  clone->datap = _unur_xmalloc(gen->s_datap);
  memcpy (clone->datap, gen->datap, gen->s_datap);

  /* set generator identifier */
  clone->genid = _unur_set_genid(type);

#ifdef UNUR_ENABLE_INFO
  /* do not copy pointer to info string */
  clone->infostr = NULL;
#endif

  /* copy distribution object into generator object */
  clone->distr_is_privatecopy = gen->distr_is_privatecopy;
  if (clone->distr_is_privatecopy) 
    clone->distr = (gen->distr) ? _unur_distr_clone(gen->distr) : NULL;
  else
    clone->distr = gen->distr;

  /* auxiliary generators */
  if (gen->gen_aux)
    clone->gen_aux = _unur_gen_clone( gen->gen_aux );
  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    clone->gen_aux_list = _unur_gen_list_clone( gen->gen_aux_list, gen->n_gen_aux_list );
    clone->n_gen_aux_list = gen->n_gen_aux_list;
  }

  /* finished clone */
  return clone;
} /* _unur_generic_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_generic_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*----------------------------------------------------------------------*/
{ 
  if (gen->gen_aux)
    _unur_free(gen->gen_aux);

  if (gen->gen_aux_list && gen->n_gen_aux_list)
    _unur_gen_list_free( gen->gen_aux_list, gen->n_gen_aux_list );

  if (gen->distr_is_privatecopy && gen->distr)
    _unur_distr_free( gen->distr );

  _unur_free_genid(gen);
  COOKIE_CLEAR(gen);
  free(gen->datap);

#ifdef UNUR_ENABLE_INFO
  /* free pointer to info string */
  if (gen->infostr) _unur_string_free(gen->infostr);  
#endif

  free(gen);
} /* end of _unur_generic_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Set and copy (clone) arrays of auxiliary generators                    **/
/**                                                                         **/
/*****************************************************************************/

struct unur_gen ** 
_unur_gen_list_set( struct unur_gen *gen, int n_gen_list )
     /*----------------------------------------------------------------------*/
     /* set all entries in list to same generator object 'gen'               */
     /*                                                                      */
     /* IMPORTANT: Be careful when using this call. When the resulting array */
     /*   is stored in some multivariate generator object then 'gen'         */
     /*   _must not_  be used any more after this call!                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_gen_list ... length of array of generator objects                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to list of generator objects                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen **gen_list;
  int i;

  /* check arguments */
  _unur_check_NULL( "gen_list_set", gen, NULL );

  if (n_gen_list < 1) {
    _unur_error("gen_list_set",UNUR_ERR_PAR_SET,"dimension < 1");
    return NULL;
  }
  
  /* allocate memory for array */
  gen_list = _unur_xmalloc (n_gen_list * sizeof(struct unur_gen *));
  
  /* copy pointer */
  for (i=0; i<n_gen_list; i++)
    gen_list[i] = gen;

  return gen_list;

} /* end of _unur_gen_list_set() */

/*---------------------------------------------------------------------------*/

struct unur_gen ** 
_unur_gen_list_clone( struct unur_gen **gen_list, int n_gen_list )
     /*----------------------------------------------------------------------*/
     /* clone list of generator objects                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen_list   ... pointer to array of generator objects               */
     /*   n_gen_list ... length of array of generator objects                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of list of generator objects                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen **clone_list;
  int i;

  /* check arguments */
  _unur_check_NULL( "gen_list_clone", gen_list, NULL );

  if (n_gen_list < 1) {
    _unur_error("gen_list_clone",UNUR_ERR_PAR_SET,"dimension < 1");
    return NULL;
  }

  for (i=0; i<n_gen_list; i++)
    _unur_check_NULL( "gen_list_clone", *(gen_list+i), NULL );

  /* allocate memory for array */
  clone_list = _unur_xmalloc (n_gen_list * sizeof(struct unur_gen *));

  /* make copy of generator objects */
  /* There are (should be) only two possibilities: 
     either all entries in the array point to the same generator object;
       (set by _unur_gen_list_set() call)
     or each entry has its own copy of some generation object.
       (set by _unur_gen_list_clone() call)
  */

  if (n_gen_list > 1 && gen_list[0] == gen_list[1]) {
      clone_list[0] = _unur_gen_clone( gen_list[0] );
      for (i=0; i<n_gen_list; i++)
	clone_list[i] = clone_list[0];
  }

  else {
    for (i=0; i<n_gen_list; i++)
      clone_list[i] = _unur_gen_clone( gen_list[i] );
  }
  
  return clone_list;

} /* end of _unur_gen_list_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_gen_list_free( struct unur_gen **gen_list, int n_gen_list )
     /*----------------------------------------------------------------------*/
     /* free list of generator objects                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen_list   ... pointer to array of generator objects               */
     /*   n_gen_list ... length of array of generator objects                */
     /*----------------------------------------------------------------------*/
{
  int i, i2, imax;

  /* check arguments */
  if (gen_list==NULL) 
    /* nothing to do */
    return; 

  if (n_gen_list < 1) {
    _unur_error("gen_list_free",UNUR_ERR_PAR_SET,"dimension < 1");
    return;
  }

  /* free memory */
  /* There are (should be) only two possibilities: 
     either all entries in the array point to the same generator object;
       (set by _unur_gen_list_set() call)
     or each entry has its own copy of some generation object.
       (set by _unur_gen_list_clone() call)
  */
  i2 = (n_gen_list>1) ? 1 : 0; 
  imax = (gen_list[0] == gen_list[i2]) ? 1 : n_gen_list;
  for (i=0; i<imax; i++)
    if (gen_list[i]) _unur_free(gen_list[i]);
  free (gen_list);

} /* end of _unur_gen_list_clone() */

/*---------------------------------------------------------------------------*/
