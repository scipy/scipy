/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      counturn.c                                                   *
 *                                                                           *
 *   Count used uniform random numbers                                       *
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
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
static char test_name[] = "CountURN";

/*---------------------------------------------------------------------------*/
/* common variables                                                          */

static long urng_counter = 0;                /* count uniform random numbers */

/*---------------------------------------------------------------------------*/
/* wrapper for uniform random number generator that performs counting        */
/*---------------------------------------------------------------------------*/
#ifdef UNUR_URNG_UNURAN
/*---------------------------------------------------------------------------*/

static double (*urng_to_use)(void*);         /* pointer to real uniform RNG  */
/*---------------------------------------------------------------------------*/

static double
_urng_with_counter(void *params)
     /*----------------------------------------------------------------------*/
     /* wrapper for uniform random number generator that performs counting   */
     /*----------------------------------------------------------------------*/
{
  ++urng_counter;
  return urng_to_use(params);
} /* end of urng_with_counter() */

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_UNURAN */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_count_urn( struct unur_gen *gen, int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* count used uniform random numbers                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of used uniform random numbers                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  long j;
  UNUR_URNG *urng_aux;

  /* check arguments */
  _unur_check_NULL(test_name,gen,-1);

  /* reset counter */
  urng_counter = 0;

  /* save auxilliary generator */
  urng_aux = gen->urng_aux;

  /* exchange pointer to uniform rng with counting wrapper */
#ifdef UNUR_URNG_UNURAN
  urng_to_use = gen->urng->sampleunif;
  gen->urng->sampleunif = _urng_with_counter;
  if (gen->urng_aux) gen->urng_aux = gen->urng;
#else
  /* no counter available */
  if (verbosity)
    fprintf(out,"\nCOUNT: ---  (cannot count URNs)\n");
  return -1;
#endif

  /* run generator */
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    for( j=0; j<samplesize; j++ )
      _unur_sample_discr(gen);
    break;

  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(gen);
    break;

  case UNUR_METH_VEC: 
    { /* we need an array for the vector */
      double *vec;
      int dim;
      dim = unur_get_dimension(gen);
      vec = _unur_xmalloc( dim * sizeof(double) );
      for( j=0; j<samplesize; j++ )
	_unur_sample_vec(gen,vec);
      free(vec);
    }
    break;

  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"method unknown!");
    return -1;
  }

  /* reset pointer to uniform rng */
#ifdef UNUR_URNG_UNURAN
  gen->urng->sampleunif = urng_to_use;
#endif

  /* restore auxilliary generator */
  gen->urng_aux = urng_aux;

  /* print result */
  if (verbosity) {
    fprintf(out,"\nCOUNT: %g urng per generated number (total = %ld)\n",
	   ((double)urng_counter)/((double) samplesize),urng_counter);
  }

  /* return total number of used uniform random numbers */
  return urng_counter;

} /* end of unur_test_count_urn() */

/*---------------------------------------------------------------------------*/
