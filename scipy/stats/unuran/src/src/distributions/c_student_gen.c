/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_student_gen.c                                              *
 *                                                                           *
 *   Special generators for Student's t distribution                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int student_trouo_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define MAX_gen_params  6      /* maximal number of parameters for generator */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define nu (DISTR.params[0])    /* shape (degrees of freedom) */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_student_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Student's distribution              */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Polar Method */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_student_tpol );
    return UNUR_SUCCESS;

  case 2:  /* Ratio of Uniforms */
    if (par!=NULL && par->distr->data.cont.params[0] < 1.) {   /* nu < 1 */
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    /* nu >= 1 !!!! */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_student_trouo );
    return student_trouo_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_chi_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Student's t Distribution: Polar Method                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from Student's t distribution with  *    
 *               parameters nu > 0.                                          *
 *                                                                           *
 * REFERENCE:  - R.W. Bailey (1994): Polar generation of random variates     *
 *               with the t-distribution,                                    *
 *               Mathematics of Computation 62, 779-781.                     *
 *                                                                           *
 * Implemented by F. Niederl, 1994                                           *
 *****************************************************************************
 *                                                                           *
 * The polar method of Box/Muller for generating Normal variates is adapted  *
 * to the Student-t distribution. The two generated variates are not         *
 * independent and the expected no. of uniforms per variate is 2.5464.       *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

double
_unur_stdgen_sample_student_tpol( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double u,v,w;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {
    u = 2. * uniform() - 1.;
    v = 2. * uniform() - 1.;
    w = u * u + v * v;
  } while (w > 1.);

  return(u * sqrt( nu * ( exp(- 2. / nu * log(w)) - 1.) / w));
  /* -X- end of generator code -X- */
} /* end of _unur_stdgen_sample_student_tpol() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Student's t Distribution: Ratio of Uniforms                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from Student's t distribution with  *    
 *               parameters nu >= 1.                                         *
 *                                                                           *
 * REFERENCE:  - A.J. Kinderman, J.F. Monahan (1980):                        *
 *               New methods for generating Student's t and gamma variables, *
 *               Computing 25, 369-377.                                      *
 *                                                                           *
 * Implemented by R. Kremer, 1990                                            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define c       (GEN->gen_param[0])
#define e       (GEN->gen_param[1])
#define p       (GEN->gen_param[2])
#define q       (GEN->gen_param[3])
#define r       (GEN->gen_param[4])
#define vm      (GEN->gen_param[5])
/*---------------------------------------------------------------------------*/

inline static int
student_trouo_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (nu < 1.) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return UNUR_ERR_GEN_CONDITION;
  }

  r = 1. / nu;
  p = 1. / (1. + r);
  q = -0.25 * (nu + 1.);
  c = 4. * pow(p, q);
  e = 16. / c;
  vm = (nu>1.0) ? sqrt(p+p) * pow( (1.-r)*p, 0.25*(nu-1.) ) : 1.;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of student_trouo_init() */

double
_unur_stdgen_sample_student_trouo( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double tru,u,v;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {

    /* step 1 */
    u = uniform();

    /* step 2 */
    v = uniform();
    v = vm * (v + v - 1.);
    tru = v / u;

    /* step 3 */
    if ( c * u <= 5. - tru * tru) 
      break;
    if (nu >= 3.) 
      if (u * (tru * tru + 3.) >= e) 
	continue;  /* goto 1 */      /* step 4 */
    if ( u <= pow(1. + tru * tru * r, q))
      break;
  }

  return tru;

} /* end of _unur_stdgen_sample_student_trouo() */

/*---------------------------------------------------------------------------*/
#undef c
#undef e
#undef p
#undef q
#undef r
#undef vm
/*---------------------------------------------------------------------------*/
