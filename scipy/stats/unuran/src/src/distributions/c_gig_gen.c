/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gig_gen.c                                                  *
 *                                                                           *
 *   Special generators for GIG distribution                                 *
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

inline static int gig_gigru_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params 10      /* maximal number of parameters for generator */

#define theta  (DISTR.params[0])    /* shape */
#define omega  (DISTR.params[1])    /* scale */
#define eta    (DISTR.params[2])    /* shape */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_gig_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for gig distribution                    */
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
  case 1:  /* Ratio of Uniforms */
    if (par!=NULL && par->distr->data.cont.params[0] <= 0.) {    /* theta <= 0 */
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    /* theta > 0 ! */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_gig_gigru );
    return gig_gigru_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_gig_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Generalized Inverse Gaussian (GIG) Distribution: Ratio of Uniforms        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the reparameterized            *
 *               Generalized Inverse Gaussian distribution with parameters   *
 *               theta > 0 and omega > 0 using Ratio of Uniforms method      *
 *               without shift for theta <= 1 and omega <= 1 and shift at    *
 *               mode m otherwise.                                           *
 *                                                                           *
 * REFERENCE:  - J.S. Dagpunar (1989): An easily implemented generalized     *
 *               inverse Gaussian generator,                                 *
 *               Commun. Statist. Simul. 18(2), 703-710.                     *
 *                                                                           *
 * Implemented by R. Kremer, 1990                                            *
 * Revised     by F. Niederl, August 1992                                    *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define m       (GEN->gen_param[0])
#define linvmax (GEN->gen_param[1])
#define vminus  (GEN->gen_param[2])
#define vdiff   (GEN->gen_param[3])
#define b2      (GEN->gen_param[4])
#define hm12    (GEN->gen_param[5])
#define a       (GEN->gen_param[6])
#define d       (GEN->gen_param[7])
#define e       (GEN->gen_param[8])
#define c       (GEN->gen_param[9])
/*---------------------------------------------------------------------------*/
static const double drittel = 0.3333333333333333;                    /* 1/3  */
static const double pdrittel = 0.037037037037037;                    /* 1/27 */
/*---------------------------------------------------------------------------*/

inline static int
gig_gigru_init( struct unur_gen *gen )
{
  double r,s,t,p,q,xeta,fi,fak,yy1,yy2,max,invy1,invy2,vplus,hm1,xm,ym;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (theta <= 0) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return UNUR_ERR_GEN_CONDITION;
  }

  if (theta<=1. && omega<=1.) {
    /* NO SHIFT m */
    e = omega * omega;

    d = theta + 1.;
    ym = (-d + sqrt(d*d + e))/omega;

    d = theta - 1.;
    xm = (d + sqrt(d*d + e))/omega;

    d = 0.5 * d;
    e = -0.25 * omega;
    r = xm + 1./xm;
    s = xm*ym;
    /* a = vplus/uplus */
    a = exp(-0.5*theta*log(s) + 0.5*log(xm/ym) - e*(r - ym - 1.0/ym));
    /* c = 1/log{sqrt[hx(xm)]} */
    c = -d * log(xm) - e * r;
    /* vminus = 0 */
  }
  else {
    /* SHIFT BY m */
    hm1 = theta - 1.;
    hm12 = hm1 * 0.5;
    b2 = omega * 0.25;
    m = (hm1 + sqrt(hm1*hm1 + omega*omega))/omega;         /* Modus      */
    max = exp(hm12 * log(m) - b2 * (m + (1./m)));          /* sqrt[hx(m)] */
    linvmax = log(1.0/max);
	 
    /* Find the points x1,x2 (-->invy1,invy2) where
       the hat function touches the density f(x)    */
    r = (6.*m + 2.*theta*m - omega*m*m + omega)/(4.*m*m);
    s = (1. + theta - omega*m)/(2.*m*m);
    t = omega/(-4.*m*m);
    p = (3.*s - r*r) * drittel;
    q = (2.*r*r*r) * pdrittel - (r*s) * drittel + t;
    xeta = sqrt(-(p*p*p)*pdrittel);
    fi = acos(-q/(2.*xeta));
    fak = 2.*exp(log(xeta)*drittel);
    yy1 = fak * cos(fi*drittel) - r*drittel;
    yy2 = fak * cos(fi*drittel + 2.*drittel*M_PI) - r*drittel;
    invy1 = 1./yy1;
    invy2 = 1./yy2;
	 
    vplus = exp(linvmax + log(invy1) + hm12*log(invy1 + m)
		- b2*(invy1 + m + 1./(invy1 + m)));
    vminus = -exp(linvmax + log(-invy2) + hm12 * log(invy2 + m)
		  - b2*(invy2 + m + 1./(invy2 + m)));
    vdiff = vplus - vminus;
    /* uplus = 1 */
  }

  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of gig_gigru_init() */

double
_unur_stdgen_sample_gig_gigru( struct unur_gen *gen )
{
  double U,V,X,Z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  /* -X- generator code -X- */
  if (theta<=1. && omega<=1.) {
    /* NO SHIFT m */
    do {
      U = uniform();                                        /* U(0,1) */
      V = uniform();                                        /* U(0,1) */
      X = a*(V/U);
    }                                         /* Acceptance/Rejection */
    while (((log(U)) > (d*log(X) + e*(X + 1./X) + c)));
  }
  else {
    /* SHIFT BY m */
    do {
      do {
	U = uniform();                                      /* U(0,1) */
	V = vminus + uniform() * vdiff;                   /* U(v-,v+) */
	Z = V/U;
      } while (Z < (-m));
      X = Z + m;
    }                                         /* Acceptance/Rejection */
    while ((log(U) > (linvmax + hm12 * log(X) - b2 * (X + 1./X))));
  }

  /* -X- end of generator code -X- */
  
  return ((DISTR.n_params==2) ? X : eta * X );

} /* end of _unur_stdgen_sample_gig_gigru() */

/*---------------------------------------------------------------------------*/
#undef m
#undef linvmax
#undef vminus
#undef vdiff
#undef b2
#undef hm12
#undef a
#undef d
#undef e
#undef c
/*---------------------------------------------------------------------------*/
