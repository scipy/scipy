/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_poisson_gen.c                                              *
 *                                                                           *
 *   Special generators for Poisson distribution                             *
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
#include <methods/cstd.h>   /* for the definition of `UNUR_STDGEN_INVERSION' */
#include <methods/dstd_struct.h>
#include <methods/x_gen_source.h>
#include <distr/distr_source.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#include "unur_distributions.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int poisson_pdtabl_init( struct unur_gen *gen );
inline static int poisson_pdac_init( struct unur_gen *gen );
inline static int poisson_pprsc_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  (39)   /* maximal number of parameters for generator */
#define MAX_gen_iparams  (5)   /* maximal number of integer param. for gen.  */

/* parameters */
#define theta  (DISTR.params[0])    /* shape */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_poisson_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Poisson distribution            */
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
  case 1:  /* Tabulated Inversion combined with Acceptance Complement */
    if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */
    if (theta < 10.) {
      /* CASE B: Tabulated Inversion */
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdtabl );
      return poisson_pdtabl_init( gen );
    }
    else { /* theta >= 10. */
      /* CASE A: acceptance complement */
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdac );
      return poisson_pdac_init( gen );
    }

  case 2:  /* Tabulated Inversion combined with Patchwork Rejection */
    if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */
    if (theta < 10.) {
      /* CASE: Tabulated Inversion --> same as case 1 !! */
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdtabl );
      return poisson_pdtabl_init( gen );
    }
    else { /* theta >= 10. */
      /* CASE: Patchwork Rejection */
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pprsc );
      return poisson_pprsc_init( gen );
    }

    /** WinRand routine `pruec' (Ratio of Uniforms/Inversion) not implemented **/

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_poisson_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Poisson Distribution: Tabulated Inversion combined with                   *
 *                       Acceptance Complement                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Poisson distribution with  *
 *               parameter theta > 0.                                        *
 *               Tabulated Inversion for  theta < 10                         *
 *               Acceptance Complement for theta >= 10.                      *
 *                                                                           *
 * REFERENCE: - J.H. Ahrens, U. Dieter (1982): Computer generation of        * 
 *              Poisson deviates from modified normal distributions,         *
 *              ACM Trans. Math. Software 8, 163-179.                        *
 *                                                                           *
 * Implemented by R. Kremer, August 1990                                     *
 * Revised by E. Stadlober, April 1992                                       *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define m    (GEN->gen_iparam[0])
#define ll   (GEN->gen_iparam[1])

#define p0   (GEN->gen_param[0])
#define q    (GEN->gen_param[1])
#define p    (GEN->gen_param[2])
#define pp   ((GEN->gen_param)+3)  /* array of length 36 */
/*---------------------------------------------------------------------------*/

int
poisson_pdtabl_init( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
    GEN->n_gen_iparam = MAX_gen_iparams;
    GEN->gen_iparam = _unur_xmalloc(GEN->n_gen_iparam * sizeof(int));
  }

  /* -X- setup code -X- */
  m = (theta > 1.) ? ((int) theta) : 1;
  ll = 0;
  p0 = q = p = exp(-theta);
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of poisson_pdtabl_init() */

int
_unur_stdgen_sample_poisson_pdtabl( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  /* -X- generator code -X- */
  double U;
  int K,i;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  
  while (1) {
    U = uniform();              /* Step U. Uniform sample */
    K = 0;
    if (U <= p0) 
      return K;

    /* Step T. Table comparison */
    if (ll != 0) {               
      i = (U > 0.458) ? _unur_min(ll,m) : 1;
      for (K = i; K <=ll; K++)
	if (U <= pp[K])
	  return K;
      if (ll == 35) continue;
    }

    /* Step C. Creation of new prob. */
    for (K = ll +1; K <= 35; K++) {
      p *= theta / (double)K;
      q += p;
      pp[K] = q;
      if (U <= q) {
	ll = K;
	return K;
      }
    }
    ll = 35;
  }
  
  /* -X- end of generator code -X- */
  
} /* end of _unur_stdgen_sample_poisson_pdtabl() */

/*---------------------------------------------------------------------------*/
#undef m 
#undef ll
#undef p0
#undef q 
#undef p 
#undef pp
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#define l     (GEN->gen_iparam[0])

#define s     (GEN->gen_param[0])
#define d     (GEN->gen_param[1])
#define omega (GEN->gen_param[2])
#define b1    (GEN->gen_param[3])
#define b2    (GEN->gen_param[4])
#define c     (GEN->gen_param[5])
#define c0    (GEN->gen_param[6])
#define c1    (GEN->gen_param[7])
#define c2    (GEN->gen_param[8])
#define c3    (GEN->gen_param[9])

#define NORMAL  gen->gen_aux    /* pointer to normal variate generator        */
/*---------------------------------------------------------------------------*/

int
poisson_pdac_init( struct unur_gen *gen )
     /* Theta >= 10: acceptance complement */
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
    GEN->n_gen_iparam = MAX_gen_iparams;
    GEN->gen_iparam = _unur_xmalloc(GEN->n_gen_iparam * sizeof(int));
  }

  /* -X- setup code -X- */

  /* make a normal variate generator (use default special generator) */
  if (NORMAL==NULL) {
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    struct unur_par *par = unur_cstd_new( distr );
    NORMAL = (par) ? _unur_init(par) : NULL;
    _unur_check_NULL( NULL, NORMAL, UNUR_ERR_NULL );
    /* need same uniform random number generator as slash generator */
    NORMAL->urng = gen->urng;
    /* copy debugging flags */
    NORMAL->debug = gen->debug;
    /* we do not need the distribution object any more */
    _unur_distr_free( distr );
  }
  /* else we are in the re-init mode 
     --> there is no necessity to make the generator object again */

  s = sqrt(theta);
  d = 6. * theta * theta;
  l = (int)(theta - 1.1484);

  /* Step P. Preparations for steps Q and H */
  omega = 0.3989423 / s;
  b1 = 0.416666666667e-1 / theta;
  b2 = 0.3 * b1 * b1;
  c3 = 0.1428571 * b1 * b2;
  c2 = b2 - 15.0 * c3;
  c1 = b1 - 6.0 * b2 + 45.0 * c3;
  c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
  c = 0.1069 / theta;

  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of poisson_pdac_init() */

/*---------------------------------------------------------------------------*/
#define  a0  -0.5000000002
#define  a1   0.3333333343
#define  a2  -0.2499998565
#define  a3   0.1999997049
#define  a4  -0.1666848753
#define  a5   0.1428833286
#define  a6  -0.1241963125
#define  a7   0.1101687109
#define  a8  -0.1142650302
#define  a9   0.1055093006
/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_poisson_pdac( struct unur_gen *gen )
     /* Theta >= 10: acceptance complement */
{
  /* -X- generator code -X- */
  /* factorial for 0 <= k <= 9 */
  static const int fac[] = {1,1,2,6,24,120,720,5040,40320,362880};

  double t,g,theta_k;
  double gx,gy,px,py,x,xx,delta,v;
  int sign;

  double E, U;
  int K;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);

  /* Step N. Normal sample */
  t = _unur_sample_cont(NORMAL);
  g = theta + s * t;

  if (g >= 0.) {
    K = (int) g;
    /* Step I. Immediate acceptance */
    if (K >= l) 
      return K;
    /* Step S. Squeeze acceptance */
    U = uniform();
    theta_k = theta - K;
    if (d * U >= theta_k * theta_k * theta_k)
      return K;

    /* FUNCTION F */
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else {  /* k >= 10 */
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8 * delta*delta*delta * (1.-1./(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    /* end FUNCTION F */

    /* Step Q. Quotient acceptance */
    if (gy * (1.0 - U)  <= py * exp(px - gx))
      return K;
  }

  /* Step E. Double exponential sample */
  while (1) {
    do {
      E = - log(uniform());
      U = uniform();
      U = U + U - 1.;
      sign = (U < 0.) ? -1 : 1;
      t = 1.8 + E * sign;
    } while (t <= -0.6744);
    K = (int)(theta + s * t);
    theta_k = theta - K;

    /* FUNCTION F */
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else { /* k >= 10 */
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    /* end FUNCTION F */

    /* Step H. Hat acceptance */
    if (c * sign * U <= py * exp(px + E) - gy * exp(gx + E)) 
      return K;
  }

  /* -X- end of generator code -X- */
  
} /* end of _unur_stdgen_sample_poisson_pdac() */

/*---------------------------------------------------------------------------*/
#undef  a0
#undef  a1
#undef  a2
#undef  a3
#undef  a4
#undef  a5
#undef  a6
#undef  a7
#undef  a8
#undef  a9

#undef l
#undef s
#undef d
#undef omega
#undef b1
#undef b2
#undef c 
#undef c0
#undef c1
#undef c2
#undef c3

#undef NORMAL
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Poisson Distribution: Tabulated Inversion combined with                   *
 *                       Patchwork Rejection                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:  - samples a random number from the Poisson distribution with   *
 *              parameter theta > 0.                                         *
 *              Tabulated Inversion for  theta < 10                          *
 *              Patchwork Rejection for theta >= 10.                         *
 *                                                                           *
 * REFERENCE: - H. Zechner (1994): Efficient sampling from continuous and    *
 *              discrete unimodal distributions,                             *
 *              Pd.D. Thesis, 156 pp., Technical University Graz, Austria.   *
 *                                                                           *
 * Implemented by H. Zechner, January 1994                                   *
 * Revised by F. Niederl, July 1994                                          *
 *****************************************************************************
 *                                                                           *
 * Patchwork Rejection:                                                      *
 * The area below the histogram function f(x) is rearranged in its body by   *
 * certain point reflections. Within a large center interval variates are    *
 * sampled efficiently by rejection from uniform hats. Rectangular immediate *
 * acceptance regions speed up the generation. The remaining tails are       *
 * covered by exponential functions.                                         *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

inline static double f(int k, double l_nu, double c_pm)
{
  return  exp(k * l_nu - _unur_SF_ln_factorial(k) - c_pm);
}

/*---------------------------------------------------------------------------*/
#define m       (GEN->gen_iparam[0])
#define k2      (GEN->gen_iparam[1])
#define k4      (GEN->gen_iparam[2])
#define k1      (GEN->gen_iparam[3])
#define k5      (GEN->gen_iparam[4])

#define dl      (GEN->gen_param[0])
#define dr      (GEN->gen_param[1])
#define r1      (GEN->gen_param[2])
#define r2      (GEN->gen_param[3])
#define r4      (GEN->gen_param[4])
#define r5      (GEN->gen_param[5])
#define ll      (GEN->gen_param[6])
#define lr      (GEN->gen_param[7])
#define l_theta (GEN->gen_param[8])
#define c_pm    (GEN->gen_param[9])
#define f2      (GEN->gen_param[10])
#define f4      (GEN->gen_param[11])
#define f1      (GEN->gen_param[12])
#define f5      (GEN->gen_param[13])
#define p1      (GEN->gen_param[14])
#define p2      (GEN->gen_param[15])
#define p3      (GEN->gen_param[16])
#define p4      (GEN->gen_param[17])
#define p5      (GEN->gen_param[18])
#define p6      (GEN->gen_param[19])
/*---------------------------------------------------------------------------*/

int
poisson_pprsc_init( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  double Ds;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
    GEN->n_gen_iparam = MAX_gen_iparams;
    GEN->gen_iparam = _unur_xmalloc(GEN->n_gen_iparam * sizeof(int));
  }

  /* -X- setup code -X- */

  /* approximate deviation of reflection points k2, k4 from theta - 1/2      */
  Ds = sqrt(theta + 0.25);

  /* mode m, reflection points k2 and k4, and points k1 and k5, which        */
  /* delimit the centre region of h(x)                                       */
  m  = (int) theta;
  k2 = (int) (theta + 0.5 - Ds);
  k4 = (int) (theta - 0.5 + Ds);
  k1 = k2 + k2 - m + 1;
  k5 = k4 + k4 - m;

  /* range width of the critical left and right centre region                */
  dl = (double) (k2 - k1);
  dr = (double) (k5 - k4);

  /* recurrence constants r(k) = p(k)/p(k-1) at k = k1, k2, k4+1, k5+1       */
  r1 = theta / (double) k1;
  r2 = theta / (double) k2;
  r4 = theta / (double)(k4 + 1);
  r5 = theta / (double)(k5 + 1);

  /* reciprocal values of the scale parameters of expon. tail envelopes      */
  ll =  log(r1);                                   /* expon. tail left */
  lr = -log(r5);                                   /* expon. tail right*/

  /* Poisson constants, necessary for computing function values f(k)         */
  l_theta = log(theta);
  c_pm = m * l_theta - _unur_SF_ln_factorial(m);

  /* function values f(k) = p(k)/p(m) at k = k2, k4, k1, k5                  */
  f2 = f(k2, l_theta, c_pm);
  f4 = f(k4, l_theta, c_pm);
  f1 = f(k1, l_theta, c_pm);
  f5 = f(k5, l_theta, c_pm);
  
  /* area of the two centre and the two exponential tail regions             */
  /* area of the two immediate acceptance regions between k2, k4             */
  p1 = f2 * (dl + 1.);                            /* immed. left      */
  p2 = f2 * dl        + p1;                       /* centre left      */
  p3 = f4 * (dr + 1.) + p2;                       /* immed. right     */
  p4 = f4 * dr        + p3;                       /* centre right     */
  p5 = f1 / ll        + p4;                       /* expon. tail left */
  p6 = f5 / lr        + p5;                       /* expon. tail right*/
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of poisson_pprsc_init() */

int
_unur_stdgen_sample_poisson_pprsc( struct unur_gen *gen )
     /* theta >= 10: Patchwork Rejection */
{
  /* -X- generator code -X- */
  int    Dk, X, Y;
  double U, V, W;

  while (1) {

    /* generate uniform number U -- U(0, p6)                                 */
    U = uniform() * p6;

    /* case distinction corresponding to U                                   */
    if (U < p2) {
      /* centre left      */

      /* immediate acceptance region R2 = [k2, m) *[0, f2),  X = k2, ... m-1 */
      V = U - p1;
      if (V < 0.)
	return (k2 + (int)(U/f2));

      /* immediate acceptance region R1 = [k1, k2)*[0, f1),  X = k1, ... k2-1 */
      W = V / dl;
      if (W < f1)
	return (k1 + (int)(V/f1));

      /* computation of candidate X < k2, and its counterpart Y > k2         */
      /* either squeeze-acceptance of X or acceptance-rejection of Y         */
      Dk = (int)(dl * uniform()) + 1;
      if (W <= f2 - Dk * (f2 - f2/r2))
	/* quick accept of  */
	return (k2 - Dk);                            /* X = k2 - Dk      */
      if ((V = f2 + f2 - W) < 1.) {
	/* quick reject of Y*/
	Y = k2 + Dk;
	if (V <= f2 + Dk * (1. - f2)/(dl + 1.))
	  /* quick accept of  */
	  return Y;                                  /* Y = k2 + Dk      */
	if (V <= f(Y, l_theta, c_pm))  
	  /* final accept of Y*/
	  return Y;
      }
      X = k2 - Dk;
    }

    else if (U < p4) {
      /* centre right     */

      /*  immediate acceptance region R3 = [m, k4+1)*[0, f4), X = m, ... k4  */
      (V = U - p3);
      if (V < 0.)
	return (k4 - (int)((U - p2)/f4));

      /* immediate acceptance region R4 = [k4+1, k5+1)*[0, f5)               */
      W = V / dr;
      if (W < f5 )
	return (k5 - (int)(V/f5));

      /* computation of candidate X > k4, and its counterpart Y < k4         */
      /* either squeeze-acceptance of X or acceptance-rejection of Y         */
      Dk = (int)(dr * uniform()) + 1;
      if (W <= f4 - Dk * (f4 - f4*r4))
	/* quick accept of  */
	return(k4 + Dk);                             /* X = k4 + Dk      */
      if ((V = f4 + f4 - W) < 1.0) {
	/* quick reject of Y*/
	Y = k4 - Dk;
	if (V <= f4 + Dk * (1.0 - f4)/ dr)
	  /* quick accept of  */
	  return Y;                                 /* Y = k4 - Dk      */
	if (V <= f(Y, l_theta, c_pm))
	  return Y;       /* final accept of Y*/
      }
      X = k4 + Dk;
    }

    else {
      W = uniform();
      if (U < p5) {
	/* expon. tail left */
	Dk = (int)(1. - log(W)/ll);
	X = k1 - Dk;
	if (X < 0)
	  continue;           /* 0 <= X <= k1 - 1 */
	W *= (U - p4) * ll;                          /* W -- U(0, h(x))      */
	if (W <= f1 - Dk * (f1 - f1/r1))
	  return X; /* quick accept of X*/
      }
      else {
	/* expon. tail right*/
	Dk = (int)(1. - log(W)/lr);
	X  = k5 + Dk;                                /* X >= k5 + 1          */
	W *= (U - p5) * lr;                          /* W -- U(0, h(x))      */
	if (W <= f5 - Dk * (f5 - f5*r5))
	  return X; /* quick accept of X*/
      }
    }

    /* acceptance-rejection test of candidate X from the original area       */
    /* test, whether  W <= f(k),    with  W = U*h(x)  and  U -- U(0, 1)      */
    /* log f(X) = (X - m)*log(theta) - log X! + log m!                       */
    if (log(W) <= X * l_theta - _unur_SF_ln_factorial(X) - c_pm)
      return X;

  }

  /* -X- end of generator code -X- */
  
} /* end of _unur_stdgen_sample_poisson_pprsc() */

/*---------------------------------------------------------------------------*/
#undef m 
#undef k2
#undef k4
#undef k1
#undef k5

#undef dl
#undef dr
#undef r1
#undef r2
#undef r4
#undef r5
#undef ll
#undef lr
#undef l_theta
#undef c_pm
#undef f2
#undef f4
#undef f1
#undef f5
#undef p1
#undef p2
#undef p3
#undef p4
#undef p5
#undef p6
/*---------------------------------------------------------------------------*/
