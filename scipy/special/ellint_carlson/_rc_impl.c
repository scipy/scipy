#ifndef _ELLINT_RC_GENERIC_GUARD
#define _ELLINT_RC_GENERIC_GUARD

#include <stddef.h>
#include <math.h>
#include "ellint_common.h"
#include "ellint_poly.h"
#include "ellint_argcheck.h"


#if ( __STDC_VERSION__ >= 199901L )	/* fma() supported by standard */
#define NPOLYNOMIALS	(4)
static inline EllInt_Num_t ELLINT_POLY_FCN(comp_horner)(EllInt_Num_t x,
	const double * restrict a, size_t n);
#define HORNER	(ELLINT_POLY_FCN(comp_horner))
#else					/* fma() not supported by standard */
static inline EllInt_Num_t ELLINT_POLY_FCN(naive_horner)(EllInt_Num_t x,
	const double * restrict a, size_t n);
#define HORNER	(ELLINT_POLY_FCN(naive_horner))
#endif


static const double ELLINT_RC_C[6] = {
    0.3000000000000000000000000,	/* 3 / 10 */
    0.1428571428571428571428571,	/* 1 / 7 */
    0.3750000000000000000000000,	/* 3 / 8 */
    0.4090909090909090909090909,	/* 9 / 22 */
    0.7644230769230769230769231,	/* 159 / 208 */
    1.1250000000000000000000000		/* 9 / 8 */
};


int
ELLINT_POLY_FCN(ellint_RC) (EllInt_Num_t x, EllInt_Num_t y,
	                    double rerr, EllInt_Num_t * restrict res)
{
    int status;

    if ( ELLINT_BAD_RERR(rerr, 2.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( ( too_small(cimag(y)) ) && ELLINT_RNEG(y) )
    {
	/* Cauchy principal value with negative real y */
	status = ELLINT_POLY_FCN(ellint_RC)(x - y, -y, rerr, res);

	if ( status == ELLINT_STATUS_SUCCESS )
	{
	    (*res) *= sqrt(x / (x - y));
	}
	return status;
    } else if ( ELLINT_RNEG(x) || ( too_small(fabs(y)) ) ) {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    unsigned int m;
    EllInt_Num_t A0;
    EllInt_Num_t xm, ym, sm, Am;
    double fterm, Q;

    A0 = (x + 2.0 * y) / 3.0;
    Q = fabs(A0 - x) / ELLINT_OCRT(3.0 * rerr);

    fterm = Q;
    Am = A0;
    xm = x;
    ym = y;
    m = 0;
    sm = y - A0;

    while ( 1 )
    {
	EllInt_Num_t xrm, yrm, lam;

        if ( fterm < fabs(Am) )
	{
            break;
	}

	if ( m > ELLINT_MAXITER )
	{
	    ELLINT_FAIL_WITH(ELLINT_STATUS_NITER);
	}

        xrm = sqrt(xm);
        yrm = sqrt(ym);
	lam = 2.0 * xrm * yrm + ym;

        Am = (Am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
	sm *= 0.25;
	fterm *= 0.25;

	m += 1;
    }
    sm /= Am;
    (*res) = (1.0 + sm * sm * HORNER(sm, ELLINT_RC_C, 5)) / sqrt(Am);

    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#if ( __STDC_VERSION__ >= 199901L )	/* If the standard supports fma() */
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
/* Ref: Gaillat, Langlois, Louvet: Compensated Horner Scheme.
 *      2006-07-24, Research Report No. RR2005-04, Univ. Perpignan Via Domitia
 *      --: Improving the Compensated Horner Scheme with a Fused Multiply and
 *      Add.
 *      2006-04-23, SAC'06 (Proc 2006 ACM Symposium Applied Computing),
 *      pp.1323--1327 <https://doi.org/10.1145/1141277.1141585> */
static inline void feft_sum(double a, double b,
                            double * restrict x, double * restrict y)
{
    double z;

    *x = a + b;
    z = (*x) - a;
    *y = (a - ((*x) - z)) + (b - z);
}


static inline void feft_prod(double a, double b,
                             double * restrict x, double * restrict y)
{
    *x = a * b;
    *y = fma(a, b, -(*x));
}


static inline double feft_horner(double x, const double * restrict a,
                                 size_t degree, 
				 double * restrict pp, double * restrict ps)
{
    double s = a[degree];
    double ptmp;

    for ( ptrdiff_t i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	feft_prod(s, x, &ptmp, pp + i);
	feft_sum(ptmp, a[i], &s, ps + i);
    }

    return s;
}


static inline double fhorner_sum(double x, const double * restrict p,
				 const double * restrict q, size_t degree)
{
    double r = p[degree] + q[degree];

    for ( ptrdiff_t i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	r = r * x + (p[i] + q[i]);
    }

    return r;
}


static inline double fcomp_horner(double x, const double * restrict a,
                                  size_t degree)
{
    double c, h;
    double pp[degree];
    double ps[degree];

    h = feft_horner(x, a, degree, pp, ps);
    c = fhorner_sum(x, pp, ps, degree - 1);

    return h + c;
}

/* Ref: Gaillat, Ménissier-Morain: Compensated Horner scheme in complex
 *      floating point arithmetic.
 *      2008-07-07, Proc 8th Conf Real Numbers and Computers, pp. 133--146
 */
#define CEQZERO(z)	( (fpclassify(creal( (z) )) == FP_ZERO) && \
                          (fpclassify(cimag( (z) )) == FP_ZERO) )
static inline void ceft_sum(double complex x, double complex y,
	                    double complex * restrict s,
			    double complex * restrict e)
{
    double s1, s2, e1, e2;

    feft_sum(creal(x), creal(y), &s1, &e1);
    feft_sum(cimag(x), cimag(y), &s2, &e2);
    *s = s1 + s2 * I;
    *e = e1 + e2 * I;
}

static inline void ceft_prod(double complex x, double complex y,
	                     double complex * restrict p,
	                     double complex * restrict e,
	                     double complex * restrict f,
	                     double complex * restrict g)
{
    double z1, h1, z2, h2, z3, h3, z4, h4, z5, h5, z6, h6;
    double a = creal(x);
    double b = cimag(x);
    double c = creal(y);
    double d = cimag(y);

    feft_prod(a, c, &z1, &h1);
    feft_prod(b, d, &z2, &h2);
    feft_prod(a, d, &z3, &h3);
    feft_prod(b, c, &z4, &h4);
    feft_sum(z1, -z2, &z5, &h5);
    feft_sum(z3, z4, &z6, &h6);
    *p = z5 + z6 * I;
    *e = h1 + h3 * I;
    *f = -h2 + h4 * I;
    *g = h5 + h6 * I;
}

static inline double complex ceft_horner(double complex x,
                                         const double * restrict a,
                                         size_t degree, 
					 double complex * * restrict ws)
{
    double complex s;

    s = a[degree] + 0.0 * I;
    for ( ptrdiff_t i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	double complex * row;
	double complex p;

	row = (double complex *)ws + NPOLYNOMIALS * i;
	ceft_prod(s, x, &p, row, row + 1, row + 2);
	ceft_sum(p, a[i] + 0.0 * I, &s, row + 3);
    }
    return s;
}

/* Ref: Rump, Ogita, Oishi: Accurate Floating-Point Summation.
 *      2005-11-13, Technical Report 05.1, Faculty Information Communication
 *      Sci, Hamburg Univ Technol
 */
typedef struct maskref
{
    bool * m;
    size_t len;
    size_t count;
} maskref_t;

static inline void mask_create(maskref_t * mask, bool * restrict mask_buf,
                               size_t n)
{
    mask->m = mask_buf;
    mask->count = 0;

    mask->len = n;
}

static inline void mask_init(maskref_t * mask,
                             const double complex * restrict arr,
                             size_t len_arr)
{
    for (size_t i = 0; i < len_arr; i++ )
    {
	if ( CEQZERO(arr[i]) )
	{
	    (mask->m)[i] = false;
	} else {
	    (mask->count)++;
	    (mask->m)[i] = true;
	}
    }
}

static inline double next_power_two(double p)
{
    double q = p / DBL_EPSILON;
    double L = fabs((q + p) - q);
    if ( fpclassify(L) == FP_ZERO )
    {
	L = fabs(p);
    }

    return L;
}

static inline void cextract_scalar(double sigma, double complex * restrict p,
                                   double complex * restrict q)
{
    *q = sigma + (*p);
    *q -= sigma;
    *p -= *q;
}

static inline void cextract_vector(double sigma,
				   size_t lenp,
				   double complex * restrict p,
                                   maskref_t * mask,
				   double complex * restrict tau)
{
    *tau = 0.0 + 0.0 * I;
    double complex q = 0.0 + 0.0 * I;
    for ( size_t i = 0; i < lenp; i++ )
    {
	if ( (mask->m)[i] )
	{
	    cextract_scalar(sigma, p + i, &q);
	    if ( CEQZERO(p[i]) )
	    {
		(mask->m)[i] = false;
		(mask->count)--;

	    }
	    *tau += q;
	}
    }
}

static inline double vmax(const double complex * restrict p, size_t lenp,
                          const maskref_t * mask,
			  double (* keyfcn)(double complex))
{
    double r = 0.0;
    for ( size_t i = 0; i < lenp; i++ )
    {
	if ( (mask->m)[i] )
	{
	    r = fmax(r, keyfcn(p[i]));
	}
    }
    return r;
}

static inline double complex cacc_sum(double complex * restrict p, size_t lenp,
                                      maskref_t * mask)
{
    if ( (lenp == 0) || (mask->count == 0) )
    {
	return 0.0;
    }

    double mu;
    mu = vmax(p, lenp, mask, cabs);
    if ( fpclassify(mu) == FP_ZERO )
    {
	return 0.0;
    }

    double twopM, sigma, phi, factor;
    double complex t, tau, tau1, tau2, res;

    twopM = next_power_two((double)(mask->count) + 2.0);
    sigma = twopM * next_power_two(mu);
    phi = DBL_EPSILON * twopM * 0.5;
    factor = DBL_EPSILON * twopM * twopM;

    t = 0.0 + 0.0 * I;
    while ( 1 )
    {
	cextract_vector(sigma, lenp, p, mask, &tau);
	tau1 = t + tau;
	if ( (cabs(tau1) >= factor * sigma) || ( sigma <= DBL_MIN ) )
	{
	    ceft_sum(t, tau, &tau1, &tau2);

	    res = 0.0 + 0.0 * I;
	    for ( size_t i = 0; i < lenp; i++ )
	    {
		res += p[i];
	    }
	    res += tau2 + tau1;

	    return res;
	}
	t = tau1;
	if ( CEQZERO(t) )
	{
	    res = cacc_sum(p, lenp, mask);
	    return res;
	}
	sigma *= phi;
    }
}

static inline double complex chorner_sum_acc(double complex x, size_t degree,
                                             double complex * * restrict ws)
{
    /* Each row of ws is an array with 4 elements, and ws has (degree + 1)
     * rows.
     */
    double complex v;
    maskref_t mask;
    double complex * row;
    bool mask_buf[NPOLYNOMIALS];

    mask_create(&mask, mask_buf, NPOLYNOMIALS);
    row = (double complex *)ws + NPOLYNOMIALS * degree;
    mask_init(&mask, row, NPOLYNOMIALS);
    v = cacc_sum(row, NPOLYNOMIALS, &mask);
    for ( ptrdiff_t i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	row = (double complex *)ws + NPOLYNOMIALS * i;
	mask_init(&mask, row, NPOLYNOMIALS);
	v = v * x + cacc_sum(row, NPOLYNOMIALS, &mask);
    }
    return v;
}

static inline double complex ccomp_horner(double complex x,
                                          const double * restrict a,
					  size_t degree)
{
    double complex ws[degree][NPOLYNOMIALS];
    double complex res, c;

    res = ceft_horner(x, a, degree, (double complex * *)ws);
    c = chorner_sum_acc(x, degree - 1, (double complex * *)ws);
    return res + c;
}

#else	/* C standard not supporting fma() */

/* Fallback, for both real and complex types. The complex version of
 * EFT multiplication is too costly. */
static inline EllInt_Num_t ELLINT_POLY_FCN(naive_horner)(EllInt_Num_t z,
			       const double * restrict a, size_t degree)
{
    EllInt_Num_t r;
    ptrdiff_t i;

    r = a[degree];
    for ( i = (ptrdiff_t)degree - 1 ; i >= 0 ; i-- )
    {
	r = r * z + a[i];
    }

    return r;
}
#endif


#endif  /* _ELLINT_RC_GENERIC_GUARD */
