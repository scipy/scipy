#ifndef _ARITHMETIC_GENERIC_GUARD
#define _ARITHMETIC_GENERIC_GUARD


#include <stddef.h>
#include <string.h>
#include <math.h>
#include <complex.h>

/* DBL_EPSILON / 2 == 2^(-53) */
#define DBL_HALF_EPSILON	(1.1102230246251565404236316680908e-16)
#define NPOLYNOMIALS		(4)


/* Basic algorithms not relying on availability of fma() */
static inline void feft_sum(double a, double b,
                            double * restrict x, double * restrict y)
{
    double z;

    *x = a + b;
    z = (*x) - a;
    *y = (a - ((*x) - z)) + (b - z);
}


static inline double fsum2(const double * restrict x, size_t n)
{
    size_t i;
    double p, s;

    p = x[0];
    s = 0.0;
    for ( i = 1; i < n; i++ )
    {
	double t, q;

	feft_sum(p, x[i], &t, &q);
	p = t;
	s += q;
    }
    return p + s;
}


#ifdef ELLINT_POLY_COMPLEX
static inline double_complex csum2(const double_complex * restrict x, size_t n)
{
    size_t i;
    double rr[n];
    double ii[n];
    for ( i = 0; i < n; i++ )
    {
	rr[i] = creal(x[i]);
	ii[i] = cimag(x[i]);
    }
    return MKCMPLX(fsum2(rr, n), fsum2(ii, n));
}


static inline void ceft_sum(double_complex x, double_complex y,
	                    double_complex * restrict s,
			    double * restrict er,
			    double * restrict ei)
{
    double sr, si;

    feft_sum(creal(x), creal(y), &sr, er);
    feft_sum(cimag(x), cimag(y), &si, ei);
    *s = MKCMPLX(sr, si);
}
#endif


#if ( __STDC_VERSION__ >= 199901L )	/* If the standard supports fma() */
static inline void feft_prod(double a, double b,
                             double * restrict x, double * restrict y)
{
    *x = a * b;
    *y = fma(a, b, -(*x));
}


#ifdef ELLINT_POLY_COMPLEX
static inline void ceft_prod(double_complex x, double_complex y,
			     double_complex * restrict p,
                             double buffer_r[3],
			     double buffer_i[3])
{
    double a, b, c, d;
    double z1, z2, z3, z4, z5, z6;
    /*  p                e    f   g
     * z5 | buffer_r : [h1, -h2, h5]
     * z6 | buffer_i : [h3,  h4, h6]
     */
    a = creal(x);
    b = cimag(x);
    c = creal(y);
    d = cimag(y);

    feft_prod(a, c, &z1, buffer_r);
    feft_prod(b, d, &z2, buffer_r + 1);
    feft_prod(a, d, &z3, buffer_i);
    feft_prod(b, c, &z4, buffer_i + 1);
    feft_sum(z1, -z2, &z5, buffer_r + 2);
    feft_sum(z3, z4, &z6, buffer_i + 2);
    *p = MKCMPLX(z5, z6);
    buffer_r[1] = -buffer_r[1];
}
#endif


#else	/* C standard not supporting fma() */


#define SPLIT_FACTOR	(134217729.0)
static inline void fsplit(double a, double * restrict x, double * restrict y)
{
    double c;

    c = SPLIT_FACTOR * a;
    *x = c - (c - a);
    *y = a - (*x);
}


static inline void feft_prod(double a, double b,
	                     double * restrict x, double * restrict y)
{
    double a1, a2, b1, b2;
    *x = a * b;
    fsplit(a, &a1, &a2);
    fsplit(b, &b1, &b2);
    *y = a2 * b2 - ((( (*x) - a1 * b1 ) - a2 * b1) - a1 * b2);
}


#ifdef ELLINT_POLY_COMPLEX
static inline void ceft_prod(double_complex x, double_complex y,
	                     double_complex * restrict p,
	                     double buffer_r[3],
	                     double buffer_i[3])
{
    double a, b, c, d;
    double a1, a2, b1, b2, c1, c2, d1, d2;
    double z1, z2, z3, z4, z5, z6;

    a = creal(x);
    b = cimag(x);
    c = creal(y);
    d = cimag(y);
    fsplit(a, &a1, &a2);
    fsplit(b, &b1, &b2);
    fsplit(c, &c1, &c2);
    fsplit(d, &d1, &d2);
    z1 = a * c;
    z2 = b * d;
    z3 = a * d;
    z4 = b * c;
    buffer_r[0] = a2 * c2 - (((z1 - a1 * c1) - a2 * c1) - a1 * c2);
    buffer_r[1] = -(b2 * d2 - (((z2 - b1 * d1) - b2 * d1) - b1 * d2));
    buffer_i[0] = a2 * d2 - (((z3 - a1 * d1) - a2 * d1) - a1 * d2);
    buffer_i[1] = b2 * c2 - (((z4 - b1 * c1) - b2 * c1) - b1 * c2);
    feft_sum(z1, -z2, &z5, buffer_r + 2);
    feft_sum(z3, z4, &z6, buffer_i + 2);
    *p = MKCMPLX(z5, z6);
}
#endif


#endif		/* End of selective compilation based on fma() availability */


/* Common algorithms utilizing ?eft_prod */
/* Ref: Gaillat, Ménissier-Morain: Compensated Horner scheme in complex
 *      floating point arithmetic.
 *      2008-07-07, Proc 8th Conf Real Numbers and Computers, pp. 133--146
 */
static inline double fdot2(const double * restrict x,
                           const double * restrict y, size_t n)
{
    size_t i;
    double p, s, h, r, t, q;

    feft_prod(x[0], y[0], &p, &s);
    for ( i = 1; i < n; i++ )
    {
	feft_prod(x[i], y[i], &h, &r);
	feft_sum(p, h, &t, &q);
	p = t;
	s += q + r;
    }
    return p + s;
}


#ifdef ELLINT_POLY_REAL
static inline double feft_horner(double x, const double * restrict a,
                                 size_t degree, 
				 double * restrict pp, double * restrict ps)
{
    ptrdiff_t i;
    double ptmp;
    double s = a[degree];

    for ( i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	feft_prod(s, x, &ptmp, pp + i);
	feft_sum(ptmp, a[i], &s, ps + i);
    }

    return s;
}


static inline double fhorner_sum(double x, const double * restrict p,
				 const double * restrict q, size_t degree)
{
    ptrdiff_t i;
    double r = p[degree] + q[degree];

    for ( i = (ptrdiff_t)degree - 1; i >= 0; i-- )
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
#endif


#ifdef ELLINT_POLY_COMPLEX
static inline double_complex cdot2(const double_complex * restrict x,
	                           const double_complex * restrict y,
                                   size_t twice_n)
{
    size_t i, j, n;
    double xx[twice_n];
    double yy[twice_n];
    double zz[twice_n];
    double rr, ii;

    n = twice_n / 2;
    for ( i = 0; i < n; i++ )
    {
	j = n + i;

	xx[i] = creal(x[i]);
	xx[j] = cimag(x[i]);

	yy[i] = zz[j] = creal(y[i]);
	yy[j] = -cimag(y[i]);

	zz[i] = cimag(y[i]);
    }
    rr = fdot2(xx, yy, twice_n);
    ii = fdot2(xx, zz, twice_n);
    return MKCMPLX(rr, ii);
}


static inline double_complex ceft_horner(double_complex x,
                                         const double * restrict a,
                                         size_t degree, 
					 double * * restrict wsr,
					 double * * restrict wsi)
{
    ptrdiff_t i;
    double_complex s;

    s = MKCR(a[degree]);
    for ( i = (ptrdiff_t)degree - 1; i >= 0; i-- )
    {
	double * restrict rowr;
	double * restrict rowi;
	double_complex p;

	rowr = (double * restrict)wsr + NPOLYNOMIALS * i;
	rowi = (double * restrict)wsi + NPOLYNOMIALS * i;
	ceft_prod(s, x, &p, rowr, rowi);
	ceft_sum(p, MKCR(a[i]), &s, rowr + 3, rowi + 3);
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
    size_t count;
} maskref_t;


static inline void mask_create(maskref_t * mask, bool * restrict mask_buf)
{
    mask->m = mask_buf;
    mask->count = 0;
}


static inline void mask_init(maskref_t * mask,
                             const double * restrict arr,
                             size_t len_arr)
{
    size_t i;
    memset(mask->m, 0, len_arr * (sizeof (bool)));
    for ( i = 0; i < len_arr; i++ )
    {
	if ( fpclassify(arr[i]) != FP_ZERO )
	{
	    (mask->count)++;
	    (mask->m)[i] = true;
	}
    }
}


static inline double next_power_two(double p)
{
    double q = p / DBL_HALF_EPSILON;
    double L = fabs((q + p) - q);
    if ( fpclassify(L) == FP_ZERO )
    {
	L = fabs(p);
    }

    return L;
}


static inline void fextract_scalar(double sigma, double * restrict p,
                                   double * restrict q)
{
    *q = (sigma + (*p)) - sigma;
    *p -= *q;
}


static inline void fextract_vector(double sigma, size_t lenp,
				   double * restrict p,
                                   maskref_t * mask,
				   double * restrict tau)
{
    double q;
    size_t i;

    q = 0.0;
    *tau = 0.0;
    for ( i = 0; i < lenp; i++ )
    {
	if ( (mask->m)[i] )
	{
	    fextract_scalar(sigma, p + i, &q);
	    if ( fpclassify(p[i]) == FP_ZERO )
	    {
		(mask->m)[i] = false;
		(mask->count)--;
	    }
	    *tau += q;
	}
    }
}


static inline double vmax(const double * restrict p, size_t lenp,
                          const maskref_t * mask,
			  double (* keyfcn)(double))
{
    size_t i;
    double r;

    r = 0.0;
    for ( i = 0; i < lenp; i++ )
    {
	if ( (mask->m)[i] )
	{
	    r = fmax(r, keyfcn(p[i]));
	}
    }
    return r;
}


static inline double facc_sum(double * restrict p, size_t lenp,
                              maskref_t * mask)
{
    double mu;
    double twopM, sigma, phi, factor;
    double t, tau, tau1, tau2, res;

    if ( (lenp == 0) || (mask->count == 0) )
    {
	return 0.0;
    }

    mu = vmax(p, lenp, mask, fabs);

    if ( fpclassify(mu) == FP_ZERO )
    {
	return 0.0;
    }

    twopM = next_power_two((double)(mask->count) + 2.0);
    sigma = twopM * next_power_two(mu);
    phi = DBL_HALF_EPSILON * twopM;
    factor = DBL_EPSILON * twopM * twopM;

    t = 0.0;
    while ( 1 )
    {
	size_t i;

	fextract_vector(sigma, lenp, p, mask, &tau);
	tau1 = t + tau;
	if ( (fabs(tau1) >= factor * sigma) || ( sigma <= DBL_MIN ) )
	{
	    /* Use Dekker's version of "eft_sum" (faster) with tau1 >= t */
	    tau2 = tau - (tau1 - t);

	    res = 0.0;
	    for ( i = 0; i < lenp; i++ )
	    {
		res += p[i];
	    }
	    res += tau2 + tau1;

	    return res;
	}
	t = tau1;
	if ( fpclassify(t) == FP_ZERO )
	{
	    res = facc_sum(p, lenp, mask);
	    return res;
	}
	sigma *= phi;
    }
}


static inline double_complex chorner_sum_acc(double_complex x, size_t degree,
                                             double * * restrict wsr,
					     double * * restrict wsi)
{
    /* Each row of ws is an array with 4 elements, and ws has (degree + 1)
     * rows.
     */
    double_complex v;
    maskref_t mask_r, mask_i;
    double * restrict rowr;
    double * restrict rowi;
    ptrdiff_t i;
    bool mask_buf_r[NPOLYNOMIALS];
    bool mask_buf_i[NPOLYNOMIALS];

    mask_create(&mask_r, mask_buf_r);
    mask_create(&mask_i, mask_buf_i);
    v = CZERO;

    for ( i = (ptrdiff_t)degree; i >= 0; i-- )
    {
	rowr = (double * restrict)wsr + NPOLYNOMIALS * i;
	rowi = (double * restrict)wsi + NPOLYNOMIALS * i;
	mask_init(&mask_r, rowr, NPOLYNOMIALS);
	mask_init(&mask_i, rowi, NPOLYNOMIALS);
	v = ADD(MULcc(v, x), MKCMPLX(facc_sum(rowr, NPOLYNOMIALS, &mask_r),
	                             facc_sum(rowi, NPOLYNOMIALS, &mask_i)));
    }
    return v;
}


static inline double_complex ccomp_horner(double_complex x,
                                          const double * restrict a,
					  size_t degree)
{
    double_complex res, c;
    double wsr[degree][NPOLYNOMIALS];
    double wsi[degree][NPOLYNOMIALS];

    res = ceft_horner(x, a, degree,
		      (double * * restrict)wsr, (double * * restrict)wsi);
    c = chorner_sum_acc(x, degree - 1,
                        (double * * restrict)wsr, (double * * restrict)wsi);
    return ADD(res, c);
}
#endif


#if 0
/* Fallback, for both real and complex types. The complex version of
 * EFT multiplication is too costly. */
static inline EllInt_Num_t ELLINT_POLY_FCN(naive_horner)(EllInt_Num_t z,
			       const double * restrict a, size_t degree)
{
    EllInt_Num_t r;
    ptrdiff_t i;

    r = MKCR(a[degree]);
    for ( i = (ptrdiff_t)degree - 1 ; i >= 0 ; i-- )
    {
	r = ADDcr(MULcc(r, z), a[i]);
    }

    return r;
}
#endif


#if ( __STDC_VERSION__ >= 199901L )	/* fma() supported by standard */
#define HORNER	(ELLINT_POLY_FCN(comp_horner))
#else					/* fma() not supported by standard */
/*
#define HORNER	(ELLINT_POLY_FCN(naive_horner))
*/
#define HORNER	(ELLINT_POLY_FCN(comp_horner))
#endif


#endif
