

/* constants */
static cdouble nc_1 = {1., 0.};
static cdouble nc_half = {0.5, 0.};
static cdouble nc_i = {0., 1.};
static cdouble nc_i2 = {0., 0.5};
/*
static cdouble nc_mi = {0., -1.};
static cdouble nc_pi2 = {M_PI/2., 0.};
*/

static void
nc_sum(cdouble *a, cdouble *b, cdouble *r)
{
	r->real = a->real + b->real;
	r->imag = a->imag + b->imag;
	return;
}

static void
nc_diff(cdouble *a, cdouble *b, cdouble *r)
{
	r->real = a->real - b->real;
	r->imag = a->imag - b->imag;
	return;
}

static void
nc_neg(cdouble *a, cdouble *r)
{
	r->real = -a->real;
	r->imag = -a->imag;
	return;
}

static void
nc_prod(cdouble *a, cdouble *b, cdouble *r)
{
	double ar=a->real, br=b->real, ai=a->imag, bi=b->imag;
	r->real = ar*br - ai*bi;
	r->imag = ar*bi + ai*br;
	return;
}

static void
nc_quot(cdouble *a, cdouble *b, cdouble *r)
{

	double ar=a->real, br=b->real, ai=a->imag, bi=b->imag;
	double d = br*br + bi*bi;
	r->real = (ar*br + ai*bi)/d;
	r->imag = (ai*br - ar*bi)/d;
	return;
}

static void
nc_floor_quot(cdouble *a, cdouble *b, cdouble *r)
{
	double ar=a->real, br=b->real, ai=a->imag, bi=b->imag;
	double d = br*br + bi*bi;
	r->real = floor((ar*br + ai*bi)/d);
	r->imag = 0;
	return;
}

static void
nc_sqrt(cdouble *x, cdouble *r)
{
	double s,d;
	if (x->real == 0. && x->imag == 0.)
		*r = *x;
	else {
		s = sqrt(0.5*(fabs(x->real) + hypot(x->real,x->imag)));
		d = 0.5*x->imag/s;
		if (x->real > 0.) {
			r->real = s;
			r->imag = d;
		}
		else if (x->imag >= 0.) {
			r->real = d;
			r->imag = s;
		}
		else {
			r->real = -d;
			r->imag = -s;
		}
	}
	return;
}

static void
nc_log(cdouble *x, cdouble *r)
{
	double l = hypot(x->real,x->imag);
	r->imag = atan2(x->imag, x->real);
	r->real = log(l);
	return;
}

static void
nc_log1p(cdouble *x, cdouble *r)
{
	double l = hypot(x->real + 1.0,x->imag);
	r->imag = atan2(x->imag, x->real + 1.0);
	r->real = log(l);
	return;
}

static void
nc_exp(cdouble *x, cdouble *r)
{
	double a = exp(x->real);
	r->real = a*cos(x->imag);
	r->imag = a*sin(x->imag);
	return;
}

static void
nc_expm1(cdouble *x, cdouble *r)
{
	double a = exp(x->real);
	r->real = a*cos(x->imag) - 1.0;
	r->imag = a*sin(x->imag);
	return;
}

static void
nc_pow(cdouble *a, cdouble *b, cdouble *r)
{
        intp n;
        double ar=a->real, br=b->real, ai=a->imag, bi=b->imag;

	if (br == 0. && bi == 0.) {
		r->real = 1.;
		r->imag = 0.;
                return;
	}
	if (ar == 0. && ai == 0.) {
		r->real = 0.;
		r->imag = 0.;
                return;
	}
        if (bi == 0 && (n=(intp)br) == br) {
            if (n > -100 && n < 100) {
                cdouble p, aa;
                intp mask = 1;
                if (n < 0) n = -n;
                aa = nc_1;
                p.real = ar; p.imag = ai;
                while (1) {
                        if (n & mask)
                            nc_prod(&aa,&p,&aa);
                        mask <<= 1;
                        if (n < mask || mask <= 0) break;
                        nc_prod(&p,&p,&p);
                }
                r->real = aa.real; r->imag = aa.imag;
                if (br < 0) nc_quot(&nc_1, r, r);
                return;
            }
        }
        /* complexobect.c uses an inline version of this formula
           investigate whether this had better performance or accuracy */
        nc_log(a, r);
        nc_prod(r, b, r);
        nc_exp(r, r);
        return;
}


static void
nc_prodi(cdouble *x, cdouble *r)
{
	r->real = -x->imag;
	r->imag = x->real;
	return;
}


static void
nc_acos(cdouble *x, cdouble *r)
{
	nc_prod(x,x,r);
	nc_diff(&nc_1, r, r);
	nc_sqrt(r, r);
	nc_prodi(r, r);
	nc_sum(x, r, r);
	nc_log(r, r);
	nc_prodi(r, r);
	nc_neg(r, r);
	return;
	/* return nc_neg(nc_prodi(nc_log(nc_sum(x,nc_prod(nc_i,
	   nc_sqrt(nc_diff(nc_1,nc_prod(x,x))))))));
	*/
}

static void
nc_acosh(cdouble *x, cdouble *r)
{
	nc_prod(x, x, r);
	nc_diff(&nc_1, r, r);
	nc_sqrt(r, r);
	nc_prodi(r, r);
	nc_sum(x, r, r);
	nc_log(r, r);
	return;
	/*
	  return nc_log(nc_sum(x,nc_prod(nc_i,
	  nc_sqrt(nc_diff(nc_1,nc_prod(x,x))))));
	*/
}

static void
nc_asin(cdouble *x, cdouble *r)
{
	cdouble a, *pa=&a;
	nc_prod(x, x, r);
	nc_diff(&nc_1, r, r);
	nc_sqrt(r, r);
	nc_prodi(x, pa);
	nc_sum(pa, r, r);
	nc_log(r, r);
	nc_prodi(r, r);
	nc_neg(r, r);
	return;
	/*
	      return nc_neg(nc_prodi(nc_log(nc_sum(nc_prod(nc_i,x),
	      nc_sqrt(nc_diff(nc_1,nc_prod(x,x)))))));
	*/
}


static void
nc_asinh(cdouble *x, cdouble *r)
{
	nc_prod(x, x, r);
	nc_sum(&nc_1, r, r);
	nc_sqrt(r, r);
	nc_diff(r, x, r);
	nc_log(r, r);
	nc_neg(r, r);
	return;
	/*
	  return nc_neg(nc_log(nc_diff(nc_sqrt(nc_sum(nc_1,nc_prod(x,x))),x)));
	*/
}

static void
nc_atan(cdouble *x, cdouble *r)
{
	cdouble a, *pa=&a;
	nc_diff(&nc_i, x, pa);
	nc_sum(&nc_i, x, r);
	nc_quot(r, pa, r);
	nc_log(r,r);
	nc_prod(&nc_i2, r, r);
	return;
	/*
	  return nc_prod(nc_i2,nc_log(nc_quot(nc_sum(nc_i,x),nc_diff(nc_i,x))));
	*/
}

static void
nc_atanh(cdouble *x, cdouble *r)
{
	cdouble a, *pa=&a;
	nc_diff(&nc_1, x, r);
	nc_sum(&nc_1, x, pa);
	nc_quot(pa, r, r);
	nc_log(r, r);
	nc_prod(&nc_half, r, r);
	return;
	/*
	  return nc_prod(nc_half,nc_log(nc_quot(nc_sum(nc_1,x),nc_diff(nc_1,x))));
	*/
}

static void
nc_cos(cdouble *x, cdouble *r)
{
	double xr=x->real, xi=x->imag;
	r->real = cos(xr)*cosh(xi);
	r->imag = -sin(xr)*sinh(xi);
	return;
}

static void
nc_cosh(cdouble *x, cdouble *r)
{
	double xr=x->real, xi=x->imag;
	r->real = cos(xi)*cosh(xr);
	r->imag = sin(xi)*sinh(xr);
	return;
}


#define M_LOG10_E 0.434294481903251827651128918916605082294397

static void
nc_log10(cdouble *x, cdouble *r)
{
	nc_log(x, r);
	r->real *= M_LOG10_E;
	r->imag *= M_LOG10_E;
	return;
}

static void
nc_sin(cdouble *x, cdouble *r)
{
	double xr=x->real, xi=x->imag;
	r->real = sin(xr)*cosh(xi);
	r->imag = cos(xr)*sinh(xi);
	return;
}

static void
nc_sinh(cdouble *x, cdouble *r)
{
	double xr=x->real, xi=x->imag;
	r->real = cos(xi)*sinh(xr);
	r->imag = sin(xi)*cosh(xr);
	return;
}

static void
nc_tan(cdouble *x, cdouble *r)
{
	double sr,cr,shi,chi;
	double rs,is,rc,ic;
	double d;
	double xr=x->real, xi=x->imag;
	sr = sin(xr);
	cr = cos(xr);
	shi = sinh(xi);
	chi = cosh(xi);
	rs = sr*chi;
	is = cr*shi;
	rc = cr*chi;
	ic = -sr*shi;
	d = rc*rc + ic*ic;
	r->real = (rs*rc+is*ic)/d;
	r->imag = (is*rc-rs*ic)/d;
	return;
}

static void
nc_tanh(cdouble *x, cdouble *r)
{
	double si,ci,shr,chr;
	double rs,is,rc,ic;
	double d;
	double xr=x->real, xi=x->imag;
	si = sin(xi);
	ci = cos(xi);
	shr = sinh(xr);
	chr = cosh(xr);
	rs = ci*shr;
	is = si*chr;
	rc = ci*chr;
	ic = si*shr;
	d = rc*rc + ic*ic;
	r->real = (rs*rc+is*ic)/d;
	r->imag = (is*rc-rs*ic)/d;
	return;
}
