/* SIGTOOLS module by Travis Oliphant

Copyright 2005 Travis Oliphant
Permission to use, copy, modify, and distribute this software without fee
is granted under the SciPy License.
*/
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_signal_ARRAY_API
#include "numpy/ndarrayobject.h"
#include "npy_2_compat.h"

#include "_sigtools.hh"
#include <stdlib.h>

#define PYERR(message) {PyErr_SetString(PyExc_ValueError, message); goto fail;}


/********************************************************
 *
 *  Code taken from remez.c by Erik Kvaleberg which was
 *    converted from an original FORTRAN by
 *
 * AUTHORS: JAMES H. MCCLELLAN
 *
 *         DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
 *         MASSACHUSETTS INSTITUTE OF TECHNOLOGY
 *         CAMBRIDGE, MASS. 02139
 *
 *         THOMAS W. PARKS
 *         DEPARTMENT OF ELECTRICAL ENGINEERING
 *         RICE UNIVERSITY
 *         HOUSTON, TEXAS 77001
 *
 *         LAWRENCE R. RABINER
 *         BELL LABORATORIES
 *         MURRAY HILL, NEW JERSEY 07974
 *
 *
 *  Adaptation to C by
 *      egil kvaleberg
 *      husebybakken 14a
 *      0379 oslo, norway
 *  Email:
 *      egil@kvaleberg.no
 *  Web:
 *      http://www.kvaleberg.com/
 *
 *
 *********************************************************/


#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define GOBACK goto
#define DOloop(a,from,to) for ( (a) = (from); (a) <= (to); ++(a))
#define PI    3.14159265358979323846
#define TWOPI (PI+PI)

/*
 *-----------------------------------------------------------------------
 * FUNCTION: lagrange_interp (d)
 *  FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
 *  COEFFICIENTS FOR USE IN THE FUNCTION gee.
 *-----------------------------------------------------------------------
 */
static double lagrange_interp(int k, int n, int m, double *x)
{
    int j, l;
    double q, retval;

    retval = 1.0;
    q = x[k];
    DOloop(l,1,m) {
	for (j = l; j <= n; j += m) {
	    if (j != k)
		retval *= 2.0 * (q - x[j]);
	}
    }
    return 1.0 / retval;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: freq_eval (gee)
 *  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
 *  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
 *-----------------------------------------------------------------------
 */
static double freq_eval(int k, int n, double *grid, double *x, double *y, double *ad)
{
    int j;
    double p,c,d,xf;

    d = 0.0;
    p = 0.0;
    xf = cos(TWOPI * grid[k]);

    DOloop(j,1,n) {
	c = ad[j] / (xf - x[j]);
	d += c;
	p += c * y[j];
    }

    return p/d;
}


/*
 *-----------------------------------------------------------------------
 * SUBROUTINE: remez
 *  THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
 *  FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
 *  FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
 *  ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
 *  DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
 *  GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
 *  EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
 *  ERROR BY DETERMINING THE BSMINEST LOCATION OF THE EXTREMAL
 *  FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
 *  THE COEFFICIENTS OF THE BEST APPROXIMATION.
 *-----------------------------------------------------------------------
 */
static int remez(double *dev, double des[], double grid[], double edge[],
	   double wt[], int ngrid, int nbands, int iext[], double alpha[],
	   int nfcns, int itrmax, double *work, int dimsize, int *niter_out)
		/* dev, iext, alpha                         are output types */
		/* des, grid, edge, wt, ngrid, nbands, nfcns are input types */
{
    int k, k1, kkk, kn, knz, klow, kup, nz, nzz, nm1;
    int cn;
    int j, jchnge, jet, jm1, jp1;
    int l, luck=0, nu, nut, nut1=0, niter;

    double ynz=0.0, comp=0.0, devl, gtemp, fsh, y1=0.0, err, dtemp, delf, dnum, dden;
    double aa=0.0, bb=0.0, ft, xe, xt;

    static double *a, *p, *q;
    static double *ad, *x, *y;

    a = work; p = a + dimsize+1; q = p + dimsize+1;
    ad = q + dimsize+1; x = ad + dimsize+1; y = x + dimsize+1;
    devl = -1.0;
    nz  = nfcns+1;
    nzz = nfcns+2;
    niter = 0;

    do {
    L100:
	iext[nzz] = ngrid + 1;
	++niter;

	if (niter > itrmax) break;

	/* printf("ITERATION %2d: ",niter); */

	DOloop(j,1,nz) {
	    x[j] = cos(grid[iext[j]]*TWOPI);
	}
	jet = (nfcns-1) / 15 + 1;

	DOloop(j,1,nz) {
	    ad[j] = lagrange_interp(j,nz,jet,x);
	}

	dnum = 0.0;
	dden = 0.0;
	k = 1;

	DOloop(j,1,nz) {
	    l = iext[j];
	    dnum += ad[j] * des[l];
	    dden += (double)k * ad[j] / wt[l];
	    k = -k;
	}
	*dev = dnum / dden;

	/* printf("DEVIATION = %lg\n",*dev); */

	nu = 1;
	if ( (*dev) > 0.0 ) nu = -1;
	(*dev) = -(double)nu * (*dev);
	k = nu;
	DOloop(j,1,nz) {
	    l = iext[j];
	    y[j] = des[l] + (double)k * (*dev) / wt[l];
	    k = -k;
	}
	if ( (*dev) <= devl ) {
	    /* finished */
	    *niter_out = niter;
	    return -1;
	}
	devl = (*dev);
	jchnge = 0;
	k1 = iext[1];
	knz = iext[nz];
	klow = 0;
	nut = -nu;
	j = 1;

    /*
     * SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
     */

    L200:
	if (j == nzz) ynz = comp;
	if (j >= nzz) goto L300;
	kup = iext[j+1];
	l = iext[j]+1;
	nut = -nut;
	if (j == 2) y1 = comp;
	comp = (*dev);
	if (l >= kup) goto L220;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L220;
	comp = (double)nut * err;
    L210:
	if (++l >= kup) goto L215;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L215;
	comp = (double)nut * err;
	GOBACK L210;

    L215:
	iext[j++] = l - 1;
	klow = l - 1;
	++jchnge;
	GOBACK L200;

    L220:
	--l;
    L225:
	if (--l <= klow) goto L250;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) > 0.0) goto L230;
	if (jchnge <= 0) goto L225;
	goto L260;

    L230:
	comp = (double)nut * err;
    L235:
	if (--l <= klow) goto L240;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L240;
	comp = (double)nut * err;
	GOBACK L235;
    L240:
	klow = iext[j];
	iext[j] = l+1;
	++j;
	++jchnge;
	GOBACK L200;

    L250:
	l = iext[j]+1;
	if (jchnge > 0) GOBACK L215;

    L255:
	if (++l >= kup) goto L260;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L255;
	comp = (double)nut * err;

	GOBACK L210;
    L260:
	klow = iext[j++];
	GOBACK L200;

    L300:
	if (j > nzz) goto L320;
	if (k1 > iext[1] ) k1 = iext[1];
	if (knz < iext[nz]) knz = iext[nz];
	nut1 = nut;
	nut = -nu;
	l = 0;
	kup = k1;
	comp = ynz*(1.00001);
	luck = 1;
    L310:
	if (++l >= kup) goto L315;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L310;
	comp = (double) nut * err;
	j = nzz;
	GOBACK L210;

    L315:
	luck = 6;
	goto L325;

    L320:
	if (luck > 9) goto L350;
	if (comp > y1) y1 = comp;
	k1 = iext[nzz];
    L325:
	l = ngrid+1;
	klow = knz;
	nut = -nut1;
	comp = y1*(1.00001);
    L330:
	if (--l <= klow) goto L340;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L330;
	j = nzz;
	comp = (double) nut * err;
	luck = luck + 10;
	GOBACK L235;
    L340:
	if (luck == 6) goto L370;
	DOloop(j,1,nfcns) {
	    iext[nzz-j] = iext[nz-j];
	}
	iext[1] = k1;
	GOBACK L100;
    L350:
	kn = iext[nzz];
	DOloop(j,1,nfcns) iext[j] = iext[j+1];
	iext[nz] = kn;

	GOBACK L100;
    L370:
	;
    } while (jchnge > 0);

/*
 *    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
 *    USING THE INVERSE DISCRETE FOURIER TRANSFORM
 */
    nm1 = nfcns - 1;
    fsh = 1.0e-06;
    gtemp = grid[1];
    x[nzz] = -2.0;
    cn  = 2*nfcns - 1;
    delf = 1.0/cn;
    l = 1;
    kkk = 0;

    if (edge[1] == 0.0 && edge[2*nbands] == 0.5) kkk = 1;

    if (nfcns <= 3) kkk = 1;
    if (kkk !=     1) {
	dtemp = cos(TWOPI*grid[1]);
	dnum  = cos(TWOPI*grid[ngrid]);
	aa    = 2.0/(dtemp-dnum);
	bb    = -(dtemp+dnum)/(dtemp-dnum);
    }

    DOloop(j,1,nfcns) {
	ft = (j - 1) * delf;
	xt = cos(TWOPI*ft);
	if (kkk != 1) {
	    xt = (xt-bb)/aa;
#if 0
	    /*XX* check up !! */
	    xt1 = sqrt(1.0-xt*xt);
	    ft = atan2(xt1,xt)/TWOPI;
#else
	    ft = acos(xt)/TWOPI;
#endif
	}
L410:
	xe = x[l];
	if (xt > xe) goto L420;
	if ((xe-xt) < fsh) goto L415;
	++l;
	GOBACK L410;
L415:
	a[j] = y[l];
	goto L425;
L420:
	if ((xt-xe) < fsh) GOBACK L415;
	grid[1] = ft;
	a[j] = freq_eval(1,nz,grid,x,y,ad);
L425:
	if (l > 1) l = l-1;
    }

    grid[1] = gtemp;
    dden = TWOPI / cn;
    DOloop (j,1,nfcns) {
	dtemp = 0.0;
	dnum = (j-1) * dden;
	if (nm1 >= 1) {
	    DOloop(k,1,nm1) {
		dtemp += a[k+1] * cos(dnum*k);
	    }
	}
	alpha[j] = 2.0 * dtemp + a[1];
    }

    DOloop(j,2,nfcns) alpha[j] *= 2.0 / cn;
    alpha[1] /= cn;

    if (kkk != 1) {
	p[1] = 2.0*alpha[nfcns]*bb+alpha[nm1];
	p[2] = 2.0*aa*alpha[nfcns];
	q[1] = alpha[nfcns-2]-alpha[nfcns];
	DOloop(j,2,nm1) {
	    if (j >= nm1) {
		aa *= 0.5;
		bb *= 0.5;
	    }
	    p[j+1] = 0.0;
	    DOloop(k,1,j) {
		a[k] = p[k];
		p[k] = 2.0 * bb * a[k];
	    }
	    p[2] += a[1] * 2.0 *aa;
	    jm1 = j - 1;
	    DOloop(k,1,jm1) p[k] += q[k] + aa * a[k+1];
	    jp1 = j + 1;
	    DOloop(k,3,jp1) p[k] += aa * a[k-1];

	    if (j != nm1) {
		DOloop(k,1,j) q[k] = -a[k];
		q[1] += alpha[nfcns - 1 - j];
	    }
	}
	DOloop(j,1,nfcns) alpha[j] = p[j];
    }

    if (nfcns <= 3) {
	  alpha[nfcns+1] = alpha[nfcns+2] = 0.0;
    }
    return 0;
}


/*
 *-----------------------------------------------------------------------
 * FUNCTION: eff
 *  FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
 *  AS A FUNCTION OF FREQUENCY.
 *  AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
 *  APPROXIMATED IF THE USER REPLACES THIS FUNCTION
 *  WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
 *  MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
 *  VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
 *-----------------------------------------------------------------------
 */
static double eff(double freq, double *fx, int lband, int jtype)
{
      if (jtype != 2) return fx[lband];
      else            return fx[lband] * freq;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: wate
 *  FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
 *  OF FREQUENCY.  SIMILAR TO THE FUNCTION eff, THIS FUNCTION CAN
 *  BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
 *  DESIRED WEIGHTING FUNCTION.
 *-----------------------------------------------------------------------
 */
static double wate(double freq, double *fx, double *wtx, int lband, int jtype)
{
      if (jtype != 2)          return wtx[lband];
      if (fx[lband] >= 0.0001) return wtx[lband] / freq;
      return                          wtx[lband];
}

/*********************************************************/

/*  This routine accepts basic input information and puts it in
 *  the form expected by remez.

 *  Adapted from main() by Travis Oliphant
 */

static int pre_remez(double *h2, int numtaps, int numbands, double *bands,
                     double *response, double *weight, int type, int maxiter,
                     int grid_density, int *niter_out) {

  int jtype, nbands, nfilt, lgrid, nz;
  int neg, nodd, nm1;
  int j, k, l, lband, dimsize;
  double delf, change, fup, temp;
  double *tempstor, *edge, *h, *fx, *wtx;
  double *des, *grid, *wt, *alpha, *work;
  double dev;
  int ngrid;
  int *iext;
  int nfcns, wrksize, total_dsize, total_isize;

  lgrid = grid_density;
  dimsize = (int) ceil(numtaps/2.0 + 2);
  wrksize = grid_density * dimsize;
  nfilt = numtaps;
  jtype = type; nbands = numbands;
  /* Note:  code assumes these arrays start at 1 */
  edge = bands-1;
  h = h2 - 1;
  fx = response - 1;
  wtx = weight - 1;

  total_dsize = (dimsize+1)*7 + 3*(wrksize+1);
  total_isize = (dimsize+1);
  /* Need space for:  (all arrays ignore the first element).

     des  (wrksize+1)
     grid (wrksize+1)
     wt   (wrksize+1)
     iext (dimsize+1)   (integer)
     alpha (dimsize+1)
     work  (dimsize+1)*6

  */
  tempstor = (double *)malloc((total_dsize)*sizeof(double)+(total_isize)*sizeof(int));
  if (tempstor == NULL) return -2;

  des = tempstor; grid = des + wrksize+1;
  wt = grid + wrksize+1; alpha = wt + wrksize+1;
  work = alpha + dimsize+1; iext = (int *)(work + (dimsize+1)*6);

  /* Set up problem on dense_grid */

  neg = 1;
  if (jtype == 1) neg = 0;
  nodd = nfilt % 2;
  nfcns = nfilt / 2;
  if (nodd == 1 && neg == 0) nfcns = nfcns + 1;

    /*
     * SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
     * IS (FILTER LENGTH + 1)*GRID DENSITY/2
     */
    grid[1] = edge[1];
    delf = lgrid * nfcns;
    delf = 0.5 / delf;
    if (neg != 0) {
	if (edge[1] < delf) grid[1] = delf;
    }
    j = 1;
    l = 1;
    lband = 1;

    /*
     * CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
     * FUNCTION ON THE GRID
     */
    for (;;) {
	fup = edge[l + 1];
	do {
	    temp = grid[j];
	    des[j] = eff(temp,fx,lband,jtype);
	    wt[j] = wate(temp,fx,wtx,lband,jtype);
	    if (++j > wrksize) {
                /* too many points, or too dense grid */
                free(tempstor);
                return -1;
            }
	    grid[j] = temp + delf;
	} while (grid[j] <= fup);

	grid[j-1] = fup;
	des[j-1] = eff(fup,fx,lband,jtype);
	wt[j-1] = wate(fup,fx,wtx,lband,jtype);
	++lband;
	l += 2;
	if (lband > nbands) break;
	grid[j] = edge[l];
    }

    ngrid = j - 1;
    if (neg == nodd) {
	if (grid[ngrid] > (0.5-delf)) --ngrid;
    }

    /*
     * SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
     * TO THE ORIGINAL PROBLEM
     */
    if (neg <= 0) {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = cos(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j] * change;
	    }
	}
    } else {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = sin(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	} else {
	    DOloop(j,1,ngrid) {
		change = sin(TWOPI * grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	}
    }

    /*XX*/
    temp = (double)(ngrid-1) / (double)nfcns;
    DOloop(j,1,nfcns) {
	iext[j] = (int)((j-1)*temp) + 1; /* round? !! */
    }
    iext[nfcns+1] = ngrid;
    nm1 = nfcns - 1;
    nz  = nfcns + 1;

    if (remez(&dev, des, grid, edge, wt, ngrid, numbands, iext, alpha, nfcns,
              maxiter, work, dimsize, niter_out) < 0) {
        free(tempstor);
        return -1;
    }

    /*
     * CALCULATE THE IMPULSE RESPONSE.
     */
    if (neg <= 0) {

	if (nodd != 0) {
	    DOloop(j,1,nm1) {
		h[j] = 0.5 * alpha[nz-j];
	    }
	    h[nfcns] = alpha[1];
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5*alpha[1] + 0.25*alpha[2];
	}
    } else {
	if (nodd != 0) {
	    h[1] = 0.25 * alpha[nfcns];
	    h[2] = 0.25 * alpha[nm1];
	    DOloop(j,3,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[3];
	    h[nz] = 0.0;
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[2];
	}
    }

    DOloop(j,1,nfcns){
        k = nfilt + 1 - j;
        if (neg == 0)
           h[k] = h[j];
        else
           h[k] = -h[j];
    }
    if (neg == 1 && nodd == 1) h[nz] = 0.0;

  free(tempstor);
  return 0;

}

/**************************************************************
 * End of remez routines
 **************************************************************/


/****************************************************/
/* End of python-independent routines               */
/****************************************************/

/******************************************/

static char doc_correlateND[] = "out = _correlateND(a,kernel,mode) \n\n   mode = 0 - 'valid', 1 - 'same', \n  2 - 'full' (default)";

/*******************************************************************/

static char doc_convolve2d[] = "out = _convolve2d(in1, in2, flip, mode, boundary, fillvalue)";

extern int pylab_convolve_2d(char*, npy_intp*, char*, npy_intp*, char*,
                             npy_intp*, npy_intp*, npy_intp*, int, char*);

static PyObject *_sigtools_convolve2d(PyObject *NPY_UNUSED(dummy), PyObject *args) {

    PyObject *in1=NULL, *in2=NULL, *fill_value=NULL;
    int mode=2, boundary=0, typenum, flag, flip=1, ret;
    npy_intp *aout_dimens=NULL;
    int i;
    PyArrayObject *ain1=NULL, *ain2=NULL, *aout=NULL;
    PyArrayObject *afill=NULL;

    if (!PyArg_ParseTuple(args, "OO|iiiO", &in1, &in2, &flip, &mode, &boundary,
                          &fill_value)) {
        return NULL;
    }

    typenum = PyArray_ObjectType(in1, 0);
    typenum = PyArray_ObjectType(in2, typenum);
    ain1 = (PyArrayObject *)PyArray_FromObject(in1, typenum, 2, 2);
    if (ain1 == NULL) goto fail;
    ain2 = (PyArrayObject *)PyArray_FromObject(in2, typenum, 2, 2);
    if (ain2 == NULL) goto fail;

    if ((boundary != PAD) && (boundary != REFLECT) && (boundary != CIRCULAR))
      PYERR("Incorrect boundary value.");

    if ((boundary == PAD) & (fill_value != NULL)) {
        afill = (PyArrayObject *)PyArray_FromObject(fill_value, typenum, 0, 0);
        if (afill == NULL) {
            /* For backwards compatibility try go via complex */
            PyArrayObject *tmp;
            PyErr_Clear();
            tmp = (PyArrayObject *)PyArray_FromObject(fill_value,
                                                      NPY_CDOUBLE, 0, 0);
            if (tmp == NULL) goto fail;
            afill = (PyArrayObject *)PyArray_Cast(tmp, typenum);
            Py_DECREF(tmp);
            if (afill == NULL) goto fail;
            PYERR("could not cast `fillvalue` directly to the output "
                  "type (it was first converted to complex).");
        }
        if (PyArray_SIZE(afill) != 1) {
            if (PyArray_SIZE(afill) == 0) {
                PyErr_SetString(PyExc_ValueError,
                                "`fillvalue` cannot be an empty array.");
                goto fail;
            }
            PYERR("`fillvalue` must be scalar or an array with "
                  "one element.");
        }
    }
    else {
        /* Create a zero filled array */
        afill = (PyArrayObject *)PyArray_ZEROS(0, NULL, typenum, 0);
        if (afill == NULL) goto fail;
    }

    aout_dimens = (npy_intp *)malloc(PyArray_NDIM(ain1)*sizeof(npy_intp));
    if (aout_dimens == NULL) goto fail;
    switch(mode & OUTSIZE_MASK) {
    case VALID:
    for (i = 0; i < PyArray_NDIM(ain1); i++) {
        aout_dimens[i] = PyArray_DIMS(ain1)[i] - PyArray_DIMS(ain2)[i] + 1;
        if (aout_dimens[i] < 0) {
        PyErr_SetString(PyExc_ValueError,
                    "no part of the output is valid, use option 1 (same) or 2 "
                    "(full) for third argument");
        goto fail;
        }
    }
    break;
    case SAME:
    for (i = 0; i < PyArray_NDIM(ain1); i++) {
            aout_dimens[i] = PyArray_DIMS(ain1)[i];
        }
    break;
    case FULL:
    for (i = 0; i < PyArray_NDIM(ain1); i++) {
            aout_dimens[i] = PyArray_DIMS(ain1)[i] + PyArray_DIMS(ain2)[i] - 1;
        }
    break;
    default:
    PyErr_SetString(PyExc_ValueError,
            "mode must be 0 (valid), 1 (same), or 2 (full)");
    goto fail;
    }

    aout = (PyArrayObject *)PyArray_SimpleNew(PyArray_NDIM(ain1), aout_dimens,
                                              typenum);
    if (aout == NULL) goto fail;

    flag = mode + boundary + (typenum << TYPE_SHIFT) + \
      (flip != 0) * FLIP_MASK;

    ret = pylab_convolve_2d((char *)PyArray_DATA(ain1),      /* Input data Ns[0] x Ns[1] */
                     PyArray_STRIDES(ain1),   /* Input strides */
                     (char *)PyArray_DATA(aout),      /* Output data */
                     PyArray_STRIDES(aout),   /* Output strides */
                     (char *)PyArray_DATA(ain2),      /* coefficients in filter */
                     (npy_intp *)PyArray_STRIDES(ain2),   /* coefficients strides */
                     PyArray_DIMS(ain2),      /* Size of kernel Nwin[2] */
                     PyArray_DIMS(ain1),      /* Size of image Ns[0] x Ns[1] */
                     flag,                    /* convolution parameters */
                     (char *)PyArray_DATA(afill));    /* fill value */


    switch (ret) {
    case 0:
      free(aout_dimens);
      Py_DECREF(ain1);
      Py_DECREF(ain2);
      Py_XDECREF(afill);
      return (PyObject *)aout;
      break;
    case -5:
    case -4:
      PyErr_SetString(PyExc_ValueError,
		      "convolve2d not available for this type.");
      goto fail;
    case -3:
      PyErr_NoMemory();
      goto fail;
    case -2:
      PyErr_SetString(PyExc_ValueError,
		      "Invalid boundary type.");
      goto fail;
    case -1:
      PyErr_SetString(PyExc_ValueError,
		      "Invalid output flag.");
      goto fail;
    }

fail:
    free(aout_dimens);
    Py_XDECREF(ain1);
    Py_XDECREF(ain2);
    Py_XDECREF(aout);
    Py_XDECREF(afill);
    return NULL;
}

/*******************************************************************/

static char doc_remez[] =
    "h = _remez(numtaps, bands, des, weight, type, fs, maxiter, grid_density)\n"
    "  returns the optimal (in the Chebyshev/minimax sense) FIR filter impulse\n"
    "  response given a set of band edges, the desired response on those bands,\n"
    "  and the weight given to the error in those bands.  Bands is a monotonic\n"
    "  vector with band edges given in frequency domain where fs is the sampling\n"
    "  frequency.";

static PyObject *_sigtools_remez(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
    PyObject *bands, *des, *weight;
    int k, numtaps, numbands, type = BANDPASS, err;
    PyArrayObject *a_bands=NULL, *a_des=NULL, *a_weight=NULL;
    PyArrayObject *h=NULL;
    npy_intp ret_dimens; int maxiter = 25, grid_density = 16;
    double oldvalue, *dptr, fs = 1.0;
    char mystr[255];
    int niter = -1;

    if (!PyArg_ParseTuple(args, "iOOO|idii", &numtaps, &bands, &des, &weight,
                          &type, &fs, &maxiter, &grid_density)) {
        return NULL;
    }

    if (type != BANDPASS && type != DIFFERENTIATOR && type != HILBERT) {
        PyErr_SetString(PyExc_ValueError,
                        "The type must be BANDPASS, DIFFERENTIATOR, or HILBERT.");
        return NULL;
	}

    if (numtaps < 2) {
        PyErr_SetString(PyExc_ValueError,
                        "The number of taps must be greater than 1.");
        return NULL;
    }

    a_bands = (PyArrayObject *)PyArray_ContiguousFromObject(bands, NPY_DOUBLE,1,1);
    if (a_bands == NULL) goto fail;
    a_des = (PyArrayObject *)PyArray_ContiguousFromObject(des, NPY_DOUBLE,1,1);
    if (a_des == NULL) goto fail;
    a_weight = (PyArrayObject *)PyArray_ContiguousFromObject(weight, NPY_DOUBLE,1,1);
    if (a_weight == NULL) goto fail;

    numbands = PyArray_DIMS(a_des)[0];
    if ((PyArray_DIMS(a_bands)[0] != 2*numbands) ||
        (PyArray_DIMS(a_weight)[0] != numbands)) {
	    PyErr_SetString(PyExc_ValueError,
                        "The inputs desired and weight must have same length.\n  "
                        "The input bands must have twice this length.");
        goto fail;
    }

    /*
     * Check the bands input to see if it is monotonic, divide by
     * fs to take from range 0 to 0.5 and check to see if in that range
     */
    dptr = (double *)PyArray_DATA(a_bands);
    oldvalue = 0;
    for (k=0; k < 2*numbands; k++) {
        if (*dptr < oldvalue) {
            PyErr_SetString(PyExc_ValueError,
                            "Bands must be monotonic starting at zero.");
            goto fail;
        }
        if (*dptr * 2 > fs) {
            PyErr_SetString(PyExc_ValueError,
                            "Band edges should be less than 1/2 the sampling frequency");
            goto fail;
        }
        oldvalue = *dptr;
        *dptr = oldvalue / fs;  /* Change so that sampling frequency is 1.0 */
        dptr++;
    }

    ret_dimens = numtaps;
    h = (PyArrayObject *)PyArray_SimpleNew(1, &ret_dimens, NPY_DOUBLE);
    if (h == NULL) goto fail;

    err = pre_remez((double *)PyArray_DATA(h), numtaps, numbands,
                    (double *)PyArray_DATA(a_bands),
                    (double *)PyArray_DATA(a_des),
                    (double *)PyArray_DATA(a_weight),
                    type, maxiter, grid_density, &niter);
    if (err < 0) {
        if (err == -1) {
            snprintf(mystr, sizeof(mystr), "Failure to converge at iteration %d, try reducing transition band width.\n", niter);
	        PyErr_SetString(PyExc_ValueError, mystr);
	        goto fail;
        }
        else if (err == -2) {
            PyErr_NoMemory();
            goto fail;
        }
    }

    Py_DECREF(a_bands);
    Py_DECREF(a_des);
    Py_DECREF(a_weight);

	return PyArray_Return(h);

fail:
    Py_XDECREF(a_bands);
    Py_XDECREF(a_des);
    Py_XDECREF(a_weight);
    Py_XDECREF(h);
    return NULL;
}

static char doc_median2d[] = "filt = _median2d(data, size)";

extern void f_medfilt2(float*,float*,npy_intp*,npy_intp*,int*);
extern void d_medfilt2(double*,double*,npy_intp*,npy_intp*,int*);
extern void b_medfilt2(unsigned char*,unsigned char*,npy_intp*,npy_intp*,int*);

static PyObject *_sigtools_median2d(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
    PyObject *image=NULL, *size=NULL;
    int typenum, errnum=-2;
    PyArrayObject *a_image=NULL, *a_size=NULL;
    PyArrayObject *a_out=NULL;
    npy_intp Nwin[2] = {3,3};

    if (!PyArg_ParseTuple(args, "O|O", &image, &size)) return NULL;

    typenum = PyArray_ObjectType(image, 0);
    a_image = (PyArrayObject *)PyArray_ContiguousFromObject(image, typenum, 2, 2);
    if (a_image == NULL) goto fail;

    if (size != NULL) {
	a_size = (PyArrayObject *)PyArray_ContiguousFromObject(size, NPY_INTP, 1, 1);
	if (a_size == NULL) goto fail;
	if ((PyArray_NDIM(a_size) != 1) || (PyArray_DIMS(a_size)[0] < 2))
	    PYERR("Size must be a length two sequence");
	Nwin[0] = ((npy_intp *)PyArray_DATA(a_size))[0];
	Nwin[1] = ((npy_intp *)PyArray_DATA(a_size))[1];
    }

    a_out = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(a_image), typenum);
    if (a_out == NULL) goto fail;

    switch (typenum) {
        case NPY_UBYTE:
            b_medfilt2((unsigned char *)PyArray_DATA(a_image),
                       (unsigned char *)PyArray_DATA(a_out),
                       Nwin, PyArray_DIMS(a_image),
                       &errnum);
            break;
        case NPY_FLOAT:
            f_medfilt2((float *)PyArray_DATA(a_image),
                       (float *)PyArray_DATA(a_out), Nwin,
                       PyArray_DIMS(a_image),
                       &errnum);
            break;
        case NPY_DOUBLE:
            d_medfilt2((double *)PyArray_DATA(a_image),
                       (double *)PyArray_DATA(a_out), Nwin,
                       PyArray_DIMS(a_image),
                       &errnum);
            break;
        default:
            PYERR("2D median filter only supports uint8, float32, and float64.");
    }
    if (errnum != 0) {
        PYERR("ERROR: unable to allocate enough memory in _medfilt2d!\n");
    }

    Py_DECREF(a_image);
    Py_XDECREF(a_size);

    return PyArray_Return(a_out);

 fail:
    Py_XDECREF(a_image);
    Py_XDECREF(a_size);
    Py_XDECREF(a_out);
    return NULL;

}

static char doc_linear_filter[] =
    "(y,Vf) = _linear_filter(b,a,X,Dim=-1,Vi=None)  " \
    "implemented using Direct Form II transposed flow " \
    "diagram. If Vi is not given, Vf is not returned.";

static struct PyMethodDef toolbox_module_methods[] = {
	{"_correlateND", scipy_signal__sigtools_correlateND, METH_VARARGS, doc_correlateND},
	{"_convolve2d", _sigtools_convolve2d, METH_VARARGS, doc_convolve2d},
	{"_linear_filter", scipy_signal__sigtools_linear_filter, METH_VARARGS, doc_linear_filter},
	{"_remez", _sigtools_remez, METH_VARARGS, doc_remez},
	{"_medfilt2d", _sigtools_median2d, METH_VARARGS, doc_median2d},
	{NULL, NULL, 0, NULL}		/* sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_sigtools",
    NULL,
    -1,
    toolbox_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__sigtools(void)
{
    PyObject *module;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    scipy_signal__sigtools_linear_filter_module_init();

    return module;
}
