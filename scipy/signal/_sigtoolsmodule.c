/* SIGTOOLS module by Travis Oliphant

Copyright 2005 Travis Oliphant
Permission to use, copy, modify, and distribute this software without fee
is granted under the SciPy License.
*/
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_signal_ARRAY_API
#include "numpy/ndarrayobject.h"
#include "npy_2_compat.h"

#include "_sigtools.h"
#include <setjmp.h>
#include <stdlib.h>

#define PYERR(message) {PyErr_SetString(PyExc_ValueError, message); goto fail;}


jmp_buf MALLOC_FAIL;

char *check_malloc(size_t size)
{
    char *the_block = malloc(size);
    if (the_block == NULL) {
        printf("\nERROR: unable to allocate %zu bytes!\n", size);
        longjmp(MALLOC_FAIL,-1);
    }
    return the_block;
}


/************************************************************************
 * Start of portable, non-python specific routines.                     *
 ************************************************************************/

/* Some core routines are written
in a portable way so that they could be used in other applications.  The
order filtering, however uses python-specific constructs in its guts
and is therefore Python dependent.  This could be changed in a
straightforward way but I haven't done it for lack of time.*/

static int index_out_of_bounds(npy_intp *indices, npy_intp *max_indices, int ndims) {
  int bad_index = 0, k = 0;

  while (!bad_index && (k++ < ndims)) {
    bad_index = ((*(indices) >= *(max_indices++)) || (*(indices) < 0));
    indices++;
  }
  return bad_index;
}

/* This maybe could be redone with stride information so it could be
 * called with non-contiguous arrays:  I think offsets is related to
 * the difference between the strides.  I'm not sure about init_offset
 * just yet.  I think it needs to be calculated because of mode_dep
 * but probably with dim1 being the size of the "original, unsliced" array
 */

static npy_intp compute_offsets (npy_uintp *offsets, npy_intp *offsets2, npy_intp *dim1,
                             npy_intp *dim2, npy_intp *dim3, npy_intp *mode_dep,
                             int nd) {
  int k,i;
  npy_intp init_offset = 0;

  for (k = 0; k < nd - 1; k++)
    {
      init_offset += mode_dep[k];
      init_offset *= dim1[k+1];
    }
  init_offset += mode_dep[k] - 2;

  k = nd;
  while(k--) {
    offsets[k] = 0;
    offsets2[k] = 0;
    for (i = k + 1; i < nd - 1; i++) {
      offsets[k] += dim1[i] - dim2[i];
      offsets[k] *= dim1[i+1];

      offsets2[k] += dim1[i] - dim3[i];
      offsets2[k] *= dim1[i+1];
    }

    if (k < nd - 1) {
      offsets[k] += dim1[i] - dim2[i];
      offsets2[k] += dim1[i] - dim3[i];
    }
    offsets[k] += 1;
    offsets2[k] += 1;
  }
  return init_offset;
}

/* increment by 1 the index into an N-D array, doing the necessary
   carrying when the index reaches the dimension along that axis */
static int increment(npy_intp *ret_ind, int nd, npy_intp *max_ind) {
    int k, incr = 1;

    k = nd - 1;
    if (++ret_ind[k] >= max_ind[k]) {
      while (k >= 0 && (ret_ind[k] >= max_ind[k]-1)) {
	incr++;
	ret_ind[k--] = 0;
      }
      if (k >= 0) ret_ind[k]++;
    }
    return incr;
}

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
  tempstor = malloc((total_dsize)*sizeof(double)+(total_isize)*sizeof(int));
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

/************************/
/* N-D Order Filtering. */


static void fill_buffer(char *ip1, PyArrayObject *ap1, PyArrayObject *ap2,
                        char *sort_buffer, int nels2, int check,
                        npy_intp *loop_ind, npy_intp *temp_ind, npy_uintp *offset){
  int i, k, incr = 1;
  int ndims = PyArray_NDIM(ap1);
  npy_intp *dims2 = PyArray_DIMS(ap2);
  npy_intp *dims1 = PyArray_DIMS(ap1);
  npy_intp is1 = PyArray_ITEMSIZE(ap1);
  npy_intp is2 = PyArray_ITEMSIZE(ap2);
  char *ip2 = PyArray_DATA(ap2);
  int elsize = PyArray_ITEMSIZE(ap1);
  char *ptr;

  i = nels2;
  ptr = PyArray_Zero(ap2);
  temp_ind[ndims-1]--;
  while (i--) {
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ip1 += offset[k]*is1;               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims)) && \
	memcmp(ip2, ptr, PyArray_ITEMSIZE(ap2))) {
      memcpy(sort_buffer, ip1, elsize);
      sort_buffer += elsize;
    }
    /* Returns number of N-D indices incremented. */
    incr = increment(loop_ind, ndims, dims2);
    ip2 += is2;

  }
  PyDataMem_FREE(ptr);
  return;
}

#define COMPARE(fname, type) \
int fname(type *ip1, type *ip2) { return *ip1 < *ip2 ? -1 : *ip1 == *ip2 ? 0 : 1; }

COMPARE(DOUBLE_compare, double)
COMPARE(FLOAT_compare, float)
COMPARE(LONGDOUBLE_compare, npy_longdouble)
COMPARE(BYTE_compare, npy_byte)
COMPARE(SHORT_compare, short)
COMPARE(INT_compare, int)
COMPARE(LONG_compare, long)
COMPARE(LONGLONG_compare, npy_longlong)
COMPARE(UBYTE_compare, npy_ubyte)
COMPARE(USHORT_compare, npy_ushort)
COMPARE(UINT_compare, npy_uint)
COMPARE(ULONG_compare, npy_ulong)
COMPARE(ULONGLONG_compare, npy_ulonglong)


int OBJECT_compare(PyObject **ip1, PyObject **ip2) {
        /* PyObject_RichCompareBool returns -1 on error; not handled here */
        if(PyObject_RichCompareBool(*ip1, *ip2, Py_LT) == 1)
          return -1;
        else if(PyObject_RichCompareBool(*ip1, *ip2, Py_EQ) == 1)
          return 0;
        else
          return 1;
}

typedef int (*CompareFunction)(const void *, const void *);

CompareFunction compare_functions[] = \
	{NULL, (CompareFunction)BYTE_compare,(CompareFunction)UBYTE_compare,\
	 (CompareFunction)SHORT_compare,(CompareFunction)USHORT_compare, \
	 (CompareFunction)INT_compare,(CompareFunction)UINT_compare, \
	 (CompareFunction)LONG_compare,(CompareFunction)ULONG_compare, \
	 (CompareFunction)LONGLONG_compare,(CompareFunction)ULONGLONG_compare,
	 (CompareFunction)FLOAT_compare,(CompareFunction)DOUBLE_compare,
	 (CompareFunction)LONGDOUBLE_compare, NULL, NULL, NULL,
	 (CompareFunction)OBJECT_compare, NULL, NULL, NULL};

PyObject *PyArray_OrderFilterND(PyObject *op1, PyObject *op2, int order) {
	PyArrayObject *ap1=NULL, *ap2=NULL, *ret=NULL;
	npy_intp *a_ind=NULL, *b_ind=NULL, *temp_ind=NULL, *mode_dep=NULL, *check_ind=NULL;
	npy_uintp *offsets=NULL;
	npy_intp *offsets2=NULL;
	npy_uintp offset1;
	int i, n2, n2_nonzero, k, check, incr = 1;
	int typenum, bytes_in_array;
	int is1, os;
	char *op, *ap1_ptr, *ap2_ptr, *sort_buffer=NULL;
	npy_intp *ret_ind=NULL;
	CompareFunction compare_func=NULL;
	char *zptr=NULL;
	PyArray_CopySwapFunc *copyswap;

	/* Get Array objects from input */
	typenum = PyArray_ObjectType(op1, 0);
	typenum = PyArray_ObjectType(op2, typenum);

	ap1 = (PyArrayObject *)PyArray_ContiguousFromObject(op1, typenum, 0, 0);
	if (ap1 == NULL) return NULL;
	ap2 = (PyArrayObject *)PyArray_ContiguousFromObject(op2, typenum, 0, 0);
	if (ap2 == NULL) goto fail;

	if (PyArray_NDIM(ap1) != PyArray_NDIM(ap2)) {
	    PyErr_SetString(PyExc_ValueError,
                "All input arrays must have the same number of dimensions.");
	  goto fail;
	}

	n2 = PyArray_Size((PyObject *)ap2);
	n2_nonzero = 0;
	ap2_ptr = PyArray_DATA(ap2);
        /*
         * Find out the number of non-zero entries in domain (allows for
         * different shapped rank-filters to be used besides just rectangles)
         */
	zptr = PyArray_Zero(ap2);
	if (zptr == NULL) goto fail;
	for (k=0; k < n2; k++) {
	  n2_nonzero += (memcmp(ap2_ptr,zptr,PyArray_ITEMSIZE(ap2)) != 0);
	  ap2_ptr += PyArray_ITEMSIZE(ap2);
	}

	if ((order >= n2_nonzero) || (order < 0)) {
	    PyErr_SetString(PyExc_ValueError,
                "Order must be non-negative and less than number of nonzero elements in domain.");
	  goto fail;
	}

	ret = (PyArrayObject *)PyArray_SimpleNew(PyArray_NDIM(ap1),
                                                 PyArray_DIMS(ap1),
                                                 typenum);
	if (ret == NULL) goto fail;

	if (PyArray_TYPE(ap1) < sizeof(compare_functions) / sizeof(compare_functions[0])) {
	    compare_func = compare_functions[PyArray_TYPE(ap1)];
	}
	if (compare_func == NULL) {
	    PyErr_SetString(PyExc_ValueError,
                        "order_filterND not available for this type");
		goto fail;
	}

	is1 = PyArray_ITEMSIZE(ap1);

	if (!(sort_buffer = malloc(n2_nonzero*is1))) goto fail;

	os = PyArray_ITEMSIZE(ret);
	op = PyArray_DATA(ret);

	copyswap = PyDataType_GetArrFuncs(PyArray_DESCR(ret))->copyswap;

	bytes_in_array = PyArray_NDIM(ap1)*sizeof(npy_intp);
	mode_dep = malloc(bytes_in_array);
	if (mode_dep == NULL) goto fail;
	for (k = 0; k < PyArray_NDIM(ap1); k++) {
	  mode_dep[k] = -((PyArray_DIMS(ap2)[k]-1) >> 1);
	}

	b_ind = (npy_intp *)malloc(bytes_in_array);  /* loop variables */
	if (b_ind == NULL) goto fail;
	memset(b_ind,0,bytes_in_array);
	a_ind = (npy_intp *)malloc(bytes_in_array);
	ret_ind = (npy_intp *)malloc(bytes_in_array);
	if (a_ind == NULL || ret_ind == NULL) goto fail;
	memset(ret_ind,0,bytes_in_array);
	temp_ind = (npy_intp *)malloc(bytes_in_array);
	check_ind = (npy_intp*)malloc(bytes_in_array);
	offsets = (npy_uintp *)malloc(PyArray_NDIM(ap1)*sizeof(npy_uintp));
	offsets2 = (npy_intp *)malloc(PyArray_NDIM(ap1)*sizeof(npy_intp));
	if (temp_ind == NULL || check_ind == NULL || offsets == NULL || offsets2 == NULL) goto fail;
	offset1 = compute_offsets(offsets, offsets2, PyArray_DIMS(ap1),
                                  PyArray_DIMS(ap2), PyArray_DIMS(ret),
                                  mode_dep, PyArray_NDIM(ap1));
	/* The filtering proceeds by looping through the output array
	   and for each value filling a buffer from the
	   element-by-element product of the two input arrays.  The buffer
	   is then sorted and the order_th element is kept as output. Index
	   counters are used for book-keeping in the area so that we
	   can tell where we are in all of the arrays and be sure that
	   we are not trying to access areas outside the arrays definition.

	   The inner loop is implemented separately but equivalently for each
	   datatype. The outer loop is similar in structure and form to
	   to the inner loop.
	*/
	/* Need to keep track of a ptr to place in big (first) input
	   array where we start the multiplication (we pass over it in the
	   inner loop (and not dereferenced)
	   if it is pointing outside dataspace)
	*/
	/* Calculate it once and the just move it around appropriately */
	PyDataMem_FREE(zptr);
	zptr = PyArray_Zero(ap1);
	if (zptr == NULL) goto fail;
	ap1_ptr = (char *)PyArray_DATA(ap1) + offset1*is1;
	for (k=0; k < PyArray_NDIM(ap1); k++) {
            a_ind[k] = mode_dep[k];
            check_ind[k] = PyArray_DIMS(ap1)[k] - PyArray_DIMS(ap2)[k] - mode_dep[k] - 1;
        }
	a_ind[PyArray_NDIM(ap1)-1]--;
	i = PyArray_Size((PyObject *)ret);
	while (i--) {
          /*
           * Zero out the sort_buffer (has effect of zero-padding
           * on boundaries). Treat object arrays right.
           */
	  ap2_ptr = sort_buffer;
	  for (k=0; k < n2_nonzero; k++) {
  	    memcpy(ap2_ptr,zptr,is1);
	    ap2_ptr += is1;
	  }

	  k = PyArray_NDIM(ap1) - 1;
	  while(--incr) {
	    a_ind[k] -= PyArray_DIMS(ret)[k] - 1;   /* Return to start */
	    k--;
	  }
	  ap1_ptr += offsets2[k]*is1;
	  a_ind[k]++;
	  memcpy(temp_ind, a_ind, bytes_in_array);

	  check = 0; k = -1;
	  while(!check && (++k < PyArray_NDIM(ap1)))
	    check = (check || (ret_ind[k] < -mode_dep[k]) ||
                     (ret_ind[k] > check_ind[k]));

	  fill_buffer(ap1_ptr,ap1,ap2,sort_buffer,n2,check,b_ind,temp_ind,offsets);
	  qsort(sort_buffer, n2_nonzero, is1, compare_func);

	  /*
	   * Use copyswap for correct refcounting with object arrays
	   * (sort_buffer has borrowed references, op owns references). Note
	   * also that os == PyArray_ITEMSIZE(ret) and we are copying a single
	   * scalar here.
	   */
	  copyswap(op, sort_buffer + order*is1, 0, NULL);

          /* increment index counter */
	  incr = increment(ret_ind,PyArray_NDIM(ret),PyArray_DIMS(ret));
          /* increment to next output index */
	  op += os;

	}
	free(b_ind); free(a_ind); free(ret_ind);
	free(offsets); free(offsets2); free(temp_ind);
	free(check_ind); free(mode_dep);
	free(sort_buffer);

	PyDataMem_FREE(zptr);
	Py_DECREF(ap1);
	Py_DECREF(ap2);

	return PyArray_Return(ret);

fail:
	if (zptr) PyDataMem_FREE(zptr);
	free(sort_buffer);
	free(mode_dep);
	free(b_ind);
	free(a_ind);
	free(ret_ind);
	free(temp_ind);
	free(check_ind);
    free(offsets);
	free(offsets2);
	Py_XDECREF(ap1);
	Py_XDECREF(ap2);
	Py_XDECREF(ret);
	return NULL;
}


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

    aout_dimens = malloc(PyArray_NDIM(ain1)*sizeof(npy_intp));
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

    ret = pylab_convolve_2d (PyArray_DATA(ain1),      /* Input data Ns[0] x Ns[1] */
		             PyArray_STRIDES(ain1),   /* Input strides */
		             PyArray_DATA(aout),      /* Output data */
		             PyArray_STRIDES(aout),   /* Output strides */
		             PyArray_DATA(ain2),      /* coefficients in filter */
		             PyArray_STRIDES(ain2),   /* coefficients strides */
		             PyArray_DIMS(ain2),      /* Size of kernel Nwin[2] */
			     PyArray_DIMS(ain1),      /* Size of image Ns[0] x Ns[1] */
		             flag,                    /* convolution parameters */
		             PyArray_DATA(afill));    /* fill value */


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

static char doc_order_filterND[] = "out = _order_filterND(a,domain,order)";

static PyObject *_sigtools_order_filterND(PyObject *NPY_UNUSED(dummy),
                                         PyObject *args) {
	PyObject *domain, *a0;
	int order=0;

	if (!PyArg_ParseTuple(args, "OO|i", &a0, &domain, &order)) return NULL;

	return PyArray_OrderFilterND(a0, domain, order);
}


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
            sprintf(mystr, "Failure to converge at iteration %d, try reducing transition band width.\n", niter);
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

extern void f_medfilt2(float*,float*,npy_intp*,npy_intp*);
extern void d_medfilt2(double*,double*,npy_intp*,npy_intp*);
extern void b_medfilt2(unsigned char*,unsigned char*,npy_intp*,npy_intp*);

static PyObject *_sigtools_median2d(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
    PyObject *image=NULL, *size=NULL;
    int typenum;
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

    if (setjmp(MALLOC_FAIL)) {
	PYERR("Memory allocation error.");
    }
    else {
	switch (typenum) {
	case NPY_UBYTE:
	    b_medfilt2((unsigned char *)PyArray_DATA(a_image),
                       (unsigned char *)PyArray_DATA(a_out),
                       Nwin, PyArray_DIMS(a_image));
	    break;
	case NPY_FLOAT:
	    f_medfilt2((float *)PyArray_DATA(a_image),
                       (float *)PyArray_DATA(a_out), Nwin,
                       PyArray_DIMS(a_image));
	    break;
	case NPY_DOUBLE:
	    d_medfilt2((double *)PyArray_DATA(a_image),
                       (double *)PyArray_DATA(a_out), Nwin,
                       PyArray_DIMS(a_image));
	    break;
	default:
	  PYERR("2D median filter only supports uint8, float32, and float64.");
	}
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
	{"_order_filterND", _sigtools_order_filterND, METH_VARARGS, doc_order_filterND},
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

    scipy_signal__sigtools_linear_filter_module_init();

    return module;
}
