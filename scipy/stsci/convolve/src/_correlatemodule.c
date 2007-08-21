#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/libnumarray.h"

typedef enum
{
	PIX_NEAREST,
	PIX_REFLECT,
	PIX_WRAP,
	PIX_CONSTANT
} PixMode;

typedef struct
{
	PixMode mode;
	long    rows, cols;
	Float64 constval;
	Float64 *data;
} PixData;

static long 
SlowCoord(long x, long maxx, PixMode m)
{
	switch(m) {
	case PIX_NEAREST:
		if (x < 0) x = 0;
		if (x >= maxx) x = maxx-1;
		return x;
	case PIX_REFLECT:
		if (x < 0) x = -x-1;
		if (x >= maxx) x = maxx - (x - maxx) - 1;
		return x;
	case PIX_WRAP:
		if (x < 0) x += maxx;
		if (x >= maxx) x -= maxx;
		return x;
	case PIX_CONSTANT:  /* handled in SlowPix, suppress warnings */
		break;
	}
	return x;
}

static Float64
SlowPix(long r, long c, PixData *p)
{
	long fr, fc;
	if (p->mode == PIX_CONSTANT) {
		if ((r <  0) || (r >= p->rows) || (c < 0) || (c >= p->cols))
			return p->constval;
		else {
			fr = r;
			fc = c;
		}
	} else {
		fr = SlowCoord(r, p->rows, p->mode);
		fc = SlowCoord(c, p->cols, p->mode);
	}
	return p->data[fr*p->cols + fc];
}

static int
_reject_complex(PyObject *a)
{
	NumarrayType t;
	if ((a == Py_None) || (a == NULL)) 
		return 0;
	t = NA_NumarrayType(a);
	if (t < 0) {
		PyErr_Clear();
		return 0;
	}
	if (t == tComplex32 || t == tComplex64) {
		PyErr_Format(PyExc_TypeError,
			     "function doesn't support complex arrays.");
		return 1;
	}
	return 0;
}

static void 
Correlate1d(long ksizex, Float64 *kernel, 
	   long dsizex, Float64 *data, 
	   Float64 *correlated)
{
	long xc;
	long halfk = ksizex/2;

	for(xc=0; xc<halfk; xc++)
		correlated[xc] = data[xc];

	for(xc=halfk; xc<dsizex-halfk; xc++) {
		long xk;
		double temp = 0;
		for (xk=0; xk<ksizex; xk++)
			temp += kernel[xk]*data[xc-halfk+xk];
		correlated[xc] = temp;
	}
		     
	for(xc=dsizex-halfk; xc<dsizex; xc++)
		correlated[xc] = data[xc];
}

static PyObject *
Py_Correlate1d(PyObject *obj, PyObject *args)
{
	PyObject   *okernel, *odata, *ocorrelated=NULL;
	PyArrayObject *kernel, *data, *correlated;

	if (!PyArg_ParseTuple(args, "OO|O:Correlate1d", 
			      &okernel, &odata, &ocorrelated))
		return NULL;

	/* Align, Byteswap, Contiguous, Typeconvert */
	kernel	= NA_InputArray(okernel, tFloat64, C_ARRAY);
	data	= NA_InputArray(odata, tFloat64, C_ARRAY);
	correlated = NA_OptionalOutputArray(ocorrelated, tFloat64, C_ARRAY, 
					    data);

	if (!kernel || !data || !correlated)
		goto _fail;

	if (_reject_complex(okernel) || _reject_complex(odata) || 
	    _reject_complex(ocorrelated))
		goto _fail;

	if ((kernel->nd != 1) || (data->nd != 1)) {
		PyErr_Format(PyExc_ValueError,
			     "Correlate1d: numarray must have exactly 1 dimension.");
		goto _fail;
	}

	if (!NA_ShapeEqual(data, correlated)) {
		PyErr_Format(PyExc_ValueError,
			     "Correlate1d: data and output must have identical length.");
		goto _fail;
	}

	Correlate1d(kernel->dimensions[0], NA_OFFSETDATA(kernel),
		    data->dimensions[0],   NA_OFFSETDATA(data),
		    NA_OFFSETDATA(correlated));

	Py_DECREF(kernel);
	Py_DECREF(data);

	/* Align, Byteswap, Contiguous, Typeconvert */
	return NA_ReturnOutput(ocorrelated, correlated);

  _fail:
	Py_XDECREF(kernel);
	Py_XDECREF(data);
	Py_XDECREF(correlated);
	return NULL;
}

/* SlowCorrelate computes 2D correlation near the boundaries of an array.
The output array shares the same dimensions as the input array, the latter
fully described by PixData.

The region defined by rmin,rmax,cmin,cmax is assumed to contain only valid
coordinates.  However, access to the input array is performed using SlowPix
because pixels reachable via "kernel offsets" may be at invalid coordinates.
*/
static void 
SlowCorrelate2d(long rmin, long rmax, long cmin, long cmax, 
	      long krows, long kcols, Float64 *kernel, 
	      PixData *pix, Float64 *output)
{
	long kr, kc, r, c;
	long halfkrows = krows/2;
	long halfkcols = kcols/2;

	for(r=rmin; r<rmax; r++) { 
		for(c=cmin; c<cmax; c++) {
			Float64 temp = 0;
			for(kr=0; kr<krows; kr++) {
				long pr = r + kr - halfkrows;
				for(kc=0; kc<kcols; kc++) {
					long pc = c + kc - halfkcols;
					temp += SlowPix(pr, pc, pix) * 
						kernel[kr*kcols+kc];
				}
			}
			output[r*pix->cols+c] = temp;
		}
	}
}

static void 
Correlate2d(long krows, long kcols, Float64 *kernel, 
	    long drows, long dcols, Float64 *data, Float64 *correlated, 
	    PixMode mode, Float64 cval)
{
	long ki, kj, di, dj;
	long halfkrows = krows/2;
	long halfkcols = kcols/2;
	
	PixData pix;
	pix.mode = mode;
	pix.data = data;
	pix.constval = cval;
	pix.rows = drows;
	pix.cols = dcols;

	/* Compute the boundaries using SlowPix */

	SlowCorrelate2d(0, halfkrows, 0, dcols, 
		      krows, kcols, kernel, &pix, correlated); /* top */
	SlowCorrelate2d(drows-halfkrows, drows, 0, dcols, 
		      krows, kcols, kernel, &pix, correlated); /* bottom */
	SlowCorrelate2d(halfkrows, drows-halfkrows, 0, halfkcols, 
		      krows, kcols, kernel, &pix, correlated); /* left */
	SlowCorrelate2d(halfkrows, drows-halfkrows, dcols-halfkcols, dcols,
		      krows, kcols, kernel, &pix, correlated); /* right */

	/* Correlate the center data using unchecked array access */
	for(di=halfkrows; di<drows-halfkrows; di++) {
		for(dj=halfkcols; dj<dcols-halfkcols; dj++) {
			Float64 temp = 0;
			for(ki=0; ki<krows; ki++) {
				long pi = di + ki - halfkrows;
				for(kj=0; kj<kcols; kj++) {
					long pj = dj + kj - halfkcols;
					temp += data[pi*dcols+pj] * 
						kernel[ki*kcols+kj];
				}
			}
			correlated[di*dcols+dj] = temp;
		}
	}
}

static PyObject *
Py_Correlate2d(PyObject *obj, PyObject *args, PyObject *kw)
{
	PyObject      *okernel, *odata, *ocorrelated=NULL;
	PyArrayObject *kernel,  *data,  *correlated;
	Float64  cval = 0;
	int      mode = PIX_NEAREST;
	char       *keywds[] = { "kernel", "data", "output", "mode", "cval", NULL };

	if (!PyArg_ParseTupleAndKeywords(
		    args, kw, "OO|Oid:Correlate2d", keywds, 
		    &okernel, &odata, &ocorrelated, &mode, &cval))
		return NULL;

	if ((mode < PIX_NEAREST) || (mode > PIX_CONSTANT))
		return PyErr_Format(PyExc_ValueError,
			    "Correlate2d: mode value not in range(%d,%d)", 
				    PIX_NEAREST, PIX_CONSTANT);

	/* Align, Byteswap, Contiguous, Typeconvert */
	kernel	= NA_InputArray(okernel, tFloat64, C_ARRAY);
	data	= NA_InputArray(odata, tFloat64, C_ARRAY);
	correlated = NA_OptionalOutputArray(ocorrelated, tFloat64, C_ARRAY, 
					    data);

	if (!kernel || !data || !correlated)
		goto _fail;

	if ((kernel->nd != 2) || (data->nd != 2) || (correlated->nd != 2)) {
		PyErr_Format(PyExc_ValueError, "Correlate2d: inputs must have 2 dimensions.");
		goto _fail;
	}

	if (!NA_ShapeEqual(data, correlated)) {
		PyErr_Format(PyExc_ValueError,
			     "Correlate2d: data and output numarray need identical shapes.");
		goto _fail;
	}	

	if (_reject_complex(okernel) || _reject_complex(odata) || 
	    _reject_complex(ocorrelated))
		goto _fail;
		
	Correlate2d(kernel->dimensions[0], kernel->dimensions[1], 
		    NA_OFFSETDATA(kernel),
		    data->dimensions[0], data->dimensions[1], 
		    NA_OFFSETDATA(data),
		    NA_OFFSETDATA(correlated), 
		    mode, cval);

	Py_DECREF(kernel);
	Py_DECREF(data);

	/* Align, Byteswap, Contiguous, Typeconvert */
	return NA_ReturnOutput(ocorrelated, correlated);

  _fail:
	Py_XDECREF(kernel);
	Py_XDECREF(data);
	Py_XDECREF(correlated);
	return NULL;
}

void Shift2d( long rows, long cols, Float64 *data, long dx, long dy, Float64 *output, int mode, Float64 cval)
{
	long r, c;
	PixData pix;
	pix.mode = mode;
	pix.constval = cval;
	pix.rows = rows;
	pix.cols = cols;
	pix.data = data;

	for(r=0; r<rows; r++)
		for(c=0; c<cols; c++)
			output[ r*cols + c] = SlowPix(r+dy, c+dx, &pix);
}

static PyObject *
Py_Shift2d(PyObject *obj, PyObject *args, PyObject *kw)
{
	PyObject      *odata, *ooutput=NULL;
	PyArrayObject *data,  *output;
	int           dx, dy;
	Float64       cval = 0;
	int           mode = PIX_NEAREST;
	char          *keywds[] = { "data", "dx", "dy", 
				    "output", "mode", "cval", NULL };

	if (!PyArg_ParseTupleAndKeywords(args, kw, "Oii|Oid:Shift2d", keywds,
			 &odata, &dx, &dy, &ooutput, &mode, &cval))
		return NULL;

	if ((mode < PIX_NEAREST) || (mode > PIX_CONSTANT))
		return PyErr_Format(PyExc_ValueError,
			    "Shift2d: mode value not in range(%d,%d)", 
				    PIX_NEAREST, PIX_CONSTANT);

	/* Align, Byteswap, Contiguous, Typeconvert */
	data	= NA_InputArray(odata, tFloat64, C_ARRAY);
	output  = NA_OptionalOutputArray(ooutput, tFloat64, C_ARRAY, 
					 data);

	if (!data || !output)
		goto _fail;

	if (_reject_complex(odata) || _reject_complex(ooutput))
		goto _fail;

	if ((data->nd != 2)) {
		PyErr_Format(PyExc_ValueError,
			     "Shift2d: numarray must have 2 dimensions.");
		goto _fail;
	}

	if (!NA_ShapeEqual(data, output)) {
		PyErr_Format(PyExc_ValueError,
			     "Shift2d: data and output numarray need identical shapes.");
		goto _fail;
	}

	/* Invert sign of deltas to match sense of 2x2 correlation. */
	Shift2d( data->dimensions[0], data->dimensions[1], NA_OFFSETDATA(data),
		 -dx, -dy, NA_OFFSETDATA(output), mode, cval);
	
	Py_XDECREF(data);

	/* Align, Byteswap, Contiguous, Typeconvert */
	return NA_ReturnOutput(ooutput, output);
  _fail:
	Py_XDECREF(data);
	Py_XDECREF(output);
	return NULL;
}

typedef struct s_BoxData BoxData;

typedef Float64 (*SumColFunc)(long,long,BoxData*);
typedef Float64 (*SumBoxFunc)(long,long,BoxData*);

struct s_BoxData {
	PixData    pix;
	long       krows, kcols;
	SumColFunc sumcol;
	SumBoxFunc sumbox;
};

static Float64
SlowSumCol(long r, long c, BoxData *D) 
{
	Float64 sum = 0;
	long i, krows = D->krows;
	for(i=0; i<krows; i++) {
		sum += SlowPix(r+i, c, &D->pix);
	}
	return sum;
}

static Float64
SlowSumBox(long r, long c, BoxData *D)
{
	long i, j;
	Float64 sum = 0;
	for(i=0; i<D->krows; i++)
		for(j=0; j<D->kcols; j++)
			sum += SlowPix(r+i, c+j, &D->pix);
	return sum;
}

static Float64
FastSumCol(long r, long c, BoxData *D) 
{
	Float64 sum = 0;
	long krows = D->krows;
	long cols = D->pix.cols;
	Float64 *data = D->pix.data;

	data += r*cols + c;
	for(; krows--; data += cols) {
		sum += *data;
	}
	return sum;
}

static Float64
FastSumBox(long r, long c, BoxData *D)
{
	long i, j;
	Float64 sum = 0;
	long cols = D->pix.cols;
	Float64 *data = D->pix.data;
	data += r*cols + c;
	for(i=0; i<D->krows; i++, data += cols-D->kcols)
		for(j=0; j<D->kcols; j++, data++)
			sum += *data;
	return sum;
}

static long bound(long x, long max)
{
	if (x < 0) return 0;
	else if (x > max) return max;
	else return x;
}

static void
BoxFunc(long rmin, long rmax, long cmin, long cmax, Float64 *output, BoxData *D)
{
	long r, c;
	long    krows2    = D->krows/2;
	long    kcols2    = D->kcols/2;
	long    kcolseven = !(D->kcols & 1);
	long    rows      = D->pix.rows;
	long    cols      = D->pix.cols;

	rmin = bound(rmin, rows);
	rmax = bound(rmax, rows);
	cmin = bound(cmin, cols);
	cmax = bound(cmax, cols);

	for(r=rmin; r<rmax; r++) {  
		Float64 sum = D->sumbox(r - krows2, cmin - kcols2, D);
		for(c=cmin; c<cmax; c++) {
			output[r*cols + c] = sum;
			sum -= D->sumcol(r - krows2, c - kcols2, D);
			sum += D->sumcol(r - krows2, c + kcols2 - kcolseven + 1, D);
		}
	}
}

/*  BoxFuncI computes a boxcar incrementally, using a formula independent of
    the size of the boxcar.  Each incremental step is based on dropping a
    whole column of the "back" of the boxcar, and adding in a new column in
    the "front".  The sums of these columns are further optimized by realizing
    they can be computed from their counterparts one element above by adding in
    bottom corners and subtracting out top corners.

    incremental pixel layout:      B  C     where S is the unknown, and A, B, C are known neighbors
                                   A  S           each of these refers to the output array

    S = A + a1 - a0            where a0 and a1 are column vectors with *bottom* elements { a, d }
    C - B = b1 - b0            where b0 and b1 are column vectors with *top* elements { b, g }
                                    column vectors and corner elements refer to the input array

    offset matrix layout:	   b      g                   where b is actually in b0
				  [b0]   [b1]                       g    "           b1
                                  [a0] S [a1]                       a    "           a0
				   a      d                         d is actually in a1
 
    a0 = b0 - b + a            column vector a0 is b0 dropping top element b and adding bottom a
    a1 = b1 - g + d            column vector a1 is b1 dropping top element g and adding bottom d

    S = A + (b1 - g + f) - (b0 - b + a)    by substitution
    S = A + (b1 - b0) - g + d + b - a      rearranging additions
    S = A + C - B - g + d + b - a          by substitution
*/
static void
BoxFuncI(long rmin, long rmax, long cmin, long cmax, Float64 *output, BoxData *D)
{
	long r, c;
	long krows2    = D->krows/2;
	long kcols2    = D->kcols/2;
	long krowseven = !(D->krows & 1);
	long kcolseven = !(D->kcols & 1);
	long rows      = D->pix.rows;
	long cols      = D->pix.cols;
	Float64 *input = D->pix.data;

	rmin = bound(rmin, rows);
	rmax = bound(rmax, rows);
	cmin = bound(cmin, cols);
	cmax = bound(cmax, cols);

	for(r=rmin; r<rmax; r++) {  
		long top    = r - krows2 - 1;
		long bottom = r + krows2 - krowseven;

		for(c=cmin; c<cmax; c++) {
			long left   = c - kcols2 - 1;
			long right  = c + kcols2 - kcolseven;

			Float64 A = output [     r  * cols + (c-1) ];
			Float64 B = output [  (r-1) * cols + (c-1) ];
			Float64 C = output [  (r-1) * cols + c     ];
			Float64 a = input  [ bottom * cols + left  ];
			Float64 b = input  [    top * cols + left  ];
			Float64 g = input  [    top * cols + right ];
			Float64 d = input  [ bottom * cols + right ];

			output[r*cols + c] = A + C - B - g + d + b - a;
		}
	}
}

static void 
Boxcar2d(long krows, long kcols, long rows, long cols, Float64 *data, 
	 Float64 *output, PixMode mode, Float64 constval)
{
	long krows2 = krows/2;
	long krowseven = !(krows&1);
	long kcols2 = kcols/2;
	long kcolseven = !(kcols&1);
	long r, c;
	Float64 karea;

	BoxData D;
	D.pix.mode = mode;
	D.pix.constval = constval;
	D.pix.rows = rows;
	D.pix.cols = cols;
	D.pix.data = data;
	D.krows = krows;
	D.kcols = kcols;
	D.sumcol = SlowSumCol;
	D.sumbox = SlowSumBox;
	
	/* The next 4 calls compute boxcars on the boundary pixels with
	   different modes detemining what data values are fetched when "out
	   of bounds".  Presumably, this is a small minority of the data and
	   therefore the implementation inefficiency is not *too* damaging.
	   It should at least beat making 2 outsized copies of the data array. */

	/* top whole plus one */
	BoxFunc(0, krows2+2, 0, cols, output, &D);
	
	/* bottom whole */
	BoxFunc(rows-krows2+krowseven, rows, 0, cols, output, &D);
	
	/* left whole plus one */
	BoxFunc(0, rows, 0, kcols2+2, output, &D);
	
	/* right whole */
	BoxFunc(0, rows, cols-kcols2+kcolseven, cols, output, &D);

	/* Do the boxcar on the "center" data, the data where the boxcar is
	   always "in bounds".  Presumably, this is the bulk of the data so
	   this should be fast.  
	*/
	D.sumcol = FastSumCol;
	D.sumbox = FastSumBox;

	BoxFuncI( krows2+2, rows-krows2+krowseven, 
		  kcols2+2, cols-kcols2+kcolseven, 
		  output, &D);

	karea = kcols * krows;
	for(r=0; r<rows; r++)
		for(c=0; c<cols; c++)
			output[r*cols + c] /= karea;
}

static PyObject *
Py_Boxcar2d(PyObject *obj, PyObject *args, PyObject *kw)
{
	PyObject   *odata, *ooutput=NULL;
	PyArrayObject *data, *output;
	int        krows, kcols, mode=PIX_NEAREST;
	Float64    cval = 0;
	char       *keywds[] = { "data", "krows", "kcols", 
				 "output", "mode", "cval", NULL };

	if (!PyArg_ParseTupleAndKeywords(args, kw, "Oii|Oid:Boxcar2d", keywds, 
			 &odata, &krows, &kcols, &ooutput, &mode, &cval))
		return NULL;

	/* Align, Byteswap, Contiguous, Typeconvert */
	data	= NA_InputArray(odata, tFloat64, C_ARRAY);
	output = NA_OptionalOutputArray(ooutput, tFloat64, C_ARRAY, data);
	if (!data || !output)
		goto _fail;

	if (_reject_complex(odata) || _reject_complex(ooutput))
		goto _fail;

	if ((krows < 0) || (kcols < 0)) {
		PyErr_Format(PyExc_ValueError, "krows and kcols must be > 0.");
		goto _fail;
	}

	if ((mode < PIX_NEAREST) || (mode > PIX_CONSTANT)) {
		PyErr_Format(PyExc_ValueError,
			     "Boxcar2d: mode value not in range(%d,%d)", 
			     PIX_NEAREST, PIX_CONSTANT);
		goto _fail;
	}

	if ((data->nd != 2)|| (output->nd != 2)) {
		PyErr_Format(PyExc_ValueError,
			     "Boxcar2d: numarray must have 2 dimensions.");
		goto _fail;
	}

	if (!NA_ShapeEqual(data, output)) {
		PyErr_Format(PyExc_ValueError,
			     "Boxcar2d: data and output numarray need identical shapes.");
		goto _fail;
	}

	if ((kcols <=0) || (krows <= 0)) {
		PyErr_Format(PyExc_ValueError,
			     "Boxcar2d: invalid data shape.");
		goto _fail;
	}
	if ((kcols > data->dimensions[1]) || (krows > data->dimensions[0])) {
		PyErr_Format(PyExc_ValueError, "Boxcar2d: boxcar shape incompatible with"
			     " data shape.");
		goto _fail;
	}

	Boxcar2d(krows, kcols, data->dimensions[0], data->dimensions[1], 
		 NA_OFFSETDATA(data), NA_OFFSETDATA(output), mode, cval);

	Py_XDECREF(data);

	/* Align, Byteswap, Contiguous, Typeconvert */
	return NA_ReturnOutput(ooutput, output);
  _fail:
	Py_XDECREF(data);
	Py_XDECREF(output);
	return NULL;
}

static PyMethodDef _correlateMethods[] = {
    {"Correlate1d", Py_Correlate1d, METH_VARARGS}, 
    {"Correlate2d", (PyCFunction) Py_Correlate2d, METH_VARARGS | METH_KEYWORDS},
    {"Shift2d", (PyCFunction) Py_Shift2d, METH_VARARGS | METH_KEYWORDS, 
     "Shift2d shifts and image by an integer number of pixels, and uses IRAF compatible modes for the boundary pixels."},
    {"Boxcar2d", (PyCFunction) Py_Boxcar2d, METH_VARARGS | METH_KEYWORDS,
    "Boxcar2d computes a sliding 2D boxcar average on a 2D array"},
    {NULL, NULL} /* Sentinel */
};

PyMODINIT_FUNC init_correlate(void)
{
	PyObject *m, *d;
	m = Py_InitModule("_correlate", _correlateMethods);
	d = PyModule_GetDict(m);
	import_libnumarray();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
