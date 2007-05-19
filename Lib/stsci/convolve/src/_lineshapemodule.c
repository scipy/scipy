/* C implementations of various lineshape functions
 *
 * Copyright (C) 2002,2003 Jochen Küpper <jochen@jochen-kuepper.de>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Python.h"

#include <math.h>

#include "numpy/libnumarray.h"


#define sqr(x) ((x)*(x))


/* These are apparently not defined in MSVC */
#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif
#if !defined(M_LN2)
#define M_LN2 0.69314718055994530942
#endif


/*** C implementation ***/

static void gauss(size_t n, double *x, double *y, double w, double xc) 
    /* Evaluate normalized Gauss profile around xc with FWHM w at all x_i,
       return in y. */
{
    int i;
    for(i=0; i<n; i++)
        y[i] = 2. * sqrt(M_LN2/M_PI) / w * exp(-4.*M_LN2 * sqr((x[i]-xc)/w));
}



static void lorentz(size_t n, double *x, double *y, double w, double xc) 
    /* Evaluate normalized Lorentz profile around xc with FWHM w at all x_i,
       return in y. */
{
    int i;
    for(i=0; i<n; i++)
        y[i] = 2.*w/M_PI / (sqr(w) + 4.*(sqr(x[i]-xc)));
}



static double humlicek_v12(double x, double y)
    /** Approximation of Voigt profile by Humlicek's 12-point formula.
     *
     * J. Humlicek, J. Quant. Spectrosc. Radiat. Transfer, 21(1978), 309.
     *
     * Voigt-Profil:
     * V(x, y) = 2/pi^(1.5) * y^2/FWHM_L * \int[-inf,+inf](e^(-y^2)/(x+y)^2+...)
     */
{
    static const double T_v12[6] = {
	0.314240376254359,    0.947788391240164,    1.597682635152605,
	2.27950708050106,     3.020637025120890,    3.889724897869782
    };
    static const double alpha_v12[6] = {
       -1.393236997981977,   -0.231152406188676,    0.155351465642094,
       -6.21836623696556e-3, -9.190829861057113e-5, 6.275259577497896e-7
    };
    static const double beta_v12[6] = {
	1.011728045548831,   -0.751971469674635,    1.255772699323164e-2,
	1.0022008145159e-2,  -2.42068134815573e-4,  5.008480613664573e-7
    };
    static const double  y0_v12 = 1.50;
    double yp, xp, xm, sum, yp2, xp2, xm2;
    int k;

    sum = 0.;
    yp = y + y0_v12;
    yp2 = yp * yp;
    if((y > 0.85) || (fabs(x) < (18.1 * y + 1.65))) {
	/* Bereich I */
	for(k=0; k<6; k++) {
	    xp = x + T_v12[k];
	    xm = x - T_v12[k];
	    sum += ((alpha_v12[k] * xm + beta_v12[k] * yp) / (xm * xm + yp2)
                    + (beta_v12[k] * yp - alpha_v12[k] * xp) / (xp * xp + yp2));
	}
    } else {
	/* Bereich II */
	for(k=0; k<6; k++) {
	    xp = x + T_v12[k];
	    xp2 = xp * xp;
	    xm = x - T_v12[k];
	    xm2 = xm * xm;
	    sum += (((beta_v12[k] * (xm2 - y0_v12 * yp) - alpha_v12[k] * xm * (yp + y0_v12))
                     / ((xm2 + yp2) * (xm2 + y0_v12 * y0_v12)))
                    + ((beta_v12[k] * (xp2 - y0_v12 * yp) + alpha_v12[k] * xp * (yp + y0_v12))
                       / ((xp2 + yp2) * (xp2 + y0_v12 * y0_v12))));
	}
	if(fabs(x) < 100.)
	    sum = y * sum + exp(-pow(x, 2));
	else
	    sum *= y;
    }
    return sum;
}


static void voigt(size_t n, double *x, double *y, double w[2], double xc) 
    /* Evaluate normalized Voigt profile at x around xc with Gaussian
     * linewidth contribution w[0] and Lorentzian linewidth
     * contribution w[1].
     */
{
    /* Transform into reduced coordinates and call Humlicek's 12 point
     * formula:
     *     x = 2 \sqrt{\ln2} \frac{\nu-\nu_0}{\Delta\nu_G}
     *     y = \sqrt{\ln2} \frac{\Delta\nu_L}{\Delta\nu_G}
     */
    int i;
    double yh = sqrt(M_LN2) * w[1] / w[0];
    for(i=0; i<n; i++) {
        double xh = 2. * sqrt(M_LN2) * (x[i]-xc) / w[0];
        y[i] = 2.*sqrt(M_LN2/M_PI)/w[0] * humlicek_v12(xh, yh);
    }
}



/*** Python interface ***/

static PyObject *_Error;


static PyObject * 
_lineshape_gauss(PyObject *self, PyObject *args, PyObject *keywds) 
{
    int f;
    double w, xc = 0.0;
    static char *kwlist[] = {"x", "w", "xc", "y", NULL};
    PyObject *ox, *oy=Py_None;
    PyArrayObject *x, *y;

    if(! PyArg_ParseTupleAndKeywords(args, keywds, "Od|dO", kwlist,
                                     &ox, &w, &xc, &oy))
        return PyErr_Format(PyExc_RuntimeError,  "gauss: invalid parameters");

    if((f = PyFloat_Check(ox)) || PyInt_Check(ox)) {
        /* scalar arguments -- always *return* Float result */
        double xa[1], ya[1];
        if(f)
            xa[0] = PyFloat_AS_DOUBLE(ox);
        else
            xa[0] = (double)PyInt_AS_LONG(ox);
        Py_BEGIN_ALLOW_THREADS;
        gauss(1, xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
        Py_DECREF(ox);
        return PyFloat_FromDouble(ya[0]);
    } else {
        /* array conversion */        
        if(! ((x = NA_InputArray(ox, tFloat64, C_ARRAY))
              && (y = NA_OptionalOutputArray(oy, tFloat64, C_ARRAY, x))))
            return 0;
        if(x->nd != 1)
            return PyErr_Format(_Error, "gauss: x must be scalar or 1d array.");
        if (!NA_ShapeEqual(x, y))
            return PyErr_Format(_Error, "gauss: x and y numarray must have same length.");

        /* calculate profile */
	{
        double *xa = NA_OFFSETDATA(x);
        double *ya = NA_OFFSETDATA(y);
        Py_BEGIN_ALLOW_THREADS;
        gauss(x->dimensions[0], xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
	}
    
        /* cleanup and return */
        Py_XDECREF(x);
        return NA_ReturnOutput(oy, y);
    }
}



static PyObject * 
_lineshape_lorentz(PyObject *self, PyObject *args, PyObject *keywds) 
{
    int f;
    double w, xc = 0.0;
    static char *kwlist[] = {"x", "w", "xc", "y", NULL};
    PyObject *ox, *oy=Py_None;
    PyArrayObject *x, *y;

    if(! PyArg_ParseTupleAndKeywords(args, keywds, "Od|dO", kwlist,
                                     &ox, &w, &xc, &oy))
        return PyErr_Format(PyExc_RuntimeError,  "lorentz: invalid parameters");

    if((f = PyFloat_Check(ox)) || PyInt_Check(ox)) {
        /* scalar arguments -- always *return* Float result */
        double xa[1], ya[1];
        if(f)
            xa[0] = PyFloat_AS_DOUBLE(ox);
        else
            xa[0] = (double)PyInt_AS_LONG(ox);
        Py_BEGIN_ALLOW_THREADS;
        lorentz(1, xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
        Py_DECREF(ox);
        return PyFloat_FromDouble(ya[0]);
    } else {
        /* array conversion */        
        if(! ((x = NA_InputArray(ox, tFloat64, C_ARRAY))
              && (y = NA_OptionalOutputArray(oy, tFloat64, C_ARRAY, x))))
            return 0;
        if(x->nd != 1)
            return PyErr_Format(_Error, "lorentz: x must be scalar or 1d array.");
        if (!NA_ShapeEqual(x, y))
            return PyErr_Format(_Error, "lorentz: x and y numarray must have same length.");

        /* calculate profile */
	{
        double *xa = NA_OFFSETDATA(x);
        double *ya = NA_OFFSETDATA(y);

        Py_BEGIN_ALLOW_THREADS;
        lorentz(x->dimensions[0], xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
	}

        /* cleanup and return */
        Py_XDECREF(x);
        return NA_ReturnOutput(oy, y);
    }
}



static PyObject * 
_lineshape_voigt(PyObject *self, PyObject *args, PyObject *keywds) 
{
    int f;
    double w[2], xc = 0.0;
    static char *kwlist[] = {"x", "w", "xc", "y", NULL};
    PyObject *wt, *ox, *oy=Py_None;
    PyArrayObject *x, *y;

    if(! PyArg_ParseTupleAndKeywords(args, keywds, "OO|dO", kwlist,
                                     &ox, &wt, &xc, &oy))
        return PyErr_Format(PyExc_RuntimeError,  "voigt: invalid parameters");

    /* parse linewidths tuple */
    if(! PyArg_ParseTuple(wt, "dd", &(w[0]), &(w[1])))
        return(0);

    if((f = PyFloat_Check(ox)) || PyInt_Check(ox)) {
        /* scalar arguments -- always *return* Float result */
        double xa[1], ya[1];
        if(f)
            xa[0] = PyFloat_AS_DOUBLE(ox);
        else
            xa[0] = (double)PyInt_AS_LONG(ox);
        Py_BEGIN_ALLOW_THREADS;
        voigt(1, xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
        Py_DECREF(ox);
        return PyFloat_FromDouble(ya[0]);
    } else {
        /* array conversion */        
        if(! ((x = NA_InputArray(ox, tFloat64, C_ARRAY))
              && (y = NA_OptionalOutputArray(oy, tFloat64, C_ARRAY, x))))
            return 0;
        if(x->nd != 1)
            return PyErr_Format(_Error, "voigt: x must be scalar or 1d array.");
        if (!NA_ShapeEqual(x, y))
            return PyErr_Format(_Error, "voigt: x and y numarray must have same length.");

        /* calculate profile */
	{
        double *xa = NA_OFFSETDATA(x);
        double *ya = NA_OFFSETDATA(y);
        Py_BEGIN_ALLOW_THREADS;
        voigt(x->dimensions[0], xa, ya, w, xc);
        Py_END_ALLOW_THREADS;
	}

        /* cleanup and return */
        Py_XDECREF(x);
        return NA_ReturnOutput(oy, y);
    }
}




/*** table of methods ***/

static PyMethodDef _lineshape_Methods[] = {
    {"gauss", (PyCFunction)_lineshape_gauss, METH_VARARGS|METH_KEYWORDS,
     "gauss(x, w, xc=0.0, y=None)\n\n"
     "Gaussian lineshape function\n\n" \
     "Calculate normalized Gaussian with full-width at half maximum |w| at |x|,\n" \
     "optionally specifying the line-center |xc|.\n" \
     "If, and only if |x| is an array an optional output array |y| can be\n" \
     "specified.  In this case |x| and |y| must be one-dimensional numarray\n" \
     "with identical shapes.\n\n" \
     "If |x| is an scalar the routine always gives the result as scalar\n" \
     "return value."
    },
    {"lorentz", (PyCFunction)_lineshape_lorentz, METH_VARARGS|METH_KEYWORDS,
     "lorentz(x, w, xc=0.0, y=None)\n\n"
     "Lorentzian lineshape function\n\n" \
     "Calculate normalized Lorentzian with full-width at half maximum |w| at |x|,\n" \
     "optionally specifying the line-center |xc|.\n" \
     "If, and only if |x| is an array an optional output array |y| can be\n" \
     "specified.  In this case |x| and |y| must be one-dimensional numarray\n" \
     "with identical shapes.\n\n" \
     "If |x| is an scalar the routine always gives the result as scalar\n" \
     "return value."
    },
    {"voigt", (PyCFunction)_lineshape_voigt, METH_VARARGS|METH_KEYWORDS,
     "voigt(x, w, xc=0.0, y=None)\n\n"
     "Voigt-lineshape function\n\n" \
     "Calculate normalized Voigt-profile with Gaussian full-width at half maximum |w[0]| and\n" \
     "Lorentzian full-width at half maximum |w[1]| at |x|, optionally specifying the line-center\n" \
     "|xc|.\n" \
     "If, and only if |x| is an array an optional output array |y| can be\n" \
     "specified.  In this case |x| and |y| must be one-dimensional numarray\n" \
     "with identical shapes.\n\n" \
     "If |x| is an scalar the routine always gives the result as scalar\n" \
     "return value.\n\n" \
     "This function uses Humlicek's 12-point formula to approximate the Voigt\n" \
     "profile (J. Humlicek, J. Quant. Spectrosc. Radiat. Transfer, 21, 309 (1978))."
    },
    {NULL, NULL, 0, ""}
};




/*** module initialization ***/

PyMODINIT_FUNC init_lineshape(void)
{
    PyObject *m, *d;
    m = Py_InitModule("_lineshape", _lineshape_Methods);
    d = PyModule_GetDict(m);
    _Error = PyErr_NewException("_lineshape.error", NULL, NULL);
    PyDict_SetItemString(d, "error", _Error);
    import_libnumarray();
}



/*
 * Local Variables:
 * mode: c
 * c-file-style: "Stroustrup"
 * fill-column: 80
 * End:
 */
