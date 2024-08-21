#include <Python.h>
#include "numpy/arrayobject.h"

#define PyInt_AsLong PyLong_AsLong

static PyObject *fitpack_error;
#include "__fitpack.h"

#ifdef HAVE_ILP64

#define F_INT npy_int64
#define F_INT_NPY NPY_INT64
#define F_INT_MAX NPY_MAX_INT64

#if NPY_BITSOF_SHORT == 64
#define F_INT_PYFMT   "h"
#elif NPY_BITSOF_INT == 64
#define F_INT_PYFMT   "i"
#elif NPY_BITSOF_LONG == 64
#define F_INT_PYFMT   "l"
#elif NPY_BITSOF_LONGLONG == 64
#define F_INT_PYFMT   "L"
#else
#error No compatible 64-bit integer size. \
       Please contact NumPy maintainers and give detailed information about your \
       compiler and platform, or set NPY_USE_BLAS64_=0
#endif

#else

#define F_INT int
#define F_INT_NPY NPY_INT
#define F_INT_MAX NPY_MAX_INT32
#if NPY_BITSOF_SHORT == 32
#define F_INT_PYFMT   "h"
#elif NPY_BITSOF_INT == 32
#define F_INT_PYFMT   "i"
#elif NPY_BITSOF_LONG == 32
#define F_INT_PYFMT   "l"
#else
#error No compatible 32-bit integer size. \
       Please contact NumPy maintainers and give detailed information about your \
       compiler and platform
#endif

#endif


/*
 * Multiply mx*my, check for integer overflow.
 * Inputs are Fortran ints, and the output is npy_intp, for use
 * in PyArray_SimpleNew et al.
 * Return -1 on overflow.
 */
npy_intp
_mul_overflow_intp(F_INT mx, F_INT my) {
    npy_intp int_max, mxy;

    /* Limit is min of (largest array size, max of Fortran int) */
    int_max = (F_INT_MAX < NPY_MAX_INTP) ? F_INT_MAX : NPY_MAX_INTP;
    /* v = int_max/my is largest integer multiple of `my` such that
       v * my <= int_max
    */
    if (my != 0 && int_max/my < mx) {
        /* Integer overflow */
        PyErr_Format(PyExc_RuntimeError,
                     "Cannot produce output of size %dx%d (size too large)",
                     mx, my);
        return -1;
    }
    mxy = (npy_intp)mx * (npy_intp)my;
    return mxy;
}


/*
 * Multiply mx*my, check for integer overflow, where both inputs and
 * the output are Fortran ints.
 * Return -1 on overflow.
 */
F_INT
_mul_overflow_f_int(F_INT mx, F_INT my) {
    F_INT mxy;

    /* v = int_max/my is largest integer multiple of `my` such that
       v * my <= F_INT_MAX
    */
    if (my != 0 && F_INT_MAX/my < mx) {
        /* Integer overflow */
        PyErr_Format(PyExc_RuntimeError,
                     "Cannot produce output of size %dx%d (size too large)",
                     mx, my);
        return -1;
    }
    mxy = mx * my;
    return mxy;
}


/*
 * Functions moved verbatim from __fitpack.h
 */


/*
 * Python-C wrapper of FITPACK (by P. Dierckx) (in netlib known as dierckx)
 * Author: Pearu Peterson <pearu@ioc.ee>
 * June 1.-4., 1999
 * June 7. 1999
 * $Revision$
 * $Date$
 */

/*  module_methods:
 * {"_parcur", fitpack_parcur, METH_VARARGS, doc_parcur},
 * {"_surfit", fitpack_surfit, METH_VARARGS, doc_surfit},
 * {"_insert", fitpack_insert, METH_VARARGS, doc_insert},
 */

/* link libraries: (one item per line)
   ddierckx
 */
/* python files: (to be imported to Multipack.py)
   fitpack.py
 */

#if defined(UPPERCASE_FORTRAN)
	#if defined(NO_APPEND_FORTRAN)
	/* nothing to do */
	#else
		#define PARCUR PARCUR_
		#define CLOCUR CLOCUR_
		#define SURFIT SURFIT_
		#define INSERT INSERT_
	#endif
#else
	#if defined(NO_APPEND_FORTRAN)
		#define PARCUR parcur
		#define CLOCUR clocur
		#define SURFIT surfit
		#define INSERT insert
	#else
		#define PARCUR parcur_
		#define CLOCUR clocur_
		#define SURFIT surfit_
		#define INSERT insert_
	#endif
#endif

void PARCUR(F_INT*,F_INT*,F_INT*,F_INT*,double*,F_INT*,double*,
        double*,double*,double*,F_INT*,double*,F_INT*,F_INT*,
        double*,F_INT*,double*,double*,double*,F_INT*,F_INT*,F_INT*);
void CLOCUR(F_INT*,F_INT*,F_INT*,F_INT*,double*,F_INT*,double*,
        double*,F_INT*,double*,F_INT*,F_INT*,double*,F_INT*,
        double*,double*,double*,F_INT*,F_INT*,F_INT*);
void SURFIT(F_INT*,F_INT*,double*,double*,double*,double*,
        double*,double*,double*,double*,F_INT*,F_INT*,double*,
        F_INT*,F_INT*,F_INT*,double*,F_INT*,double*,F_INT*,double*,
        double*,double*,double*,F_INT*,double*,F_INT*,F_INT*,F_INT*,F_INT*);
void INSERT(F_INT*,double*,F_INT*,double*,F_INT*,double*,double*,
        F_INT*,double*,F_INT*,F_INT*);

/* Note that curev, cualde need no interface. */

static char doc_surfit[] = " [tx,ty,c,o] = _surfit(x, y, z, w, xb, xe, yb, ye,"\
      " kx,ky,iopt,s,eps,tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)";
static PyObject *
fitpack_surfit(PyObject *dummy, PyObject *args)
{
    F_INT iopt, m, kx, ky, nxest, nyest, lwrk1, lwrk2, *iwrk, kwrk, ier;
    F_INT lwa, nxo, nyo, i, lcest, nmax, nx, ny, lc;
    npy_intp dims[1], lc_intp;
    double *x, *y, *z, *w, xb, xe, yb, ye, s, *tx, *ty, *c, fp;
    double *wrk1, *wrk2, *wa = NULL, eps;
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z, *ap_w = NULL;
    PyArrayObject *ap_tx = NULL,*ap_ty = NULL,*ap_c = NULL, *ap_wrk = NULL;
    PyObject *x_py = NULL, *y_py = NULL, *z_py = NULL, *w_py = NULL;
    PyObject *tx_py = NULL, *ty_py = NULL, *wrk_py = NULL;

    nx = ny = ier = nxo = nyo = 0;
    if (!PyArg_ParseTuple(args, ("OOOOdddd" F_INT_PYFMT F_INT_PYFMT F_INT_PYFMT
                                 "ddOO" F_INT_PYFMT F_INT_PYFMT "O"
                                 F_INT_PYFMT F_INT_PYFMT),
                &x_py, &y_py, &z_py, &w_py, &xb, &xe, &yb, &ye,
                &kx, &ky, &iopt, &s, &eps, &tx_py, &ty_py, &nxest,
                &nyest, &wrk_py, &lwrk1, &lwrk2)) {
        return NULL;
    }
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_DOUBLE, 0, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, NPY_DOUBLE, 0, 1);
    ap_z = (PyArrayObject *)PyArray_ContiguousFromObject(z_py, NPY_DOUBLE, 0, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, NPY_DOUBLE, 0, 1);
    ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, NPY_DOUBLE, 0, 1);
    /*ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, F_INT_NPY, 0, 1);*/
    if (ap_x == NULL
            || ap_y == NULL
            || ap_z == NULL
            || ap_w == NULL
            || ap_wrk == NULL) {
        goto fail;
    }
    x = (double *) PyArray_DATA(ap_x);
    y = (double *) PyArray_DATA(ap_y);
    z = (double *) PyArray_DATA(ap_z);
    w = (double *) PyArray_DATA(ap_w);
    m = PyArray_DIMS(ap_x)[0];
    nmax = nxest;
    if (nmax < nyest) {
        nmax = nyest;
    }
    /* lcest = (nxest - kx - 1)*(nyest - ky - 1); */
    lcest = _mul_overflow_f_int(nxest - kx - 1, nyest - ky - 1);
    if (lcest < 0) {
        PyErr_NoMemory();
        goto fail;
    }
    /* kwrk computation is unlikely to overflow if lcest above did not.*/
    kwrk = m + (nxest - 2*kx - 1)*(nyest - 2*ky - 1);
    lwa = 2*nmax + lcest + lwrk1 + lwrk2 + kwrk;
    if ((wa = malloc(lwa*sizeof(double))) == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    /*
     * NOTE: the work arrays MUST be aligned on double boundaries, as Fortran
     *       compilers (e.g. gfortran) may assume that.  Malloc gives suitable
     *       alignment, so we are here careful not to mess it up.
     */
    tx = wa;
    ty = tx + nmax;
    c = ty + nmax;
    wrk1 = c + lcest;
    iwrk = (F_INT *)(wrk1 + lwrk1);
    wrk2 = ((double *)iwrk) + kwrk;
    if (iopt) {
        ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, NPY_DOUBLE, 0, 1);
        ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, NPY_DOUBLE, 0, 1);
        if (ap_tx == NULL || ap_ty == NULL) {
            goto fail;
        }
        nx = nxo = PyArray_DIMS(ap_tx)[0];
        ny = nyo = PyArray_DIMS(ap_ty)[0];
        memcpy(tx, PyArray_DATA(ap_tx), nx*sizeof(double));
        memcpy(ty, PyArray_DATA(ap_ty), ny*sizeof(double));
    }
    if (iopt==1) {
        /* lc = (nx - kx - 1)*(ny - ky - 1); */
        lc = _mul_overflow_f_int(nx - kx - 1, ny - ky -1);
        if (lc < 0) {
            goto fail;
        }

        memcpy(wrk1, PyArray_DATA(ap_wrk), lc*sizeof(double));
        /*memcpy(iwrk,PyArray_DATA(ap_iwrk),n*sizeof(int));*/
    }
    SURFIT(&iopt, &m, x, y, z, w, &xb, &xe, &yb, &ye, &kx, &ky,
            &s, &nxest, &nyest, &nmax, &eps, &nx, tx, &ny, ty,
            c, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
    i = 0;
    while ((ier > 10) && (i++ < 5)) {
        lwrk2 = ier;
        if ((wrk2 = malloc(lwrk2*sizeof(double))) == NULL) {
            PyErr_NoMemory();
            goto fail;
        }
        SURFIT(&iopt, &m, x, y, z, w, &xb, &xe, &yb, &ye, &kx, &ky,
                &s, &nxest, &nyest, &nmax, &eps, &nx, tx, &ny, ty,
                c, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
        free(wrk2);
    }
    if (ier == 10) {
        PyErr_SetString(PyExc_ValueError, "Invalid inputs.");
        goto fail;
    }

    lc_intp = _mul_overflow_intp(nx - kx - 1, ny - ky -1);
    if (lc_intp < 0) {
        goto fail;
    }

    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    dims[0] = nx;
    ap_tx = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    dims[0] = ny;
    ap_ty = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    dims[0] = lc_intp;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (ap_tx == NULL
            || ap_ty == NULL
            || ap_c == NULL) {
        goto fail;
    }
    if ((iopt == 0)||(nx > nxo)||(ny > nyo)) {
        Py_XDECREF(ap_wrk);
        dims[0] = lc_intp;
        ap_wrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (ap_wrk == NULL) {
            goto fail;
        }
        /*ap_iwrk = (PyArrayObject *)PyArray_SimpleNew(1,&n,F_INT_NPY);*/
    }
    if (PyArray_DIMS(ap_wrk)[0] < lc_intp) {
        Py_XDECREF(ap_wrk);
        dims[0] = lc_intp;
        ap_wrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (ap_wrk == NULL) {
            goto fail;
        }
    }
    memcpy(PyArray_DATA(ap_tx), tx, nx*sizeof(double));
    memcpy(PyArray_DATA(ap_ty), ty, ny*sizeof(double));
    memcpy(PyArray_DATA(ap_c), c, lc_intp*sizeof(double));
    memcpy(PyArray_DATA(ap_wrk), wrk1, lc_intp*sizeof(double));
    /*memcpy(PyArray_DATA(ap_iwrk),iwrk,n*sizeof(int));*/
    free(wa);
    Py_DECREF(ap_x);
    Py_DECREF(ap_y);
    Py_DECREF(ap_z);
    Py_DECREF(ap_w);
    return Py_BuildValue("NNN{s:N,s:i,s:d}",PyArray_Return(ap_tx),
            PyArray_Return(ap_ty),PyArray_Return(ap_c),
            "wrk",PyArray_Return(ap_wrk),
            "ier",ier,"fp",fp);

fail:
    free(wa);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_wrk);
    /*Py_XDECREF(ap_iwrk);*/
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError,
                "An error occurred.");
    }
    return NULL;
}


static char doc_parcur[] = " [t,c,o] = _parcur(x,w,u,ub,ue,k,iopt,ipar,s,t,nest,wrk,iwrk,per)";
static PyObject *
fitpack_parcur(PyObject *dummy, PyObject *args)
{
    F_INT k, iopt, ipar, nest, *iwrk, idim, m, mx, no=0, nc, ier, lwa, lwrk, i, per;
    F_INT n=0,  lc;
    npy_intp dims[1];
    double *x, *w, *u, *c, *t, *wrk, *wa=NULL, ub, ue, fp, s;
    PyObject *x_py = NULL, *u_py = NULL, *w_py = NULL, *t_py = NULL;
    PyObject *wrk_py=NULL, *iwrk_py=NULL;
    PyArrayObject *ap_x = NULL, *ap_u = NULL, *ap_w = NULL, *ap_t = NULL, *ap_c = NULL;
    PyArrayObject *ap_wrk = NULL, *ap_iwrk = NULL;

    if (!PyArg_ParseTuple(args, ("OOOdd" F_INT_PYFMT F_INT_PYFMT F_INT_PYFMT
                                 "dO" F_INT_PYFMT "OO" F_INT_PYFMT),
                          &x_py, &w_py, &u_py, &ub, &ue, &k, &iopt, &ipar,
                          &s, &t_py, &nest, &wrk_py, &iwrk_py, &per)) {
        return NULL;
    }
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_DOUBLE, 0, 1);
    ap_u = (PyArrayObject *)PyArray_ContiguousFromObject(u_py, NPY_DOUBLE, 0, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, NPY_DOUBLE, 0, 1);
    ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, NPY_DOUBLE, 0, 1);
    ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, F_INT_NPY, 0, 1);
    if (ap_x == NULL
            || ap_u == NULL
            || ap_w == NULL
            || ap_wrk == NULL
            || ap_iwrk == NULL) {
        goto fail;
    }
    x = (double *) PyArray_DATA(ap_x);
    u = (double *) PyArray_DATA(ap_u);
    w = (double *) PyArray_DATA(ap_w);
    m = PyArray_DIMS(ap_w)[0];
    mx = PyArray_DIMS(ap_x)[0];
    idim = mx/m;
    if (per) {
        lwrk = m*(k + 1) + nest*(7 + idim + 5*k);
    }
    else {
        lwrk = m*(k + 1) + nest*(6 + idim + 3*k);
    }
    nc = idim*nest;
    lwa = nc + 2*nest + lwrk;
    if ((wa = malloc(lwa*sizeof(double))) == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    t = wa;
    c = t + nest;
    wrk = c + nc;
    iwrk = (F_INT *)(wrk + lwrk);
    if (iopt) {
        ap_t=(PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_DOUBLE, 0, 1);
        if (ap_t == NULL) {
            goto fail;
        }
        n = no = PyArray_DIMS(ap_t)[0];
        memcpy(t, PyArray_DATA(ap_t), n*sizeof(double));
        Py_DECREF(ap_t);
        ap_t = NULL;
    }
    if (iopt == 1) {
        memcpy(wrk, PyArray_DATA(ap_wrk), n*sizeof(double));
        memcpy(iwrk, PyArray_DATA(ap_iwrk), n*sizeof(F_INT));
    }
    if (per) {
        CLOCUR(&iopt, &ipar, &idim, &m, u, &mx, x, w, &k, &s, &nest,
                &n, t, &nc, c, &fp, wrk, &lwrk, iwrk, &ier);
    }
    else {
        PARCUR(&iopt, &ipar, &idim, &m, u, &mx, x, w, &ub, &ue, &k,
                &s, &nest, &n, t, &nc, c, &fp, wrk, &lwrk, iwrk, &ier);
    }
    if (ier == 10) {
        PyErr_SetString(PyExc_ValueError, "Invalid inputs.");
        goto fail;
    }
    if (ier > 0 && n == 0) {
        n = 1;
    }
    lc = (n - k - 1)*idim;
    dims[0] = n;
    ap_t = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    dims[0] = lc;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (ap_t == NULL || ap_c == NULL) {
        goto fail;
    }
    if (iopt != 1|| n > no) {
        Py_XDECREF(ap_wrk);
        ap_wrk = NULL;
        Py_XDECREF(ap_iwrk);
        ap_iwrk = NULL;

        dims[0] = n;
        ap_wrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (ap_wrk == NULL) {
            goto fail;
        }
        ap_iwrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, F_INT_NPY);
        if (ap_iwrk == NULL) {
            goto fail;
        }
    }
    memcpy(PyArray_DATA(ap_t), t, n*sizeof(double));
    for (i = 0; i < idim; i++)
        memcpy((double *)PyArray_DATA(ap_c) + i*(n - k - 1), c + i*n, (n - k - 1)*sizeof(double));
    memcpy(PyArray_DATA(ap_wrk), wrk, n*sizeof(double));
    memcpy(PyArray_DATA(ap_iwrk), iwrk, n*sizeof(F_INT));
    free(wa);
    Py_DECREF(ap_x);
    Py_DECREF(ap_w);
    return Py_BuildValue(("NN{s:N,s:d,s:d,s:N,s:N,s:" F_INT_PYFMT ",s:d}"), PyArray_Return(ap_t),
            PyArray_Return(ap_c), "u", PyArray_Return(ap_u), "ub", ub, "ue", ue,
            "wrk", PyArray_Return(ap_wrk), "iwrk", PyArray_Return(ap_iwrk),
            "ier", ier, "fp",fp);
fail:
    free(wa);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_u);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;
}

static char doc_insert[] = " [tt,cc,ier] = _insert(iopt,t,c,k,x,m)";
static PyObject *
fitpack_insert(PyObject *dummy, PyObject*args)
{
    F_INT iopt, n, nn, k, ier, m, nest;
    npy_intp dims[1];
    double x;
    double *t_in, *c_in, *t_out, *c_out, *t_buf = NULL, *c_buf = NULL, *p;
    double *t1, *t2, *c1, *c2;
    PyArrayObject *ap_t_in = NULL, *ap_c_in = NULL, *ap_t_out = NULL, *ap_c_out = NULL;
    PyObject *t_py = NULL, *c_py = NULL;
    PyObject *ret = NULL;

    if (!PyArg_ParseTuple(args,(F_INT_PYFMT "OO" F_INT_PYFMT "d" F_INT_PYFMT),
                          &iopt,&t_py,&c_py,&k, &x, &m)) {
        return NULL;
    }
    ap_t_in = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_DOUBLE, 0, 1);
    ap_c_in = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_DOUBLE, 0, 1);
    if (ap_t_in == NULL || ap_c_in == NULL) {
        goto fail;
    }
    t_in = (double *)PyArray_DATA(ap_t_in);
    c_in = (double *)PyArray_DATA(ap_c_in);
    n = PyArray_DIMS(ap_t_in)[0];
    nest = n + m;
    dims[0] = nest;
    ap_t_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    ap_c_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (ap_t_out == NULL || ap_c_out == NULL) {
        goto fail;
    }
    t_out = (double *)PyArray_DATA(ap_t_out);
    c_out = (double *)PyArray_DATA(ap_c_out);

    /*
     * Call the INSERT routine m times to insert m-multiplicity knot, ie.:
     *
     *     for _ in range(n, nest):
     *         t, c = INSERT(t, c)
     *     return t, c
     *
     * We need to ensure that input and output buffers given to INSERT routine
     * do not point to same memory, which is not allowed by Fortran. For this,
     * we use temporary storage, and cycle between it and the output buffers.
     */
    t2 = t_in;
    c2 = c_in;
    t1 = t_out;
    c1 = c_out;

    for ( ; n < nest; n++) {
        /* Swap buffers */
        p = t2; t2 = t1; t1 = p;
        p = c2; c2 = c1; c1 = p;

        /* Allocate temporary buffer (needed for m > 1) */
        if (t2 == t_in) {
            if (t_buf == NULL) {
                t_buf = calloc(nest, sizeof(double));
                c_buf = calloc(nest, sizeof(double));
                if (t_buf == NULL || c_buf == NULL) {
                    PyErr_NoMemory();
                    goto fail;
                }
            }
            t2 = t_buf;
            c2 = c_buf;
        }

        INSERT(&iopt, t1, &n, c1, &k, &x, t2, &nn, c2, &nest, &ier);

        if (ier) {
            break;
        }
    }

    /* Ensure output ends up in output buffers */
    if (t2 != t_out) {
        memcpy(t_out, t2, nest * sizeof(double));
        memcpy(c_out, c2, nest * sizeof(double));
    }

    Py_DECREF(ap_c_in);
    Py_DECREF(ap_t_in);
    free(t_buf);
    free(c_buf);
    ret = Py_BuildValue(("NN" F_INT_PYFMT), PyArray_Return(ap_t_out), PyArray_Return(ap_c_out), ier);
    return ret;

fail:
    Py_XDECREF(ap_c_out);
    Py_XDECREF(ap_t_out);
    Py_XDECREF(ap_c_in);
    Py_XDECREF(ap_t_in);
    free(t_buf);
    free(c_buf);
    return NULL;
}

static struct PyMethodDef fitpack_module_methods[] = {
{"_parcur",
    fitpack_parcur,
    METH_VARARGS, doc_parcur},
{"_surfit",
    fitpack_surfit,
    METH_VARARGS, doc_surfit},
{"_insert",
    fitpack_insert,
    METH_VARARGS, doc_insert},
{NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fitpack",
    NULL,
    -1,
    fitpack_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__fitpack(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);

    fitpack_error = PyErr_NewException ("_fitpack.error", NULL, NULL);
    if (fitpack_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", fitpack_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}
