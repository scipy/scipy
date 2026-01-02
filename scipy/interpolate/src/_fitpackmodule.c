#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/arrayobject.h"
#include "dfitpack.h"

#define PyInt_AsLong PyLong_AsLong

static PyObject *fitpack_error;


static char doc_bispeu[] = " z,ier = bispeu(tx,ty,c,kx,ky,x,y)";
static PyObject *
fitpack_bispeu(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int nx, ny, kx, ky, m, ier, lwrk;
    npy_intp dims[1];
    double *tx, *ty, *c, *x, *y, *z, *wrk;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z = NULL;
    PyObject *tx_py = NULL, *ty_py = NULL, *c_py = NULL;
    PyObject *x_py = NULL, *y_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOOiiOO"),
                          &tx_py, &ty_py, &c_py, &kx, &ky, &x_py, &y_py)) {
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, NPY_FLOAT64, 1, 1);

    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL ||
        ap_x == NULL || ap_y == NULL) {
        goto fail;
    }

    tx = (double *)PyArray_DATA(ap_tx);
    ty = (double *)PyArray_DATA(ap_ty);
    c = (double *)PyArray_DATA(ap_c);
    x = (double *)PyArray_DATA(ap_x);
    y = (double *)PyArray_DATA(ap_y);

    nx = (int)PyArray_DIMS(ap_tx)[0];
    ny = (int)PyArray_DIMS(ap_ty)[0];
    m = (int)PyArray_DIMS(ap_x)[0];

    /* Check dimensions match pyf specification */
    if (PyArray_DIMS(ap_y)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x and y must have same length");
        goto fail;
    }
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError,
            "c must have length (nx-kx-1)*(ny-ky-1)");
        goto fail;
    }

    /* Allocate output array */
    dims[0] = m;
    ap_z = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_z == NULL) {
        goto fail;
    }
    z = (double *)PyArray_DATA(ap_z);

    /* Allocate work array: lwrk = kx + ky + 2 */
    lwrk = kx + ky + 2;
    wrk = (double *)malloc(lwrk * sizeof(double));
    if (wrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Call the C function */
    bispeu(tx, nx, ty, ny, c, kx, ky, x, y, z, m, wrk, lwrk, &ier);

    free(wrk);
    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);
    Py_DECREF(ap_y);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_z), ier);

fail:
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    return NULL;
}

static char doc_bispev[] = " z,ier = bispev(tx,ty,c,kx,ky,x,y)";
static PyObject *
fitpack_bispev(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int nx, ny, kx, ky, mx, my, ier, lwrk, kwrk;
    npy_intp dims[2];
    double *tx, *ty, *c, *x, *y, *z, *wrk;
    int *iwrk;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z = NULL;
    PyObject *tx_py = NULL, *ty_py = NULL, *c_py = NULL;
    PyObject *x_py = NULL, *y_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOOiiOO"),
                          &tx_py, &ty_py, &c_py, &kx, &ky, &x_py, &y_py)) {
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, NPY_FLOAT64, 1, 1);

    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL || ap_x == NULL || ap_y == NULL) {
        goto fail;
    }

    tx = (double *)PyArray_DATA(ap_tx);
    ty = (double *)PyArray_DATA(ap_ty);
    c = (double *)PyArray_DATA(ap_c);
    x = (double *)PyArray_DATA(ap_x);
    y = (double *)PyArray_DATA(ap_y);

    nx = (int)PyArray_DIMS(ap_tx)[0];
    ny = (int)PyArray_DIMS(ap_ty)[0];
    mx = (int)PyArray_DIMS(ap_x)[0];
    my = (int)PyArray_DIMS(ap_y)[0];

    /* Check dimensions match pyf specification */
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError,
            "c must have length (nx-kx-1)*(ny-ky-1)");
        goto fail;
    }

    /* Allocate output array z(mx, my) - note Fortran order */
    dims[0] = mx;
    dims[1] = my;
    ap_z = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (ap_z == NULL) {
        goto fail;
    }
    z = (double *)PyArray_DATA(ap_z);

    /* Allocate work arrays */
    lwrk = mx * (kx + 1) + my * (ky + 1);
    kwrk = mx + my;
    wrk = (double *)malloc(lwrk * sizeof(double));
    iwrk = (int *)malloc(kwrk * sizeof(int));
    if (wrk == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        free(wrk);
        free(iwrk);
        goto fail;
    }

    /* Call the C function */
    bispev(tx, nx, ty, ny, c, kx, ky, x, mx, y, my, z, wrk, lwrk, iwrk, kwrk, &ier);

    free(wrk);
    free(iwrk);
    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);
    Py_DECREF(ap_y);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_z), ier);

fail:
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    return NULL;
}

static char doc_splev[] = " y,ier = splev(t,c,k,x,[e])";
static PyObject *
fitpack_splev(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int n, nc, k, m, e, ier;
    npy_intp dims[1];
    double *t, *c, *x, *y;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_x = NULL, *ap_y = NULL;
    PyObject *t_py = NULL, *c_py = NULL, *x_py = NULL;

    e = 0;  /* Default value */
    if (!PyArg_ParseTuple(args, ("OOiO|i"),
                          &t_py, &c_py, &k, &x_py, &e)) {
        return NULL;
    }

    /* Validate e parameter */
    if (e < 0 || e > 3) {
        PyErr_SetString(PyExc_ValueError, "e must be between 0 and 3");
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);

    if (ap_t == NULL || ap_c == NULL || ap_x == NULL) {
        goto fail;
    }

    t = (double *)PyArray_DATA(ap_t);
    c = (double *)PyArray_DATA(ap_c);
    x = (double *)PyArray_DATA(ap_x);

    n = (int)PyArray_DIMS(ap_t)[0];
    nc = (int)PyArray_DIMS(ap_c)[0];
    m = (int)PyArray_DIMS(ap_x)[0];

    /* Check c dimension: len(c) >= n-k-1 */
    if (nc < n - k - 1) {
        PyErr_SetString(PyExc_ValueError, "c array is too small");
        goto fail;
    }

    /* Allocate output array */
    dims[0] = m;
    ap_y = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_y == NULL) {
        goto fail;
    }
    y = (double *)PyArray_DATA(ap_y);

    // Call the C function
    splev(t, n, c, nc, k, x, y, m, e, &ier);

    Py_DECREF(ap_t);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_y), ier);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    return NULL;
}

static char doc_splder[] = " dy,ier = splder(t,c,k,nu,x,e)";
static PyObject *
fitpack_splder(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int n, nc, k, nu, m, e, ier;
    npy_intp dims[1];
    double *t, *c, *x, *y, *wrk;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_x = NULL, *ap_y = NULL;
    PyObject *t_py = NULL, *c_py = NULL, *x_py = NULL;

    nu = 1;  /* Default value */
    e = 0;   /* Default value */
    if (!PyArg_ParseTuple(args, ("OOiiO|i"),
                          &t_py, &c_py, &k, &nu, &x_py, &e)) {
        return NULL;
    }

    /* Validate parameters */
    if (nu < 0 || nu > k) {
        PyErr_SetString(PyExc_ValueError, "nu must be between 0 and k");
        return NULL;
    }
    if (e < 0 || e > 3) {
        PyErr_SetString(PyExc_ValueError, "e must be between 0 and 3");
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);

    if (ap_t == NULL || ap_c == NULL || ap_x == NULL) {
        goto fail;
    }

    t = (double *)PyArray_DATA(ap_t);
    c = (double *)PyArray_DATA(ap_c);
    x = (double *)PyArray_DATA(ap_x);

    n = (int)PyArray_DIMS(ap_t)[0];
    nc = (int)PyArray_DIMS(ap_c)[0];
    m = (int)PyArray_DIMS(ap_x)[0];

    /* Check c dimension: len(c) >= n-k-1 */
    if (nc < n - k - 1) {
        PyErr_SetString(PyExc_ValueError, "c array is too small");
        goto fail;
    }

    /* Allocate output and work arrays */
    dims[0] = m;
    ap_y = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_y == NULL) {
        goto fail;
    }
    y = (double *)PyArray_DATA(ap_y);

    wrk = (double *)malloc(n * sizeof(double));
    if (wrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Call the C function */
    splder(t, n, c, nc, k, nu, x, y, m, e, wrk, &ier);

    free(wrk);
    Py_DECREF(ap_t);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_y), ier);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    return NULL;
}

static char doc_splint[] = " aint = splint(t,c,k,a,b)";
static PyObject *
fitpack_splint(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int n, nc, k;
    double a, b, aint;
    double *t, *c, *wrk;
    PyArrayObject *ap_t = NULL, *ap_c = NULL;
    PyObject *t_py = NULL, *c_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOidd"),
                          &t_py, &c_py, &k, &a, &b)) {
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);

    if (ap_t == NULL || ap_c == NULL) {
        goto fail;
    }

    t = (double *)PyArray_DATA(ap_t);
    c = (double *)PyArray_DATA(ap_c);

    n = (int)PyArray_DIMS(ap_t)[0];
    nc = (int)PyArray_DIMS(ap_c)[0];

    /* Check c dimension: len(c) >= n-k-1 */
    if (nc < n - k - 1) {
        PyErr_SetString(PyExc_ValueError, "The length of c must be >=n-k-1");
        goto fail;
    }

    /* Allocate work array */
    wrk = (double *)malloc(n * sizeof(double));
    if (wrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Call the C function */
    aint = splint(t, n, c, nc, k, a, b, wrk);

    free(wrk);
    Py_DECREF(ap_t);
    Py_DECREF(ap_c);

    return Py_BuildValue("d", aint);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    return NULL;
}

static char doc_sproot[] = " zero,m,ier = sproot(t,c,[mest])";
static PyObject *
fitpack_sproot(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int n, nc, mest, m, ier;
    npy_intp dims[1];
    double *t, *c, *zero;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_zero = NULL;
    PyObject *t_py = NULL, *c_py = NULL;

    mest = -1;  /* Will be computed from n if not provided */
    if (!PyArg_ParseTuple(args, ("OO|i"),
                          &t_py, &c_py, &mest)) {
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);

    if (ap_t == NULL || ap_c == NULL) {
        goto fail;
    }

    t = (double *)PyArray_DATA(ap_t);
    c = (double *)PyArray_DATA(ap_c);

    n = (int)PyArray_DIMS(ap_t)[0];
    nc = (int)PyArray_DIMS(ap_c)[0];

    /* Check n >= 8 (required for sproot) */
    if (n < 8) {
        PyErr_SetString(PyExc_ValueError, "t array must have at least 8 elements");
        goto fail;
    }

    /* sproot only works with cubic splines (k=3) */
    /* Check c dimension: len(c) >= n-3-1 = n-4 */
    if (nc < n - 4) {
        PyErr_SetString(PyExc_ValueError, "c array is too small");
        goto fail;
    }

    /* Default mest = 3*(n-7) */
    if (mest < 0) {
        mest = 3 * (n - 7);
    }

    /* Allocate output array */
    dims[0] = mest;
    ap_zero = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_zero == NULL) {
        goto fail;
    }
    zero = (double *)PyArray_DATA(ap_zero);

    /* Call the C function */
    sproot(t, n, c, nc, zero, mest, &m, &ier);

    Py_DECREF(ap_t);
    Py_DECREF(ap_c);

    return Py_BuildValue(("Nii"), PyArray_Return(ap_zero), m, ier);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_zero);
    return NULL;
}

static char doc_spalde[] = " d,ier = spalde(t,c,k1,x)";
static PyObject *
fitpack_spalde(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int n, nc, k1, ier;
    npy_intp dims[1];
    double x;
    double *t, *c, *d;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_d = NULL;
    PyObject *t_py = NULL, *c_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOid"),
                          &t_py, &c_py, &k1, &x)) {
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);

    if (ap_t == NULL || ap_c == NULL) {
        goto fail;
    }

    t = (double *)PyArray_DATA(ap_t);
    c = (double *)PyArray_DATA(ap_c);

    n = (int)PyArray_DIMS(ap_t)[0];
    nc = (int)PyArray_DIMS(ap_c)[0];

    /* Check c dimension: len(c) >= n-k1 */
    if (nc < n - k1) {
        PyErr_SetString(PyExc_ValueError, "c array is too small");
        goto fail;
    }

    /* Allocate output array d of size k1 */
    dims[0] = k1;
    ap_d = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_d == NULL) {
        goto fail;
    }
    d = (double *)PyArray_DATA(ap_d);

    /* Call the C function */
    spalde(t, n, c, nc, k1, x, d, &ier);

    Py_DECREF(ap_t);
    Py_DECREF(ap_c);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_d), ier);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_d);
    return NULL;
}

static char doc_fpchec[] = " ier = fpchec(x,t,k)";
static PyObject *
fitpack_fpchec(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int m, n, k, ier;
    double *x, *t;
    PyArrayObject *ap_x = NULL, *ap_t = NULL;
    PyObject *x_py = NULL, *t_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOi"),
                          &x_py, &t_py, &k)) {
        return NULL;
    }

    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);
    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, NPY_FLOAT64, 1, 1);

    if (ap_x == NULL || ap_t == NULL) {
        goto fail;
    }

    x = (double *)PyArray_DATA(ap_x);
    t = (double *)PyArray_DATA(ap_t);

    m = (int)PyArray_DIMS(ap_x)[0];
    n = (int)PyArray_DIMS(ap_t)[0];

    /* Call the C function */
    fpchec(x, m, t, n, k, &ier);

    Py_DECREF(ap_x);
    Py_DECREF(ap_t);

    return Py_BuildValue("i", ier);

fail:
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_t);
    return NULL;
}

static char doc_parder[] = " z,ier = parder(tx,ty,c,kx,ky,nux,nuy,x,y)";
static PyObject *
fitpack_parder(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int nx, ny, kx, ky, nux, nuy, mx, my, ier, lwrk, kwrk;
    npy_intp dims[2];
    double *tx, *ty, *c, *x, *y, *z, *wrk;
    int *iwrk;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z = NULL;
    PyObject *tx_py = NULL, *ty_py = NULL, *c_py = NULL;
    PyObject *x_py = NULL, *y_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOOiiiiOO"),
                          &tx_py, &ty_py, &c_py, &kx, &ky, &nux, &nuy, &x_py, &y_py)) {
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, NPY_FLOAT64, 1, 1);

    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL ||
        ap_x == NULL || ap_y == NULL) {
        goto fail;
    }

    tx = (double *)PyArray_DATA(ap_tx);
    ty = (double *)PyArray_DATA(ap_ty);
    c = (double *)PyArray_DATA(ap_c);
    x = (double *)PyArray_DATA(ap_x);
    y = (double *)PyArray_DATA(ap_y);

    nx = (int)PyArray_DIMS(ap_tx)[0];
    ny = (int)PyArray_DIMS(ap_ty)[0];
    mx = (int)PyArray_DIMS(ap_x)[0];
    my = (int)PyArray_DIMS(ap_y)[0];

    /* Check dimensions match pyf specification */
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError,
            "c must have length (nx-kx-1)*(ny-ky-1)");
        goto fail;
    }

    /* Allocate output array z(mx, my) */
    dims[0] = mx;
    dims[1] = my;
    ap_z = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (ap_z == NULL) {
        goto fail;
    }
    z = (double *)PyArray_DATA(ap_z);

    /* Allocate work arrays */
    lwrk = (nx * ny) + (kx + 1) * mx + (ky + 1) * my;
    kwrk = mx + my;
    wrk = (double *)malloc(lwrk * sizeof(double));
    iwrk = (int *)malloc(kwrk * sizeof(int));
    if (wrk == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        free(wrk);
        free(iwrk);
        goto fail;
    }

    /* Call the C function */
    parder(tx, nx, ty, ny, c, kx, ky, nux, nuy, x, mx, y, my, z, wrk, lwrk, iwrk, kwrk, &ier);

    free(wrk);
    free(iwrk);
    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);
    Py_DECREF(ap_y);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_z), ier);

fail:
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    return NULL;
}

static char doc_pardtc[] = " newc,ier = pardtc(tx,ty,c,kx,ky,nux,nuy)";
static PyObject *
fitpack_pardtc(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int nx, ny, kx, ky, nux, nuy, ier;
    npy_intp dims[1];
    double *tx, *ty, *c, *newc;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL, *ap_newc = NULL;
    PyObject *tx_py = NULL, *ty_py = NULL, *c_py = NULL;

    if (!PyArg_ParseTuple(args, ("OOOiiii"),
                          &tx_py, &ty_py, &c_py, &kx, &ky, &nux, &nuy)) {
        return NULL;
    }

    /* Validate kx, ky */
    if (kx < 1 || kx > 5) {
        PyErr_SetString(PyExc_ValueError, "kx must be between 1 and 5");
        return NULL;
    }
    if (ky < 1 || ky > 5) {
        PyErr_SetString(PyExc_ValueError, "ky must be between 1 and 5");
        return NULL;
    }
    /* Validate nux, nuy */
    if (nux < 0 || nux >= kx) {
        PyErr_SetString(PyExc_ValueError, "nux must be between 0 and kx-1");
        return NULL;
    }
    if (nuy < 0 || nuy >= ky) {
        PyErr_SetString(PyExc_ValueError, "nuy must be between 0 and ky-1");
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, NPY_FLOAT64, 1, 1);

    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL) {
        goto fail;
    }

    tx = (double *)PyArray_DATA(ap_tx);
    ty = (double *)PyArray_DATA(ap_ty);
    c = (double *)PyArray_DATA(ap_c);

    nx = (int)PyArray_DIMS(ap_tx)[0];
    ny = (int)PyArray_DIMS(ap_ty)[0];

    /* Check dimensions match pyf specification */
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError,
            "c must have length (nx-kx-1)*(ny-ky-1)");
        goto fail;
    }

    /* Allocate output array newc - same size as c */
    dims[0] = (nx - kx - 1) * (ny - ky - 1);
    ap_newc = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_newc == NULL) {
        goto fail;
    }
    newc = (double *)PyArray_DATA(ap_newc);

    /* Call the C function */
    pardtc(tx, nx, ty, ny, c, kx, ky, nux, nuy, newc, &ier);

    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_newc), ier);

fail:
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_newc);
    return NULL;
}

static char doc_pardeu[] = " [z,ier] = pardeu(tx,ty,c,kx,ky,nux,nuy,x,y)";
static PyObject *
fitpack_pardeu(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z = NULL;
    PyObject *tx, *ty, *c, *x, *y;
    int kx, ky, nux, nuy, m, nx, ny, lwrk, kwrk, ier;
    npy_intp dims[1];
    double *wrk = NULL;
    int *iwrk = NULL;

    if (!PyArg_ParseTuple(args, "OOOiiiiOO:pardeu",
                          &tx, &ty, &c, &kx, &ky, &nux, &nuy, &x, &y)) {
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y, NPY_FLOAT64, 1, 1);
    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL ||
        ap_x == NULL || ap_y == NULL) {
        goto fail;
    }

    nx = PyArray_DIMS(ap_tx)[0];
    ny = PyArray_DIMS(ap_ty)[0];
    m = PyArray_DIMS(ap_x)[0];

    /* Check c array dimensions */
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError, "Invalid c array dimensions");
        goto fail;
    }

    /* Check that x and y have same length */
    if (PyArray_DIMS(ap_y)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x and y arrays must have same length");
        goto fail;
    }

    /* Allocate output array */
    dims[0] = m;
    ap_z = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_z == NULL) {
        goto fail;
    }

    /* Allocate work arrays */
    lwrk = (nx * ny) + (kx + 1) * m + (ky + 1) * m;
    kwrk = m + m;
    wrk = malloc(lwrk * sizeof(double));
    iwrk = malloc(kwrk * sizeof(int));
    if (wrk == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Call the C function */
    pardeu((double *)PyArray_DATA(ap_tx), nx,
           (double *)PyArray_DATA(ap_ty), ny,
           (double *)PyArray_DATA(ap_c), kx, ky, nux, nuy,
           (double *)PyArray_DATA(ap_x),
           (double *)PyArray_DATA(ap_y),
           (double *)PyArray_DATA(ap_z), m,
           wrk, lwrk, iwrk, kwrk, &ier);

    free(wrk);
    free(iwrk);
    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);
    Py_DECREF(ap_x);
    Py_DECREF(ap_y);

    return Py_BuildValue(("Ni"), PyArray_Return(ap_z), ier);

fail:
    free(wrk);
    free(iwrk);
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    return NULL;
}

static char doc_dblint[] = " [iy] = dblint(tx,ty,c,kx,ky,xb,xe,yb,ye)";
static PyObject *
fitpack_dblint(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyObject *tx, *ty, *c;
    int kx, ky, nx, ny;
    double xb, xe, yb, ye, result;
    double *wrk = NULL;
    int wrk_size;

    if (!PyArg_ParseTuple(args, "OOOiidddd:dblint",
                          &tx, &ty, &c, &kx, &ky, &xb, &xe, &yb, &ye)) {
        return NULL;
    }

    ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx, NPY_FLOAT64, 1, 1);
    ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c, NPY_FLOAT64, 1, 1);
    if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL) {
        goto fail;
    }

    nx = PyArray_DIMS(ap_tx)[0];
    ny = PyArray_DIMS(ap_ty)[0];

    /* Check c array dimensions */
    if (PyArray_DIMS(ap_c)[0] != (nx - kx - 1) * (ny - ky - 1)) {
        PyErr_SetString(PyExc_ValueError, "Invalid c array dimensions");
        goto fail;
    }

    /* Allocate work array */
    wrk_size = nx + ny - kx - ky - 2;
    wrk = malloc(wrk_size * sizeof(double));
    if (wrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Call the C function */
    result = dblint((double *)PyArray_DATA(ap_tx), nx,
           (double *)PyArray_DATA(ap_ty), ny,
           (double *)PyArray_DATA(ap_c), kx, ky,
           xb, xe, yb, ye, wrk);

    free(wrk);
    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);
    Py_DECREF(ap_c);

    return Py_BuildValue("d", result);

fail:
    free(wrk);
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    return NULL;
}

static char doc_curfit[] = " [n,t,c,fp,ier] = curfit(iopt,x,y,w,xb,xe,k,s,t,wrk,iwrk)";
static PyObject *
fitpack_curfit(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_w = NULL;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_wrk = NULL, *ap_iwrk = NULL;
    PyObject *x, *y, *w, *t_in, *wrk_in, *iwrk_in;
    int iopt, m, k, nest, n, lwrk, ier;
    double xb, xe, s, fp;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "iOOOddidOOO:curfit",
                          &iopt, &x, &y, &w, &xb, &xe, &k, &s,
                          &t_in, &wrk_in, &iwrk_in)) {
        return NULL;
    }

    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y, NPY_FLOAT64, 1, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w, NPY_FLOAT64, 1, 1);
    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_in, NPY_FLOAT64, 1, 1);
    ap_wrk = (PyArrayObject *)PyArray_ContiguousFromObject(wrk_in, NPY_FLOAT64, 1, 1);
    ap_iwrk = (PyArrayObject *)PyArray_ContiguousFromObject(iwrk_in, NPY_INT32, 1, 1);
    if (ap_x == NULL || ap_y == NULL || ap_w == NULL ||
        ap_t == NULL || ap_wrk == NULL || ap_iwrk == NULL) {
        goto fail;
    }

    m = PyArray_DIMS(ap_x)[0];
    nest = PyArray_DIMS(ap_t)[0];
    lwrk = PyArray_DIMS(ap_wrk)[0];

    m = PyArray_DIMS(ap_x)[0];

    /* Validation checks */
    if (m <= k) {
        PyErr_SetString(PyExc_ValueError, "m must be > k");
        goto fail;
    }
    if (PyArray_DIMS(ap_y)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x and y arrays must have same length");
        goto fail;
    }
    if (PyArray_DIMS(ap_w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w array must have same length as x");
        goto fail;
    }
    if (PyArray_DIMS(ap_t)[0] != nest) {
        PyErr_SetString(PyExc_ValueError, "t array length must equal nest");
        goto fail;
    }

    lwrk = PyArray_DIMS(ap_wrk)[0];

    /* Allocate output c array */
    dims[0] = nest;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Initialize n */
    n = nest;

    /* Call the C function */
    curfit(iopt, m,
           (double *)PyArray_DATA(ap_x),
           (double *)PyArray_DATA(ap_y),
           (double *)PyArray_DATA(ap_w),
           xb, xe, k, s, nest, &n,
           (double *)PyArray_DATA(ap_t),
           (double *)PyArray_DATA(ap_c),
           &fp,
           (double *)PyArray_DATA(ap_wrk), lwrk,
           (int *)PyArray_DATA(ap_iwrk), &ier);

    /* Zero out the padding in c array (FITPACK convention: len(c) = len(t), but only c[0:n-k-1] are meaningful) */
    {
        double *c_data = (double *)PyArray_DATA(ap_c);
        int nc = n - k - 1;
        for (int i = nc; i < n; i++) {
            c_data[i] = 0.0;
        }
    }

    Py_DECREF(ap_x);
    Py_DECREF(ap_y);
    Py_DECREF(ap_w);
    /* Don't decref ap_t, ap_wrk, ap_iwrk - they are modified in-place */

    return Py_BuildValue(("iNNdi"),
                         n, PyArray_Return(ap_t), PyArray_Return(ap_c), fp, ier);

fail:
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;
}

static char doc_percur[] = " [n,t,c,fp,ier] = percur(iopt,x,y,w,k,s,nest,t,wrk,iwrk)";
static PyObject *
fitpack_percur(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_w = NULL;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_wrk = NULL, *ap_iwrk = NULL;
    PyObject *x, *y, *w, *t_in, *wrk_in, *iwrk_in;
    int iopt, m, k, nest, n, lwrk, ier;
    double s, fp;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "iOOOidOOO:percur",
                          &iopt, &x, &y, &w, &k, &s,
                          &t_in, &wrk_in, &iwrk_in)) {
        return NULL;
    }

    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y, NPY_FLOAT64, 1, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w, NPY_FLOAT64, 1, 1);
    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_in, NPY_FLOAT64, 1, 1);
    ap_wrk = (PyArrayObject *)PyArray_ContiguousFromObject(wrk_in, NPY_FLOAT64, 1, 1);
    ap_iwrk = (PyArrayObject *)PyArray_ContiguousFromObject(iwrk_in, NPY_INT32, 1, 1);
    if (ap_x == NULL || ap_y == NULL || ap_w == NULL ||
        ap_t == NULL || ap_wrk == NULL || ap_iwrk == NULL) {
        goto fail;
    }

    m = PyArray_DIMS(ap_x)[0];
    nest = PyArray_DIMS(ap_t)[0];
    lwrk = PyArray_DIMS(ap_wrk)[0];

    m = PyArray_DIMS(ap_x)[0];

    /* Validation checks */
    if (m <= k) {
        PyErr_SetString(PyExc_ValueError, "m must be > k");
        goto fail;
    }
    if (PyArray_DIMS(ap_y)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x and y arrays must have same length");
        goto fail;
    }
    if (PyArray_DIMS(ap_w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w array must have same length as x");
        goto fail;
    }
    if (PyArray_DIMS(ap_t)[0] != nest) {
        PyErr_SetString(PyExc_ValueError, "t array length must equal nest");
        goto fail;
    }

    lwrk = PyArray_DIMS(ap_wrk)[0];

    /* Allocate output c array */
    dims[0] = nest;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Initialize n */
    n = nest;

    /* Call the C function */
    percur(iopt, m,
           (double *)PyArray_DATA(ap_x),
           (double *)PyArray_DATA(ap_y),
           (double *)PyArray_DATA(ap_w),
           k, s, nest, &n,
           (double *)PyArray_DATA(ap_t),
           (double *)PyArray_DATA(ap_c),
           &fp,
           (double *)PyArray_DATA(ap_wrk), lwrk,
           (int *)PyArray_DATA(ap_iwrk), &ier);

    /* Zero out the padding in c array (FITPACK convention: len(c) = len(t), but only c[0:n-k-1] are meaningful) */
    {
        double *c_data = (double *)PyArray_DATA(ap_c);
        int nc = n - k - 1;
        for (int i = nc; i < n; i++) {
            c_data[i] = 0.0;
        }
    }

    Py_DECREF(ap_x);
    Py_DECREF(ap_y);
    Py_DECREF(ap_w);
    /* Don't decref ap_t, ap_wrk, ap_iwrk - they are modified in-place */

    return Py_BuildValue(("iNNdi"),
                         n, PyArray_Return(ap_t), PyArray_Return(ap_c), fp, ier);

fail:
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;
}

/**
  * There are three wrappers for surfit below, for different callsites in SciPy as legacy code.
  * Two are from the Fortran interface (_surfit_lsq and _surfit_smth in fitpack2.py) that is
  * used in LSQBivariateSpline and SmoothBivariateSpline.
  * The remaining is a a C wrapper for use in the legacy function "scipy.interpolate.bisplrep".
  * None of them call surfit as is but assume different parameters as optional or fixed values.
  * Worse is that bisplrep is stateful and keeps previous results in a global dict. Thus there
  * is some care needed to handle all that back and forth correctly.
*/
static char doc_surfit_lsq[] = " [c,fp,ier] = surfit_lsq(x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps,tx,ty)\n\nNotes\n-----\nThis LSQ-only wrapper modifies tx and ty in-place (both must be writable float64 ndarrays of length nmax).";
static PyObject *
fitpack_surfit_lsq(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    double *ap_x = NULL, *ap_y = NULL, *ap_z = NULL, *ap_w = NULL, *ap_tx = NULL, *ap_ty = NULL;
    PyArrayObject *ap_c = NULL;
    PyArrayObject *x = NULL, *y = NULL, *z = NULL, *w = NULL;
    PyArrayObject *tx = NULL, *ty = NULL;

    int m, kx, ky, nxest, nyest, nmax;
    int nx, ny, lwrk1, lwrk2, kwrk, ier;
    double xb, xe, yb, ye, s, eps, fp;

    double *wrk1 = NULL;
    double *wrk2 = NULL;
    int *iwrk = NULL;

    npy_intp dims[1];

    /* LSQ-only wrapper: iopt is fixed at -1; no external workspace/state. */
    if (!PyArg_ParseTuple(args, "O!O!O!O!ddddiidiiidO!O!:surfit_lsq",
                          &PyArray_Type, (PyObject **)&x,
                          &PyArray_Type, (PyObject **)&y,
                          &PyArray_Type, (PyObject **)&z,
                          &PyArray_Type, (PyObject **)&w,
                          &xb, &xe, &yb, &ye,
                          &kx, &ky,
                          &s,
                          &nxest, &nyest, &nmax,
                          &eps,
                          &PyArray_Type, (PyObject **)&tx,
                          &PyArray_Type, (PyObject **)&ty)) {
        return NULL;
    }

    ap_x = (double*)PyArray_DATA(x);
    ap_y = (double*)PyArray_DATA(y);
    ap_z = (double*)PyArray_DATA(z);
    ap_w = (double*)PyArray_DATA(w);
    ap_tx = (double*)PyArray_DATA(tx);
    ap_ty = (double*)PyArray_DATA(ty);

    m = (int)PyArray_DIMS(x)[0];

    if (PyArray_DIMS(y)[0] != m || PyArray_DIMS(z)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x, y, z arrays must have the same length");
        goto fail;
    }
    if (m < (kx + 1) * (ky + 1)) {
        PyErr_SetString(PyExc_ValueError, "m must be >= (kx+1)*(ky+1)");
        goto fail;
    }
    if (nxest <= 0 || nyest <= 0) {
        PyErr_SetString(PyExc_ValueError, "nxest and nyest must be positive");
        goto fail;
    }
    if (nmax <= 0) {
        PyErr_SetString(PyExc_ValueError, "nmax must be positive");
        goto fail;
    }
    if (nmax < nxest || nmax < nyest) {
        /* This is the key safety condition: surfit writes into tx/ty buffers of length nmax. */
        PyErr_SetString(PyExc_ValueError, "nmax must be >= max(nxest, nyest)");
        goto fail;
    }
    if (PyArray_NDIM(tx) != 1 || PyArray_NDIM(ty) != 1) {
        PyErr_SetString(PyExc_ValueError, "tx and ty must be 1-D arrays");
        goto fail;
    }
    if ((int)PyArray_DIMS(tx)[0] != nmax) {
        PyErr_SetString(PyExc_ValueError, "tx must have length == nmax");
        goto fail;
    }
    if ((int)PyArray_DIMS(ty)[0] != nmax) {
        PyErr_SetString(PyExc_ValueError, "ty must have length == nmax");
        goto fail;
    }
    if ((int)PyArray_DIMS(w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w must have the same length as x, y, z");
        goto fail;
    }

    nx = nxest;
    ny = nyest;

    /* Workspace sizing: same formulas as legacy f2py wrapper, but allocated via malloc. */
    {
        int u = nxest - kx - 1;
        int v = nyest - ky - 1;
        int km = (kx > ky ? kx : ky) + 1;
        int ne = (nxest > nyest ? nxest : nyest);
        int bx = kx * v + ky + 1;
        int by = ky * u + kx + 1;
        int b1 = bx;
        int b2 = bx + v - ky;
        if (bx > by) {
            b1 = by;
            b2 = by + u - kx;
        }

        lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
        lwrk2 = u * v * (b2 + 1) + b2;
        kwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);

        if (lwrk1 <= 0 || lwrk2 <= 0 || kwrk <= 0) {
            PyErr_SetString(PyExc_ValueError, "invalid workspace size computed");
            goto fail;
        }
    }

    wrk1 = (double *)malloc((size_t)lwrk1 * sizeof(double));
    wrk2 = (double *)malloc((size_t)lwrk2 * sizeof(double));
    iwrk = (int *)malloc((size_t)kwrk * sizeof(int));
    if (wrk1 == NULL || wrk2 == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }


    npy_intp nc = (npy_intp)(nxest - kx - 1) * (npy_intp)(nyest - ky - 1);
    if (nc < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid coefficient size");
        goto fail;
    }
    dims[0] = nc;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Call underlying routine: iopt=-1 (fixed knots LSQ). */
    surfit(-1, m, ap_x, ap_y, ap_z, ap_w, xb, xe, yb, ye, kx, ky, s,
           nxest, nyest, nmax, eps,
           &nx, ap_tx, &ny, ap_ty,
           (double *)PyArray_DATA(ap_c), &fp,
           wrk1, lwrk1,
           wrk2, lwrk2,
           iwrk, kwrk, &ier);

    free(wrk1);
    free(wrk2);
    free(iwrk);

    return Py_BuildValue("Ndi", PyArray_Return(ap_c), fp, ier);

fail:
    if (wrk1 != NULL) {
        free(wrk1);
    }
    if (wrk2 != NULL) {
        free(wrk2);
    }
    if (iwrk != NULL) {
        free(iwrk);
    }
    return NULL;
}


static char doc_surfit_smth[] = " [nx,tx,ny,ty,c,fp,ier] = surfit_smth(x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps);";
static PyObject *
fitpack_surfit_smth(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    double *ap_x = NULL, *ap_y = NULL, *ap_z = NULL, *ap_w = NULL;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *x = NULL, *y = NULL, *z = NULL, *w = NULL;

    int m, kx, ky, nxest, nyest, nmax;
    int nx, ny, lwrk1, lwrk2, kwrk, ier;
    double xb, xe, yb, ye, s, eps, fp;

    double *wrk1 = NULL;
    double *wrk2 = NULL;
    int *iwrk = NULL;

    npy_intp dims[1];

    /* Smoothing-only wrapper: iopt is fixed at 0; allocates tx, ty, c internally. */
    if (!PyArg_ParseTuple(args, "O!O!O!O!ddddiidiiid:surfit_smth",
                          &PyArray_Type, (PyObject **)&x,
                          &PyArray_Type, (PyObject **)&y,
                          &PyArray_Type, (PyObject **)&z,
                          &PyArray_Type, (PyObject **)&w,
                          &xb, &xe, &yb, &ye,
                          &kx, &ky,
                          &s,
                          &nxest, &nyest, &nmax,
                          &eps)) {
        return NULL;
    }

    ap_x = (double*)PyArray_DATA(x);
    ap_y = (double*)PyArray_DATA(y);
    ap_z = (double*)PyArray_DATA(z);
    ap_w = (double*)PyArray_DATA(w);

    m = (int)PyArray_DIMS(x)[0];

    if (PyArray_DIMS(y)[0] != m || PyArray_DIMS(z)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x, y, z arrays must have the same length");
        goto fail;
    }
    if (PyArray_DIMS(w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w must have the same length as x, y, z");
        goto fail;
    }
    if (m < (kx + 1) * (ky + 1)) {
        PyErr_SetString(PyExc_ValueError, "m must be >= (kx+1)*(ky+1)");
        goto fail;
    }
    if (nxest <= 0 || nyest <= 0) {
        PyErr_SetString(PyExc_ValueError, "nxest and nyest must be positive");
        goto fail;
    }
    if (nmax <= 0) {
        PyErr_SetString(PyExc_ValueError, "nmax must be positive");
        goto fail;
    }
    if (nmax < nxest || nmax < nyest) {
        PyErr_SetString(PyExc_ValueError, "nmax must be >= max(nxest, nyest)");
        goto fail;
    }

    /* Workspace sizing: same formulas as LSQ mode. */
    {
        int u = nxest - kx - 1;
        int v = nyest - ky - 1;
        int km = (kx > ky ? kx : ky) + 1;
        int ne = (nxest > nyest ? nxest : nyest);
        int bx = kx * v + ky + 1;
        int by = ky * u + kx + 1;
        int b1 = bx;
        int b2 = bx + v - ky;
        if (bx > by) {
            b1 = by;
            b2 = by + u - kx;
        }

        lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
        lwrk2 = u * v * (b2 + 1) + b2;
        kwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);

        if (lwrk1 <= 0 || lwrk2 <= 0 || kwrk <= 0) {
            PyErr_SetString(PyExc_ValueError, "invalid workspace size computed");
            goto fail;
        }
    }

    wrk1 = (double *)malloc((size_t)lwrk1 * sizeof(double));
    wrk2 = (double *)malloc((size_t)lwrk2 * sizeof(double));
    iwrk = (int *)malloc((size_t)kwrk * sizeof(int));
    if (wrk1 == NULL || wrk2 == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Allocate output arrays: tx, ty of length nmax, c of length (nxest-kx-1)*(nyest-ky-1). */
    dims[0] = nmax;
    ap_tx = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tx == NULL) {
        goto fail;
    }

    dims[0] = nmax;
    ap_ty = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_ty == NULL) {
        goto fail;
    }

    npy_intp nc = (npy_intp)(nxest - kx - 1) * (npy_intp)(nyest - ky - 1);
    if (nc < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid coefficient size");
        goto fail;
    }
    dims[0] = nc;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Call underlying routine: iopt=0 (smoothing spline). */
    surfit(0, m, ap_x, ap_y, ap_z, ap_w, xb, xe, yb, ye, kx, ky, s,
           nxest, nyest, nmax, eps,
           &nx, (double *)PyArray_DATA(ap_tx),
           &ny, (double *)PyArray_DATA(ap_ty),
           (double *)PyArray_DATA(ap_c), &fp,
           wrk1, lwrk1,
           wrk2, lwrk2,
           iwrk, kwrk, &ier);

    free(wrk1);
    free(wrk2);
    free(iwrk);

    return Py_BuildValue("iNiNNdi", nx, PyArray_Return(ap_tx), ny, PyArray_Return(ap_ty), PyArray_Return(ap_c), fp, ier);

fail:
    if (wrk1 != NULL) {
        free(wrk1);
    }
    if (wrk2 != NULL) {
        free(wrk2);
    }
    if (iwrk != NULL) {
        free(iwrk);
    }
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    return NULL;
}


static char doc_surfit[] = " [nx,tx,ny,ty,c,fp,ier,wrk] = surfit(iopt,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps,wrk)\n\nNotes\n-----\nThis is the legacy bisplrep-compatible wrapper that handles iopt=0,1,-1 and maintains wrk state.";
static PyObject *
fitpack_surfit(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    double *ap_x = NULL, *ap_y = NULL, *ap_z = NULL, *ap_w = NULL;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_wrk_in = NULL, *ap_wrk_out = NULL;
    PyArrayObject *x = NULL, *y = NULL, *z = NULL, *w = NULL;
    PyObject *wrk_in = NULL;

    int iopt, m, kx, ky, nxest, nyest, nmax;
    int nx, ny, lwrk1, lwrk2, kwrk, ier;
    double xb, xe, yb, ye, s, eps, fp;

    double *wrk1 = NULL;
    double *wrk2 = NULL;
    int *iwrk = NULL;

    npy_intp dims[1];

    /* Legacy stateful wrapper for bisplrep: handles iopt=0,1,-1 and returns wrk for caching. */
    if (!PyArg_ParseTuple(args, "iO!O!O!O!ddddiidiiidO!:surfit",
                          &iopt,
                          &PyArray_Type, (PyObject **)&x,
                          &PyArray_Type, (PyObject **)&y,
                          &PyArray_Type, (PyObject **)&z,
                          &PyArray_Type, (PyObject **)&w,
                          &xb, &xe, &yb, &ye,
                          &kx, &ky,
                          &s,
                          &nxest, &nyest, &nmax,
                          &eps,
                          &PyArray_Type, (PyObject **)&wrk_in)) {
        return NULL;
    }

    ap_x = (double*)PyArray_DATA(x);
    ap_y = (double*)PyArray_DATA(y);
    ap_z = (double*)PyArray_DATA(z);
    ap_w = (double*)PyArray_DATA(w);

    m = (int)PyArray_DIMS(x)[0];

    if (PyArray_DIMS(y)[0] != m || PyArray_DIMS(z)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "x, y, z arrays must have the same length");
        goto fail;
    }
    if (PyArray_DIMS(w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w must have the same length as x, y, z");
        goto fail;
    }
    if (m < (kx + 1) * (ky + 1)) {
        PyErr_SetString(PyExc_ValueError, "m must be >= (kx+1)*(ky+1)");
        goto fail;
    }
    if (nxest <= 0 || nyest <= 0) {
        PyErr_SetString(PyExc_ValueError, "nxest and nyest must be positive");
        goto fail;
    }
    if (nmax <= 0) {
        PyErr_SetString(PyExc_ValueError, "nmax must be positive");
        goto fail;
    }
    if (nmax < nxest || nmax < nyest) {
        PyErr_SetString(PyExc_ValueError, "nmax must be >= max(nxest, nyest)");
        goto fail;
    }

    /* Workspace sizing: same formulas as other wrappers. */
    {
        int u = nxest - kx - 1;
        int v = nyest - ky - 1;
        int km = (kx > ky ? kx : ky) + 1;
        int ne = (nxest > nyest ? nxest : nyest);
        int bx = kx * v + ky + 1;
        int by = ky * u + kx + 1;
        int b1 = bx;
        int b2 = bx + v - ky;
        if (bx > by) {
            b1 = by;
            b2 = by + u - kx;
        }

        lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
        lwrk2 = u * v * (b2 + 1) + b2;
        kwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);

        if (lwrk1 <= 0 || lwrk2 <= 0 || kwrk <= 0) {
            PyErr_SetString(PyExc_ValueError, "invalid workspace size computed");
            goto fail;
        }
    }

    wrk1 = (double *)malloc((size_t)lwrk1 * sizeof(double));
    wrk2 = (double *)malloc((size_t)lwrk2 * sizeof(double));
    iwrk = (int *)malloc((size_t)kwrk * sizeof(int));
    if (wrk1 == NULL || wrk2 == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* For iopt=1, copy wrk_in into wrk1 for continuation. */
    if (iopt == 1) {
        ap_wrk_in = (PyArrayObject *)wrk_in;
        npy_intp wrk_in_size = PyArray_DIMS(ap_wrk_in)[0];
        npy_intp lc = (npy_intp)(nxest - kx - 1) * (npy_intp)(nyest - ky - 1);

        if (wrk_in_size < lc) {
            PyErr_SetString(PyExc_ValueError, "wrk array too small for iopt=1");
            goto fail;
        }

        memcpy(wrk1, PyArray_DATA(ap_wrk_in), (size_t)lc * sizeof(double));
    }

    /* Allocate output arrays: tx, ty of length nmax, c of length (nxest-kx-1)*(nyest-ky-1). */
    dims[0] = nmax;
    ap_tx = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tx == NULL) {
        goto fail;
    }

    dims[0] = nmax;
    ap_ty = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_ty == NULL) {
        goto fail;
    }

    npy_intp nc = (npy_intp)(nxest - kx - 1) * (npy_intp)(nyest - ky - 1);
    if (nc < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid coefficient size");
        goto fail;
    }
    dims[0] = nc;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Call underlying routine with appropriate iopt. */
    surfit(iopt, m, ap_x, ap_y, ap_z, ap_w, xb, xe, yb, ye, kx, ky, s,
           nxest, nyest, nmax, eps,
           &nx, (double *)PyArray_DATA(ap_tx),
           &ny, (double *)PyArray_DATA(ap_ty),
           (double *)PyArray_DATA(ap_c), &fp,
           wrk1, lwrk1,
           wrk2, lwrk2,
           iwrk, kwrk, &ier);

    /* Create wrk_out array containing wrk1 for state preservation. */
    dims[0] = nc;  /* wrk size should match coefficient array size */
    ap_wrk_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_wrk_out == NULL) {
        goto fail;
    }
    memcpy(PyArray_DATA(ap_wrk_out), wrk1, (size_t)nc * sizeof(double));

    free(wrk1);
    free(wrk2);
    free(iwrk);

    /* Slice tx, ty to actual sizes nx, ny for return. */
    PyArrayObject *ap_tx_sliced = NULL;
    PyArrayObject *ap_ty_sliced = NULL;

    dims[0] = nx;
    ap_tx_sliced = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tx_sliced == NULL) {
        Py_DECREF(ap_wrk_out);
        goto fail;
    }
    memcpy(PyArray_DATA(ap_tx_sliced), PyArray_DATA(ap_tx), (size_t)nx * sizeof(double));

    dims[0] = ny;
    ap_ty_sliced = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_ty_sliced == NULL) {
        Py_DECREF(ap_tx_sliced);
        Py_DECREF(ap_wrk_out);
        goto fail;
    }
    memcpy(PyArray_DATA(ap_ty_sliced), PyArray_DATA(ap_ty), (size_t)ny * sizeof(double));

    Py_DECREF(ap_tx);
    Py_DECREF(ap_ty);

    return Py_BuildValue("iNiNNdiN", nx, PyArray_Return(ap_tx_sliced), ny, PyArray_Return(ap_ty_sliced),
                         PyArray_Return(ap_c), fp, ier, PyArray_Return(ap_wrk_out));

fail:
    if (wrk1 != NULL) {
        free(wrk1);
    }
    if (wrk2 != NULL) {
        free(wrk2);
    }
    if (iwrk != NULL) {
        free(iwrk);
    }
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    return NULL;
}


static char doc_regrid[] = " [nx,tx,ny,ty,c,fp,ier] = regrid(iopt,x,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,maxit,wrk,iwrk)";
static PyObject *
fitpack_regrid(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_x = NULL, *ap_y = NULL, *ap_z = NULL;
    PyArrayObject *ap_tx = NULL, *ap_ty = NULL, *ap_c = NULL;
    PyArrayObject *ap_wrk = NULL, *ap_iwrk = NULL;
    PyObject *x, *y, *z, *wrk_in, *iwrk_in;
    int iopt, mx, my, kx, ky, nxest, nyest, maxit, nx, ny, lwrk, kwrk, ier;
    double xb, xe, yb, ye, s, fp;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "iOOOddddiidiiiOO:regrid",
                          &iopt, &x, &y, &z, &xb, &xe, &yb, &ye,
                          &kx, &ky, &s, &nxest, &nyest, &maxit, &wrk_in, &iwrk_in)) {
        return NULL;
    }

    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x, NPY_FLOAT64, 1, 1);
    ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y, NPY_FLOAT64, 1, 1);
    ap_z = (PyArrayObject *)PyArray_ContiguousFromObject(z, NPY_FLOAT64, 1, 1);
    ap_wrk = (PyArrayObject *)PyArray_ContiguousFromObject(wrk_in, NPY_FLOAT64, 1, 1);
    ap_iwrk = (PyArrayObject *)PyArray_ContiguousFromObject(iwrk_in, NPY_INT32, 1, 1);
    if (ap_x == NULL || ap_y == NULL || ap_z == NULL ||
        ap_wrk == NULL || ap_iwrk == NULL) {
        goto fail;
    }

    mx = PyArray_DIMS(ap_x)[0];
    my = PyArray_DIMS(ap_y)[0];

    /* Validation checks */
    if (mx <= kx) {
        PyErr_SetString(PyExc_ValueError, "mx must be > kx");
        goto fail;
    }
    if (my <= ky) {
        PyErr_SetString(PyExc_ValueError, "my must be > ky");
        goto fail;
    }
    if (PyArray_DIMS(ap_z)[0] != mx * my) {
        PyErr_SetString(PyExc_ValueError, "z array length must equal mx*my");
        goto fail;
    }

    lwrk = PyArray_DIMS(ap_wrk)[0];
    kwrk = PyArray_DIMS(ap_iwrk)[0];

    /* Allocate output arrays */
    dims[0] = nxest;
    ap_tx = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tx == NULL) {
        goto fail;
    }

    dims[0] = nyest;
    ap_ty = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_ty == NULL) {
        goto fail;
    }

    dims[0] = (nxest - kx - 1) * (nyest - ky - 1);
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Call the C function */
    regrid(iopt, mx,
           (double *)PyArray_DATA(ap_x), my,
           (double *)PyArray_DATA(ap_y),
           (double *)PyArray_DATA(ap_z),
           xb, xe, yb, ye, kx, ky, s,
           nxest, nyest, maxit,
           &nx, (double *)PyArray_DATA(ap_tx),
           &ny, (double *)PyArray_DATA(ap_ty),
           (double *)PyArray_DATA(ap_c), &fp,
           (double *)PyArray_DATA(ap_wrk), lwrk,
           (int *)PyArray_DATA(ap_iwrk), kwrk, &ier);

    Py_DECREF(ap_x);
    Py_DECREF(ap_y);
    Py_DECREF(ap_z);
    Py_DECREF(ap_wrk);
    Py_DECREF(ap_iwrk);

    return Py_BuildValue(("iNiNNdi"),
                         nx, PyArray_Return(ap_tx), ny, PyArray_Return(ap_ty),
                         PyArray_Return(ap_c), fp, ier);

fail:
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_z);
    Py_XDECREF(ap_tx);
    Py_XDECREF(ap_ty);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;
}

static char doc_sphere[] = " [nt,tt,np,tp,c,fp,ier] = sphere(iopt,teta,phi,r,w,s,ntest,npest,tt,tp,eps,wrk1,wrk2,iwrk)";
static PyObject *
fitpack_sphere(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_teta = NULL, *ap_phi = NULL, *ap_r = NULL, *ap_w = NULL;
    PyArrayObject *ap_tt = NULL, *ap_tp = NULL, *ap_c = NULL;
    double *ap_wrk1 = NULL, *ap_wrk2 = NULL;
    int *ap_iwrk = NULL;
    PyArrayObject *tt = NULL, *tp =  NULL, *wrk1 = NULL, *wrk2 = NULL, *iwrk = NULL;
    PyObject *teta, *phi, *r, *w;
    int iopt, m, ntest, npest, nt, np, lwrk1, lwrk2, kwrk, ier;
    double s, eps, fp;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "iOOOOdiiO!O!dO!O!O!:sphere",
                          &iopt,                             // i
                          &teta,                             // O
                          &phi,                              // O
                          &r,                                // O
                          &w,                                // O
                          &s,                                // d
                          &ntest,                            // i
                          &npest,                            // i
                          &PyArray_Type, (PyObject **)&tt,   // O!
                          &PyArray_Type, (PyObject **)&tp,   // O!
                          &eps,                              // d
                          &PyArray_Type, (PyObject **)&wrk1, // O!
                          &PyArray_Type, (PyObject **)&wrk2, // O!
                          &PyArray_Type, (PyObject **)&iwrk  // O!
    )) {
        return NULL;
    }

    ap_teta = (PyArrayObject *)PyArray_ContiguousFromObject(teta, NPY_FLOAT64, 1, 1);
    ap_phi = (PyArrayObject *)PyArray_ContiguousFromObject(phi, NPY_FLOAT64, 1, 1);
    ap_r = (PyArrayObject *)PyArray_ContiguousFromObject(r, NPY_FLOAT64, 1, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w, NPY_FLOAT64, 1, 1);
    ap_wrk1 = (double*)PyArray_DATA(wrk1);
    ap_wrk2 = (double*)PyArray_DATA(wrk2);
    ap_iwrk = (int*)PyArray_DATA(iwrk);
    if (ap_teta == NULL || ap_phi == NULL || ap_r == NULL || ap_w == NULL) {
        goto fail;
    }

    // Currently SciPy FITPACK invocations are inconsistent due to legacy reasons
    // _spherefit_lsq is meant for iopt=-1 workflow and hence passes tt and tp of
    // size ntest and npest respectively. However, _spherefit_smooth does not
    // and generates tt and tp as output arguments.
    // Hence we have to accomodate for these cases here to avoid multiple wrappers
    // as was the case previously with the Fortran FITPACK bindings.

    if (iopt == -1){
        // LSQ
        // tt and tp are input arrays that will be modified in-place
        // Acquire pointers to the input arrays without copy
        ap_tt = tt;
        ap_tp = tp;

    } else if (iopt == 0) {
        // SMOOTH
        // tt and tp are output arrays to be allocated
        dims[0] = ntest;
        ap_tt = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_tt == NULL) {
            goto fail;
        }

        dims[0] = npest;
        ap_tp = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_tp == NULL) {
            goto fail;
        }
    }
    m = PyArray_DIMS(ap_teta)[0];
    // Initialize nt and np, In the smooth case they are ignored and populated
    nt = ntest;
    np = npest;
    /* Validation checks */
    if (m < 2) {
        PyErr_SetString(PyExc_ValueError, "m must be >= 2");
        goto fail;
    }
    if (PyArray_DIMS(ap_phi)[0] != m || PyArray_DIMS(ap_r)[0] != m ||
        PyArray_DIMS(ap_w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "teta, phi, r, w arrays must have same length");
        goto fail;
    }

    lwrk1 = PyArray_DIMS(wrk1)[0];
    lwrk2 = PyArray_DIMS(wrk2)[0];
    kwrk = PyArray_DIMS(iwrk)[0];

    dims[0] = (ntest - 4) * (npest - 4);
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        goto fail;
    }

    /* Call the C function */
    sphere(iopt, m,
           (double *)PyArray_DATA(ap_teta),
           (double *)PyArray_DATA(ap_phi),
           (double *)PyArray_DATA(ap_r),
           (double *)PyArray_DATA(ap_w),
           s, ntest, npest, eps,
           &nt, (double *)PyArray_DATA(ap_tt),
           &np, (double *)PyArray_DATA(ap_tp),
           (double *)PyArray_DATA(ap_c), &fp,
           ap_wrk1, lwrk1,
           ap_wrk2, lwrk2,
           ap_iwrk, kwrk, &ier);

    Py_DECREF(ap_teta);
    Py_DECREF(ap_phi);
    Py_DECREF(ap_r);
    Py_DECREF(ap_w);

    if (iopt == -1) {
        // LSQ: tt and tp were borrowed, INCREF before returning.
        Py_INCREF(ap_tt);
        Py_INCREF(ap_tp);
    }
    return Py_BuildValue(("iNiNNdi"),
                         nt, PyArray_Return(ap_tt), np, PyArray_Return(ap_tp),
                         PyArray_Return(ap_c), fp, ier);

fail:
    Py_XDECREF(ap_teta);
    Py_XDECREF(ap_phi);
    Py_XDECREF(ap_r);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_tt);
    Py_XDECREF(ap_tp);
    Py_XDECREF(ap_c);
    return NULL;
}

static char doc_spgrid[] = " [nu,tu,nv,tv,c,fp,ier] = spgrid(iopt,ider,u,v,r,r0,r1,s,nuest,nvest)";
static PyObject *
fitpack_spgrid(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    int *ap_iopt = NULL, *ap_ider = NULL;
    double *ap_u = NULL, *ap_v = NULL, *ap_r = NULL;
    PyArrayObject *iopt = NULL, *ider = NULL, *u = NULL, *v = NULL, *r = NULL;
    PyArrayObject *ap_tu = NULL, *ap_tv = NULL, *ap_c = NULL;
    PyArrayObject *ap_tu_full = NULL, *ap_tv_full = NULL;
    PyArrayObject *ap_c_full = NULL;

    int mu, mv, nuest, nvest, nu, nv, lwrk, kwrk, ier;
    double r0, r1, s, fp;
    double *wrk = NULL;
    int *iwrk = NULL;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!dddii:spgrid",
                          &PyArray_Type, (PyObject **)&iopt,
                          &PyArray_Type, (PyObject **)&ider,
                          &PyArray_Type, (PyObject **)&u,
                          &PyArray_Type, (PyObject **)&v,
                          &PyArray_Type, (PyObject **)&r,
                          &r0, &r1, &s, &nuest, &nvest)) {
        return NULL;
    }

    ap_iopt = (int *)PyArray_DATA(iopt);
    ap_ider = (int *)PyArray_DATA(ider);
    ap_u = (double *)PyArray_DATA(u);
    ap_v = (double *)PyArray_DATA(v);
    ap_r = (double *)PyArray_DATA(r);

    mu = (int)PyArray_DIMS(u)[0];
    mv = (int)PyArray_DIMS(v)[0];

    /* Validation checks */
    if (PyArray_DIMS(iopt)[0] != 3) {
        PyErr_SetString(PyExc_ValueError, "iopt array must have length 3");
        goto fail;
    }
    if (PyArray_DIMS(ider)[0] != 4) {
        PyErr_SetString(PyExc_ValueError, "ider array must have length 4");
        goto fail;
    }
    if (PyArray_DIMS(r)[0] != mu * mv) {
        PyErr_SetString(PyExc_ValueError, "r array length must equal mu*mv");
        goto fail;
    }

    /* Calculate workspace requirements */
    lwrk = 12 + nuest * (mv + nvest + 3) + nvest * 24 + 4 * mu + 8 * mv +
           (mv + nvest > nuest ? mv + nvest : nuest);
    kwrk = 5 + mu + mv + nuest + nvest;

    /* Allocate workspace arrays */
    wrk = malloc(lwrk * sizeof(double));
    iwrk = malloc(kwrk * sizeof(int));
    if (wrk == NULL || iwrk == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Allocate full-size output arrays (will slice later) */
    dims[0] = nuest;
    ap_tu_full = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tu_full == NULL) {
        goto fail;
    }

    dims[0] = nvest;
    ap_tv_full = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tv_full == NULL) {
        goto fail;
    }

    dims[0] = (nuest - 4) * (nvest - 4);
    ap_c_full = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c_full == NULL) {
        goto fail;
    }

    /* Call the C function */
    spgrid(ap_iopt, ap_ider,
           mu, ap_u,
           mv, ap_v,
           ap_r,
           r0, r1, s, nuest, nvest,
           &nu, (double *)PyArray_DATA(ap_tu_full),
           &nv, (double *)PyArray_DATA(ap_tv_full),
            (double *)PyArray_DATA(ap_c_full), &fp,
           wrk, lwrk, iwrk, kwrk, &ier);

    free(wrk);
    free(iwrk);

    /* Check for errors before creating sliced arrays */
    if (ier == 10) {
        /* Input validation error - nu and nv are not initialized */
        Py_DECREF(ap_tu_full);
        Py_DECREF(ap_tv_full);
        Py_DECREF(ap_c_full);
        PyErr_SetString(PyExc_ValueError, "Invalid input parameters to spgrid");
        return NULL;
    }

    /* Slice output arrays to actual sizes */
    dims[0] = nu;
    ap_tu = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tu == NULL) {
        Py_DECREF(ap_tu_full);
        Py_DECREF(ap_tv_full);
        Py_DECREF(ap_c);
        return NULL;
    }
    memcpy(PyArray_DATA(ap_tu), PyArray_DATA(ap_tu_full), (size_t)nu * sizeof(double));

    dims[0] = nv;
    ap_tv = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_tv == NULL) {
        Py_DECREF(ap_tu);
        Py_DECREF(ap_tu_full);
        Py_DECREF(ap_tv_full);
        Py_DECREF(ap_c_full);
        return NULL;
    }
    memcpy(PyArray_DATA(ap_tv), PyArray_DATA(ap_tv_full), (size_t)nv * sizeof(double));

    Py_DECREF(ap_tu_full);
    Py_DECREF(ap_tv_full);

    /* Slice coefficients to actual size (nu-4)*(nv-4). */
    {
        npy_intp ncof = (npy_intp)(nu - 4) * (npy_intp)(nv - 4);
        dims[0] = ncof;
        ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_c == NULL) {
            Py_DECREF(ap_tu);
            Py_DECREF(ap_tv);
            Py_DECREF(ap_c_full);
            return NULL;
        }
        memcpy(PyArray_DATA(ap_c), PyArray_DATA(ap_c_full), (size_t)ncof * sizeof(double));
    }

    Py_DECREF(ap_c_full);

    return Py_BuildValue("iNiNNdi",
                         nu, PyArray_Return(ap_tu),
                         nv, PyArray_Return(ap_tv),
                         PyArray_Return(ap_c), fp, ier);

fail:
    free(wrk);
    free(iwrk);
    Py_XDECREF(ap_tu_full);
    Py_XDECREF(ap_tv_full);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_c_full);
    return NULL;
}

static char doc_parcur[] = " [n,t,c,fp,ier] = parcur(iopt,ipar,idim,u,x,w,ub,ue,k,s,nest,t,wrk,iwrk)";
static PyObject *
fitpack_parcur(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_u = NULL, *ap_x = NULL, *ap_w = NULL;
    PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_wrk = NULL, *ap_iwrk = NULL;
    PyObject *u, *x, *w, *t_in, *wrk_in, *iwrk_in;
    int iopt, ipar, idim, m, mx, k, nest, n, nc, lwrk, ier, per;
    double ub, ue, s, fp;
    npy_intp dims[1];

    /* Signature: x, w, u, ub, ue, k, task, ipar, s, t, nest, wrk, iwrk, per */
    if (!PyArg_ParseTuple(args, "OOOddiiidOiOOi:parcur",
                          &x, &w, &u, &ub, &ue, &k, &iopt, &ipar, &s,
                          &t_in, &nest, &wrk_in, &iwrk_in, &per)) {
        return NULL;
    }

    ap_u = (PyArrayObject *)PyArray_ContiguousFromObject(u, NPY_FLOAT64, 1, 1);
    ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x, NPY_FLOAT64, 1, 1);
    ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w, NPY_FLOAT64, 1, 1);
    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_in, NPY_FLOAT64, 1, 1);
    ap_wrk = (PyArrayObject *)PyArray_ContiguousFromObject(wrk_in, NPY_FLOAT64, 1, 1);
    ap_iwrk = (PyArrayObject *)PyArray_ContiguousFromObject(iwrk_in, NPY_INT32, 1, 1);
    if (ap_u == NULL || ap_x == NULL || ap_w == NULL ||
        ap_t == NULL || ap_wrk == NULL || ap_iwrk == NULL) {
        goto fail;
    }

    m = PyArray_DIMS(ap_u)[0];
    mx = PyArray_DIMS(ap_x)[0];

    /* Derive idim from x array size: idim = mx / m */
    if (mx % m != 0) {
        PyErr_SetString(PyExc_ValueError, "len(x) must be divisible by len(u)");
        goto fail;
    }
    idim = mx / m;

    /* Validation checks */
    if (m <= k) {
        PyErr_SetString(PyExc_ValueError, "m must be > k");
        goto fail;
    }
    if (iopt < -1 || iopt > 1) {
        PyErr_SetString(PyExc_ValueError, "iopt must be -1, 0, or 1");
        goto fail;
    }
    if (ipar != 0 && ipar != 1) {
        PyErr_SetString(PyExc_ValueError, "ipar must be 0 or 1");
        goto fail;
    }
    if (idim < 1 || idim > 10) {
        PyErr_SetString(PyExc_ValueError, "idim must be between 1 and 10");
        goto fail;
    }
    if (mx < idim * m) {
        PyErr_SetString(PyExc_ValueError, "x array length must be >= idim*m");
        goto fail;
    }
    if (PyArray_DIMS(ap_w)[0] != m) {
        PyErr_SetString(PyExc_ValueError, "w array must have same length as u");
        goto fail;
    }

    /* For task >= 0, t can be empty initially, so don't check its size */
    /* The t array will be resized appropriately by parcur */

    /* Calculate required work array sizes */
    int k1 = k + 1;
    /* Calculate lwest based on which function will be called:
     * parcur: m*k1 + nest*(6 + idim + 3*k)
     * clocur: m*k1 + nest*(7 + idim + 5*k)
     */
    int lwest;
    if (per == 0) {
        lwest = m * k1 + nest * (6 + idim + 3*k);
    } else {
        lwest = m * k1 + nest * (7 + idim + 5*k);
    }

    lwrk = PyArray_DIMS(ap_wrk)[0];
    if (lwrk < lwest) {
        Py_DECREF(ap_wrk);
        dims[0] = lwest;
        ap_wrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_wrk == NULL) {
            goto fail;
        }
        lwrk = lwest;
    }

    /* Ensure iwrk is at least nest elements */
    int kwrk = PyArray_DIMS(ap_iwrk)[0];
    if (kwrk < nest) {
        Py_DECREF(ap_iwrk);
        dims[0] = nest;
        ap_iwrk = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT32);
        if (ap_iwrk == NULL) {
            goto fail;
        }
    }

    nc = idim * nest;

    /* Pad t array to size nest */
    PyArrayObject *ap_t_pad = NULL;
    dims[0] = nest;
    ap_t_pad = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
    if (ap_t_pad == NULL) {
        goto fail;
    }
    int t_len = PyArray_DIMS(ap_t)[0];
    if (t_len > 0) {
        memcpy(PyArray_DATA(ap_t_pad), PyArray_DATA(ap_t),
               (t_len < nest ? t_len : nest) * sizeof(double));
    }

    /* Allocate output c array */
    dims[0] = nc;
    ap_c = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (ap_c == NULL) {
        Py_DECREF(ap_t_pad);
        goto fail;
    }

    /* Initialize n based on the input t array length
     * For iopt >= 0 (compute new knots), if t is empty, n starts at 0
     * For iopt == 1 (continue from previous), n is the length of t from cache
     * For iopt == -1 (user-provided knots), n is the length of the input t
     */
    n = t_len;

    /* Call the appropriate C function based on per parameter */
    if (per == 0) {
        /* parcur for non-periodic curves - has ub, ue parameters */
        parcur(iopt, ipar, idim, m,
               (double *)PyArray_DATA(ap_u), mx,
               (double *)PyArray_DATA(ap_x),
               (double *)PyArray_DATA(ap_w),
               &ub, &ue,
               k, s, nest, &n,
               (double *)PyArray_DATA(ap_t_pad), nc,
               (double *)PyArray_DATA(ap_c), &fp,
               (double *)PyArray_DATA(ap_wrk), lwrk,
               (int *)PyArray_DATA(ap_iwrk), &ier);
    } else {
        /* clocur for periodic/closed curves - no ub, ue parameters */
        clocur(iopt, ipar, idim, m,
               (double *)PyArray_DATA(ap_u), mx,
               (double *)PyArray_DATA(ap_x),
               (double *)PyArray_DATA(ap_w),
               k, s, nest, &n,
               (double *)PyArray_DATA(ap_t_pad), nc,
               (double *)PyArray_DATA(ap_c), &fp,
               (double *)PyArray_DATA(ap_wrk), lwrk,
               (int *)PyArray_DATA(ap_iwrk), &ier);
    }

    Py_DECREF(ap_x);
    Py_DECREF(ap_w);
    Py_DECREF(ap_t);

    /* Resize t and c to actual output size n */
    PyArrayObject *ap_t_out = NULL;
    PyArrayObject *ap_c_out = NULL;
    if (ier <= 3) {
        /* ier <= 0: normal return, ier=1,2,3: warnings but with valid output */
        dims[0] = n;
        ap_t_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_t_out == NULL) {
            Py_DECREF(ap_t_pad);
            goto fail_after_call;
        }
        memcpy(PyArray_DATA(ap_t_out), PyArray_DATA(ap_t_pad), n * sizeof(double));

        dims[0] = idim * (n - k - 1);
        ap_c_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_c_out == NULL) {
            Py_DECREF(ap_t_out);
            Py_DECREF(ap_t_pad);
            goto fail_after_call;
        }
        /* Copy coefficients dimension by dimension.
         * fpclos stores coefficients with spacing: c[n*(j-1) + i-1] for dimension j, coef i
         * But Python expects them packed contiguously.
         */
        double *c_data = (double *)PyArray_DATA(ap_c);
        double *c_out_data = (double *)PyArray_DATA(ap_c_out);
        int coefs_per_dim = n - k - 1;
        for (int j = 0; j < idim; j++) {
            memcpy(&c_out_data[j * coefs_per_dim],
                   &c_data[j * n],
                   coefs_per_dim * sizeof(double));
        }
    } else {
        /* Severe error - return empty arrays */
        dims[0] = 0;
        ap_t_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        ap_c_out = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        if (ap_t_out == NULL || ap_c_out == NULL) {
            Py_XDECREF(ap_t_out);
            Py_XDECREF(ap_c_out);
            Py_DECREF(ap_t_pad);
            goto fail_after_call;
        }
    }

    Py_DECREF(ap_t_pad);
    Py_DECREF(ap_c);

    /* Build output dict o = {'u': u, 'ub': ub, 'ue': ue, 'wrk': wrk, 'iwrk': iwrk, 'ier': ier, 'fp': fp} */
    PyObject *o = PyDict_New();
    if (o == NULL) {
        Py_DECREF(ap_t_out);
        Py_DECREF(ap_c_out);
        goto fail_after_call;
    }

    PyDict_SetItemString(o, "u", (PyObject *)ap_u);
    PyDict_SetItemString(o, "ub", PyFloat_FromDouble(ub));
    PyDict_SetItemString(o, "ue", PyFloat_FromDouble(ue));
    PyDict_SetItemString(o, "wrk", (PyObject *)ap_wrk);
    PyDict_SetItemString(o, "iwrk", (PyObject *)ap_iwrk);
    PyDict_SetItemString(o, "ier", PyLong_FromLong(ier));
    PyDict_SetItemString(o, "fp", PyFloat_FromDouble(fp));

    Py_DECREF(ap_u);
    Py_DECREF(ap_wrk);
    Py_DECREF(ap_iwrk);

    return Py_BuildValue(("NNO"), PyArray_Return(ap_t_out), PyArray_Return(ap_c_out), o);

fail_after_call:
    Py_XDECREF(ap_u);
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;

fail:
    Py_XDECREF(ap_u);
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    return NULL;
}

static char doc_insert[] = " [tt,cc,ier] = insert(iopt,t,c,k,x,nest)";
static PyObject *
fitpack_insert(PyObject* Py_UNUSED(dummy), PyObject *args)
{
    PyArrayObject *ap_t = NULL, *ap_c = NULL;
    PyArrayObject *ap_tt = NULL, *ap_cc = NULL;
    PyArrayObject *ap_t_pad = NULL, *ap_c_pad = NULL;
    PyObject *t, *c;
    int iopt, n, k, m, nest, nn, ier;
    double x;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "iOOidi:insert",
                          &iopt, &t, &c, &k, &x, &m)) {
        return NULL;
    }

    ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t, NPY_FLOAT64, 1, 1);
    ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c, NPY_FLOAT64, 1, 1);
    if (ap_t == NULL || ap_c == NULL) {
        goto fail;
    }

    n = PyArray_DIMS(ap_t)[0];
    /* For periodic splines (iopt != 0), need extra space for boundary conditions */
    nest = n + m + (iopt != 0 ? k : 0);

    /* Validation checks */
    if (m < 1) {
        PyErr_SetString(PyExc_ValueError, "m must be >= 1");
        goto fail;
    }

    /* Pad input arrays to size nest */
    dims[0] = nest;
    ap_t_pad = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
    ap_c_pad = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
    if (ap_t_pad == NULL || ap_c_pad == NULL) {
        goto fail;
    }
    memcpy(PyArray_DATA(ap_t_pad), PyArray_DATA(ap_t), n * sizeof(double));
    memcpy(PyArray_DATA(ap_c_pad), PyArray_DATA(ap_c), PyArray_DIMS(ap_c)[0] * sizeof(double));

    /* Allocate output arrays of size nest */
    ap_tt = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
    ap_cc = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
    if (ap_tt == NULL || ap_cc == NULL) {
        goto fail;
    }

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
    double *t_in = (double *)PyArray_DATA(ap_t_pad);
    double *c_in = (double *)PyArray_DATA(ap_c_pad);
    double *t_out = (double *)PyArray_DATA(ap_tt);
    double *c_out = (double *)PyArray_DATA(ap_cc);
    double *t_buf = NULL, *c_buf = NULL;
    double *t1, *t2, *c1, *c2, *p;

    t2 = t_in;
    c2 = c_in;
    t1 = t_out;
    c1 = c_out;

    /* Loop exactly m times to insert m knots */
    for (int iter = 0; iter < m; iter++) {
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
                    Py_DECREF(ap_t);
                    Py_DECREF(ap_c);
                    Py_DECREF(ap_t_pad);
                    Py_DECREF(ap_c_pad);
                    Py_XDECREF(ap_tt);
                    Py_XDECREF(ap_cc);
                    free(t_buf);
                    free(c_buf);
                    return NULL;
                }
            }
            t2 = t_buf;
            c2 = c_buf;
        }

        /* Use n + iter as the current array size */
        insert(iopt, t1, n + iter, c1, k, x, t2, &nn, c2, nest, &ier);

        if (ier) {
            break;
        }
    }

    /* Ensure output ends up in output buffers */
    if (t2 != t_out) {
        memcpy(t_out, t2, nest * sizeof(double));
        memcpy(c_out, c2, nest * sizeof(double));
    }

    free(t_buf);
    free(c_buf);
    Py_DECREF(ap_t);
    Py_DECREF(ap_c);
    Py_DECREF(ap_t_pad);
    Py_DECREF(ap_c_pad);

    /* Resize output arrays to actual output size nn */
    PyArrayObject *ap_tt_out = NULL;
    PyArrayObject *ap_cc_out = NULL;
    if (ier == 0) {
        dims[0] = nn;
        ap_tt_out = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
        ap_cc_out = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
        if (ap_tt_out == NULL || ap_cc_out == NULL) {
            Py_XDECREF(ap_tt);
            Py_XDECREF(ap_cc);
            Py_XDECREF(ap_tt_out);
            Py_XDECREF(ap_cc_out);
            return NULL;
        }
        memcpy(PyArray_DATA(ap_tt_out), PyArray_DATA(ap_tt), nn * sizeof(double));
        memcpy(PyArray_DATA(ap_cc_out), PyArray_DATA(ap_cc), nn * sizeof(double));
        Py_DECREF(ap_tt);
        Py_DECREF(ap_cc);
    } else {
        /* On error, return empty arrays */
        ap_tt_out = ap_tt;
        ap_cc_out = ap_cc;
    }

    /* Return arrays of actual output size nn */
    return Py_BuildValue(("NNi"),
                         PyArray_Return(ap_tt_out), PyArray_Return(ap_cc_out), ier);

fail:
    Py_XDECREF(ap_t);
    Py_XDECREF(ap_c);
    Py_XDECREF(ap_t_pad);
    Py_XDECREF(ap_c_pad);
    Py_XDECREF(ap_tt);
    Py_XDECREF(ap_cc);
    return NULL;
}


static struct PyMethodDef fitpack_module_methods[] = {
    {"bispeu",  fitpack_bispeu, METH_VARARGS, doc_bispeu},
    {"bispev",  fitpack_bispev, METH_VARARGS, doc_bispev},
    {"curfit",  fitpack_curfit, METH_VARARGS, doc_curfit},
    {"dblint",  fitpack_dblint, METH_VARARGS, doc_dblint},
    {"fpchec",  fitpack_fpchec, METH_VARARGS, doc_fpchec},
    {"insert",  fitpack_insert, METH_VARARGS, doc_insert},
    {"parder",  fitpack_parder, METH_VARARGS, doc_parder},
    {"pardtc",  fitpack_pardtc, METH_VARARGS, doc_pardtc},
    {"pardeu",  fitpack_pardeu, METH_VARARGS, doc_pardeu},
    {"percur",  fitpack_percur, METH_VARARGS, doc_percur},
    {"parcur",  fitpack_parcur, METH_VARARGS, doc_parcur},
    {"surfit",  fitpack_surfit, METH_VARARGS, doc_surfit},
    {"surfit_lsq", fitpack_surfit_lsq, METH_VARARGS, doc_surfit_lsq},
    {"surfit_smth", fitpack_surfit_smth, METH_VARARGS, doc_surfit_smth},
    {"regrid",  fitpack_regrid, METH_VARARGS, doc_regrid},
    {"spalde",  fitpack_spalde, METH_VARARGS, doc_spalde},
    {"spgrid",  fitpack_spgrid, METH_VARARGS, doc_spgrid},
    {"sphere",  fitpack_sphere, METH_VARARGS, doc_sphere},
    {"splder",  fitpack_splder, METH_VARARGS, doc_splder},
    {"splev",   fitpack_splev,  METH_VARARGS, doc_splev},
    {"splint",  fitpack_splint, METH_VARARGS, doc_splint},
    {"sproot",  fitpack_sproot, METH_VARARGS, doc_sproot},
    {NULL, NULL, 0, NULL}
};


static int fitpack_module_exec(PyObject *module) {
    if (_import_array() < 0) { return -1; }

    fitpack_error = PyErr_NewException("scipy.interpolate._fitpack.error", NULL, NULL);
    if (fitpack_error == NULL) {
        return -1;
    }
    Py_INCREF(fitpack_error);
    if (PyModule_AddObject(module, "error", fitpack_error) < 0) {
        Py_DECREF(fitpack_error);
        return -1;
    }

    return 0;
}


static struct PyModuleDef_Slot fitpack_module_slots[] = {
    {Py_mod_exec, fitpack_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_fitpack",
    .m_doc = "C-translations of algorithms from Fortran77 FITPACK collection",
    .m_size = 0,
    .m_methods = fitpack_module_methods,
    .m_slots = fitpack_module_slots
};


PyMODINIT_FUNC PyInit__fitpack(void) {
    return PyModuleDef_Init(&moduledef);
}
