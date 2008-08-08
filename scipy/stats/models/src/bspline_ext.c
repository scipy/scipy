#include "Python.h"
#include "numpy/arrayobject.h"

/*  function prototypes */

double *bspline(double**, double*, int, double *, int, int, int, int, int); 
void bspline_gram(double **, double *, int, int, int, int);
void invband_compute(double **, double *, int, int);


static PyObject *BSpline_Invband(PyObject *self, PyObject *args)
{

    double *data;
    double *L_data;
    npy_intp *dims_invband;
    npy_intp *dims_L;
    PyArrayObject *L       = NULL;
    PyArrayObject *invband = NULL;

    if(!PyArg_ParseTuple(args, "O", &L)) 
	    goto exit;

    dims_L = PyArray_DIMS(L);
    L_data = (double *)PyArray_DATA(L);

    dims_invband = calloc(2, sizeof(npy_intp));
    dims_invband[0] = dims_L[0];
    dims_invband[1] = dims_L[1];

    invband = (PyArrayObject*)PyArray_SimpleNew(2, dims_invband, PyArray_DOUBLE);
    data    = (double *)PyArray_DATA(invband);
    free(dims_invband);

    invband_compute(&data, L_data, (int)dims_L[0], (int)dims_L[1]);

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("O", invband);

}



static PyObject *BSpline_Gram(PyObject *self, PyObject *args)
{

    int m;
    int dl;
    int dr;
    double *knots;
    double *data;
    npy_intp *nknots;
    npy_intp *dims_gram;
    PyArrayObject *knots_array = NULL;
    PyArrayObject *gram_array  = NULL;

    if(!PyArg_ParseTuple(args, "Oiii", &knots_array, &m, &dl, &dr)) 
	    goto exit;

    nknots = PyArray_DIMS(knots_array);
    knots  = (double *)PyArray_DATA(knots_array);

    dims_gram = calloc(2, sizeof(npy_intp));
    dims_gram[0] = (int)nknots[0] - m;
    dims_gram[1] = m; 

    gram_array  = (PyArrayObject*)PyArray_SimpleNew(2, dims_gram, PyArray_DOUBLE);
    data        = (double *)PyArray_DATA(gram_array);
    free(dims_gram);

    bspline_gram(&data, knots, (int)nknots[0], m, dl, dr);

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("O", gram_array);

}


static PyObject *BSpline_Evaluate(PyObject *self, PyObject *args)
{

    int i;
    int upper;
    int lower;
    int m;
    int d;
    double *knots;
    double *x;
    double *data;
    npy_intp *nknots;
    npy_intp *nx;
    npy_intp dims_basis[2];
    PyArrayObject *knots_array = NULL;
    PyArrayObject *x_array     = NULL;
    PyArrayObject *basis_array = NULL;

    if(!PyArg_ParseTuple(args, "OOiiii", &x_array, &knots_array, &m, &d, &lower, &upper)) 
	    goto exit;

    nknots = PyArray_DIMS(knots_array);
    nx     = PyArray_DIMS(x_array);

    knots  = (double *)PyArray_DATA(knots_array);
    x      = (double *)PyArray_DATA(x_array);

    dims_basis[0] = upper-lower;
    dims_basis[1] = (int)nx[0];
    basis_array   = (PyArrayObject*)PyArray_SimpleNew(2, dims_basis, PyArray_DOUBLE);
    data          = (double *)PyArray_DATA(basis_array);

    bspline(&data, x, (int)nx[0], knots, (int)nknots[0], m, d, lower, upper); 

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("O", basis_array);

}


static PyMethodDef BSplineMethods[] =
{
    { "evaluate",     BSpline_Evaluate,  METH_VARARGS, NULL },
    { "gram",         BSpline_Gram,      METH_VARARGS, NULL },
    { "invband",      BSpline_Invband,   METH_VARARGS, NULL },
    {  NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC init_hbspline(void)
{
    Py_InitModule("_hbspline", BSplineMethods);
    import_array();
}

