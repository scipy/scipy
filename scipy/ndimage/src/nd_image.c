/* Copyright (C) 2003-2005 Peter J. Verveer
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 *    products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * The order of these first two imports should not be changed, see the note
 * in ni_support.h for details.
 */
#include "nd_image.h"
#include "ni_support.h"

#include "ni_filters.h"
#include "ni_fourier.h"
#include "ni_morphology.h"
#include "ni_interpolation.h"
#include "ni_measure.h"

#include "ccallback.h"


typedef struct {
    PyObject *extra_arguments;
    PyObject *extra_keywords;
} NI_PythonCallbackData;

/* Numarray Helper Functions */

/*
 * Creates a new numpy array of the requested type and shape, and either
 * copies into it the contents of buffer, or sets it to all zeros if
 * buffer is NULL.
 */
static PyArrayObject *
NA_NewArray(void *buffer, enum NPY_TYPES type, int ndim, npy_intp *shape)
{
    PyArrayObject *result;

    if (type == NPY_NOTYPE) {
        type = NPY_DOUBLE;
    }

    result = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, type);
    if (result == NULL) {
        return NULL;
    }

    if (buffer == NULL) {
        memset(PyArray_DATA(result), 0, PyArray_NBYTES(result));
    }
    else {
        memcpy(PyArray_DATA(result), buffer, PyArray_NBYTES(result));
    }

    return result;
}

/* Converts a Python array-like object into a behaved input array. */
static int
NI_ObjectToInputArray(PyObject *object, PyArrayObject **array)
{
    int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_NOTSWAPPED;
    *array = (PyArrayObject *)PyArray_CheckFromAny(object, NULL, 0, 0, flags,
                                                   NULL);
    return *array != NULL;
}

/* Like NI_ObjectToInputArray, but with special handling for Py_None. */
static int
NI_ObjectToOptionalInputArray(PyObject *object, PyArrayObject **array)
{
    if (object == Py_None) {
        *array = NULL;
        return 1;
    }
    return NI_ObjectToInputArray(object, array);
}

/* Converts a Python array-like object into a behaved output array. */
static int
NI_ObjectToOutputArray(PyObject *object, PyArrayObject **array)
{
    int flags = NPY_ARRAY_BEHAVED_NS | NPY_ARRAY_WRITEBACKIFCOPY;
    /*
     * This would also be caught by the PyArray_CheckFromAny call, but
     * we check it explicitly here to provide a saner error message.
     */
    if (PyArray_Check(object) &&
            !PyArray_ISWRITEABLE((PyArrayObject *)object)) {
        PyErr_SetString(PyExc_ValueError, "output array is read-only.");
    return 0;
    }
    /*
     * If the input array is not aligned or is byteswapped, this call
     * will create a new aligned, native byte order array, and copy the
     * contents of object into it. For an output array, the copy is
     * unnecessary, so this could be optimized. It is very easy to not
     * do NPY_ARRAY_UPDATEIFCOPY right, so we let NumPy do it for us
     * and pay the performance price.
     */
    *array = (PyArrayObject *)PyArray_CheckFromAny(object, NULL, 0, 0, flags,
                                                   NULL);
    return *array != NULL;
}

/* Like NI_ObjectToOutputArray, but with special handling for Py_None. */
static int
NI_ObjectToOptionalOutputArray(PyObject *object, PyArrayObject **array)
{
    if (object == Py_None) {
        *array = NULL;
        return 1;
    }
    return NI_ObjectToOutputArray(object, array);
}

/* Converts a Python array-like object into a behaved input/output array. */
static int
NI_ObjectToInputOutputArray(PyObject *object, PyArrayObject **array)
{
    /*
     * This is also done in NI_ObjectToOutputArray, double checking here
     * to provide a more specific error message.
     */
    if (PyArray_Check(object) &&
            !PyArray_ISWRITEABLE((PyArrayObject *)object)) {
        PyErr_SetString(PyExc_ValueError, "input/output array is read-only.");
        return 0;
    }
    return NI_ObjectToOutputArray(object, array);
}

/* Checks that an origin value was received for each array dimension. */
static int
_validate_origin(PyArrayObject *array, PyArray_Dims origin)
{
    if (origin.len != PyArray_NDIM(array)) {
        PyErr_Format(PyExc_ValueError,
                     "Invalid %d element 'origin' sequence for "
                     "%d-dimensional input array.",
                     origin.len, PyArray_NDIM(array));
        return 0;
    }
    return 1;
}

/*********************************************************************/
/* wrapper functions: */
/*********************************************************************/

static PyObject *Py_Correlate1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *weights = NULL;
    int axis, mode;
    double cval;
    npy_intp origin;

    if (!PyArg_ParseTuple(args, "O&O&iO&idn" ,
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &weights, &axis,
                          NI_ObjectToOutputArray, &output, &mode, &cval,
                          &origin))
        goto exit;

    NI_Correlate1D(input, weights, axis, output, (NI_ExtendMode)mode, cval,
                   origin);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(weights);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_Correlate(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *weights = NULL;
    PyArray_Dims origin = {NULL, 0};
    int mode;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&idO&", NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &weights,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval,
                          PyArray_IntpConverter, &origin)) {
        goto exit;
    }
    if (!_validate_origin(input, origin)) {
        goto exit;
    }

    NI_Correlate(input, weights, output, (NI_ExtendMode)mode, cval,
                 origin.ptr);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(weights);
    Py_XDECREF(output);
    PyDimMem_FREE(origin.ptr);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_UniformFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, mode;
    npy_intp filter_size, origin;
    double cval;

    if (!PyArg_ParseTuple(args, "O&niO&idn",
                          NI_ObjectToInputArray, &input,
                          &filter_size, &axis,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval, &origin))
        goto exit;

    NI_UniformFilter1D(input, filter_size, axis, output, (NI_ExtendMode)mode,
                       cval, origin);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_MinOrMaxFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, mode, minimum;
    npy_intp filter_size, origin;
    double cval;

    if (!PyArg_ParseTuple(args, "O&niO&idni",
                          NI_ObjectToInputArray, &input,
                          &filter_size, &axis,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval, &origin, &minimum))
        goto exit;

    NI_MinOrMaxFilter1D(input, filter_size, axis, output, (NI_ExtendMode)mode,
                        cval, origin, minimum);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_MinOrMaxFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *footprint = NULL;
    PyArrayObject *structure = NULL;
    PyArray_Dims origin = {NULL, 0};
    int mode, minimum;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&O&idO&i",
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &footprint,
                          NI_ObjectToOptionalInputArray, &structure,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval,
                          PyArray_IntpConverter, &origin,
                          &minimum)) {
        goto exit;
    }
    if (!_validate_origin(input, origin)) {
        goto exit;
    }

    NI_MinOrMaxFilter(input, footprint, structure, output, (NI_ExtendMode)mode,
                      cval, origin.ptr, minimum);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(footprint);
    Py_XDECREF(structure);
    Py_XDECREF(output);
    PyDimMem_FREE(origin.ptr);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_RankFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *footprint = NULL;
    PyArray_Dims origin = {NULL, 0};
    int mode, rank;
    double cval;

    if (!PyArg_ParseTuple(args, "O&iO&O&idO&",
                          NI_ObjectToInputArray, &input, &rank,
                          NI_ObjectToInputArray, &footprint,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval,
                          PyArray_IntpConverter, &origin)) {
        goto exit;
    }
    if (!_validate_origin(input, origin)) {
        goto exit;
    }

    NI_RankFilter(input, rank, footprint, output, (NI_ExtendMode)mode, cval,
                  origin.ptr);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(footprint);
    Py_XDECREF(output);
    PyDimMem_FREE(origin.ptr);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_Filter1DFunc(double *iline, npy_intp ilen,
                           double *oline, npy_intp olen, void *data)
{
    PyArrayObject *py_ibuffer = NULL, *py_obuffer = NULL;
    PyObject *rv = NULL, *args = NULL, *tmp = NULL;
    npy_intp ii;
    double *po = NULL;
    ccallback_t *callback = (ccallback_t *)data;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)callback->info_p;

    py_ibuffer = NA_NewArray(iline, NPY_DOUBLE, 1, &ilen);
    py_obuffer = NA_NewArray(NULL, NPY_DOUBLE, 1, &olen);
    if (!py_ibuffer || !py_obuffer)
        goto exit;
    tmp = Py_BuildValue("(OO)", py_ibuffer, py_obuffer);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rv = PyObject_Call(callback->py_function, args, cbdata->extra_keywords);
    if (!rv)
        goto exit;
    po = (double*)PyArray_DATA(py_obuffer);
    for(ii = 0; ii < olen; ii++)
        oline[ii] = po[ii];
exit:
    Py_XDECREF(py_ibuffer);
    Py_XDECREF(py_obuffer);
    Py_XDECREF(rv);
    Py_XDECREF(args);
    Py_XDECREF(tmp);
    return PyErr_Occurred() ? 0 : 1;
}

static PyObject *Py_GenericFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    PyObject *fnc = NULL, *extra_arguments = NULL, *extra_keywords = NULL;
    void *func = NULL, *data = NULL;
    NI_PythonCallbackData cbdata;
    int axis, mode;
    npy_intp origin, filter_size;
    double cval;
    ccallback_t callback;
    static ccallback_signature_t callback_signatures[] = {
        {"int (double *, intptr_t, double *, intptr_t, void *)"},
        {"int (double *, npy_intp, double *, npy_intp, void *)"},
#if NPY_SIZEOF_INTP == NPY_SIZEOF_SHORT
        {"int (double *, short, double *, short, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_INT
        {"int (double *, int, double *, int, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONG
        {"int (double *, long, double *, long, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONGLONG
        {"int (double *, long long, double *, long long, void *)"},
#endif
        {NULL}
    };

    callback.py_function = NULL;
    callback.c_function = NULL;

    if (!PyArg_ParseTuple(args, "O&OniO&idnOO",
                          NI_ObjectToInputArray, &input,
                          &fnc, &filter_size, &axis,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval, &origin,
                          &extra_arguments, &extra_keywords))
        goto exit;

    if (!PyTuple_Check(extra_arguments)) {
        PyErr_SetString(PyExc_RuntimeError, "extra_arguments must be a tuple");
        goto exit;
    }
    if (!PyDict_Check(extra_keywords)) {
        PyErr_SetString(PyExc_RuntimeError,
                                        "extra_keywords must be a dictionary");
        goto exit;
    }
    if (PyCapsule_CheckExact(fnc) && PyCapsule_GetName(fnc) == NULL) {
        /* 'Legacy' low-level callable */
        func = PyCapsule_GetPointer(fnc, NULL);
        data = PyCapsule_GetContext(fnc);
    } else {
        int ret;

        ret = ccallback_prepare(&callback, callback_signatures, fnc, CCALLBACK_DEFAULTS);
        if (ret == -1) {
            goto exit;
        }

        if (callback.py_function != NULL) {
            cbdata.extra_arguments = extra_arguments;
            cbdata.extra_keywords = extra_keywords;
            callback.info_p = (void*)&cbdata;
            func = Py_Filter1DFunc;
            data = (void*)&callback;
        }
        else {
            func = callback.c_function;
            data = callback.user_data;
        }
    }

    NI_GenericFilter1D(input, func, data, filter_size, axis, output,
                       (NI_ExtendMode)mode, cval, origin);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    if (callback.py_function != NULL || callback.c_function != NULL) {
        ccallback_release(&callback);
    }
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_FilterFunc(double *buffer, npy_intp filter_size,
                                                 double *output, void *data)
{
    PyArrayObject *py_buffer = NULL;
    PyObject *rv = NULL, *args = NULL, *tmp = NULL;
    ccallback_t *callback = (ccallback_t *)data;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)callback->info_p;

    py_buffer = NA_NewArray(buffer, NPY_DOUBLE, 1, &filter_size);
    if (!py_buffer)
        goto exit;
    tmp = Py_BuildValue("(O)", py_buffer);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rv = PyObject_Call(callback->py_function, args, cbdata->extra_keywords);
    if (!rv)
        goto exit;
    *output = PyFloat_AsDouble(rv);
exit:
    Py_XDECREF(py_buffer);
    Py_XDECREF(rv);
    Py_XDECREF(args);
    Py_XDECREF(tmp);
    return PyErr_Occurred() ? 0 : 1;
}

static PyObject *Py_GenericFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *footprint = NULL;
    PyObject *fnc = NULL, *extra_arguments = NULL, *extra_keywords = NULL;
    void *func = NULL, *data = NULL;
    NI_PythonCallbackData cbdata;
    int mode;
    PyArray_Dims origin = {NULL, 0};
    double cval;
    ccallback_t callback;
    static ccallback_signature_t callback_signatures[] = {
        {"int (double *, intptr_t, double *, void *)"},
        {"int (double *, npy_intp, double *, void *)"},
#if NPY_SIZEOF_INTP == NPY_SIZEOF_SHORT
        {"int (double *, short, double *, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_INT
        {"int (double *, int, double *, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONG
        {"int (double *, long, double *, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONGLONG
        {"int (double *, long long, double *, void *)"},
#endif
        {NULL}
    };

    callback.py_function = NULL;
    callback.c_function = NULL;

    if (!PyArg_ParseTuple(args, "O&OO&O&idO&OO",
                          NI_ObjectToInputArray, &input,
                          &fnc,
                          NI_ObjectToInputArray, &footprint,
                          NI_ObjectToOutputArray, &output,
                          &mode, &cval,
                          PyArray_IntpConverter, &origin,
                          &extra_arguments, &extra_keywords)) {
        goto exit;
    }
    if (!_validate_origin(input, origin)) {
        goto exit;
    }
    if (!PyTuple_Check(extra_arguments)) {
        PyErr_SetString(PyExc_RuntimeError, "extra_arguments must be a tuple");
        goto exit;
    }
    if (!PyDict_Check(extra_keywords)) {
        PyErr_SetString(PyExc_RuntimeError,
                                        "extra_keywords must be a dictionary");
        goto exit;
    }
    if (PyCapsule_CheckExact(fnc) && PyCapsule_GetName(fnc) == NULL) {
        func = PyCapsule_GetPointer(fnc, NULL);
        data = PyCapsule_GetContext(fnc);
    } else {
        int ret;

        ret = ccallback_prepare(&callback, callback_signatures, fnc, CCALLBACK_DEFAULTS);
        if (ret == -1) {
            goto exit;
        }

        if (callback.py_function != NULL) {
            cbdata.extra_arguments = extra_arguments;
            cbdata.extra_keywords = extra_keywords;
            callback.info_p = (void*)&cbdata;
            func = Py_FilterFunc;
            data = (void*)&callback;
        }
        else {
            func = callback.c_function;
            data = callback.user_data;
        }
    }

    NI_GenericFilter(input, func, data, footprint, output, (NI_ExtendMode)mode,
                     cval, origin.ptr);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    if (callback.py_function != NULL || callback.c_function != NULL) {
        ccallback_release(&callback);
    }
    Py_XDECREF(input);
    Py_XDECREF(output);
    Py_XDECREF(footprint);
    PyDimMem_FREE(origin.ptr);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_FourierFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *parameters = NULL;
    int axis, filter_type;
    npy_intp n;

    if (!PyArg_ParseTuple(args, "O&O&niO&i",
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &parameters,
                          &n, &axis,
                          NI_ObjectToOutputArray, &output,
                          &filter_type))
        goto exit;

    NI_FourierFilter(input, parameters, n, axis, output, filter_type);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(parameters);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_FourierShift(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *shifts = NULL;
    int axis;
    npy_intp n;

    if (!PyArg_ParseTuple(args, "O&O&niO&",
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &shifts,
                          &n, &axis,
                          NI_ObjectToOutputArray, &output))
        goto exit;

    NI_FourierShift(input, shifts, n, axis, output);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(shifts);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_SplineFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, order, mode;

    if (!PyArg_ParseTuple(args, "O&iiO&i",
                          NI_ObjectToInputArray, &input, &order, &axis,
                          NI_ObjectToOutputArray, &output, &mode))
        goto exit;

    NI_SplineFilter1D(input, order, axis, mode, output);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_Map(npy_intp *ocoor, double* icoor, int orank, int irank,
                                    void *data)
{
    PyObject *coors = NULL, *rets = NULL, *args = NULL, *tmp = NULL;
    npy_intp ii;
    ccallback_t *callback = (ccallback_t *)data;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)callback->info_p;

    coors = PyTuple_New(orank);
    if (!coors)
        goto exit;
    for(ii = 0; ii < orank; ii++) {
        PyTuple_SetItem(coors, ii, PyLong_FromSsize_t(ocoor[ii]));
        if (PyErr_Occurred())
            goto exit;
    }
    tmp = Py_BuildValue("(O)", coors);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rets = PyObject_Call(callback->py_function, args, cbdata->extra_keywords);
    if (!rets)
        goto exit;
    for(ii = 0; ii < irank; ii++) {
        icoor[ii] = PyFloat_AsDouble(PyTuple_GetItem(rets, ii));
        if (PyErr_Occurred())
            goto exit;
    }
exit:
    Py_XDECREF(coors);
    Py_XDECREF(tmp);
    Py_XDECREF(rets);
    Py_XDECREF(args);
    return PyErr_Occurred() ? 0 : 1;
}


static PyObject *Py_GeometricTransform(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    PyArrayObject *coordinates = NULL, *matrix = NULL, *shift = NULL;
    PyObject *fnc = NULL, *extra_arguments = NULL, *extra_keywords = NULL;
    int mode, order, nprepad;
    double cval;
    void *func = NULL, *data = NULL;
    NI_PythonCallbackData cbdata;
    ccallback_t callback;
    static ccallback_signature_t callback_signatures[] = {
        {"int (intptr_t *, double *, int, int, void *)"},
        {"int (npy_intp *, double *, int, int, void *)"},
#if NPY_SIZEOF_INTP == NPY_SIZEOF_SHORT
        {"int (short *, double *, int, int, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_INT
        {"int (int *, double *, int, int, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONG
        {"int (long *, double *, int, int, void *)"},
#endif
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONGLONG
        {"int (long long *, double *, int, int, void *)"},
#endif
        {NULL}
    };

    callback.py_function = NULL;
    callback.c_function = NULL;

    if (!PyArg_ParseTuple(args, "O&OO&O&O&O&iidiOO",
                          NI_ObjectToInputArray, &input,
                          &fnc,
                          NI_ObjectToOptionalInputArray, &coordinates,
                          NI_ObjectToOptionalInputArray, &matrix,
                          NI_ObjectToOptionalInputArray, &shift,
                          NI_ObjectToOutputArray, &output,
                          &order, &mode, &cval, &nprepad,
                          &extra_arguments, &extra_keywords))
        goto exit;

    if (fnc != Py_None) {
        if (!PyTuple_Check(extra_arguments)) {
            PyErr_SetString(PyExc_RuntimeError,
                                            "extra_arguments must be a tuple");
            goto exit;
        }
        if (!PyDict_Check(extra_keywords)) {
            PyErr_SetString(PyExc_RuntimeError,
                                            "extra_keywords must be a dictionary");
            goto exit;
        }
        if (PyCapsule_CheckExact(fnc) && PyCapsule_GetName(fnc) == NULL) {
            func = PyCapsule_GetPointer(fnc, NULL);
            data = PyCapsule_GetContext(fnc);
        } else {
            int ret;

            ret = ccallback_prepare(&callback, callback_signatures, fnc, CCALLBACK_DEFAULTS);
            if (ret == -1) {
                goto exit;
            }

            if (callback.py_function != NULL) {
                cbdata.extra_arguments = extra_arguments;
                cbdata.extra_keywords = extra_keywords;
                callback.info_p = (void*)&cbdata;
                func = Py_Map;
                data = (void*)&callback;
            }
            else {
                func = callback.c_function;
                data = callback.user_data;
            }
        }
    }

    NI_GeometricTransform(input, func, data, matrix, shift, coordinates,
                          output, order, (NI_ExtendMode)mode, cval, nprepad);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    if (callback.py_function != NULL || callback.c_function != NULL) {
        ccallback_release(&callback);
    }
    Py_XDECREF(input);
    Py_XDECREF(output);
    Py_XDECREF(coordinates);
    Py_XDECREF(matrix);
    Py_XDECREF(shift);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_ZoomShift(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *shift = NULL;
    PyArrayObject *zoom = NULL;
    int mode, order, nprepad, grid_mode;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&O&iidii",
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToOptionalInputArray, &zoom,
                          NI_ObjectToOptionalInputArray, &shift,
                          NI_ObjectToOutputArray, &output,
                          &order, &mode, &cval, &nprepad, &grid_mode))
        goto exit;

    NI_ZoomShift(input, zoom, shift, output, order, (NI_ExtendMode)mode, cval,
                 nprepad, grid_mode);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(shift);
    Py_XDECREF(zoom);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_FindObjects(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL;
    PyObject *result = NULL, *tuple = NULL, *start = NULL, *end = NULL;
    PyObject *slc = NULL;
    int jj;
    npy_intp max_label;
    npy_intp ii, *regions = NULL;

    if (!PyArg_ParseTuple(args, "O&n",
                          NI_ObjectToInputArray, &input, &max_label))
        goto exit;

    if (max_label < 0)
        max_label = 0;
    if (max_label > 0) {
        if (PyArray_NDIM(input) > 0) {
            regions = (npy_intp*)malloc(2 * max_label * PyArray_NDIM(input) *
                                        sizeof(npy_intp));
        } else {
            regions = (npy_intp*)malloc(max_label * sizeof(npy_intp));
        }
        if (!regions) {
            PyErr_NoMemory();
            goto exit;
        }
    }

    if (!NI_FindObjects(input, max_label, regions))
        goto exit;

    result = PyList_New(max_label);
    if (!result) {
        PyErr_NoMemory();
        goto exit;
    }

    for(ii = 0; ii < max_label; ii++) {
        npy_intp idx =
                PyArray_NDIM(input) > 0 ? 2 * PyArray_NDIM(input) * ii : ii;
        if (regions[idx] >= 0) {
            tuple = PyTuple_New(PyArray_NDIM(input));
            if (!tuple) {
                PyErr_NoMemory();
                goto exit;
            }
            for(jj = 0; jj < PyArray_NDIM(input); jj++) {
                start = PyLong_FromSsize_t(regions[idx + jj]);
                end = PyLong_FromSsize_t(regions[idx + jj +
                                             PyArray_NDIM(input)]);
                if (!start || !end) {
                    PyErr_NoMemory();
                    goto exit;
                }
                slc = PySlice_New(start, end, NULL);
                if (!slc) {
                    PyErr_NoMemory();
                    goto exit;
                }
                Py_DECREF(start);
                Py_DECREF(end);
                start = end = NULL;
                PyTuple_SetItem(tuple, jj, slc);
                slc = NULL;
            }
            PyList_SetItem(result, ii, tuple);
            tuple = NULL;
        } else {
            Py_INCREF(Py_None);
            PyList_SetItem(result, ii, Py_None);
        }
    }

    Py_INCREF(result);

 exit:
    Py_XDECREF(input);
    Py_XDECREF(result);
    Py_XDECREF(tuple);
    Py_XDECREF(start);
    Py_XDECREF(end);
    Py_XDECREF(slc);
    free(regions);
    if (PyErr_Occurred()) {
        return NULL;
    } else {
        return result;
    }
}

/*
   Implement the ndimage.value_indices() function.
   Makes 3 passes through the array data. We use ndimage's NI_Iterator to
   iterate through all elements of the input array.

   To support all numpy datatypes, we define several macros which take a
   parameter for the array datatype in question. By using these, all comparisons
   are carried out using the same datatype as the input array.
*/
#define VALUEINDICES_MINVAL(valType) (*((valType *)PyArray_GETPTR1(minMaxArr, 0)))
#define VALUEINDICES_MAXVAL(valType) (*((valType *)PyArray_GETPTR1(minMaxArr, 1)))
#define VALUEINDICES_IGNOREVAL(valType) (*((valType *)PyArray_GETPTR1(ignorevalArr, 0)))
#define CASE_VALUEINDICES_SET_MINMAX(valType) {\
    valType val = *((valType*)arrData); \
    if (ignoreValIsNone || (val != VALUEINDICES_IGNOREVAL(valType))) {\
        if (minMaxUnset) { \
            VALUEINDICES_MINVAL(valType) = val;  \
            VALUEINDICES_MAXVAL(valType) = val; \
            minMaxUnset = 0; \
        } else { \
            if (val < VALUEINDICES_MINVAL(valType)) VALUEINDICES_MINVAL(valType) = val; \
            if (val > VALUEINDICES_MAXVAL(valType)) VALUEINDICES_MAXVAL(valType) = val; \
        } \
    }\
}
#define CASE_VALUEINDICES_MAKEHISTOGRAM(valType) {\
    numPossibleVals = (VALUEINDICES_MAXVAL(valType) - VALUEINDICES_MINVAL(valType) + 1); \
    hist = (npy_intp *)calloc(numPossibleVals, sizeof(npy_intp)); \
    if (hist != NULL) { \
        NI_InitPointIterator(arr, &ndiIter); \
        arrData = (char *)PyArray_DATA(arr); \
        for (iterIndex=0; iterIndex<arrSize; iterIndex++) { \
            valType val = *((valType*)arrData); \
            if (ignoreValIsNone || (val != VALUEINDICES_IGNOREVAL(valType))) { \
                ii = val - VALUEINDICES_MINVAL(valType); \
                hist[ii] += 1; \
            } \
            NI_ITERATOR_NEXT(ndiIter, arrData); \
        } \
    } else { \
        PyErr_SetString(PyExc_MemoryError, "Couldn't allocate hist"); \
    } \
}
#define CASE_VALUEINDICES_GET_VALUEOFFSET(valType) { \
    valType val = *((valType*)arrData); \
    ii = val - VALUEINDICES_MINVAL(valType); \
    valueIsIgnore = (val == VALUEINDICES_IGNOREVAL(valType)); \
}
#define CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(valType, ii) { \
    valType val = ii + VALUEINDICES_MINVAL(valType); \
    valObj = PyArray_ToScalar(&val, minMaxArr); \
}
static PyObject *NI_ValueIndices(PyObject *self, PyObject *args)
{
    PyArrayObject *arr, *ndxArr, *minMaxArr, *ignorevalArr;
    PyObject *t=NULL, *valObj=NULL, **ndxPtr=NULL, *ndxTuple, *valDict;
    int ignoreValIsNone, valueIsIgnore=0, ndim, j, arrType, minMaxUnset=1;
    NI_Iterator ndiIter;
    char *arrData;
    npy_intp *hist=NULL, *valCtr=NULL, ii, numPossibleVals=0;
    npy_intp arrSize, iterIndex, dims[1];

    /* Get arguments passed in */
    if (!PyArg_ParseTuple(args, "O!iO!", &PyArray_Type, &arr, &ignoreValIsNone,
            &PyArray_Type, &ignorevalArr))
        return NULL;

    arrSize = PyArray_SIZE(arr);
    arrType = PyArray_TYPE(arr);
    arrData = (char *)PyArray_DATA(arr);
    ndim = PyArray_NDIM(arr);
    if (!PyTypeNum_ISINTEGER(arrType)) {
        PyErr_SetString(PyExc_ValueError, "Parameter 'arr' must be an integer array");
        return NULL;
    }

    /* This dictionary is the final return value */
    valDict = PyDict_New();
    if (valDict == NULL) return PyErr_NoMemory();

    /* We use a small numpy array for the min and max values, as this will
       take the same datatype as the input array */
    dims[0] = 2;
    minMaxArr = (PyArrayObject *)PyArray_SimpleNew(1, dims, arrType);
    if (minMaxArr == NULL) return PyErr_NoMemory();

    /* First pass. Just set the min & max */
    NI_InitPointIterator(arr, &ndiIter);
    for (iterIndex=0; iterIndex<arrSize; iterIndex++) {
        switch(arrType) {
        case NPY_INT8:   CASE_VALUEINDICES_SET_MINMAX(npy_int8); break;
        case NPY_UINT8:  CASE_VALUEINDICES_SET_MINMAX(npy_uint8); break;
        case NPY_INT16:  CASE_VALUEINDICES_SET_MINMAX(npy_int16); break;
        case NPY_UINT16: CASE_VALUEINDICES_SET_MINMAX(npy_uint16); break;
        case NPY_INT32:  CASE_VALUEINDICES_SET_MINMAX(npy_int32); break;
        case NPY_UINT32: CASE_VALUEINDICES_SET_MINMAX(npy_uint32); break;
        case NPY_INT64:  CASE_VALUEINDICES_SET_MINMAX(npy_int64); break;
        case NPY_UINT64: CASE_VALUEINDICES_SET_MINMAX(npy_uint64); break;
        default:
            switch(arrType) {
            case NPY_UINT: CASE_VALUEINDICES_SET_MINMAX(npy_uint); break;
            case NPY_INT: CASE_VALUEINDICES_SET_MINMAX(npy_int); break;
            }
        }
        NI_ITERATOR_NEXT(ndiIter, arrData);
    }

    /* Second pass, creates a histogram of all the possible values between
       min and max. If min/max were not set, then the array was all ignore
       values. */
    if (!minMaxUnset) {
        switch(arrType) {
        case NPY_INT8:   CASE_VALUEINDICES_MAKEHISTOGRAM(npy_int8); break;
        case NPY_UINT8:  CASE_VALUEINDICES_MAKEHISTOGRAM(npy_uint8); break;
        case NPY_INT16:  CASE_VALUEINDICES_MAKEHISTOGRAM(npy_int16); break;
        case NPY_UINT16: CASE_VALUEINDICES_MAKEHISTOGRAM(npy_uint16); break;
        case NPY_INT32:  CASE_VALUEINDICES_MAKEHISTOGRAM(npy_int32); break;
        case NPY_UINT32: CASE_VALUEINDICES_MAKEHISTOGRAM(npy_uint32); break;
        case NPY_INT64:  CASE_VALUEINDICES_MAKEHISTOGRAM(npy_int64); break;
        case NPY_UINT64: CASE_VALUEINDICES_MAKEHISTOGRAM(npy_uint64); break;
        default:
            switch(arrType) {
            case NPY_INT:  CASE_VALUEINDICES_MAKEHISTOGRAM(npy_int); break;
            case NPY_UINT: CASE_VALUEINDICES_MAKEHISTOGRAM(npy_uint); break;
            }
        }
    }

    if (hist != NULL) {
        /* Allocate local data structures to track where we are up to while
           assigning index values */
        valCtr = (npy_intp *)calloc(numPossibleVals, sizeof(npy_intp));
        ndxPtr = (PyObject **)calloc(numPossibleVals, sizeof(PyObject  *));
        if (valCtr == NULL)
            PyErr_SetString(PyExc_MemoryError, "Couldn't allocate valCtr");
        else if (ndxPtr == NULL)
            PyErr_SetString(PyExc_MemoryError, "Couldn't allocate ndxPtr");
    }

    if ((valCtr != NULL) && (ndxPtr != NULL)) {
        /* For each value with non-zero histogram count, allocate an array index
           tuple of arrays of shape (count,), in the corresponding entry of
           the dictionary */
        for (ii=0; ii<numPossibleVals; ii++) {
            if (hist[ii] > 0) {
                /* Create a numpy scalar <valObj> for the value corresponding to
                   offset <ii> */
                switch(arrType) {
                case NPY_INT8:   CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_int8, ii); break;
                case NPY_UINT8:  CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_uint8, ii); break;
                case NPY_INT16:  CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_int16, ii); break;
                case NPY_UINT16: CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_uint16, ii); break;
                case NPY_INT32:  CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_int32, ii); break;
                case NPY_UINT32: CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_uint32, ii); break;
                case NPY_INT64:  CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_int64, ii); break;
                case NPY_UINT64: CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_uint64, ii); break;
                default:
                    switch(arrType) {
                    case NPY_INT:  CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_int, ii); break;
                    case NPY_UINT: CASE_VALUEINDICES_MAKE_VALUEOBJ_FROMOFFSET(npy_uint, ii); break;
                    }
                }
                /* Create a tuple of <ndim> index arrays */
                t = PyTuple_New(ndim);
                if ((t != NULL) && (valObj != NULL)) {
                    for (j=0; j<ndim; j++) {
                        dims[0] = hist[ii];
                        ndxArr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INTP);
                        if (ndxArr == NULL) {
                            PyErr_SetString(PyExc_MemoryError, "Couldn't allocate ndxArr");
                        } else {
                            PyTuple_SetItem(t, j, (PyObject *)ndxArr);
                        }
                    }
                    /* Set the tuple as the dict entry for the value */
                    PyDict_SetItem(valDict, valObj, t);
                }
                if (valObj != NULL) Py_DECREF(valObj);
                if (t != NULL) Py_DECREF(t);
                /* Save a pointer to the tuple, to avoid later dict lookups */
                ndxPtr[ii] = t;

                if ((t == NULL) || (valObj == NULL)) {
                    PyErr_SetString(PyExc_MemoryError, "Couldn't allocate dict entry");
                    break;
                }
            }
        }

        if (!PyErr_Occurred()) {
            /* Third pass. Loop over all array elements, adding indices in the
               right place. */
            NI_InitPointIterator(arr, &ndiIter);
            arrData = (char *)PyArray_DATA(arr);

            for (iterIndex=0; iterIndex<arrSize; iterIndex++) {
                /* Get the offset <ii> for the current array value. Also
                   set valueIsIgnore flag. */
                switch(arrType) {
                case NPY_INT8:   CASE_VALUEINDICES_GET_VALUEOFFSET(npy_int8); break;
                case NPY_UINT8:  CASE_VALUEINDICES_GET_VALUEOFFSET(npy_uint8); break;
                case NPY_INT16:  CASE_VALUEINDICES_GET_VALUEOFFSET(npy_int16); break;
                case NPY_UINT16: CASE_VALUEINDICES_GET_VALUEOFFSET(npy_uint16); break;
                case NPY_INT32:  CASE_VALUEINDICES_GET_VALUEOFFSET(npy_int32); break;
                case NPY_UINT32: CASE_VALUEINDICES_GET_VALUEOFFSET(npy_uint32); break;
                case NPY_INT64:  CASE_VALUEINDICES_GET_VALUEOFFSET(npy_int64); break;
                case NPY_UINT64: CASE_VALUEINDICES_GET_VALUEOFFSET(npy_uint64); break;
                default:
                    switch(arrType) {
                    case NPY_INT:  CASE_VALUEINDICES_GET_VALUEOFFSET(npy_int); break;
                    case NPY_UINT: CASE_VALUEINDICES_GET_VALUEOFFSET(npy_uint); break;
                    }
                }

                if (ignoreValIsNone || (!valueIsIgnore)) {
                    ndxTuple = ndxPtr[ii];
                    for (j=0; j<ndim; j++) {
                        ndxArr = (PyArrayObject *)PyTuple_GetItem(ndxTuple, j);
                        *(npy_intp *)PyArray_GETPTR1(ndxArr, (npy_intp)valCtr[ii]) =
                            ndiIter.coordinates[j];
                    }
                    valCtr[ii] += 1;
                }

                NI_ITERATOR_NEXT(ndiIter, arrData);
            }
        }
    }

    /* Clean up everything */
    if (hist != NULL) free(hist);
    if (valCtr != NULL) free(valCtr);
    if (ndxPtr != NULL) free(ndxPtr);
    Py_DECREF(minMaxArr);

    if (PyErr_Occurred()) {
        Py_DECREF(valDict);
        return NULL;
    } else
        return valDict;
}

static PyObject *Py_WatershedIFT(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *markers = NULL;
    PyArrayObject *strct = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&O&", NI_ObjectToInputArray, &input,
                    NI_ObjectToInputArray, &markers, NI_ObjectToInputArray,
                    &strct, NI_ObjectToOutputArray, &output))
        goto exit;

    NI_WatershedIFT(input, markers, strct, output);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(markers);
    Py_XDECREF(strct);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_DistanceTransformBruteForce(PyObject *obj,
                                                                                                PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *features = NULL;
    PyArrayObject *sampling = NULL;
    int metric;

    if (!PyArg_ParseTuple(args, "O&iO&O&O&",
                          NI_ObjectToInputArray, &input,
                          &metric,
                          NI_ObjectToOptionalInputArray, &sampling,
                                                NI_ObjectToOptionalOutputArray, &output,
                                                NI_ObjectToOptionalOutputArray, &features))
        goto exit;

    NI_DistanceTransformBruteForce(input, metric, sampling, output, features);
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(sampling);
    Py_XDECREF(output);
    Py_XDECREF(features);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_DistanceTransformOnePass(PyObject *obj, PyObject *args)
{
    PyArrayObject *strct = NULL, *distances = NULL, *features = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&",
                          NI_ObjectToInputArray, &strct,
                          NI_ObjectToInputOutputArray, &distances,
                          NI_ObjectToOptionalOutputArray, &features))
        goto exit;

    NI_DistanceTransformOnePass(strct, distances, features);

exit:
    Py_XDECREF(strct);
    Py_XDECREF(distances);
    Py_XDECREF(features);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_EuclideanFeatureTransform(PyObject *obj,
                                                                                            PyObject *args)
{
    PyArrayObject *input = NULL, *features = NULL, *sampling = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&",
                          NI_ObjectToInputArray, &input,
                                                NI_ObjectToOptionalInputArray, &sampling,
                                                NI_ObjectToOutputArray, &features))
        goto exit;

    NI_EuclideanFeatureTransform(input, sampling, features);

exit:
    Py_XDECREF(input);
    Py_XDECREF(sampling);
    Py_XDECREF(features);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static void _FreeCoordinateList(PyObject *obj)
{
    NI_FreeCoordinateList((NI_CoordinateList*)PyCapsule_GetPointer(obj, NULL));
}

static PyObject *Py_BinaryErosion(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *strct = NULL;
    PyArrayObject *mask = NULL;
    PyObject *cobj = NULL;
    int border_value, invert, center_is_true;
    int changed = 0, return_coordinates;
    NI_CoordinateList *coordinate_list = NULL;
    PyArray_Dims origin = {NULL, 0};

    if (!PyArg_ParseTuple(args, "O&O&O&O&iO&iii",
                          NI_ObjectToInputArray, &input,
                          NI_ObjectToInputArray, &strct,
                          NI_ObjectToOptionalInputArray, &mask,
                          NI_ObjectToOutputArray, &output,
                          &border_value,
                          PyArray_IntpConverter, &origin,
                          &invert, &center_is_true, &return_coordinates)) {
        goto exit;
    }
    if (!_validate_origin(input, origin)) {
        goto exit;
    }
    if (!NI_BinaryErosion(input, strct, mask, output, border_value,
                          origin.ptr, invert, center_is_true, &changed,
                          return_coordinates ? &coordinate_list : NULL)) {
        goto exit;
    }
    if (return_coordinates) {
        cobj = PyCapsule_New(coordinate_list, NULL, _FreeCoordinateList);
    }
    PyArray_ResolveWritebackIfCopy(output);

exit:
    Py_XDECREF(input);
    Py_XDECREF(strct);
    Py_XDECREF(mask);
    Py_XDECREF(output);
    PyDimMem_FREE(origin.ptr);
    if (PyErr_Occurred()) {
        Py_XDECREF(cobj);
        return NULL;
    } else {
        if (return_coordinates) {
            return Py_BuildValue("iN", changed, cobj);
        } else {
            return Py_BuildValue("i", changed);
        }
    }
}

static PyObject *Py_BinaryErosion2(PyObject *obj, PyObject *args)
{
    PyArrayObject *array = NULL, *strct = NULL, *mask = NULL;
    PyObject *cobj = NULL;
    int invert, niter;
    PyArray_Dims origin = {NULL, 0};

    if (!PyArg_ParseTuple(args, "O&O&O&iO&iO",
                          NI_ObjectToInputOutputArray, &array,
                          NI_ObjectToInputArray, &strct,
                          NI_ObjectToOptionalInputArray,
                          &mask, &niter,
                          PyArray_IntpConverter, &origin,
                          &invert, &cobj)) {
        goto exit;
    }
    if (!_validate_origin(array, origin)) {
        goto exit;
    }
    if (PyCapsule_CheckExact(cobj)) {
        NI_CoordinateList *cobj_data = PyCapsule_GetPointer(cobj, NULL);
        if (!NI_BinaryErosion2(array, strct, mask, niter, origin.ptr, invert,
                               &cobj_data)) {
            goto exit;
        }
    }
    else {
        PyErr_SetString(PyExc_RuntimeError, "cannot convert CObject");
    }
exit:
    Py_XDECREF(array);
    Py_XDECREF(strct);
    Py_XDECREF(mask);
    PyDimMem_FREE(origin.ptr);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyMethodDef methods[] = {
    {"correlate1d",           (PyCFunction)Py_Correlate1D,
     METH_VARARGS, NULL},
    {"correlate",             (PyCFunction)Py_Correlate,
     METH_VARARGS, NULL},
    {"uniform_filter1d",      (PyCFunction)Py_UniformFilter1D,
     METH_VARARGS, NULL},
    {"min_or_max_filter1d",   (PyCFunction)Py_MinOrMaxFilter1D,
        METH_VARARGS, NULL},
    {"min_or_max_filter",     (PyCFunction)Py_MinOrMaxFilter,
        METH_VARARGS, NULL},
    {"rank_filter",           (PyCFunction)Py_RankFilter,
     METH_VARARGS, NULL},
    {"generic_filter",        (PyCFunction)Py_GenericFilter,
     METH_VARARGS, NULL},
    {"generic_filter1d",      (PyCFunction)Py_GenericFilter1D,
     METH_VARARGS, NULL},
    {"fourier_filter",        (PyCFunction)Py_FourierFilter,
     METH_VARARGS, NULL},
    {"fourier_shift",         (PyCFunction)Py_FourierShift,
     METH_VARARGS, NULL},
    {"spline_filter1d",       (PyCFunction)Py_SplineFilter1D,
     METH_VARARGS, NULL},
    {"geometric_transform",   (PyCFunction)Py_GeometricTransform,
        METH_VARARGS, NULL},
    {"zoom_shift",            (PyCFunction)Py_ZoomShift,
     METH_VARARGS, NULL},
    {"find_objects",          (PyCFunction)Py_FindObjects,
     METH_VARARGS, NULL},
    {"value_indices",         (PyCFunction)NI_ValueIndices,
     METH_VARARGS, NULL},
    {"watershed_ift",         (PyCFunction)Py_WatershedIFT,
     METH_VARARGS, NULL},
    {"distance_transform_bf", (PyCFunction)Py_DistanceTransformBruteForce,
     METH_VARARGS, NULL},
    {"distance_transform_op", (PyCFunction)Py_DistanceTransformOnePass,
     METH_VARARGS, NULL},
    {"euclidean_feature_transform",
     (PyCFunction)Py_EuclideanFeatureTransform,
     METH_VARARGS, NULL},
    {"binary_erosion",        (PyCFunction)Py_BinaryErosion,
     METH_VARARGS, NULL},
    {"binary_erosion2",       (PyCFunction)Py_BinaryErosion2,
     METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_nd_image",
    NULL,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__nd_image(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
