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

#define ND_IMPORT_ARRAY
#include "nd_image.h"
#undef ND_IMPORT_ARRAY
#include "ni_support.h"
#include "ni_filters.h"
#include "ni_fourier.h"
#include "ni_morphology.h"
#include "ni_interpolation.h"
#include "ni_measure.h"

typedef struct {
    PyObject *function;
    PyObject *extra_arguments;
    PyObject *extra_keywords;
} NI_PythonCallbackData;

/* Convert an input array of any type, not necessarily contiguous */
static int
NI_ObjectToInputArray(PyObject *object, PyArrayObject **array)
{
    *array = NA_InputArray(object, tAny, NPY_ALIGNED|NPY_NOTSWAPPED);
    return *array ? 1 : 0;
}

/* Convert an input array of any type, not necessarily contiguous */
static int
NI_ObjectToOptionalInputArray(PyObject *object, PyArrayObject **array)
{
    if (object == Py_None) {
        *array = NULL;
        return 1;
    } else {
        *array = NA_InputArray(object, tAny, NPY_ALIGNED|NPY_NOTSWAPPED);
        return *array ? 1 : 0;
    }
}

/* Convert an output array of any type, not necessarily contiguous */
static int
NI_ObjectToOutputArray(PyObject *object, PyArrayObject **array)
{
    *array = NA_OutputArray(object, tAny, NPY_ALIGNED|NPY_NOTSWAPPED);
    return *array ? 1 : 0;
}

/* Convert an output array of any type, not necessarily contiguous */
static int
NI_ObjectToOptionalOutputArray(PyObject *object, PyArrayObject **array)
{
    if (object == Py_None) {
        *array = NULL;
        return 1;
    } else {
        *array = NA_OutputArray(object, tAny, NPY_ALIGNED|NPY_NOTSWAPPED);
        return *array ? 1 : 0;
    }
}

/* Convert an input/output array of any type, not necessarily contiguous */
static int
NI_ObjectToIoArray(PyObject *object, PyArrayObject **array)
{
    *array = NA_IoArray(object, tAny, NPY_ALIGNED|NPY_NOTSWAPPED);
    return *array ? 1 : 0;
}

/* Convert an Long sequence */
static maybelong
NI_ObjectToLongSequenceAndLength(PyObject *object, maybelong **sequence)
{
    long *pa, ii;
    PyArrayObject *array = NA_InputArray(object, PyArray_LONG, NPY_CARRAY);
    maybelong length = PyArray_SIZE(array);

    *sequence = (maybelong*)malloc(length * sizeof(maybelong));
    if (!*sequence) {
        PyErr_NoMemory();
        Py_XDECREF(array);
        return -1;
    }
    pa = (long*)PyArray_DATA(array);
    for(ii = 0; ii < length; ii++)
        (*sequence)[ii] = pa[ii];
    Py_XDECREF(array);
    return length;
}

static int
NI_ObjectToLongSequence(PyObject *object, maybelong **sequence)
{
    return NI_ObjectToLongSequenceAndLength(object, sequence) >= 0;
}

/*********************************************************************/
/* wrapper functions: */
/*********************************************************************/

static PyObject *Py_Correlate1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *weights = NULL;
    int axis, mode;
    long origin;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&iO&idl", NI_ObjectToInputArray, &input,
                                    NI_ObjectToInputArray, &weights, &axis,
                                    NI_ObjectToOutputArray, &output, &mode, &cval, &origin))
        goto exit;
    if (!NI_Correlate1D(input, weights, axis, output,
                                            (NI_ExtendMode)mode, cval, origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(weights);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_Correlate(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *weights = NULL;
    maybelong *origin = NULL;
    int mode;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&idO&", NI_ObjectToInputArray, &input,
                    NI_ObjectToInputArray, &weights, NI_ObjectToOutputArray, &output,
                    &mode, &cval, NI_ObjectToLongSequence, &origin))
        goto exit;
    if (!NI_Correlate(input, weights, output, (NI_ExtendMode)mode, cval,
                                        origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(weights);
    Py_XDECREF(output);
    if (origin)
        free(origin);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_UniformFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, mode;
    long filter_size, origin;
    double cval;

    if (!PyArg_ParseTuple(args, "O&liO&idl", NI_ObjectToInputArray, &input,
                                        &filter_size, &axis, NI_ObjectToOutputArray, &output,
                                        &mode, &cval, &origin))
        goto exit;
    if (!NI_UniformFilter1D(input, filter_size, axis, output,
                                                                             (NI_ExtendMode)mode, cval, origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_MinOrMaxFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, mode, minimum;
    long filter_size, origin;
    double cval;

    if (!PyArg_ParseTuple(args, "O&liO&idli", NI_ObjectToInputArray, &input,
                                        &filter_size, &axis, NI_ObjectToOutputArray, &output,
                                        &mode, &cval, &origin, &minimum))
        goto exit;
    if (!NI_MinOrMaxFilter1D(input, filter_size, axis, output,
                                                            (NI_ExtendMode)mode, cval, origin, minimum))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_MinOrMaxFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *footprint = NULL;
    PyArrayObject *structure = NULL;
    maybelong *origin = NULL;
    int mode, minimum;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&O&idO&i", NI_ObjectToInputArray,
                                        &input, NI_ObjectToInputArray, &footprint,
                                        NI_ObjectToOptionalInputArray, &structure,
                                        NI_ObjectToOutputArray, &output, &mode, &cval,
                                        NI_ObjectToLongSequence, &origin, &minimum))
        goto exit;
    if (!NI_MinOrMaxFilter(input, footprint, structure, output,
                                                (NI_ExtendMode)mode, cval, origin, minimum))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(footprint);
    Py_XDECREF(structure);
    Py_XDECREF(output);
    if (origin)
        free(origin);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_RankFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *footprint = NULL;
    maybelong *origin = NULL;
    int mode, rank;
    double cval;

    if (!PyArg_ParseTuple(args, "O&iO&O&idO&", NI_ObjectToInputArray,
                                        &input, &rank, NI_ObjectToInputArray, &footprint,
                                        NI_ObjectToOutputArray, &output, &mode, &cval,
                                        NI_ObjectToLongSequence, &origin))
        goto exit;
    if (!NI_RankFilter(input, rank, footprint, output, (NI_ExtendMode)mode,
                                         cval, origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(footprint);
    Py_XDECREF(output);
    if (origin)
        free(origin);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_Filter1DFunc(double *iline, maybelong ilen,
                                                     double *oline, maybelong olen, void *data)
{
    PyArrayObject *py_ibuffer = NULL, *py_obuffer = NULL;
    PyObject *rv = NULL, *args = NULL, *tmp = NULL;
    maybelong ii;
    double *po = NULL;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)data;

    py_ibuffer = NA_NewArray(iline, PyArray_DOUBLE, 1, &ilen);
    py_obuffer = NA_NewArray(NULL, PyArray_DOUBLE, 1, &olen);
    if (!py_ibuffer || !py_obuffer)
        goto exit;
    tmp = Py_BuildValue("(OO)", py_ibuffer, py_obuffer);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rv = PyObject_Call(cbdata->function, args, cbdata->extra_keywords);
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
    void *func = Py_Filter1DFunc, *data = NULL;
    NI_PythonCallbackData cbdata;
    int axis, mode;
    long origin, filter_size;
    double cval;

    if (!PyArg_ParseTuple(args, "O&OliO&idlOO", NI_ObjectToInputArray,
                &input, &fnc, &filter_size, &axis, NI_ObjectToOutputArray,
                &output, &mode, &cval, &origin, &extra_arguments, &extra_keywords))
        goto exit;
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
    if (PyCObject_Check(fnc)) {
        func = PyCObject_AsVoidPtr(fnc);
        data = PyCObject_GetDesc(fnc);
    } else if (PyCallable_Check(fnc)) {
        cbdata.function = fnc;
        cbdata.extra_arguments = extra_arguments;
        cbdata.extra_keywords = extra_keywords;
        data = (void*)&cbdata;
    } else {
        PyErr_SetString(PyExc_RuntimeError,
                                        "function parameter is not callable");
        goto exit;
    }
    if (!NI_GenericFilter1D(input, func, data, filter_size, axis, output,
                                                                                (NI_ExtendMode)mode, cval, origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_FilterFunc(double *buffer, maybelong filter_size,
                                                 double *output, void *data)
{
    PyArrayObject *py_buffer = NULL;
    PyObject *rv = NULL, *args = NULL, *tmp = NULL;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)data;

    py_buffer = NA_NewArray(buffer, PyArray_DOUBLE, 1, &filter_size);
    if (!py_buffer)
        goto exit;
    tmp = Py_BuildValue("(O)", py_buffer);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rv = PyObject_Call(cbdata->function, args, cbdata->extra_keywords);
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
    void *func = Py_FilterFunc, *data = NULL;
    NI_PythonCallbackData cbdata;
    int mode;
    maybelong *origin = NULL;
    double cval;

    if (!PyArg_ParseTuple(args, "O&OO&O&idO&OO", NI_ObjectToInputArray,
                                                &input, &fnc, NI_ObjectToInputArray, &footprint,
                                                NI_ObjectToOutputArray, &output, &mode, &cval,
                                                NI_ObjectToLongSequence, &origin,
                                                &extra_arguments, &extra_keywords))
        goto exit;
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
    if (PyCObject_Check(fnc)) {
        func = PyCObject_AsVoidPtr(fnc);
        data = PyCObject_GetDesc(fnc);
    } else if (PyCallable_Check(fnc)) {
        cbdata.function = fnc;
        cbdata.extra_arguments = extra_arguments;
        cbdata.extra_keywords = extra_keywords;
        data = (void*)&cbdata;
    } else {
        PyErr_SetString(PyExc_RuntimeError,
                                        "function parameter is not callable");
        goto exit;
    }
    if (!NI_GenericFilter(input, func, data, footprint, output,
                                                (NI_ExtendMode)mode, cval, origin))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    Py_XDECREF(footprint);
    if (origin)
        free(origin);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_FourierFilter(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *parameters = NULL;
    int axis, filter_type;
    long n;

    if (!PyArg_ParseTuple(args, "O&O&liO&i", NI_ObjectToInputArray, &input,
                                        NI_ObjectToInputArray, &parameters, &n, &axis,
                                        NI_ObjectToOutputArray, &output, &filter_type))
        goto exit;

    if (!NI_FourierFilter(input, parameters, n, axis, output, filter_type))
        goto exit;

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
    long n;

    if (!PyArg_ParseTuple(args, "O&O&liO&", NI_ObjectToInputArray, &input,
                                        NI_ObjectToInputArray, &shifts, &n, &axis,
                                        NI_ObjectToOutputArray, &output))
        goto exit;

    if (!NI_FourierShift(input, shifts, n, axis, output))
        goto exit;

exit:
    Py_XDECREF(input);
    Py_XDECREF(shifts);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_SplineFilter1D(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL;
    int axis, order;

    if (!PyArg_ParseTuple(args, "O&iiO&", NI_ObjectToInputArray, &input,
                    &order, &axis, NI_ObjectToOutputArray, &output))
        goto exit;

    if (!NI_SplineFilter1D(input, order, axis, output))
        goto exit;

exit:
    Py_XDECREF(input);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int Py_Map(maybelong *ocoor, double* icoor, int orank, int irank,
                                    void *data)
{
    PyObject *coors = NULL, *rets = NULL, *args = NULL, *tmp = NULL;
    maybelong ii;
    NI_PythonCallbackData *cbdata = (NI_PythonCallbackData*)data;

    coors = PyTuple_New(orank);
    if (!coors)
        goto exit;
    for(ii = 0; ii < orank; ii++) {
        PyTuple_SetItem(coors, ii, PyInt_FromLong(ocoor[ii]));
        if (PyErr_Occurred())
            goto exit;
    }
    tmp = Py_BuildValue("(O)", coors);
    if (!tmp)
        goto exit;
    args = PySequence_Concat(tmp, cbdata->extra_arguments);
    if (!args)
        goto exit;
    rets = PyObject_Call(cbdata->function, args, cbdata->extra_keywords);
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
    int mode, order;
    double cval;
    void *func = NULL, *data = NULL;
    NI_PythonCallbackData cbdata;

    if (!PyArg_ParseTuple(args, "O&OO&O&O&O&iidOO", NI_ObjectToInputArray,
                                                &input, &fnc, NI_ObjectToOptionalInputArray,
                                                &coordinates, NI_ObjectToOptionalInputArray,
                                                &matrix, NI_ObjectToOptionalInputArray, &shift,
                                                NI_ObjectToOutputArray, &output, &order, &mode,
                                                &cval, &extra_arguments, &extra_keywords))
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
        if (PyCObject_Check(fnc)) {
            func = PyCObject_AsVoidPtr(fnc);
            data = PyCObject_GetDesc(fnc);
        } else if (PyCallable_Check(fnc)) {
            func = Py_Map;
            cbdata.function = fnc;
            cbdata.extra_arguments = extra_arguments;
            cbdata.extra_keywords = extra_keywords;
            data = (void*)&cbdata;
        } else {
            PyErr_SetString(PyExc_RuntimeError,
                                            "function parameter is not callable");
            goto exit;
        }
    }

    if (!NI_GeometricTransform(input, func, data, matrix, shift, coordinates,
                                                    output, order, (NI_ExtendMode)mode, cval))
        goto exit;

exit:
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
    int mode, order;
    double cval;

    if (!PyArg_ParseTuple(args, "O&O&O&O&iid", NI_ObjectToInputArray,
                    &input, NI_ObjectToOptionalInputArray, &zoom,
                    NI_ObjectToOptionalInputArray, &shift, NI_ObjectToOutputArray,
                    &output, &order, &mode, &cval))
        goto exit;

    if (!NI_ZoomShift(input, zoom, shift, output, order, (NI_ExtendMode)mode,
                                        cval))
        goto exit;

exit:
    Py_XDECREF(input);
    Py_XDECREF(shift);
    Py_XDECREF(zoom);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static PyObject *Py_Label(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *strct = NULL;
    maybelong max_label;

    if (!PyArg_ParseTuple(args, "O&O&O&", NI_ObjectToInputArray, &input,
                    NI_ObjectToInputArray, &strct, NI_ObjectToOutputArray, &output))
        goto exit;

    if (!NI_Label(input, strct, &max_label, output))
        goto exit;

exit:
    Py_XDECREF(input);
    Py_XDECREF(strct);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("l", (long)max_label);
}

static PyObject *Py_FindObjects(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL;
    PyObject *result = NULL, *tuple = NULL, *start = NULL, *end = NULL;
    PyObject *slc = NULL;
    int jj;
    long max_label;
    maybelong ii, *regions = NULL;

    if (!PyArg_ParseTuple(args, "O&l", NI_ObjectToInputArray, &input,
                                                &max_label))
        goto exit;

    if (max_label < 0)
        max_label = 0;
    if (max_label > 0) {
        if (input->nd > 0) {
            regions = (maybelong*)malloc(2 * max_label * input->nd *
                                                                                                        sizeof(maybelong));
        } else {
            regions = (maybelong*)malloc(max_label * sizeof(maybelong));
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
        maybelong idx = input->nd > 0 ? 2 * input->nd * ii : ii;
        if (regions[idx] >= 0) {
            PyObject *tuple = PyTuple_New(input->nd);
            if (!tuple) {
                PyErr_NoMemory();
                goto exit;
            }
            for(jj = 0; jj < input->nd; jj++) {
                start = PyInt_FromLong(regions[idx + jj]);
                end = PyInt_FromLong(regions[idx + jj + input->nd]);
                if (!start || !end) {
                    PyErr_NoMemory();
                    goto exit;
                }
                slc = PySlice_New(start, end, NULL);
                if (!slc) {
                    PyErr_NoMemory();
                    goto exit;
                }
                Py_XDECREF(start);
                Py_XDECREF(end);
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
    if (regions)
        free(regions);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        return NULL;
    } else {
        return result;
    }
}

static PyObject *Py_WatershedIFT(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *markers = NULL;
    PyArrayObject *strct = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&O&", NI_ObjectToInputArray, &input,
                    NI_ObjectToInputArray, &markers, NI_ObjectToInputArray,
                    &strct, NI_ObjectToOutputArray, &output))
        goto exit;

    if (!NI_WatershedIFT(input, markers, strct, output))
        goto exit;

exit:
    Py_XDECREF(input);
    Py_XDECREF(markers);
    Py_XDECREF(strct);
    Py_XDECREF(output);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static int _NI_GetIndices(PyObject* indices_object,
                                                    maybelong** result_indices, maybelong* min_label,
                                                    maybelong* max_label, maybelong* n_results)
{
    maybelong *indices = NULL, n_indices, ii;

    if (indices_object == Py_None) {
        *min_label = -1;
        *n_results = 1;
    } else {
        n_indices = NI_ObjectToLongSequenceAndLength(indices_object, &indices);
        if (n_indices < 0)
            goto exit;
        if (n_indices < 1) {
            PyErr_SetString(PyExc_RuntimeError, "no correct indices provided");
            goto exit;
        } else {
            *min_label = *max_label = indices[0];
            if (*min_label < 0) {
                PyErr_SetString(PyExc_RuntimeError,
                                                "negative indices not allowed");
                goto exit;
            }
            for(ii = 1; ii < n_indices; ii++) {
                if (indices[ii] < 0) {
                    PyErr_SetString(PyExc_RuntimeError,
                                                    "negative indices not allowed");
                    goto exit;
                }
                if (indices[ii] < *min_label)
                    *min_label = indices[ii];
                if (indices[ii] > *max_label)
                    *max_label = indices[ii];
            }
            *result_indices = (maybelong*)malloc((*max_label - *min_label + 1) *
                                                                                     sizeof(maybelong));
            if (!*result_indices) {
                PyErr_NoMemory();
                goto exit;
            }
            for(ii = 0; ii < *max_label - *min_label + 1; ii++)
                (*result_indices)[ii] = -1;
            *n_results = 0;
            for(ii = 0; ii < n_indices; ii++) {
                if ((*result_indices)[indices[ii] - *min_label] >= 0) {
                    PyErr_SetString(PyExc_RuntimeError, "duplicate index");
                    goto exit;
                }
                (*result_indices)[indices[ii] - *min_label] = ii;
                ++(*n_results);
            }
        }
    }
 exit:
    if (indices)
        free(indices);
    return PyErr_Occurred() == NULL;
}


PyObject* _NI_BuildMeasurementResultArrayObject(maybelong n_results,
                                                                                                PyArrayObject** values)
{
    PyObject *result = NULL;
    if (n_results > 1) {
        result = PyList_New(n_results);
        if (result) {
            maybelong ii;
            for(ii = 0; ii < n_results; ii++) {
                PyList_SET_ITEM(result, ii, (PyObject*)values[ii]);
                Py_XINCREF(values[ii]);
            }
        }
    } else {
        result = (PyObject*)values[0];
        Py_XINCREF(values[0]);
    }
    return result;
}


PyObject* _NI_BuildMeasurementResultDouble(maybelong n_results,
                                                                                     double* values)
{
    PyObject *result = NULL;
    if (n_results > 1) {
        result = PyList_New(n_results);
        if (result) {
            int ii;
            for(ii = 0; ii < n_results; ii++) {
                PyObject* val = PyFloat_FromDouble(values[ii]);
                if (!val) {
                    Py_XDECREF(result);
                    return NULL;
                }
                PyList_SET_ITEM(result, ii, val);
            }
        }
    } else {
        result = Py_BuildValue("d", values[0]);
    }
    return result;
}


PyObject* _NI_BuildMeasurementResultDoubleTuple(maybelong n_results,
                                                                                        int tuple_size, double* values)
{
    PyObject *result = NULL;
    maybelong ii;
    int jj;

    if (n_results > 1) {
        result = PyList_New(n_results);
        if (result) {
            for(ii = 0; ii < n_results; ii++) {
                PyObject* val = PyTuple_New(tuple_size);
                if (!val) {
                    Py_XDECREF(result);
                    return NULL;
                }
                for(jj = 0; jj < tuple_size; jj++) {
                    maybelong idx = jj + ii * tuple_size;
                    PyTuple_SetItem(val, jj, PyFloat_FromDouble(values[idx]));
                    if (PyErr_Occurred()) {
                        Py_XDECREF(result);
                        return NULL;
                    }
                }
                PyList_SET_ITEM(result, ii, val);
            }
        }
    } else {
        result = PyTuple_New(tuple_size);
        if (result) {
            for(ii = 0; ii < tuple_size; ii++) {
                PyTuple_SetItem(result, ii, PyFloat_FromDouble(values[ii]));
                if (PyErr_Occurred()) {
                    Py_XDECREF(result);
                    return NULL;
                }
            }
        }
    }
    return result;
}


PyObject* _NI_BuildMeasurementResultInt(maybelong n_results,
                                                                                maybelong* values)
{
    PyObject *result = NULL;
    if (n_results > 1) {
        result = PyList_New(n_results);
        if (result) {
            maybelong ii;
            for(ii = 0; ii < n_results; ii++) {
                PyObject* val = PyInt_FromLong(values[ii]);
                if (!val) {
                    Py_XDECREF(result);
                    return NULL;
                }
                PyList_SET_ITEM(result, ii, val);
            }
        }
    } else {
        result = Py_BuildValue("l", values[0]);
    }
    return result;
}


static PyObject *Py_Statistics(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *labels = NULL;
    PyObject *indices_object, *result = NULL;
    PyObject *res1 = NULL, *res2 = NULL, *res3 = NULL, *res4 = NULL;
    double *dresult1 = NULL, *dresult2 = NULL;
    maybelong *lresult1 = NULL, *lresult2 = NULL;
    maybelong min_label, max_label, *result_indices = NULL, n_results, ii;
    int type;

    if (!PyArg_ParseTuple(args, "O&O&Oi", NI_ObjectToInputArray, &input,
                    NI_ObjectToOptionalInputArray, &labels, &indices_object, &type))
        goto exit;

    if (!_NI_GetIndices(indices_object, &result_indices, &min_label,
                                            &max_label, &n_results))
        goto exit;

    if (type >= 0 && type <= 7) {
        dresult1 = (double*)malloc(n_results * sizeof(double));
        if (!dresult1) {
            PyErr_NoMemory();
            goto exit;
        }
    }
    if (type == 2 || type == 7) {
        dresult2 = (double*)malloc(n_results * sizeof(double));
        if (!dresult2) {
            PyErr_NoMemory();
            goto exit;
        }
    }
    if (type == 1 || type == 2 || (type >= 5 && type <= 7)) {
        lresult1 = (maybelong*)malloc(n_results * sizeof(maybelong));
        if (!lresult1) {
            PyErr_NoMemory();
            goto exit;
        }
    }
    if (type == 7) {
        lresult2 = (maybelong*)malloc(n_results * sizeof(maybelong));
        if (!lresult2) {
            PyErr_NoMemory();
            goto exit;
        }
    }
    switch(type) {
    case 0:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                                    n_results, dresult1, NULL, NULL, NULL, NULL, NULL, NULL))
            goto exit;
        result = _NI_BuildMeasurementResultDouble(n_results, dresult1);
        break;
    case 1:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                            n_results, dresult1, lresult1, NULL, NULL, NULL, NULL, NULL))
            goto exit;
        for(ii = 0; ii < n_results; ii++)
            dresult1[ii] = lresult1[ii] > 0 ? dresult1[ii] / lresult1[ii] : 0.0;

        result = _NI_BuildMeasurementResultDouble(n_results, dresult1);
        break;
    case 2:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                    n_results, dresult1, lresult1, dresult2, NULL, NULL, NULL, NULL))
            goto exit;
        result = _NI_BuildMeasurementResultDouble(n_results, dresult2);
        break;
    case 3:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                                    n_results, NULL, NULL, NULL, dresult1, NULL, NULL, NULL))
            goto exit;
        result = _NI_BuildMeasurementResultDouble(n_results, dresult1);
        break;
    case 4:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                                    n_results, NULL, NULL, NULL, NULL, dresult1, NULL, NULL))
            goto exit;
        result = _NI_BuildMeasurementResultDouble(n_results, dresult1);
        break;
    case 5:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                            n_results, NULL, NULL, NULL, dresult1, NULL, lresult1, NULL))
            goto exit;
        result = _NI_BuildMeasurementResultInt(n_results, lresult1);
        break;
    case 6:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                            n_results, NULL, NULL, NULL, NULL, dresult1, NULL, lresult1))
            goto exit;
        result = _NI_BuildMeasurementResultInt(n_results, lresult1);
        break;
    case 7:
        if (!NI_Statistics(input, labels, min_label, max_label, result_indices,
                                             n_results, NULL, NULL, NULL, dresult1, dresult2,
                                             lresult1, lresult2))
            goto exit;
        res1 = _NI_BuildMeasurementResultDouble(n_results, dresult1);
        res2 = _NI_BuildMeasurementResultDouble(n_results, dresult2);
        res3 = _NI_BuildMeasurementResultInt(n_results, lresult1);
        res4 = _NI_BuildMeasurementResultInt(n_results, lresult2);
        if (!res1 || !res2 || !res3 || !res4)
            goto exit;
        result = Py_BuildValue("OOOO", res1, res2, res3, res4);
        break;
    default:
        PyErr_SetString(PyExc_RuntimeError, "operation not supported");
        goto exit;
    }

 exit:
    Py_XDECREF(input);
    Py_XDECREF(labels);
    if (result_indices)
        free(result_indices);
    if (dresult1)
        free(dresult1);
    if (dresult2)
        free(dresult2);
    if (lresult1)
        free(lresult1);
    if (lresult2)
        free(lresult2);
    return result;
}


static PyObject *Py_CenterOfMass(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *labels = NULL;
    PyObject *indices_object, *result = NULL;
    double *center_of_mass = NULL;
    maybelong min_label, max_label, *result_indices = NULL, n_results;

    if (!PyArg_ParseTuple(args, "O&O&O", NI_ObjectToInputArray, &input,
                                    NI_ObjectToOptionalInputArray, &labels, &indices_object))
        goto exit;

    if (!_NI_GetIndices(indices_object, &result_indices, &min_label,
                                            &max_label, &n_results))
        goto exit;

    center_of_mass = (double*)malloc(input->nd * n_results *
                                                                     sizeof(double));
    if (!center_of_mass) {
        PyErr_NoMemory();
        goto exit;
    }

    if (!NI_CenterOfMass(input, labels, min_label, max_label,
                                             result_indices, n_results, center_of_mass))
        goto exit;

    result = _NI_BuildMeasurementResultDoubleTuple(n_results, input->nd,
                                                                                                 center_of_mass);

 exit:
    Py_XDECREF(input);
    Py_XDECREF(labels);
    if (result_indices)
        free(result_indices);
    if (center_of_mass)
        free(center_of_mass);
    return result;
}

static PyObject *Py_Histogram(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *labels = NULL, **histograms = NULL;
    PyObject *indices_object, *result = NULL;
    maybelong min_label, max_label, *result_indices = NULL, n_results;
    maybelong jj, nbins;
    long nbins_in;
    double min, max;

    if (!PyArg_ParseTuple(args, "O&ddlO&O", NI_ObjectToInputArray, &input,
                                                &min, &max, &nbins_in, NI_ObjectToOptionalInputArray,
                                                &labels, &indices_object))
        goto exit;
    nbins = nbins_in;

    if (!_NI_GetIndices(indices_object, &result_indices, &min_label,
                                            &max_label, &n_results))
        goto exit;

    /* Set all pointers to NULL, so that freeing the memory */
    /* doesn't cause problems. */
    histograms = (PyArrayObject**)calloc(input->nd * n_results,
                                                                             sizeof(PyArrayObject*));
    if (!histograms) {
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < n_results; jj++) {
        histograms[jj] = NA_NewArray(NULL, tInt32, 1, &nbins);
        if (!histograms[jj]) {
            PyErr_NoMemory();
            goto exit;
        }
    }

    if (!NI_Histogram(input, labels, min_label, max_label, result_indices,
                                        n_results, histograms, min, max, nbins))
        goto exit;

    result = _NI_BuildMeasurementResultArrayObject(n_results, histograms);

 exit:
    Py_XDECREF(input);
    Py_XDECREF(labels);
    if (result_indices)
        free(result_indices);
    if (histograms) {
        for(jj = 0; jj < n_results; jj++) {
            Py_XDECREF(histograms[jj]);
        }
        free(histograms);
    }
    return result;
}

static PyObject *Py_DistanceTransformBruteForce(PyObject *obj,
                                                                                                PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *features = NULL;
    PyArrayObject *sampling = NULL;
    int metric;

    if (!PyArg_ParseTuple(args, "O&iO&O&O&", NI_ObjectToInputArray, &input,
                                                &metric, NI_ObjectToOptionalInputArray, &sampling,
                                                NI_ObjectToOptionalOutputArray, &output,
                                                NI_ObjectToOptionalOutputArray, &features))
        goto exit;
    if (!NI_DistanceTransformBruteForce(input, metric, sampling,
                                                                            output, features))
        goto exit;
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

    if (!PyArg_ParseTuple(args, "O&O&O&", NI_ObjectToInputArray, &strct,
                                                NI_ObjectToIoArray, &distances,
                                                NI_ObjectToOptionalOutputArray, &features))
        goto exit;
    if (!NI_DistanceTransformOnePass(strct, distances, features))
        goto exit;
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

    if (!PyArg_ParseTuple(args, "O&O&O&", NI_ObjectToInputArray, &input,
                                                NI_ObjectToOptionalInputArray, &sampling,
                                                NI_ObjectToOutputArray, &features))
        goto exit;
    if (!NI_EuclideanFeatureTransform(input, sampling, features))
        goto exit;
exit:
    Py_XDECREF(input);
    Py_XDECREF(sampling);
    Py_XDECREF(features);
    return PyErr_Occurred() ? NULL : Py_BuildValue("");
}

static void _FreeCoordinateList(void* ptr)
{
    NI_FreeCoordinateList((NI_CoordinateList*)ptr);
}

static PyObject *Py_BinaryErosion(PyObject *obj, PyObject *args)
{
    PyArrayObject *input = NULL, *output = NULL, *strct = NULL;
    PyArrayObject *mask = NULL;
    PyObject *cobj = NULL;
    int border_value, invert, center_is_true;
    int changed = 0, return_coordinates;
    NI_CoordinateList *coordinate_list = NULL;
    maybelong *origins = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&O&iO&iii", NI_ObjectToInputArray,
                                                &input, NI_ObjectToInputArray, &strct,
                                                NI_ObjectToOptionalInputArray, &mask,
                                                NI_ObjectToOutputArray, &output, &border_value,
                                                NI_ObjectToLongSequence, &origins,  &invert,
                                                &center_is_true, &return_coordinates))
        goto exit;
    if (!NI_BinaryErosion(input, strct, mask, output, border_value,
                                                origins, invert, center_is_true, &changed,
                                                return_coordinates ? &coordinate_list : NULL))
        goto exit;
    if (return_coordinates) {
        cobj = PyCObject_FromVoidPtr(coordinate_list, _FreeCoordinateList);
    }
exit:
    Py_XDECREF(input);
    Py_XDECREF(strct);
    Py_XDECREF(mask);
    Py_XDECREF(output);
    if (origins)
        free(origins);
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
    maybelong *origins = NULL;

    if (!PyArg_ParseTuple(args, "O&O&O&iO&iO", NI_ObjectToIoArray, &array,
                    NI_ObjectToInputArray, &strct, NI_ObjectToOptionalInputArray,
                    &mask, &niter, NI_ObjectToLongSequence, &origins, &invert,
                    &cobj))
        goto exit;

    if (PyCObject_Check(cobj)) {
        NI_CoordinateList *cobj_data = PyCObject_AsVoidPtr(cobj);
        if (!NI_BinaryErosion2(array, strct, mask, niter, origins, invert,
                                                     &cobj_data))
            goto exit;
    } else {
        PyErr_SetString(PyExc_RuntimeError, "cannot convert CObject");
        goto exit;
    }
exit:
    Py_XDECREF(array);
    Py_XDECREF(strct);
    Py_XDECREF(mask);
    if (origins) free(origins);
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
    {"label",                 (PyCFunction)Py_Label,
     METH_VARARGS, NULL},
    {"find_objects",          (PyCFunction)Py_FindObjects,
     METH_VARARGS, NULL},
    {"watershed_ift",         (PyCFunction)Py_WatershedIFT,
     METH_VARARGS, NULL},
    {"statistics",            (PyCFunction)Py_Statistics,
     METH_VARARGS, NULL},
    {"center_of_mass",        (PyCFunction)Py_CenterOfMass,
     METH_VARARGS, NULL},
    {"histogram",             (PyCFunction)Py_Histogram,
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

PyMODINIT_FUNC init_nd_image(void)
{
    Py_InitModule("_nd_image", methods);
    import_array();
}
