#include <Python.h>
#include <numpy/npy_common.h>


static void
_destructor(PyObject *obj)
{
    void *callback_data = PyCapsule_GetContext(obj);
    PyMem_Free(callback_data);
}


static int
_filter1d(double *input_line, npy_intp input_length, double *output_line,
	  npy_intp output_length, void *callback_data)
{
    npy_intp i, j;
    npy_intp filter_size = *(npy_intp *)callback_data;

    for (i = 0; i < output_length; i++) {
	output_line[i] = 0;
	for (j = 0; j < filter_size; j++) {
	    output_line[i] += input_line[i+j];
	}
	output_line[i] /= filter_size;
    }
    return 1;
}


static PyObject *
py_filter1d(PyObject *obj, PyObject *args)
{
    npy_intp *callback_data = NULL;
    PyObject *capsule = NULL;

    callback_data = PyMem_Malloc(sizeof(npy_intp));
    if (!callback_data) {
	PyErr_NoMemory();
	goto error;
    }
    if (!PyArg_ParseTuple(args, "n", callback_data)) goto error;

    capsule = PyCapsule_New(_filter1d, NULL, _destructor);
    if (!capsule) goto error;
    if (PyCapsule_SetContext(capsule, callback_data) != 0) {
	Py_DECREF(capsule);
	goto error;
    }
    return capsule;
 error:
    PyMem_Free(callback_data);
    return NULL;
}


static int
_filter2d(double *buffer, npy_intp filter_size, double *res,
	  void *callback_data)
{
    npy_intp i;
    double *weights = (double *)callback_data;

    *res = 0;
    for (i = 0; i < filter_size; i++) {
	*res += weights[i]*buffer[i];
    }
    return 1;
}


static PyObject *
py_filter2d(PyObject *obj, PyObject *args)
{
    Py_ssize_t i, size;
    double *callback_data = NULL;
    PyObject *seq = NULL, *item = NULL, *capsule = NULL;

    if (!PyArg_ParseTuple(args, "O", &seq)) goto error;

    size = PySequence_Length(seq);
    if (size == -1) goto error;
    callback_data = PyMem_Malloc(size*sizeof(double));
    if (!callback_data) {
	PyErr_NoMemory();
	goto error;
    }

    for (i = 0; i < size; i++) {
	item = PySequence_GetItem(seq, i);
	if (!item) {
	    PyErr_SetString(PyExc_IndexError, "failed to get item");
	    goto error;
	}
	callback_data[i] = PyFloat_AsDouble(item);
        Py_DECREF(item);
        item = NULL;
	if (PyErr_Occurred()) goto error;
    }

    capsule = PyCapsule_New(_filter2d, NULL, _destructor);
    if (!capsule) goto error;
    if (PyCapsule_SetContext(capsule, callback_data) != 0) {
	Py_DECREF(capsule);
	goto error;
    }
    return capsule;
 error:
    PyMem_Free(callback_data);
    return NULL;
}


static int
_transform(npy_intp *output_coordinates, double *input_coordinates,
	   npy_intp output_rank, npy_intp input_rank, void *callback_data)
{
    npy_intp i;
    double shift = *(double *)callback_data;

    for (i = 0; i < input_rank; i++) {
	input_coordinates[i] = output_coordinates[i] - shift;
    }
    return 1;
}


static PyObject *
py_transform(PyObject *obj, PyObject *args)
{
    double *callback_data = PyMem_Malloc(sizeof(double));
    PyObject *capsule = NULL;

    if (!callback_data) {
	PyErr_NoMemory();
	goto error;
    }
    if (!PyArg_ParseTuple(args, "d", callback_data)) goto error;

    capsule = PyCapsule_New(_transform, NULL, _destructor);
    if (!capsule) goto error;
    if (PyCapsule_SetContext(capsule, callback_data) != 0) {
	Py_DECREF(capsule);
	goto error;
    }
    return capsule;
 error:
    PyMem_Free(callback_data);
    return NULL;
}


static PyMethodDef _CTestMethods[] = {
    {"transform", (PyCFunction)py_transform, METH_VARARGS, ""},
    {"filter1d", (PyCFunction)py_filter1d, METH_VARARGS, ""},
    {"filter2d", (PyCFunction)py_filter2d, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};


/* Initialize the module */
static struct PyModuleDef _ctest = {
    PyModuleDef_HEAD_INIT,
    "_ctest",
    NULL,
    -1,
    _CTestMethods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC
PyInit__ctest(void)
{
    PyObject *module;

    module = PyModule_Create(&_ctest);
    if (module == NULL) {
        return module;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}
