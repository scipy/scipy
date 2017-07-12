#ifndef NI_CONVERTERS_H
#define NI_CONVERTERS_H

/*
 * Creates a new numpy array of the requested type and shape, and either
 * copies into it the contents of buffer, or sets it to all zeros if
 * buffer is NULL.
 */
static PyArrayObject *
NI_NewArray(void *buffer, enum NPY_TYPES type, int ndim, npy_intp *shape)
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
    int flags = NPY_ARRAY_BEHAVED_NS | NPY_ARRAY_UPDATEIFCOPY;
    /*
     * These would also be caught by the PyArray_CheckFromAny call, but
     * we check them explicitly here to provide a saner error message.
     */
    if (!PyArray_Check(object)) {
        PyErr_SetString(PyExc_TypeError, "output array must be an array");
        return 0;
    }
    if (!PyArray_ISWRITEABLE((PyArrayObject *)object)) {
        PyErr_SetString(PyExc_ValueError, "output array is read-only");
        return 0;
    }
    /*
     * If the output array is not aligned or is byteswapped, this call
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
    int flags = NPY_ARRAY_BEHAVED_NS | NPY_ARRAY_UPDATEIFCOPY;
    /*
     * These would also be caught by the PyArray_CheckFromAny call, but
     * we check them explicitly here to provide a saner error message.
     */
    if (!PyArray_Check(object)) {
        PyErr_SetString(PyExc_TypeError,
                        "input-output array must be an array");
        return 0;
    }
    if (!PyArray_ISWRITEABLE((PyArrayObject *)object)) {
        PyErr_SetString(PyExc_ValueError, "input-output array is read-only");
        return 0;
    }
    *array = (PyArrayObject *)PyArray_CheckFromAny(object, NULL, 0, 0, flags,
                                                   NULL);
    return *array != NULL;}

/* Checks that an origin value was received for each array dimension. */
static int
_validate_origin(PyArrayObject *array, PyArray_Dims origin)
{
    if (origin.len != PyArray_NDIM(array)) {
        PyErr_Format(PyExc_ValueError,
                     "invalid %d element 'origin' sequence for "
                     "%d-dimensional input array",
                     origin.len, PyArray_NDIM(array));
        return 0;
    }
    return 1;
}

#endif  /* NI_CONVERTERS_H */
