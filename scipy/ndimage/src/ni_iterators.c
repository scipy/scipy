#include "ni_iterators.h"
#include "ni_iterator_utils.h"


/*
 ***********************************************************************
 ***                           NDIterator                            ***
 ***********************************************************************
 */

struct NDIterator_Internal {
    /* Number of dimensions being iterated minus one. */
    int ndim_1;
    /* Total iterator size */
    npy_intp size;
    /* Sizes of the dimensions being iterated, minus 1 */
    npy_intp dims_1[NPY_MAXDIMS];
    /* Strides over the dimensions being iterated */
    npy_intp strides[NPY_MAXDIMS];
    /* How far to jump back for the dimensions being iterated */
    npy_intp backstrides[NPY_MAXDIMS];
    /* Current iterator position  */
    npy_intp index;
    /* Current iterator position for multidimensional iteration */
    npy_intp coords[NPY_MAXDIMS];
    /* Pointer to the current iteration item */
    char *data;
    /* The array being iterated */
    PyArrayObject *array;
};


static NDIterator*
_NDI_NewExceptAxes(PyArrayObject *array, npy_uint32 axes)
{
    NDIterator *nditer = malloc(sizeof(NDIterator));
    int axis;

    if (nditer == NULL) {
        return (NDIterator *)PyErr_NoMemory();
    }

    Py_INCREF(array);
    nditer->array = array;
    nditer->data = PyArray_BYTES(array);
    nditer->ndim_1 = -1;
    nditer->size = 1;
    for (axis = 0; axis < PyArray_NDIM(array); ++axis) {
        if (!(axes & (npy_uint32)0x1)) {
            const npy_intp dim = PyArray_DIM(array, axis);
            const npy_intp stride = PyArray_STRIDE(array, axis);
            const int iter_axis = ++(nditer->ndim_1);

            nditer->size *= dim;
            nditer->dims_1[iter_axis] = dim - 1;
            nditer->strides[iter_axis] = stride;
            nditer->backstrides[iter_axis] = stride * (dim - 1);
            nditer->coords[iter_axis] = 0;
        }
        axes >>= 1;
    }
    nditer->index = 0;

    return nditer;
}

NDIterator*
NDI_New(PyArrayObject *array)
{
    if (!_is_input_array(array)) {
        return NULL;
    }

    return _NDI_NewExceptAxes(array, (npy_uint32)0x0);
}


NDIterator*
NDI_NewExceptAxis(PyArrayObject *array, int axis)
{
    if (!_is_input_array(array) || !_is_valid_axis(array, &axis)) {
        return NULL;
    }

    return _NDI_NewExceptAxes(array, ((npy_uint32)0x1) << axis);
}


int
NDI_Next(NDIterator *nditer)
{
    int axis;

    if (++(nditer->index) >= nditer->size) {
        return 0;
    }

    for (axis = nditer->ndim_1; axis >= 0; --axis) {
        if (nditer->coords[axis] < nditer->dims_1[axis]) {
            ++(nditer->coords[axis]);
            nditer->data += nditer->strides[axis];
            break;
        }
        else {
            nditer->coords[axis] = 0;
            nditer->data -= nditer->backstrides[axis];
        }
    }

    return 1;
}

static int
_NDI_NextForTwo(NDIterator *nditer1, NDIterator *nditer2)
{
    int axis;

    ++(nditer1->index);
    ++(nditer2->index);
    if (nditer1->index >= nditer1->size) {
        return 0;
    }

    for (axis = nditer1->ndim_1; axis >= 0; --axis) {
        if (nditer1->coords[axis] < nditer1->dims_1[axis]) {
            ++(nditer1->coords[axis]);
            ++(nditer2->coords[axis]);
            nditer1->data += nditer1->strides[axis];
            nditer2->data += nditer2->strides[axis];
            break;
        }
        else {
            nditer1->coords[axis] = 0;
            nditer2->coords[axis] = 0;
            nditer1->data -= nditer1->backstrides[axis];
            nditer2->data -= nditer2->backstrides[axis];
        }
    }

    return 1;
}


void
NDI_Delete(NDIterator *nditer)
{
    if (nditer != NULL) {
        Py_DECREF(nditer->array);
        free(nditer);
    }
}


/*
 ***********************************************************************
 ***                       LineBufferIterator                        ***
 ***********************************************************************
 */

struct LineBufferIterator_Internal {
    /* Input and output buffers. */
    npy_double *in_buffer;
    npy_double *out_buffer;
    /* Size of filter. */
    npy_intp filter_size;
    /* Size of filter before and after origin, size = before + 1 + after. */
    npy_intp size_before;
    npy_intp size_after;
    /* Length of array lines. */
    npy_intp line_length;
    /* Iterators over input and output array lines. */
    NDIterator *in_nditer;
    NDIterator *out_nditer;
    /* Line strides of the input and output arrays. */
    npy_intp in_stride;
    npy_intp out_stride;
    /* Functions to read and write array lines from/to buffers. */
    read_line_func *read_func;
    write_line_func *write_func;
    /* Extend function and value for constant mode. */
    extend_line_func *extend_func;
    npy_double extend_value;
};


static NPY_INLINE void
_LBI_ExtendLine(LineBufferIterator *lbiter)
{
    lbiter->extend_func(lbiter->in_buffer, lbiter->line_length,
                        lbiter->size_before, lbiter->size_after,
                        lbiter->extend_value);
}


static NPY_INLINE void
_LBI_ReadLine(LineBufferIterator *lbiter)
{
    lbiter->read_func(lbiter->in_buffer + lbiter->size_before,
                      lbiter->in_nditer->data,
                      lbiter->line_length, lbiter->in_stride);
    _LBI_ExtendLine(lbiter);
}


LineBufferIterator*
LBI_New(PyArrayObject *input, int axis, PyArrayObject *output,
        npy_intp filter_size, npy_intp filter_offset,
        NI_ExtendMode extend_mode, npy_double extend_value)
{
    LineBufferIterator *lbiter = NULL;

    if (!_is_input_array(input) || !_is_output_array(output) ||
            !_arrays_are_compatible(input, output) ||
            !_is_valid_axis(input, &axis)) {
        return NULL;
    }

    lbiter = malloc(sizeof(LineBufferIterator));
    if (lbiter == NULL) {
        return (LineBufferIterator *)PyErr_NoMemory();
    }

    lbiter->in_buffer = NULL;
    lbiter->out_buffer = NULL;
    lbiter->in_nditer = NULL;
    lbiter->out_nditer = NULL;

    if (!_split_filter_size(filter_size, filter_offset,
                            &lbiter->size_before, &lbiter->size_after)) {
        return LBI_Delete(lbiter);
    }
    lbiter->filter_size = filter_size;

    lbiter->read_func = _get_read_line_func(PyArray_TYPE(input));
    if (lbiter->read_func == NULL) {
        return LBI_Delete(lbiter);
    }
    lbiter->in_nditer = NDI_NewExceptAxis(input, axis);
    if (lbiter->in_nditer == NULL) {
        return LBI_Delete(lbiter);
    }

    lbiter->write_func = _get_write_line_func(PyArray_TYPE(output));
    if (lbiter->write_func == NULL) {
        return LBI_Delete(lbiter);
    }
    lbiter->out_nditer = NDI_NewExceptAxis(output, axis);
    if (lbiter->out_nditer == NULL) {
        return LBI_Delete(lbiter);
    }

    lbiter->line_length = PyArray_DIM(input, axis);
    lbiter->in_stride = PyArray_STRIDE(input, axis);
    lbiter->out_stride = PyArray_STRIDE(output, axis);

    lbiter->in_buffer = malloc(sizeof(npy_double) *
            (lbiter->line_length + lbiter->size_before + lbiter->size_after));
    if (lbiter->in_buffer == NULL) {
        PyErr_NoMemory();
        return LBI_Delete(lbiter);
    }
    lbiter->out_buffer = malloc(sizeof(npy_double) * lbiter->line_length);
    if (lbiter->out_buffer == NULL) {
        PyErr_NoMemory();
        return LBI_Delete(lbiter);
    }

    lbiter->extend_func = _get_extend_line_func(extend_mode);
    if (lbiter->extend_func == NULL) {
        return LBI_Delete(lbiter);
    }
    lbiter->extend_value = extend_value;

    _LBI_ReadLine(lbiter);

    return lbiter;
}


npy_double*
LBI_GetInputBuffer(LineBufferIterator *lbiter)
{
    return lbiter->in_buffer;
}


npy_double*
LBI_GetOutputBuffer(LineBufferIterator *lbiter)
{
    return lbiter->out_buffer;
}

npy_intp
LBI_GetLineLength(LineBufferIterator *lbiter)
{
    return lbiter->line_length;
}


static NPY_INLINE void
_LBI_WriteLine(LineBufferIterator *lbiter)
{
    lbiter->write_func(lbiter->out_buffer, lbiter->out_nditer->data,
                       lbiter->line_length, lbiter->out_stride);
}


int
LBI_Next(LineBufferIterator *lbiter)
{
    _LBI_WriteLine(lbiter);
    if (_NDI_NextForTwo(lbiter->in_nditer, lbiter->out_nditer)) {
        _LBI_ReadLine(lbiter);
        return 1;
    }
    return 0;
}


LineBufferIterator*
LBI_Delete(LineBufferIterator *lbiter)
{
    if (lbiter != NULL) {
        NDI_Delete(lbiter->in_nditer);
        NDI_Delete(lbiter->out_nditer);
        free(lbiter->in_buffer);
        free(lbiter->out_buffer);
        free(lbiter);
    }
    return NULL;
}
