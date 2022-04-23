from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.pycapsule cimport (
    PyCapsule_New, PyCapsule_SetContext, PyCapsule_GetContext, PyCapsule_GetPointer
)

cimport numpy as np
from numpy cimport npy_intp as intp

np.import_array()

cdef void _destructor(obj):
    cdef void *callback_data = PyCapsule_GetContext(obj)
    PyMem_Free(callback_data)


cdef void _destructor_data(obj):
    cdef void *callback_data = PyCapsule_GetPointer(obj, NULL)
    PyMem_Free(callback_data)


cdef int _filter1d(double *input_line, intp input_length, double *output_line,
	           intp output_length, void *callback_data):
    cdef intp i, j
    cdef intp filter_size = (<intp *>callback_data)[0]

    for i in range(output_length):
        output_line[i] = 0
        for j in range(filter_size):
            output_line[i] += input_line[i+j]
        output_line[i] /= filter_size
    return 1


def filter1d(intp filter_size, with_signature=False):
    cdef intp *callback_data = <intp *>PyMem_Malloc(sizeof(intp))
    cdef char *signature = NULL
    if not callback_data:
        raise MemoryError()
    callback_data[0] = filter_size

    if with_signature:
        signature = "int (double *, npy_intp, double *, npy_intp, void *)"

    try:
        capsule = PyCapsule_New(<void *>_filter1d, signature, _destructor)
        res = PyCapsule_SetContext(capsule, callback_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule


def filter1d_capsule(intp filter_size):
    cdef intp *callback_data = <intp *>PyMem_Malloc(sizeof(intp))
    if not callback_data:
        raise MemoryError()
    callback_data[0] = filter_size

    try:
        capsule = PyCapsule_New(<void *>callback_data, NULL, _destructor_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule


cdef int _filter2d(double *buffer, intp filter_size, double *res,
	           void *callback_data):
    cdef intp i
    cdef double *weights = <double *>callback_data

    res[0] = 0
    for i in range(filter_size):
        res[0] += weights[i]*buffer[i]
    return 1


def filter2d(seq, with_signature=False):
    cdef double *callback_data = <double *>PyMem_Malloc(len(seq)*sizeof(double))
    cdef char *signature = NULL
    if not callback_data:
        raise MemoryError()
    for i, item in enumerate(seq):
        callback_data[i] = float(item)

    if with_signature:
        signature = "int (double *, npy_intp, double *, void *)"

    try:
        capsule = PyCapsule_New(<void *>_filter2d, signature, _destructor)
        PyCapsule_SetContext(capsule, callback_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule


def filter2d_capsule(seq):
    cdef double *callback_data = <double *>PyMem_Malloc(len(seq)*sizeof(double))
    if not callback_data:
        raise MemoryError()
    for i, item in enumerate(seq):
        callback_data[i] = float(item)

    try:
        capsule = PyCapsule_New(<void *>callback_data, NULL, _destructor_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule


cdef int _transform(intp *output_coordinates, double *input_coordinates,
	            int output_rank, int input_rank, void *callback_data):
    cdef intp i
    cdef double shift = (<double *>callback_data)[0]

    for i in range(input_rank):
        input_coordinates[i] = output_coordinates[i] - shift
    return 1


def transform(double shift, with_signature=False):
    cdef double *callback_data = <double *>PyMem_Malloc(sizeof(double))
    cdef char *signature = NULL
    if not callback_data:
        raise MemoryError()
    callback_data[0] = shift

    if with_signature:
        signature = "int (npy_intp *, double *, int, int, void *)"

    try:
        capsule = PyCapsule_New(<void *>_transform, signature, _destructor)
        PyCapsule_SetContext(capsule, callback_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule


def transform_capsule(double shift):
    cdef double *callback_data = <double *>PyMem_Malloc(sizeof(double))
    if not callback_data:
        raise MemoryError()
    callback_data[0] = shift

    try:
        capsule = PyCapsule_New(<void *>callback_data, NULL, _destructor_data)
    except:  # noqa: E722
        PyMem_Free(callback_data)
        raise
    return capsule
