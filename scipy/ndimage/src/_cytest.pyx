from libc.stdlib cimport malloc, free
from cpython.pycapsule cimport (
    PyCapsule_New, PyCapsule_SetContext, PyCapsule_GetContext
)

cimport numpy as np
from numpy cimport npy_intp as intp


cdef void _destructor(obj):
    cdef void *callback_data = PyCapsule_GetContext(obj)
    free(callback_data)


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
    

def filter1d(intp filter_size):
    cdef intp *callback_data = <intp *>malloc(sizeof(intp))
    callback_data[0] = filter_size
    capsule = PyCapsule_New(<void *>_filter1d, NULL, _destructor)
    PyCapsule_SetContext(capsule, callback_data)
    return capsule


cdef int _filter2d(double *buffer, intp filter_size, double *res,
	           void *callback_data):
    cdef intp i
    cdef double *weights = <double *>callback_data

    res[0] = 0
    for i in range(filter_size):
        res[0] += weights[i]*buffer[i]
    return 1
    
    
def filter2d(seq):
    cdef double *callback_data = <double *>malloc(len(seq)*sizeof(double))
    for i, item in enumerate(seq):
        callback_data[i] = float(item)

    capsule = PyCapsule_New(<void *>_filter2d, NULL, _destructor)
    PyCapsule_SetContext(capsule, callback_data)
    return capsule


cdef int _transform(intp *output_coordinates, double *input_coordinates,
	            intp output_rank, intp input_rank, void *callback_data):
    cdef intp i
    cdef double shift = (<double *>callback_data)[0]

    for i in range(input_rank):
        input_coordinates[i] = output_coordinates[i] - shift
    return 1


def transform(double shift):
    cdef double *callback_data = <double *>malloc(sizeof(double))
    callback_data[0] = shift
    capsule = PyCapsule_New(<void *>_transform, NULL, _destructor)
    PyCapsule_SetContext(capsule, callback_data)
    return capsule
