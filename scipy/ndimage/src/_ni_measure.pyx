######################################################################
# Cython version of scipy.ndimage.measurement.find_objects() function.
######################################################################

cimport cython
from cython cimport sizeof
import numpy as np
cimport numpy as np

np.import_array()

cdef extern from *:
   ctypedef int Py_intptr_t

ctypedef void (*PyArray_CopySwapFunc)(void *, void *, int, void *)

cdef extern from "numpy/arrayobject.h" nogil:
    ctypedef struct PyArrayIterObject:
        np.intp_t *coordinates
        char *dataptr
        np.npy_bool contiguous

    ctypedef struct PyArrayObject:
        int nd

    void *PyDataMem_NEW(size_t)
    void PyDataMem_FREE(void *)


######################################################################
# Use of Cython's type templates for type specialization
######################################################################

# Only integer values are allowed.
ctypedef fused data_t:
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.uint64_t

#####################################################################
# Function Specializers and asociate function for using fused type
#####################################################################

ctypedef void (*func_p)(void *data, np.flatiter iti, np.intp_t max_label, 
               int *regions, int rank) nogil

def get_funcs(np.ndarray[data_t] input):
    return (<Py_intptr_t> findObjectsPoint[data_t])

######################################################################
# Update Regions According to Input Data Type
######################################################################

cdef inline findObjectsPoint(data_t *data, np.flatiter _iti, np.intp_t max_label, 
            int *regions, int rank):
    cdef int kk =0
    cdef np.intp_t cc

    # only integer or boolean values are allowed, since s_index is being used in indexing
    cdef np.intp_t s_index = <np.intp_t> (<data_t *> data)[0] - 1

    if s_index >=0  and s_index < max_label:
        if rank > 0:
            s_index *= 2 * rank
            if regions[s_index] < 0:
                for kk in range(rank):
                    cc = (<PyArrayIterObject *>_iti).coordinates[kk]
                    regions[s_index + kk] = cc
                    regions[s_index + kk + rank] = cc + 1

            else:
                for kk in range(rank):
                    cc = (<PyArrayIterObject *>_iti).coordinates[kk]
                    if cc < regions[s_index + kk]:
                        regions[s_index + kk] = cc
                    if cc + 1 > regions[s_index + kk + rank]:
                        regions[s_index + kk + rank] = cc + 1
        else:
            regions[s_index] = 1

    return 1

######################################################################
# Implementaion of find_Objects function:-
######################################################################

cpdef _findObjects(np.ndarray input_raw, np.intp_t max_label):
    cdef:
        funcs = get_funcs(input_raw.take([0]))

        int ii, rank, size_regions
        int start, jj, idx, end
        int *regions = NULL

        np.flatiter _iti
        PyArrayIterObject *iti

        func_p findObjectsPoint = <func_p> <void *> <Py_intptr_t> funcs

        int flags = np.NPY_CONTIGUOUS | np.NPY_NOTSWAPPED | np.NPY_ALIGNED

    input = np.PyArray_FROM_OF(input_raw, flags)

    rank = input.ndim

    # Array Iterator defining and Initialization:
    _iti = np.PyArray_IterNew(input)
    iti = <PyArrayIterObject *> _iti

    # If iterator is contiguous, PyArray_ITER_NEXT will treat it as 1D Array
    iti.contiguous = 0
 
    # Declaring output array
    size_regions = 0

    if max_label < 0:
        max_label = 0
    elif max_label > 0:
        if rank > 0:
            size_regions = 2 * max_label * rank

        else:
            size_regions = max_label
        
        regions = <int *> PyDataMem_NEW(size_regions * sizeof(int))
        if regions == NULL:
            raise MemoryError()

    else:
        regions = NULL
    
    try:
        with nogil:
            if rank > 0:
                for jj in range(size_regions):
                    regions[jj] = -1

            # Iteration over array:
            while np.PyArray_ITER_NOTDONE(_iti):
                findObjectsPoint(np.PyArray_ITER_DATA(_iti), _iti, max_label, regions, rank)
                np.PyArray_ITER_NEXT(_iti)

        result = []

        for ii in range(max_label):
            if rank > 0:
                idx = 2 * rank * ii

            else:
                idx = ii

            if regions[idx] >= 0:
                slc = ()
                for jj in range(rank):
                    start = regions[idx + jj]
                    end = regions[idx + jj + rank]

                    slc += (slice(start, end),)
                result.append(slc)

            else:
                result.append(None)
    except:
        # clean up and re-raise
        PyDataMem_FREE(regions)
        raise

    return result
