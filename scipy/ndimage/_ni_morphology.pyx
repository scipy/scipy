#Input arrays nrrd to be matched from the original function
#Byteswapping needs to be used.

######################################################################
# Cython version of ndimage.morphology file
######################################################################
cimport cython
from cython cimport sizeof
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

np.import_array()

DEF DISTANCE_EUCLIDIAN = 1
DEF DISTANCE_CITY_BLOCK = 2
DEF NI_DISTANCE_CHESSBOARD = 3

cdef extern from 'float.h':
    double DBL_MAX

cdef extern from "limits.h":
    long long LONG_MAX
    unsigned int UINT_MAX

cdef extern from *:
   ctypedef int Py_intptr_t

cdef extern from "numpy/arrayobject.h" nogil:
    ctypedef struct PyArrayIterObject:
        np.intp_t *coordinates
        char *dataptr
        np.npy_bool contiguous

    void *PyDataMem_NEW(size_t)
    void PyDataMem_FREE(void *)

##########################################################################
# Cython version of DistanceTransformBruteForce function
##########################################################################

ctypedef struct border_element:
    np.intp_t *coordinates
    np.intp_t index
    void *next

cpdef _distance_transform_bf(np.ndarray input, 	int metric, np.ndarray sampling_arr, 
                                  np.ndarray distances, np.ndarray features):
    cdef:
        np.intp_t size, jj, min_index = 0
        int kk
        border_element *border_elements = NULL, *temp
        np.flatiter _ii, _di, _fi
        PyArrayIterObject *ii, *di, *fi
        np.float64_t *sampling = <np.float64_t *> np.PyArray_DATA(sampling_arr) if sampling_arr is not None else NULL

    size = input.size

    # Iterator Intialization
    # Making the contiguous is important as without it the arrays 
    # behave as flat arrays.

    if distances:
        _di = np.PyArray_IterNew(distances)
        di = <PyArrayIterObject *> _di
        di.contiguous = 0
        
    if features:
        _fi = np.PyArray_IterNew(features)
        fi = <PyArrayIterObject *> _fi
        fi.contiguous = 0

    if input:
        _ii = np.PyArray_IterNew(input)
        ii = <PyArrayIterObject *> _ii
        ii.contiguous = 0

    for jj in range(size):
    	# take this input from a functions
        if (<np.int8_t *> np.PyArray_ITER_DATA(_ii))[0] < 0:
            temp = <border_element *> PyDataMem_NEW(sizeof(border_element))
            if temp == NULL:
                raise MemoryError()

            temp.next = border_elements
            border_elements = temp
            temp.index = jj
            temp.coordinates = <np.intp_t *> PyDataMem_NEW(input.ndim * sizeof(np.intp_t))
            for kk in range(input.ndim):
                temp.coordinates[kk] = ii.coordinates[kk]
        np.PyArray_ITER_NEXT(_ii)

    np.PyArray_ITER_RESET(_ii)
    
    cdef double distance, d, t
    if metric == DISTANCE_EUCLIDIAN:
        for jj in range(size):
            if (<np.int8_t *> np.PyArray_ITER_DATA(_ii))[0] > 0:
                distance = DBL_MAX
                temp = border_elements
                while temp:
                    d = 0.0
                    for kk in range(input.ndim):
                        t = ii.coordinates[kk] - temp.coordinates[kk]
                        if sampling is not NULL:
                            t *= sampling[kk]
                        d += t * t

                    if d < distance:
                        distance = d
                        if features is not None:
                            min_index = temp.index

                    temp = <border_element *>temp.next
                    # print 'reading data', (<np.int8_t *> np.PyArray_ITER_DATA(_ii))[0]

                if distances:
                    (<np.float64_t *>np.PyArray_ITER_DATA(_di) )[0] = sqrt(distance)
                if features:
                    (<np.intp_t *>np.PyArray_ITER_DATA(_fi))[0]= min_index
            else:
                if distances:
                    (<np.float64_t*>np.PyArray_ITER_DATA(_di) )[0] = 0.0
                if features:
                    (<np.int32_t *> np.PyArray_ITER_DATA(_fi))[0] = jj
            
            np.PyArray_ITER_NEXT(_ii)
            if features is not None and distances is not None:
                # NI_ITERATOR_NEXT3ii, di, fi, pi, np.PyArray_ITER_DATA(_di), pf
                np.PyArray_ITER_NEXT(_di)
                np.PyArray_ITER_NEXT(_fi)
            elif distances is not None:
                np.PyArray_ITER_NEXT(_di)
            else:
                np.PyArray_ITER_NEXT(_fi)

    
    cdef unsigned int distance_n, d_n
    cdef np.intp_t t_n
    if metric == NI_DISTANCE_CHESSBOARD or metric == DISTANCE_CITY_BLOCK:
        for jj in range(size):
            if (<np.int8_t *> np.PyArray_ITER_DATA(_ii))[0] > 0:
                distance_n = UINT_MAX
                temp = border_elements
                while temp:
                    d_n = 0
                    for kk in range(input.ndim):
                        t_n = ii.coordinates[kk] - temp.coordinates[kk]
                        if t_n < 0:
                            t_n = -t_n
                        if metric == DISTANCE_CITY_BLOCK:
                            d_n += t_n
                        else:
                            if <np.uintp_t> t_n > d_n:
                                d_n = t_n
                    # print 'd_n = ', d_n, distance_n
                    if d_n < distance_n:
                        distance_n = d_n
                        if features is not None:
                            min_index = temp.index

                    temp = <border_element *> temp.next

                if distances is not None:
                   (<np.uint32_t *> np.PyArray_ITER_DATA(_di))[0] = distance_n
                if features is not None:
                    (<np.int32_t *> np.PyArray_ITER_DATA(_fi))[0] = min_index
            else:
                if distances is not None:
                    (<np.uint32_t *> np.PyArray_ITER_DATA(_di))[0] = 0
                if features is not None:
                    (<np.int32_t *> np.PyArray_ITER_DATA(_fi))[0] = jj

            # print 'Hi FInally', (<np.uint32_t *> np.PyArray_ITER_DATA(_di))[0]

            if features is not None and distances is not None:                # NI_ITERATOR_NEXT3ii, di, fi, pi, pd, pf
                np.PyArray_ITER_NEXT(_ii)
                np.PyArray_ITER_NEXT(_di)
                np.PyArray_ITER_NEXT(_fi)	
            elif distances is not None:
                # NI_ITERATOR_NEXT2ii, di, pi, pd
                np.PyArray_ITER_NEXT(_ii)
                np.PyArray_ITER_NEXT(_di)
            
            else:
                # NI_ITERATOR_NEXT2ii, fi, pi, pf
                np.PyArray_ITER_NEXT(_ii)
                np.PyArray_ITER_NEXT(_fi)

    # else:
    #     PyErr_SetString(PyExc_RuntimeError,  "distance metric not supported")

    while border_elements:
        temp = border_elements
        border_elements = <border_element *> border_elements.next
        PyDataMem_FREE(temp.coordinates)
        PyDataMem_FREE(temp)
