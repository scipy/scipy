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

cdef extern from "numpy/arrayobject.h" nogil:
    ctypedef struct PyArrayIterObject:
        np.intp_t *coordinates
        char *dataptr
        np.npy_bool contiguous

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

cdef inline np.ndarray arr_corr(np.ndarray raw_array, int flags):
    cdef np.ndarray array
    if raw_array is not None:
        array = np.PyArray_FROM_OF(raw_array, flags)
    else:
        array = raw_array
    return array

cdef inline findObjectsPoint(data_t *data, np.flatiter _iti, np.intp_t max_label, 
                             int *regions, int rank):
    cdef int kk =0
    cdef np.intp_t cc

    # only integer or boolean values are allowed, since s_index is being used in indexing
    cdef np.intp_t s_index = <np.intp_t> (<data_t *> data)[0] - 1
    cdef np.intp_t val
    if s_index >= 0  and s_index < max_label:
        if rank > 0:
            s_index *= 2 * rank

            val = regions[s_index]
            for kk in range(rank):
                cc = (<PyArrayIterObject *>_iti).coordinates[kk]
                if val < 0 or cc < regions[s_index + kk]:
                    regions[s_index + kk] = cc
                if val < 0 or cc + 1 > regions[s_index + kk + rank]:
                    regions[s_index + kk + rank] = cc + 1

        else:
            regions[s_index] = 1

######################################################################
# Implementaion of find_Objects function:-
######################################################################

cpdef _findObjects(input_raw, m_label):
    cdef:
        funcs = get_funcs(input_raw.take([0]))

        np.intp_t ii, rank, size_regions, max_label = m_label
        np.intp_t start, jj, idx, end
        int *regions = NULL

        np.flatiter _iti
        PyArrayIterObject *iti

        func_p findObjectsPoint = <func_p> <void *> <Py_intptr_t> funcs

        int flags = np.NPY_NOTSWAPPED | np.NPY_ALIGNED
        np.ndarray input

    input = arr_corr(input_raw, flags)

    rank = input.ndim

    # Array Iterator defining and Initialization:
    _iti = np.PyArray_IterNew(input)
    iti = <PyArrayIterObject *> _iti

    # If iterator is contiguous, PyArray_ITER_NEXT will treat it as 1D Array
    iti.contiguous = 0
 
    # Declaring output array
    size_regions = 0

    if max_label > 0:
        if rank > 0:
            size_regions = 2 * max_label * rank
        else:
            size_regions = max_label
        
        regions = <int *> PyDataMem_NEW(size_regions * sizeof(int))
        if regions == NULL:
            raise MemoryError()
    else:
        max_label = 0
    
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


# ##############################################################################
# ##############################################################################
# #   Py_WatershedIFT function in cython
# ##############################################################################
# ##############################################################################


# DEF WS_MAXDIM = 7
# #############################################################################
# # Fused type delarations
# #############################################################################

# ctypedef fused data_markers:
#     np.int8_t
#     np.int16_t
#     np.int32_t
#     np.int64_t
#     np.uint8_t
#     np.uint16_t
#     np.uint32_t
#     np.uint64_t


# ctypedef fused data_output:
#     np.int8_t
#     np.int16_t
#     np.int32_t
#     np.int64_t
#     np.uint8_t
#     np.uint16_t
#     np.uint32_t
#     np.uint64_t

# def get_funcp_watershed(np.ndarray[data_markers] markers, np.ndarray[data_output] output):
#     return (<Py_intptr_t> case_windex1[data_output],
#             <Py_intptr_t> markers_to_output[data_markers, data_output],
#             <Py_intptr_t> case_windex2[data_output])

# ################################################################################
# # Pointer Declaration:
# #################################################################################

# ctypedef np.intp_t (*funcp_markers_to_output)(void * data_m, void * data_o, np.flatiter _mi, 
#                     np.flatiter _li)

# ctypedef np.intp_t (*funcp_windex1)(void *data, np.intp_t v_index, np.intp_t p_index, np.intp_t *strides, int i_contiguous, 
#                             np.intp_t p_idx, np.intp_t v_idx, np.flatiter _ii, np.intp_t *vval, np.intp_t *pval, np.ndarray input)

# ctypedef np.intp_t (*funcp_windex2) (void *data, np.intp_t v_index, np.intp_t p_index, np.intp_t *strides, np.ndarray output, np.ndarray input, 
#     np.intp_t idx, int o_contiguous, int *label)

# ################################################################################
# # Asociate function Declaration
# #################################################################################

# cdef void get_index(data_output *data, np.intp_t index, np.intp_t *strides, np.ndarray array, 
#                     np.intp_t rank, np.intp_t *out, np.intp_t contiguous):
#     cdef int qq
#     cdef np.intp_t cc, idx
    
    
#     if contiguous:
#         out[0] = index * sizeof(data_output)
#     else:
#         idx = index
#         out[0] = 0
#         for qq in range(rank):
#             cc = idx / strides[qq]
#             idx -= cc * strides[qq]
#             out[0] += array.strides[qq] * cc


# cdef np.intp_t markers_to_output(data_markers *data_m, data_output *data_p, np.flatiter _mi, 
#                                 np.flatiter _li):
#     cdef np.intp_t temp = 4
#     temp = (<data_markers *> np.PyArray_ITER_DATA(_mi))[0]
#     (<data_output *> np.PyArray_ITER_DATA(_mi))[0] = < data_output > temp
#     return temp
    
# cdef np.intp_t case_windex1(data_output *data, np.intp_t v_index, np.intp_t p_index, np.intp_t *strides, int i_contiguous, 
#                             np.intp_t p_idx, np.intp_t v_idx, np.flatiter _ii, np.intp_t *vval, np.intp_t *pval, np.ndarray input):
    
#     # Another Function defination
#     vval[0] = (<data_output *>(np.PyArray_ITER_DATA(_ii) + v_idx))[0]
#     pval[0] = (<data_output *>(np.PyArray_ITER_DATA(_ii) + p_idx))[0]

# cdef np.intp_t case_windex2(data_output *data, np.intp_t v_index, np.intp_t p_index, np.intp_t *strides, np.ndarray output, np.ndarray input, 
#     np.intp_t idx, int o_contiguous, int *label):
#     # Another function definations
#     get_index(data, v_index, strides, output, dim, idx, o_contiguous)
#     label[0] = (<data_output *> data)[0]
#     #same function call
#     get_index(data, p_index, strides, output.strides, out)
#     (<data_output *> data)[0] = label[0]
#     # cdef int del = 0
#     return 1

# cdef inline get_value(np.ndarray array, np.flatiter iter):
#     cdef np.intp_t ival
#     if input.dtype == np.uint8_t:
#         ival[0] = (<np.uint8_t *>iter)[0]
#     elif input.dtype == np.uint16_t:
#         ival[0] = (<np.uint16_t *>iter)[0]
#     return  ival

# #############################################################################
# # Basic function
# #############################################################################

# cdef struct WatershedElement:
#     np.uintp_t index
#     np.uint16_t cost
#     void *next, *prev
#     np.uint8_t done

# cpdef int watershed_ift(np.ndarray input_raw, np.ndarray markers, np.ndarray structure, 
#                         np.ndarray output):
#     cdef int flags = np.NPY_NOTSWAPPED | np.NPY_ALIGNED

#     arr_corr(input_raw, flags, input)
#     arr_corr(markers_raw, flags, markers)
#     arr_corr(structure_arr, flags, structure)
#     arr_corr(output_arr, flags, output)
    
#     cdef:
#         funcs = get_funcp_watershed(markers.take([0]), output.take([0]))
#         int ll, jj, hh, kk, i_contiguous, o_contiguous, label
#         np.intp_t size, maxval, nneigh, ssize, ival
#         np.intp_t strides[WS_MAXDIM], coordinates[WS_MAXDIM]
#         np.intp_t *nstrides = NULL
#         bint *ps = NULL
#         np.flatiter _ii, _li, _mi
#         PyArrayIterObject *ii, *li, *mi
#         WatershedElement *temp = NULL, **first = NULL, **last = NULL
#         dim = input.ndim
    
    
#     if dim > WS_MAXDIM:
#         raise RuntimeError("Too many dimensions")

#     ssize = structure.ndim
#     size = dim

#     temp = <WatershedElement *> PyDataMem_NEW(size * sizeof(WatershedElement))
#     # Error condition

#     # Iterator inititalization
#     _ii = np.PyArray_IterNew(input)
#     _li = np.PyArray_IterNew(output)
#     _mi = np.PyArray_IterNew(markers)

#     ii = <PyArrayIterObject *> _ii
#     li = <PyArrayIterObject *> _li
#     mi = <PyArrayIterObject *> _mi

#     cdef funcp_windex1 case_windex1 = <funcp_windex1> <void *> <Py_intptr_t> funcs[0]
#     cdef funcp_markers_to_output markers_to_output = <funcp_markers_to_output> <void *> <Py_intptr_t> funcs[1]
#     cdef funcp_windex2 case_windex2 = <funcp_windex2> <void *> <Py_intptr_t> funcs[2]
    
#     maxval = 0
#     for jj in range(size):
#         # Need value in function in ival from pi using fused_typeget_value(np.PyArray_ITER_DATA(_ii), _ii, ii, input)
#         ival = get_value(_ii, input)
        
#         temp[jj].index = jj
#         temp[jj].done = 0
#         if ival > maxval:
#             maxval = ival

#         np.PyArray_ITER_NEXT(_ii)
    
#     # Allocate and initialize the storage for the queue
#     first = <WatershedElement ** >  PyDataMem_NEW((maxval + 1) * sizeof(WatershedElement *))
#     last =  <WatershedElement ** > PyDataMem_NEW((maxval + 1) * sizeof(WatershedElement *))
#     # error in allocations
    
#     for hh in range(maxval + 1):
#         first[hh] = NULL
#         last[hh] = NULL

#     for ll in range(dim):
#         coordinates[ll] = 0

#     # use memset in here

#     for jj in range(size):
#         # check funtions
#         label = markers_to_output(np.PyArray_ITER_DATA(_mi), np.PyArray_ITER_DATA(_li), _mi, _li)
#         np.PyArray_ITER_NEXT(_mi)
#         np.PyArray_ITER_NEXT(_li)
#         if label:
#             # This node is a marker
#             temp[jj].cost = 0
#             # check whether vython count NULL as 0
#             # if !first[0]:
#             if first[0] is NULL:
#                 first[0] = &(temp[jj])
#                 # beware here.. could get erreors
#                 first[0].prev = NULL
#                 first[0].next = NULL
#                 last[0] = first[0]

#             else:
#                 if label > 0:
#                     # object markers are enqueued at the beginning, so they
#                     # are processed first.
#                     temp[jj].next = first[0]
#                     temp[jj].prev = NULL
#                     first[0].prev = &(temp[jj])
#                     first[0] = &(temp[jj])

#                 else:
#                     # background markers are enqueued at the end, so they are
#                     # processed after the object markers. 
#                     temp[jj].next = NULL
#                     temp[jj].prev = last[0]
#                     last[0].next = &(temp[jj])
#                     last[0] = &(temp[jj])

#         else:
#             # This node is not a marker 
#             temp[jj].cost = maxval + 1
#             temp[jj].next = NULL
#             temp[jj].prev = NULL

#         ll = dim - 1
#         while ll >=0:
#             if coordinates[ll] < input.dimensions[ll] - 1:
#                 coordinates[ll] += 1
#                 break

#             else:
#                 coordinates[ll] = 0
#             ll -= 1

#     nneigh = 0
#     for kk in range(ssize):
#         if ps[kk] and kk != ssize/2:
#             nneigh += 1

#     nstrides = <np.intp_t *> PyDataMem_NEW(nneigh * sizeof(np.intp_t))
#     #Error in allocation

#     strides[dim -1] = 1

#     for ll in range(dim -1) :
#         strides[ll] = input.dimensions[ll + 1] * strides[ll + 1]

#     for ll in range(dim):
#         coordinates[ll] = -1

#     for kk in range(nneigh):
#         nstrides[kk] = 0

#     jj = 0
#     cdef int offset = 0

#     for kk in range(ssize):
#         if ps[kk]:
#             offset = 0
#             for ll in range(dim):
#                 offset += coordinates[ll] * strides[ll]
#             if offset:
#                 nstrides[jj] += offset
#                 jj += 1

#         ll = dim -1
#         while ll >= 0:
#             if coordinates[ll] < 1:
#                 coordinates[ll] += 1

#             else:
#                 coordinates[ll] = -1
#             ll -= 1


#     # Propogation Phase
#     cdef:
#         WatershedElement *v, *p, *prev, *next
#         np.intp_t v_index, p_index, idx, cc
#         int qq, outside
#         np.intp_t max, pval, vval, wvp, pcost, p_idx, v_idx
    
#     for jj in range(maxval + 1):
#         while first[jj] is not NULL:
#             v = first[jj]
#             first[jj] = <WatershedElement *>first[jj].next
#             if first[jj] is not NULL:
#                 (first[jj].prev)[0] = NULL

#             (v[0]).prev = NULL
#             (v[0]).next = NULL
#             # Mark elements as Done:
#             (v[0]).done = 1

#             for hh in range(nneigh):
#                 v_index = v.index
#                 p_index = v.index + nstrides[hh]
#                 outside = 0
#                 # Check if the neighbour is within the extent of the array
#                 idx = p_index
#                 for qq in range(dim):
#                     cc = idx / strides[qq]
#                     if cc < 0 or cc >=input.dimensions[qq]:
#                         outside = 1
#                         break
#                     idx -= cc * strides[qq]

#                 if outside:
#                     p = &(temp[p_index])
#                     # If the neighbour is not processed Yet
#                     if not ((p[0]).done):

#                         case_windex1(np.PyArray_ITER_DATA(_ii), v_index, p_index, strides,
#                                      i_contiguous, p_idx, v_idx, _ii, &vval, &pval, input)


# # IMPLEMENT THIS FUNCTION USING IF and ELSE condition.



#                         # Calculate Cost
#                         wvp = pval - vval
#                         if wvp < 0:
#                             wvp = -wvp
#                         # Find the maximum of this cost and the current element cost
#                         pcost = (p[0]).cost
#                         if (v[0]).cost > wvp:
#                             max = (v[0]).cost
#                         else:
#                             max = wvp

#                         if max < pcost:
#                             # If this maximum is less than the neighbors cost,
#                             # adapt the cost and the label of the neighbor: 
#                             (p[0]).cost = max
#                             case_windex2(np.PyArray_ITER_DATA(_li), v_index, p_index, strides, output, input, idx,
#                                         o_contiguous, &label)
                            
#                             # If the neighbor is in a queue, remove it:
#                             if p.next is not NULL or p.prev is not NULL:
#                                 prev = <WatershedElement *> p.prev
#                                 next = <WatershedElement *> p.next
#                                 if first[pcost] == p:
#                                     first[pcost] = next

#                                 if last[pcost] == p:
#                                     last[pcost] = prev

#                                 if prev is not NULL:
#                                     prev.next = next

#                                 if next is not NULL:
#                                     next.prev = prev

#                             # Insert the neighbor in the appropiate queue:
#                             if label < 0:
#                                 p[0].prev = last[max]
#                                 p[0].next = NULL
#                                 if last[max] is not NULL:
#                                     ((last[max])[0]).next = p
#                                 last[max] = p
#                                 if first[max] is NULL:
#                                     first[max] = p

#                             else:
#                                 (p[0]).next = first[max]
#                                 p[0].prev = NULL
#                                 if first[max] is not NULL:
#                                     ((first[max])[0]).prev = p
#                                 first[max] = p
#                                 if last[max] is NULL:
#                                     last[max] = p
    
#     PyDataMem_FREE(temp)
#     PyDataMem_FREE(first)
#     PyDataMem_FREE(last)
#     PyDataMem_FREE(nstrides)

#     return 1
