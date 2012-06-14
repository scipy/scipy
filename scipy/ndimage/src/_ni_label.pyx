######################################################################
# Cython version of scipy.ndimage.measurements.label().
# Requires Cython version 0.17 or greater due to type templating.
######################################################################

cimport cython
from cython cimport sizeof
import numpy as np
cimport numpy as np

np.import_array()

cdef extern from *:
   ctypedef int Py_intptr_t

cdef enum:
    BACKGROUND = 0
    FOREGROUND = 1

cdef extern from "numpy/arrayobject.h" nogil:
    ctypedef struct PyArrayIterObject:
        np.npy_intp *coordinates

    void PyArray_ITER_NEXT(PyArrayIterObject *it)
    int PyArray_ITER_NOTDONE(PyArrayIterObject *it)
    void PyArray_ITER_RESET(PyArrayIterObject *it)
    void *PyArray_ITER_DATA(PyArrayIterObject *it)

    void *PyDataMem_NEW(size_t)
    void PyDataMem_FREE(void *)
    void *PyDataMem_RENEW(void *, size_t)

class NeedMoreBits(Exception):
    pass

######################################################################
# Use Cython's type templates for type specialization
######################################################################
ctypedef fused data_t:
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.uint64_t
    np.float32_t
    np.float64_t

# needed for double-specialization of copy_nonzero_line
ctypedef fused data2_t:
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.uint64_t
    np.float32_t
    np.float64_t

######################################################################
# Copy one line of the input array to the output, setting output to FOREGROUND
# wherever the input is nonzero
######################################################################
cdef void copy_nonzero_line(data_t *inl, data2_t *out,
                            np.intp_t L, np.intp_t si, np.intp_t so):
    cdef np.uint_t i
    for i in range(L):
        out[0] = FOREGROUND if (inl[0] != 0) else BACKGROUND
        inl = <data_t *> ((<void *> inl) + si)
        out = <data2_t *> ((<void *> out) + so)


######################################################################
# Mark two labels to be merged
######################################################################
cdef inline np.uintp_t mark_for_merge(np.uintp_t a,
                                      np.uintp_t b,
                                      np.uintp_t *mergetable) nogil:
    # we keep the mergetable such that merged labels always point to the
    # smallest value in the merge.
    if mergetable[a] < mergetable[b]:
        mergetable[b] = mergetable[a]
    else:
        mergetable[a] = mergetable[b]
    return mergetable[a]


######################################################################
# Take the label of a neighbor, or mark them for merging
######################################################################
cdef inline np.uintp_t take_label_or_merge(np.uintp_t neighbor_label,
                                           np.uintp_t cur_label,
                                           np.uintp_t *mergetable) nogil:
    if neighbor_label == BACKGROUND:
        return cur_label
    if cur_label == FOREGROUND:
        return neighbor_label
    if neighbor_label:
        if cur_label != neighbor_label:
            cur_label = mark_for_merge(neighbor_label, cur_label, mergetable)
    return cur_label

######################################################################
# Label one line of input, using a neighbor line that has already been labeled.
######################################################################
cdef void label_line_with_neighbor(data_t *line,
                                   data_t *neighbor,
                                   np.intp_t stride,
                                   int use_previous,
                                   int use_adjacent,
                                   int use_next,
                                   np.intp_t L,
                                   np.uintp_t *mergetable) nogil:
    cdef:
        np.uintp_t cur_label, i
        data_t *prev, *next

    # special case the first element
    next = <data_t *> ((<void *> neighbor) + stride)
    if line[0] != BACKGROUND:
        cur_label = <np.uintp_t> line[0]
        if use_adjacent:
            cur_label = take_label_or_merge(<np.uintp_t> neighbor[0], cur_label, mergetable)
        if use_next:
            cur_label = take_label_or_merge(<np.uintp_t> next[0], cur_label, mergetable)
        line[0] = cur_label
    for i in range(1, L - 1):
        # increment pointers
        line = <data_t *> ((<void *> line) + stride)
        prev = neighbor
        neighbor = next
        next = <data_t *> ((<void *> next) + stride)
        if line[0] == BACKGROUND:
            continue
        cur_label = <np.uintp_t> line[0]
        if use_previous:
            cur_label = take_label_or_merge(<np.uintp_t> prev[0], cur_label, mergetable)
        if use_adjacent:
            cur_label = take_label_or_merge(<np.uintp_t> neighbor[0], cur_label, mergetable)
        if use_next:
            cur_label = take_label_or_merge(<np.uintp_t> next[0], cur_label, mergetable)
        line[0] = cur_label
    # special case the last element
    # increment pointers
    line = <data_t *> ((<void *> line) + stride)
    prev = neighbor
    neighbor = next
    if line[0] != BACKGROUND:
        cur_label = <np.uintp_t> line[0]
        if use_previous:
            cur_label = take_label_or_merge(<np.uintp_t> prev[0], cur_label, mergetable)
        if use_adjacent:
            cur_label = take_label_or_merge(<np.uintp_t> neighbor[0], cur_label, mergetable)
        line[0] = cur_label

######################################################################
# Label one line of input after labeling it with any neighbors.
# returns updated count of regions.
######################################################################
cdef np.uintp_t label_line(data_t *line,
                           np.intp_t stride,
                           bint use_previous,
                           np.uintp_t L,
                           np.uintp_t next_label,
                           np.uintp_t *mergetable) nogil except 0:
    cdef:
        np.uintp_t cur_label, i
        data_t *prev

    # special case first element
    if line[0] != BACKGROUND:
        if <np.uintp_t> line[0] == FOREGROUND:
            # create a new region
            line[0] = <data_t> next_label
            mergetable[next_label] = next_label
            next_label += 1
            # check that we aren't going to overflow
            if next_label != <np.intp_t> <data_t> next_label:
                with gil:
                    raise NeedMoreBits()
    for i in range(1, L):
        prev = line
        line = <data_t *> ((<void *> line) + stride)
        if line[0] == BACKGROUND:
            continue
        elif line[0] == FOREGROUND:
            if use_previous and (prev[0] != BACKGROUND):
                line[0] = prev[0]
            else:
                # create a new region
                line[0] = <data_t> next_label
                mergetable[next_label] = next_label
                next_label += 1
                # check that we aren't going to overflow
                if next_label != <np.intp_t> <data_t> next_label:
                    with gil:
                        raise NeedMoreBits()
        elif use_previous and (prev[0] != BACKGROUND) and \
                line[0] != prev[0]:
            line[0] = <data_t> mark_for_merge(<np.uintp_t> line[0],
                                              <np.uintp_t> prev[0],
                                              mergetable)
    return next_label

######################################################################
# Remap region labels through the mergetable
######################################################################
cdef void remap_line(data_t *line,
                     np.intp_t stride,
                     np.uintp_t L,
                     np.uintp_t *mergetable) nogil:
    for i in range(0, L):
        line[0] = <data_t> mergetable[<np.uintp_t> line[0]]
        line = <data_t *> ((<void *> line) + stride)


######################################################################
# Function specializers
######################################################################
ctypedef void (*copy_nonzero_line_func_t)(void *inl, void *out,
                                          np.intp_t L,
                                          np.intp_t si, np.intp_t so) nogil
ctypedef void (*label_line_with_neighbor_func_t)(void *line,
                                                 void *neighbor,
                                                 np.intp_t line_stride,
                                                 int use_previous,
                                                 int use_adjacent,
                                                 int use_next,
                                                 np.intp_t L,
                                                 np.uintp_t *mergetable) nogil
ctypedef np.uintp_t (*label_line_func_t)(void *line,
                                         np.intp_t line_stride,
                                         int use_previous,
                                         np.uintp_t L,
                                         np.uintp_t next_label,
                                         np.uintp_t *mergetable) nogil
ctypedef void (*remap_line_func_t)(void *line,
                                   np.intp_t line_stride,
                                   np.uintp_t L,
                                   np.uintp_t *mergetable) nogil

def get_funcs(np.ndarray[data_t] input, np.ndarray[data2_t] output):
    return (<Py_intptr_t> copy_nonzero_line[data_t, data2_t],
            <Py_intptr_t> label_line_with_neighbor[data2_t],
            <Py_intptr_t> label_line[data2_t],
            <Py_intptr_t> remap_line[data2_t])

cpdef _label(np.ndarray input,
             np.ndarray structure,
             np.ndarray output):
    # check dimensions
    # To understand the need for the casts to object, see
    # http://trac.cython.org/cython_trac/ticket/302
    assert (<object> input).shape == (<object> output).shape, \
        ("Shapes must match for input and output,"
         "{} != {}".format((<object> input).shape, (<object> output).shape))

    structure = np.asanyarray(structure, dtype=np.int).copy()
    assert input.ndim == structure.ndim, \
        ("Structuring element must have same "
         "# of dimensions as input, "
         "{:d} != {:d}".format(input.ndim, structure.ndim))

    # Check that structuring element is of size 3 in every dimension
    assert set((<object> structure).shape) <= set([3]), \
        ("Structuring element must be size 3 in every dimension, "
         "was {}".format((<object> structure).shape))

    # check structuring element for symmetry
    assert np.all(structure == structure[(np.s_[::-1],) * structure.ndim]), \
        "Structuring element is not symmetric"

    # make sure we're dealing with a non-empty, non-scalar array
    assert input.ndim > 0 and input.size > 0, "Cannot label scalars or empty arrays"

    # if we're handed booleans, we treat them as uint8s
    if input.dtype == np.bool:
        input = input.view(dtype=np.uint8)
    if output.dtype == np.bool:
        # XXX - trigger special check for bit depth?
        output = output.view(dtype=np.uint8)

    funcs = get_funcs(input.take([0]), output.take([0]))
    cdef:
        copy_nonzero_line_func_t copy_nonzero_line = \
            <copy_nonzero_line_func_t> <void *> <Py_intptr_t> funcs[0]
        label_line_with_neighbor_func_t label_line_with_neighbor = \
            <label_line_with_neighbor_func_t> <void *> <Py_intptr_t> funcs[1]
        label_line_func_t label_line = \
            <label_line_func_t> <void *> <Py_intptr_t> funcs[2]
        remap_line_func_t remap_line = \
            <remap_line_func_t> <void *> <Py_intptr_t> funcs[3]
        np.uintp_t next_label, src_label, dest_label
        np.flatiter _iti, _ito, _itstruct
        PyArrayIterObject *iti, *ito, *itstruct
        int axis, L, delta
        np.intp_t si, so, ss
        np.intp_t total_offset
        bint valid, center, use_previous
        int mergetable_size
        np.uintp_t *mergetable

    # make sure center of structuring element is empty, and everything past it,
    # which are lines we haven't processed, yet.
    structure.flat[structure.size // 2] = 0

    axis = -1  # choose best axis based on output
    _ito = np.PyArray_IterAllButAxis(output, &axis)
    _iti = np.PyArray_IterAllButAxis(input, &axis)
    _itstruct = np.PyArray_IterAllButAxis(structure, &axis)

    ito = <PyArrayIterObject *> _ito
    iti = <PyArrayIterObject *> _iti
    itstruct = <PyArrayIterObject *> _itstruct

    mergetable_size = output.shape[axis]
    mergetable = <np.uintp_t *> PyDataMem_NEW(mergetable_size * sizeof(np.uintp_t))
    if mergetable == NULL:
        raise MemoryError()

    try:
        L = input.shape[axis]
        # strides
        si = input.strides[axis]
        so = output.strides[axis]
        ss = structure.strides[axis]

        # 0 = background
        # 1 = foreground, needs label
        # 2... = working labels, will be compacted on output
        next_label = 2

        # Used for labeling single lines
        temp = [1] * structure.ndim
        temp[axis] = 0
        use_previous = (structure[tuple(temp)] != 0)

        with nogil:
            while PyArray_ITER_NOTDONE(iti):
                # nonzero input -> output
                copy_nonzero_line(PyArray_ITER_DATA(iti), PyArray_ITER_DATA(ito), L, si, so)

                # Take neighbor labels
                PyArray_ITER_RESET(itstruct)
                while PyArray_ITER_NOTDONE(itstruct):
                    if not ((<np.int_t *> PyArray_ITER_DATA(itstruct))[0] or
                            (<np.int_t *> (PyArray_ITER_DATA(itstruct) + ss))[0] or
                            (<np.int_t *> (PyArray_ITER_DATA(itstruct) + 2 * ss))[0]):
                        PyArray_ITER_NEXT(itstruct)
                        continue
                    valid = True
                    for idim in range(structure.ndim):
                        if idim == axis:
                            continue
                        delta = (itstruct.coordinates[idim] - 1)  # 1,1,1... is center
                        if delta != 0:
                            if not (0 <= (ito.coordinates[idim] + delta) < output.shape[idim]):
                                valid = False
                                break
                    if valid:
                        total_offset = 0
                        for idim in range(structure.ndim):
                            if idim == axis:
                                continue
                            delta = (itstruct.coordinates[idim] - 1)
                            total_offset += delta * output.strides[idim]
                        if total_offset == 0:  # don't treat ourselves as a neighbor
                            # this indicates we've hit the center, and should stop
                            break
                        else:
                            label_line_with_neighbor(PyArray_ITER_DATA(ito),
                                                     PyArray_ITER_DATA(ito) + total_offset,
                                                     so,
                                                     (<np.int_t *> (PyArray_ITER_DATA(itstruct)))[0],
                                                     (<np.int_t *> (PyArray_ITER_DATA(itstruct) + ss))[0],
                                                     (<np.int_t *> (PyArray_ITER_DATA(itstruct) + 2 * ss))[0],
                                                     L,
                                                     mergetable)
                    PyArray_ITER_NEXT(itstruct)

                # Label unlabeled pixels
                # be conservative about how much space we may need
                while mergetable_size < (next_label + L):
                    mergetable_size *= 2
                    mergetable = <np.uintp_t *> PyDataMem_RENEW(<void *> mergetable,
                                                                 mergetable_size * sizeof(np.uintp_t))
                next_label = label_line(PyArray_ITER_DATA(ito),
                                        so, use_previous, L,
                                        next_label, mergetable)
                PyArray_ITER_NEXT(iti)
                PyArray_ITER_NEXT(ito)

            # compact the mergetable, mapping each value to its destination
            mergetable[BACKGROUND] = BACKGROUND
            mergetable[FOREGROUND] = -1  # should never be encountered
            mergetable[2] = 1  # labels started here
            dest_label = 2
            for src_label in range(3, next_label):
                # labels that map to themselves are new regions
                if mergetable[src_label] == src_label:
                    mergetable[src_label] = dest_label
                    dest_label += 1
                else:
                    # we've compacted every label below this, and the mergetable has an
                    # invariant (from mark_for_merge()) that it always points downward.
                    # Therefore, we can fetch the final lable by two steps of
                    # indirection.
                    mergetable[src_label] = mergetable[mergetable[src_label]]

            PyArray_ITER_RESET(ito)
            while PyArray_ITER_NOTDONE(ito):
                remap_line(PyArray_ITER_DATA(ito), so, L, mergetable)
                PyArray_ITER_NEXT(ito)
    except:
        # clean up and re-raise
        PyDataMem_FREE(<void *> mergetable)
        raise

    return dest_label - 1
