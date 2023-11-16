######################################################################
# Cython version of scipy.ndimage.measurements.label().
# Requires Cython version 0.17 or greater due to type templating.
######################################################################

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


######################################################################
# Load a line from a fused data array, setting the line to FOREGROUND wherever
# the fused data is nonzero, BACKGROUND elsewhere
######################################################################
cdef void fused_nonzero_line(data_t *p, np.intp_t stride,
                             np.uintp_t *line, np.intp_t L) noexcept nogil:
    cdef np.intp_t i
    for i in range(L):
        line[i] = FOREGROUND if \
            (<data_t *> ((<char *> p) + i * stride))[0] \
            else BACKGROUND


######################################################################
# Load a line from a fused data array to a np.uintp_t array
######################################################################
cdef void fused_read_line(data_t *p, np.intp_t stride,
                          np.uintp_t *line, np.intp_t L) noexcept nogil:
    cdef np.intp_t i
    for i in range(L):
        line[i] = <np.uintp_t> (<data_t *> ((<char *> p) + i * stride))[0]


######################################################################
# Store a line from a np.uintp_t array to a fused data array if possible,
# returning True if overflowed
######################################################################
cdef bint fused_write_line(data_t *p, np.intp_t stride,
                           np.uintp_t *line, np.intp_t L) noexcept nogil:
    cdef np.intp_t i
    for i in range(L):
        # Check before overwrite, as this prevents us accidentally writing a 0
        # in the foreground, which allows us to retry even when operating
        # in-place.
        if line[i] != <np.uintp_t> <data_t> line[i]:
            return True
        (<data_t *> ((<char *> p) + i * stride))[0] = <data_t> line[i]
    return False


######################################################################
# Function specializers
######################################################################
def get_nonzero_line(np.ndarray[data_t] a):
    return <Py_intptr_t> fused_nonzero_line[data_t]

def get_read_line(np.ndarray[data_t] a):
    return <Py_intptr_t> fused_read_line[data_t]

def get_write_line(np.ndarray[data_t] a):
    return <Py_intptr_t> fused_write_line[data_t]


######################################################################
# Typedefs for referring to specialized instances of fused functions
######################################################################
ctypedef void (*nonzero_line_func_t)(void *p, np.intp_t stride,
                                     np.uintp_t *line, np.intp_t L) noexcept nogil
ctypedef void (*read_line_func_t)(void *p, np.intp_t stride,
                                  np.uintp_t *line, np.intp_t L) noexcept nogil
ctypedef bint (*write_line_func_t)(void *p, np.intp_t stride,
                                   np.uintp_t *line, np.intp_t L) noexcept nogil


######################################################################
# Mark two labels to be merged
######################################################################
cdef inline np.uintp_t mark_for_merge(np.uintp_t a,
                                      np.uintp_t b,
                                      np.uintp_t *mergetable) noexcept nogil:

    cdef:
        np.uintp_t orig_a, orig_b, minlabel

    orig_a = a
    orig_b = b
    # find smallest root for each of a and b
    while a != mergetable[a]:
        a = mergetable[a]
    while b != mergetable[b]:
        b = mergetable[b]
    minlabel = a if (a < b) else b

    # merge roots
    mergetable[a] = mergetable[b] = minlabel

    # merge every step to minlabel
    a = orig_a
    b = orig_b
    while a != minlabel:
        a, mergetable[a] = mergetable[a], minlabel
    while b != minlabel:
        b, mergetable[b] = mergetable[b], minlabel

    return minlabel


######################################################################
# Take the label of a neighbor, or mark them for merging
######################################################################
cdef inline np.uintp_t take_label_or_merge(np.uintp_t cur_label,
                                           np.uintp_t neighbor_label,
                                           np.uintp_t *mergetable) noexcept nogil:
    if neighbor_label == BACKGROUND:
        return cur_label
    if cur_label == FOREGROUND:
        return neighbor_label  # neighbor is not BACKGROUND
    if neighbor_label:
        if cur_label != neighbor_label:
            cur_label = mark_for_merge(neighbor_label, cur_label, mergetable)
    return cur_label


######################################################################
# Label one line of input, using a neighbor line that has already been labeled.
######################################################################
cdef np.uintp_t label_line_with_neighbor(np.uintp_t *line,
                                         np.uintp_t *neighbor,
                                         int neighbor_use_previous,
                                         int neighbor_use_adjacent,
                                         int neighbor_use_next,
                                         np.intp_t L,
                                         bint label_unlabeled,
                                         bint use_previous,
                                         np.uintp_t next_region,
                                         np.uintp_t *mergetable) noexcept nogil:
    cdef:
        np.intp_t i

    for i in range(L):
        if line[i] != BACKGROUND:
            # See allocation of line_buffer for why this is valid when i = 0
            if neighbor_use_previous:
                line[i] = take_label_or_merge(line[i], neighbor[i - 1], mergetable)
            if neighbor_use_adjacent:
                line[i] = take_label_or_merge(line[i], neighbor[i], mergetable)
            if neighbor_use_next:
                line[i] = take_label_or_merge(line[i], neighbor[i + 1], mergetable)
            if label_unlabeled:
                if use_previous:
                    line[i] = take_label_or_merge(line[i], line[i - 1], mergetable)
                if line[i] == FOREGROUND:  # still needs a label
                    line[i] = next_region
                    mergetable[next_region] = next_region
                    next_region += 1
    return next_region

######################################################################
# Label regions
######################################################################
cpdef _label(np.ndarray input,
             np.ndarray structure,
             np.ndarray output) noexcept:
    # check dimensions
    # To understand the need for the casts to object, see
    # http://trac.cython.org/cython_trac/ticket/302
    assert (<object> input).shape == (<object> output).shape, \
        ("Shapes must match for input and output,"
         "{} != {}".format((<object> input).shape, (<object> output).shape))

    structure = np.asanyarray(structure, dtype=np.bool_).copy()
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
    if input.dtype == np.bool_:
        input = input.view(dtype=np.uint8)
    if output.dtype == np.bool_:
        # XXX - trigger special check for bit depth?
        output = output.view(dtype=np.uint8)

    cdef:
        nonzero_line_func_t nonzero_line = \
            <nonzero_line_func_t> <void *> <Py_intptr_t> get_nonzero_line(input.take([0]))
        read_line_func_t read_line = \
            <read_line_func_t> <void *> <Py_intptr_t> get_read_line(output.take([0]))
        write_line_func_t write_line = \
            <write_line_func_t> <void *> <Py_intptr_t> get_write_line(output.take([0]))
        np.flatiter _iti, _ito, _itstruct
        PyArrayIterObject *iti
        PyArrayIterObject *ito
        PyArrayIterObject *itstruct
        int axis, idim, num_neighbors, ni
        np.intp_t L, delta, i
        np.intp_t si, so, ss
        np.intp_t total_offset
        np.intp_t output_ndim, structure_ndim
        bint needs_self_labeling, valid, use_previous, overflowed
        np.ndarray _line_buffer, _neighbor_buffer
        np.uintp_t *line_buffer
        np.uintp_t *neighbor_buffer
        np.uintp_t *tmp
        np.uintp_t next_region, src_label, dest_label
        np.uintp_t mergetable_size
        np.uintp_t *mergetable

    axis = -1  # choose best axis based on output
    _ito = np.PyArray_IterAllButAxis(output, &axis)
    _iti = np.PyArray_IterAllButAxis(input, &axis)
    _itstruct = np.PyArray_IterAllButAxis(structure, &axis)

    ito = <PyArrayIterObject *> _ito
    iti = <PyArrayIterObject *> _iti
    itstruct = <PyArrayIterObject *> _itstruct

    # we only process this many neighbors from the itstruct iterator before
    # reaching the center, where we stop
    num_neighbors = structure.size // (3 * 2)

    # Create two buffer arrays for reading/writing labels.
    # Add an entry at the end and beginning to simplify some bounds checks.
    L = input.shape[axis]
    _line_buffer = np.empty(L + 2, dtype=np.uintp)
    _neighbor_buffer = np.empty(L + 2, dtype=np.uintp)
    line_buffer = <np.uintp_t *> _line_buffer.data
    neighbor_buffer = <np.uintp_t *> _neighbor_buffer.data

    # Add fenceposts with background values
    line_buffer[0] = neighbor_buffer[0] = BACKGROUND
    line_buffer[L + 1] = neighbor_buffer[L + 1] = BACKGROUND
    line_buffer = line_buffer + 1
    neighbor_buffer = neighbor_buffer + 1

    mergetable_size = 2 * output.shape[axis]
    mergetable = <np.uintp_t *> PyDataMem_NEW(mergetable_size * sizeof(np.uintp_t))
    if mergetable == NULL:
        raise MemoryError()

    try:
        # strides
        si = input.strides[axis]
        so = output.strides[axis]
        ss = structure.strides[axis]

        # 0 = background
        # 1 = foreground, needs label
        # 2... = working labels, will be compacted on output
        next_region = 2

        structure_ndim = structure.ndim
        temp = [1] * structure_ndim
        temp[axis] = 0
        use_previous = (structure[tuple(temp)] != 0)
        output_ndim = output.ndim
        output_shape = output.shape
        output_strides = output.strides

        with nogil:
            while PyArray_ITER_NOTDONE(iti):
                # Optimization - for 2D, line_buffer becomes next iteration's
                # neighbor buffer
                if output_ndim == 2:
                    tmp = line_buffer
                    line_buffer = neighbor_buffer
                    neighbor_buffer = tmp

                # copy nonzero values in input to line_buffer as FOREGROUND
                nonzero_line(PyArray_ITER_DATA(iti), si, line_buffer, L)

                # Used for labeling single lines
                needs_self_labeling = True

                # Take neighbor labels
                PyArray_ITER_RESET(itstruct)
                for ni in range(num_neighbors):
                    neighbor_use_prev = (<np.npy_bool *> PyArray_ITER_DATA(itstruct))[0]
                    neighbor_use_adjacent = (<np.npy_bool *> (<char *> PyArray_ITER_DATA(itstruct) + ss))[0]
                    neighbor_use_next = (<np.npy_bool *> (<char *> PyArray_ITER_DATA(itstruct) + 2 * ss))[0]
                    if not (neighbor_use_prev or
                            neighbor_use_adjacent or
                            neighbor_use_next):
                        PyArray_ITER_NEXT(itstruct)
                        continue

                    # Check that the neighbor line is in bounds
                    valid = True
                    total_offset = 0
                    for idim in range(structure_ndim):
                        if idim == axis:
                            continue
                        delta = (itstruct.coordinates[idim] - 1)  # 1,1,1... is center
                        if not (0 <= (ito.coordinates[idim] + delta) < output_shape[idim]):
                            valid = False
                            break
                        total_offset += delta * output_strides[idim]

                    if valid:
                        # Optimization (see above) - for 2D, line_buffer
                        # becomes next iteration's neighbor buffer, so no
                        # need to read it here.
                        if output_ndim != 2:
                            read_line(<char *> PyArray_ITER_DATA(ito) + total_offset, so,
                                      neighbor_buffer, L)

                        # be conservative about how much space we may need
                        while mergetable_size < (next_region + L):
                            mergetable_size *= 2
                            mergetable = <np.uintp_t *> \
                                PyDataMem_RENEW(<void *> mergetable,
                                                 mergetable_size * sizeof(np.uintp_t))

                        next_region = label_line_with_neighbor(line_buffer,
                                                              neighbor_buffer,
                                                              neighbor_use_prev,
                                                              neighbor_use_adjacent,
                                                              neighbor_use_next,
                                                              L,
                                                              ni == (num_neighbors - 1),
                                                              use_previous,
                                                              next_region,
                                                              mergetable)
                        if ni == (num_neighbors - 1):
                            needs_self_labeling = False
                    PyArray_ITER_NEXT(itstruct)

                if needs_self_labeling:
                    # We didn't call label_line_with_neighbor above with
                    # label_unlabeled=True, so call it now in such a way as to
                    # cause unlabeled regions to get a label.
                    while mergetable_size < (next_region + L):
                            mergetable_size *= 2
                            mergetable = <np.uintp_t *> \
                                PyDataMem_RENEW(<void *> mergetable,
                                                 mergetable_size * sizeof(np.uintp_t))

                    next_region = label_line_with_neighbor(line_buffer,
                                                          neighbor_buffer,
                                                          False, False, False,  # no neighbors
                                                          L,
                                                          True,
                                                          use_previous,
                                                          next_region,
                                                          mergetable)

                overflowed = write_line(PyArray_ITER_DATA(ito), so,
                                        line_buffer, L)
                if overflowed:
                    with gil:
                        raise NeedMoreBits()

                PyArray_ITER_NEXT(iti)
                PyArray_ITER_NEXT(ito)

            # compact the mergetable, mapping each value to its destination
            mergetable[BACKGROUND] = BACKGROUND
            mergetable[FOREGROUND] = -1  # should never be encountered
            mergetable[2] = 1  # labels started here
            dest_label = 2
            # next_region is still the original value -> we found no regions
            # set dest_label to 1 so we return 0
            if next_region < 3:
                dest_label = 1
            for src_label in range(3, next_region):
                # labels that map to themselves are new regions
                if mergetable[src_label] == src_label:
                    mergetable[src_label] = dest_label
                    dest_label += 1
                else:
                    # we've compacted every label below this, and the
                    # mergetable has an invariant (from mark_for_merge()) that
                    # it always points downward.  Therefore, we can fetch the
                    # final label by two steps of indirection.
                    mergetable[src_label] = mergetable[mergetable[src_label]]

            PyArray_ITER_RESET(ito)
            while PyArray_ITER_NOTDONE(ito):
                read_line(PyArray_ITER_DATA(ito), so, line_buffer, L)
                for i in range(L):
                    line_buffer[i] = mergetable[line_buffer[i]]
                write_line(PyArray_ITER_DATA(ito), so, line_buffer, L)
                PyArray_ITER_NEXT(ito)
    except:  # noqa: E722
        # clean up and re-raise
        PyDataMem_FREE(<void *> mergetable)
        raise

    PyDataMem_FREE(<void *> mergetable)
    return dest_label - 1
