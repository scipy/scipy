# -*- coding: utf-8 -*-

# Code adapted from "upfirdn" python library with permission:
#
# Copyright (c) 2009, Motorola, Inc
#
# All Rights Reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# * Neither the name of Motorola nor the names of its contributors may be
# used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cimport cython
cimport numpy as np
import numpy as np
from cython import bint  # boolean integer type
from libc.stdlib cimport malloc, free
from libc.string cimport memset


ctypedef double complex double_complex
ctypedef float complex float_complex

ctypedef fused DTYPE_t:
    # Eventually we could add "object", too, but then we'd lose the "nogil"
    # on the _apply_impl function.
    float
    float_complex
    double
    double_complex

cdef struct ArrayInfo:
    np.intp_t * shape
    np.intp_t * strides
    np.intp_t ndim


def _output_len(np.intp_t len_h,
                np.intp_t in_len,
                np.intp_t up,
                np.intp_t down):
    """The output length that results from a given input"""
    # ceil(((in_len - 1) * up + len_h) / down), but using integer arithmetic
    return (((in_len - 1) * up + len_h) - 1) // down + 1


# Signal extension modes
ctypedef enum MODE:
    MODE_CONSTANT = 0
    MODE_SYMMETRIC = 1       #  3 2 1 | 1 2 3 | 3 2 1
    MODE_CONSTANT_EDGE = 2
    MODE_SMOOTH = 3
    MODE_PERIODIC = 4
    MODE_REFLECT = 5         #  3 2 | 1 2 3 | 2 1
    MODE_ANTISYMMETRIC = 6
    MODE_ANTIREFLECT = 7
    MODE_LINE = 8  # slope determined by first and last entries of the array


cpdef MODE mode_enum(mode):
    if mode == 'constant':
        return MODE_CONSTANT
    elif mode == 'symmetric':
        return MODE_SYMMETRIC
    elif mode == 'edge':
        return MODE_CONSTANT_EDGE
    elif mode == 'smooth':
        return MODE_SMOOTH
    elif mode == 'wrap':
        return MODE_PERIODIC
    elif mode == 'reflect':
        return MODE_REFLECT
    elif mode == 'antisymmetric':
        return MODE_ANTISYMMETRIC
    elif mode == 'antireflect':
        return MODE_ANTIREFLECT
    elif mode == 'line':
        return MODE_LINE
    else:
        raise ValueError("Unknown mode: {}".format(mode))


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef DTYPE_t _extend_left(DTYPE_t *x, np.intp_t idx, np.intp_t len_x,
                          MODE mode, DTYPE_t cval) nogil:
    cdef DTYPE_t le = 0.
    cdef DTYPE_t lin_slope = 0.

    # note: idx will be < 0
    if mode == MODE_SYMMETRIC:
        if (-idx) < len_x:
            return x[-idx - 1]
        else:
            # general case for multiple reflections:
            # the pattern repeats with periodicity 2*len_x
            idx = (-idx - 1) % (2 * len_x)
            if idx < len_x:
                return x[idx]
            else:
                return x[len_x - 1 - (idx - len_x)]
    elif mode == MODE_REFLECT:
        if (-idx) < (len_x - 1):
            return x[-idx]
        else:
            # general case for multiple reflections:
            # the pattern repeats with periodicity 2*(len_x - 1)
            idx = (-idx - 1) % (2 * (len_x - 1))
            if idx < (len_x - 1):
                return x[idx + 1]
            else:
                return x[len_x - 2 - (idx - (len_x - 1))]
    elif mode == MODE_PERIODIC:
        idx = (-idx - 1) % len_x
        return x[len_x - idx - 1]
    elif mode == MODE_SMOOTH:
        return x[0] + idx * (x[1] - x[0])
    elif mode == MODE_LINE:
        lin_slope = (x[len_x - 1] - x[0]) / (len_x - 1)
        return x[0] + idx * lin_slope
    elif mode == MODE_ANTISYMMETRIC:
        if (-idx) < len_x:
            return -x[-idx - 1]
        else:
            idx = (-idx - 1) % (2 * len_x)
            if idx < len_x:
                return -x[idx]
            else:
                return x[len_x - 1 - (idx - len_x)]
    elif mode == MODE_ANTIREFLECT:
        if (-idx) < len_x:
            return x[0] - (x[-idx] - x[0])
        else:
            le = x[0] + (x[0] - x[len_x - 1]) * ((-(idx) - 1) // (len_x - 1))
            idx = (-idx - 1) % (2 * (len_x - 1))
            if idx < (len_x - 1):
                return le - (x[idx + 1] - x[0])
            else:
                return le - (x[len_x - 1] - x[len_x - 2 - (idx - (len_x - 1))])
    elif mode == MODE_CONSTANT_EDGE:
        return x[0]
    elif mode == MODE_CONSTANT:
        return cval
    else:
        return -1.


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef DTYPE_t _extend_right(DTYPE_t *x, np.intp_t idx, np.intp_t len_x,
                           MODE mode, DTYPE_t cval) nogil:
    # note: idx will be >= len_x
    cdef DTYPE_t re = 0.
    cdef DTYPE_t lin_slope = 0.

    if mode == MODE_SYMMETRIC:
        if idx < (2 * len_x):
            return x[len_x - 1 - (idx - len_x)]
        else:
            idx = idx % (2 * len_x)
            if idx < len_x:
                return x[idx]
            else:
                return x[len_x - 1 - (idx - len_x)]
    elif mode == MODE_REFLECT:
        if idx < (2 * len_x - 1):
            return x[len_x - 2 - (idx - len_x)]
        else:
            idx = idx % (2 * (len_x - 1))
            if idx < (len_x - 1):
                return x[idx]
            else:
                return x[len_x - 1 - (idx - (len_x - 1))]
    elif mode == MODE_PERIODIC:
        return x[idx % len_x]
    elif mode == MODE_SMOOTH:
        return x[len_x - 1] + (idx - len_x + 1) * (x[len_x - 1] - x[len_x - 2])
    elif mode == MODE_LINE:
        lin_slope = (x[len_x - 1] - x[0]) / (len_x - 1)
        return x[len_x - 1] + (idx - len_x + 1) * lin_slope
    elif mode == MODE_CONSTANT_EDGE:
        return x[len_x - 1]
    elif mode == MODE_ANTISYMMETRIC:
        if idx < (2 * len_x):
            return -x[len_x - 1 - (idx - len_x)]
        else:
            idx = idx % (2 * len_x)
            if idx < len_x:
                return x[idx]
            else:
                return -x[len_x - 1 - (idx - len_x)]
    elif mode == MODE_ANTIREFLECT:
        if idx < (2 * len_x - 1):
            return x[len_x - 1] - (x[len_x - 2 - (idx - len_x)] - x[len_x - 1])
        else:
            re = x[len_x - 1] + (x[len_x - 1] - x[0]) * (idx // (len_x - 1) - 1)
            idx = idx % (2 * (len_x - 1))
            if idx < (len_x - 1):
                return re + (x[idx] - x[0])
            else:
                return re + (x[len_x - 1] - x[len_x - 1 - (idx - (len_x - 1))])
    elif mode == MODE_CONSTANT:
        return cval
    else:
        return -1.


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef _pad_test(np.ndarray[DTYPE_t] data, np.intp_t npre=0, np.intp_t npost=0,
                object mode=0, DTYPE_t cval=0):
    """1D test function for signal extension modes.

    Returns ``data extended by ``npre``, ``npost`` at the beginning, end.
    """
    cdef np.intp_t idx
    cdef np.intp_t cnt = 0
    cdef np.intp_t len_x = data.size
    cdef np.intp_t len_out = npre + len_x + npost
    cdef DTYPE_t xval
    cdef DTYPE_t [::1] out
    cdef DTYPE_t* data_ptr
    cdef MODE _mode
    _mode = mode_enum(mode)

    if DTYPE_t is float:
        out = np.zeros((len_out,), dtype=np.float32)
    elif DTYPE_t is float_complex:
        out = np.zeros((len_out,), dtype=np.complex64)
    elif DTYPE_t is double:
        out = np.zeros((len_out,), dtype=np.float64)
    elif DTYPE_t is double_complex:
        out = np.zeros((len_out,), dtype=np.complex128)
    else:
        raise ValueError("unsupported dtype")

    data_ptr = <DTYPE_t*> data.data
    with nogil:
        for idx in range(-npre, len_x + npost, 1):
            if idx < 0:
                xval = _extend_left(data_ptr, idx, len_x, _mode, cval)
            elif idx >= len_x:
                xval = _extend_right(data_ptr, idx, len_x, _mode, cval)
            else:
                xval = data_ptr[idx]
            out[cnt] = xval
            cnt += 1
    return np.asarray(out)


def _apply(np.ndarray data, DTYPE_t [::1] h_trans_flip, np.ndarray out,
           np.intp_t up, np.intp_t down, np.intp_t axis, np.intp_t mode,
           DTYPE_t cval):
    cdef ArrayInfo data_info, output_info
    cdef np.intp_t len_h = h_trans_flip.size
    cdef DTYPE_t *data_ptr
    cdef DTYPE_t *filter_ptr
    cdef DTYPE_t *out_ptr
    cdef int retval
    cdef np.intp_t len_out = out.shape[axis]

    data_info.ndim = data.ndim
    data_info.strides = <np.intp_t *> data.strides
    data_info.shape = <np.intp_t *> data.shape

    output_info.ndim = out.ndim
    output_info.strides = <np.intp_t *> out.strides
    output_info.shape = <np.intp_t *> out.shape

    data_ptr = <DTYPE_t*> data.data
    filter_ptr = <DTYPE_t*> &h_trans_flip[0]
    out_ptr = <DTYPE_t*> out.data

    with nogil:
        retval = _apply_axis_inner(data_ptr, data_info,
                                   filter_ptr, len_h,
                                   out_ptr, output_info,
                                   up, down, axis, <MODE>mode, cval, len_out)
    if retval == 1:
        raise ValueError("failure in _apply_axis_inner: data and output arrays"
                         " must have the same number of dimensions.")
    elif retval == 2:
        raise ValueError(
            ("failure in _apply_axis_inner: axis = {}, ".format(axis) +
             "but data_info.ndim is only {}.".format(data_info.ndim)))
    elif retval == 3 or retval == 4:
        raise MemoryError()
    return out


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _apply_axis_inner(DTYPE_t* data, ArrayInfo data_info,
                           DTYPE_t* h_trans_flip, np.intp_t len_h,
                           DTYPE_t* output, ArrayInfo output_info,
                           np.intp_t up, np.intp_t down,
                           np.intp_t axis, MODE mode, DTYPE_t cval,
                           np.intp_t len_out) nogil:
    cdef np.intp_t i
    cdef np.intp_t num_loops = 1
    cdef bint make_temp_data, make_temp_output
    cdef DTYPE_t* temp_data = NULL
    cdef DTYPE_t* temp_output = NULL
    cdef size_t row_size_bytes

    if data_info.ndim != output_info.ndim:
        return 1
    if axis >= data_info.ndim:
        return 2

    make_temp_data = data_info.strides[axis] != sizeof(DTYPE_t);
    make_temp_output = output_info.strides[axis] != sizeof(DTYPE_t);
    if make_temp_data:
        temp_data = <DTYPE_t*>malloc(data_info.shape[axis] * sizeof(DTYPE_t))
        if not temp_data:
            free(temp_data)
            return 3
    if make_temp_output:
        row_size_bytes = output_info.shape[axis] * sizeof(DTYPE_t)
        temp_output = <DTYPE_t*>malloc(row_size_bytes)
        if not temp_output:
            free(temp_data)
            free(temp_output)
            return 4

    for i in range(output_info.ndim):
        if i != axis:
            num_loops *= output_info.shape[i]

    # strides in number of elements rather than number of bytes
    cdef np.intp_t idx_stride = data_info.strides[axis] / sizeof(DTYPE_t)
    cdef np.intp_t idx_stride_out = output_info.strides[axis] / sizeof(DTYPE_t)

    cdef np.intp_t j
    cdef np.intp_t data_offset
    cdef np.intp_t output_offset
    cdef DTYPE_t* data_row
    cdef DTYPE_t* output_row
    cdef np.intp_t reduced_idx
    cdef np.intp_t j_rev
    cdef np.intp_t axis_idx
    cdef DTYPE_t* tmp_ptr = NULL
    for i in range(num_loops):
        data_offset = 0
        output_offset = 0
        # Calculate offset into linear buffer
        reduced_idx = i
        for j in range(output_info.ndim):
            j_rev = output_info.ndim - 1 - j
            if j_rev != axis:
                axis_idx = reduced_idx % output_info.shape[j_rev]
                reduced_idx /= output_info.shape[j_rev]
                data_offset += (axis_idx * data_info.strides[j_rev])
                output_offset += (axis_idx * output_info.strides[j_rev])

        # Copy to temporary data if necessary
        if make_temp_data:
            # Offsets are byte offsets, so need to cast to char and back
            tmp_ptr = <DTYPE_t *>((<char *> data) + data_offset)
            for j in range(data_info.shape[axis]):
                temp_data[j] = tmp_ptr[idx_stride*j]

        # Select temporary or direct output and data
        if make_temp_data:
            data_row = temp_data
        else:
            data_row = <DTYPE_t *>((<char *>data) + data_offset)
        if make_temp_output:
            output_row = temp_output
            memset(output_row, 0, row_size_bytes)
        else:
            output_row = <DTYPE_t *>((<char *>output) + output_offset)

        # call 1D upfirdn
        _apply_impl(data_row, data_info.shape[axis],
                    h_trans_flip, len_h, output_row, up, down, mode, cval,
                    len_out)

        # Copy from temporary output if necessary
        if make_temp_output:
            tmp_ptr = <DTYPE_t *>((<char *>output) + output_offset)
            for j in range(output_info.shape[axis]):
                tmp_ptr[idx_stride_out*j] = output_row[j]

    # cleanup
    free(temp_data)
    free(temp_output)
    return 0


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef void _apply_impl(DTYPE_t *x, np.intp_t len_x, DTYPE_t *h_trans_flip,
                      np.intp_t len_h, DTYPE_t *out,
                      np.intp_t up, np.intp_t down, MODE mode,
                      DTYPE_t cval, np.intp_t len_out) nogil:
    cdef np.intp_t h_per_phase = len_h // up
    cdef np.intp_t padded_len = len_x + h_per_phase - 1
    cdef np.intp_t x_idx = 0
    cdef np.intp_t y_idx = 0
    cdef np.intp_t h_idx = 0
    cdef np.intp_t t = 0
    cdef np.intp_t x_conv_idx = 0
    cdef DTYPE_t xval
    cdef bint zpad

    zpad = (mode == MODE_CONSTANT and cval == 0)
    if len_out == 0:
        return

    while x_idx < len_x:
        h_idx = t * h_per_phase
        x_conv_idx = x_idx - h_per_phase + 1
        if x_conv_idx < 0:
            if zpad:
                h_idx -= x_conv_idx
            else:
                for x_conv_idx in range(x_conv_idx, 0):
                    xval = _extend_left(x, x_conv_idx, len_x, mode, cval)
                    out[y_idx] += xval * h_trans_flip[h_idx]
                    h_idx += 1
            x_conv_idx = 0
        for x_conv_idx in range(x_conv_idx, x_idx + 1):
            out[y_idx] = out[y_idx] + x[x_conv_idx] * h_trans_flip[h_idx]
            h_idx += 1
        # store and increment
        y_idx += 1
        if y_idx >= len_out:
            return
        t += down
        x_idx += t // up
        # which phase of the filter to use
        t = t % up

    # Use a second simplified loop to flush out the last bits
    while x_idx < padded_len:
        h_idx = t * h_per_phase
        x_conv_idx = x_idx - h_per_phase + 1
        for x_conv_idx in range(x_conv_idx, x_idx + 1):
            if x_conv_idx >= len_x:
                xval = _extend_right(x, x_conv_idx, len_x, mode, cval)
            elif x_conv_idx < 0:
                xval = _extend_left(x, x_conv_idx, len_x, mode, cval)
            else:
                xval = x[x_conv_idx]
            out[y_idx] += xval * h_trans_flip[h_idx]
            h_idx += 1
        y_idx += 1
        if y_idx >= len_out:
            return
        t += down
        x_idx += t // up
        t = t % up
