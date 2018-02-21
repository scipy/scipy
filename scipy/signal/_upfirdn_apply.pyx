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

from __future__ import absolute_import

cimport cython
cimport numpy as np
import numpy as np
from cython import bint  # boolean integer type
from libc.stdlib cimport malloc, free


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
    cdef np.intp_t nt
    cdef np.intp_t in_len_copy
    in_len_copy = in_len + (len_h + (-len_h % up)) // up - 1
    nt = in_len_copy * up
    cdef np.intp_t need = nt // down
    if nt % down > 0:
        need += 1
    return need


def _apply(np.ndarray data, DTYPE_t [::1] h_trans_flip, np.ndarray out,
           np.intp_t up, np.intp_t down, np.intp_t axis):
    cdef ArrayInfo data_info, output_info
    cdef np.intp_t len_h = h_trans_flip.size
    cdef DTYPE_t *data_ptr
    cdef DTYPE_t *filter_ptr
    cdef DTYPE_t *out_ptr
    cdef int retval

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
                                   up, down, axis)
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
cdef int _apply_axis_inner(DTYPE_t* data, ArrayInfo data_info,
                           DTYPE_t* h_trans_flip, np.intp_t len_h,
                           DTYPE_t* output, ArrayInfo output_info,
                           np.intp_t up, np.intp_t down,
                           np.intp_t axis) nogil:
    cdef np.intp_t i
    cdef np.intp_t num_loops = 1
    cdef bint make_temp_data, make_temp_output
    cdef DTYPE_t* temp_data = NULL
    cdef DTYPE_t* temp_output = NULL

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
        temp_output = <DTYPE_t*>malloc(output_info.shape[axis] * sizeof(DTYPE_t))
        if not temp_output:
            free(temp_data)
            free(temp_output)
            return 4

    for i in range(output_info.ndim):
        if i != axis:
            num_loops *= output_info.shape[i]

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
            for j in range(data_info.shape[axis]):
                # Offsets are byte offsets, to need to cast to char and back
                tmp_ptr = <DTYPE_t *>((<char *> data) + data_offset +
                    j * data_info.strides[axis])
                temp_data[j] = tmp_ptr[0]

        # Select temporary or direct output and data
        if make_temp_data:
            data_row = temp_data
        else:
            data_row = <DTYPE_t *>((<char *>data) + data_offset)
        if make_temp_output:
            output_row = temp_output
            for j in range(output_info.shape[axis]):
                output_row[j] = 0.0  # initialize as zeros
        else:
            output_row = <DTYPE_t *>((<char *>output) + output_offset)

        # call 1D upfirdn
        _apply_impl(data_row, data_info.shape[axis],
                    h_trans_flip, len_h, output_row, up, down)

        # Copy from temporary output if necessary
        if make_temp_output:
            for j in range(output_info.shape[axis]):
                # Offsets are byte offsets, to need to cast to char and back
                tmp_ptr = <DTYPE_t *>((<char *>output) + output_offset +
                    j * output_info.strides[axis])
                tmp_ptr[0] = output_row[j]

    # cleanup
    free(temp_data)
    free(temp_output)
    return 0


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef void _apply_impl(DTYPE_t *x, np.intp_t len_x, DTYPE_t *h_trans_flip,
                      np.intp_t len_h, DTYPE_t *out,
                      np.intp_t up, np.intp_t down) nogil:
    cdef np.intp_t h_per_phase = len_h / up
    cdef np.intp_t padded_len = len_x + h_per_phase - 1
    cdef np.intp_t x_idx = 0
    cdef np.intp_t y_idx = 0
    cdef np.intp_t h_idx = 0
    cdef np.intp_t t = 0
    cdef np.intp_t x_conv_idx = 0

    while x_idx < len_x:
        h_idx = t * h_per_phase
        x_conv_idx = x_idx - h_per_phase + 1
        if x_conv_idx < 0:
            h_idx -= x_conv_idx
            x_conv_idx = 0
        for x_conv_idx in range(x_conv_idx, x_idx + 1):
            out[y_idx] = out[y_idx] + x[x_conv_idx] * h_trans_flip[h_idx]
            h_idx += 1
        # store and increment
        y_idx += 1
        t += down
        x_idx += t / up  # integer div
        # which phase of the filter to use
        t = t % up

    # Use a second simplified loop to flush out the last bits
    while x_idx < padded_len:
        h_idx = t * h_per_phase
        x_conv_idx = x_idx - h_per_phase + 1
        for x_conv_idx in range(x_conv_idx, x_idx + 1):
            if x_conv_idx < len_x and x_conv_idx >= 0:
                out[y_idx] = out[y_idx] + x[x_conv_idx] * h_trans_flip[h_idx]
            h_idx += 1
        y_idx += 1
        t += down
        x_idx += t / up  # integer div
        t = t % up
