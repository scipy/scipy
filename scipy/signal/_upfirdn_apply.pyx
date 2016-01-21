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


ctypedef double complex double_complex
ctypedef float complex float_complex

ctypedef fused DTYPE_t:
    # Eventually we could add "object", too, but then we'd lose the "nogil"
    # on the _apply_impl function.
    float
    float_complex
    double
    double_complex


def _output_len(Py_ssize_t len_h,
                Py_ssize_t in_len,
                Py_ssize_t up,
                Py_ssize_t down):
    """The output length that results from a given input"""
    cdef Py_ssize_t np
    cdef Py_ssize_t in_len_copy
    in_len_copy = in_len + (len_h + (-len_h % up)) // up - 1
    np = in_len_copy * up
    cdef Py_ssize_t need = np // down
    if np % down > 0:
        need += 1
    return need


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
def _rat(x, tol=None):
    """Calculate the continued fraction representation of x"""
    # Adapted from https://en.wikipedia.org/wiki/Continued_fraction
    x = float(x)
    if not np.isfinite(x):
        raise RuntimeError('x must be finite')
    if tol is None:
        tol = np.abs(x) * 1e-12
    else:
        tol = float(tol)
    cdef double use_tol = tol
    cdef double sign = np.sign(x)
    cdef double use_x = x * sign
    cdef double next_ = use_x
    # compute the convergents (ck) for each step
    cdef int n = 0
    cdef int hks[2][2]
    cdef int next_hk[2]
    cdef double frac_part
    cdef int integer_part
    cdef double min_frac = 1. / 2147483647

    n = 0
    while True:
        integer_part = int(next_)
        if n > 1:
            next_hk[0] = integer_part * hks[0][0] + hks[1][0]
            next_hk[1] = integer_part * hks[0][1] + hks[1][1]
        elif n == 1:
            next_hk[0] = integer_part * hks[0][0] + 1
            next_hk[1] = integer_part
        elif n == 0:
            next_hk[0] = integer_part
            next_hk[1] = 1
        hks[1][0] = hks[0][0]
        hks[1][1] = hks[0][1]
        hks[0][0] = next_hk[0]
        hks[0][1] = next_hk[1]
        frac_part = next_ - integer_part
        if (frac_part == 0 or
            abs(hks[0][0] / float(hks[0][1]) - use_x) < use_tol or
            frac_part <= min_frac):
            break
        if n > 10000:
            break
        next_ = 1. / frac_part
        n += 1
    if sign < 0:
        hks[0][0] *= -1
    return np.array(hks[0])


def _apply(DTYPE_t [:] x, DTYPE_t [:] h_trans_flip, DTYPE_t [:] out,
                 Py_ssize_t up, Py_ssize_t down):
    _apply_impl(x, h_trans_flip, out, up, down)


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef void _apply_impl(DTYPE_t [:] x, DTYPE_t [:] h_trans_flip, DTYPE_t [:] out,
                 Py_ssize_t up, Py_ssize_t down) nogil:
    cdef Py_ssize_t len_x = x.shape[0]
    cdef Py_ssize_t h_per_phase = h_trans_flip.shape[0] / up
    cdef Py_ssize_t padded_len = len_x + h_per_phase - 1
    cdef Py_ssize_t x_idx = 0
    cdef Py_ssize_t y_idx = 0
    cdef Py_ssize_t h_idx = 0
    cdef Py_ssize_t t = 0
    cdef Py_ssize_t x_conv_idx = 0

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
            if x_conv_idx < len_x and x_conv_idx > 0:
                out[y_idx] = out[y_idx] + x[x_conv_idx] * h_trans_flip[h_idx]
            h_idx += 1
        y_idx += 1
        t += down
        x_idx += t / up  # integer div
        t = t % up
