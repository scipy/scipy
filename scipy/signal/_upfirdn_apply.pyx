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


def _pad_h(h, up):
    """Store coefficients in a transposed, flipped arrangement.

    For example, suppose upRate is 3, and the
    input number of coefficients is 10, represented as h[0], ..., h[9].

    Then the internal buffer will look like this::

       h[9], h[6], h[3], h[0],   // flipped phase 0 coefs
       0,    h[7], h[4], h[1],   // flipped phase 1 coefs (zero-padded)
       0,    h[8], h[5], h[2],   // flipped phase 2 coefs (zero-padded)

    """
    h_padlen = len(h) + (-len(h) % up)
    h_full = np.zeros(h_padlen)
    h_full[:len(h)] = h
    h_full = h_full.reshape(-1, up).T[:, ::-1].ravel()
    return h_full


class UpFIRDown(object):
    def __init__(self, h, up, down):
        """Helper for resampling"""
        h = np.asarray(h, np.float64)
        if h.ndim != 1 or h.size == 0:
            raise ValueError('h must be 1D with non-zero length')
        self._up = int(up)
        self._down = int(down)
        if self._up < 1 or self._down < 1:
            raise RuntimeError('Both up and down must be >= 1')
        # This both transposes, and "flips" each phase for filtering
        self._h_trans_flip = _pad_h(h, self._up)

    def _output_len(self, len_x):
        """Helper to get the output length given an input length"""
        return _output_len(len_x + len(self._h_trans_flip) // self._up - 1,
                           self._up, self._down, 0)

    def apply(self, x):
        """Apply the prepared filter to a 1D signal x"""
        out = np.zeros(self._output_len(len(x)), dtype=np.float64)
        _apply(np.asarray(x, np.float64), self._h_trans_flip, out,
               self._up, self._down)
        return out


@cython.cdivision(True)
def _output_len(Py_ssize_t in_len,
                Py_ssize_t up,
                Py_ssize_t down,
                Py_ssize_t len_h):
    """The output length that results from a given input"""
    cdef Py_ssize_t np
    if len_h != 0:  # figure out what the padding needs to be
        # the modulo below is the C-equivalent of a Python modulo,
        # i.e. -1 % 10 = 9
        in_len = in_len + (len_h + ((-len_h % up) + up) % up) // up - 1
    np = in_len * up
    cdef Py_ssize_t need = np // down
    if np % down > 0:
        need += 1
    return need


ctypedef fused DTYPE_t:
    np.float64_t


@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
cdef void _apply(DTYPE_t [:] x, DTYPE_t [:] h_trans_flip, DTYPE_t [:] out,
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
            out[y_idx] += x[x_conv_idx] * h_trans_flip[h_idx]
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
                out[y_idx] += x[x_conv_idx] * h_trans_flip[h_idx]
            h_idx += 1
        y_idx += 1
        t += down
        x_idx += t / up  # integer div
        t = t % up
