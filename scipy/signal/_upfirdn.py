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

import numpy as np

from ._upfirdn_apply import _output_len, _apply


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



def upfirdn(x, h, up=1, down=1, axis=-1):
    """Upsample, FIR filter, and downsample

    Parameters
    ----------
    x : array-like
        Input signal array.
    h : array-like
        1-dimensional FIR (finite-impulse response) filter coefficients.
    up : int
        Upsampling rate.
    down : int
        Downsampling rate.
    axis : int
        The axis of ``x`` to operate along.

    Returns
    -------
    y : ndarray
        The output signal array. Dimensions will be the same as ``x`` except
        for along ``axis``, which will change size according to the ``h``,
        ``up``,  and ``down`` parameters.

    Notes
    -----
    The algorithm is an implementation of the block diagram shown on page 129
    of the Vaidyanathan text [1]_ (Figure 4.3-8d).

    .. [1] P. P. Vaidyanathan, Multirate Systems and Filter Banks,
       Prentice Hall, 1993.

    The direct approach of upsampling by factor of P with zero insertion,
    FIR filtering of length ``N``, and downsampling by factor of Q is
    O(N*Q) per output sample. The polyphase implementation used here is
    O(N/P).

    .. versionadded:: 0.17

    Examples
    --------
    Simple operations:

    >>> from scipy.signal import upfirdn
    >>> upfirdn([1,1,1], [1,1,1])   # FIR filter
    array([ 1.,  2.,  3.,  2.,  1.])
    >>> upfirdn([1, 2, 3], [1], 3)  # upsampling with zeros insertion
    array([ 1.,  0.,  0.,  2.,  0.,  0.,  3.,  0.,  0.])
    >>> upfirdn([1,2,3], [1,1,1], 3)  # upsampling with sample-and-hold
    array([ 1.,  1.,  1.,  2.,  2.,  2.,  3.,  3.,  3.])
    >>> upfirdn([1,1,1], [.5,1,.5], 2)  # linear interpolation
    array([ 0.5,  1. ,  1. ,  1. ,  1. ,  1. ,  0.5,  0. ])
    >>> upfirdn(np.arange(10), [1], 1, 3)  # decimation by 3
    array([ 0.,  3.,  6.,  9.])
    >>> upfirdn(np.arange(10), [.5,1,.5], 2, 3)  # linear interp, rate 2/3
    array([ 0. ,  1. ,  2.5,  4. ,  5.5,  7. ,  8.5,  0. ])

    Apply a single filter to multiple signals:

    >>> x = np.reshape(np.arange(8), (4,2))
    >>> x
    array([[0, 1],
           [2, 3],
           [4, 5],
           [6, 7]])

    Apply along the last dimension of ``x``:

    >>> h = [1, 1]
    >>> upfirdn(x, h, 2)
    array([[ 0.,  0.,  1.,  1.],
           [ 2.,  2.,  3.,  3.],
           [ 4.,  4.,  5.,  5.],
           [ 6.,  6.,  7.,  7.]])

    Apply along the 0th dimension of ``x``:

    >>> upfirdn(x, h, 2, axis=0)
    array([[ 0.,  1.],
           [ 0.,  1.],
           [ 2.,  3.],
           [ 2.,  3.],
           [ 4.,  5.],
           [ 4.,  5.],
           [ 6.,  7.],
           [ 6.,  7.]])

    """
    ufd = UpFIRDown(h, up, down)
    return np.apply_along_axis(ufd.apply, axis, x)
