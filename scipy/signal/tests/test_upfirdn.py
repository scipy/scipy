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
from numpy.testing import assert_allclose, assert_equal, assert_raises

from scipy.signal import upfirdn
from scipy.signal._upfirdn import _output_len


def upfirdn_naive(x, h, up=1, down=1):
    """Naive upfirdn processing in Python

    Note: arg order (x, h) differs to facilitate apply_along_axis use.
    """
    h = np.asarray(h)
    out = np.zeros(len(x) * up, x.dtype)
    out[::up] = x
    out = np.convolve(h, out)[::down][:_output_len(len(h), len(x), up, down)]
    return out


class _UpFIRDnCase(object):
    """Test _UpFIRDn object"""
    def __init__(self, up, down, h, x_dtype):
        self.up = up
        self.down = down
        self.h = np.atleast_1d(h)
        self.x_dtype = x_dtype
        self.rng = np.random.RandomState(17)

    def __call__(self):
        # tiny signal
        self.scrub(np.ones(1, self.x_dtype))
        # ones
        self.scrub(np.ones(10, self.x_dtype))  # ones
        # randn
        x = self.rng.randn(10).astype(self.x_dtype)
        if self.x_dtype in (np.complex64, np.complex128):
            x += 1j * self.rng.randn(10)
        self.scrub(x)
        # ramp
        self.scrub(np.arange(10).astype(self.x_dtype))
        # 3D, random
        size = (2, 3, 5)
        x = self.rng.randn(*size).astype(self.x_dtype)
        if self.x_dtype in (np.complex64, np.complex128):
            x += 1j * self.rng.randn(*size)
        for axis in range(len(size)):
            self.scrub(x, axis=axis)
        x = x[:, ::2, 1::3].T
        for axis in range(len(size)):
            self.scrub(x, axis=axis)

    def scrub(self, x, axis=-1):
        yr = np.apply_along_axis(upfirdn_naive, axis, x,
                                 self.h, self.up, self.down)
        y = upfirdn(self.h, x, self.up, self.down, axis=axis)
        dtypes = (self.h.dtype, x.dtype)
        if all(d == np.complex64 for d in dtypes):
            assert_equal(y.dtype, np.complex64)
        elif np.complex64 in dtypes and np.float32 in dtypes:
            assert_equal(y.dtype, np.complex64)
        elif all(d == np.float32 for d in dtypes):
            assert_equal(y.dtype, np.float32)
        elif np.complex128 in dtypes or np.complex64 in dtypes:
            assert_equal(y.dtype, np.complex128)
        else:
            assert_equal(y.dtype, np.float64)
        assert_allclose(yr, y)


def test_upfirdn():
    """Test upfirdn"""
    big, small = 100, 10  # up/down factors
    longest_h = 25
    n_rep = 3
    try_types = (int, np.float32, np.complex64, float, complex)
    random_state = np.random.RandomState(17)
    for x_dtype in try_types:
        # some trivial cases
        for h in (1., 1j):
            yield _UpFIRDnCase(1, 1, h, x_dtype)
            yield _UpFIRDnCase(2, 2, h, x_dtype)
            yield _UpFIRDnCase(3, 2, h, x_dtype)
            yield _UpFIRDnCase(2, 3, h, x_dtype)
        for h_dtype in try_types:
            # mixture of big, small, and both directions (net up and net down)
            for p_max, q_max in ((big, big),
                                 (small, big),
                                 (big, small),
                                 (small, small)):
                for _ in range(n_rep):
                    p_add = q_max if p_max > q_max else 1
                    q_add = p_max if q_max > p_max else 1
                    p = random_state.randint(p_max) + p_add
                    q = random_state.randint(q_max) + q_add
                    len_h = random_state.randint(longest_h) + 1
                    h = np.atleast_1d(random_state.randint(len_h))
                    h = h.astype(h_dtype)
                    if h_dtype == complex:
                        h += 1j * random_state.randint(len_h)
                    yield _UpFIRDnCase(p, q, h, x_dtype)
    # Degenerate cases
    assert_raises(ValueError, upfirdn, [1], [1], 1, 0)  # up or down < 1
    assert_raises(ValueError, upfirdn, [], [1], 1, 1)  # h.ndim != 1
    assert_raises(ValueError, upfirdn, [[1]], [1], 1, 1)


if __name__ == '__main__':
    # Execute the test suite
    for t in test_upfirdn():
        t()
