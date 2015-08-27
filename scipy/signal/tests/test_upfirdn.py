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
from numpy.testing import assert_allclose

from scipy.signal._upfirdn import upfirdn, _output_len


def upfirdn_naive(x, h, up=1, down=1):
    """Naive upfirdn processing in Python"""
    out = np.zeros(len(x) * up)
    out[::up] = x
    out = np.convolve(h, out)[::down][:_output_len(len(x), up, down, len(h))]
    return out


# Try some trivial cases first
random_state = np.random.RandomState(17)


class UpFIRDownCase(object):
    """Test UpFIRDown object"""
    def __init__(self, up, down, h):
        self.up = up
        self.down = down
        self.h = h

    def __call__(self):
        self.scrub(np.ones(1))  # tiny signal
        self.scrub(np.ones(100))  # ones
        self.scrub(random_state.randn(100))  # randn
        self.scrub(np.arange(100))  # ramp
        x = random_state.randn(25, 50)
        self.scrub(x, axis=0)
        self.scrub(x, axis=1)

    def scrub(self, x, axis=-1):
        yr = np.apply_along_axis(upfirdn_naive, axis, x,
                                 self.h, self.up, self.down)
        y = upfirdn(x, self.h, self.up, self.down, axis=axis)
        assert_allclose(yr, y)


def test_upfirdn():
    # some trivial cases
    yield UpFIRDownCase(1, 1, [1.])
    yield UpFIRDownCase(3, 2, [1.])
    yield UpFIRDownCase(2, 3, [1.])
    # mixture of big, small, and both directions (net up and net down)
    for p_max, q_max in ((100, 100), (10, 100), (100, 10), (10, 10)):
        for i in range(10):
            p_add = q_max if p_max > q_max else 1
            q_add = p_max if q_max > p_max else 1
            p = random_state.randint(p_max) + p_add
            q = random_state.randint(q_max) + q_add
            len_h = random_state.randint(200) + 1
            h = np.atleast_1d(random_state.randint(len_h))  # can be singleton
            yield UpFIRDownCase(p, q, h)


if __name__ == '__main__':
    # Execute the test suite
    for t in test_upfirdn():
        t()
