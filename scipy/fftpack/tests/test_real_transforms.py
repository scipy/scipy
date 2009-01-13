#!/usr/bin/env python
from os.path import join, dirname

import numpy as np
from numpy.fft import fft as numfft
from numpy.testing import assert_array_almost_equal

from scipy.io import loadmat

TDATA = loadmat(join(dirname(__file__), 'test.mat'), 
                squeeze_me=True,  struct_as_record=True)
X = [TDATA['x%d' % i] for i in range(8)]
Y = [TDATA['y%d' % i] for i in range(8)]

def direct_dct(x):
    """Compute a Discrete Cosine Transform, type II.

    The DCT type II is defined as:

        \forall u \in 0...N-1, 
        dct(u) = a(u) sum_{i=0}^{N-1}{f(i)cos((i + 0.5)\pi u}

    Where a(0) = sqrt(1/(4N)), a(u) = sqrt(1/(2N)) for u > 0
    """
    x = np.asarray(x)
    if not np.isrealobj(x):
        raise ValueError("Complex input not supported")
    n = x.size
    y = np.zeros(n * 4, x.dtype)
    y[1:2*n:2] = x
    y[2*n+1::2] = x[-1::-1]
    y = np.real(numfft(y))[:n]
    y[0] *= np.sqrt(.25 / n)
    y[1:] *= np.sqrt(.5 / n)
    return y

def test_ref():
    for i in range(len(X)):
        assert_array_almost_equal(direct_dct(X[i]), Y[i])

if __name__ == "__main__":
    np.testing.run_module_suite()
