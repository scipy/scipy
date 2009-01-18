#!/usr/bin/env python
from os.path import join, dirname

import numpy as np
from numpy.fft import fft as numfft
from numpy.testing import assert_array_almost_equal, TestCase

from scipy.io import loadmat
from scipy.fftpack.realtransforms import dct1, dct2, dct3

TDATA = loadmat(join(dirname(__file__), 'test.mat'),
                squeeze_me=True,  struct_as_record=True, mat_dtype=True)
X = [TDATA['x%d' % i] for i in range(8)]
Y = [TDATA['y%d' % i] for i in range(8)]

def direct_fft_dct2(x, matlab=False):
    """Compute a Discrete Cosine Transform, type II.

    The DCT type II is defined as (matlab=False):

        \forall u \in 0...N-1,
        dct(u) = 2 * sum_{i=0}^{N-1}{f(i)cos((i + 0.5)\pi u/N}

    Or (matlab=True)

        \forall u \in 0...N-1,
        dct(u) = a(u) sum_{i=0}^{N-1}{f(i)cos((i + 0.5)\pi u/N}

    Where a(0) = sqrt(1/N), a(u) = sqrt(2/N) for u > 0

    This is the exact same formula as Matlab.
    """
    x = np.asarray(x)
    if not np.isrealobj(x):
        raise ValueError("Complex input not supported")
    n = x.size
    y = np.zeros(n * 4, x.dtype)
    y[1:2*n:2] = x
    y[2*n+1::2] = x[-1::-1]
    y = np.real(numfft(y))[:n]
    if matlab:
        y[0] *= np.sqrt(.25 / n)
        y[1:] *= np.sqrt(.5 / n)
    return y

def direct_dct2(x):
    """Direct implementation (O(n^2)) of dct II.

    dct(u) = 2 * sum_{i=0}^{N-1}{f(i)cos((i + 0.5)\pi u/N}

    Note that it is not 'normalized'
    """
    n = x.size
    #a = np.empty((n, n), dtype = x.dtype)
    #for i in xrange(n):
    #    for j in xrange(n):
    #        a[i, j] = x[j] * np.cos(np.pi * (0.5 + j) * i / n)
    grd = np.outer(np.linspace(0, n - 1, n),  np.linspace(0.5, 0.5 + n - 1, n))
    a = np.cos(np.pi / n * grd) * x

    return 2 * a.sum(axis = 1)

def fdct2(x):
    """Compute a 'Fast' Discrete Cosine Transform, type II, using a N point fft
    instead of a direct 4n point DFT

        \forall u \in 0...N-1,
        dct(u) = sum_{i=0}^{N-1}{f(i)cos((i + 0.5)\pi u/N}

    See 'A Fast Cosine Transform in One and Two Dimensions', by J. Makhoul, in
    IEEE Transactions on acoustics, speech and signal processing.

    Note that it is not 'normalized'
    """
    x = np.asarray(x)
    n = x.size
    v = np.empty(x.size, x.dtype)
    if (n/2) * 2  == n:
        iseven = True
    else:
        iseven = False
    cut = (n-1)/2 + 1
    v[:cut] = x[::2]
    if iseven:
        v[cut:] = x[-1:0:-2]
    else:
        v[cut:] = x[-2::-2]
    t = 2 *  numfft(v) *  np.exp(-1j * np.pi * 0.5 / n * np.linspace(0, n-1, n))
    v[:n/2+1] = np.real(t)[:n/2+1]
    if iseven:
        v[n/2+1:] = -np.imag(t)[n/2-1:0:-1]
    else:
        v[n/2+1:] = -np.imag(t)[n/2:0:-1]
    return v

def test_refs():
    for i in range(len(X)):
        assert_array_almost_equal(direct_fft_dct2(X[i], matlab=True), Y[i])
        assert_array_almost_equal(direct_fft_dct2(X[i], matlab=False), direct_dct2(X[i]))

    for i in range(len(X)):
        x = X[i]
        y = direct_fft_dct2(x, matlab=True)
        y[0] *= np.sqrt(x.size*4)
        y[1:] *= np.sqrt(x.size*2)
        assert_array_almost_equal(y, direct_dct2(x))

def test_fdct2():
    for i in range(len(X)):
        assert_array_almost_equal(direct_dct2(X[i]), fdct2(X[i]))

class _TestDCTIIBase(TestCase):
    def setUp(self):
        self.rdt = None

    def test_definition(self):
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            yr = direct_dct2(x)
            y = dct2(x)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(y, yr)

    def test_definition_ortho(self):
        """Test orthornomal mode."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            yr = direct_fft_dct2(x, matlab=True)
            y = dct2(x, norm="ortho")
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(y, yr)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct2(x)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct2(x[j]))

            x = x.T
            y = dct2(x, axis=0)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct2(x[:,j]))

class TestDCTIIDouble(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.double

class TestDCTIIFloat(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.double

class _TestDCTIIIBase(TestCase):
    def setUp(self):
        self.rdt = None

    def test_definition(self):
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            y = dct3(x)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(dct2(y) / (2*x.size), x)

    def test_definition_ortho(self):
        """Test orthornomal mode."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            y = dct2(x, norm='ortho')
            xi = dct3(y, norm="ortho")
            self.failUnless(xi.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (xi.dtype, self.rdt))
            assert_array_almost_equal(xi, x)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct3(x)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct3(x[j]))

            x = x.T
            y = dct3(x, axis=0)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct3(x[:,j]))

class TestDCTIIIDouble(_TestDCTIIIBase):
    def setUp(self):
        self.rdt = np.double

if __name__ == "__main__":
    np.testing.run_module_suite()
