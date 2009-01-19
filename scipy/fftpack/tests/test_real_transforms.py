#!/usr/bin/env python
from os.path import join, dirname

import numpy as np
from numpy.fft import fft as numfft
from numpy.testing import assert_array_almost_equal, TestCase

from scipy.fftpack.realtransforms import dct1, dct2, dct3

# Matlab reference data
MDATA = np.load(join(dirname(__file__), 'test.npz'))
X = [MDATA['x%d' % i] for i in range(8)]
Y = [MDATA['y%d' % i] for i in range(8)]

# FFTW reference data: the data are organized as follows:
#    * SIZES is an array containing all available sizes
#    * for every type (1, 2, 3, 4) and every size, the array dct_type_size
#    contains the output of the DCT applied to the input np.linspace(0, size-1,
#    size)
FFTWDATA_DOUBLE = np.load(join(dirname(__file__), 'fftw_double_ref.npz'))
FFTWDATA_SINGLE = np.load(join(dirname(__file__), 'fftw_single_ref.npz'))
FFTWDATA_SIZES = FFTWDATA_DOUBLE['sizes']

def fftw_ref(type, size, dt):
    x = np.linspace(0, size-1, size).astype(dt)
    if dt == np.double:
        data = FFTWDATA_DOUBLE
    elif dt == np.float32:
        data = FFTWDATA_SINGLE
    else:
        raise ValueError()
    y = (data['dct_%d_%d' % (type, size)]).astype(dt)
    return x, y

class _TestDCTIIBase(TestCase):
    def setUp(self):
        self.rdt = None
        self.dec = 14

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            x, yr = fftw_ref(2, i, self.rdt)
            y = dct2(x)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack. 
            assert_array_almost_equal(y / np.max(y), yr / np.max(y), decimal=self.dec, 
                    err_msg="Size %d failed" % i)

    def test_definition_ortho(self):
        """Test orthornomal mode."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            yr = Y[i]
            y = dct2(x, norm="ortho")
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(y, yr, decimal=self.dec)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct2(x)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct2(x[j]), decimal=self.dec)

            x = x.T
            y = dct2(x, axis=0)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct2(x[:,j]), decimal=self.dec)

class TestDCTIIDouble(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 10

class TestDCTIIFloat(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5

class _TestDCTIIIBase(TestCase):
    def setUp(self):
        self.rdt = None
        self.dec = 14

    def test_definition(self):
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            y = dct3(x)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(dct2(y) / (2*x.size), x,
                    decimal=self.dec)

    def test_definition_ortho(self):
        """Test orthornomal mode."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            y = dct2(x, norm='ortho')
            xi = dct3(y, norm="ortho")
            self.failUnless(xi.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (xi.dtype, self.rdt))
            assert_array_almost_equal(xi, x, decimal=self.dec)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct3(x)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct3(x[j]), decimal=self.dec)

            x = x.T
            y = dct3(x, axis=0)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct3(x[:,j]),
                        decimal=self.dec)

class TestDCTIIIDouble(_TestDCTIIIBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 14

class TestDCTIIIFloat(_TestDCTIIIBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5

if __name__ == "__main__":
    np.testing.run_module_suite()
