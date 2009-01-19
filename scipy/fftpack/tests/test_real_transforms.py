#!/usr/bin/env python
from os.path import join, dirname

import numpy as np
from numpy.fft import fft as numfft
from numpy.testing import assert_array_almost_equal, TestCase

from scipy.fftpack.realtransforms import dct, idct

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

class _TestDCTBase(TestCase):
    def setUp(self):
        self.rdt = None
        self.dec = 14
        self.type = None

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            x, yr = fftw_ref(self.type, i, self.rdt)
            y = dct(x, type=self.type)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack. 
            assert_array_almost_equal(y / np.max(y), yr / np.max(y), decimal=self.dec, 
                    err_msg="Size %d failed" % i)

    def test_axis(self):
        nt = 2
        for i in [7, 8, 9, 16, 32, 64]:
            x = np.random.randn(nt, i)
            y = dct(x, type=self.type)
            for j in range(nt):
                assert_array_almost_equal(y[j], dct(x[j], type=self.type),
                        decimal=self.dec)

            x = x.T
            y = dct(x, axis=0, type=self.type)
            for j in range(nt):
                assert_array_almost_equal(y[:,j], dct(x[:,j], type=self.type),
                        decimal=self.dec)

class _TestDCTIIBase(_TestDCTBase):
    def test_definition_matlab(self):
        """Test correspondance with matlab (orthornomal mode)."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            yr = Y[i]
            y = dct(x, norm="ortho", type=2)
            self.failUnless(y.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (y.dtype, self.rdt))
            assert_array_almost_equal(y, yr, decimal=self.dec)

class _TestDCTIIIBase(_TestDCTBase):
    def test_definition_ortho(self):
        """Test orthornomal mode."""
        for i in range(len(X)):
            x = np.array(X[i], dtype=self.rdt)
            y = dct(x, norm='ortho', type=2)
            xi = dct(y, norm="ortho", type=3)
            self.failUnless(xi.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (xi.dtype, self.rdt))
            assert_array_almost_equal(xi, x, decimal=self.dec)

class TestDCTIDouble(_TestDCTBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 1

class TestDCTIFloat(_TestDCTBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 1

class TestDCTIIDouble(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 2

class TestDCTIIFloat(_TestDCTIIBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 2

class TestDCTIIIDouble(_TestDCTIIIBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3

class TestDCTIIIFloat(_TestDCTIIIBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 3

class _TestIDCTBase(TestCase):
    def setUp(self):
        self.rdt = None
        self.dec = 14
        self.type = None

    def test_definition(self):
        for i in FFTWDATA_SIZES:
            xr, yr = fftw_ref(self.type, i, self.rdt)
            y = dct(xr, type=self.type)
            x = idct(yr, type=self.type)
            if self.type == 1:
                x /= 2 * (i-1)
            else:
                x /= 2 * i
            self.failUnless(x.dtype == self.rdt,
                    "Output dtype is %s, expected %s" % (x.dtype, self.rdt))
            # XXX: we divide by np.max(y) because the tests fail otherwise. We
            # should really use something like assert_array_approx_equal. The
            # difference is due to fftw using a better algorithm w.r.t error
            # propagation compared to the ones from fftpack. 
            assert_array_almost_equal(x / np.max(x), xr / np.max(x), decimal=self.dec, 
                    err_msg="Size %d failed" % i)

class TestIDCTIDouble(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 1

class TestIDCTIFloat(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 4
        self.type = 1

class TestIDCTIIDouble(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 10
        self.type = 2

class TestIDCTIIFloat(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 2

class TestIDCTIIIDouble(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.double
        self.dec = 14
        self.type = 3

class TestIDCTIIIFloat(_TestIDCTBase):
    def setUp(self):
        self.rdt = np.float32
        self.dec = 5
        self.type = 3

if __name__ == "__main__":
    np.testing.run_module_suite()
