#!/usr/bin/env python
# Created by Pearu Peterson, September 2002
""" Test functions for fftpack.basic module
"""
__usage__ = """
Build fftpack:
  python setup_fftpack.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.fftpack.test()'
Run tests if fftpack is not installed:
  python tests/test_basic.py
"""

from numpy.testing import *
from scipy.fftpack import ifft,fft,fftn,ifftn,rfft,irfft, fft2
from scipy.fftpack import _fftpack as fftpack

from numpy import arange, add, array, asarray, zeros, dot, exp, pi,\
     swapaxes, double, cdouble
import numpy.fft

from numpy.random import rand
def random(size):
    return rand(*size)

def get_mat(n):
    data = arange(n)
    data = add.outer(data,data)
    return data

def direct_dft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n,dtype=cdouble)
    w = -arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w),x)
    return y

def direct_idft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n,dtype=cdouble)
    w = arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w),x)/n
    return y

def direct_dftn(x):
    x = asarray(x)
    for axis in range(len(x.shape)):
        x = fft(x,axis=axis)
    return x

def direct_idftn(x):
    x = asarray(x)
    for axis in range(len(x.shape)):
        x = ifft(x,axis=axis)
    return x

def direct_rdft(x):
    x = asarray(x)
    n = len(x)
    w = -arange(n)*(2j*pi/n)
    r = zeros(n,dtype=double)
    for i in range(n/2+1):
        y = dot(exp(i*w),x)
        if i:
            r[2*i-1] = y.real
            if 2*i<n:
                r[2*i] = y.imag
        else:
            r[0] = y.real
    return r

def direct_irdft(x):
    x = asarray(x)
    n = len(x)
    x1 = zeros(n,dtype=cdouble)
    for i in range(n/2+1):
        if i:
            if 2*i<n:
                x1[i] = x[2*i-1] + 1j* x[2*i]
                x1[n-i] = x[2*i-1] - 1j* x[2*i]
            else:
                x1[i] = x[2*i-1]
        else:
            x1[0] = x[0]
    return direct_idft(x1).real

class TestFft(TestCase):

    def test_definition(self):
        x = [1,2,3,4+1j,1,2,3,4+2j]
        y = fft(x)
        y1 = direct_dft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4+0j,5]
        assert_array_almost_equal(fft(x),direct_dft(x))

    def test_n_argument_real(self):
        x1 = [1,2,3,4]
        x2 =  [1,2,3,4]
        y = fft([x1,x2],n=4)
        assert_equal(y.shape,(2,4))
        assert_array_almost_equal(y[0],direct_dft(x1))
        assert_array_almost_equal(y[1],direct_dft(x2))

    def _test_n_argument_complex(self):
        x1 = [1,2,3,4+1j]
        x2 =  [1,2,3,4+1j]
        y = fft([x1,x2],n=4)
        assert_equal(y.shape,(2,4))
        assert_array_almost_equal(y[0],direct_dft(x1))
        assert_array_almost_equal(y[1],direct_dft(x2))

    def test_djbfft(self):
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y = fftpack.zfft(x)
            y2 = numpy.fft.fft(x)
            assert_array_almost_equal(y,y2)
            y = fftpack.zrfft(x)
            assert_array_almost_equal(y,y2)


class TestIfft(TestCase):

    def test_definition(self):
        x = [1,2,3,4+1j,1,2,3,4+2j]
        y = ifft(x)
        y1 = direct_idft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4+0j,5]
        assert_array_almost_equal(ifft(x),direct_idft(x))
        x = [1,2,3,4,1,2,3,4]
        y = ifft(x)
        y1 = direct_idft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4,5]
        assert_array_almost_equal(ifft(x),direct_idft(x))

    def test_djbfft(self):
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y = fftpack.zfft(x,direction=-1)
            y2 = numpy.fft.ifft(x)
            assert_array_almost_equal(y,y2)
            y = fftpack.zrfft(x,direction=-1)
            assert_array_almost_equal(y,y2)

    def test_random_complex(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype(cdouble)
            x = random([size]).astype(cdouble) +1j*x
            assert_array_almost_equal (ifft(fft(x)),x)
            assert_array_almost_equal (fft(ifft(x)),x)

    def test_random_real(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype(double)
            assert_array_almost_equal (ifft(fft(x)),x)
            assert_array_almost_equal (fft(ifft(x)),x)


class TestRfft(TestCase):

    def test_definition(self):
        x = [1,2,3,4,1,2,3,4]
        y = rfft(x)
        y1 = direct_rdft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4,1,2,3,4,5]
        y = rfft(x)
        y1 = direct_rdft(x)
        assert_array_almost_equal(y,y1)

    def test_djbfft(self):
        from numpy.fft import fft as numpy_fft
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y2 = numpy_fft(x)
            y1 = zeros((n,),dtype=double)
            y1[0] = y2[0].real
            y1[-1] = y2[n/2].real
            for k in range(1,n/2):
                y1[2*k-1] = y2[k].real
                y1[2*k] = y2[k].imag
            y = fftpack.drfft(x)
            assert_array_almost_equal(y,y1)


class TestIrfft(TestCase):

    def test_definition(self):
        x = [1,2,3,4,1,2,3,4]
        x1 = [1,2+3j,4+1j,2+3j,4,2-3j,4-1j,2-3j]
        y = irfft(x)
        y1 = direct_irdft(x)
        assert_array_almost_equal(y,y1)
        assert_array_almost_equal(y,ifft(x1))
        x = [1,2,3,4,1,2,3,4,5]
        x1 = [1,2+3j,4+1j,2+3j,4+5j,4-5j,2-3j,4-1j,2-3j]
        y = irfft(x)
        y1 = direct_irdft(x)
        assert_array_almost_equal(y,y1)
        assert_array_almost_equal(y,ifft(x1))

    def test_djbfft(self):
        from numpy.fft import ifft as numpy_ifft
        for i in range(2,14):
            n = 2**i
            x = range(n)
            x1 = zeros((n,),dtype=cdouble)
            x1[0] = x[0]
            for k in range(1,n/2):
                x1[k] = x[2*k-1]+1j*x[2*k]
                x1[n-k] = x[2*k-1]-1j*x[2*k]
            x1[n/2] = x[-1]
            y1 = numpy_ifft(x1)
            y = fftpack.drfft(x,direction=-1)
            assert_array_almost_equal(y,y1)

    def test_random_real(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype(double)
            assert_array_almost_equal (irfft(rfft(x)),x)
            assert_array_almost_equal (rfft(irfft(x)),x)

class Testfft2(TestCase):
    def test_regression_244(self):
        """fft returns wrong result with axes parameter."""
        # fftn (and hence fft2) used to break when both axes and shape were
        # used
        x = numpy.ones((4,4,2))
        y = fft2(x, shape=(8,8), axes=(-3,-2))
        y_r = numpy.fft.fftn(x, s=(8, 8), axes=(-3,  -2))
        assert_array_almost_equal(y, y_r)

class TestFftn(TestCase):

    def test_definition(self):
        x = [[1,2,3],[4,5,6],[7,8,9]]
        y = fftn(x)
        assert_array_almost_equal(y,direct_dftn(x))
        x = random((20,26))
        assert_array_almost_equal(fftn(x),direct_dftn(x))
        x = random((5,4,3,20))
        assert_array_almost_equal(fftn(x),direct_dftn(x))

    def test_axes_argument(self):
        #plane == ji_plane, x== kji_space
        plane1 = [[1,2,3],[4,5,6],[7,8,9]]
        plane2 = [[10,11,12],[13,14,15],[16,17,18]]
        plane3 = [[19,20,21],[22,23,24],[25,26,27]]
        ki_plane1 = [[1,2,3],[10,11,12],[19,20,21]]
        ki_plane2 = [[4,5,6],[13,14,15],[22,23,24]]
        ki_plane3 = [[7,8,9],[16,17,18],[25,26,27]]
        jk_plane1 = [[1,10,19],[4,13,22],[7,16,25]]
        jk_plane2 = [[2,11,20],[5,14,23],[8,17,26]]
        jk_plane3 = [[3,12,21],[6,15,24],[9,18,27]]
        kj_plane1 = [[1,4,7],[10,13,16],[19,22,25]]
        kj_plane2 = [[2,5,8],[11,14,17],[20,23,26]]
        kj_plane3 = [[3,6,9],[12,15,18],[21,24,27]]
        ij_plane1 = [[1,4,7],[2,5,8],[3,6,9]]
        ij_plane2 = [[10,13,16],[11,14,17],[12,15,18]]
        ij_plane3 = [[19,22,25],[20,23,26],[21,24,27]]
        ik_plane1 = [[1,10,19],[2,11,20],[3,12,21]]
        ik_plane2 = [[4,13,22],[5,14,23],[6,15,24]]
        ik_plane3 = [[7,16,25],[8,17,26],[9,18,27]]
        ijk_space = [jk_plane1,jk_plane2,jk_plane3]
        ikj_space = [kj_plane1,kj_plane2,kj_plane3]
        jik_space = [ik_plane1,ik_plane2,ik_plane3]
        jki_space = [ki_plane1,ki_plane2,ki_plane3]
        kij_space = [ij_plane1,ij_plane2,ij_plane3]
        x = array([plane1,plane2,plane3])

        assert_array_almost_equal(fftn(x),fftn(x,axes=(-3,-2,-1))) # kji_space
        assert_array_almost_equal(fftn(x),fftn(x,axes=(0,1,2)))
        y = fftn(x,axes=(2,1,0)) # ijk_space
        assert_array_almost_equal(swapaxes(y,-1,-3),fftn(ijk_space))
        y = fftn(x,axes=(2,0,1)) # ikj_space
        assert_array_almost_equal(swapaxes(swapaxes(y,-1,-3),
                                                   -1,-2)
                                  ,fftn(ikj_space))
        y = fftn(x,axes=(1,2,0)) # jik_space
        assert_array_almost_equal(swapaxes(swapaxes(y,-1,-3),
                                                   -3,-2)
                                  ,fftn(jik_space))
        y = fftn(x,axes=(1,0,2)) # jki_space
        assert_array_almost_equal(swapaxes(y,-2,-3),fftn(jki_space))
        y = fftn(x,axes=(0,2,1)) # kij_space
        assert_array_almost_equal(swapaxes(y,-2,-1),
                                  fftn(kij_space))

        y = fftn(x,axes=(-2,-1)) # ji_plane
        assert_array_almost_equal(fftn(plane1),y[0])
        assert_array_almost_equal(fftn(plane2),y[1])
        assert_array_almost_equal(fftn(plane3),y[2])
        y = fftn(x,axes=(1,2)) # ji_plane
        assert_array_almost_equal(fftn(plane1),y[0])
        assert_array_almost_equal(fftn(plane2),y[1])
        assert_array_almost_equal(fftn(plane3),y[2])
        y = fftn(x,axes=(-3,-2)) # kj_plane
        assert_array_almost_equal(fftn(x[:,:,0]),y[:,:,0])
        assert_array_almost_equal(fftn(x[:,:,1]),y[:,:,1])
        assert_array_almost_equal(fftn(x[:,:,2]),y[:,:,2])
        y = fftn(x,axes=(-3,-1)) # ki_plane
        assert_array_almost_equal(fftn(x[:,0,:]),y[:,0,:])
        assert_array_almost_equal(fftn(x[:,1,:]),y[:,1,:])
        assert_array_almost_equal(fftn(x[:,2,:]),y[:,2,:])
        y = fftn(x,axes=(-1,-2)) # ij_plane
        assert_array_almost_equal(fftn(ij_plane1),swapaxes(y[0],-2,-1))
        assert_array_almost_equal(fftn(ij_plane2),swapaxes(y[1],-2,-1))
        assert_array_almost_equal(fftn(ij_plane3),swapaxes(y[2],-2,-1))
        y = fftn(x,axes=(-1,-3)) # ik_plane
        assert_array_almost_equal(fftn(ik_plane1),swapaxes(y[:,0,:],-1,-2))
        assert_array_almost_equal(fftn(ik_plane2),swapaxes(y[:,1,:],-1,-2))
        assert_array_almost_equal(fftn(ik_plane3),swapaxes(y[:,2,:],-1,-2))
        y = fftn(x,axes=(-2,-3)) # jk_plane
        assert_array_almost_equal(fftn(jk_plane1),swapaxes(y[:,:,0],-1,-2))
        assert_array_almost_equal(fftn(jk_plane2),swapaxes(y[:,:,1],-1,-2))
        assert_array_almost_equal(fftn(jk_plane3),swapaxes(y[:,:,2],-1,-2))

        y = fftn(x,axes=(-1,)) # i_line
        for i in range(3):
            for j in range(3):
                assert_array_almost_equal(fft(x[i,j,:]),y[i,j,:])
        y = fftn(x,axes=(-2,)) # j_line
        for i in range(3):
            for j in range(3):
                assert_array_almost_equal(fft(x[i,:,j]),y[i,:,j])
        y = fftn(x,axes=(0,)) # k_line
        for i in range(3):
            for j in range(3):
                assert_array_almost_equal(fft(x[:,i,j]),y[:,i,j])

        y = fftn(x,axes=()) # point
        assert_array_almost_equal(y,x)

    def test_shape_argument(self):
        small_x = [[1,2,3],[4,5,6]]
        large_x1 = [[1,2,3,0],[4,5,6,0],[0,0,0,0],[0,0,0,0]]
        y = fftn(small_x,shape=(4,4))
        assert_array_almost_equal (y,fftn(large_x1))
        y = fftn(small_x,shape=(3,4))
        assert_array_almost_equal (y,fftn(large_x1[:-1]))

    def test_shape_axes_argument(self):
        small_x = [[1,2,3],[4,5,6],[7,8,9]]
        large_x1 = array([[1,2,3,0],
                                  [4,5,6,0],
                                  [7,8,9,0],
                                  [0,0,0,0]])
        # Disable tests with shape and axes of different lengths
        #y = fftn(small_x,shape=(4,4),axes=(-1,))
        #for i in range(4):
        #    assert_array_almost_equal (y[i],fft(large_x1[i]))
        #y = fftn(small_x,shape=(4,4),axes=(-2,))
        #for i in range(4):
        #    assert_array_almost_equal (y[:,i],fft(large_x1[:,i]))
        y = fftn(small_x,shape=(4,4),axes=(-2,-1))
        assert_array_almost_equal (y,fftn(large_x1))
        y = fftn(small_x,shape=(4,4),axes=(-1,-2))
        assert_array_almost_equal (y,swapaxes(\
            fftn(swapaxes(large_x1,-1,-2)),-1,-2))

    def test_shape_axes_argument2(self):
        # Change shape of the last axis
        x = numpy.random.random((10, 5, 3, 7))
        y = fftn(x, axes=(-1,), shape=(8,))
        assert_array_almost_equal(y, fft(x, axis=-1, n=8))

        # Change shape of an arbitrary axis which is not the last one
        x = numpy.random.random((10, 5, 3, 7))
        y = fftn(x, axes=(-2,), shape=(8,))
        assert_array_almost_equal(y, fft(x, axis=-2, n=8))

        # Change shape of axes: cf #244, where shape and axes were mixed up
        x = numpy.random.random((4,4,2))
        y = fftn(x, axes=(-3,-2), shape=(8,8))
        assert_array_almost_equal(y, numpy.fft.fftn(x, axes=(-3, -2), s=(8, 8)))

    def test_shape_argument_more(self):
        # Test that fftn raise a value error exception when s.shape is longer
        # than x.shape
        x = zeros((4, 4, 2))
        try:
            fx = fftn(x, shape = (8, 8, 2, 1))
            raise AssertionError("s.shape longer than x.shape succeded, "\
                                 "but should not have.")
        except ValueError:
            pass

class TestIfftn(TestCase):

    def test_definition(self):
        x = [[1,2,3],[4,5,6],[7,8,9]]
        y = ifftn(x)
        assert_array_almost_equal(y,direct_idftn(x))
        x = random((20,26))
        assert_array_almost_equal(ifftn(x),direct_idftn(x))
        x = random((5,4,3,20))
        assert_array_almost_equal(ifftn(x),direct_idftn(x))

    def test_random_complex(self):
        for size in [1,2,51,32,64,92]:
            x = random([size,size]) + 1j*random([size,size])
            assert_array_almost_equal (ifftn(fftn(x)),x)
            assert_array_almost_equal (fftn(ifftn(x)),x)

if __name__ == "__main__":
    run_module_suite()
