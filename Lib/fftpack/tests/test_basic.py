#!/usr/bin/env python
# Created by Pearu Peterson, September 2002
""" Test functions for fftpack.basic module
"""
__usage__ = """
Build fftpack:
  python setup_fftpack.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.fftpack.test(<level>)'
Run tests if fftpack is not installed:
  python tests/test_basic.py [<level>]
"""
import sys
from scipy_test.testing import *
set_package_path()
from fftpack import ifft,fft,fftn,ifftn,rfft,irfft
from fftpack import _fftpack as fftpack
del sys.path[0]


import Numeric
from Numeric import arange, add, array

from scipy_test.testing import rand
def random(size):
    return rand(*size)

import unittest

def get_mat(n):
    data = arange(n)
    data = add.outer(data,data)
    return data

def direct_dft(x):
    x = Numeric.asarray(x)
    n = len(x)
    y = Numeric.zeros(n,'D')
    w = -Numeric.arange(n)*(2j*Numeric.pi/n)
    for i in range(n):
        y[i] = Numeric.dot(Numeric.exp(i*w),x)
    return y

def direct_idft(x):
    x = Numeric.asarray(x)
    n = len(x)
    y = Numeric.zeros(n,'D')
    w = Numeric.arange(n)*(2j*Numeric.pi/n)
    for i in range(n):
        y[i] = Numeric.dot(Numeric.exp(i*w),x)/n
    return y

def direct_dftn(x):
    x = Numeric.asarray(x)
    for axis in range(len(x.shape)):
        x = fft(x,axis=axis)
    return x

def direct_idftn(x):
    x = Numeric.asarray(x)
    for axis in range(len(x.shape)):
        x = ifft(x,axis=axis)
    return x

def direct_rdft(x):
    x = Numeric.asarray(x)
    n = len(x)
    w = -Numeric.arange(n)*(2j*Numeric.pi/n)
    r = Numeric.zeros(n,'d')
    for i in range(n/2+1):
        y = Numeric.dot(Numeric.exp(i*w),x)
        if i:
            r[2*i-1] = y.real
            if 2*i<n:
                r[2*i] = y.imag
        else:
            r[0] = y.real
    return r

def direct_irdft(x):
    x = Numeric.asarray(x)
    n = len(x)
    x1 = Numeric.zeros(n,'D')
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

class test_fft(ScipyTestCase):

    def check_definition(self):
        x = [1,2,3,4+1j,1,2,3,4+2j]
        y = fft(x)
        y1 = direct_dft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4+0j,5]
        assert_array_almost_equal(fft(x),direct_dft(x))

    def check_n_argument_real(self):
        x1 = [1,2,3,4]
        x2 =  [1,2,3,4]
        y = fft([x1,x2],n=4)
        assert_equal(y.shape,(2,4))
        assert_array_almost_equal(y[0],direct_dft(x1))
        assert_array_almost_equal(y[1],direct_dft(x2))

    def check_n_argument_complex(self):
        x1 = [1,2,3,4+1j]
        x2 =  [1,2,3,4+1j]
        y = fft([x1,x2],n=4)
        assert_equal(y.shape,(2,4))
        assert_array_almost_equal(y[0],direct_dft(x1))
        assert_array_almost_equal(y[1],direct_dft(x2))

    def check_djbfft(self):
        try:
            from FFT import fft as fft2
        except ImportError:
            print 'Skipping check_djbfft (failed to import FFT)'
            return
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y = fftpack.zfft(x)
            y2 = fft2(x)
            assert_array_almost_equal(y,y2)
            y = fftpack.zrfft(x)
            assert_array_almost_equal(y,y2)

    def bench_random(self,level=5):
        try:
            from FFT import fft as FFT_fft
        except ImportError:
            FFT_fft = None
        print
        print '                 Fast Fourier Transform'
        print '================================================='
        print '      |    real input     |   complex input    '
        print '-------------------------------------------------'
        print ' size |  scipy  | Numeric |  scipy  | Numeric '
        print '-------------------------------------------------'
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print '%5s' % size,
            sys.stdout.flush()

            for x in [random([size]).astype('d'),
                      random([size]).astype('D')+random([size]).astype('D')*1j
                      ]:
                if size > 500: y = fft(x)
                else: y = direct_dft(x)                
                assert_array_almost_equal(fft(x),y)
                print '|%8.2f' % self.measure('fft(x)',repeat),
                sys.stdout.flush()
                
                if FFT_fft is not None:
                    assert_array_almost_equal(FFT_fft(x),y)
                    print '|%8.2f' % self.measure('FFT_fft(x)',repeat),
                else:
                    print '|  N/A  ',
                sys.stdout.flush()

            print ' (secs for %s calls)' % (repeat)


class test_ifft(ScipyTestCase):

    def check_definition(self):
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

    def check_djbfft(self):
        try:
            from FFT import inverse_fft as ifft2
        except ImportError:
            return
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y = fftpack.zfft(x,direction=-1)
            y2 = ifft2(x)
            assert_array_almost_equal(y,y2)
            y = fftpack.zrfft(x,direction=-1)
            assert_array_almost_equal(y,y2)

    def check_random_complex(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype('D')
            x = random([size]).astype('D') +1j*x
            assert_array_almost_equal (ifft(fft(x)),x)
            assert_array_almost_equal (fft(ifft(x)),x)

    def check_random_real(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype('d')
            assert_array_almost_equal (ifft(fft(x)),x)
            assert_array_almost_equal (fft(ifft(x)),x)

    def bench_random(self,level=5):
        try:
            from FFT import inverse_fft as FFT_ifft
        except ImportError:
            FFT_ifft = None

        print
        print '       Inverse Fast Fourier Transform'
        print '==============================================='
        print '      |     real input    |    complex input   '
        print '-----------------------------------------------'
        print ' size |  scipy  | Numeric |  scipy  | Numeric  '
        print '-----------------------------------------------'
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print '%5s' % size,
            sys.stdout.flush()

            for x in [random([size]).astype('d'),
                      random([size]).astype('D')+random([size]).astype('D')*1j
                      ]:
                if size > 500: y = ifft(x)
                else: y = direct_idft(x)                
                assert_array_almost_equal(ifft(x),y)
                print '|%8.2f' % self.measure('ifft(x)',repeat),
                sys.stdout.flush()

                if FFT_ifft is not None:
                    assert_array_almost_equal(FFT_ifft(x),y)
                    print '|%8.2f' % self.measure('FFT_ifft(x)',repeat),
                else:
                    print '|  N/A  ',
                sys.stdout.flush()

            print ' (secs for %s calls)' % (repeat)


class test_rfft(ScipyTestCase):

    def check_definition(self):
        x = [1,2,3,4,1,2,3,4]
        y = rfft(x)
        y1 = direct_rdft(x)
        assert_array_almost_equal(y,y1)
        x = [1,2,3,4,1,2,3,4,5]
        y = rfft(x)
        y1 = direct_rdft(x)
        assert_array_almost_equal(y,y1)

    def check_djbfft(self):
        try:
            from FFT import fft as fft2
        except ImportError,errmess:
            print 'Skipping check_djbfft (failed to import FFT: %s)' % errmess
            return
        for i in range(2,14):
            n = 2**i
            x = range(n)
            y2 = fft2(x)
            y1 = Numeric.zeros((n,),'d')
            y1[0] = y2[0].real
            y1[-1] = y2[n/2].real
            for k in range(1,n/2):
                y1[2*k-1] = y2[k].real
                y1[2*k] = y2[k].imag
            y = fftpack.drfft(x)
            assert_array_almost_equal(y,y1)

    def bench_random(self,level=5):
        try:
            from FFT import real_fft as FFT_rfft
        except ImportError:
            FFT_rfft = None

        print
        print 'Fast Fourier Transform (real data)'
        print '=================================='
        print ' size |  scipy  | Numeric  '
        print '----------------------------------'
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print '%5s' % size,
            sys.stdout.flush()

            x = random([size]).astype('d')
            print '|%8.2f' % self.measure('rfft(x)',repeat),
            sys.stdout.flush()
            
            if FFT_rfft is not None:
                print '|%8.2f' % self.measure('FFT_rfft(x)',repeat),
            else:
                print '|  N/A  ',
            sys.stdout.flush()

            print ' (secs for %s calls)' % (repeat)


class test_irfft(ScipyTestCase):

    def check_definition(self):
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

    def check_djbfft(self):
        try:
            from FFT import inverse_fft as ifft2
        except ImportError:
            print 'Skipping check_djbfft (failed to import FFT)'
            return
        for i in range(2,14):
            n = 2**i
            x = range(n)
            x1 = Numeric.zeros((n,),'D')
            x1[0] = x[0]
            for k in range(1,n/2):
                x1[k] = x[2*k-1]+1j*x[2*k]
                x1[n-k] = x[2*k-1]-1j*x[2*k]
            x1[n/2] = x[-1]
            y1 = ifft2(x1)
            y = fftpack.drfft(x,direction=-1)
            assert_array_almost_equal(y,y1)

    def check_random_real(self):
        for size in [1,51,111,100,200,64,128,256,1024]:
            x = random([size]).astype('d')
            assert_array_almost_equal (irfft(rfft(x)),x)
            assert_array_almost_equal (rfft(irfft(x)),x)

    def bench_random(self,level=5):
        try:
            from FFT import inverse_real_fft as FFT_irfft
        except ImportError:
            FFT_irfft = None

        print
        print 'Inverse Fast Fourier Transform (real data)'
        print '=================================='
        print ' size |  scipy  | Numeric  '
        print '----------------------------------'
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print '%5s' % size,
            sys.stdout.flush()

            x = random([size]).astype('d')
            x1 = Numeric.zeros(size/2+1,'D')
            x1[0] = x[0]
            for i in range(1,size/2):
                x1[i] = x[2*i-1] + 1j * x[2*i]
            if not size%2:
                x1[-1] = x[-1]
            y = irfft(x)

            print '|%8.2f' % self.measure('irfft(x)',repeat),
            sys.stdout.flush()
            
            if FFT_irfft is not None:
                assert_array_almost_equal(FFT_irfft(x1,size),y)
                print '|%8.2f' % self.measure('FFT_irfft(x1,size)',repeat),
            else:
                print '|  N/A  ',
            sys.stdout.flush()

            print ' (secs for %s calls)' % (repeat)

class test_fftn(ScipyTestCase):

    def check_definition(self):
        x = [[1,2,3],[4,5,6],[7,8,9]]
        y = fftn(x)
        assert_array_almost_equal(y,direct_dftn(x))
        x = random((20,26))
        assert_array_almost_equal(fftn(x),direct_dftn(x))
        x = random((5,4,3,20))
        assert_array_almost_equal(fftn(x),direct_dftn(x))

    def check_axes_argument(self):
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
        x = Numeric.array([plane1,plane2,plane3])

        assert_array_almost_equal(fftn(x),fftn(x,axes=(-3,-2,-1))) # kji_space
        assert_array_almost_equal(fftn(x),fftn(x,axes=(0,1,2)))
        y = fftn(x,axes=(2,1,0)) # ijk_space
        assert_array_almost_equal(Numeric.swapaxes(y,-1,-3),fftn(ijk_space))
        y = fftn(x,axes=(2,0,1)) # ikj_space
        assert_array_almost_equal(Numeric.swapaxes(Numeric.swapaxes(y,-1,-3),
                                                   -1,-2)
                                  ,fftn(ikj_space))
        y = fftn(x,axes=(1,2,0)) # jik_space
        assert_array_almost_equal(Numeric.swapaxes(Numeric.swapaxes(y,-1,-3),
                                                   -3,-2)
                                  ,fftn(jik_space))
        y = fftn(x,axes=(1,0,2)) # jki_space
        assert_array_almost_equal(Numeric.swapaxes(y,-2,-3),fftn(jki_space))
        y = fftn(x,axes=(0,2,1)) # kij_space
        assert_array_almost_equal(Numeric.swapaxes(y,-2,-1),
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
        assert_array_almost_equal(fftn(ij_plane1),Numeric.swapaxes(y[0],-2,-1))
        assert_array_almost_equal(fftn(ij_plane2),Numeric.swapaxes(y[1],-2,-1))
        assert_array_almost_equal(fftn(ij_plane3),Numeric.swapaxes(y[2],-2,-1))
        y = fftn(x,axes=(-1,-3)) # ik_plane
        assert_array_almost_equal(fftn(ik_plane1),Numeric.swapaxes(y[:,0,:],-1,-2))
        assert_array_almost_equal(fftn(ik_plane2),Numeric.swapaxes(y[:,1,:],-1,-2))
        assert_array_almost_equal(fftn(ik_plane3),Numeric.swapaxes(y[:,2,:],-1,-2))
        y = fftn(x,axes=(-2,-3)) # jk_plane
        assert_array_almost_equal(fftn(jk_plane1),Numeric.swapaxes(y[:,:,0],-1,-2))
        assert_array_almost_equal(fftn(jk_plane2),Numeric.swapaxes(y[:,:,1],-1,-2))
        assert_array_almost_equal(fftn(jk_plane3),Numeric.swapaxes(y[:,:,2],-1,-2))

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

    def check_shape_argument(self):
        small_x = [[1,2,3],[4,5,6]]
        large_x1 = [[1,2,3,0],[4,5,6,0],[0,0,0,0],[0,0,0,0]]
        y = fftn(small_x,shape=(4,4))
        assert_array_almost_equal (y,fftn(large_x1))
        y = fftn(small_x,shape=(3,4))
        assert_array_almost_equal (y,fftn(large_x1[:-1]))

    def check_shape_axes_argument(self):
        small_x = [[1,2,3],[4,5,6],[7,8,9]]
        large_x1 = Numeric.array([[1,2,3,0],
                                  [4,5,6,0],
                                  [7,8,9,0],
                                  [0,0,0,0]])
        y = fftn(small_x,shape=(4,4),axes=(-1,))
        for i in range(4):
            assert_array_almost_equal (y[i],fft(large_x1[i]))
        y = fftn(small_x,shape=(4,4),axes=(-2,))
        for i in range(4):
            assert_array_almost_equal (y[:,i],fft(large_x1[:,i]))
        y = fftn(small_x,shape=(4,4),axes=(-2,-1))
        assert_array_almost_equal (y,fftn(large_x1))
        y = fftn(small_x,shape=(4,4),axes=(-1,-2))
        assert_array_almost_equal (y,Numeric.swapaxes(\
            fftn(Numeric.swapaxes(large_x1,-1,-2)),-1,-2))

    def bench_random(self,level=5):
        try:
            from FFT import fftnd as FFT_fftn
        except ImportError:
            FFT_fftn = None

        print
        print '    Multi-dimensional Fast Fourier Transform'
        print '==================================================='
        print '          |    real input     |   complex input    '
        print '---------------------------------------------------'
        print '   size   |  scipy  | Numeric |  scipy  |  Numeric '
        print '---------------------------------------------------'
        for size,repeat in [((100,100),100),((1000,100),7),
                            ((256,256),10),
                            ((512,512),3),
                            ]:
            print '%9s' % ('%sx%s'%size),
            sys.stdout.flush()

            for x in [random(size).astype('d'),
                      random(size).astype('D')+random(size).astype('D')*1j
                      ]:
                y = fftn(x)
                #if size > 500: y = fftn(x)
                #else: y = direct_dft(x)                
                assert_array_almost_equal(fftn(x),y)
                print '|%8.2f' % self.measure('fftn(x)',repeat),
                sys.stdout.flush()
                
                if FFT_fftn is not None:
                    assert_array_almost_equal(FFT_fftn(x),y)
                    print '|%8.2f' % self.measure('FFT_fftn(x)',repeat),
                else:
                    print '|  N/A  ',
                sys.stdout.flush()

            print ' (secs for %s calls)' % (repeat)


class test_ifftn(ScipyTestCase):

    def check_definition(self):
        x = [[1,2,3],[4,5,6],[7,8,9]]
        y = ifftn(x)
        assert_array_almost_equal(y,direct_idftn(x))
        x = random((20,26))
        assert_array_almost_equal(ifftn(x),direct_idftn(x))
        x = random((5,4,3,20))
        assert_array_almost_equal(ifftn(x),direct_idftn(x))

    def check_random_complex(self):
        for size in [1,2,51,32,64,92]:
            x = random([size,size]) + 1j*random([size,size])
            assert_array_almost_equal (ifftn(fftn(x)),x)
            assert_array_almost_equal (fftn(ifftn(x)),x)

if __name__ == "__main__":
    ScipyTest('fftpack.basic').run()
