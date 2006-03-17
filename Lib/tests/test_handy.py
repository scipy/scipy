""" Test functions for handy module

"""

import unittest
import scipy.misc.limits as limits
from numpy.testing import assert_array_equal, assert_equal
from numpy.testing import assert_almost_equal, assert_array_almost_equal
from scipy import *

##################################################

class test_grid(unittest.TestCase):
    def check_basic(self):
        a = grid[-1:1:10j]
        b = grid[-1:1:0.1]
        assert(a.shape == (10,))
        assert(b.shape == (20,))
        assert(a[0] == -1)
        assert_almost_equal(a[-1],1)
        assert(b[0] == -1)
        assert_almost_equal(b[1]-b[0],0.1,11)
        assert_almost_equal(b[-1],b[0]+19*0.1,11)
        assert_almost_equal(a[1]-a[0],2.0/9.0,11)

    def check_nd(self):
        c = grid[-1:1:10j,-2:2:10j]
        d = grid[-1:1:0.1,-2:2:0.2]
        assert(c.shape == (2,10,10))
        assert(d.shape == (2,20,20))
        assert_array_equal(c[0][0,:],-ones(10,'d'))
        assert_array_equal(c[1][:,0],-2*ones(10,'d'))
        assert_array_almost_equal(c[0][-1,:],ones(10,'d'),11)
        assert_array_almost_equal(c[1][:,-1],2*ones(10,'d'),11)
        assert_array_almost_equal(d[0,1,:]-d[0,0,:], 0.1*ones(20,'d'),11)
        assert_array_almost_equal(d[1,:,1]-d[1,:,0], 0.2*ones(20,'d'),11)

class test_concatenator(unittest.TestCase):
    def check_1d(self):
        assert_array_equal(r_[1,2,3,4,5,6],array([1,2,3,4,5,6]))
        b = ones(5)
        c = r_[b,0,0,b]
        assert_array_equal(c,[1,1,1,1,1,0,0,1,1,1,1,1])
        c = c_[b,0,0,b]
        assert_array_equal(c,[1,1,1,1,1,0,0,1,1,1,1,1])

    def check_2d(self):
        b = rand(5,5)
        c = rand(5,5)
        d = c_[b,c]  # append columns
        assert(d.shape == (5,10))
        assert_array_equal(d[:,:5],b)
        assert_array_equal(d[:,5:],c)
        d = r_[b,c]
        assert(d.shape == (10,5))
        assert_array_equal(d[:5,:],b)
        assert_array_equal(d[5:,:],c)

class test_logspace(unittest.TestCase):
    def check_basic(self):
        y = logspace(0,6)
        assert(len(y)==50)
        y = logspace(0,6,num=100)
        assert(y[-1] == 10**6)
        y = logspace(0,6,endpoint=0)
        assert(y[-1] < 10**6)
        y = logspace(0,6,num=7)
        assert_array_equal(y,[1,10,100,1e3,1e4,1e5,1e6])

class test_linspace(unittest.TestCase):
    def check_basic(self):
        y = linspace(0,10)
        assert(len(y)==50)
        y = linspace(2,10,num=100)
        assert(y[-1] == 10)
        y = linspace(2,10,endpoint=0)
        assert(y[-1] < 10)
        y,st = linspace(2,10,retstep=1)
        assert_almost_equal(st,8/49.0)
        assert_array_almost_equal(y,grid[2:10:50j],13)

class test_any(unittest.TestCase):
    def check_basic(self):
        y1 = [0,0,1,0]
        y2 = [0,0,0,0]
        y3 = [1,0,1,0]
        assert(any(y1))
        assert(any(y3))
        assert(not any(y2))

    def check_nd(self):
        y1 = [[0,0,0],[0,1,0],[1,1,0]]
        assert_array_equal(any(y1),[1,1,0])
        assert_array_equal(any(y1,axis=1),[0,1,1])

class test_all(unittest.TestCase):
    def check_basic(self):
        y1 = [0,1,1,0]
        y2 = [0,0,0,0]
        y3 = [1,1,1,1]
        assert(not all(y1))
        assert(all(y3))
        assert(not all(y2))
        assert(all(~array(y2)))

    def check_nd(self):
        y1 = [[0,0,1],[0,1,1],[1,1,1]]
        assert_array_equal(all(y1),[0,0,1])
        assert_array_equal(all(y1,axis=1),[0,0,1])

class test_fix(unittest.TestCase):
    def check_basic(self):
        y1 = [-3.5,4.2,-0.0,0.0,5.5,-1.7]
        assert_array_almost_equal(fix(y1),[-3,4,0.0,0.0,5,-1.0],13)

class test_fftshift(unittest.TestCase):
    def check_basic(self):
        y1 = grid[:10]
        y2 = grid[:11]
        y1s = fftshift(y1)
        y2s = fftshift(y2)
        assert_array_equal(y1s[:5],y1[5:])
        assert_array_equal(y1s[5:],y1[:5])
        assert_array_equal(y2s[:5],y2[6:])
        assert_array_equal(y2s[5:],y2[:6])

    def check_nd(self):
        y1 = grid[0:10,0:11,0:13][0]
        y1s = fftshift(y1,(0,2))
        assert_array_equal(y1s[0],5*ones((11,13),'d'))
        assert_array_equal(y1s[:,0,0],fftshift(grid[0:10]))
        assert_array_equal(y1,ifftshift(y1s,(0,2)))


class test_ifftshift(unittest.TestCase):
    def check_basic(self):
        for k in range(5,11):
            y1 = rand(k)
            assert_array_equal(y1,ifftshift(fftshift(y1)))

    def check_nd(self):
        for k in range(5,11):
            y1 = rand(k-1,k,k+1)
            assert_array_equal(y1,ifftshift(fftshift(y1,(0,1)),(0,1)))
            assert_array_equal(y1,ifftshift(fftshift(y1,(1,2)),(1,2)))

class test_fftfreq(unittest.TestCase):
    def check_basic(self):
        for k in range(4,11):
            y = fftfreq(k)*k
            lower = -(k/2)
            upper = (k-1)/2
            assert_array_equal(fftshift(y),grid[lower:upper+1])

class test_cont_ft(unittest.TestCase):
    def check_basic(self):
        n = arange(100) - 49.5
        x = grid[-2:2:100j]
        dx = x[2] - x[1]
        gn = exp(-pi*x**2)
        fr = grid[-2:2:20j]
        actual = exp(-pi*fr**2)
        comp = cont_ft(gn,fr,delta=dx,n=n).real
        assert_array_almost_equal(comp,actual,6)

class test_atleast_1d(unittest.TestCase):
    def check_basic(self):
        assert(atleast_1d(3).shape == (1,))
        assert(atleast_1d(3j).shape == (1,))
        assert(atleast_1d(3L).shape == (1,))
        assert(atleast_1d(3.0).shape == (1,))
        assert(atleast_1d([[2,3],[4,5]]).shape == (2,2))

class test_r2array(unittest.TestCase):
    def check_basic(self):
        assert(r2array(3).shape == (1,1))
        assert(r2array([3j,1]).shape == (1,2))
        assert(r2array([[[3,1],[4,5]],[[3,5],[1,2]]]).shape == (2,2,2))

class test_isscalar(unittest.TestCase):
    def check_basic(self):
        assert(isscalar(3))
        assert(not isscalar([3]))
        assert(not isscalar((3,)))
        assert(isscalar(3j))
        assert(isscalar(10L))
        assert(isscalar(4.0))

class test_toeplitz(unittest.TestCase):
    def check_basic(self):
        y = toeplitz([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,1,2],[3,2,1]])
        y = toeplitz([1,2,3],[1,4,5])
        assert_array_equal(y,[[1,4,5],[2,1,4],[3,2,1]])

class test_hankel(unittest.TestCase):
    def check_basic(self):
        y = hankel([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,3,0],[3,0,0]])
        y = hankel([1,2,3],[3,4,5])
        assert_array_equal(y,[[1,2,3],[2,3,4],[3,4,5]])


class test_real_if_close(unittest.TestCase):
    def check_basic(self):
        a = randn(10)
        b = real_if_close(a+1e-15j)
        assert(isrealobj(b))
        assert_array_equal(a,b)
        b = real_if_close(a+1e-7j)
        assert(iscomplexobj(b))
        b = real_if_close(a+1e-7j,tol=1e-6)
        assert(isrealobj(b))

class test_polyint(unittest.TestCase):
    def check_order1_noconstant(self):
        a = rand(10)*4


##################################################

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_grid,'check_') )
        suites.append( unittest.makeSuite(test_concatenator,'check_') )
        suites.append( unittest.makeSuite(test_logspace,'check_') )
        suites.append( unittest.makeSuite(test_linspace,'check_') )
        suites.append( unittest.makeSuite(test_any,'check_') )
        suites.append( unittest.makeSuite(test_all,'check_') )
        suites.append( unittest.makeSuite(test_fix,'check_') )
        suites.append( unittest.makeSuite(test_fftshift,'check_') )
        suites.append( unittest.makeSuite(test_ifftshift,'check_') )
        suites.append( unittest.makeSuite(test_fftfreq,'check_') )
        suites.append( unittest.makeSuite(test_cont_ft,'check_') )
        suites.append( unittest.makeSuite(test_atleast_1d,'check_') )
        suites.append( unittest.makeSuite(test_r2array,'check_') )
        suites.append( unittest.makeSuite(test_isscalar,'check_') )
        suites.append( unittest.makeSuite(test_toeplitz,'check_') )
        suites.append( unittest.makeSuite(test_hankel,'check_') )
        suites.append( unittest.makeSuite(test_real_if_close,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
