""" Test functions for basic module

"""

import unittest
import scipy_base.limits as limits
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_almost_equal
from scipy import *

##################################################

### Test for rand

class test_rand(unittest.TestCase):
    def __init__(self,*args,**kwds):
        unittest.TestCase.__init__(self,*args,**kwds)
        self.z = rand(100,10,20)
        
    def check_shape(self):
        assert_equal(self.z.shape,(100,10,20))

    def check_mean(self):
        mn = add.reduce(ravel(self.z)) / 200 / 100.0
        assert_almost_equal(mn,0.5,2)

    def check_std(self):
        std = add.reduce((ravel(self.z)-0.5)**2) / 200.0 / 100.0
        assert_almost_equal(std,1.0/12.0,2)

class test_randn(unittest.TestCase):
    def __init__(self,*args,**kwds):
        unittest.TestCase.__init__(self,*args,**kwds)
        self.z = randn(100,10,20)
        
    def check_shape(self):
        assert_equal(self.z.shape,(100,10,20))

    def check_mean(self):
        mn = add.reduce(ravel(self.z)) / 200 / 100.0
        assert_almost_equal(mn,0.0,1)

    def check_std(self):
        std = sqrt(add.reduce(ravel(self.z)**2) / 200.0 / 100.0)
        assert_almost_equal(std,1.0,1)


val = limits.double_resolution

class test_eye(unittest.TestCase):
    def check_basic(self):
        assert_equal(eye(4),array([[1,0,0,0],
                                   [0,1,0,0],
                                   [0,0,1,0],
                                   [0,0,0,1]]))
        assert_equal(eye(4,typecode='f'),array([[1,0,0,0],
                                                [0,1,0,0],
                                                [0,0,1,0],
                                                [0,0,0,1]],'f'))
    def check_diag(self):
        assert_equal(eye(4,k=1),array([[0,1,0,0],
                                       [0,0,1,0],
                                       [0,0,0,1],
                                       [0,0,0,0]]))
        assert_equal(eye(4,k=-1),array([[0,0,0,0],
                                        [1,0,0,0],
                                        [0,1,0,0],
                                        [0,0,1,0]]))
    def check_2d(self):
        assert_equal(eye(4,3),array([[1,0,0],
                                     [0,1,0],
                                     [0,0,1],
                                     [0,0,0]]))
        assert_equal(eye(3,4),array([[1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,1,0]]))        
    def check_diag2d(self):
        assert_equal(eye(3,4,k=2),array([[0,0,1,0],
                                         [0,0,0,1],
                                         [0,0,0,0]]))
        assert_equal(eye(4,3,k=-2),array([[0,0,0],
                                          [0,0,0],
                                          [1,0,0],
                                          [0,1,0]]))

class test_tri(unittest.TestCase):
    def check_basic(self):
        assert_equal(tri(4),array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,1,1,0],
                                   [1,1,1,1]]))
        assert_equal(tri(4,typecode='f'),array([[1,0,0,0],
                                                [1,1,0,0],
                                                [1,1,1,0],
                                                [1,1,1,1]],'f'))
    def check_diag(self):
        assert_equal(tri(4,k=1),array([[1,1,0,0],
                                       [1,1,1,0],
                                       [0,1,1,1],
                                       [1,1,1,1]]))
        assert_equal(tri(4,k=-1),array([[0,0,0,0],
                                        [1,0,0,0],
                                        [1,1,0,0],
                                        [1,1,1,0]]))
    def check_2d(self):
        assert_equal(tri(4,3),array([[1,0,0],
                                     [1,1,0],
                                     [1,1,1],
                                     [1,1,1]]))
        assert_equal(tri(3,4),array([[1,0,0,0],
                                     [1,1,0,0],
                                     [1,1,1,0]]))        
    def check_diag2d(self):
        assert_equal(tri(3,4,k=2),array([[1,1,1,0],
                                         [1,1,1,1],
                                         [1,1,1,1]]))
        assert_equal(tri(4,3,k=-2),array([[0,0,0],
                                          [0,0,0],
                                          [1,0,0],
                                          [1,1,0]]))

class test_diag(unittest.TestCase):
    def check_vector(self):
        vals = (100*rand(5)).astype('l')
        b = zeros((5,5))
        for k in range(5):
            b[k,k] = vals[k]
        assert_equal(diag(vals),b)
        b = zeros((7,7))
        c = b.copy()
        for k in range(5):
            b[k,k+2] = vals[k]
            c[k+2,k] = vals[k]
        assert_equal(diag(vals,k=2), b)
        assert_equal(diag(vals,k=-2), c)

    def check_matrix(self):
        vals = (100*rand(5,5)+1).astype('l')
        b = zeros((5,))
        for k in range(5):
            b[k] = vals[k,k]
        assert_equal(diag(vals),b)
        b = b*0
        for k in range(3):
            b[k] = vals[k,k+2]
        assert_equal(diag(vals,2),b[:3])
        for k in range(3):
            b[k] = vals[k+2,k]
        assert_equal(diag(vals,-2),b[:3])

class test_fliplr(unittest.TestCase):
    def check_basic(self):
        self.failUnlessRaises(ValueError, fliplr, ones(4))        
        self.failUnlessRaises(ValueError, fliplr, ones((4,3,2)))
        a = rand(4,4)
        b = a[:,::-1]
        assert_equal(fliplr(a),b)
        a = [[0,1,2],
             [3,4,5]]
        b = [[2,1,0],
             [5,4,3]]
        assert_equal(fliplr(a),b)

class test_flipud(unittest.TestCase):
    def check_basic(self):
        self.failUnlessRaises(ValueError, flipud, ones(4))
        self.failUnlessRaises(ValueError, flipud, ones((4,3,2)))
        a = rand(4,4)
        b = a[::-1,:]
        assert_equal(flipud(a),b)
        a = [[0,1,2],
             [3,4,5]]
        b = [[3,4,5],
             [0,1,2]]
        assert_equal(flipud(a),b)

class test_rot90(unittest.TestCase):
    def check_basic(self):
        self.failUnlessRaises(ValueError, rot90, ones(4))
        self.failUnlessRaises(ValueError, rot90, ones((4,3,2)))

        a = [[0,1,2],
             [3,4,5]]
        b1 = [[2,5],
              [1,4],
              [0,3]]
        b2 = [[5,4,3],
              [2,1,0]]
        b3 = [[3,0],
              [4,1],
              [5,2]]
        b4 = [[0,1,2],
              [3,4,5]]

        for k in range(-3,13,4):
            assert_equal(rot90(a,k=k),b1)
        for k in range(-2,13,4):
            assert_equal(rot90(a,k=k),b2)
        for k in range(-1,13,4):
            assert_equal(rot90(a,k=k),b3)
        for k in range(0,13,4):
            assert_equal(rot90(a,k=k),b4)

class test_tril(unittest.TestCase):
    def check_basic(self):
        a = (100*rand(5,5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[k,l] = 0
        assert_equal(tril(a),b)

    def check_diag(self):        
        a = (100*rand(5,5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[k,l] = 0
        assert_equal(tril(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[k,l] = 0
        assert_equal(tril(a,k=-2),b)

class test_triu(unittest.TestCase):
    def check_basic(self):
        a = (100*rand(5,5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[l,k] = 0
        assert_equal(triu(a),b)

    def check_diag(self):        
        a = (100*rand(5,5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[l,k] = 0
        assert_equal(triu(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[l,k] = 0
        assert_equal(tril(a,k=-2),b)

class test_amax(unittest.TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(amax(a),10.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(amax(b),[8.0,10.0,9.0])
        assert_equal(amax(b,axis=1),[9.0,10.0,8.0])
        
class test_amin(unittest.TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(amin(a),-5.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(amin(b),[3.0,3.0,2.0])
        assert_equal(amin(b,axis=1),[3.0,4.0,2.0])


class test_ptp(unittest.TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(ptp(a),15.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(ptp(b),[5.0,7.0,7.0])
        assert_equal(ptp(b,axis=1),[6.0,6.0,6.0])


class test_cumsum(unittest.TestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in ['1','b','s','i','l','f','d','F','D']:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            assert_array_equal(cumsum(a), array([1,3,13,24,30,35,39],ctype))
            assert_array_equal(cumsum(a2), array([[1,2,3,4],[6,8,10,13],
                                                  [16,11,14,18]],ctype))
            assert_array_equal(cumsum(a2,axis=1),
                               array([[1,3,6,10],
                                      [5,11,18,27],
                                      [10,13,17,22]],ctype))

class test_prod(unittest.TestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in ['1','b','s','i','l','f','d','F','D']:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            if ctype in ['1', 'b']:
                self.failUnlessRaises(ArithmeticError, prod, a)
                self.failUnlessRaises(ArithmeticError, prod, a2, 1)
                self.failUnlessRaises(ArithmeticError, prod, a)
            else:                
                assert_equal(prod(a),26400)
                assert_array_equal(prod(a2), array([50,36,84,180],ctype))
                assert_array_equal(prod(a2,axis=1),
                                   array([24, 1890, 600],ctype))

class test_cumprod(unittest.TestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in ['1','b','s','i','l','f','d','F','D']:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            if ctype in ['1', 'b']:
                self.failUnlessRaises(ArithmeticError, cumprod, a)
                self.failUnlessRaises(ArithmeticError, cumprod, a2, 1)
                self.failUnlessRaises(ArithmeticError, cumprod, a)
            else:                
                assert_array_equal(cumprod(a),
                                   array([1, 2, 20, 220,
                                          1320, 6600, 26400],ctype))
                assert_array_equal(cumprod(a2),
                                   array([[ 1,  2,  3,   4],
                                          [ 5, 12, 21,  36],
                                          [50, 36, 84, 180]],ctype))
                assert_array_equal(cumprod(a2,axis=1),
                                   array([[ 1,  2,   6,   24],
                                          [ 5, 30, 210, 1890],
                                          [10, 30, 120,  600]],ctype))


class test_diff(unittest.TestCase):
    def check_basic(self):
        x = [1,4,6,7,12]
        out = array([3,2,1,5])
        out2 = array([-1,-1,4])
        out3 = array([0,5])
        assert_array_equal(diff(x),out)
        assert_array_equal(diff(x,n=2),out2)
        assert_array_equal(diff(x,n=3),out3)

    def check_nd(self):
        x = 20*rand(10,20,30)
        out1 = x[:,:,1:] - x[:,:,:-1]
        out2 = out1[:,:,1:] - out1[:,:,:-1]
        out3 = x[1:,:,:] - x[:-1,:,:]
        out4 = out3[1:,:,:] - out3[:-1,:,:]
        assert_array_equal(diff(x),out1)
        assert_array_equal(diff(x,n=2),out2)
        assert_array_equal(diff(x,axis=0),out3)
        assert_array_equal(diff(x,n=2,axis=0),out4)


class test_cov(unittest.TestCase):
    def check_basic(self):
        x = [[1,1.3],
             [1.1,1.5],
             [0.9,1.8]]
        c = cov(x)
        tmp = [[0.0100, -0.0150],[-0.0150, 0.0633]]        
        assert_array_almost_equal(c,tmp,4)

    def check_twoargs(self):
        x = [1,1.1,0.9]
        y = [1.3,1.5,1.8]
        c = cov(x,y)
        z = [[1,1.3],
             [1.1,1.5],
             [0.9,1.8]]
        d = cov(z)
        e = cov([x,y],rowvar=1)
        tmp = [[0.0100, -0.0150],[-0.0150, 0.0633]]
        assert_array_almost_equal(c,tmp,4)
        assert_array_almost_equal(c,d,11)
        assert_array_almost_equal(e,d,11)

    def check_bias(self):
        z = [[1,1.3],
             [1.1,1.5],
             [0.9,1.8]]
        d = cov(z,bias=1)
        tmp = [[0.0067,-0.0100],[-0.0100, 0.0422]]
        assert_array_almost_equal(d,tmp,4)

class test_corrcoef(unittest.TestCase):
    def check_basic(self):
        x = [[1,1.3],
             [1.1,1.5],
             [0.9,1.8]]
        c = corrcoef(x)
        tmp = [[1.0000, -0.5960],[-0.5960, 1.0000]]
        assert_array_almost_equal(c,tmp,4)

    def check_twoargs(self):
        x = [1,1.1,0.9]
        y = [1.3,1.5,1.8]
        c = corrcoef(x,y)
        z = [[1,1.3],
             [1.1,1.5],
             [0.9,1.8]]
        d = corrcoef(z)
        e = corrcoef([x,y],rowvar=1)
        tmp = [[1.0000, -0.5960],[-0.5960, 1.0000]]
        assert_array_almost_equal(c,tmp,4)
        assert_array_almost_equal(c,d,11)
        assert_array_almost_equal(e,d,11)

class test_squeeze(unittest.TestCase):
    def check_basic(self):
        a = rand(20,10,10,1,1)
        b = rand(20,1,10,1,20)
        c = rand(1,1,20,10)
        assert_array_equal(squeeze(a),reshape(a,(20,10,10)))
        assert_array_equal(squeeze(b),reshape(b,(20,10,20)))
        assert_array_equal(squeeze(c),reshape(c,(20,10)))

class test_sinc(unittest.TestCase):
    def check_basic(self):
        x = r_[grid[-10:0:50j],grid[0:10:50j]]
        y = sinc(x)
        tmp = x*pi
        tmp = sin(tmp)/tmp
        tmp[50] = 1.0
        tmp[49] = 1.0
        assert_array_almost_equal(y,tmp,11)
        
class test_angle(unittest.TestCase):
    def check_basic(self):
        x = [1+3j,sqrt(2)/2.0+1j*sqrt(2)/2,1,1j,-1,-1j,1-3j,-1+3j]
        y = angle(x)
        yo = [arctan(3.0/1.0),arctan(1.0),0,pi/2,pi,-pi/2.0,-arctan(3.0/1.0),pi-arctan(3.0/1.0)]
        z = angle(x,deg=1)
        zo = array(yo)*180/pi
        assert_array_almost_equal(y,yo,11)
        assert_array_almost_equal(z,zo,11)
        
##################################################

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_rand,'check_') )
        suites.append( unittest.makeSuite(test_randn,'check_') )
        suites.append( unittest.makeSuite(test_eye,'check_') )
        suites.append( unittest.makeSuite(test_tri,'check_') )
        suites.append( unittest.makeSuite(test_diag,'check_') )
        suites.append( unittest.makeSuite(test_fliplr,'check_') )
        suites.append( unittest.makeSuite(test_flipud,'check_') )
        suites.append( unittest.makeSuite(test_rot90,'check_') )
        suites.append( unittest.makeSuite(test_triu,'check_') )
        suites.append( unittest.makeSuite(test_tril,'check_') )
        suites.append( unittest.makeSuite(test_amax,'check_') )
        suites.append( unittest.makeSuite(test_amin,'check_') )
        suites.append( unittest.makeSuite(test_ptp,'check_') )
        suites.append( unittest.makeSuite(test_mean,'check_') )
        suites.append( unittest.makeSuite(test_median,'check_') )
        suites.append( unittest.makeSuite(test_std,'check_') )
        suites.append( unittest.makeSuite(test_cumsum,'check_') )
        suites.append( unittest.makeSuite(test_prod,'check_') )
        suites.append( unittest.makeSuite(test_cumprod,'check_') )
        suites.append( unittest.makeSuite(test_trapz,'check_') )
        suites.append( unittest.makeSuite(test_diff,'check_') )
        suites.append( unittest.makeSuite(test_cov,'check_') )
        suites.append( unittest.makeSuite(test_corrcoef,'check_') )
        suites.append( unittest.makeSuite(test_squeeze,'check_') )
        suites.append( unittest.makeSuite(test_sinc,'check_') )
        suites.append( unittest.makeSuite(test_angle,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
