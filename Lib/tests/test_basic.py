""" Test functions for basic module

"""

import unittest
import scipy.limits as limits
from scipy.scipy_test import assert_array_equal, assert_equal
from scipy.scipy_test import assert_almost_equal
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
            
        
##################################################

def test_suite():
    suites = []
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
    
    
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test():
    all_tests = test_suite()
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
