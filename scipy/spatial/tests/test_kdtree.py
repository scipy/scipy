# Copyright Anne M. Archibald 2008
# Released under the scipy license
from numpy.testing import *

import numpy as np
from scipy.spatial import KDTree, distance

class CheckSmall(NumpyTestCase):
    def setUp(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = KDTree(self.data)

    def test_nearest(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 1),
                (0.1,0))
    def test_nearest_two(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 2),
                ([0.1,0.9],[0,1]))
class CheckSmallNonLeaf(NumpyTestCase):
    def setUp(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = KDTree(self.data,leafsize=1)

    def test_nearest(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 1),
                (0.1,0))
    def test_nearest_two(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 2),
                ([0.1,0.9],[0,1]))

class CheckRandom(NumpyTestCase):
    def setUp(self):
        self.n = 1000
        self.k = 4
        self.data = np.random.randn(self.n, self.k)
        self.kdtree = KDTree(self.data)

    def test_nearest(self):
        x = np.random.randn(self.k)
        d, i = self.kdtree.query(x, 1)
        assert_almost_equal(d**2,np.sum((x-self.data[i])**2))
        eps = 1e-8
        assert np.all(np.sum((self.data-x[np.newaxis,:])**2,axis=1)>d**2-eps)
        
    def test_m_nearest(self):
        x = np.random.randn(self.k)
        m = 10
        dd, ii = self.kdtree.query(x, m)
        d = np.amax(dd)
        i = ii[np.argmax(dd)]
        assert_almost_equal(d**2,np.sum((x-self.data[i])**2))
        eps = 1e-8
        assert_equal(np.sum(np.sum((self.data-x[np.newaxis,:])**2,axis=1)<d**2+eps),m)

    def test_points_near(self):
        x = np.random.randn(self.k)
        d = 0.2
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d==np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d**2,np.sum((x-self.data[near_i])**2))
            assert near_d<d+eps, "near_d=%g should be less than %g" % (near_d,d)
        assert_equal(np.sum(np.sum((self.data-x[np.newaxis,:])**2,axis=1)<d**2+eps),hits)

    def test_points_near_l1(self):
        x = np.random.randn(self.k)
        d = 0.2
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=1, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d==np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d,distance(x,self.data[near_i],1))
            assert near_d<d+eps, "near_d=%g should be less than %g" % (near_d,d)
        assert_equal(np.sum(distance(self.data,x,1)<d+eps),hits)
    def test_points_near_linf(self):
        x = np.random.randn(self.k)
        d = 0.2
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=np.inf, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d==np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d,distance(x,self.data[near_i],np.inf))
            assert near_d<d+eps, "near_d=%g should be less than %g" % (near_d,d)
        assert_equal(np.sum(distance(self.data,x,np.inf)<d+eps),hits)

    def test_approx(self):
        x = np.random.randn(self.k)
        m = 10
        eps = 0.1
        d_real, i_real = self.kdtree.query(x, m)
        d, i = self.kdtree.query(x, m, eps=eps)
        assert np.all(d<=d_real*(1+eps))

class CheckVectorization(NumpyTestCase):
    def setUp(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = KDTree(self.data)

    def test_single_query(self):
        d, i = self.kdtree.query([0,0,0])
        assert isinstance(d,float)
        assert isinstance(i,int)

    def test_vectorized_query(self):
        d, i = self.kdtree.query(np.zeros((2,4,3)))
        assert_equal(np.shape(d),(2,4))
        assert_equal(np.shape(i),(2,4))

    def test_single_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query([0,0,0],k=kk)
        assert_equal(np.shape(d),(kk,))
        assert_equal(np.shape(i),(kk,))
        assert np.all(~np.isfinite(d[-s:]))
        assert np.all(i[-s:]==self.kdtree.n)
    def test_vectorized_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query(np.zeros((2,4,3)),k=kk)
        assert_equal(np.shape(d),(2,4,kk))
        assert_equal(np.shape(i),(2,4,kk))
        assert np.all(~np.isfinite(d[:,:,-s:]))
        assert np.all(i[:,:,-s:]==self.kdtree.n)
    def test_single_query_all_neighbors(self):
        d, i = self.kdtree.query([0,0,0],k=None,distance_upper_bound=1.1)
        assert isinstance(d,list)
        assert isinstance(i,list)
    def test_vectorized_query_all_neighbors(self):
        d, i = self.kdtree.query(np.zeros((2,4,3)),k=None,distance_upper_bound=1.1)
        assert_equal(np.shape(d),(2,4))
        assert_equal(np.shape(i),(2,4))

        assert isinstance(d[0,0],list)
        assert isinstance(i[0,0],list)



    
if __name__=='__main__':
    import unittest
    unittest.main()


