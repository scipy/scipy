#!/usr/bin/env python

import sys
from tempfile import mktemp
from scipy_test.testing import *
from scipy_base import array,transpose

set_package_path()
from io.mmio import mminfo,mmread,mmwrite
restore_path()

class test_mmio_array(ScipyTestCase):

    def check_simple(self):
        a = [[1,2],[3,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_rectangular(self):
        a = [[1,2,3],[4,5,6]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,3,6,'array','integer','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,3,6))
        assert_array_almost_equal(a,b)

    def check_simple_rectangular_real(self):
        a = [[1,2],[3.5,4],[5,6]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(3,2,6,'array','real','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(3,2,6))
        assert_array_almost_equal(a,b)

    def check_simple_real(self):
        a = [[1,2],[3,4.0]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_complex(self):
        a = [[1,2],[3,4j]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_symmetric(self):
        a = [[1,2],[2,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','symmetric'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_skew_symmetric(self):
        a = [[1,2],[-2,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','skew-symmetric'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_skew_symmetric_float(self):
        a = array([[1,2],[-2.0,4]],'f')
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','skew-symmetric'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_simple_hermitian(self):
        a = [[1,2+3j],[2-3j,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','hermitian'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(2,2,4))
        assert_array_almost_equal(a,b)

    def check_random_symmetric_real(self):
        sz = (20,20)
        a = rand(*sz)
        a = a + transpose(a)
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,20,400,'array','real','symmetric'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(20,20,400))
        assert_array_almost_equal(a,b)

    def check_random_rect_real(self):
        sz = (20,15)
        a = rand(*sz)
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,15,300,'array','real','general'))
        b,m,n,s = mmread(fn)
        assert_equal((m,n,s),(20,15,300))
        assert_array_almost_equal(a,b)

if __name__ == "__main__":
    ScipyTest('io.mmio').run()

