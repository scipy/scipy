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
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_rectangular(self):
        a = [[1,2,3],[4,5,6]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,3,6,'array','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_rectangular_real(self):
        a = [[1,2],[3.5,4],[5,6]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(3,2,6,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_real(self):
        a = [[1,2],[3,4.0]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_complex(self):
        a = [[1,2],[3,4j]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_symmetric(self):
        a = [[1,2],[2,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_skew_symmetric(self):
        a = [[1,2],[-2,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_skew_symmetric_float(self):
        a = array([[1,2],[-2.0,4]],'f')
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_simple_hermitian(self):
        a = [[1,2+3j],[2-3j,4]]
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','hermitian'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_random_symmetric_real(self):
        sz = (20,20)
        a = rand(*sz)
        a = a + transpose(a)
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,20,400,'array','real','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def check_random_rect_real(self):
        sz = (20,15)
        a = rand(*sz)
        fn = mktemp()
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,15,300,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

_exmpl_mtx = '''\
%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This ASCII file represents a sparse MxN matrix with L 
% nonzeros in the following Matrix Market format:
%
% +----------------------------------------------+
% |%%MatrixMarket matrix coordinate real general | <--- header line
% |%                                             | <--+
% |% comments                                    |    |-- 0 or more comment lines
% |%                                             | <--+         
% |    M  N  L                                   | <--- rows, columns, entries
% |    I1  J1  A(I1, J1)                         | <--+
% |    I2  J2  A(I2, J2)                         |    |
% |    I3  J3  A(I3, J3)                         |    |-- L lines
% |        . . .                                 |    |
% |    IL JL  A(IL, JL)                          | <--+
% +----------------------------------------------+   
%
% Indices are 1-based, i.e. A(1,1) is the first element.
%
%=================================================================================
  5  5  8
    1     1   1.000e+00
    2     2   1.050e+01
    3     3   1.500e-02
    1     4   6.000e+00
    4     2   2.505e+02
    4     4  -2.800e+02
    4     5   3.332e+01
    5     5   1.200e+01
'''

class test_mmio_coordinate(ScipyTestCase):

    def check_simple_todense(self):
        fn = mktemp()
        f = open(fn,'w')
        f.write(_exmpl_mtx)
        f.close()
        assert_equal(mminfo(fn),(5,5,8,'coordinate','real','general'))
        a = [[1,    0,      0,       6,      0],
             [0,   10.5,    0,       0,      0],
             [0,    0,    .015,      0,      0],
             [0,  250.5,    0,     -280,    33.32],
             [0,    0,      0,       0,     12]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)
        
if __name__ == "__main__":
    ScipyTest('io.mmio').run()

