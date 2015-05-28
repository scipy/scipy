#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from tempfile import mkdtemp, mktemp
import os
import shutil
from numpy import array,transpose, pi
from numpy.testing import (TestCase, run_module_suite, assert_equal,
    assert_array_equal, assert_array_almost_equal, rand)

import scipy.sparse
from scipy.io.mmio import mminfo,mmread,mmwrite


class TestMMIOArray(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.fn = os.path.join(self.tmpdir, 'testfile.mtx')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_simple(self):
        a = [[1,2],[3,4]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_upper_triangle(self):
        a = [[0,1],[0,0]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_lower_triangle(self):
        a = [[0,0],[1,0]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_rectangular(self):
        a = [[1,2,3],[4,5,6]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,3,6,'array','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_rectangular_real(self):
        a = [[1,2],[3.5,4],[5,6]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(3,2,6,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_real(self):
        a = [[1,2],[3,4.0]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_complex(self):
        a = [[1,2],[3,4j]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_symmetric(self):
        a = [[1,2],[2,4]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_skew_symmetric(self):
        a = [[1,2],[-2,4]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','integer','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_skew_symmetric_float(self):
        a = array([[1,2],[-2.0,4]],'f')
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','real','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_simple_hermitian(self):
        a = [[1,2+3j],[2-3j,4]]
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'array','complex','hermitian'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_random_symmetric_real(self):
        sz = (20,20)
        a = rand(*sz)
        a = a + transpose(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,20,400,'array','real','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)

    def test_random_rect_real(self):
        sz = (20,15)
        a = rand(*sz)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,15,300,'array','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a,b)


class TestMMIOSparseCSR(TestMMIOArray):

    def setUp(self):
        self.tmpdir = mkdtemp()
        self.fn = os.path.join(self.tmpdir, 'testfile.mtx')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
    
    def test_simple(self):
        a = [[1,2],[3,4]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'coordinate','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_upper_triangle(self):
        a = [[0,1],[0,0]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,1,'coordinate','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_lower_triangle(self):
        a = [[0,0],[1,0]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,1,'coordinate','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_rectangular(self):
        a = [[1,2,3],[4,5,6]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,3,6,'coordinate','integer','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_rectangular_real(self):
        a = [[1,2],[3.5,4],[5,6]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(3,2,6,'coordinate','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_real(self):
        a = [[1,2],[3,4.0]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'coordinate','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_complex(self):
        a = [[1,2],[3,4j]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,4,'coordinate','complex','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_symmetric(self):
        a = [[1,2],[2,4]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,3,'coordinate','integer','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_skew_symmetric(self):
        a = [[1,2],[-2,4]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,3,'coordinate','integer','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_skew_symmetric_float(self):
        a = array([[1,2],[-2.0,4]],'f')
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,3,'coordinate','real','skew-symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_simple_hermitian(self):
        a = [[1,2+3j],[2-3j,4]]
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(2,2,3,'coordinate','complex','hermitian'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_random_symmetric_real(self):
        sz = (20,20)
        a = rand(*sz)
        a = a + transpose(a)
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,20,210,'coordinate','real','symmetric'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())

    def test_random_rect_real(self):
        sz = (20,15)
        a = rand(*sz)
        a = scipy.sparse.csr_matrix(a)
        fn = self.fn
        mmwrite(fn,a)
        assert_equal(mminfo(fn),(20,15,300,'coordinate','real','general'))
        b = mmread(fn)
        assert_array_almost_equal(a.todense(),b.todense())


_general_example = '''\
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

_hermitian_example = '''\
%%MatrixMarket matrix coordinate complex hermitian
  5  5  7
    1     1     1.0      0
    2     2    10.5      0
    4     2   250.5     22.22
    3     3     1.5e-2   0
    4     4    -2.8e2    0
    5     5    12.       0
    5     4     0       33.32
'''

_skew_example = '''\
%%MatrixMarket matrix coordinate real skew-symmetric
  5  5  7
    1     1     1.0
    2     2    10.5
    4     2   250.5
    3     3     1.5e-2
    4     4    -2.8e2
    5     5    12.
    5     4     0
'''

_symmetric_example = '''\
%%MatrixMarket matrix coordinate real symmetric
  5  5  7
    1     1     1.0
    2     2    10.5
    4     2   250.5
    3     3     1.5e-2
    4     4    -2.8e2
    5     5    12.
    5     4     8
'''

_symmetric_pattern_example = '''\
%%MatrixMarket matrix coordinate pattern symmetric
  5  5  7
    1     1
    2     2
    4     2
    3     3
    4     4
    5     5
    5     4
'''


class TestMMIOCoordinate(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.fn = os.path.join(self.tmpdir, 'testfile.mtx')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_read_general(self):
        fn = self.fn
        f = open(fn,'w')
        f.write(_general_example)
        f.close()
        assert_equal(mminfo(fn),(5,5,8,'coordinate','real','general'))
        a = [[1, 0, 0, 6, 0],
             [0, 10.5, 0, 0, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 33.32],
             [0, 0, 0, 0, 12]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_read_hermitian(self):
        fn = self.fn
        f = open(fn,'w')
        f.write(_hermitian_example)
        f.close()
        assert_equal(mminfo(fn),(5,5,7,'coordinate','complex','hermitian'))
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, 250.5 - 22.22j, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5 + 22.22j, 0, -280, -33.32j],
             [0, 0, 0, 33.32j, 12]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_read_skew(self):
        fn = self.fn
        f = open(fn,'w')
        f.write(_skew_example)
        f.close()
        assert_equal(mminfo(fn),(5,5,7,'coordinate','real','skew-symmetric'))
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, -250.5, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 0],
             [0, 0, 0, 0, 12]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_read_symmetric(self):
        fn = self.fn
        f = open(fn,'w')
        f.write(_symmetric_example)
        f.close()
        assert_equal(mminfo(fn),(5,5,7,'coordinate','real','symmetric'))
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, 250.5, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 8],
             [0, 0, 0, 8, 12]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_read_symmetric_pattern(self):
        fn = self.fn
        f = open(fn,'w')
        f.write(_symmetric_pattern_example)
        f.close()
        assert_equal(mminfo(fn),(5,5,7,'coordinate','pattern','symmetric'))
        a = [[1, 0, 0, 0, 0],
             [0, 1, 0, 1, 0],
             [0, 0, 1, 0, 0],
             [0, 1, 0, 1, 1],
             [0, 0, 0, 1, 1]]
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_empty_write_read(self):
        #http://projects.scipy.org/scipy/ticket/883

        b = scipy.sparse.coo_matrix((10,10))
        fn = self.fn
        mmwrite(fn,b)

        assert_equal(mminfo(fn),(10,10,0,'coordinate','real','symmetric'))
        a = b.todense()
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_bzip2_py3(self):
        # test if fix for #2152 works
        try:
            # bz2 module isn't always built when building Python.
            import bz2
        except:
            return
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])

        b = scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5))

        fn = self.fn
        mmwrite(fn, b)

        fn_bzip2 = "%s.bz2" % fn
        with open(fn, 'rb') as f_in:
            f_out = bz2.BZ2File(fn_bzip2, 'wb')
            f_out.write(f_in.read())
            f_out.close()

        a = mmread(fn_bzip2).todense()
        assert_array_almost_equal(a, b.todense())

    def test_gzip_py3(self):
        # test if fix for #2152 works
        try:
            # gzip module can be missing from Python installation
            import gzip
        except:
            return
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])

        b = scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5))

        fn = self.fn
        mmwrite(fn, b)

        fn_gzip = "%s.gz" % fn
        with open(fn, 'rb') as f_in:
            f_out = gzip.open(fn_gzip, 'wb')
            f_out.write(f_in.read())
            f_out.close()

        a = mmread(fn_gzip).todense()
        assert_array_almost_equal(a, b.todense())

    def test_real_write_read(self):
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])

        b = scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5))

        fn = self.fn
        mmwrite(fn,b)

        assert_equal(mminfo(fn),(5,5,8,'coordinate','real','general'))
        a = b.todense()
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_complex_write_read(self):
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0 + 3j, 6.0 + 2j, 10.50 + 0.9j, 0.015 + -4.4j,
                   250.5 + 0j, -280.0 + 5j, 33.32 + 6.4j, 12.00 + 0.8j])

        b = scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5))

        fn = self.fn
        mmwrite(fn,b)

        assert_equal(mminfo(fn),(5,5,8,'coordinate','complex','general'))
        a = b.todense()
        b = mmread(fn).todense()
        assert_array_almost_equal(a,b)

    def test_sparse_formats(self):
        mats = []

        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])

        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])
        mats.append(scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5)))

        V = array([1.0 + 3j, 6.0 + 2j, 10.50 + 0.9j, 0.015 + -4.4j,
                   250.5 + 0j, -280.0 + 5j, 33.32 + 6.4j, 12.00 + 0.8j])
        mats.append(scipy.sparse.coo_matrix((V,(I,J)),shape=(5,5)))

        for mat in mats:
            expected = mat.todense()
            for fmt in ['csr','csc','coo']:
                fn = mktemp(dir=self.tmpdir)  # safe, we own tmpdir
                mmwrite(fn, mat.asformat(fmt))

                result = mmread(fn).todense()
                assert_array_almost_equal(result, expected)

    def test_precision(self):
        test_values = [pi] + [10**(i) for i in range(0, -10, -1)]
        test_precisions = range(1, 10)
        for value in test_values:
            for precision in test_precisions:
                # construct sparse matrix with test value at last main diagonal
                n = 10**precision + 1
                A = scipy.sparse.dok_matrix((n, n))
                A[n-1, n-1] = value
                # write matrix with test precision and read again
                mmwrite(self.fn, A, precision=precision)
                A = scipy.io.mmread(self.fn)
                # check for right entries in matrix
                assert_array_equal(A.row, [n-1])
                assert_array_equal(A.col, [n-1])
                assert_array_almost_equal(A.data,
                    [float('%%.%dg' % precision % value)])


if __name__ == "__main__":
    run_module_suite()
