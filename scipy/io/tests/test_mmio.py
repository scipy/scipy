#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from tempfile import mkdtemp, mktemp
import os
import shutil
from numpy import array,transpose
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal, \
            assert_equal, rand

import scipy.sparse
from scipy.io.mmio import mminfo,mmread,mmwrite


class TestMMIOArray(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.fn = os.path.join(self.tmpdir, 'testfile.mtx')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
    
    def test_all(self):
        def single_test(a, info):
            mmwrite(self.fn, a)
            assert_equal(mminfo(self.fn), info)
            b = mmread(self.fn)
            assert_array_almost_equal(a, b)
        
        # simple (integer)
        single_test([[1, 2], [3, 4]], (2, 2, 4, 'array', 'integer', 'general'))
        # simple upper triangle (integer)
        single_test([[0, 1], [0, 0]], (2, 2, 4, 'array', 'integer', 'general'))
        # simple lower triangle (integer)
        single_test([[0, 0], [1, 0]], (2, 2, 4, 'array', 'integer', 'general'))
        # simple rectangular (integer, 2x3)
        single_test([[1, 2, 3], [4, 5, 6]], (2, 3, 6, 'array', 'integer', 'general'))
        # simple rectangular (float, 3x2)
        single_test([[1, 2], [3.5, 4], [5, 6]], (3, 2, 6, 'array', 'real', 'general'))
        # simple (float)
        single_test([[1, 2], [3, 4.0]], (2, 2, 4, 'array', 'real', 'general'))
        # simple (complex)
        single_test([[1, 2], [3, 4j]], (2, 2, 4, 'array', 'complex', 'general'))
        # simple symmetric (integer)
        single_test([[1, 2], [2, 4]], (2, 2, 4, 'array', 'integer', 'symmetric'))
        # simple skew symmetric (integer)
        single_test([[1, 2], [-2, 4]], (2, 2, 4, 'array', 'integer', 'skew-symmetric'))
        # simple skew symmetric (float)
        single_test(array([[1,2], [-2.0,4]], 'f'), (2, 2, 4, 'array', 'real', 'skew-symmetric'))
        # simple hermitian (complex)
        single_test([[1, 2+3j], [2-3j, 4]], (2, 2, 4, 'array', 'complex', 'hermitian'))
        # random symmetric (float)
        sz = (20, 20)
        a = rand(*sz)
        a = a + transpose(a)
        single_test(a, (20, 20, 400, 'array', 'real', 'symmetric'))
        # random rectangular (float)
        sz = (20, 15)
        a = rand(*sz)
        single_test(a, (20, 15, 300, 'array', 'real', 'general'))


class TestMMIOSparseCSR(TestMMIOArray):

    def setUp(self):
        self.tmpdir = mkdtemp()
        self.fn = os.path.join(self.tmpdir, 'testfile.mtx')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
    
    def test_all(self):
        def single_test(a, info):
            mmwrite(self.fn, a)
            assert_equal(mminfo(self.fn), info)
            b = mmread(self.fn)
            assert_array_almost_equal(a.todense(), b.todense())
        
        # simple (integer)
        single_test(scipy.sparse.csr_matrix([[1, 2], [3, 4]]), (2, 2, 4, 'coordinate', 'integer', 'general'))
        # simple upper triangle (integer)
        single_test(scipy.sparse.csr_matrix([[0, 1], [0, 0]]), (2, 2, 1, 'coordinate', 'integer', 'general'))
        # simple lower triangle (integer)
        single_test(scipy.sparse.csr_matrix([[0, 0], [1, 0]]), (2, 2, 1, 'coordinate', 'integer', 'general'))
        # simple rectangular (integer, 2x3)
        single_test(scipy.sparse.csr_matrix([[1, 2, 3], [4, 5, 6]]), (2, 3, 6, 'coordinate', 'integer', 'general'))
        # simple rectangular (float, 3x2)
        single_test(scipy.sparse.csr_matrix([[1, 2], [3.5, 4], [5, 6]]), (3, 2, 6, 'coordinate', 'real', 'general'))
        # simple (float)
        single_test(scipy.sparse.csr_matrix([[1, 2], [3, 4.0]]), (2, 2, 4, 'coordinate', 'real', 'general'))
        # simple (complex)
        single_test(scipy.sparse.csr_matrix([[1, 2], [3, 4j]]), (2, 2, 4, 'coordinate', 'complex', 'general'))
        # simple symmetric (integer)
        single_test(scipy.sparse.csr_matrix([[1, 2], [2, 4]]), (2, 2, 3, 'coordinate', 'integer', 'symmetric'))
        # simple skew symmetric (integer)
        single_test(scipy.sparse.csr_matrix([[1, 2], [-2, 4]]), (2, 2, 3, 'coordinate', 'integer', 'skew-symmetric'))
        # simple skew symmetric (float)
        single_test(scipy.sparse.csr_matrix(array([[1,2], [-2.0,4]], 'f')), (2, 2, 3, 'coordinate', 'real', 'skew-symmetric'))
        # simple hermitian (complex)
        single_test(scipy.sparse.csr_matrix([[1, 2+3j], [2-3j, 4]]), (2, 2, 3, 'coordinate', 'complex', 'hermitian'))
        # random symmetric (float)
        sz = (20, 20)
        a = rand(*sz)
        a = a + transpose(a)
        a = scipy.sparse.csr_matrix(a)
        single_test(scipy.sparse.csr_matrix(a), (20, 20, 210, 'coordinate', 'real', 'symmetric'))
        # random rectangular (float)
        sz = (20, 15)
        a = rand(*sz)
        a = scipy.sparse.csr_matrix(a)
        single_test(scipy.sparse.csr_matrix(a), (20, 15, 300, 'coordinate', 'real', 'general'))


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
    
    def test_read(self):
        def single_test(example, a, info):
            f = open(self.fn, 'w')
            f.write(example)
            f.close()
            assert_equal(mminfo(self.fn), info)
            b = mmread(self.fn).todense()
            assert_array_almost_equal(a, b)
        
        # general
        a = [[1, 0, 0, 6, 0],
             [0, 10.5, 0, 0, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 33.32],
             [0, 0, 0, 0, 12]]
        single_test(_general_example, a, (5, 5, 8, 'coordinate', 'real', 'general'))
        # hermitian
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, 250.5 - 22.22j, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5 + 22.22j, 0, -280, -33.32j],
             [0, 0, 0, 33.32j, 12]]
        single_test(_hermitian_example, a, (5, 5, 7, 'coordinate', 'complex', 'hermitian'))
        # skew
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, -250.5, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 0],
             [0, 0, 0, 0, 12]]
        single_test(_skew_example, a, (5, 5, 7, 'coordinate', 'real', 'skew-symmetric'))
        # symmetric
        a = [[1, 0, 0, 0, 0],
             [0, 10.5, 0, 250.5, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 8],
             [0, 0, 0, 8, 12]]
        single_test(_symmetric_example, a, (5, 5, 7, 'coordinate', 'real', 'symmetric'))
        # symmetric pattern
        a = [[1, 0, 0, 6, 0],
             [0, 10.5, 0, 0, 0],
             [0, 0, .015, 0, 0],
             [0, 250.5, 0, -280, 33.32],
             [0, 0, 0, 0, 12]]
        single_test(_symmetric_pattern_example, a, (5, 5, 7, 'coordinate', 'pattern', 'symmetric'))

    def test_empty_write_read(self):
        #http://projects.scipy.org/scipy/ticket/883

        b = scipy.sparse.coo_matrix((10, 10))
        mmwrite(self.fn, b)

        assert_equal(mminfo(self.fn), (10, 10, 0, 'coordinate', 'real', 'symmetric'))
        a = b.todense()
        b = mmread(self.fn).todense()
        assert_array_almost_equal(a, b)

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

        b = scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5))

        mmwrite(self.fn, b)

        fn_bzip2 = "%s.bz2" % self.fn
        with open(self.fn, 'rb') as f_in:
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

        b = scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5))

        mmwrite(self.fn, b)

        fn_gzip = "%s.gz" % self.fn
        with open(self.fn, 'rb') as f_in:
            f_out = gzip.open(fn_gzip, 'wb')
            f_out.write(f_in.read())
            f_out.close()

        a = mmread(fn_gzip).todense()
        assert_array_almost_equal(a, b.todense())

    def test_real_write_read(self):
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])

        b = scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5))

        mmwrite(self.fn, b)

        assert_equal(mminfo(self.fn),(5,5,8,'coordinate','real','general'))
        a = b.todense()
        b = mmread(self.fn).todense()
        assert_array_almost_equal(a, b)

    def test_complex_write_read(self):
        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])
        V = array([1.0 + 3j, 6.0 + 2j, 10.50 + 0.9j, 0.015 + -4.4j,
                   250.5 + 0j, -280.0 + 5j, 33.32 + 6.4j, 12.00 + 0.8j])

        b = scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5))

        mmwrite(self.fn, b)

        assert_equal(mminfo(self.fn),(5, 5, 8, 'coordinate', 'complex', 'general'))
        a = b.todense()
        b = mmread(self.fn).todense()
        assert_array_almost_equal(a, b)

    def test_sparse_formats(self):
        mats = []

        I = array([0, 0, 1, 2, 3, 3, 3, 4])
        J = array([0, 3, 1, 2, 1, 3, 4, 4])

        V = array([1.0, 6.0, 10.5, 0.015, 250.5, -280.0, 33.32, 12.0])
        mats.append(scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5)))

        V = array([1.0 + 3j, 6.0 + 2j, 10.50 + 0.9j, 0.015 + -4.4j,
                   250.5 + 0j, -280.0 + 5j, 33.32 + 6.4j, 12.00 + 0.8j])
        mats.append(scipy.sparse.coo_matrix((V, (I, J)), shape=(5, 5)))

        for mat in mats:
            expected = mat.todense()
            for fmt in ['csr','csc','coo']:
                fn = mktemp(dir=self.tmpdir)  # safe, we own tmpdir
                mmwrite(fn, mat.asformat(fmt))

                result = mmread(fn).todense()
                assert_array_almost_equal(result, expected)


if __name__ == "__main__":
    run_module_suite()
