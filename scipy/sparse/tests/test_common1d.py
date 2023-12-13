"""Test of 1D aspects of sparse array classes"""

import contextlib
import operator
import pytest

import numpy as np
from numpy import array

from numpy.testing import (assert_equal, assert_array_equal,
        assert_array_almost_equal, assert_almost_equal, assert_,
        assert_allclose,suppress_warnings)

from scipy.sparse import (coo_array, csr_array, csc_array,
                          dok_array, sparray, issparse,
                          SparseEfficiencyWarning)
from scipy.sparse._sputils import supported_dtypes, isscalarlike, matrix


sup_complex = suppress_warnings()
sup_complex.filter(np.ComplexWarning)


def assert_array_equal_dtype(x, y, **kwargs):
    assert_(x.dtype == y.dtype)
    assert_array_equal(x, y, **kwargs)


def toarray(a):
    if isinstance(a, np.ndarray) or isscalarlike(a):
        return a
    return a.toarray()


class _Common1D:
    """test common functionality shared by 1D sparse formats"""

    # Canonical data.
    dat1d = np.array([3, 0, 1, 0], 'd')

    # Some sparse and dense matrices with data for every supported dtype.
    # This set union is a workaround for numpy#6295, which means that
    # two np.int64 dtypes don't hash to the same value.
    math_dtypes = [np.int_, np.float_, np.complex_]
    checked_dtypes = set(supported_dtypes).union(math_dtypes)
    dat_dtypes = {}
    for dtype in checked_dtypes:
        dat_dtypes[dtype] = dat1d.astype(dtype)

    def test_class_vars(self):
        # Check that the original data is equivalent to the
        # corresponding dat_dtypes & datsp_dtypes.
        assert_equal(self.dat1d, self.dat_dtypes[np.float64])
        assert_equal(self.datsp.toarray(),
                     self.datsp_dtypes[np.float64].toarray())

    def test_empty(self):
        # create empty matrices
        assert_equal(self.spcreator((3,)).toarray(), np.zeros(3))
        assert_equal(self.spcreator((3,)).nnz, 0)
        assert_equal(self.spcreator((3,)).count_nonzero(), 0)

    def test_invalid_shapes(self):
        pytest.raises(ValueError, self.spcreator, (-3,))

    def test_repr(self):
        repr(self.datsp)

    def test_str(self):
        str(self.datsp)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_empty_arithmetic(self):
        # Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        shape = (5,)
        for mytype in [np.dtype('int32'), np.dtype('float32'),
                np.dtype('float64'), np.dtype('complex64'),
                np.dtype('complex128')]:
            a = self.spcreator(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a @ a.tocsc()
            e = a @ a.tocsr()
            f = a @ a.tocoo()
            for m in [a,b,c,d,e,f]:
                assert_equal(toarray(m), a.toarray()@a.toarray())
                assert_equal(m.dtype, mytype)
                assert_equal(toarray(m).dtype, mytype)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_abs(self):
        A = np.array([-1, 0, 17, 0, -5, 0, 1, -4, 0, 0, 0, 0], 'd')
        assert_equal(abs(A), abs(self.spcreator(A)).toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_round(self):
        decimal = 1
        A = np.array([-1.35, 0.56, 17.25, -5.98], 'd')
        assert_equal(np.around(A, decimals=decimal),
                     round(self.spcreator(A), ndigits=decimal).toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_elementwise_power(self):
        A = np.array([-4, -3, -2, -1, 0, 1, 2, 3, 4], 'd')
        assert_equal(np.power(A, 2), self.spcreator(A).power(2).toarray())

        #it's element-wise power function, input has to be a scalar
        pytest.raises(NotImplementedError, self.spcreator(A).power, A)

    def test_neg(self):
        A = np.array([-1, 0, 17, 0, -5, 0, 1, -4, 0, 0, 0, 0], 'd')
        assert_equal(-A, (-self.spcreator(A)).toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_real(self):
        D = np.array([1 + 3j, 2 - 4j])
        A = self.spcreator(D)
        assert_equal(A.real.toarray(), D.real)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_imag(self):
        D = np.array([1 + 3j, 2 - 4j])
        A = self.spcreator(D)
        assert_equal(A.imag.toarray(), D.imag)

    def test_reshape_1d_tofrom_row_or_column(self):
        # add a dimension
        x = self.spcreator([1, 0, 7, 0, 0, 0, 0, -3, 0, 0, 0, 5])
        y = x.reshape(1, 12)
        desired = [[1, 0, 7, 0, 0, 0, 0, -3, 0, 0, 0, 5]]
        assert_array_equal(y.toarray(), desired)

        # remove a size-1 dimension
        x = self.spcreator(desired)
        y = x.reshape(12)
        assert_array_equal(y.toarray(), desired[0])
        y2 = x.reshape((12,))
        assert_equal(y.shape, y2.shape)

        # make a column 1d
        y = x.T.reshape(12)
        assert_array_equal(y.toarray(), desired[0])

    def test_reshape(self):
        x = self.spcreator([1, 0, 7, 0, 0, 0, 0, -3, 0, 0, 0, 5])
        y = x.reshape((4, 3))
        desired = [[1, 0, 7], [0, 0, 0], [0, -3, 0], [0, 0, 5]]
        assert_array_equal(y.toarray(), desired)

        y = x.reshape((12,))
        assert_(y is x)

        y = x.reshape(12)
        assert_array_equal(y.toarray(), x.toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_getrowcol(self):
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            assert_array_equal(self.datsp.getrow(0).toarray(), self.dat1d[None, :])
            assert_array_equal(self.datsp.getrow(-1).toarray(), self.dat1d[None, :])
            assert_array_equal(self.datsp.getcol(2).toarray(), self.dat1d[None, [2]])
            assert_array_equal(self.datsp.getcol(-2).toarray(), self.dat1d[None, [-2]])

    def test_sum(self):
        np.random.seed(1234)
        dat_1 = np.array([0, 1, 2, 3, -4, 5, -6, 7, 9])
        dat_2 = np.random.rand(5)
        dat_3 = np.array([])
        dat_4 = np.zeros((40, ))
        matrices = [dat_1, dat_2, dat_3, dat_4]

        for m in matrices:
            dat = np.array(m)
            datsp = self.spcreator(dat)
            with np.errstate(over='ignore'):
                assert_(np.isscalar(datsp.sum()))
                assert_array_almost_equal(dat.sum(), datsp.sum())
                assert_array_almost_equal(dat.sum(axis=None), datsp.sum(axis=None))
                assert_array_almost_equal(dat.sum(axis=0), datsp.sum(axis=0))
                assert_array_almost_equal(dat.sum(axis=-1), datsp.sum(axis=-1))

        # test out parameter
        datsp.sum(axis=0, out=np.zeros(()))

    def test_sum_invalid_params(self):
        out = np.zeros((3,))  # wrong size for out
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        pytest.raises(ValueError, datsp.sum, axis=1)
        pytest.raises(TypeError, datsp.sum, axis=(0, 1))
        pytest.raises(TypeError, datsp.sum, axis=1.5)
        pytest.raises(ValueError, datsp.sum, axis=0, out=out)

    def test_numpy_sum(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        dat_sum = np.sum(dat)
        datsp_sum = np.sum(datsp)

        assert_array_almost_equal(dat_sum, datsp_sum)

    def test_mean(self):
        keepdims = not isinstance(self, sparray)
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        assert_array_almost_equal(dat.mean(), datsp.mean())
        assert_(np.isscalar(datsp.mean(axis=None)))
        assert_array_almost_equal(
            dat.mean(axis=None, keepdims=keepdims), datsp.mean(axis=None)
        )
        assert_array_almost_equal(
            dat.mean(axis=0, keepdims=keepdims), datsp.mean(axis=0)
        )
        assert_array_almost_equal(
            dat.mean(axis=-1, keepdims=keepdims), datsp.mean(axis=-1)
        )

        with pytest.raises(ValueError, match='axis'):
            datsp.mean(axis=1)
        with pytest.raises(ValueError, match='axis'):
            datsp.mean(axis=-2)

    def test_mean_invalid_params(self):
        out = np.asarray(np.zeros((1, 3)))
        dat = np.array([[0, 1, 2],
                     [3, -4, 5],
                     [-6, 7, 9]])

        if self.spcreator._format == 'uni':
            pytest.raises(ValueError, self.spcreator, dat)
            return

        datsp = self.spcreator(dat)
        pytest.raises(ValueError, datsp.mean, axis=3)
        pytest.raises(TypeError, datsp.mean, axis=(0, 1))
        pytest.raises(TypeError, datsp.mean, axis=1.5)
        pytest.raises(ValueError, datsp.mean, axis=1, out=out)

    def test_sum_dtype(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        for dtype in self.checked_dtypes:
            dat_sum = dat.sum(dtype=dtype)
            datsp_sum = datsp.sum(dtype=dtype)

            assert_array_almost_equal(dat_sum, datsp_sum)
            assert_equal(dat_sum.dtype, datsp_sum.dtype)

    def test_mean_dtype(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        for dtype in self.checked_dtypes:
            dat_mean = dat.mean(dtype=dtype)
            datsp_mean = datsp.mean(dtype=dtype)

            assert_array_almost_equal(dat_mean, datsp_mean)
            assert_equal(dat_mean.dtype, datsp_mean.dtype)

    def test_mean_out(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        dat_out = np.array([0])
        datsp_out = np.array([0])

        dat.mean(out=dat_out, keepdims=True)
        datsp.mean(out=datsp_out)
        assert_array_almost_equal(dat_out, datsp_out)

        dat.mean(axis=0, out=dat_out, keepdims=True)
        datsp.mean(axis=0, out=datsp_out)
        assert_array_almost_equal(dat_out, datsp_out)

    def test_numpy_mean(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)

        dat_mean = np.mean(dat)
        datsp_mean = np.mean(datsp)

        assert_array_almost_equal(dat_mean, datsp_mean)
        assert_equal(dat_mean.dtype, datsp_mean.dtype)

    @sup_complex
    def test_from_array(self):
        A = np.array([2,3,4])
        assert_array_equal(self.spcreator(A).toarray(), A)

        A = np.array([1.0 + 3j, 0, -1])
        assert_array_equal(self.spcreator(A).toarray(), A)
        assert_array_equal(self.spcreator(A, dtype='int16').toarray(),A.astype('int16'))

    @sup_complex
    def test_from_list(self):
        A = [2,3,4]
        assert_array_equal(self.spcreator(A).toarray(), A)

        A = [1.0 + 3j, 0, -1]
        assert_array_equal(self.spcreator(A).toarray(), np.array(A))
        assert_array_equal(
            self.spcreator(A, dtype='int16').toarray(), np.array(A).astype('int16')
        )

    @sup_complex
    def test_from_sparse(self):
        D = np.array([1,0,0])
        S = coo_array(D)
        assert_array_equal(self.spcreator(S).toarray(), D)
        S = self.spcreator(D)
        assert_array_equal(self.spcreator(S).toarray(), D)

        D = np.array([1.0 + 3j, 0, -1])
        S = coo_array(D)
        assert_array_equal(self.spcreator(S).toarray(), D)
        assert_array_equal(
            self.spcreator(S, dtype='int16').toarray(), D.astype('int16')
        )
        S = self.spcreator(D)
        assert_array_equal(self.spcreator(S).toarray(), D)
        assert_array_equal(
            self.spcreator(S, dtype='int16').toarray(), D.astype('int16')
        )

    def test_toarray(self):
        # Check C- or F-contiguous (default).
        dat = np.asarray(self.dat1d)
        chk = self.datsp.toarray()
        assert_array_equal(chk, dat)
        assert_(chk.flags.c_contiguous == chk.flags.f_contiguous)
        # Check C-contiguous (with arg).
        chk = self.datsp.toarray(order='C')
        assert_array_equal(chk, dat)
        assert_(chk.flags.c_contiguous)
        assert_(chk.flags.f_contiguous)
        # Check F-contiguous (with arg).
        chk = self.datsp.toarray(order='F')
        assert_array_equal(chk, dat)
        assert_(chk.flags.c_contiguous)
        assert_(chk.flags.f_contiguous)

        # Check with output arg.
        out = np.zeros(self.datsp.shape, dtype=self.datsp.dtype)
        self.datsp.toarray(out=out)
        assert_array_equal(chk, dat)
        # Check that things are fine when we don't initialize with zeros.
        out[...] = 1.
        self.datsp.toarray(out=out)
        assert_array_equal(chk, dat)
        a = np.array([1., 2., 3., 4.])
        dense_dot_dense = np.dot(a, dat)
        check = np.dot(a, self.datsp.toarray())
        assert_array_equal(dense_dot_dense, check)
        b = np.array([1., 2., 3., 4.])
        dense_dot_dense = np.dot(dat, b)
        check2 = np.dot(self.datsp.toarray(), b)
        assert_array_equal(dense_dot_dense, check2)

        # Check bool data works.
        spbool = self.spcreator(dat, dtype=bool)
        arrbool = dat.astype(bool)
        assert_array_equal(spbool.toarray(), arrbool)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_mul_scalar(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            assert_array_equal(dat*2, (datsp*2).toarray())
            assert_array_equal(dat*17.3, (datsp*17.3).toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_rmul_scalar(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            assert_array_equal(2*dat, (2*datsp).toarray())
            assert_array_equal(17.3*dat, (17.3*datsp).toarray())

    def test_add(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            a = dat.copy()
            a[0] = 2.0
            b = datsp
            c = b + a
            assert_array_equal(c, b.toarray() + a)

            # test broadcasting
            # Note: cdnt add nonzero scalar to sparray. Can add len 1 array
            c = b + a[0:1]
            assert_array_equal(c, b.toarray() + a[0])

    def test_radd(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            a = dat.copy()
            a[0] = 2.0
            b = datsp
            c = a + b
            assert_array_equal(c, a + b.toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_sub(self):
        for dtype in self.math_dtypes:
            if dtype == np.dtype('bool'):
                # boolean array subtraction deprecated in 1.9.0
                continue

            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            assert_array_equal((datsp - datsp).toarray(), np.zeros(4))
            assert_array_equal((datsp - 0).toarray(), dat)

            A = self.spcreator([1, -4, 0, 2], dtype='d')
            assert_array_equal((datsp - A).toarray(), dat - A.toarray())
            assert_array_equal((A - datsp).toarray(), A.toarray() - dat)

            # test broadcasting
            assert_array_equal(datsp.toarray() - dat[0], dat - dat[0])

    def test_rsub(self):
        for dtype in self.math_dtypes:
            if dtype == np.dtype('bool'):
                # boolean array subtraction deprecated in 1.9.0
                continue
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            assert_array_equal((dat - datsp),[0,0,0,0])
            assert_array_equal((datsp - dat),[0,0,0,0])
            assert_array_equal((0 - datsp).toarray(), -dat)

            A = self.spcreator([1, -4, 0, 2], dtype='d')
            assert_array_equal((dat - A), dat - A.toarray())
            assert_array_equal((A - dat), A.toarray() - dat)
            assert_array_equal(A.toarray() - datsp, A.toarray() - dat)
            assert_array_equal(datsp - A.toarray(), dat - A.toarray())

            # test broadcasting
            assert_array_equal(dat[:1] - datsp, dat[:1] - dat)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_add0(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            # Adding 0 to a sparse matrix
            assert_array_equal((datsp + 0).toarray(), dat)
            # use sum (which takes 0 as a starting value)
            sumS = sum([k * datsp for k in range(1, 3)])
            sumD = sum([k * dat for k in range(1, 3)])
            assert_almost_equal(sumS.toarray(), sumD)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_elementwise_multiply(self):
        # real/real
        A = np.array([4,0,9])
        B = np.array([0,7,-1])
        Asp = self.spcreator(A)
        Bsp = self.spcreator(B)
        assert_almost_equal(Asp.multiply(Bsp).toarray(), A*B)  # sparse/sparse
        assert_almost_equal(Asp.multiply(B).toarray(), A*B)  # sparse/dense

        # complex/complex
        C = np.array([1-2j,0+5j,-1+0j])
        D = np.array([5+2j,7-3j,-2+1j])
        Csp = self.spcreator(C)
        Dsp = self.spcreator(D)
        assert_almost_equal(Csp.multiply(Dsp).toarray(), C*D)  # sparse/sparse
        assert_almost_equal(Csp.multiply(D).toarray(), C*D)  # sparse/dense

        # real/complex
        assert_almost_equal(Asp.multiply(Dsp).toarray(), A*D)  # sparse/sparse
        assert_almost_equal(Asp.multiply(D).toarray(), A*D)  # sparse/dense

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_elementwise_multiply_broadcast(self):
        A = np.array([4])
        B = np.array([[-9]])
        C = np.array([1,-1,0])
        D = np.array([[7,9,-9]])
        E = np.array([[3],[2],[1]])
        F = np.array([[8,6,3],[-4,3,2],[6,6,6]])
        G = [1, 2, 3]
        H = np.ones((3, 4))
        J = H.T
        K = np.array([[0]])
        L = np.array([[[1,2],[0,1]]])

        # Some arrays can't be cast as spmatrices (A,C,L) so leave
        # them out.
        Asp = self.spcreator(A)
        Csp = self.spcreator(C)
        Gsp = self.spcreator(G)
        # 2d arrays
        creator = self.spcreator if Asp.format != 'uni' else csr_array
        Bsp = creator(B)
        Dsp = creator(D)
        Esp = creator(E)
        Fsp = creator(F)
        Hsp = creator(H)
        Hspp = creator(H[0,None])
        Jsp = creator(J)
        Jspp = creator(J[:,0,None])
        Ksp = creator(K)

        matrices = [A, B, C, D, E, F, G, H, J, K, L]
        spmatrices = [Asp, Bsp, Csp, Dsp, Esp, Fsp, Gsp, Hsp, Hspp, Jsp, Jspp, Ksp]
        sp1dmatrices = [Asp, Csp, Gsp]

        # sparse/sparse
        for i in sp1dmatrices:
            for j in spmatrices:
                try:
                    dense_mult = i.toarray() * j.toarray()
                except ValueError:
                    pytest.raises(ValueError, i.multiply, j)
                    continue
                sp_mult = i.multiply(j)
                assert_almost_equal(sp_mult.toarray(), dense_mult)

        # sparse/dense
        for i in sp1dmatrices:
            for j in matrices:
                try:
                    dense_mult = i.toarray() * j
                except TypeError:
                    continue
                except ValueError:
                    pytest.raises(ValueError, i.multiply, j)
                    continue
                sp_mult = i.multiply(j)
                if issparse(sp_mult):
                    assert_almost_equal(sp_mult.toarray(), dense_mult)
                else:
                    assert_almost_equal(sp_mult, dense_mult)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_elementwise_divide(self):
        expected = [1, np.nan, 1, np.nan]
        assert_array_equal(self.datsp / self.datsp, expected)

        denom = self.spcreator([1, 0, 0, 4],dtype='d')
        expected = [3, np.nan, np.inf, 0]
        assert_array_equal(toarray(self.datsp / denom), expected)

        # complex
        A = np.array([1-2j, 0+5j, -1+0j])
        B = np.array([5+2j, 7-3j, -2+1j])
        Asp = self.spcreator(A)
        Bsp = self.spcreator(B)
        assert_almost_equal(toarray(Asp / Bsp), A/B)

        # integer
        A = np.array([1, 2, 3])
        B = np.array([0, 1, 2])
        Asp = self.spcreator(A)
        Bsp = self.spcreator(B)
        with np.errstate(divide='ignore'):
            assert_array_equal(toarray(Asp / Bsp), A / B)

        # mismatching sparsity patterns
        A = np.array([0, 1])
        B = np.array([1, 0])
        Asp = self.spcreator(A)
        Bsp = self.spcreator(B)
        with np.errstate(divide='ignore', invalid='ignore'):
            assert_array_equal(array(toarray(Asp / Bsp)), A / B)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_pow(self):
        A = np.array([1, 0, 2, 0])
        B = self.spcreator(A)

        # unusual exponents
        with pytest.raises(ValueError, match='negative integer powers'):
            B ** -1
        with pytest.raises(NotImplementedError, match='zero power'):
            B ** 0

        for exponent in [1, 2, 3, 2.2]:
            ret_sp = B**exponent
            ret_np = A**exponent
            assert_array_equal(ret_sp.toarray(), ret_np)
            assert_equal(ret_sp.dtype, ret_np.dtype)

    def test_rmatvec(self):
        M = self.datsp
        assert_array_almost_equal([1,2,3,4] @ M, np.dot([1,2,3,4], M.toarray()))
        row = np.array([[1,2,3,4]])
        assert_array_almost_equal(row @ M, row @ M.toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_dot_scalar(self):
        M = self.datsp
        scalar = 10
        actual = M.dot(scalar)
        expected = M * scalar

        assert_allclose(actual.toarray(), expected.toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_matmul(self):
        M = self.spcreator([2,0,3.0])
        creator = self.spcreator if M.format != 'uni' else csr_array
        B = creator(array([[0,1],[1,0],[0,2]],'d'))
        col = np.array([[1,2,3]]).T

        matmul = operator.matmul
        # check matrix-vector
        assert_array_almost_equal(matmul(M, col), M.toarray() @ col)

        # check matrix-matrix
        assert_array_almost_equal(matmul(M, B).toarray(), (M @ B).toarray())
        assert_array_almost_equal(matmul(M.toarray(), B), (M @ B).toarray())
        assert_array_almost_equal(matmul(M, B.toarray()), (M @ B).toarray())

        # check error on matrix-scalar
        pytest.raises(ValueError, matmul, M, 1)
        pytest.raises(ValueError, matmul, 1, M)

    def test_matvec(self):
        A = np.array([2,0,3.0])
        Asp = self.spcreator(A)
        col = np.array([[1,2,3]]).T

        assert_array_almost_equal(Asp @ col, Asp.toarray() @ col)

        assert_equal((A @ np.array([1,2,3])).shape,())
        assert Asp @ np.array([1,2,3]) == 11
        assert_equal((Asp @ np.array([1,2,3])).shape,())
        assert_equal((Asp @ np.array([[1],[2],[3]])).shape,())
        # check result type
        assert isinstance(Asp @ matrix([[1,2,3]]).T, np.ndarray)

        # ensure exception is raised for improper dimensions
        bad_vecs = [array([1,2]), np.array([1,2,3,4]), array([[1],[2]]),
                    matrix([1,2,3]), matrix([[1],[2]])]
        for x in bad_vecs:
            pytest.raises(ValueError, Asp.__matmul__, x)

        # The current relationship between sparse matrix products and array
        # products is as follows:
        dot_result = np.dot(Asp.toarray(),[1,2,3])
        assert_array_almost_equal(Asp@array([1,2,3]), dot_result)
        assert_array_almost_equal(Asp@[[1],[2],[3]], dot_result.T)
        # Note that the result of Asp * x is dense if x has a singleton dimension.

    def test_transpose(self):
        dat_1 = self.dat1d
        dat_2 = np.array([])
        matrices = [dat_1, dat_2]

        for dtype in self.checked_dtypes:
            for j in range(len(matrices)):
                dat = np.array(matrices[j], dtype=dtype)
                datsp = self.spcreator(dat)

                a = datsp.transpose()
                b = dat.transpose()

                assert_array_equal(a.toarray(), b)
                assert_array_equal(a.transpose().toarray(), dat)
                assert_equal(a.dtype, b.dtype)

    def test_add_dense(self):
        for dtype in self.math_dtypes:
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            # adding a dense matrix to a sparse matrix
            sum1 = dat + datsp
            assert_array_equal(sum1, dat + dat)
            sum2 = datsp + dat
            assert_array_equal(sum2, dat + dat)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_sub_dense(self):
        # subtracting a dense matrix to/from a sparse matrix
        for dtype in self.math_dtypes:
            if dtype == np.dtype('bool'):
                # boolean array subtraction deprecated in 1.9.0
                continue
            dat = self.dat_dtypes[dtype]
            datsp = self.datsp_dtypes[dtype]

            # Manually add to avoid upcasting from scalar
            # multiplication.
            sum1 = (dat + dat + dat) - datsp
            assert_array_equal(sum1, dat + dat)
            sum2 = (datsp + datsp + datsp) - dat
            assert_array_equal(sum2, dat + dat)

    # test that __iter__ is compatible with NumPy matrix
    def test_iterator(self):
        B = np.arange(5)
        A = self.spcreator(B)

        if A.format not in ['coo', 'dia', 'bsr']:
            for x, y in zip(A, B):
                assert_equal(x, y)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_size_zero_matrix_arithmetic(self):
        # Test basic matrix arithmetic with shapes like (0,0), (10,0),
        # (0, 3), etc.
        mat = np.array([])
        a = mat.reshape(0)
        d = mat.reshape((1, 0))
        f = np.ones([5, 5])

        asp = self.spcreator(a)
        creator = self.spcreator if asp.format != 'uni' else csr_array
        dsp = creator(d)
        # bad shape for addition
        pytest.raises(ValueError, asp.__add__, dsp)

        # matrix product.
        assert_equal(asp.dot(asp), np.dot(a, a))

        # bad matrix products
        pytest.raises(ValueError, asp.dot, f)

        # elemente-wise multiplication
        assert_array_equal(asp.multiply(asp).toarray(), np.multiply(a, a))

        assert_array_equal(asp.multiply(a).toarray(), np.multiply(a, a))

        assert_array_equal(asp.multiply(6).toarray(), np.multiply(a, 6))

        # bad element-wise multiplication
        pytest.raises(ValueError, asp.multiply, f)

        # Addition
        assert_array_equal(asp.__add__(asp).toarray(), a.__add__(a))

    def test_resize(self):
        # resize(shape) resizes the matrix in-place
        D = np.array([1, 0, 3, 4])
        S = self.spcreator(D)
        assert_(S.resize((3,)) is None)
        assert_array_equal(S.toarray(), [1, 0, 3])
        S.resize((5,))
        assert_array_equal(S.toarray(), [1, 0, 3, 0, 0])


class _GetSet1D:
    def test_getelement(self):
        D = np.array([4,3,0])
        A = self.spcreator(D)
        N = D.shape[0]

        for j in range(-N, N):
            assert_equal(A[j], D[j])
        assert_array_equal(A[()].toarray(), D[()])
        for ij in [(0,3),(3,),3,-4]:
            pytest.raises((IndexError, TypeError), A.__getitem__, ij)


    def test_setelement(self):
        dtype = np.float64
        A = self.spcreator((12,), dtype=dtype)
        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
                "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            A[0] = dtype(0)
            A[np.array(6)] = dtype(4.0)  # scalar index
            A[1] = dtype(3)
            A[8] = dtype(9.0)
            A[-9,] = dtype(8)
            A[-2] = dtype(7)
            A[np.array(8)] = dtype(2.0)  # overwrite with scalar index
            A[1,] = dtype(5)  # overwrite with 1-tuple index 
            A[5] = 9

            assert_array_equal(A.toarray(), [0, 5, 0, 8, 0, 9, 4, 0, 2, 0, 7, 0])

        for ij in [(13,),13,-14]:
            pytest.raises(IndexError, A.__setitem__, ij, 123.0)

        for v in [(), (0, 3), [1,2,3], np.array([1,2,3])]:
            pytest.raises(ValueError, A.__setitem__, 0, v)


@contextlib.contextmanager
def check_remains_sorted(X):
    """Checks that sorted indices property is retained through an operation
    """
    if not hasattr(X, 'has_sorted_indices') or not X.has_sorted_indices:
        yield
        return
    yield
    indices = X.indices.copy()
    X.has_sorted_indices = False
    X.sort_indices()
    assert_array_equal(indices, X.indices,
                       'Expected sorted indices, found unsorted')


class _SlicingAndFancy1D:
    ####################
    #  1d Slice as index
    ####################
    def test_dtype_preservation(self):
        assert_equal(self.spcreator((10,), dtype=np.int16)[1:5].dtype, np.int16)
        assert_equal(self.spcreator((6,), dtype=np.int32)[0:0:2].dtype, np.int32)
        assert_equal(self.spcreator((6,), dtype=np.int64)[:].dtype, np.int64)

    def test_get_1d_slice(self):
        B = np.arange(50.)
        A = self.spcreator(B)
        assert_array_equal(B[:], A[:].toarray())
        assert_array_equal(B[2:5], A[2:5].toarray())

        C = np.array([4, 0, 6, 0, 0, 0, 0, 0, 1])
        D = self.spcreator(C)
        assert_array_equal(C[1:3], D[1:3].toarray())

        # Now test slicing when a row contains only zeros
        E = np.array([0, 0, 0, 0, 0])
        F = self.spcreator(E)
        assert_array_equal(E[1:3], F[1:3].toarray())
        assert_array_equal(E[-2:], F[-2:].toarray())
        assert_array_equal(E[:], F[:].toarray())
        assert_array_equal(E[slice(None)], F[slice(None)].toarray())

    def test_slicing_idx_slice(self):
        B = np.arange(50)
        A = self.spcreator(B)

        # [i]
        assert_equal(A[2], B[2])
        assert_equal(A[-1], B[-1])
        assert_equal(A[array(-2)], B[-2])

        # [1:2]
        assert_equal(A[:].toarray(), B[:])
        assert_equal(A[5:-2].toarray(), B[5:-2])
        assert_equal(A[5:12:3].toarray(), B[5:12:3])

        # int8 slice
        s = slice(np.int8(2), np.int8(4), None)
        assert_equal(A[s].toarray(), B[2:4])

        # np.s_
        s_ = np.s_
        slices = [s_[:2], s_[1:2], s_[3:], s_[3::2],
                  s_[15:20], s_[3:2],
                  s_[8:3:-1], s_[4::-2], s_[:5:-1],
                  0, 1, s_[:], s_[1:5], -1, -2, -5,
                  np.array(-1), np.int8(-3)]

        for j, a in enumerate(slices):
            x = A[a]
            y = B[a]
            if y.shape == ():
                assert_equal(x, y, repr(a))
            else:
                if x.size == 0 and y.size == 0:
                    pass
                else:
                    assert_array_equal(x.toarray(), y, repr(a))

    def test_ellipsis_1d_slicing(self):
        B = np.arange(50)
        A = self.spcreator(B)
        assert_array_equal(A[...].toarray(), B[...])
        assert_array_equal(A[...,].toarray(), B[...,])

    ##########################
    #  Assignment with Slicing
    ##########################
    def test_slice_scalar_assign(self):
        A = self.spcreator((5,))
        B = np.zeros((5,))
        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
               "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            for C in [A, B]:
                C[0:1] = 1
                C[2:0] = 4
                C[2:3] = 9
                C[3:] = 1
                C[3::-1] = 9
        assert_array_equal(A.toarray(), B)

    def test_slice_assign_2(self):
        shape = (10,)

        for idx in [slice(3), slice(None, 10, 4), slice(5, -2)]:
            A = self.spcreator(shape)
            with suppress_warnings() as sup:
                sup.filter(
                    SparseEfficiencyWarning,
                    "Changing the sparsity structure of a cs[cr]_matrix is expensive"
                )
                A[idx] = 1
            B = np.zeros(shape)
            B[idx] = 1
            msg = f"idx={idx!r}"
            assert_array_almost_equal(A.toarray(), B, err_msg=msg)

    def test_self_self_assignment(self):
        # Tests whether a row of one lil_matrix can be assigned to another.
        B = self.spcreator((5,))
        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
               "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            B[0] = 2
            B[1] = 0
            B[2] = 3
            B[3] = 10

            A = B / 10
            B[:] = A[:]
            assert_array_equal(A[:].toarray(), B[:].toarray())
            B.eliminate_zeros()

            A = B / 10
            B[:] = A[:1]
            assert_array_equal(np.zeros((5,)) + A[0], B.toarray())

            A = B / 10
            B[:-1] = A[1:]
            assert_array_equal(A[1:].toarray(), B[:-1].toarray())

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_slice_assignment(self):
        B = self.spcreator((4,))
        expected = np.array([10, 0, 14, 0])
        block = [2, 1]

        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
               "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            B[0] = 5
            B[2] = 7
            B[:] = B+B
            assert_array_equal(B.toarray(), expected)

            B[:2] = csc_array(block)
            assert_array_equal(B.toarray()[:2], block)

    def test_set_slice(self):
        A = self.spcreator((5,))
        B = np.zeros(5, float)
        s_ = np.s_
        slices = [s_[:2], s_[1:2], s_[3:], s_[3::2],
                  s_[8:3:-1], s_[4::-2], s_[:5:-1],
                  0, 1, s_[:], s_[1:5], -1, -2, -5,
                  np.array(-1), np.int8(-3)]

        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
               "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            for j, a in enumerate(slices):
                A[a] = j
                B[a] = j
                assert_array_equal(A.toarray(), B, repr(a))

            A[1:10:2] = range(1, 5, 2)
            B[1:10:2] = range(1, 5, 2)
            assert_array_equal(A.toarray(), B)

        # The next commands should raise exceptions
        toobig = list(range(100))
        pytest.raises(ValueError, A.__setitem__, 0, toobig)
        pytest.raises(ValueError, A.__setitem__, slice(None), toobig)

    def test_assign_empty(self):
        A = self.spcreator(np.ones(3))
        B = self.spcreator((2,))
        A[:2] = B
        assert_array_equal(A.toarray(), [0, 0, 1])

    ####################
    #  1d Fancy Indexing
    ####################
    def test_dtype_preservation_empty_index(self):
        A = self.spcreator((2,), dtype=np.int16)
        assert_equal(A[[False, False]].dtype, np.int16)
        assert_equal(A[[]].dtype, np.int16)

    def test_bad_index(self):
        A = self.spcreator(np.zeros(5))
        pytest.raises((IndexError, ValueError, TypeError), A.__getitem__, "foo")
        pytest.raises((IndexError, ValueError, TypeError), A.__getitem__, (2, "foo"))

    def test_fancy_indexing(self):
        B = np.arange(50)
        A = self.spcreator(B)

        # [i]
        assert_equal(A[[3]].toarray(), B[[3]])

        # [array]
        assert_equal(A[[1, 3]].toarray(), B[[1, 3]])
        assert_equal(A[[2, -5]].toarray(), B[[2, -5]])
        assert_equal(A[array(-1)], B[-1])
        assert_equal(A[array([-1, 2])].toarray(), B[[-1, 2]])
        assert_equal(A[array(5)], B[array(5)])


        # [[[1],[2]]]
        assert_equal(A[[[1], [3]]].toarray(), B[[[1], [3]]])
        assert_equal(
            A[[[-1], [-3], [-2]]].toarray(),
            B[[[-1], [-3], [-2]]]
        )

        # [[1,2]]
        assert_equal(A[[1, 3]].toarray(), B[[1, 3]])
        assert_equal(A[[-1, -3]].toarray(), B[[-1, -3]])
        assert_equal(A[array([-1, -3])].toarray(), B[[-1, -3]])

        # [[1,2]][[1,2]]
        assert_equal(A[[1, 5, 2, 8]][[1, 3]].toarray(),
                     B[[1, 5, 2, 8]][[1, 3]])
        assert_equal(A[[-1, -5, 2, 8]][[1, -4]].toarray(),
                     B[[-1, -5, 2, 8]][[1, -4]])

    def test_fancy_indexing_boolean(self):
        np.random.seed(1234)  # make runs repeatable

        B = np.arange(50)
        A = self.spcreator(B)

        I = np.array(np.random.randint(0, 2, size=50), dtype=bool)

        assert_equal(toarray(A[I]), B[I])
        assert_equal(toarray(A[B > 9]), B[B > 9])

        Z1 = np.zeros(51, dtype=bool)
        Z2 = np.zeros(51, dtype=bool)
        Z2[-1] = True
        Z3 = np.zeros(51, dtype=bool)
        Z3[0] = True

        pytest.raises(IndexError, A.__getitem__, Z1)
        pytest.raises(IndexError, A.__getitem__, Z2)
        pytest.raises(IndexError, A.__getitem__, Z3)

    @pytest.mark.skip(reason="tocsr() not valid for 1d sparse array")
    def test_fancy_indexing_sparse_boolean(self):
        np.random.seed(1234)  # make runs repeatable

        B = np.arange(20)
        A = self.spcreator(B)

        X = np.array(np.random.randint(0, 2, size=20), dtype=bool)
        Xsp = csr_array(X)

        assert_equal(toarray(A[Xsp]), B[X])
        assert_equal(toarray(A[A > 9]), B[B > 9])

        Y = np.array(np.random.randint(0, 2, size=60), dtype=bool)

        Ysp = csr_array(Y)

        pytest.raises(IndexError, A.__getitem__, Ysp)
        pytest.raises((IndexError, ValueError), A.__getitem__, (Xsp, 1))

    def test_fancy_indexing_seq_assign(self):
        mat = self.spcreator(np.array([1, 0]))
        pytest.raises(ValueError, mat.__setitem__, 0, np.array([1,2]))

    def test_fancy_indexing_empty(self):
        B = np.arange(50)
        B[3:9] = 0
        B[30] = 0
        A = self.spcreator(B)

        K = np.array([False] * 50)
        assert_equal(toarray(A[K]), B[K])
        K = np.array([], dtype=int)
        assert_equal(toarray(A[K]), B[K])
        J = np.array([0, 1, 2, 3, 4], dtype=int)
        assert_equal(toarray(A[J]), B[J])

    ############################
    #  1d Fancy Index Assignment
    ############################
    def test_bad_index_assign(self):
        A = self.spcreator(np.zeros(5))
        pytest.raises(
            (IndexError, ValueError, TypeError), A.__setitem__, "foo", 2
        )

    def test_fancy_indexing_set(self):
        M = (5,)

        # [1:2]
        for j in [[2, 3, 4], slice(None, 10, 4), np.arange(3),
                     slice(5, -2), slice(2, 5)]:
            A = self.spcreator(M)
            B = np.zeros(M)
            with suppress_warnings() as sup:
                sup.filter(
                    SparseEfficiencyWarning,
                   "Changing the sparsity structure of a cs[cr]_matrix is expensive"
                )
                B[j] = 1
                with check_remains_sorted(A):
                    A[j] = 1
            assert_array_almost_equal(A.toarray(), B)

    def test_sequence_assignment(self):
        A = self.spcreator((4,))
        B = self.spcreator((3,))

        i0 = [0,1,2]
        i1 = (0,1,2)
        i2 = np.array(i0)

        with suppress_warnings() as sup:
            sup.filter(
                SparseEfficiencyWarning,
               "Changing the sparsity structure of a cs[cr]_matrix is expensive"
            )
            with check_remains_sorted(A):
                A[i0] = B[i0]
                pytest.raises(IndexError, B.__getitem__, i1)
                A[i2] = B[i2]
            assert_array_equal(A[:3].toarray(), B.toarray())
            assert A.shape == (4,)

            # slice
            A = self.spcreator((4,))
            with check_remains_sorted(A):
                A[1:3] = [10,20]
            assert_array_equal(A.toarray(), [0, 10, 20, 0])

            # array
            A = self.spcreator((4,))
            B = np.zeros(4)
            with check_remains_sorted(A):
                for C in [A, B]:
                    C[[0,1,2]] = [4,5,6]
            assert_array_equal(A.toarray(), B)

    def test_fancy_assign_empty(self):
        B = np.arange(50)
        B[2] = 0
        B[[3, 6]] = 0
        A = self.spcreator(B)

        K = np.array([False] * 50)
        A[K] = 42
        assert_equal(A.toarray(), B)

        K = np.array([], dtype=int)
        A[K] = 42
        assert_equal(A.toarray(), B)


class _MinMaxMixin1D:
    def test_minmax(self):
        D = np.arange(5)
        X = self.spcreator(D)

        assert_equal(X.min(), 0)
        assert_equal(X.max(), 4)
        assert_equal((-X).min(), -4)
        assert_equal((-X).max(), 0)


    def test_minmax_axis(self):
        D = np.arange(50)
        X = self.spcreator(D)

        for axis in [0, -1]:
            assert_array_equal(
                toarray(X.max(axis=axis)), D.max(axis=axis, keepdims=True)
            )
            assert_array_equal(
                toarray(X.min(axis=axis)), D.min(axis=axis, keepdims=True)
            )
        for axis in [-2, 1]:
            pytest.raises(ValueError, X.min, axis=axis)
            pytest.raises(ValueError, X.max, axis=axis)


    def test_numpy_minmax(self):
        dat = np.array([0, 1, 2])
        datsp = self.spcreator(dat)
        assert_array_equal(np.min(datsp), np.min(dat))
        assert_array_equal(np.max(datsp), np.max(dat))


    def test_argmax(self):
        D1 = np.array([-1, 5, 2, 3])
        D2 = np.array([0, 0, -1, -2])
        D3 = np.array([-1, -2, -3, -4])
        D4 = np.array([1, 2, 3, 4])
        D5 = np.array([1, 2, 0, 0])

        for D in [D1, D2, D3, D4, D5]:
            mat = self.spcreator(D)

            assert_equal(mat.argmax(), np.argmax(D))
            assert_equal(mat.argmin(), np.argmin(D))

            assert_equal(mat.argmax(axis=0), np.argmax(D, axis=0))
            assert_equal(mat.argmin(axis=0), np.argmin(D, axis=0))

        D6 = np.empty((0,))

        for axis in [None, 0]:
            mat = self.spcreator(D6)
            pytest.raises(ValueError, mat.argmax, axis=axis)
            pytest.raises(ValueError, mat.argmin, axis=axis)


class TestCOO1D(_Common1D, _MinMaxMixin1D):
    spcreator = coo_array
    datsp = spcreator(_Common1D.dat1d)
    datsp_dtypes = {}
    for dtype in _Common1D.checked_dtypes:
        datsp_dtypes[dtype] = spcreator(_Common1D.dat1d.astype(dtype))


#class TestUNI1D(_Common1D, _MinMaxMixin1D):
#    spcreator = uni_array
#    datsp = spcreator(_Common1D.dat1d)
#    datsp_dtypes = {}
#    for dtype in _Common1D.checked_dtypes:
#        datsp_dtypes[dtype] = spcreator(_Common1D.dat1d.astype(dtype))


#class TestCSC1D(_Common1D, _MinMaxMixin1D, _GetSet1D, _SlicingAndFancy1D):
#    spcreator = csc_array
#    datsp = spcreator(_Common1D.dat1d)
#    datsp_dtypes = {}
#    for dtype in _Common1D.checked_dtypes:
#        datsp_dtypes[dtype] = spcreator(_Common1D.dat1d.astype(dtype))


#class TestCSR1D(_Common1D, _MinMaxMixin1D, _SlicingAndFancy1D):
#    spcreator = csr_array
#    datsp = spcreator(_Common1D.dat1d)
#    datsp_dtypes = {}
#    for dtype in _Common1D.checked_dtypes:
#        datsp_dtypes[dtype] = spcreator(_Common1D.dat1d.astype(dtype))


class TestDOK1D(_Common1D, _SlicingAndFancy1D):
    spcreator = dok_array
    datsp = spcreator(_Common1D.dat1d)
    datsp_dtypes = {}
    for dtype in _Common1D.checked_dtypes:
        datsp_dtypes[dtype] = spcreator(_Common1D.dat1d.astype(dtype))

