"""test sparse matrix construction functions"""

import numpy as np
from numpy import array
from numpy.testing import (assert_equal, assert_,
        assert_array_equal, assert_array_almost_equal_nulp)
import pytest
from pytest import raises as assert_raises
from scipy._lib._testutils import check_free_memory

from scipy.sparse import (csr_array, coo_array,
                          csc_array, bsr_array,
                          dia_array, dok_array,
                          lil_array, sparray,
                          _construct as construct)
from scipy.sparse._construct import random_array as sprand

sparse_formats = ['csr','csc','coo','bsr','dia','lil','dok']

#TODO check whether format=XXX is respected


def _sprandn(shape, density=0.01, format="coo", dtype=None, rng=None):
    # Helper function for testing.
    rng = np.random.default_rng(rng)
    data_sampler = rng.standard_normal
    return construct.random_array(shape, density=density, format=format,
                                  dtype=dtype, rng=rng, data_sampler=data_sampler)


class TestConstructUtils:

    @pytest.mark.parametrize("cls", [
        csc_array, csr_array, coo_array, bsr_array,
        dia_array, dok_array, lil_array
    ])
    def test_singleton_array_constructor(self, cls):
        with pytest.raises(
            ValueError,
            match=(
                'scipy sparse array classes do not support '
                'instantiation from a scalar'
            )
        ):
            cls(0)

    @pytest.mark.filterwarnings("ignore:.* output has been cast to:FutureWarning")
    def test_spdiags(self):
        diags1 = array([[1, 2, 3, 4, 5]])
        diags2 = array([[1, 2, 3, 4, 5],
                         [6, 7, 8, 9,10]])
        diags3 = array([[1, 2, 3, 4, 5],
                         [6, 7, 8, 9,10],
                         [11,12,13,14,15]])

        # case => (diags, offsets, m, n, result)
        cases = []
        cases.append((diags1, 0, 1, 1, [[1]]))
        cases.append((diags1, [0], 1, 1, [[1]]))
        cases.append((diags1, [0], 2, 1, [[1],[0]]))
        cases.append((diags1, [0], 1, 2, [[1,0]]))
        cases.append((diags1, [1], 1, 2, [[0,2]]))
        cases.append((diags1,[-1], 1, 2, [[0,0]]))
        cases.append((diags1, [0], 2, 2, [[1,0],[0,2]]))
        cases.append((diags1,[-1], 2, 2, [[0,0],[1,0]]))
        cases.append((diags1, [3], 2, 2, [[0,0],[0,0]]))
        cases.append((diags1, [0], 3, 4, [[1,0,0,0],[0,2,0,0],[0,0,3,0]]))
        cases.append((diags1, [1], 3, 4, [[0,2,0,0],[0,0,3,0],[0,0,0,4]]))
        cases.append((diags1, [2], 3, 5, [[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]]))

        cases.append((diags2, [0,2], 3, 3, [[1,0,8],[0,2,0],[0,0,3]]))
        cases.append((diags2, [-1,0], 3, 4, [[6,0,0,0],[1,7,0,0],[0,2,8,0]]))
        cases.append((diags2, [2,-3], 6, 6, [[0,0,3,0,0,0],
                                              [0,0,0,4,0,0],
                                              [0,0,0,0,5,0],
                                              [6,0,0,0,0,0],
                                              [0,7,0,0,0,0],
                                              [0,0,8,0,0,0]]))

        cases.append((diags3, [-1,0,1], 6, 6, [[6,12, 0, 0, 0, 0],
                                                [1, 7,13, 0, 0, 0],
                                                [0, 2, 8,14, 0, 0],
                                                [0, 0, 3, 9,15, 0],
                                                [0, 0, 0, 4,10, 0],
                                                [0, 0, 0, 0, 5, 0]]))
        cases.append((diags3, [-4,2,-1], 6, 5, [[0, 0, 8, 0, 0],
                                                 [11, 0, 0, 9, 0],
                                                 [0,12, 0, 0,10],
                                                 [0, 0,13, 0, 0],
                                                 [1, 0, 0,14, 0],
                                                 [0, 2, 0, 0,15]]))
        cases.append((diags3, [-1, 1, 2], len(diags3[0]), len(diags3[0]),
                      [[0, 7, 13, 0, 0],
                       [1, 0, 8, 14, 0],
                       [0, 2, 0, 9, 15],
                       [0, 0, 3, 0, 10],
                       [0, 0, 0, 4, 0]]))

        for d, o, m, n, result in cases:
            assert_equal(dia_array((d, o), shape=(m, n)).toarray(), result)

    def test_diags(self):
        a = array([1.0, 2.0, 3.0, 4.0, 5.0])
        b = array([6.0, 7.0, 8.0, 9.0, 10.0])
        c = array([11.0, 12.0, 13.0, 14.0, 15.0])

        cases = []
        cases.append((a[:1], 0, (1, 1), [[1]]))
        cases.append(([a[:1]], [0], (1, 1), [[1]]))
        cases.append(([a[:1]], [0], (2, 1), [[1],[0]]))
        cases.append(([a[:1]], [0], (1, 2), [[1,0]]))
        cases.append(([a[:1]], [1], (1, 2), [[0,1]]))
        cases.append(([a[:2]], [0], (2, 2), [[1,0],[0,2]]))
        cases.append(([a[:1]],[-1], (2, 2), [[0,0],[1,0]]))
        cases.append(([a[:3]], [0], (3, 4), [[1,0,0,0],[0,2,0,0],[0,0,3,0]]))
        cases.append(([a[:3]], [1], (3, 4), [[0,1,0,0],[0,0,2,0],[0,0,0,3]]))
        cases.append(([a[:1]], [-2], (3, 5), [[0,0,0,0,0],[0,0,0,0,0],[1,0,0,0,0]]))
        cases.append(([a[:2]], [-1], (3, 5), [[0,0,0,0,0],[1,0,0,0,0],[0,2,0,0,0]]))
        cases.append(([a[:3]], [0], (3, 5), [[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0]]))
        cases.append(([a[:3]], [1], (3, 5), [[0,1,0,0,0],[0,0,2,0,0],[0,0,0,3,0]]))
        cases.append(([a[:3]], [2], (3, 5), [[0,0,1,0,0],[0,0,0,2,0],[0,0,0,0,3]]))
        cases.append(([a[:2]], [3], (3, 5), [[0,0,0,1,0],[0,0,0,0,2],[0,0,0,0,0]]))
        cases.append(([a[:1]], [4], (3, 5), [[0,0,0,0,1],[0,0,0,0,0],[0,0,0,0,0]]))
        cases.append(([a[:1]], [-4], (5, 3), [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[1,0,0]]))
        cases.append(([a[:2]], [-3], (5, 3), [[0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,2,0]]))
        cases.append(([a[:3]], [-2], (5, 3), [[0,0,0],[0,0,0],[1,0,0],[0,2,0],[0,0,3]]))
        cases.append(([a[:3]], [-1], (5, 3), [[0,0,0],[1,0,0],[0,2,0],[0,0,3],[0,0,0]]))
        cases.append(([a[:3]], [0], (5, 3), [[1,0,0],[0,2,0],[0,0,3],[0,0,0],[0,0,0]]))
        cases.append(([a[:2]], [1], (5, 3), [[0,1,0],[0,0,2],[0,0,0],[0,0,0],[0,0,0]]))
        cases.append(([a[:1]], [2], (5, 3), [[0,0,1],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]))

        cases.append(([a[:3],b[:1]], [0,2], (3, 3), [[1,0,6],[0,2,0],[0,0,3]]))
        cases.append(([a[:2],b[:3]], [-1,0], (3, 4), [[6,0,0,0],[1,7,0,0],[0,2,8,0]]))
        cases.append(([a[:4],b[:3]], [2,-3], (6, 6), [[0,0,1,0,0,0],
                                                     [0,0,0,2,0,0],
                                                     [0,0,0,0,3,0],
                                                     [6,0,0,0,0,4],
                                                     [0,7,0,0,0,0],
                                                     [0,0,8,0,0,0]]))

        cases.append(([a[:4],b,c[:4]], [-1,0,1], (5, 5), [[6,11, 0, 0, 0],
                                                            [1, 7,12, 0, 0],
                                                            [0, 2, 8,13, 0],
                                                            [0, 0, 3, 9,14],
                                                            [0, 0, 0, 4,10]]))
        cases.append(([a[:2],b[:3],c], [-4,2,-1], (6, 5), [[0, 0, 6, 0, 0],
                                                          [11, 0, 0, 7, 0],
                                                          [0,12, 0, 0, 8],
                                                          [0, 0,13, 0, 0],
                                                          [1, 0, 0,14, 0],
                                                          [0, 2, 0, 0,15]]))

        # too long arrays are OK
        cases.append(([a], [0], (1, 1), [[1]]))
        cases.append(([a[:3],b], [0,2], (3, 3), [[1, 0, 6], [0, 2, 0], [0, 0, 3]]))
        cases.append((
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            [0,-1],
            (3, 3),
            [[1, 0, 0], [4, 2, 0], [0, 5, 3]]
        ))

        # scalar case: broadcasting
        cases.append(([1.0,-2.0,1.0], [1,0,-1], (3, 3), [[-2, 1, 0],
                                                         [1, -2, 1],
                                                         [0, 1, -2]]))

        for d, o, shape, result in cases:
            err_msg = f"{d!r} {o!r} {shape!r} {result!r}"
            assert_equal(construct.diags_array(d, offsets=o, shape=shape).toarray(),
                         result, err_msg=err_msg)

            if (shape[0] == shape[1]
                and hasattr(d[0], '__len__')
                and len(d[0]) <= max(shape)):
                # should be able to find the shape automatically
                assert_equal(construct.diags_array(d, offsets=o).toarray(),
                             result, err_msg=err_msg)

    def test_diags_default(self):
        a = array([1.0, 2.0, 3.0, 4.0, 5.0])
        assert_equal(construct.diags_array(a).toarray(), np.diag(a))

    def test_diags_default_bad(self):
        a = array([[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]])
        assert_raises(ValueError, construct.diags_array, a)

    def test_diags_bad(self):
        a = array([1.0, 2.0, 3.0, 4.0, 5.0])
        b = array([6.0, 7.0, 8.0, 9.0, 10.0])
        c = array([11.0, 12.0, 13.0, 14.0, 15.0])

        cases = []
        cases.append(([a[:0]], 0, (1, 1)))
        cases.append(([a[:4],b,c[:3]], [-1,0,1], (5, 5)))
        cases.append(([a[:2],c,b[:3]], [-4,2,-1], (6, 5)))
        cases.append(([a[:2],c,b[:3]], [-4,2,-1], None))
        cases.append(([], [-4,2,-1], None))
        cases.append(([1.0], [-5], (4, 4)))
        cases.append(([a], 0, None))

        for d, o, shape in cases:
            assert_raises(ValueError, construct.diags_array, d, offsets=o, shape=shape)

        assert_raises(TypeError, construct.diags_array, [[None]], offsets=[0])

    def test_diags_vs_diag(self):
        # Check that
        #
        #    diags([a, b, ...], [i, j, ...]) == diag(a, i) + diag(b, j) + ...
        #

        rng = np.random.RandomState(1234)

        for n_diags in [1, 2, 3, 4, 5, 10]:
            n = 1 + n_diags//2 + rng.randint(0, 10)

            offsets = np.arange(-n+1, n-1)
            rng.shuffle(offsets)
            offsets = offsets[:n_diags]

            diagonals = [rng.rand(n - abs(q)) for q in offsets]

            mat = construct.diags_array(diagonals, offsets=offsets)
            dense_mat = sum([np.diag(x, j) for x, j in zip(diagonals, offsets)])

            assert_array_almost_equal_nulp(mat.toarray(), dense_mat)

            if len(offsets) == 1:
                mat = construct.diags_array(diagonals[0], offsets=offsets[0])
                dense_mat = np.diag(diagonals[0], offsets[0])
                assert_array_almost_equal_nulp(mat.toarray(), dense_mat)

    def test_diags_dtype(self):
        x = construct.diags_array([2.2], offsets=[0], shape=(2, 2), dtype=int)
        assert_equal(x.dtype, int)
        assert_equal(x.toarray(), [[2, 0], [0, 2]])

    def test_diags_one_diagonal(self):
        d = [0.0, 1.0, 2.0, 3.0, 4.0]
        for k in range(-5, 6):
            assert_equal(construct.diags_array(d, offsets=k).toarray(),
                         construct.diags_array([d], offsets=[k]).toarray())

    def test_diags_empty(self):
        x = construct.diags_array([])
        assert_equal(x.shape, (0, 0))

    def test_identity(self):
        self.check_identity(construct.eye_array)

#    @pytest.mark.parametrize("identity", [construct.identity, construct.eye_array])
    def check_identity(self, identity):
        assert_equal(identity(1).toarray(), [[1]])
        assert_equal(identity(2).toarray(), [[1,0],[0,1]])

        I = identity(3, dtype='int8', format='dia')
        assert_equal(I.dtype, np.dtype('int8'))
        assert_equal(I.format, 'dia')

        for fmt in sparse_formats:
            I = identity(3, format=fmt)
            assert_equal(I.format, fmt)
            assert_equal(I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

#    @pytest.mark.parametrize("eye", [construct.eye, construct.eye_array])
    def test_eye(self):
        self.check_eye(construct.eye_array)

    def check_eye(self, eye):
        assert_equal(eye(1,1).toarray(), [[1]])
        assert_equal(eye(2,3).toarray(), [[1,0,0],[0,1,0]])
        assert_equal(eye(3,2).toarray(), [[1,0],[0,1],[0,0]])
        assert_equal(eye(3,3).toarray(), [[1,0,0],[0,1,0],[0,0,1]])

        assert_equal(eye(3,3,dtype='int16').dtype, np.dtype('int16'))

        for m in [3, 5]:
            for n in [3, 5]:
                for k in range(-5,6):
                    # scipy.sparse.eye deviates from np.eye here. np.eye will
                    # create arrays of all 0's when the diagonal offset is
                    # greater than the size of the array. For sparse arrays
                    # this makes less sense, especially as it results in dia
                    # arrays with negative diagonals. Therefore sp.sparse.eye
                    # validates that diagonal offsets fall within the shape of
                    # the array. See gh-18555.
                    if (k > 0 and k > n) or (k < 0 and abs(k) > m):
                        with pytest.raises(
                            ValueError, match="Offset.*out of bounds"
                        ):
                            eye(m, n, k=k)

                    else:
                        assert_equal(
                            eye(m, n, k=k).toarray(),
                            np.eye(m, n, k=k)
                        )
                        if m == n:
                            assert_equal(
                                eye(m, k=k).toarray(),
                                np.eye(m, n, k=k)
                            )

#    @pytest.mark.parametrize("eye", [construct.eye, construct.eye_array])
    def test_eye_one(self):
        self.check_eye_one(construct.eye_array)

    def check_eye_one(self, eye):
        assert_equal(eye(1).toarray(), [[1]])
        assert_equal(eye(2).toarray(), [[1,0],[0,1]])

        I = eye(3, dtype='int8', format='dia')
        assert_equal(I.dtype, np.dtype('int8'))
        assert_equal(I.format, 'dia')

        for fmt in sparse_formats:
            I = eye(3, format=fmt)
            assert_equal(I.format, fmt)
            assert_equal(I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

    @pytest.mark.parametrize("arr,kw_format,out_format", [
        ([[0, 0], [0, 1]], None, 'coo'),  # 2D sparse
        ([[1, 0], [1, 1]], None, 'bsr'),  # 2D dense
        ([[[1, 0], [1, 1]]], None, 'coo'),  # 3D dense
    ])
    def test_kron_output_format(self, arr, kw_format, out_format):
        sparr = coo_array(arr)
        assert construct.kron(sparr, sparr, format=kw_format).format == out_format
        assert construct.kron(sparr, arr, format=kw_format).format == out_format
        assert construct.kron(arr, sparr, format=kw_format).format == out_format

    def test_kron(self):
        cases = []

        cases.append(array([[0]]))
        cases.append(array([[-1]]))
        cases.append(array([[4]]))
        cases.append(array([[10]]))
        cases.append(array([[0], [0]]))
        cases.append(array([[0, 0]]))
        cases.append(array([[1, 2], [3, 4]]))
        cases.append(array([[0, 2], [5, 0]]))
        cases.append(array([[0, 2, -6], [8, 0, 14]]))
        cases.append(array([[5, 4], [0, 0], [6, 0]]))
        cases.append(array([[5, 4, 4], [1, 0, 0], [6, 0, 8]]))
        cases.append(array([[0, 1, 0, 2, 0, 5, 8]]))
        cases.append(array([[0.5, 0.125, 0, 3.25], [0, 2.5, 0, 0]]))

        # test all cases with some formats
        for a in cases:
            ca = csr_array(a)
            for b in cases:
                cb = csr_array(b)
                expected = np.kron(a, b)
                for fmt in sparse_formats[1:4]:
                    result = construct.kron(ca, cb, format=fmt)
                    assert_equal(result.format, fmt)
                    assert_array_equal(result.toarray(), expected)
                    assert isinstance(result, sparray)

        # nD cases
        cases.append(array([0, 1, 2]))
        cases.append(array([[[0, 1, 2], [0, 1, 0]]]))
        cases.append(array([[[0, 1]], [[2, 2]], [[1, 0]], [[2, 0]]]))

        for a in cases:
            ca = coo_array(a)
            for b in cases:
                cb = coo_array(b)
                expected = np.kron(a, b)
                result = construct.kron(ca, cb, format="coo")
                assert_array_equal(result.toarray(), expected)

        # test one case with all formats
        a = array([[0.5, 0.125, 0, 3.25], [0, 2.5, 0, 0]])
        b = array([[5, 4, 4], [1, 0, 0], [6, 0, 8]])
        ca = csr_array(a)
        cb = csr_array(b)

        expected = np.kron(a, b)
        for fmt in sparse_formats:
            result = construct.kron(ca, cb, format=fmt)
            assert_equal(result.format, fmt)
            assert_array_equal(result.toarray(), expected)
            assert isinstance(result, sparray)

    def test_kron_ndim_exceptions(self):
        # check that 3D input is ok
        construct.kron([[0], [1]], [[[0, 1]]])
        construct.kron([[[1, 1]]], [[1], [1]])

        # no exception for 3D if any sparse arrays input
        construct.kron(coo_array([[[0, 1]]]), [[[0], [1]]])
        construct.kron([[[0, 1]]], coo_array([[[0], [1]]]))

        # no exception for 1D if either sparray or spmatrix
        construct.kron([[0], [1]], [0, 1, 0])  # spmatrix b/c lists; 1d-list -> 2d
        construct.kron([1, 1], [[1], [1]])
        construct.kron([[0], [1]], coo_array([0, 1, 0]))  # sparray 1d-list -> 1d
        construct.kron(coo_array([1, 1]), [[1], [1]])

    def test_kron_large(self):
        n = 2**16
        a = construct.diags_array([1], shape=(1, n), offsets=n-1, dtype=None)
        b = construct.diags_array([1], shape=(n, 1), offsets=1-n, dtype=None)

        construct.kron(a, a)
        construct.kron(b, b)

    def test_kronsum(self):
        cases = []

        cases.append(array([[0]]))
        cases.append(array([[-1]]))
        cases.append(array([[4]]))
        cases.append(array([[10]]))
        cases.append(array([[1,2],[3,4]]))
        cases.append(array([[0,2],[5,0]]))
        cases.append(array([[0,2,-6],[8,0,14],[0,3,0]]))
        cases.append(array([[1,0,0],[0,5,-1],[4,-2,8]]))

        # test all cases with default format
        for a in cases:
            for b in cases:
                result = construct.kronsum(csr_array(a), csr_array(b)).toarray()
                expected = (np.kron(np.eye(b.shape[0]), a)
                            + np.kron(b, np.eye(a.shape[0])))
                assert_array_equal(result, expected)

    def test_kronsum_ndim_exceptions(self):
        with pytest.raises(ValueError, match='requires 2D input'):
            construct.kronsum([[0], [1]], csr_array([0, 1]))
        with pytest.raises(ValueError, match='requires 2D input'):
            construct.kronsum(csr_array([0, 1]), [[0], [1]])
        with pytest.raises(ValueError, match='requires 2D input'):
            construct.kronsum([[0, 1], [1, 0]], [2])

    def test_vstack(self):
        coo_cls = coo_array
        A = coo_cls([[1,2],[3,4]])
        B = coo_cls([[5,6]])

        expected = array([[1, 2],
                          [3, 4],
                          [5, 6]])
        assert_equal(construct.vstack([A, B]).toarray(), expected)
        assert_equal(construct.vstack([A, B], dtype=np.float32).dtype,
                     np.float32)

        assert_equal(construct.vstack([A.todok(), B.todok()]).toarray(), expected)

        assert_equal(construct.vstack([A.tocsr(), B.tocsr()]).toarray(),
                     expected)
        result = construct.vstack([A.tocsr(), B.tocsr()],
                                  format="csr", dtype=np.float32)
        assert_equal(result.dtype, np.float32)
        assert_equal(result.indices.dtype, np.int32)
        assert_equal(result.indptr.dtype, np.int32)

        assert_equal(construct.vstack([A.tocsc(), B.tocsc()]).toarray(),
                     expected)
        result = construct.vstack([A.tocsc(), B.tocsc()],
                                  format="csc", dtype=np.float32)
        assert_equal(result.dtype, np.float32)
        assert_equal(result.indices.dtype, np.int32)
        assert_equal(result.indptr.dtype, np.int32)

    def test_vstack_maintain64bit_idx_dtype(self):
        # see gh-20389 v/hstack returns int32 idx_dtype with input int64 idx_dtype
        X = csr_array([[1, 0, 0], [0, 1, 0], [0, 1, 0]])
        X.indptr = X.indptr.astype(np.int64)
        X.indices = X.indices.astype(np.int64)
        assert construct.vstack([X, X]).indptr.dtype == np.int64
        assert construct.hstack([X, X]).indptr.dtype == np.int64

        X = csc_array([[1, 0, 0], [0, 1, 0], [0, 1, 0]])
        X.indptr = X.indptr.astype(np.int64)
        X.indices = X.indices.astype(np.int64)
        assert construct.vstack([X, X]).indptr.dtype == np.int64
        assert construct.hstack([X, X]).indptr.dtype == np.int64

        X = coo_array([[1, 0, 0], [0, 1, 0], [0, 1, 0]])
        X.coords = tuple(co.astype(np.int64) for co in X.coords)
        assert construct.vstack([X, X]).coords[0].dtype == np.int64
        assert construct.hstack([X, X]).coords[0].dtype == np.int64

    def test_vstack_1d_with_2d(self):
        # fixes gh-21064
        arr = csr_array([[1, 0, 0], [0, 1, 0]])
        arr1d = csr_array([1, 0, 0])
        arr1dcoo = coo_array([1, 0, 0])
        assert construct.vstack([arr, np.array([0, 0, 0])]).shape == (3, 3)
        assert construct.hstack([arr1d, np.array([[0]])]).shape == (1, 4)
        assert construct.hstack([arr1d, arr1d]).shape == (1, 6)
        assert construct.vstack([arr1d, arr1d]).shape == (2, 3)

        # check csr specialty stacking code like _stack_along_minor_axis
        assert construct.hstack([arr, arr]).shape == (2, 6)
        assert construct.hstack([arr1d, arr1d]).shape == (1, 6)

        assert construct.hstack([arr1d, arr1dcoo]).shape == (1, 6)
        assert construct.vstack([arr, arr1dcoo]).shape == (3, 3)
        assert construct.vstack([arr1d, arr1dcoo]).shape == (2, 3)

        with pytest.raises(ValueError, match="incompatible row dimensions"):
            construct.hstack([arr, np.array([0, 0])])
        with pytest.raises(ValueError, match="incompatible column dimensions"):
            construct.vstack([arr, np.array([0, 0])])

    def test_hstack(self):
        coo_cls = coo_array
        A = coo_cls([[1,2],[3,4]])
        B = coo_cls([[5],[6]])

        expected = array([[1, 2, 5],
                          [3, 4, 6]])
        assert_equal(construct.hstack([A, B]).toarray(), expected)
        assert_equal(construct.hstack([A, B], dtype=np.float32).dtype,
                     np.float32)

        assert_equal(construct.hstack([A.todok(), B.todok()]).toarray(), expected)

        assert_equal(construct.hstack([A.tocsc(), B.tocsc()]).toarray(),
                     expected)
        assert_equal(construct.hstack([A.tocsc(), B.tocsc()],
                                      dtype=np.float32).dtype,
                     np.float32)
        assert_equal(construct.hstack([A.tocsr(), B.tocsr()]).toarray(),
                     expected)
        assert_equal(construct.hstack([A.tocsr(), B.tocsr()],
                                      dtype=np.float32).dtype,
                     np.float32)

    def test_block_creation(self):
        block_array = construct.block_array

        A = coo_array([[1, 2], [3, 4]])
        B = coo_array([[5],[6]])
        C = coo_array([[7]])
        D = coo_array((0, 0))

        expected = array([[1, 2, 5],
                          [3, 4, 6],
                          [0, 0, 7]])
        assert_equal(block_array([[A, B], [None, C]]).toarray(), expected)
        E = csr_array((1, 2), dtype=np.int32)
        assert_equal(block_array([[A.tocsr(), B.tocsr()],
                                  [E, C.tocsr()]]).toarray(),
                     expected)
        assert_equal(block_array([[A.tocsc(), B.tocsc()],
                                  [E.tocsc(), C.tocsc()]]).toarray(),
                     expected)

        expected = array([[1, 2, 0],
                          [3, 4, 0],
                          [0, 0, 7]])
        assert_equal(block_array([[A, None], [None, C]]).toarray(), expected)
        assert_equal(block_array([[A.tocsr(), E.T.tocsr()],
                                  [E, C.tocsr()]]).toarray(),
                     expected)
        assert_equal(block_array([[A.tocsc(), E.T.tocsc()],
                                  [E.tocsc(), C.tocsc()]]).toarray(),
                     expected)

        Z = csr_array((1, 1), dtype=np.int32)
        expected = array([[0, 5],
                          [0, 6],
                          [7, 0]])
        assert_equal(block_array([[None, B], [C, None]]).toarray(), expected)
        assert_equal(block_array([[E.T.tocsr(), B.tocsr()],
                                  [C.tocsr(), Z]]).toarray(),
                     expected)
        assert_equal(block_array([[E.T.tocsc(), B.tocsc()],
                                  [C.tocsc(), Z.tocsc()]]).toarray(),
                     expected)

        expected = np.empty((0, 0))
        assert_equal(block_array([[None, None]]).toarray(), expected)
        assert_equal(block_array([[None, D], [D, None]]).toarray(),
                     expected)

        # test bug reported in gh-5976
        expected = array([[7]])
        assert_equal(block_array([[None, D], [C, None]]).toarray(),
                     expected)

        # test failure cases
        with assert_raises(ValueError) as excinfo:
            block_array([[A], [B]])
        excinfo.match(r'Got blocks\[1,0\]\.shape\[1\] == 1, expected 2')

        with assert_raises(ValueError) as excinfo:
            block_array([[A.tocsr()], [B.tocsr()]])
        excinfo.match(r'incompatible dimensions for axis 1')

        with assert_raises(ValueError) as excinfo:
            block_array([[A.tocsc()], [B.tocsc()]])
        excinfo.match(r'Mismatching dimensions along axis 1: ({1, 2}|{2, 1})')

        with assert_raises(ValueError) as excinfo:
            block_array([[A, C]])
        excinfo.match(r'Got blocks\[0,1\]\.shape\[0\] == 1, expected 2')

        with assert_raises(ValueError) as excinfo:
            block_array([[A.tocsr(), C.tocsr()]])
        excinfo.match(r'Mismatching dimensions along axis 0: ({1, 2}|{2, 1})')

        with assert_raises(ValueError) as excinfo:
            block_array([[A.tocsc(), C.tocsc()]])
        excinfo.match(r'incompatible dimensions for axis 0')

    def test_block_return_type(self):
        block = construct.block_array

        # csr format ensures we hit _compressed_sparse_stack
        # shape of F,G ensure we hit _stack_along_minor_axis
        # list version ensure we hit the path with neither helper function
        Fl, Gl = [[1, 2],[3, 4]], [[7], [5]]
        Fm, Gm = csr_array(Fl), csr_array(Gl)
        assert isinstance(block([[None, Fl], [Gl, None]], format="csr"), sparray)
        assert isinstance(block([[None, Fm], [Gm, None]], format="csr"), sparray)
        assert isinstance(block([[Fm, Gm]], format="csr"), sparray)


    @pytest.mark.xslow
    @pytest.mark.xfail_on_32bit("Can't create large array for test")
    def test_concatenate_int32_overflow(self):
        """ test for indptr overflow when concatenating matrices """
        check_free_memory(30000)

        n = 33000
        A = csr_array(np.ones((n, n), dtype=bool))
        B = A.copy()
        C = construct._compressed_sparse_stack((A, B), axis=0,
                                               return_spmatrix=False)

        assert_(np.all(np.equal(np.diff(C.indptr), n)))
        assert_equal(C.indices.dtype, np.int64)
        assert_equal(C.indptr.dtype, np.int64)

    def test_block_diag_basic(self):
        """ basic test for block_diag """
        A = coo_array([[1,2],[3,4]])
        B = coo_array([[5],[6]])
        C = coo_array([[7]])

        expected = array([[1, 2, 0, 0],
                          [3, 4, 0, 0],
                          [0, 0, 5, 0],
                          [0, 0, 6, 0],
                          [0, 0, 0, 7]])

        ABC = construct.block_diag((A, B, C))
        assert_equal(ABC.toarray(), expected)
        assert ABC.coords[0].dtype == np.int32

    def test_block_diag_idx_dtype(self):
        X = coo_array([[1, 0, 0], [0, 1, 0], [0, 1, 0]])
        X.coords = tuple(co.astype(np.int64) for co in X.coords)
        assert construct.block_diag([X, X]).coords[0].dtype == np.int64

    def test_block_diag_scalar_1d_args(self):
        """ block_diag with scalar and 1d arguments """
        # one 1d matrix and a scalar
        assert_array_equal(construct.block_diag([[2, 3], 4]).toarray(),
                           [[2, 3, 0], [0, 0, 4]])
        # 1d sparse arrays
        A = coo_array([1,0,3])
        B = coo_array([0,4])
        assert_array_equal(construct.block_diag([A, B]).toarray(),
                           [[1, 0, 3, 0, 0], [0, 0, 0, 0, 4]])

    def test_block_diag_1(self):
        """ block_diag with one matrix """
        assert_equal(construct.block_diag([[1, 0]]).toarray(),
                     array([[1, 0]]))
        assert_equal(construct.block_diag([[[1, 0]]]).toarray(),
                     array([[1, 0]]))
        assert_equal(construct.block_diag([[[1], [0]]]).toarray(),
                     array([[1], [0]]))
        # just on scalar
        assert_equal(construct.block_diag([1]).toarray(),
                     array([[1]]))

    def test_block_diag_sparse_arrays(self):
        """ block_diag with sparse arrays """

        A = coo_array([[1, 2, 3]], shape=(1, 3))
        B = coo_array([[4, 5]], shape=(1, 2))
        assert_equal(construct.block_diag([A, B]).toarray(),
                     array([[1, 2, 3, 0, 0], [0, 0, 0, 4, 5]]))

        A = coo_array([[1], [2], [3]], shape=(3, 1))
        B = coo_array([[4], [5]], shape=(2, 1))
        assert_equal(construct.block_diag([A, B]).toarray(),
                     array([[1, 0], [2, 0], [3, 0], [0, 4], [0, 5]]))

    def test_random_sampling(self):
        # Simple sanity checks for sparse random sampling.
        for f in sprand, _sprandn:
            for t in [np.float32, np.float64, np.longdouble,
                      np.int32, np.int64, np.complex64, np.complex128]:
                x = f((5, 10), density=0.1, dtype=t)
                assert_equal(x.dtype, t)
                assert_equal(x.shape, (5, 10))
                assert_equal(x.nnz, 5)

            x1 = f((5, 10), density=0.1, rng=4321)
            assert_equal(x1.dtype, np.float64)

            x2 = f((5, 10), density=0.1, rng=np.random.default_rng(4321))

            assert_array_equal(x1.data, x2.data)
            assert_array_equal(x1.row, x2.row)
            assert_array_equal(x1.col, x2.col)

            for density in [0.0, 0.1, 0.5, 1.0]:
                x = f((5, 10), density=density)
                assert_equal(x.nnz, int(density * np.prod(x.shape)))

            for fmt in ['coo', 'csc', 'csr', 'lil']:
                x = f((5, 10), format=fmt)
                assert_equal(x.format, fmt)

            assert_raises(ValueError, lambda: f((5, 10), density=1.1))
            assert_raises(ValueError, lambda: f((5, 10), density=-0.1))

    @pytest.mark.parametrize("rng", [None, 4321, np.random.default_rng(4321)])
    def test_rand(self, rng):
        # Simple distributional checks for sparse.rand.
        x = sprand((10, 20), density=0.5, dtype=np.float64, rng=rng)
        assert_(np.all(np.less_equal(0, x.data)))
        assert_(np.all(np.less_equal(x.data, 1)))

    @pytest.mark.parametrize("rng", [None, 4321, np.random.default_rng(4321)])
    def test_randn(self, rng):
        # Simple distributional checks for sparse.randn.
        # Statistically, some of these should be negative
        # and some should be greater than 1.
        x = _sprandn((10, 20), density=0.5, dtype=np.float64, rng=rng)
        assert_(np.any(np.less(x.data, 0)))
        assert_(np.any(np.less(1, x.data)))

    def test_random_accept_str_dtype(self):
        # anything that np.dtype can convert to a dtype should be accepted
        # for the dtype
        construct.random_array((10, 10), dtype='d')
        construct.random_array((10, 10, 10), dtype='d')
        construct.random_array((10, 10, 10, 10, 10), dtype='d')

    def test_random_array_maintains_array_shape(self):
        # preserve use of old random_state during SPEC 7 transition
        arr = construct.random_array((0, 4), density=0.3, dtype=int, random_state=0)
        assert arr.shape == (0, 4)

        arr = construct.random_array((10, 10, 10), density=0.3, dtype=int, rng=0)
        assert arr.shape == (10, 10, 10)

        arr = construct.random_array((10, 10, 10, 10, 10), density=0.3, dtype=int,
                                     rng=0)
        assert arr.shape == (10, 10, 10, 10, 10)

    def test_random_array_idx_dtype(self):
        A = construct.random_array((10, 10))
        assert A.coords[0].dtype == np.int32

    def test_random_sparse_matrix_returns_correct_number_of_non_zero_elements(self):
        # A 10 x 10 matrix, with density of 12.65%, should have 13 nonzero elements.
        # 10 x 10 x 0.1265 = 12.65, which should be rounded up to 13, not 12.
        sparse_array = construct.random_array((10, 10), density=0.1265)
        assert_equal(sparse_array.count_nonzero(),13)
        assert isinstance(sparse_array, sparray)
        # check big size
        shape = (2**33, 2**33)
        sparse_array = construct.random_array(shape, density=2.7105e-17)
        assert_equal(sparse_array.count_nonzero(),2000)

        # for n-D
        # check random_array
        sparse_array = construct.random_array((10, 10, 10, 10), density=0.12658)
        assert_equal(sparse_array.count_nonzero(),1266)
        assert isinstance(sparse_array, sparray)
        # check big size
        shape = (2**33, 2**33, 2**33)
        sparse_array = construct.random_array(shape, density=2.7105e-28)
        assert_equal(sparse_array.count_nonzero(),172)


def test_diags_array():
    """Tests of diags_array that do not rely on diags wrapper."""
    diag = np.arange(1.0, 5.0)

    assert_array_equal(construct.diags_array(diag, dtype=None).toarray(), np.diag(diag))

    assert_array_equal(
        construct.diags_array(diag, offsets=2, dtype=None).toarray(), np.diag(diag, k=2)
    )

    assert_array_equal(
        construct.diags_array(diag, offsets=2, shape=(4, 4), dtype=None).toarray(),
        np.diag(diag, k=2)[:4, :4]
    )

    # Offset outside bounds when shape specified
    with pytest.raises(ValueError, match=".*out of bounds"):
        construct.diags_array(np.arange(1.0, 5.0), offsets=5, shape=(4, 4))


def test_diags_int():
    d = [[3], [1, 2], [4]]
    offsets = [-1, 0, 1]
    # Until the deprecation period is over, `dtype=None` must be given
    # explicitly to avoid the warning and the cast to an inexact type
    # in diags_array() (gh-23102).
    arr = construct.diags_array(d, offsets=offsets, dtype=None)
    expected = np.array([[1, 4], [3, 2]])
    assert_array_equal(arr.toarray(), expected, strict=True)


def test_diags_int_to_float64():
    d = [[3], [1, 2], [4]]
    offsets = [-1, 0, 1]
    # Until the deprecation period is over, diags and diag_array will cast
    # integer inputs to float64 by default.  A warning will be generated
    # that indicates this behavior is deprecated.
    # See gh-23102.
    with pytest.warns(FutureWarning, match="output has been cast to"):
        arr = construct.diags_array(d, offsets=offsets)
    expected = np.array([[1.0, 4.0], [3.0, 2.0]])
    assert_array_equal(arr.toarray(), expected, strict=True)


def test_swapaxes():
    # Borrowed from Numpy swapaxes tests
    x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                  8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                  3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                  6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                  7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                  7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
    sX = coo_array(x).reshape(6, 6)
    sXswapped = construct.swapaxes(sX, 0, 1)
    assert_equal(sXswapped[-1].toarray(), sX[:, -1].toarray())

    sXX = sX.reshape(3, 2, 2, 3)
    sXXswapped = construct.swapaxes(sXX, 0, 2)
    assert_equal(sXXswapped.shape, (2, 2, 3, 3))


def test_3d_swapaxes():
    tgt = [[[0, 0], [2, 6]], [[1, 5], [0, 7]]]
    x = np.array([[[0, 1], [2, 0]], [[0, 5], [6, 7]]])
    A = coo_array(x) #[[[0, 1], [2, 0]], [[0, 5], [6, 7]]])
    out = construct.swapaxes(A, 0, 2)
    assert_equal(out.toarray(), tgt)
    assert_equal(out.toarray(), np.swapaxes(x, 0, 2))


@pytest.mark.parametrize("format", sparse_formats)
def test_sparse_format_swapaxes(format):
    A = np.array([[2, 0, 1], [3, 5, 0]])
    SA = coo_array(A).asformat(format)

    out = construct.swapaxes(SA, 1, 0)
    assert out.format == "coo"
    assert out.shape == (3, 2)
    assert_equal(out.toarray(), np.swapaxes(A, 1, 0))
    assert not out.has_canonical_format


def test_axis_swapaxes():
    A = coo_array([[2, 0], [3, 5]])
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, -4, 0)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, 0, -4)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, 3, 0)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, 0, 3)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, 1.2, 1)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.swapaxes(A, 1, 1.2)
    assert_equal(construct.swapaxes(A, 0, 0).toarray(), A.toarray())
    for i in range(2):
        assert_equal(
            construct.swapaxes(A, i, 1 - i).toarray(),
            construct.swapaxes(A, i - 2, -1 - i).toarray()
        )


def test_permute_dims():
    # Borrowed from Numpy tests.
    x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                  8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                  3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                  6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                  7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                  7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
    npx = x.reshape(6, 6)
    sX = coo_array(x).reshape(6, 6)

    sXpermuted = construct.permute_dims(sX, axes=(1, 0), copy=True)
    sXtransposed = sX.transpose(axes=(1, 0))
    assert_equal(sXpermuted.toarray(), sXtransposed.toarray())
    assert_equal(sXpermuted[-1].toarray(), sX[:, -1].toarray())

    npxx = npx.reshape(3, 2, 2, 3)
    sXX = sX.reshape(3, 2, 2, 3)
    sXXpermuted = construct.permute_dims(sXX, axes=(0, 2, 1, 3), copy=True)
    assert_equal(sXXpermuted.shape, (3, 2, 2, 3))
    sXXtransposed = sXX.transpose(axes=(0, 2, 1, 3))
    assert_equal(sXXtransposed.shape, (3, 2, 2, 3))
    assert_equal(sXXpermuted.toarray(), sXXtransposed.toarray())
    # TODO change np.transpose to np.permute_dims when numpy 2 is min supported version
    assert_equal(sXXpermuted.toarray(), np.transpose(npxx, axes=(0, 2, 1, 3)))


def test_3d_permute_dims():
    tgt = [[[0], [2], [0], [6]], [[1], [0], [5], [7]]]
    x = np.array([[[0, 1], [2, 0], [0, 5], [6, 7]]])
    A = coo_array(x)

    out = construct.permute_dims(A, axes=(2, 1, 0))
    assert_equal(out.shape, (2, 4, 1))
    assert_equal(out.toarray(), tgt)
    # TODO change np.transpose to np.permute_dims when numpy 2 is min supported version
    assert_equal(out.toarray(), np.transpose(x, axes=(2, 1, 0)))


def test_canonical_format_permute_dims():
    A = coo_array([[2, 0, 1], [3, 5, 0]])
    # identity axes keep has_canoncial_format True after permute_dims.
    assert construct.permute_dims(A, axes=(0, 1)).has_canonical_format is True
    assert construct.permute_dims(A, axes=[0, 1]).has_canonical_format is True
    # order changes set has_canonical_format to False
    assert construct.permute_dims(A, axes=[1, 0]).has_canonical_format is False


def test_axis_permute_dims():
    A = coo_array([[2, 0, 1], [3, 5, 0]])

    with assert_raises(ValueError, match="Incorrect number of axes"):
        construct.permute_dims(A, axes=(2, 0, 1))
    with assert_raises(ValueError, match="duplicate value in axis"):
        construct.permute_dims(A, axes=(0, 0))
    with assert_raises(TypeError, match="axis must be an integer/tuple"):
        construct.permute_dims(A, axes={1, 0})

    with assert_raises(ValueError, match="axis out of range"):
        construct.permute_dims(A, axes=(-3, 0))
    with assert_raises(ValueError, match="axis out of range"):
        construct.permute_dims(A, axes=(0, -3))
    with assert_raises(ValueError, match="axis out of range"):
        construct.permute_dims(A, axes=(2, 0))
    with assert_raises(ValueError, match="axis out of range"):
        construct.permute_dims(A, axes=(0, 2))
    with assert_raises(TypeError, match="axis must be an integer"):
        construct.permute_dims(A, axes=(1.2, 0))

    assert_equal(
        construct.permute_dims(A, axes=(1, 0), copy=True).toarray(),
        A.transpose(axes=(1, 0), copy=True).toarray()
    )
    # use lists for axes
    assert_equal(
        construct.permute_dims(A, axes=[1, 0], copy=True).toarray(),
        A.transpose(axes=[1, 0], copy=True).toarray()
    )
    assert_equal(
        construct.permute_dims(A, axes=None, copy=True).toarray(),
        A.transpose(axes=(1, 0), copy=True).toarray()
    )
    assert_equal(
        construct.permute_dims(A, axes=(0, 1), copy=True).toarray(), A.toarray()
    )


@pytest.mark.parametrize("format", sparse_formats)
def test_sparse_format_permute_dims(format):
    A = np.array([[2, 0, 1], [3, 5, 0]])
    SA = coo_array(A).asformat(format)

    out = construct.permute_dims(SA, axes=(1, 0))
    assert out.format == "coo"
    assert out.shape == (3, 2)
    # TODO change np.transpose to np.permute_dims when numpy 2 is min supported version
    assert_equal(out.toarray(), np.transpose(A, axes=(1, 0)))
    assert not out.has_canonical_format


def test_expand_dims():
    # Borrowed from Numpy tests.
    x = np.array([8.375, 7.545, 8.828, 8.5, 1.757, 5.928,
                  8.43, 7.78, 9.865, 5.878, 8.979, 4.732,
                  3.012, 6.022, 5.095, 3.116, 5.238, 3.957,
                  6.04, 9.63, 7.712, 3.382, 4.489, 6.479,
                  7.189, 9.645, 5.395, 4.961, 9.894, 2.893,
                  7.357, 9.828, 6.272, 3.758, 6.693, 0.993])
    npx = x.reshape(6, 6)
    sX = coo_array(npx)

    npx_expanded = np.expand_dims(npx, axis=1)
    sXexpanded = construct.expand_dims(sX, axis=1)
    assert_equal(sXexpanded[-1].toarray(), sX[-1, np.newaxis, :].toarray())
    assert_equal(sXexpanded.toarray(), npx_expanded)

    npxx = npx.reshape(3, 2, 2, 3)
    sXX = sX.reshape(3, 2, 2, 3)

    npxx_expanded = np.expand_dims(npxx, axis=2)
    sXXexpanded = construct.expand_dims(sXX, axis=2)
    assert_equal(sXXexpanded.shape, (3, 2, 1, 2, 3))
    assert_equal(sXXexpanded.toarray(), npxx_expanded)

    npxx_expanded = np.expand_dims(npxx, axis=-2)
    sXXexpanded = construct.expand_dims(sXX, axis=-2)
    assert_equal(sXXexpanded.shape, (3, 2, 2, 1, 3))
    assert_equal(sXXexpanded.toarray(), npxx_expanded)


def test_3d_expand_dims():
    tgt = [[[[0, 0], [2, 6]]], [[[1, 5], [0, 7]]]]
    A = coo_array([[[0, 0], [2, 6]], [[1, 5], [0, 7]]])
    out = construct.expand_dims(A, axis=1)
    assert_equal(out.toarray(), tgt)


@pytest.mark.parametrize("format", sparse_formats)
def test_sparse_format_expand_dims(format):
    A = np.array([[2, 0], [3, 5]])
    SA = coo_array(A).asformat(format)

    out = construct.expand_dims(SA, axis=1)
    assert out.format == "coo"
    assert out.shape == (2, 1, 2)
    assert_equal(out.toarray(), np.expand_dims(A, axis=1))
    assert SA.tocoo().has_canonical_format == out.has_canonical_format


def test_axis_expand_dims():
    A = coo_array([[2, 0], [3, 5]])
    with assert_raises(ValueError, match="Invalid axis"):
        construct.expand_dims(A, axis=-4)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.expand_dims(A, axis=3)
    with assert_raises(ValueError, match="Invalid axis"):
        construct.expand_dims(A, axis=1.2)
    for i in range(3):
        assert_equal(
            construct.expand_dims(A, axis=i).toarray(),
            construct.expand_dims(A, axis=i - 3).toarray()
        )
