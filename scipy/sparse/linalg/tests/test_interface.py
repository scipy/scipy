"""Test functions for the sparse.linalg._interface module
"""

from functools import partial
from itertools import product
import operator
from typing import NamedTuple
import pytest
from pytest import raises as assert_raises, warns
from numpy.testing import assert_, assert_equal, assert_allclose

import numpy as np
import scipy.sparse as sparse

import scipy.sparse.linalg._interface as interface
from scipy.sparse._sputils import matrix
from scipy._lib._gcutils import assert_deallocated, IS_PYPY
from scipy._lib._util import np_vecdot


class TestLinearOperator:
    def setup_method(self):
        self.A = np.array([[1,2,3],
                           [4,5,6]])
        self.B = np.array([[1,2],
                           [3,4],
                           [5,6]])
        self.C = np.array([[1,2],
                           [3,4]])

    def test_matvec(self):
        def get_matvecs(A):
            return [{
                        'shape': A.shape,
                        'matvec': lambda x: np.dot(A, x).reshape(A.shape[0]),
                        'rmatvec': lambda x: np.dot(A.T.conj(),
                                                    x).reshape(A.shape[1])
                    },
                    {
                        'shape': A.shape,
                        'matvec': lambda x: np.dot(A, x),
                        'rmatvec': lambda x: np.dot(A.T.conj(), x),
                        'rmatmat': lambda x: np.dot(A.T.conj(), x),
                        'matmat': lambda x: np.dot(A, x)
                    }]

        for matvecs in get_matvecs(self.A):
            A = interface.LinearOperator(**matvecs)

            assert_(A.args == ())

            assert_equal(A.matvec(np.array([1,2,3])), [14,32])
            assert_equal(A.matvec(np.array([[1],[2],[3]])), [[14],[32]])
            assert_equal(A @ np.array([1,2,3]), [14,32])
            assert_equal(A @ np.array([[1],[2],[3]]), [[14],[32]])
            assert_equal(A.dot(np.array([1,2,3])), [14,32])
            assert_equal(A.dot(np.array([[1],[2],[3]])), [[14],[32]])

            assert_equal(A.matvec(matrix([[1],[2],[3]])), [[14],[32]])
            assert_equal(A @ matrix([[1],[2],[3]]), [[14],[32]])
            assert_equal(A.dot(matrix([[1],[2],[3]])), [[14],[32]])

            assert_equal((2*A)@[1,1,1], [12,30])
            assert_equal((2 * A).rmatvec([1, 1]), [10, 14, 18])
            assert_equal((2*A).H.matvec([1,1]), [10, 14, 18])
            assert_equal((2*A).adjoint().matvec([1,1]), [10, 14, 18])
            assert_equal((2*A)@[[1],[1],[1]], [[12],[30]])
            assert_equal((2 * A).matmat([[1], [1], [1]]), [[12], [30]])
            assert_equal((A*2)@[1,1,1], [12,30])
            assert_equal((A*2)@[[1],[1],[1]], [[12],[30]])
            assert_equal((2j*A)@[1,1,1], [12j,30j])
            assert_equal((A+A)@[1,1,1], [12, 30])
            assert_equal((A + A).rmatvec([1, 1]), [10, 14, 18])
            assert_equal((A+A).H.matvec([1,1]), [10, 14, 18])
            assert_equal((A+A).adjoint().matvec([1,1]), [10, 14, 18])
            assert_equal((A+A)@[[1],[1],[1]], [[12], [30]])
            assert_equal((A+A).matmat([[1],[1],[1]]), [[12], [30]])
            assert_equal((-A)@[1,1,1], [-6,-15])
            assert_equal((-A)@[[1],[1],[1]], [[-6],[-15]])
            assert_equal((A-A)@[1,1,1], [0,0])
            assert_equal((A - A) @ [[1], [1], [1]], [[0], [0]])

            X = np.array([[1, 2], [3, 4]])
            # A_asarray = np.array([[1, 2, 3], [4, 5, 6]])
            assert_equal((2 * A).rmatmat(X), np.dot((2 * self.A).T, X))
            assert_equal((A * 2).rmatmat(X), np.dot((self.A * 2).T, X))
            assert_equal((2j * A).rmatmat(X),
                         np.dot((2j * self.A).T.conj(), X))
            assert_equal((A * 2j).rmatmat(X),
                         np.dot((self.A * 2j).T.conj(), X))
            assert_equal((A + A).rmatmat(X),
                         np.dot((self.A + self.A).T, X))
            assert_equal((A + 2j * A).rmatmat(X),
                         np.dot((self.A + 2j * self.A).T.conj(), X))
            assert_equal((-A).rmatmat(X), np.dot((-self.A).T, X))
            assert_equal((A - A).rmatmat(X),
                         np.dot((self.A - self.A).T, X))
            assert_equal((2j * A).rmatmat(2j * X),
                         np.dot((2j * self.A).T.conj(), 2j * X))

            z = A+A
            assert_(len(z.args) == 2 and z.args[0] is A and z.args[1] is A)
            z = 2*A
            assert_(len(z.args) == 2 and z.args[0] is A and z.args[1] == 2)

            assert_(isinstance(A.matvec([1, 2, 3]), np.ndarray))
            assert_(isinstance(A.matvec(np.array([[1],[2],[3]])), np.ndarray))
            assert_(isinstance(A @ np.array([1,2,3]), np.ndarray))
            assert_(isinstance(A @ np.array([[1],[2],[3]]), np.ndarray))
            assert_(isinstance(A.dot(np.array([1,2,3])), np.ndarray))
            assert_(isinstance(A.dot(np.array([[1],[2],[3]])), np.ndarray))

            assert_(isinstance(A.matvec(matrix([[1],[2],[3]])), np.ndarray))
            assert_(isinstance(A @ matrix([[1],[2],[3]]), np.ndarray))
            assert_(isinstance(A.dot(matrix([[1],[2],[3]])), np.ndarray))

            assert_(isinstance(2*A, interface._ScaledLinearOperator))
            assert_(isinstance(2j*A, interface._ScaledLinearOperator))
            assert_(isinstance(A+A, interface._SumLinearOperator))
            assert_(isinstance(-A, interface._ScaledLinearOperator))
            assert_(isinstance(A-A, interface._SumLinearOperator))
            assert_(isinstance(A/2, interface._ScaledLinearOperator))
            assert_(isinstance(A/2j, interface._ScaledLinearOperator))
            assert_(((A * 3) / 3).args[0] is A)  # check for simplification

            # Test that prefactor is of _ScaledLinearOperator is not mutated
            # when the operator is multiplied by a number
            result = A @ np.array([1, 2, 3])
            B = A * 3
            C = A / 5
            assert_equal(A @ np.array([1, 2, 3]), result)

            assert_((2j*A).dtype == np.complex128)

            # Test division by non-scalar
            msg = "Can only divide a linear operator by a scalar."
            with assert_raises(ValueError, match=msg):
                A / np.array([1, 2])

            assert_raises(ValueError, A.matvec, np.array([1,2]))
            assert_raises(ValueError, A.matvec, np.array([1,2,3,4]))
            assert_raises(ValueError, A.matvec, np.array([[1],[2]]))
            assert_raises(ValueError, A.matvec, np.array([[1],[2],[3],[4]]))

            assert_raises(ValueError, lambda: A@A)
            assert_raises(ValueError, lambda: A**2)

        for matvecsA, matvecsB in product(get_matvecs(self.A),
                                          get_matvecs(self.B)):
            A = interface.LinearOperator(**matvecsA)
            B = interface.LinearOperator(**matvecsB)
            # AtimesB = np.array([[22, 28], [49, 64]])
            AtimesB = self.A.dot(self.B)
            X = np.array([[1, 2], [3, 4]])

            assert_equal((A @ B).rmatmat(X), np.dot((AtimesB).T, X))
            assert_equal((2j * A @ B).rmatmat(X),
                         np.dot((2j * AtimesB).T.conj(), X))

            assert_equal((A@B)@[1,1], [50,113])
            assert_equal((A@B)@[[1],[1]], [[50],[113]])
            assert_equal((A@B).matmat([[1],[1]]), [[50],[113]])

            assert_equal((A @ B).rmatvec([1, 1]), [71, 92])
            assert_equal((A @ B).H.matvec([1, 1]), [71, 92])
            assert_equal((A @ B).adjoint().matvec([1, 1]), [71, 92])

            assert_(isinstance(A@B, interface._ProductLinearOperator))

            assert_raises(ValueError, lambda: A+B)
            assert_raises(ValueError, lambda: A**2)

            z = A@B
            assert_(len(z.args) == 2 and z.args[0] is A and z.args[1] is B)

        for matvecsC in get_matvecs(self.C):
            C = interface.LinearOperator(**matvecsC)
            X = np.array([[1, 2], [3, 4]])

            assert_equal(C.rmatmat(X), np.dot((self.C).T, X))
            assert_equal((C**2).rmatmat(X),
                         np.dot((np.dot(self.C, self.C)).T, X))

            assert_equal((C**2)@[1,1], [17,37])
            assert_equal((C**2).rmatvec([1, 1]), [22, 32])
            assert_equal((C**2).H.matvec([1, 1]), [22, 32])
            assert_equal((C**2).adjoint().matvec([1, 1]), [22, 32])
            assert_equal((C**2).matmat([[1],[1]]), [[17],[37]])

            assert_(isinstance(C**2, interface._PowerLinearOperator))

    def test_matmul(self):
        D = {'shape': self.A.shape,
             'matvec': lambda x: np.dot(self.A, x).reshape(self.A.shape[0]),
             'rmatvec': lambda x: np.dot(self.A.T.conj(),
                                         x).reshape(self.A.shape[1]),
             'rmatmat': lambda x: np.dot(self.A.T.conj(), x),
             'matmat': lambda x: np.dot(self.A, x)}
        A = interface.LinearOperator(**D)
        B = np.array([[1 + 1j, 2, 3],
                      [4, 5, 6],
                      [7, 8, 9]])
        b = B[0]

        assert_equal(operator.matmul(A, b), A * b)
        assert_equal(operator.matmul(A, b.reshape(-1, 1)), A * b.reshape(-1, 1))
        assert_equal(operator.matmul(A, B), A @ B)
        assert_equal(operator.matmul(b, A.H), b * A.H)
        assert_equal(operator.matmul(b, A.adjoint()), b * A.adjoint())
        assert_equal(operator.matmul(b.reshape(1, -1), A.H), b.reshape(1, -1) * A.H)
        assert_equal(operator.matmul(b.reshape(1, -1), A.adjoint()),
                     b.reshape(1, -1) * A.adjoint())
        assert_equal(operator.matmul(B, A.H), B @ A.H)
        assert_equal(operator.matmul(B, A.adjoint()), B @ A.adjoint())
        assert_raises(ValueError, operator.matmul, A, 2)
        assert_raises(ValueError, operator.matmul, 2, A)


class TestDotTests:
    """
    This class aims to help ensure correctness of the LinearOperator
    interface, by verifying correctness properties based on equivalent
    computations using 'forward' and 'adjoint' modes.
    """
    class OperatorArgs(NamedTuple):
        """
        shape: shape of the operator
        op_dtype: dtype of the operator
        data_dtype: real dtype corresponding to op_dtype for data generation
        complex: the operator has a complex dtype
        """
        shape: tuple[int, ...]
        op_dtype: str
        data_dtype: str
        complex: bool

    real_square_args: OperatorArgs = OperatorArgs(
        (12, 12), "float64", "float64", False
    )
    integer_square_args: OperatorArgs = OperatorArgs(
        (9, 9), "int32", "float32", False
    )
    complex_square_args: OperatorArgs = OperatorArgs(
        (13, 13), "complex64", "float32", True
    )
    real_overdetermined_args: OperatorArgs = OperatorArgs(
        (17, 11), "float64", "float64", False
    )
    complex_overdetermined_args: OperatorArgs = OperatorArgs(
        (17, 11), "complex128", "float64", True
    )
    real_underdetermined_args: OperatorArgs = OperatorArgs(
        (5, 9), "float64", "float64", False
    )

    square_args_list: list[OperatorArgs] = [
        real_square_args, integer_square_args, complex_square_args
    ]
    all_args_list: list[OperatorArgs]  = square_args_list + [
        real_overdetermined_args, complex_overdetermined_args, real_underdetermined_args
    ]

    def check_matvec(
        self, op: interface.LinearOperator, data_dtype: str, complex_data: bool = False,
        check_operators: bool = False, check_dot: bool = False
    ):
        """
        This check verifies the equivalence of the forward and adjoint computation,
        using `matvec` and `rmatvec` respectively, on randomised data.

        Data is generated with the real dtype `data_dtype` and operated on by the
        linear operator `op`.

        If `complex_data` is set to `True`, complex data is instead generated
        by combining randomised real and imaginary components, each generated
        with `data_dtype`.

        If `check_operators` is set to `True`, equivalence is checked between
        `matvec` and `*` and `@`,
        and between `rmatvec` and the composition of `.H` and `*`.

        If `check_dot` is set to `True`, equivalence is checked between
        `matvec` and `.dot`,
        and between `rmatvec` and the composition of `.H` and `.dot`.
        """
        rng = np.random.default_rng(42)

        u = rng.standard_normal(op.shape[-1], dtype=data_dtype)
        v = rng.standard_normal(op.shape[-2], dtype=data_dtype)
        if complex_data:
            u = u + (1j * rng.standard_normal(op.shape[-1], dtype=data_dtype))
            v = v + (1j * rng.standard_normal(op.shape[-2], dtype=data_dtype))

        op_u = op.matvec(u)
        opH_v = op.rmatvec(v)

        if check_operators:
            assert_allclose(op_u, op * u)
            assert_allclose(op_u, op @ u)
            assert_allclose(opH_v, op.H * v)
            assert_allclose(opH_v, op.H @ v)

        if check_dot:
            assert_allclose(op_u, op.dot(u))
            assert_allclose(opH_v, op.H.dot(v))

        op_u_H_v = np_vecdot(op_u, v, axis=-1)
        uH_opH_v = np_vecdot(u, opH_v, axis=-1)

        rtol = 1e-12 if np.finfo(data_dtype).eps < 1e-8 else 1e-5
        assert_allclose(op_u_H_v, uH_opH_v, rtol=rtol)

    def check_matmat(
        self, op: interface.LinearOperator, data_dtype: str, complex_data: bool = False,
        check_operators: bool = False, check_dot: bool = False
    ):
        """
        This check verifies the equivalence of the forward and adjoint computation,
        using `matmat` and `rmatmat` respectively, on randomised data.

        Data is generated with the real dtype `data_dtype` and operated on by the
        linear operator `op`.

        If `complex_data` is set to `True`, complex data is instead generated
        by combining randomised real and imaginary components, each generated
        with `data_dtype`.

        If `check_operators` is set to `True`, equivalence is checked between
        `matmat` and `*` and `@`,
        and between `rmatvec` and the composition of `.H` and `@`.

        If `check_dot` is set to `True`, equivalence is checked between
        `matmat` and `.dot`,
        and between `rmatmat` and the composition of `.H` and `.dot`.
        """
        rng = np.random.default_rng(42)
        k = rng.integers(2, 100)

        U = rng.standard_normal(size=(op.shape[-1], k), dtype=data_dtype)
        V = rng.standard_normal(size=(op.shape[-2], k), dtype=data_dtype)
        if complex_data:
            U = U + (1j * rng.standard_normal(size=(op.shape[-1], k), dtype=data_dtype))
            V = V + (1j * rng.standard_normal(size=(op.shape[-2], k), dtype=data_dtype))

        op_U = op.matmat(U)
        opH_V = op.rmatmat(V)

        if check_operators:
            assert_allclose(op_U, op * U)
            assert_allclose(op_U, op @ U)
            assert_allclose(opH_V, op.H * V)
            assert_allclose(opH_V, op.H @ V)

        if check_dot:
            assert_allclose(op_U, op.dot(U))
            assert_allclose(opH_V, op.H.dot(V))

        op_U_H = np.conj(op_U).T
        UH = np.conj(U).T

        op_U_H_V = np.matmul(op_U_H, V)
        UH_opH_V = np.matmul(UH, opH_V)

        rtol = 3e-12 if np.finfo(data_dtype).eps < 1e-8 else 6e-4
        assert_allclose(op_U_H_V, UH_opH_V, rtol=rtol)

    @pytest.mark.parametrize("args", square_args_list)
    def test_identity_square(self, args):
        """Simple identity operator on square matrices"""
        def identity(x):
            return x

        op = interface.LinearOperator(
            shape=args.shape, dtype=args.op_dtype,
            matvec=identity, rmatvec=identity
        )

        self.check_matvec(op, data_dtype=args.data_dtype, complex_data=args.complex)
        self.check_matmat(op, data_dtype=args.data_dtype, complex_data=args.complex)
    
    @pytest.mark.parametrize("args", all_args_list)
    def test_identity_nonsquare(self, args):
        """Identity operator with zero-padding on non-square matrices"""
        def mv(x):
            # handle column vectors too
            # (`LinearOperator` handles reshape in post-processing)
            x = x.flatten()

            match np.sign(x.shape[0] - args.shape[-2]):
                case 0:  # square
                    return x
                case 1:  # crop x to size
                    return x[:args.shape[-2]]
                case -1:  # pad with zeros
                    pad_width = (0, args.shape[-2] - x.shape[0])
                    return np.pad(x, pad_width, mode='constant', constant_values=0)

        def rmv(x):
            # handle column vectors too
            # (`LinearOperator` handles reshape in post-processing)
            x = x.flatten()
            
            match np.sign(args.shape[-1] - x.shape[0]):
                case 0:  # square
                    return x
                case 1:  # pad with zeros
                    pad_width = (0, args.shape[-1] - x.shape[0])
                    return np.pad(x, pad_width, mode='constant', constant_values=0)
                case -1:  # crop x to size
                    return x[:args.shape[-1]]

        op = interface.LinearOperator(
            shape=args.shape, dtype=args.op_dtype, matvec=mv, rmatvec=rmv
        )
        
        self.check_matvec(op, data_dtype=args.data_dtype, complex_data=args.complex)
        self.check_matmat(op, data_dtype=args.data_dtype, complex_data=args.complex)
        
    @pytest.mark.parametrize("args", square_args_list)
    def test_scaling_square(self, args):
        """Simple (complex) scaling operator on square matrices"""
        def scale(x):
            return (3 + 2j) * x

        def r_scale(x):
            return (3 - 2j) * x

        op = interface.LinearOperator(
            shape=args.shape, dtype=args.op_dtype, matvec=scale, rmatvec=r_scale
        )
        self.check_matvec(
            op, data_dtype=args.data_dtype, complex_data=args.complex,
            check_operators=True, check_dot=True
        )
        self.check_matmat(
            op, data_dtype=args.data_dtype, complex_data=args.complex,
            check_operators=True, check_dot=True
        )

    def test_subclass_matmat(self):
        """
        Simple rotation operator defined by `matmat` and `adjoint`,
        subclassing `LinearOperator`.
        """
        def rmatmat(X):
            theta = np.pi / 2
            R_inv = np.array([
                [np.cos(theta),  np.sin(theta)],
                [-np.sin(theta), np.cos(theta)]
            ])
            return R_inv @ X

        class RotOp(interface.LinearOperator):
            
            def __init__(self, dtype, shape, theta):
                self._theta = theta
                super().__init__(dtype, shape)
            
            def _matmat(self, X):
                theta = self._theta
                R = np.array([
                    [np.cos(theta), -np.sin(theta)],
                    [np.sin(theta),  np.cos(theta)]
                ])
                return R @ X

            def _adjoint(self):
                negative_theta = -self._theta
                return RotOp(self.dtype, self.shape, negative_theta)
            
        theta = np.pi / 2
        dtype = "float64"
        op = RotOp(shape=(2, 2), dtype=dtype, theta=theta)

        self.check_matvec(
            op, data_dtype=dtype, complex_data=False,
            check_operators=True, check_dot=True
        )
        self.check_matmat(
            op, data_dtype=dtype, complex_data=False,
            check_operators=True, check_dot=True
        )
    
    @pytest.mark.parametrize(
        "matrix", [
            np.asarray([[1, 2j, 3j], [4j, 5j, 6]]),
            sparse.random_array((5, 5))
        ]
    )
    def test_aslinearop(self, matrix):
        op = interface.aslinearoperator(matrix)
        data_dtype = "float64"
        self.check_matvec(op, data_dtype=data_dtype, complex_data=True)
        self.check_matmat(op, data_dtype=data_dtype, complex_data=True)


class TestAsLinearOperator:
    def setup_method(self):
        self.cases = []

        def make_cases(original, dtype):
            cases = []

            cases.append((matrix(original, dtype=dtype), original))
            cases.append((np.array(original, dtype=dtype), original))
            cases.append((sparse.csr_array(original, dtype=dtype), original))

            # Test default implementations of _adjoint and _rmatvec, which
            # refer to each other.
            def mv(x, dtype):
                y = original.dot(x)
                if len(x.shape) == 2:
                    y = y.reshape(-1, 1)
                return y

            def rmv(x, dtype):
                return original.T.conj().dot(x)

            class BaseMatlike(interface.LinearOperator):
                args = ()

                def __init__(self, dtype):
                    self.dtype = np.dtype(dtype)
                    self.shape = original.shape

                def _matvec(self, x):
                    return mv(x, self.dtype)

            class HasRmatvec(BaseMatlike):
                args = ()

                def _rmatvec(self,x):
                    return rmv(x, self.dtype)

            class HasAdjoint(BaseMatlike):
                args = ()

                def _adjoint(self):
                    shape = self.shape[1], self.shape[0]
                    matvec = partial(rmv, dtype=self.dtype)
                    rmatvec = partial(mv, dtype=self.dtype)
                    return interface.LinearOperator(matvec=matvec,
                                                    rmatvec=rmatvec,
                                                    dtype=self.dtype,
                                                    shape=shape)

            class HasRmatmat(HasRmatvec):
                def _matmat(self, x):
                    return original.dot(x)

                def _rmatmat(self, x):
                    return original.T.conj().dot(x)

            cases.append((HasRmatvec(dtype), original))
            cases.append((HasAdjoint(dtype), original))
            cases.append((HasRmatmat(dtype), original))
            return cases

        original = np.array([[1,2,3], [4,5,6]])
        self.cases += make_cases(original, np.int32)
        self.cases += make_cases(original, np.float32)
        self.cases += make_cases(original, np.float64)
        self.cases += [(interface.aslinearoperator(M).T, A.T)
                       for M, A in make_cases(original.T, np.float64)]
        self.cases += [(interface.aslinearoperator(M).H, A.T.conj())
                       for M, A in make_cases(original.T, np.float64)]
        self.cases += [(interface.aslinearoperator(M).adjoint(), A.T.conj())
                       for M, A in make_cases(original.T, np.float64)]

        original = np.array([[1, 2j, 3j], [4j, 5j, 6]])
        self.cases += make_cases(original, np.complex128)
        self.cases += [(interface.aslinearoperator(M).T, A.T)
                       for M, A in make_cases(original.T, np.complex128)]
        self.cases += [(interface.aslinearoperator(M).H, A.T.conj())
                       for M, A in make_cases(original.T, np.complex128)]
        self.cases += [(interface.aslinearoperator(M).adjoint(), A.T.conj())
                       for M, A in make_cases(original.T, np.complex128)]

    def test_basic(self):

        for M, A_array in self.cases:
            A = interface.aslinearoperator(M)
            M,N = A.shape

            xs = [np.array([1, 2, 3]),
                  np.array([[1], [2], [3]])]
            ys = [np.array([1, 2]), np.array([[1], [2]])]

            if A.dtype == np.complex128:
                xs += [np.array([1, 2j, 3j]),
                       np.array([[1], [2j], [3j]])]
                ys += [np.array([1, 2j]), np.array([[1], [2j]])]

            x2 = np.array([[1, 4], [2, 5], [3, 6]])

            for x in xs:
                assert_equal(A.matvec(x), A_array.dot(x))
                assert_equal(A @ x, A_array.dot(x))

            assert_equal(A.matmat(x2), A_array.dot(x2))
            assert_equal(A @ x2, A_array.dot(x2))

            for y in ys:
                assert_equal(A.rmatvec(y), A_array.T.conj().dot(y))
                assert_equal(A.T.matvec(y), A_array.T.dot(y))
                assert_equal(A.H.matvec(y), A_array.T.conj().dot(y))
                assert_equal(A.adjoint().matvec(y), A_array.T.conj().dot(y))

            for y in ys:
                if y.ndim < 2:
                    continue
                assert_equal(A.rmatmat(y), A_array.T.conj().dot(y))
                assert_equal(A.T.matmat(y), A_array.T.dot(y))
                assert_equal(A.H.matmat(y), A_array.T.conj().dot(y))
                assert_equal(A.adjoint().matmat(y), A_array.T.conj().dot(y))

            if hasattr(M,'dtype'):
                assert_equal(A.dtype, M.dtype)

            assert_(hasattr(A, 'args'))

    def test_dot(self):

        for M, A_array in self.cases:
            A = interface.aslinearoperator(M)
            M,N = A.shape

            x0 = np.array([1, 2, 3])
            x1 = np.array([[1], [2], [3]])
            x2 = np.array([[1, 4], [2, 5], [3, 6]])

            assert_equal(A.dot(x0), A_array.dot(x0))
            assert_equal(A.dot(x1), A_array.dot(x1))
            assert_equal(A.dot(x2), A_array.dot(x2))


def test_repr():
    A = interface.LinearOperator(shape=(1, 1), matvec=lambda x: 1)
    repr_A = repr(A)
    assert_('unspecified dtype' not in repr_A, repr_A)


def test_identity():
    ident = interface.IdentityOperator((3, 3))
    assert_equal(ident @ [1, 2, 3], [1, 2, 3])
    assert_equal(ident.dot(np.arange(9).reshape(3, 3)).ravel(), np.arange(9))

    assert_raises(ValueError, ident.matvec, [1, 2, 3, 4])


def test_attributes():
    A = interface.aslinearoperator(np.arange(16).reshape(4, 4))

    def always_four_ones(x):
        x = np.asarray(x)
        assert_(x.shape == (3,) or x.shape == (3, 1))
        return np.ones(4)

    B = interface.LinearOperator(shape=(4, 3), matvec=always_four_ones)

    ops = [A, B, A * B, A @ B, A.H, A.adjoint(), A + A, B + B, A**4]
    for op in ops:
        assert_(hasattr(op, "dtype"))
        assert_(hasattr(op, "shape"))
        assert_(hasattr(op, "_matvec"))

def matvec(x):
    """ Needed for test_pickle as local functions are not pickleable """
    return np.zeros(3)

def test_pickle():
    import pickle

    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        A = interface.LinearOperator((3, 3), matvec)
        s = pickle.dumps(A, protocol=protocol)
        B = pickle.loads(s)

        for k in A.__dict__:
            assert_equal(getattr(A, k), getattr(B, k))


def test_inheritance():
    class Empty(interface.LinearOperator):
        pass

    with warns(RuntimeWarning, match="should implement at least"):
        assert_raises(TypeError, Empty)

    class Identity(interface.LinearOperator):
        def __init__(self, n):
            super().__init__(dtype=None, shape=(n, n))

        def _matvec(self, x):
            return x

    id3 = Identity(3)
    assert_equal(id3.matvec([1, 2, 3]), [1, 2, 3])
    assert_raises(NotImplementedError, id3.rmatvec, [4, 5, 6])

    class MatmatOnly(interface.LinearOperator):
        def __init__(self, A):
            super().__init__(A.dtype, A.shape)
            self.A = A

        def _matmat(self, x):
            return self.A.dot(x)

    mm = MatmatOnly(np.random.randn(5, 3))
    assert_equal(mm.matvec(np.random.randn(3)).shape, (5,))

def test_dtypes_of_operator_sum():
    # gh-6078

    mat_complex = np.random.rand(2,2) + 1j * np.random.rand(2,2)
    mat_real = np.random.rand(2,2)

    complex_operator = interface.aslinearoperator(mat_complex)
    real_operator = interface.aslinearoperator(mat_real)

    sum_complex = complex_operator + complex_operator
    sum_real = real_operator + real_operator

    assert_equal(sum_real.dtype, np.float64)
    assert_equal(sum_complex.dtype, np.complex128)

def test_no_double_init():
    call_count = [0]

    def matvec(v):
        call_count[0] += 1
        return v

    # It should call matvec exactly once (in order to determine the
    # operator dtype)
    interface.LinearOperator((2, 2), matvec=matvec)
    assert_equal(call_count[0], 1)

INT_DTYPES = (np.int8, np.int16, np.int32, np.int64)
REAL_DTYPES = (np.float32, np.float64, np.longdouble)
COMPLEX_DTYPES = (np.complex64, np.complex128, np.clongdouble)
INEXACTDTYPES = REAL_DTYPES + COMPLEX_DTYPES
ALLDTYPES = INT_DTYPES + INEXACTDTYPES


@pytest.mark.parametrize("test_dtype", ALLDTYPES)
def test_determine_lo_dtype_from_matvec(test_dtype):
    # gh-19209
    scalar = np.array(1, dtype=test_dtype)
    def mv(v):
        return np.array([scalar * v[0], v[1]])

    lo = interface.LinearOperator((2, 2), matvec=mv)
    assert lo.dtype == np.dtype(test_dtype)

def test_determine_lo_dtype_for_int():
    # gh-19209
    # test Python int larger than int8 max cast to some int
    def mv(v):
        return np.array([128 * v[0], v[1]])

    lo = interface.LinearOperator((2, 2), matvec=mv)
    assert lo.dtype in INT_DTYPES

def test_adjoint_conjugate():
    X = np.array([[1j]])
    A = interface.aslinearoperator(X)

    B = 1j * A
    Y = 1j * X

    v = np.array([1])

    assert_equal(B.dot(v), Y.dot(v))
    assert_equal(B.H.dot(v), Y.T.conj().dot(v))
    assert_equal(B.adjoint().dot(v), Y.T.conj().dot(v))

def test_ndim():
    X = np.array([[1]])
    A = interface.aslinearoperator(X)
    assert_equal(A.ndim, 2)

def test_transpose_noconjugate():
    X = np.array([[1j]])
    A = interface.aslinearoperator(X)

    B = 1j * A
    Y = 1j * X

    v = np.array([1])

    assert_equal(B.dot(v), Y.dot(v))
    assert_equal(B.T.dot(v), Y.T.dot(v))

def test_transpose_multiplication():
    class MyMatrix(interface.LinearOperator):
        def __init__(self, A):
            super().__init__(A.dtype, A.shape)
            self.A = A
        def _matmat(self, other): return self.A @ other
        def _rmatmat(self, other): return self.A.T @ other

    A = MyMatrix(np.array([[1, 2], [3, 4]]))
    X = np.array([1, 2])
    B = np.array([[10, 20], [30, 40]])
    X2 = X.reshape(-1, 1)
    Y = np.array([[1, 2], [3, 4]])

    assert_equal(A @ B, Y @ B)
    assert_equal(B.T @ A, B.T @ Y)
    assert_equal(A.T @ B, Y.T @ B)
    assert_equal(A @ X, Y @ X)
    assert_equal(X.T @ A, X.T @ Y)
    assert_equal(A.T @ X, Y.T @ X)
    assert_equal(A @ X2, Y @ X2)
    assert_equal(X2.T @ A, X2.T @ Y)
    assert_equal(A.T @ X2, Y.T @ X2)

def test_sparse_matmat_exception():
    A = interface.LinearOperator((2, 2), matvec=lambda x: x)
    B = sparse.eye_array(2)
    msg = "Unable to multiply a LinearOperator with a sparse matrix."
    with assert_raises(TypeError, match=msg):
        A @ B
    with assert_raises(TypeError, match=msg):
        B @ A
    with assert_raises(ValueError):
        A @ np.identity(4)
    with assert_raises(ValueError):
        np.identity(4) @ A


@pytest.mark.skipif(IS_PYPY, reason="Test not meaningful on PyPy")
def test_MatrixLinearOperator_refcycle():
    # gh-10634
    # Test that MatrixLinearOperator can be automatically garbage collected
    A = np.eye(2)
    with assert_deallocated(interface.MatrixLinearOperator, A) as op:
        op.adjoint()
        del op
