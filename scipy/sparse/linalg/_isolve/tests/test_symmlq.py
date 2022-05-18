import numpy as np
from numpy.testing import assert_equal, assert_
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg._isolve import symmlq


def get_sample_problem():
    matrix = csr_matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], dtype=float)
    rhs = np.ones(matrix.shape[0])
    return matrix, rhs


count = [0]
residuals = []


def cb(x):
    count[0] += 1
    residuals.append(x)


class TestSYMMLQ():
    def test_happy_breakdown(self, capsys):
        # Test initial happy breakdown
        A, b = get_sample_problem()
        eps = 1e-10
        x0 = np.array([1.5, 2., 1.5+eps])

        x, info = symmlq(A, b, x0=x0, tol=eps, callback=cb, show=True)
        out, err = capsys.readouterr()
        assert_equal(residuals[0], np.linalg.norm(b - A@x0))
        assert_equal(out, "SYMMLQ: Linear solve converged due to HAPPY "
                          "BREAKDOWN iterations 0\n")
        assert_equal(err, '')

        # Test happy breakdown in iterations
        eps = 1e-4
        x0[2] = 1.5 + eps
        count[0] = 0
        residuals.clear()
        x, info = symmlq(A, b, x0=x0, callback=cb, show=True)
        out, err = capsys.readouterr()
        assert_(count[0] > 0)
        assert_equal(residuals[len(residuals)-1], np.linalg.norm(b - A@x))
        assert_equal(out, "SYMMLQ: Linear solve converged due to reach TOL "
                          f"iterations {count[0]}\n")
        assert_equal(err, '')

    def test_diverged_indefinite_pc(self, capsys):
        # Test diverged indefinite preconditioner under initial step
        A, b = get_sample_problem()
        eps = 1e-5
        x0 = np.array([1.5, 2., 1.5+eps])

        def matvec(x):
            ret = np.array([x[0], -x[1], -x[2]])
            return ret
        M = LinearOperator(A.shape, matvec=matvec, dtype=A.dtype)

        count[0] = 0
        residuals.clear()
        x, info = symmlq(A, b, x0=x0, M=M, callback=cb, show=True)
        out, err = capsys.readouterr()
        assert_equal(residuals[0], np.linalg.norm(M@(b - A@x0)))
        assert_equal(out, "SYMMLQ: Linear solve not converged due to DIVERGED "
                          "INDEFINITE PC iterations 0\n")
        assert_equal(err, '')

        # Test diverged indefinite preconditioner in iterations
        def matvec(x):
            ret = np.array([x[0], -x[1], x[2]])
            return ret
        M = LinearOperator(A.shape, matvec=matvec, dtype=A.dtype)

        count[0] = 0
        residuals.clear()
        x, info = symmlq(A, b, x0=x0, M=M, callback=cb, show=True)
        out, err = capsys.readouterr()
        assert_(count[0] > 0)
        assert_equal(residuals[len(residuals)-1], np.linalg.norm(M@(b - A@x)))
        assert_equal(out, "SYMMLQ: Linear solve not converged due to DIVERGED "
                          f"INDEFINITE PC iterations {count[0]}\n")
        assert_equal(err, '')
