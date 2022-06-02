import pytest
import numpy as np
from numpy.testing import assert_equal, assert_array_equal, assert_
import scipy.io as io
from scipy.sparse.linalg import LinearOperator, splu
from scipy.sparse.linalg._isolve import symmlq
from .test_iterative import assert_normclose


# Build CFD problem (real symmetric indefinite system)
def get_sample_real_problem():
    problem = "Harwell-Boeing/sherman/sherman1"
    mm = np.lib._datasource.Repository('https://math.nist.gov/pub/MatrixMarket2/')
    f = mm.open(f"{problem}.mtx.gz")
    matrix = io.mmread(f).tocsc()
    f.close()

    f = mm.open(f'{problem}_rhs1.mtx.gz')
    rhs = np.array(io.mmread(f)).ravel()
    f.close()
    return matrix, rhs


# Build MHD problem (Hermitian indefinite system)
def get_sample_complex_problem():
    problem = "NEP/mhd/mhd1280b"
    mm = np.lib._datasource.Repository('https://math.nist.gov/pub/MatrixMarket2/')
    f = mm.open(f"{problem}.mtx.gz")
    matrix = io.mmread(f).tocsc()
    f.close()
    rng = np.random.default_rng(1234)
    n = matrix.shape[0]
    rhs = rng.standard_normal(n) + rng.standard_normal(n)*1j
    return matrix, rhs


A, b = get_sample_real_problem()
Ac, bc = get_sample_complex_problem()
count = [0]
residuals = []


class TestSYMMLQ():
    def test_happy_breakdown(self, capsys):
        def cb(x):
            residuals.append(np.linalg.norm(b - A@x))

        eps = 1e-8
        rtol = 1e-12
        x0 = splu(A).solve(b)
        x0[0] += eps
        x, info = symmlq(A, b, x0=x0, tol=rtol, callback=cb, verbose=True)
        out, err = capsys.readouterr()
        assert_equal(residuals[0], np.linalg.norm(b - A@x))
        assert_equal(out, "SYMMLQ: Linear solve converged due to HAPPY "
                          "BREAKDOWN iterations 0\n")
        assert_equal(err, '')

    def test_diverged_indefinite_pc(self, capsys):
        def cb(x):
            count[0] += 1
            residuals.append(np.linalg.norm(b - A@x))

        # Test diverged indefinite preconditioner under initial step
        residuals.clear()
        pc = splu(A)
        M = LinearOperator(shape=A.shape, matvec=pc.solve, dtype=A.dtype)
        x, info = symmlq(A, b, M=M, callback=cb, verbose=True)
        out, err = capsys.readouterr()
        assert_equal(residuals[len(residuals)-1], np.linalg.norm(b - A@x))
        assert_equal(out, "SYMMLQ: Linear solve not converged due to DIVERGED "
                          "INDEFINITE PC iterations 0\n")
        assert_equal(err, '')

        # Test diverged indefinite preconditioner in iterations
        count[0] = 0
        residuals.clear()
        def matvec(x):
            D = np.ones(A.shape[0])
            D[1] = -1.
            return D * x
        M = LinearOperator(A.shape, matvec=matvec, dtype=A.dtype)
        x, info = symmlq(A, b, M=M, callback=cb, verbose=True)
        out, err = capsys.readouterr()
        assert_equal(residuals[len(residuals)-1], np.linalg.norm(b - A@x))
        assert_equal(out, "SYMMLQ: Linear solve not converged due to DIVERGED "
                          f"INDEFINITE PC iterations {count[0]}\n")
        assert_equal(err, '')

    @pytest.mark.parametrize(('mat', 'rhs'), [(A, b), (Ac, bc)])
    def test_convergence(self, mat, rhs):
        def cb(x):
            count[0] += 1

        # Unpreconditioned SYMMLQ
        rtol = 1e-8
        count[0] = 0
        # Converged only for CFD, appropriate preconditioner need to
        # be constructed for MHD
        if not np.iscomplex(rhs).any():
            x, info = symmlq(mat, rhs, tol=rtol, callback=cb)
            assert_equal(info, 0)
            assert_normclose(mat@x, rhs, tol=rtol)

        # Preconditioned SYMMLQ
        countWithoutPre = count[0]
        count[0] = 0
        def matvec(x):
            invD = 1. / abs(mat.diagonal())
            return invD * x
        M = LinearOperator(mat.shape, matvec=matvec, dtype=mat.dtype)
        x, info = symmlq(mat, rhs, tol=rtol, M=M, callback=cb)
        assert_equal(info, 0)
        assert_normclose(mat@x, rhs, tol=rtol)
        # Number of iterations of SYMMLQ with preconditioner should not
        # be more than without preconditioner
        if not np.iscomplex(rhs).any():
            assert_(count[0] <= countWithoutPre)

    @pytest.mark.parametrize(('mat', 'rhs'), [(A, b), (Ac, bc)])
    def test_non_default_x0(self, capsys, mat, rhs):
        # Unpreconditioned SYMMLQ
        rtol = 1e-8
        rng = np.random.default_rng(12345)
        x0 = rng.standard_normal(mat.shape[0])
        # Converged only for CFD, appropriate preconditioner need to
        # be constructed for MHD
        if not np.iscomplex(rhs).any():
            x, info = symmlq(mat, rhs, x0=x0, tol=rtol)
            assert_equal(info, 0)
            assert_normclose(mat@x, rhs, tol=rtol)

        # Preconditioned SYMMLQ
        def matvec(x):
            invD = 1. / abs(mat.diagonal())
            return invD * x
        M = LinearOperator(mat.shape, matvec=matvec, dtype=mat.dtype)
        # Use random initial guess
        x, info = symmlq(mat, rhs, x0=x0, tol=rtol, M=M)
        assert_equal(info, 0)
        assert_normclose(mat@x, rhs, tol=rtol)
        # Use Knoll's initial guess (preconditioned RHS)
        x, info = symmlq(mat, rhs, x0='Mb', tol=rtol, M=M)
        assert_equal(info, 0)
        assert_normclose(mat@x, rhs, tol=rtol)

    @pytest.mark.parametrize(('mat', 'rhs'), [(A, b), (Ac, bc)])
    def test_exact_x0(self, capsys, mat, rhs):
        xe = splu(mat).solve(rhs)
        rtole = np.linalg.norm(mat@xe - rhs) / np.linalg.norm(rhs)
        x, info = symmlq(mat, rhs, x0=xe, tol=rtole, verbose=True)
        assert_equal(info, 0)
        assert_array_equal(x, xe)
        out, err = capsys.readouterr()
        assert_equal(out, "SYMMLQ: Linear solve converged due to reach TOL "
                          "iterations 0\n")
        assert_equal(err, '')

    @pytest.mark.parametrize(('mat', 'rhs'), [(A, b), (Ac, bc)])
    def test_b_shape_Nx1(self, mat, rhs):
        def matvec(x):
            invD = 1. / abs(mat.diagonal())
            return invD * x
        M = LinearOperator(mat.shape, matvec=matvec, dtype=mat.dtype)
        bblock = rhs.reshape((rhs.shape[0], 1))
        x, info = symmlq(mat, bblock, M=M)
        assert_equal(info, 0)
        assert_normclose(mat@x, rhs, tol=1e-5)

    @pytest.mark.parametrize(('mat', 'rhs'), [(A, b), (Ac, bc)])
    def test_atol(self, mat, rhs):
        def matvec(x):
            invD = 1. / abs(mat.diagonal())
            return invD * x
        M = LinearOperator(mat.shape, matvec=matvec, dtype=mat.dtype)
        rng = np.random.default_rng(12345)
        randtol = rng.beta(.1, 1.2, size=(2,))
        rtol = randtol[0]
        atol = np.linalg.norm(rhs) * randtol[1]
        x, info = symmlq(mat, rhs, tol=rtol, atol=atol, M=M)
        assert_equal(info, 0)
        assert_normclose(mat@x, rhs, tol=max(randtol))
