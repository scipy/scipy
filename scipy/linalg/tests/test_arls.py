# coverage test for scipy arls.py
import numpy as np
from scipy.linalg import hilbert
from scipy.linalg.misc import norm
from numpy.testing import assert_, assert_array_almost_equal_nulp

from scipy.linalg.arls import (arls, strange, arlsqr, arlsusv, prepeq, arlseq,
    arlsgt, arlsnn, splita, splitb, decide_width)


def myrandom(m,n):
    A = np.zeros((m,n))
    for i in range(0, m):
        for j in range(0, n):
            A[i,j] = abs(np.sin(float(2*m+3*n) + 2.0*float(i) + 2.5*float(j)))
    return A
    

# TEST LOW LEVEL UTILTIES
def test_decide_width():
    k = decide_width(2)
    assert_(k == 1, "Algorithm in decide_width has changed.")
    k = decide_width(5)
    assert_(k == 2, "Algorithm in decide_width has changed.")
    k = decide_width(12)
    assert_(k == 4, "Algorithm in decide_width has changed.")
    k = decide_width(30)
    assert_(k == 6, "Algorithm in decide_width has changed.")
    k = decide_width(50)
    assert_(k == 8, "Algorithm in decide_width has changed.")
    k = decide_width(80)
    assert_(k == 10, "Algorithm in decide_width has changed.")
    k = decide_width(150)
    assert_(k == 14, "Algorithm in decide_width has changed.")
    return


def test_splita():
    g = np.array([1.0, 1.0, 0.0, 20.0])
    ans = np.array([1, 2, 3, 3])
    r = np.zeros(4)
    for i in range(1, 5):
        r[i - 1] = splita(g, i)
        assert_(r[i-1] == ans[i-1], "Splita is miscalibrated.")
    return


def test_splita_in_arlseq():
    n = 14
    E = np.eye(n)
    f = np.array([21., 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
                 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 20.])
    EE, ff = prepeq(E, f, 0.0000001)
    assert_(EE.shape[0] < n, "splita failed inside prepeq.")
    return


def test_splitb():
    g = np.array([1.0, 0.1, 0.01, 0.1, 1.0, 10.0])
    ans = np.array([1, 2, 3, 4, 4, 4])
    r = np.zeros(6)
    for i in range(1, 7):
        r[i - 1] = splitb(g, i)
        assert_(r[i-1] == ans[i-1], "Splitb is miscalibrated.")
    return


def test_zero_rhs():
    # test solvers with zero right hand side
    A = np.eye(3)
    Z = A
    b = np.zeros(3)
    x, nr, ur, sigma, lambdah = arls(A, b)
    assert_(norm(x) == 0.0, "Solution of arls() is incorrectly non-zero.")
    x, nr, ur, sigma, lambdah = arlsqr(A, b)
    assert_(norm(x) == 0.0, "Solution of arlsqr() is incorrectly non-zero.")
    x, nr, ur, sigma, lambdah = arlsusv(A, b, Z, Z, Z)
    assert_(norm(x) == 0.0, "Solution of arlsusv() is incorrectly non-zero.")
    x, nr, ur, sigma, lambdah = arlseq(A, b, A, b)
    assert_(norm(x) == 0.0, "Solution of arlseq() is incorrectly non-zero.")
    x, nr, ur, sigma, lambdah = arlsgt(A, b, A, b)
    assert_(norm(x) == 0.0, "Solution of arlsgt() is incorrectly non-zero.")
    x, nr, ur, sigma, lambdah = arlsnn(A, b)
    assert_(norm(x) == 0.0, "Solution of arlsnn() is incorrectly non-zero.")
    return


# TEST ARLS
def test_multiple_columns():
    A = np.eye(3)
    b = np.ones((3, 2))
    b[2, 1] = 0.0
    ans = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 0.0]])
    x = arls(A, b)[0]
    assert_array_almost_equal_nulp(x, ans, 10)
    x = arlsqr(A, b)[0]
    assert_array_almost_equal_nulp(x, ans, 10)
    return


def test_regularization():
    A = hilbert(12)
    xx = np.ones(12)
    b = A @ xx
    for i in range(0, 12):
        b[i] += 0.00001 * np.sin(float(i + i))
    ans = np.array([0.998635, 1.013942, 0.980540, 0.986143,
                    1.000395, 1.011578, 1.016739, 1.015970,
                    1.010249, 1.000679, 0.988241, 0.973741])
    x = arls(A, b)[0]
    d = norm(ans - x)
    assert_(d < 1.0e-5, "arls() residual too large.")

    x = arlsqr(A, b)[0]
    d = norm(ans - x)
    assert_(d < 1.0e-5, "arlsqr() residual too large.")
    return


# TEST ARLSQR
def test_arlsqr_overdetermined():
    A = myrandom(6, 4)
    x = np.zeros(4)
    for i in range(0, 4):
        x[i] = float(i)
    b = A @ x
    xx = arlsqr(A, b)[0]
    res = norm(A @ xx - b)
    assert_(res < 1.0e-8, "arls() residual too large.")


def test_arlsqr_underdetermined():
    A = myrandom(4, 6)
    x = np.zeros(6)
    for i in range(0, 6):
        x[i] = float(i)
    b = A @ x
    xx = arlsqr(A, b)[0]
    res = norm(A @ xx - b)
    assert_(res < 1.0e-8, "arlsqr() residual too large.")


def test_strange():
    A = myrandom(4, 6)
    A[2,:] = 0.0
    x = np.zeros(6)
    for i in range(0, 6):
        x[i] = float(i)
    b = A @ x
    odd, maxnorm = strange(A,b)
    assert_(odd, "strange() should report TRUE")


# TEST ARLSEQ
def test_row_interchange():
    A = np.eye(3)
    b = np.ones(3)
    E = np.array([[1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    f = np.array([2.0, 1.0, 2.0])
    ans = np.array([1.0, 1.0, 2.0])
    x = arlseq(A, b, E, f)[0]
    assert_array_almost_equal_nulp(x, ans, 10)
    return


def test_row_deletion():
    A = np.eye(3)
    b = np.ones(3)
    E = A.copy()
    E[2, 2] = 0.0
    f = b.copy()
    f[2] = 0.0
    ans = np.array([1.0, 1.0, 1.0])
    x = arlseq(A, b, E, f)[0]
    assert_array_almost_equal_nulp(x, ans, 10)
    return


def test_arlseq():
    n = 14
    A = np.eye(n)
    x = np.ones(n)
    b = A @ x

    E = hilbert(n)
    E[7, 7] = 3.
    x[7] = 5.0
    f = E @ x

    xx = arlseq(A, b, E, f)[0]
    assert_(abs(xx[7] - 5.0) < 1.0e-4, "Constraint not obeyed in arlseq.")

    d = norm(x - xx)
    assert_(d < 0.01, "Residual too large in arsleq.")
    return


# TEST ARLSGT
def test_arlsgt():
    x = np.array([6.0, 5.0, 4.0, 3.0])
    A = hilbert(5)
    A = np.delete(A, 4, 1)   # delete last column
    b = A @ x
    G = np.array([0.0, 0.0, 0.0, 1.0])
    h = 5
    x = arlsgt(A, b, G, h)[0]
    res = A @ x - b
    assert_(norm(res) < 0.002, "Residual too large in arlsgt hilbert test.")
    return


def test_forced_zero_solution():
    A = np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
    b = np.array([0.0, 1.0])
    G = np.eye(3)
    h = np.zeros(3)
    x = arlsgt(A, b, G, h)[0]
    assert_(norm(x) == 0.0, "Solution by arlsgt should be zero.")
    return


# TEST ARLSNN
def test_arlsnn_single_column():
    A = np.ones((3, 1))
    b = [1.0, 1.0, 1.0]
    x, nr, ur, sigma, lambdah = arlsnn(A, b)
    assert_(x[0] == 1.0, "arlsnn not handling single column right.")
    return


def test_arlsnn_with_impossible():
    A = np.eye(3)
    b = np.ones(3)
    b = -b
    x = arlsnn(A, b)[0]
    assert_(norm(x) == 0.0, "Solution of arlsnn is incorrectly non-zero.")
    return
