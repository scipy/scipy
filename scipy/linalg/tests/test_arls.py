# coverage test for scipy arls.py
import numpy as np
from scipy.linalg import hilbert
from scipy.linalg.misc import norm
from numpy.testing import assert_, assert_array_almost_equal_nulp

from scipy.linalg.arls import (arls, arlsusv, cull, prepeq, arlseq,
           arlsgt, arlsnn, splita, splitb, decide_width)


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
        r[i - 1] = splita(i, g)
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
        r[i - 1] = splitb(i, g)
        assert_(r[i-1] == ans[i-1], "Splitb is miscalibrated.")  
    return

        
def test_cull():
    E = np.eye(4)
    f = np.ones(4)
    E[2, 2] = 0.0000001
    E[3, 3] = 0.0000001
    E, f = cull(E, f, 0.00001)
    assert_(E.shape[0] == 2, "cull is not deleting row properly.")
    return


def test_zero_rhs():
    # test solvers with zero right hand side
    A = np.eye(3)
    Z = A
    b = np.zeros(3)
    x = arlsusv(A, b, Z, Z, Z)[0]
    assert_(norm(x) == 0.0, "Solution of arls() is incorrectly non-zero.")
    x = arls(A, b)[0]
    assert_(norm(x) == 0.0, "Solution of arls() is incorrectly non-zero.")
    x = arlseq(A, b, A, b)[0]
    assert_(norm(x) == 0.0, "Solution of arlseq() is incorrectly non-zero.")
    x = arlsgt(A, b, A, b)[0]
    assert_(norm(x) == 0.0, "Solution of arlsgt() is incorrectly non-zero.")
    x = arlsnn(A, b)[0]
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
    return


def test_regularization():
    A = hilbert(12)
    xx = np.ones(12)
    b = A @ xx
    for i in range(0, 12):
        b[i] += 0.00001 * np.sin(float(i + i))
    ans =np.array([0.998635, 1.013942, 0.980540, 0.986143,
                   1.000395, 1.011578, 1.016739, 1.015970,
                   1.010249, 1.000679, 0.988241, 0.973741])
    x = arls(A, b)[0]
    d = norm(ans - x)
    assert_(d < 1.0e-5, "arls() residual too large.")
    return


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
    E[2,2] = 0.0
    f = b.copy()
    f[2] = 0.0
    ans = np.array([1.0, 1.0, 1.0])
    x = arlseq(A, b, E, f)[0]
    assert_array_almost_equal_nulp(x, ans, 10)
    return


def test_arlseq_with_poor_E():
    n=14
    x = np.ones(n)
    x[1] = -1.0
    x[5] = -1.0
    x[9] = -1.0
    A = np.eye(n)
    b = A @ x
    E = hilbert(n)
    f = E @ x
    xx=arlseq(A,b,E,f)[0]
    ans =np.array([ 1.00090446, -0.99817657,  1.00236951,  1.00279655,
                    1.00221899, -0.99679898,  1.00738824,  0.9862784,
                    1.02311775, -0.99162213,  0.95629113,  1.05844299,
                     0.97236068,  1.0087402 ])
    d = norm(xx-ans)
    assert_(d < 1.0e-6,"Residual too large in arlseq.")
    return


# TEST ARLSGT
def test_arlsgt():
    x = np.array([6.0, 5.0, 4.0, 3.0])
    A = hilbert(5)
    A = np.delete(A,4,1) #delete last column
    b = A @ x
    G = np.array([0.0, 0.0, 0.0, 1.0])
    h = 5
    x = arlsgt(A, b, G, h)[0]
    ans = [5.90761758, 6.18916743, 0.99658155, 5.0]
    d = norm(ans - x)
    assert_(d < 1.0e-5,"Residual too large in arlsgt hilbert test.")
    return    

def test_forced_zero_solution():
    A = np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
    b = np.array([0.0, 1.0])
    G = np.eye(3)
    h = np.zeros(3)
    x = arlsgt(A, b, G, h)[0]
    d = norm(x)
    assert_(norm(x) == 0.0, "Solution by arlsgt should be zero.")
    return


# TEST ARLSNN
def test_arlsnn_single_column():
    A = np.ones((3, 1))
    b = [1.0, 1.0, 1.0]
    x = arlsnn(A, b)[0]
    assert_(x[0] == 1.0, "arlsnn not handling single column right.")
    return

def test_arlsnn_with_impossible():
    A = np.eye(3)
    b = np.ones(3)
    b = -b
    x = arlsnn(A, b)[0]
    assert_(norm(x) == 0.0, "Solution of arlsnn is incorrectly non-zero.")
    return


