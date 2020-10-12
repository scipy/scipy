# coverage test for scipy arls.py
import numpy as np
from scipy.linalg import hilbert
from scipy.linalg.misc import norm
from scipy.linalg.autosol import (arls, arlsusv, cull, prepeq, arlseq, arlsgt,
                      arlsnn, splita, splitb, decide_width)
from scipy.linalg import autosol
from numpy.testing import assert_, assert_array_almost_equal_nulp


def hilb(m, n):
    A = np.zeros((m, n))
    for i in range(0, m):
        for j in range(0, n):
            A[i, j] = 1 / float(i + j + 1)
    return A


# test solvers with zero b
A = np.eye(3)
Z = A
b = np.zeros(3)
x = arlsusv(A, b, Z, Z, Z)
assert_(norm(x) == 0.0, "Solution of arls() is incorrectly non-zero.")
x = arls(A, b)
assert_(norm(x) == 0.0, "Solution of arls() is incorrectly non-zero.")
x = arlseq(A, b, A, b)
assert_(norm(x) == 0.0, "Solution of arlseq() is incorrectly non-zero.")
x = arlsgt(A, b, A, b)
assert_(norm(x) == 0.0, "Solution of arlsgt() is incorrectly non-zero.")
x = arlsnn(A, b)
assert_(norm(x) == 0.0, "Solution of arlsnn() is incorrectly non-zero.")

# TEST ARLS UTILITIES

E = np.eye(4)
f = np.ones(4)
E[2, 2] = 0.0000001
E[3, 3] = 0.0000001
E, f = cull(E, f, 0.00001)
assert_(E.shape[0] == 2, "cull is not deleting row properly.")

# TEST ARLS

# test arls() with mulriple columns
A = np.eye(3)
b = np.ones((3, 2))
b[2, 1] = 0.0
ans = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 0.0]])
x = arls(A, b)
assert_array_almost_equal_nulp(x, ans, 10)

# test arls() with ill-conditioned problem
A = hilbert(12)
xx = np.ones(12)
b = A @ xx
for i in range(0, 12):
    b[i] += 0.00001 * np.sin(float(i + i))
ans =np.array([0.998635, 1.013942, 0.980540, 0.986143,
               1.000395, 1.011578, 1.016739, 1.015970,
               1.010249, 1.000679, 0.988241, 0.973741])
x = arls(A, b)
d = norm(ans - x)
assert_(d < 1.0e-5,
"arls() residual too large in test of ill-conditioned problem.")

# TEST ARLSEQ UTILITIES

# assure splita is working inside arsleq
n = 14
E = np.eye(n)
f = np.array([21., 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 20.])
EE, ff = prepeq(E, f, 0.0000001)
assert_(EE.shape[0] < n, "splita failed inside prepeq.")

# TEST ARLSEQ

# test arlseq() row interchange
A = np.eye(3)
b = np.ones(3)
E = np.array([[1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
f = np.array([2.0, 1.0, 2.0])
ans = np.array([1.0, 1.0, 2.0])
x = arlseq(A, b, E, f)
assert_array_almost_equal_nulp(x, ans, 10)

# test arlseq() row deletion
A = np.eye(3)
b = np.ones(3)
E = A.copy()
E[2,2] = 0.0
f = b.copy()
f[2] = 0.0
ans = np.array([1.0, 1.0, 1.0])
x = arlseq(A, b, E, f)
assert_array_almost_equal_nulp(x, ans, 10)

# TEST ARLSGT

# test of arlsgt
x = np.array([6.0, 5.0, 4.0, 3.0])
A = hilb(5, 4)
b = A @ x
G = np.array([0.0, 0.0, 0.0, 1.0])
h = 5
x = arlsgt(A, b, G, h)
ans = [5.90761758, 6.18916743, 0.99658155, 5.0]
d = norm(ans - x)
assert_(d < 1.0e-6, "Residual too large in arlsgt hilbert test.")

# force solution of arlsgt to be zero
A = np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
b = np.array([0.0, 1.0])
G = np.eye(3)
h = np.zeros(3)
x = arlsgt(A, b, G, h)
d = norm(x)
assert_(norm(x) == 0.0, "Solution of arlsgt() should be zero.")

# TEST ARLSNN

# test too few columns
A = np.ones((3, 1))
b = [1.0, 1.0, 1.0]
x = arlsnn(A, b)
assert_(x[0] == 1.0, "arlsnn not handling single column right.")

# test arlsnn with unnecessary call
A = np.eye(4)
b = np.ones(4)
x = arlsnn(A, b)
assert_(abs(norm(x) - 2.0) < 1.0e-8, "Solution of arls() is wrong.")

# test arlsnn with impossible problem
A = np.eye(3)
b = np.ones(3)
b = -b
x = arlsnn(A, b)
assert_(norm(x) == 0.0, "Solution of arls() is  incorrectly non-zero.")

# TEST UTILITY ROUTINES

# assure decide_width is working for large n
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

# test splita()
g = np.array([1.0, 1.0, 0.0, 20.0])
ans = np.array([1, 2, 3, 3])
r = np.zeros(4)
for i in range(1, 5):
    r[i - 1] = splita(i, g)
assert_array_almost_equal_nulp(r, ans, 10)

# test splitb()
g = np.array([1.0, 0.1, 0.01, 0.1, 1.0, 10.0])
ans = np.array([1, 2, 3, 4, 4, 4])
r = np.zeros(6)
for i in range(1, 7):
    r[i - 1] = splitb(i, g)
assert_array_almost_equal_nulp(r, ans, 10)
