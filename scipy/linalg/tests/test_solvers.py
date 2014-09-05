from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import inv

from numpy.testing import TestCase, rand, run_module_suite, assert_raises, \
    assert_equal, assert_almost_equal, assert_array_almost_equal, assert_, \
    assert_allclose

from scipy.linalg import solve_sylvester, solve_lyapunov, \
    solve_discrete_lyapunov, solve_continuous_are, solve_discrete_are


class TestSolveLyapunov(TestCase):

    cases = [
        (np.array([[1, 2], [3, 4]]),
         np.array([[9, 10], [11, 12]])),
        # a, q all complex.
        (np.array([[1.0+1j, 2.0], [3.0-4.0j, 5.0]]),
         np.array([[2.0-2j, 2.0+2j],[-1.0-1j, 2.0]])),
        # a real; q complex.
        (np.array([[1.0, 2.0], [3.0, 5.0]]),
         np.array([[2.0-2j, 2.0+2j],[-1.0-1j, 2.0]])),
        # a complex; q real.
        (np.array([[1.0+1j, 2.0], [3.0-4.0j, 5.0]]),
         np.array([[2.0, 2.0],[-1.0, 2.0]])),
        # An example from Kitagawa, 1977
        (np.array([[3, 9, 5, 1, 4], [1, 2, 3, 8, 4], [4, 6, 6, 6, 3],
                   [1, 5, 2, 0, 7], [5, 3, 3, 1, 5]]),
         np.array([[2, 4, 1, 0, 1], [4, 1, 0, 2, 0], [1, 0, 3, 0, 3],
                   [0, 2, 0, 1, 0], [1, 0, 3, 0, 4]])),
        # Companion matrix example. a complex; q real; a.shape[0] = 11
        (np.array([[0.100+0.j, 0.091+0.j, 0.082+0.j, 0.073+0.j, 0.064+0.j,
                    0.055+0.j, 0.046+0.j, 0.037+0.j, 0.028+0.j, 0.019+0.j,
                    0.010+0.j],
                   [1.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 1.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 1.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 1.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 1.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    1.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 1.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 1.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 1.000+0.j, 0.000+0.j,
                    0.000+0.j],
                   [0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j,
                    0.000+0.j, 0.000+0.j, 0.000+0.j, 0.000+0.j, 1.000+0.j,
                    0.000+0.j]]),
         np.eye(11)),
        ]

    def check_continuous_case(self, a, q):
        x = solve_lyapunov(a, q)
        assert_array_almost_equal(np.dot(a, x) + np.dot(x, a.conj().transpose()), q)

    def check_discrete_case(self, a, q, method=None):
        x = solve_discrete_lyapunov(a, q, method=method)
        assert_array_almost_equal(np.dot(np.dot(a, x),a.conj().transpose()) - x, -1.0*q)

    def test_cases(self):
        for case in self.cases:
            self.check_continuous_case(case[0], case[1])
            self.check_discrete_case(case[0], case[1])
            self.check_discrete_case(case[0], case[1], method='direct')
            self.check_discrete_case(case[0], case[1], method='bilinear')


class TestSolveContinuousARE(TestCase):

    cases = [
        # An example from Laub, A. J.
        # (http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf)
        (np.matrix([[0, 1], [0, 0]]),
         np.matrix([[0,], [1,]]),
         np.matrix([[1, 0], [0, 2]]),
         np.matrix([[1,],])),
        # Difficult from a numerical standpoint, again from Laub, A. J.
        (np.matrix([[4, 3], [-9.0/2.0, -7.0/2.0]]),
         np.matrix([[1,], [-1,]]),
         np.matrix([[9, 6], [6, 4]]),
         np.matrix([[1,],])),
        # Complex a; real b, q, r
        (np.matrix([[0, 1-2j], [0, -3j]]),
         np.matrix([[0,], [1,]]),
         np.matrix([[1, 0], [0, 2]]),
         np.matrix([[1,],])),
        # Real a, q, r; complex b
        (np.matrix([[0, 1], [0, -1]]),
         np.matrix([[-2j,], [1j,]]),
         np.matrix([[1, 0], [0, 2]]),
         np.matrix([[1,],])),
        # Real a, b; complex q, r
        (np.matrix([[0, 1], [0, -1]]),
         np.matrix([[1, 2], [1, 3]]),
         np.matrix([[1, -3j], [1-1j, 2]]),
         np.matrix([[-2j, 2], [1j, 3]])),
        ]

    def check_case(self, a, b, q, r):
        """Checks if (A'X + XA - XBR^-1B'X+Q=0) is true"""

        x = solve_continuous_are(a, b, q, r)
        assert_array_almost_equal(
            a.getH()*x + x*a - x*b*inv(r)*b.getH()*x + q, 0.0)

    def test_cases(self):
        for case in self.cases:
            self.check_case(case[0], case[1], case[2], case[3])


class TestSolveDiscreteARE(TestCase):

    cases = [
        # Difficult from a numerical standpoint, again from Laub, A. J.
        # (http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf)
        (np.matrix([[4, 3], [-9.0/2.0, -7.0/2.0]]),
         np.matrix([[1,], [-1,]]),
         np.matrix([[9, 6], [6, 4]]),
         np.matrix([[1,],])),
        # Another example from Laub
        (np.matrix([[0.9512, 0], [0, 0.9048]]),
         np.matrix([[4.877, 4.877], [-1.1895, 3.569]]),
         np.matrix([[0.005, 0],[0, 0.02]]),
         np.matrix([[1.0/3.0, 0],[0, 3]])),
        # Complex a; real b, q, r
        (np.matrix([[2, 1-2j], [0, -3j]]),
         np.matrix([[0,], [1,]]),
         np.matrix([[1, 0], [0, 2]]),
         np.matrix([[1,],])),
        # Real a, q, r; complex b
        (np.matrix([[2, 1], [0, -1]]),
         np.matrix([[-2j,], [1j,]]),
         np.matrix([[1, 0], [0, 2]]),
         np.matrix([[1,],])),
        # Real a, b; complex q, r
        (np.matrix([[3, 1], [0, -1]]),
         np.matrix([[1, 2], [1, 3]]),
         np.matrix([[1, -3j], [1-1j, 2]]),
         np.matrix([[-2j, 2], [1j, 3]])),
        ]

    def check_case(self, a, b, q, r):
        """Checks if X = A'XA-(A'XB)(R+B'XB)^-1(B'XA)+Q) is true"""

        x = solve_discrete_are(a, b, q, r)
        assert_array_almost_equal(
            a.getH()*x*a-(a.getH()*x*b)*inv(r+b.getH()*x*b)*(b.getH()*x*a)+q-x, 0.0)

    def test_cases(self):
        for case in self.cases:
            self.check_case(case[0], case[1], case[2], case[3])


class TestSolveSylvester(TestCase):

    cases = [
        # a, b, c all real.
        (np.array([[1, 2], [0, 4]]),
         np.array([[5, 6], [0, 8]]),
         np.array([[9, 10], [11, 12]])),
        # a, b, c all real, 4x4. a and b have non-trival 2x2 blocks in their
        # quasi-triangular form.
        (np.array([[1.0, 0, 0, 0], [0, 1.0, 2.0, 0.0], [0, 0, 3.0, -4], [0, 0, 2, 5]]),
         np.array([[2.0, 0, 0,1.0], [0, 1.0, 0.0, 0.0], [0, 0, 1.0, -1], [0, 0, 1, 1]]),
         np.array([[1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]])),
        # a, b, c all complex.
        (np.array([[1.0+1j, 2.0], [3.0-4.0j, 5.0]]),
         np.array([[-1.0, 2j], [3.0, 4.0]]),
         np.array([[2.0-2j, 2.0+2j],[-1.0-1j, 2.0]])),
        # a and b real; c complex.
        (np.array([[1.0, 2.0], [3.0, 5.0]]),
         np.array([[-1.0, 0], [3.0, 4.0]]),
         np.array([[2.0-2j, 2.0+2j],[-1.0-1j, 2.0]])),
        # a and c complex; b real.
        (np.array([[1.0+1j, 2.0], [3.0-4.0j, 5.0]]),
         np.array([[-1.0, 0], [3.0, 4.0]]),
         np.array([[2.0-2j, 2.0+2j],[-1.0-1j, 2.0]])),
        # a complex; b and c real.
        (np.array([[1.0+1j, 2.0], [3.0-4.0j, 5.0]]),
         np.array([[-1.0, 0], [3.0, 4.0]]),
         np.array([[2.0, 2.0],[-1.0, 2.0]])),
        # not square matrices, real
        (np.array([[8, 1, 6], [3, 5, 7], [4, 9, 2]]),
         np.array([[2, 3], [4, 5]]),
         np.array([[1, 2], [3, 4], [5, 6]])),
        # not square matrices, complex
        (np.array([[8, 1j, 6+2j], [3, 5, 7], [4, 9, 2]]),
         np.array([[2, 3], [4, 5-1j]]),
         np.array([[1, 2j], [3, 4j], [5j, 6+7j]])),
    ]

    def check_case(self, a, b, c):
        x = solve_sylvester(a, b, c)
        assert_array_almost_equal(np.dot(a, x) + np.dot(x, b), c)

    def test_cases(self):
        for case in self.cases:
            self.check_case(case[0], case[1], case[2])

    def test_trivial(self):
        a = np.array([[1.0, 0.0], [0.0, 1.0]])
        b = np.array([[1.0]])
        c = np.array([2.0, 2.0]).reshape(-1,1)
        x = solve_sylvester(a, b, c)
        assert_array_almost_equal(x, np.array([1.0, 1.0]).reshape(-1,1))

if __name__ == "__main__":
    run_module_suite()
