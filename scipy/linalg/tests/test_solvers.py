from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import inv, solve

from numpy.testing import TestCase, rand, run_module_suite, assert_raises, \
    assert_equal, assert_almost_equal, assert_array_almost_equal, assert_, \
    assert_allclose

from numpy.testing.noseclasses import KnownFailureTest
    
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
        # https://github.com/scipy/scipy/issues/4176
        (np.matrix([[0, 1], [-1/2, -1]]),
         (np.matrix([0, 3]).T * np.matrix([0, 3]).T.T)),
        # https://github.com/scipy/scipy/issues/4176
        (np.matrix([[0, 1], [-1/2, -1]]),
         (np.array(np.matrix([0, 3]).T * np.matrix([0, 3]).T.T))),
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

def test_solve_discrete_are():

    cases = [
        # Darex examples taken from (with default parameters):
        # [1] P.BENNER, A.J. LAUB, V. MEHRMANN: 'A Collection of Benchmark 
        #     Examples for the Numerical Solution of Algebraic Riccati 
        #     Equations II: Discrete-Time Case', Tech. Report SPC 95_23, 
        #     Fak. f. Mathematik, TU Chemnitz-Zwickau (Germany), 1995.
        # [2] T. GUDMUNDSSON, C. KENNEY, A.J. LAUB: 'Scaling of the 
        #     Discrete-Time Algebraic Riccati Equation to Enhance Stability 
        #     of the Schur Solution Method', IEEE Trans.Aut.Cont., vol.37(4) 
        #
        # The format of the data is (a, b, q, r, knownfailure), where
        # knownfailure is None if the test passes or a string
        # indicating the reason for failure.
        #
        # Complex a; real b, q, r
        (np.array([[2, 1-2j], [0, -3j]]),
         np.array([[0,], [1,]]),
         np.array([[1, 0], [0, 2]]),
         np.array([[1,],]),
         None),
        # Real a, q, r; complex b
        (np.array([[2, 1], [0, -1]]),
         np.array([[-2j,], [1j,]]),
         np.array([[1, 0], [0, 2]]),
         np.array([[1,],]),
         None),
        # Real a, b; complex q, r 
        (np.array([[3, 1], [0, -1]]),
         np.array([[1, 2], [1, 3]]),
         np.array([[1, 1+1j], [1-1j, 2]]),
         np.array([[2, -2j], [2j, 3]]),
         None),
        # User-reported gh-2251 (Trac #1732)
        (np.array([[0.63399379, 0.54906824, 0.76253406],
                   [0.5404729, 0.53745766, 0.08731853],
                   [0.27524045, 0.84922129, 0.4681622]]),
         np.array([[0.96861695],[0.05532739],[0.78934047]]),
         np.eye(3),
         np.eye(1),
         None),
        # darex #1
        (np.array([[4, 3],[-4.5, -3.5]]),
         np.array([[1],[-1]]),
         np.array([[9, 6],[6, 4]]),
         np.array([[1]]),
         None),
        # darex #2
        (np.array([[0.9512, 0],[0, 0.9048]]),
         np.array([[4.877, 4.877],[-1.1895, 3.569]]),
         np.array([[0.005, 0],[0, 0.02]]),
         np.array([[1/3, 0],[0, 3]]),
         None),
        # darex #3
        (np.array([[2, -1],[1, 0]]),
         np.array([[1],[0]]),
         np.array([[0, 0],[0, 1]]),
         np.array([[0]]),
         None),
        # darex #4 (skipped the gen. Ric. term S)
        (np.array([[0, 1],[0, -1]]),
         np.array([[1, 0],[2, 1]]),
         np.array([[-4, -4],[-4, 7]]) * (1/11),
         np.array([[9, 3],[3, 1]]),
         None),
        # darex #5
        (np.array([[0, 1],[0, 0]]),
         np.array([[0],[1]]),
         np.array([[1, 2],[2, 4]]),
         np.array([[1]]),
         None),
        # darex #6
        (np.array([[0.998, 0.067, 0, 0],
                   [-.067, 0.998, 0, 0],
                   [0, 0, 0.998, 0.153],
                   [0, 0, -.153, 0.998]]),
         np.array([[0.0033, 0.0200],
                   [0.1000, -.0007],
                   [0.0400, 0.0073],
                   [-.0028, 0.1000]]),
         np.array([[1.87, 0, 0, -0.244],
                   [0, 0.744, 0.205, 0],
                   [0, 0.205, 0.589, 0],
                   [-0.244, 0, 0, 1.048]]),
         np.eye(2),
         None),
        # darex #7
        (np.array([[0.984750, -.079903, 0.0009054, -.0010765],
                   [0.041588, 0.998990, -.0358550, 0.0126840],
                   [-.546620, 0.044916, -.3299100, 0.1931800],
                   [2.662400, -.100450, -.9245500, -.2632500]]),
         np.array([[0.0037112, 0.0007361],
                   [-.0870510, 9.3411e-6],
                   [-1.198440, -4.1378e-4],
                   [-3.192700, 9.2535e-4]]),
         np.eye(4)*1e-2,
         np.eye(2),
         None),
        # darex #8
        (np.array([[-0.6000000, -2.2000000, -3.6000000, -5.4000180],
                   [1.0000000, 0.6000000, 0.8000000, 3.3999820],
                   [0.0000000, 1.0000000, 1.8000000, 3.7999820],
                   [0.0000000, 0.0000000, 0.0000000, -0.9999820]]),
         np.array([[1.0, -1.0, -1.0, -1.0],
                   [0.0, 1.0, -1.0, -1.0],
                   [0.0, 0.0, 1.0, -1.0],
                   [0.0, 0.0, 0.0, 1.0]]),
         np.array([[2, 1, 3, 6],
                   [1, 2, 2, 5],
                   [3, 2, 6, 11],
                   [6, 5, 11, 22]]),
         np.eye(4),
         None),
        # darex #9
        (np.array([[95.4070, 1.9643, 0.3597, 0.0673, 0.0190],
                    [40.8490, 41.3170, 16.0840, 4.4679, 1.1971],
                   [12.2170, 26.3260, 36.1490, 15.9300, 12.3830],
                   [4.1118, 12.8580, 27.2090, 21.4420, 40.9760],
                   [0.1305, 0.5808, 1.8750, 3.6162, 94.2800]]) * 0.01,
         np.array([[0.0434, -0.0122],
                   [2.6606, -1.0453],
                   [3.7530, -5.5100],
                   [3.6076, -6.6000],
                   [0.4617, -0.9148]]) * 0.01,
         np.eye(5),
         np.eye(2),
         None),
        # darex #10
        (np.kron(np.eye(2),np.diag([1,1],k=1)),
         np.kron(np.eye(2),np.array([[0],[0],[1]])),
         np.array([[1, 1, 0, 0, 0, 0],
               [1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, -1, 0],
               [0, 0, 0, -1, 1, 0],
               [0, 0, 0, 0, 0, 0]]),
         np.array([[3, 0],[0, 1]]),
         None),
        # darex #11
        (0.001 * np.array(
         [[870.1, 135.0, 11.59, .5014, -37.22, .3484, 0, 4.242, 7.249],
         [76.55, 897.4, 12.72, 0.5504, -40.16, .3743, 0, 4.53, 7.499],
         [-127.2, 357.5, 817, 1.455, -102.8, .987, 0, 11.85, 18.72],
         [-363.5, 633.9, 74.91, 796.6, -273.5, 2.653, 0, 31.72, 48.82],
         [-960, 1645.9, -128.9, -5.597, 71.42, 7.108, 0, 84.52, 125.9],
         [-664.4, 112.96, -88.89, -3.854, 84.47, 13.6, 0, 144.3, 101.6],
         [-410.2, 693, -54.71, -2.371, 66.49, 12.49, .1063, 99.97, 69.67],
         [-179.9, 301.7, -23.93, -1.035, 60.59, 22.16, 0, 213.9, 35.54],
         [-345.1, 580.4, -45.96, -1.989, 105.6, 19.86, 0, 219.1, 215.2]]),
         np.array([[4.7600, -0.5701, -83.6800],
               [0.8790, -4.7730, -2.7300],
               [1.4820, -13.1200, 8.8760],
               [3.8920, -35.1300, 24.8000],
               [10.3400, -92.7500, 66.8000],
               [7.2030, -61.5900, 38.3400],
               [4.4540, -36.8300, 20.2900],
               [1.9710, -15.5400, 6.9370],
               [3.7730, -30.2800, 14.6900]]) * 0.001,
         np.diag([50, 0, 0, 0, 50, 0, 0, 0, 0]),
         np.eye(3),
         None),
        # darex #12
        (np.array([[0, 1e6],[0, 0]]),
         np.array([[0],[1]]),
         np.eye(2),
         np.array([[1]]),
         "Bad absolute accuracy"),
        # darex #13
        (np.array([[16, 10, -2],
                  [10, 13, -8],
                  [-2, -8, 7]]) * (1/9),
         np.eye(3),
         1e6 * np.eye(3),
         1e6 * np.eye(3),
         "Fails to find a valid solution"),
        # darex #14
        (np.array([[1 - 1/1e8, 0, 0, 0],
                  [1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0]]),
         np.array([[1e-08],[0],[0],[0]]),
         np.diag([0, 0, 0, 1]),
         np.array([[0.25]]),
         "Bad absolute accuracy"),
        # darex #15
        (np.eye(100, k=1),
         np.flipud(np.eye(100, 1)),
         np.eye(100),
         np.array([[1]]),
         None)
        ]

    def _test_factory(case):
        """Checks if X = A'XA-(A'XB)(R+B'XB)^-1(B'XA)+Q) is true"""
        a, b, q, r, knownfailure = case
        if knownfailure:
            raise KnownFailureTest(knownfailure)

        x = solve_discrete_are(a, b, q, r)
        res = a.conj().T.dot(x.dot(a)) - x + q
        res -= a.conj().T.dot(x.dot(b)).dot(
                    solve(r+b.conj().T.dot(x.dot(b)),b.conj().T).dot(x.dot(a))
                    )
        assert_array_almost_equal(res,np.zeros_like(res))

    for case in cases:
        yield _test_factory, case


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
