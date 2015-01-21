from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy.testing import (assert_almost_equal, assert_equal, assert_allclose,
                           assert_, assert_raises, TestCase, run_module_suite)
from scipy.signal.ltisys import (ss2tf, tf2ss, lsim2, impulse2, step2, lti,
                                 bode, freqresp, impulse, step,
                                 abcd_normalize, place_poles)
from scipy.signal.filter_design import BadCoefficients
import scipy.linalg as linalg


def _assert_poles_close(P1,P2, rtol=1e-8, atol=1e-8):
    """
    Check each pole in P1 is close to a pole in P2 with a 1e-8
    relative tolerance or 1e-8 absolute tolerance (useful for zero poles). 
    These tolerance are very scrict but the systems tested are known to 
    accept these poles so we should not be far from what is requested.    
    """
    P2 = P2.copy()
    for p1 in P1:
        found = False
        for p2_idx in range(P2.shape[0]):
            if np.allclose([np.real(p1), np.imag(p1)],
                           [np.real(P2[p2_idx]), np.imag(P2[p2_idx])],
                           rtol, atol):
                found = True
                np.delete(P2, p2_idx)
                break
        if not found:
            raise ValueError("Can't find pole " + str(p1) + " in " + str(P2))


class TestPlacePoles(TestCase):

    def _check(self, A, B, P, method=None):
        """
        Perform the most common tests on the poles computed by place_poles
        and return the Bunch object for further specific tests
        """  
        if method is not None:
            fsf = place_poles(A, B, P, method=method)
        else:
            fsf = place_poles(A, B, P)

        expected, _ = np.linalg.eig(A - np.dot(B, fsf.gain_matrix))
        _assert_poles_close(expected,fsf.requested_poles)
        _assert_poles_close(expected,fsf.computed_poles)
        _assert_poles_close(P,fsf.requested_poles)
            
        return fsf

    def test_real(self):
        # Test real pole placement using KNV and YT0 algorithm and example 1 in
        # section 4 of the reference publication (see place_poles docstring)
        A = np.array([1.380, -0.2077, 6.715, -5.676, -0.5814, -4.290, 0,
                      0.6750, 1.067, 4.273, -6.654, 5.893, 0.0480, 4.273,
                      1.343, -2.104]).reshape(4, 4)
        B = np.array([0, 5.679, 1.136, 1.136, 0, 0, -3.146,0]).reshape(4, 2)
        P = np.array([-0.2, -0.5, -5.0566, -8.6659])

        # Check that both KNV and YT compute correct K matrix
        self._check(A, B, P, method='KNV0')
        self._check(A, B, P, method='YT')

    def test_complex1(self):
        # Test complex pole placement on a linearized car model, taken from L.
        # Jaulin, Automatique pour la robotique, Cours et Exercices, iSTE
        # editions p 184/185
        A = np.array([0,7,0,0,0,0,0,7/3.,0,0,0,0,0,0,0,0]).reshape(4,4)
        B = np.array([0,0,0,0,1,0,0,1]).reshape(4,2)
        # Test complex poles on YT
        P = np.array([-3, -1, -2-1j, -2+1j])
        self._check(A, B, P)

    def test_complex2(self):
        # need a 5x5 array to ensure YT handles properly when there
        # is only one real pole and several complex
        A = np.array([0,7,0,0,0,0,0,7/3.,0,0,0,0,0,0,0,0,
                      0,0,0,5,0,0,0,0,9]).reshape(5,5)
        B = np.array([0,0,0,0,1,0,0,1,2,3]).reshape(5,2)
        P = np.array([-2, -3+1j, -3-1j, -1+1j, -1-1j])
        place_poles(A,B,P)

        # same test with an odd number of real poles > 1
        # this is another specific case of YT
        P = np.array([-2, -3, -4, -1+1j, -1-1j])
        self._check(A, B, P)

    def test_tricky_B(self):
        # check we handle as we should the 1 column B matrices and
        # n column B matrices (with n such as shape(A)=(n, n))
        A = np.array([1.380, -0.2077, 6.715, -5.676, -0.5814, -4.290, 0,
                      0.6750, 1.067, 4.273, -6.654, 5.893, 0.0480, 4.273,
                      1.343, -2.104]).reshape(4, 4)
        B = np.array([0, 5.679, 1.136, 1.136, 0, 0, -3.146, 0, 1, 2, 3, 4,
                      5, 6, 7, 8]).reshape(4, 4)

        # KNV or YT are not called here, it's a specific case with only
        # one unique solution
        P = np.array([-0.2, -0.5, -5.0566, -8.6659])
        fsf = self._check(A, B, P)
        #rtol and nb_iter should be set to 0 as the solution is unique
        assert_equal(fsf.rtol,0)        
        assert_equal(fsf.nb_iter,0)        

        # check with complex poles too as they trigger a specific case in
        # the specific case :-)
        P = np.array((-2+1j,-2-1j,-3,-2))
        fsf = self._check(A, B, P)
        assert_equal(fsf.rtol,0)        
        assert_equal(fsf.nb_iter,0)       
        
        #now test with a B matrix with only one column (no optimisation)
        B = B[:,0].reshape(4,1)
        P = np.array((-2+1j,-2-1j,-3,-2))
        fsf = self._check(A, B, P)
        
        #rtol and nb_iter are meaningless here as we can't optimize anything,
        #check they are set to NaN as expected         
        assert_equal(fsf.rtol,np.nan)        
        assert_equal(fsf.nb_iter,np.nan)   
                    
    def test_errors(self):
        # Test input mistakes from user
        A = np.array([0,7,0,0,0,0,0,7/3.,0,0,0,0,0,0,0,0]).reshape(4,4)
        B = np.array([0,0,0,0,1,0,0,1]).reshape(4,2)
        
        #should fail as the method keyword is invalid
        assert_raises(ValueError, place_poles, A, B, (-2,-2,-2,-2),
                      method="foo")

        #should fail as the rtol is greater than 1
        assert_raises(ValueError, place_poles, A, B, (-2,-2,-2,-2),
                      rtol=42)

        #should fail as maxiter is smaller than 1
        assert_raises(ValueError, place_poles, A, B, (-2,-2,-2,-2),
                      maxiter=-42)

        # should fail as rank(B) is two
        assert_raises(ValueError, place_poles, A, B, (-2,-2,-2,-2))

        # Should not raise ValueError as the poles can be placed but should
        # raise a warning as the convergence is not reached
        with warnings.catch_warnings(record=True) as w:
            fsf = place_poles(A, B, (-1,-2,-3,-4), rtol=1e-16, maxiter=42)
            assert_(len(w) == 1)
            assert_(issubclass(w[-1].category, UserWarning))
            assert_("Convergence was not reached after maxiter iterations"
                    in str(w[-1].message))
            assert_equal(fsf.nb_iter,42)   

        # should fail as a complex misses its conjugate
        assert_raises(ValueError, place_poles, A, B, (-2+1j,-2-1j,-2+3j,-2))

        # should fail as poles are not in a 1D array-like
        assert_raises(ValueError, place_poles, A, B, ((-2+1j,-2-1j,-2+3j,-2)))

        # should fail as A is not square
        assert_raises(ValueError, place_poles, A[:,:3], B, (-2,-3,-4,-5))

        # should fail as B has not the same number of lines as A
        assert_raises(ValueError, place_poles, A, B[:3,:], (-2,-3,-4,-5))

        # should fail as KNV0 does not support complex poles
        assert_raises(ValueError, place_poles, A, B,
                      (-2+1j,-2-1j,-2+3j,-2-3j), method="KNV0")


class TestSS2TF:

    def tst_matrix_shapes(self, p, q, r):
        ss2tf(np.zeros((p, p)),
              np.zeros((p, q)),
              np.zeros((r, p)),
              np.zeros((r, q)), 0)

    def test_shapes(self):
        # Each tuple holds:
        #   number of states, number of inputs, number of outputs
        for p, q, r in [(3, 3, 3), (1, 3, 3), (1, 1, 1)]:
            yield self.tst_matrix_shapes, p, q, r

    def test_basic(self):
        # Test a round trip through tf2ss and sst2f.
        b = np.array([1.0, 3.0, 5.0])
        a = np.array([1.0, 2.0, 3.0])

        A, B, C, D = tf2ss(b, a)
        assert_allclose(A, [[-2, -3], [1, 0]], rtol=1e-13)
        assert_allclose(B, [[1], [0]], rtol=1e-13)
        assert_allclose(C, [[1, 2]], rtol=1e-13)
        # N.b. the shape of D returned by tf2ss is (1,), not (1, 1). Sigh.
        assert_allclose(D, [1], rtol=1e-14)

        bb, aa = ss2tf(A, B, C, D)
        assert_allclose(bb[0], b, rtol=1e-13)
        assert_allclose(aa, a, rtol=1e-13)

    def test_multioutput(self):
        # Regression test for gh-2669.

        # 4 states
        A = np.array([[-1.0, 0.0, 1.0, 0.0],
                      [-1.0, 0.0, 2.0, 0.0],
                      [-4.0, 0.0, 3.0, 0.0],
                      [-8.0, 8.0, 0.0, 4.0]])

        # 1 input
        B = np.array([[0.3],
                      [0.0],
                      [7.0],
                      [0.0]])

        # 3 outputs
        C = np.array([[0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [8.0, 8.0, 0.0, 0.0]])

        D = np.array([[0.0],
                      [0.0],
                      [1.0]])

        # Get the transfer functions for all the outputs in one call.
        b_all, a = ss2tf(A, B, C, D)

        # Get the transfer functions for each output separately.
        b0, a0 = ss2tf(A, B, C[0], D[0])
        b1, a1 = ss2tf(A, B, C[1], D[1])
        b2, a2 = ss2tf(A, B, C[2], D[2])

        # Check that we got the same results.
        assert_allclose(a0, a, rtol=1e-13)
        assert_allclose(a1, a, rtol=1e-13)
        assert_allclose(a2, a, rtol=1e-13)
        assert_allclose(b_all, np.vstack((b0, b1, b2)), rtol=1e-13, atol=1e-14)


class Test_lsim2(object):

    def test_01(self):
        t = np.linspace(0,10,1001)
        u = np.zeros_like(t)
        # First order system: x'(t) + x(t) = u(t), x(0) = 1.
        # Exact solution is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0])
        expected_x = np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_02(self):
        t = np.array([0.0, 1.0, 1.0, 3.0])
        u = np.array([0.0, 0.0, 1.0, 1.0])
        # Simple integrator: x'(t) = u(t)
        system = ([1.0],[1.0,0.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0])
        expected_x = np.maximum(1.0, tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_03(self):
        t = np.array([0.0, 1.0, 1.0, 1.1, 1.1, 2.0])
        u = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        # Simple integrator:  x'(t) = u(t)
        system = ([1.0],[1.0, 0.0])
        tout, y, x = lsim2(system, u, t, hmax=0.01)
        expected_x = np.array([0.0, 0.0, 0.0, 0.1, 0.1, 0.1])
        assert_almost_equal(x[:,0], expected_x)

    def test_04(self):
        t = np.linspace(0, 10, 1001)
        u = np.zeros_like(t)
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = 0.
        # With initial conditions x(0)=1.0 and x'(t)=0.0, the exact solution
        # is (1-t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0, 0.0])
        expected_x = (1.0 - tout) * np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_05(self):
        # The call to lsim2 triggers a "BadCoefficients" warning from
        # scipy.signal.filter_design, but the test passes.  I think the warning
        # is related to the incomplete handling of multi-input systems in
        # scipy.signal.

        # A system with two state variables, two inputs, and one output.
        A = np.array([[-1.0, 0.0], [0.0, -2.0]])
        B = np.array([[1.0, 0.0], [0.0, 1.0]])
        C = np.array([1.0, 0.0])
        D = np.zeros((1,2))

        t = np.linspace(0, 10.0, 101)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BadCoefficients)
            tout, y, x = lsim2((A,B,C,D), T=t, X0=[1.0, 1.0])
        expected_y = np.exp(-tout)
        expected_x0 = np.exp(-tout)
        expected_x1 = np.exp(-2.0*tout)
        assert_almost_equal(y, expected_y)
        assert_almost_equal(x[:,0], expected_x0)
        assert_almost_equal(x[:,1], expected_x1)

    def test_06(self):
        """Test use of the default values of the arguments `T` and `U`."""
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = 0.
        # With initial conditions x(0)=1.0 and x'(t)=0.0, the exact solution
        # is (1-t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y, x = lsim2(system, X0=[1.0, 0.0])
        expected_x = (1.0 - tout) * np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)


class _TestImpulseFuncs(object):
    # Common tests for impulse/impulse2 (= self.func)

    def test_01(self):
        # First order system: x'(t) + x(t) = u(t)
        # Exact impulse response is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system)
        expected_y = np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_02(self):
        # Specify the desired time values for the output.

        # First order system: x'(t) + x(t) = u(t)
        # Exact impulse response is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        n = 21
        t = np.linspace(0, 2.0, n)
        tout, y = self.func(system, T=t)
        assert_equal(tout.shape, (n,))
        assert_almost_equal(tout, t)
        expected_y = np.exp(-t)
        assert_almost_equal(y, expected_y)

    def test_03(self):
        # Specify an initial condition as a scalar.

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact impulse response is x(t) = 4*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system, X0=3.0)
        expected_y = 4.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_04(self):
        # Specify an initial condition as a list.

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact impulse response is x(t) = 4*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system, X0=[3.0])
        expected_y = 4.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_05(self):
        # Simple integrator: x'(t) = u(t)
        system = ([1.0],[1.0,0.0])
        tout, y = self.func(system)
        expected_y = np.ones_like(tout)
        assert_almost_equal(y, expected_y)

    def test_array_like(self):
        # Test that function can accept sequences, scalars.
        system = ([1.0], [1.0, 2.0, 1.0])
        # TODO: add meaningful test where X0 is a list
        tout, y = self.func(system, X0=[3], T=[5, 6])
        tout, y = self.func(system, X0=[3], T=[5])


class TestImpulse2(_TestImpulseFuncs):
    def setup(self):
        self.func = impulse2

    def test_array_like2(self):
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y = self.func(system, X0=3, T=5)

    def test_06(self):
        # Second order system with a repeated root:
        #     x''(t) + 2*x(t) + x(t) = u(t)
        # The exact impulse response is t*exp(-t).
        # Doesn't pass for `impulse` (on some systems, see gh-2654)
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y = self.func(system)
        expected_y = tout * np.exp(-tout)
        assert_almost_equal(y, expected_y)


class TestImpulse(_TestImpulseFuncs):
    def setup(self):
        self.func = impulse


class _TestStepFuncs(object):
    def test_01(self):
        # First order system: x'(t) + x(t) = u(t)
        # Exact step response is x(t) = 1 - exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system)
        expected_y = 1.0 - np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_02(self):
        # Specify the desired time values for the output.

        # First order system: x'(t) + x(t) = u(t)
        # Exact step response is x(t) = 1 - exp(-t).
        system = ([1.0],[1.0,1.0])
        n = 21
        t = np.linspace(0, 2.0, n)
        tout, y = self.func(system, T=t)
        assert_equal(tout.shape, (n,))
        assert_almost_equal(tout, t)
        expected_y = 1 - np.exp(-t)
        assert_almost_equal(y, expected_y)

    def test_03(self):
        # Specify an initial condition as a scalar.

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact step response is x(t) = 1 + 2*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system, X0=3.0)
        expected_y = 1 + 2.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_04(self):
        # Specify an initial condition as a list.

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact step response is x(t) = 1 + 2*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = self.func(system, X0=[3.0])
        expected_y = 1 + 2.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_array_like(self):
        # Test that function can accept sequences, scalars.
        system = ([1.0], [1.0, 2.0, 1.0])
        # TODO: add meaningful test where X0 is a list
        tout, y = self.func(system, T=[5, 6])


class TestStep2(_TestStepFuncs):
    def setup(self):
        self.func = step2

    def test_05(self):
        # Simple integrator: x'(t) = u(t)
        # Exact step response is x(t) = t.
        system = ([1.0],[1.0,0.0])
        tout, y = self.func(system, atol=1e-10, rtol=1e-8)
        expected_y = tout
        assert_almost_equal(y, expected_y)

    def test_06(self):
        # Second order system with a repeated root:
        #     x''(t) + 2*x(t) + x(t) = u(t)
        # The exact step response is 1 - (1 + t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y = self.func(system, atol=1e-10, rtol=1e-8)
        expected_y = 1 - (1 + tout) * np.exp(-tout)
        assert_almost_equal(y, expected_y)


class TestStep(_TestStepFuncs):
    def setup(self):
        self.func = step

    def test_complex_input(self):
        # Test that complex input doesn't raise an error.
        # `step` doesn't seem to have been designed for complex input, but this
        # works and may be used, so add regression test.  See gh-2654.
        step(([], [-1], 1+0j))


def test_lti_instantiation():
    # Test that lti can be instantiated with sequences, scalars.  See PR-225.
    s = lti([1], [-1])
    s = lti(np.array([]), np.array([-1]), 1)
    s = lti([], [-1], 1)
    s = lti([1], [-1], 1, 3)


class Test_abcd_normalize(object):
    def setup(self):
        self.A = np.array([[1.0, 2.0], [3.0, 4.0]])
        self.B = np.array([[-1.0], [5.0]])
        self.C = np.array([[4.0, 5.0]])
        self.D = np.array([[2.5]])

    def test_no_matrix_fails(self):
        assert_raises(ValueError, abcd_normalize)

    def test_A_nosquare_fails(self):
        assert_raises(ValueError, abcd_normalize, [1, -1],
                      self.B, self.C, self.D)

    def test_AB_mismatch_fails(self):
        assert_raises(ValueError, abcd_normalize, self.A, [-1, 5],
                      self.C, self.D)

    def test_AC_mismatch_fails(self):
        assert_raises(ValueError, abcd_normalize, self.A, self.B,
                      [[4.0], [5.0]], self.D)

    def test_CD_mismatch_fails(self):
        assert_raises(ValueError, abcd_normalize, self.A, self.B,
                      self.C, [2.5, 0])

    def test_BD_mismatch_fails(self):
        assert_raises(ValueError, abcd_normalize, self.A, [-1, 5],
                      self.C, self.D)

    def test_normalized_matrices_unchanged(self):
        A, B, C, D = abcd_normalize(self.A, self.B, self.C, self.D)
        assert_equal(A, self.A)
        assert_equal(B, self.B)
        assert_equal(C, self.C)
        assert_equal(D, self.D)

    def test_shapes(self):
        A, B, C, D = abcd_normalize(self.A, self.B, [1, 0], 0)
        assert_equal(A.shape[0], A.shape[1])
        assert_equal(A.shape[0], B.shape[0])
        assert_equal(A.shape[0], C.shape[1])
        assert_equal(C.shape[0], D.shape[0])
        assert_equal(B.shape[1], D.shape[1])

    def test_zero_dimension_is_not_none1(self):
        B_ = np.zeros((2, 0))
        D_ = np.zeros((0, 0))
        A, B, C, D = abcd_normalize(A=self.A, B=B_, D=D_)
        assert_equal(A, self.A)
        assert_equal(B, B_)
        assert_equal(D, D_)
        assert_equal(C.shape[0], D_.shape[0])
        assert_equal(C.shape[1], self.A.shape[0])

    def test_zero_dimension_is_not_none2(self):
        B_ = np.zeros((2, 0))
        C_ = np.zeros((0, 2))
        A, B, C, D = abcd_normalize(A=self.A, B=B_, C=C_)
        assert_equal(A, self.A)
        assert_equal(B, B_)
        assert_equal(C, C_)
        assert_equal(D.shape[0], C_.shape[0])
        assert_equal(D.shape[1], B_.shape[1])

    def test_missing_A(self):
        A, B, C, D = abcd_normalize(B=self.B, C=self.C, D=self.D)
        assert_equal(A.shape[0], A.shape[1])
        assert_equal(A.shape[0], B.shape[0])
        assert_equal(A.shape, (self.B.shape[0], self.B.shape[0]))

    def test_missing_B(self):
        A, B, C, D = abcd_normalize(A=self.A, C=self.C, D=self.D)
        assert_equal(B.shape[0], A.shape[0])
        assert_equal(B.shape[1], D.shape[1])
        assert_equal(B.shape, (self.A.shape[0], self.D.shape[1]))

    def test_missing_C(self):
        A, B, C, D = abcd_normalize(A=self.A, B=self.B, D=self.D)
        assert_equal(C.shape[0], D.shape[0])
        assert_equal(C.shape[1], A.shape[0])
        assert_equal(C.shape, (self.D.shape[0], self.A.shape[0]))

    def test_missing_D(self):
        A, B, C, D = abcd_normalize(A=self.A, B=self.B, C=self.C)
        assert_equal(D.shape[0], C.shape[0])
        assert_equal(D.shape[1], B.shape[1])
        assert_equal(D.shape, (self.C.shape[0], self.B.shape[1]))

    def test_missing_AB(self):
        A, B, C, D = abcd_normalize(C=self.C, D=self.D)
        assert_equal(A.shape[0], A.shape[1])
        assert_equal(A.shape[0], B.shape[0])
        assert_equal(B.shape[1], D.shape[1])
        assert_equal(A.shape, (self.C.shape[1], self.C.shape[1]))
        assert_equal(B.shape, (self.C.shape[1], self.D.shape[1]))

    def test_missing_AC(self):
        A, B, C, D = abcd_normalize(B=self.B, D=self.D)
        assert_equal(A.shape[0], A.shape[1])
        assert_equal(A.shape[0], B.shape[0])
        assert_equal(C.shape[0], D.shape[0])
        assert_equal(C.shape[1], A.shape[0])
        assert_equal(A.shape, (self.B.shape[0], self.B.shape[0]))
        assert_equal(C.shape, (self.D.shape[0], self.B.shape[0]))

    def test_missing_AD(self):
        A, B, C, D = abcd_normalize(B=self.B, C=self.C)
        assert_equal(A.shape[0], A.shape[1])
        assert_equal(A.shape[0], B.shape[0])
        assert_equal(D.shape[0], C.shape[0])
        assert_equal(D.shape[1], B.shape[1])
        assert_equal(A.shape, (self.B.shape[0], self.B.shape[0]))
        assert_equal(D.shape, (self.C.shape[0], self.B.shape[1]))

    def test_missing_BC(self):
        A, B, C, D = abcd_normalize(A=self.A, D=self.D)
        assert_equal(B.shape[0], A.shape[0])
        assert_equal(B.shape[1], D.shape[1])
        assert_equal(C.shape[0], D.shape[0])
        assert_equal(C.shape[1], A.shape[0])
        assert_equal(B.shape, (self.A.shape[0], self.D.shape[1]))
        assert_equal(C.shape, (self.D.shape[0], self.A.shape[0]))

    def test_missing_ABC_fails(self):
        assert_raises(ValueError, abcd_normalize, D=self.D)

    def test_missing_BD_fails(self):
        assert_raises(ValueError, abcd_normalize, A=self.A, C=self.C)

    def test_missing_CD_fails(self):
        assert_raises(ValueError, abcd_normalize, A=self.A, B=self.B)


class Test_bode(object):

    def test_01(self):
        """Test bode() magnitude calculation (manual sanity check)."""
        # 1st order low-pass filter: H(s) = 1 / (s + 1),
        # cutoff: 1 rad/s, slope: -20 dB/decade
        #   H(s=0.1) ~= 0 dB
        #   H(s=1) ~= -3 dB
        #   H(s=10) ~= -20 dB
        #   H(s=100) ~= -40 dB
        system = lti([1], [1, 1])
        w = [0.1, 1, 10, 100]
        w, mag, phase = bode(system, w=w)
        expected_mag = [0, -3, -20, -40]
        assert_almost_equal(mag, expected_mag, decimal=1)

    def test_02(self):
        """Test bode() phase calculation (manual sanity check)."""
        # 1st order low-pass filter: H(s) = 1 / (s + 1),
        #   angle(H(s=0.1)) ~= -5.7 deg
        #   angle(H(s=1)) ~= -45 deg
        #   angle(H(s=10)) ~= -84.3 deg
        system = lti([1], [1, 1])
        w = [0.1, 1, 10]
        w, mag, phase = bode(system, w=w)
        expected_phase = [-5.7, -45, -84.3]
        assert_almost_equal(phase, expected_phase, decimal=1)

    def test_03(self):
        """Test bode() magnitude calculation."""
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        system = lti([1], [1, 1])
        w = [0.1, 1, 10, 100]
        w, mag, phase = bode(system, w=w)
        jw = w * 1j
        y = np.polyval(system.num, jw) / np.polyval(system.den, jw)
        expected_mag = 20.0 * np.log10(abs(y))
        assert_almost_equal(mag, expected_mag)

    def test_04(self):
        """Test bode() phase calculation."""
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        system = lti([1], [1, 1])
        w = [0.1, 1, 10, 100]
        w, mag, phase = bode(system, w=w)
        jw = w * 1j
        y = np.polyval(system.num, jw) / np.polyval(system.den, jw)
        expected_phase = np.arctan2(y.imag, y.real) * 180.0 / np.pi
        assert_almost_equal(phase, expected_phase)

    def test_05(self):
        """Test that bode() finds a reasonable frequency range."""
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        system = lti([1], [1, 1])
        n = 10
        # Expected range is from 0.01 to 10.
        expected_w = np.logspace(-2, 1, n)
        w, mag, phase = bode(system, n=n)
        assert_almost_equal(w, expected_w)

    def test_06(self):
        """Test that bode() doesn't fail on a system with a pole at 0."""
        # integrator, pole at zero: H(s) = 1 / s
        system = lti([1], [1, 0])
        w, mag, phase = bode(system, n=2)
        assert_equal(w[0], 0.01)  # a fail would give not-a-number

    def test_07(self):
        """bode() should not fail on a system with pure imaginary poles."""
        # The test passes if bode doesn't raise an exception.
        system = lti([1], [1, 0, 100])
        w, mag, phase = bode(system, n=2)

    def test_08(self):
        """Test that bode() return continuous phase, issues/2331."""
        system = lti([], [-10, -30, -40, -60, -70], 1)
        w, mag, phase = system.bode(w=np.logspace(-3, 40, 100))
        assert_almost_equal(min(phase), -450, decimal=15)

    def test_from_state_space(self):
        # Ensure that bode works with a system that was created from the
        # state space representation matrices A, B, C, D.  In this case,
        # system.num will be a 2-D array with shape (1, n+1), where (n,n)
        # is the shape of A.
        # A Butterworth lowpass filter is used, so we know the exact
        # frequency response.
        a = np.array([1.0, 2.0, 2.0, 1.0])
        A = linalg.companion(a).T
        B = np.array([[0.0],[0.0],[1.0]])
        C = np.array([[1.0, 0.0, 0.0]])
        D = np.array([[0.0]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BadCoefficients)
            system = lti(A, B, C, D)
        w, mag, phase = bode(system, n=100)
        expected_magnitude = 20 * np.log10(np.sqrt(1.0 / (1.0 + w**6)))
        assert_almost_equal(mag, expected_magnitude)


class Test_freqresp(object):

    def test_real_part_manual(self):
        # Test freqresp() real part calculation (manual sanity check).
        # 1st order low-pass filter: H(s) = 1 / (s + 1),
        #   re(H(s=0.1)) ~= 0.99
        #   re(H(s=1)) ~= 0.5
        #   re(H(s=10)) ~= 0.0099
        system = lti([1], [1, 1])
        w = [0.1, 1, 10]
        w, H = freqresp(system, w=w)
        expected_re = [0.99, 0.5, 0.0099]
        assert_almost_equal(H.real, expected_re, decimal=1)

    def test_imag_part_manual(self):
        # Test freqresp() imaginary part calculation (manual sanity check).
        # 1st order low-pass filter: H(s) = 1 / (s + 1),
        #   im(H(s=0.1)) ~= -0.099
        #   im(H(s=1)) ~= -0.5
        #   im(H(s=10)) ~= -0.099
        system = lti([1], [1, 1])
        w = [0.1, 1, 10]
        w, H = freqresp(system, w=w)
        expected_im = [-0.099, -0.5, -0.099]
        assert_almost_equal(H.imag, expected_im, decimal=1)

    def test_real_part(self):
        # Test freqresp() real part calculation.
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        system = lti([1], [1, 1])
        w = [0.1, 1, 10, 100]
        w, H = freqresp(system, w=w)
        jw = w * 1j
        y = np.polyval(system.num, jw) / np.polyval(system.den, jw)
        expected_re = y.real
        assert_almost_equal(H.real, expected_re)

    def test_imag_part(self):
        # Test freqresp() imaginary part calculation.
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        system = lti([1], [1, 1])
        w = [0.1, 1, 10, 100]
        w, H = freqresp(system, w=w)
        jw = w * 1j
        y = np.polyval(system.num, jw) / np.polyval(system.den, jw)
        expected_im = y.imag
        assert_almost_equal(H.imag, expected_im)

    def test_freq_range(self):
        # Test that freqresp() finds a reasonable frequency range.
        # 1st order low-pass filter: H(s) = 1 / (s + 1)
        # Expected range is from 0.01 to 10.
        system = lti([1], [1, 1])
        n = 10
        expected_w = np.logspace(-2, 1, n)
        w, H = freqresp(system, n=n)
        assert_almost_equal(w, expected_w)

    def test_pole_zero(self):
        # Test that freqresp() doesn't fail on a system with a pole at 0.
        # integrator, pole at zero: H(s) = 1 / s
        system = lti([1], [1, 0])
        w, H = freqresp(system, n=2)
        assert_equal(w[0], 0.01)  # a fail would give not-a-number

    def test_from_state_space(self):
        # Ensure that freqresp works with a system that was created from the
        # state space representation matrices A, B, C, D.  In this case,
        # system.num will be a 2-D array with shape (1, n+1), where (n,n) is
        # the shape of A.
        # A Butterworth lowpass filter is used, so we know the exact
        # frequency response.
        a = np.array([1.0, 2.0, 2.0, 1.0])
        A = linalg.companion(a).T
        B = np.array([[0.0],[0.0],[1.0]])
        C = np.array([[1.0, 0.0, 0.0]])
        D = np.array([[0.0]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BadCoefficients)
            system = lti(A, B, C, D)
        w, H = freqresp(system, n=100)
        expected_magnitude = np.sqrt(1.0 / (1.0 + w**6))
        assert_almost_equal(np.abs(H), expected_magnitude)


if __name__ == "__main__":
    run_module_suite()
