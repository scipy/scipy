"""
Unit tests for trust-region subproblem solver trlib.

To run it in its simplest form::
  nosetests test_optimize.py

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.optimize._trlib import (TRLIBQuadraticSubproblem)
from numpy.testing import (TestCase, assert_, assert_array_equal,
                           assert_almost_equal,
                           assert_equal, assert_array_almost_equal,
                           assert_array_less, run_module_suite)

TRLIBQP = TRLIBQuadraticSubproblem(tol_rel_i=1e-8, tol_rel_b=1e-6)

class TestTRLIBQuadraticSubproblem(TestCase):

    def test_for_the_easy_case(self):

        # `H` is chosen such that `g` is not orthogonal to the
        # eigenvector associated with the smallest eigenvalue.
        H = np.array([[1.0, 0.0, 4.0],
             [0.0, 2.0, 0.0],
             [4.0, 0.0, 3.0]])
        g = np.array([5.0, 0.0, 4.0])

        # Trust Radius
        trust_radius = 1.0

        # Solve Subproblem
        subprob = TRLIBQP(x=0,
                          fun=lambda x: 0,
                          jac=lambda x: g,
                          hess=lambda x: None,
                          hessp=lambda x, y: H.dot(y))
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(p, np.array([-1.0, 0.0, 0.0]))
        assert_equal(hits_boundary, True)
        # check kkt satisfaction
        assert_almost_equal(
                np.linalg.norm(H.dot(p) + subprob.lam * p + g),
                0.0)
        # check trust region constraint
        assert_almost_equal(np.linalg.norm(p), trust_radius)

        trust_radius = 0.5
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(p,
                np.array([-0.46125446, 0., -0.19298788]))
        assert_equal(hits_boundary, True)
        # check kkt satisfaction
        assert_almost_equal(
                np.linalg.norm(H.dot(p) + subprob.lam * p + g),
                0.0)
        # check trust region constraint
        assert_almost_equal(np.linalg.norm(p), trust_radius)

    def test_for_the_hard_case(self):

        # `H` is chosen such that `g` is orthogonal to the
        # eigenvector associated with the smallest eigenvalue.
        H = np.array([[1.0, 0.0, 4.0],
             [0.0, 2.0, 0.0],
             [4.0, 0.0, 3.0]])
        g = np.array([0.0, 2.0, 0.0])

        # Trust Radius
        trust_radius = 1.0

        # Solve Subproblem
        subprob = TRLIBQP(x=0,
                          fun=lambda x: 0,
                          jac=lambda x: g,
                          hess=lambda x: None,
                          hessp=lambda x, y: H.dot(y))
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(p, np.array([0.0, -1.0, 0.0]))
        # check kkt satisfaction
        assert_almost_equal(
                np.linalg.norm(H.dot(p) + subprob.lam * p + g),
                0.0)
        # check trust region constraint
        assert_almost_equal(np.linalg.norm(p), trust_radius)

        trust_radius = 0.5
        p, hits_boundary = subprob.solve(trust_radius)

        assert_array_almost_equal(p, np.array([0.0, -0.5, 0.0]))
        # check kkt satisfaction
        assert_almost_equal(
                np.linalg.norm(H.dot(p) + subprob.lam * p + g),
                0.0)
        # check trust region constraint
        assert_almost_equal(np.linalg.norm(p), trust_radius)

if __name__ == '__main__':
    run_module_suite()
