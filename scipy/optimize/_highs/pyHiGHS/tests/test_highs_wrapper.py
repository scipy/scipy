'''
To run tests from project root directory:

    python -m pyHiGHS.tests.test_highs_wrapper
'''

import unittest

import numpy as np
from scipy.sparse import csc_matrix

from _highs.highs_wrapper import highs_wrapper


class TestHighsWrapperCoverage(unittest.TestCase):
    '''
    Smoke tests to try to cover as many lines of code
    as possible.
    '''

    def setUp(self):
        '''Default values for simple testing'''
        self.m, self.n = 2, 3
        self.c = np.zeros(self.n)
        self.A = csc_matrix(np.ones((self.m, self.n)))
        self.lhs = np.zeros(self.m)
        self.rhs = np.zeros(self.m)
        self.lb = np.zeros(self.n)
        self.ub = np.zeros(self.n)
        self.options = {'message_level': 0}

    def test_yes_rows_cols(self):
        '''
        Test numrow > 0 branch and numcol > 0 branch.
        Also tests nonempty A.
        '''
        assert(self.A.shape[0] > 0)  # numrow > 0
        assert(self.A.shape[1] > 0)  # numcol > 0
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            self.options)

    def test_no_rows(self):
        '''Test numrow == 0 branch'''
        A = csc_matrix(np.ones((0, self.n)))
        lhs, rhs = np.zeros(0), np.zeros(0)
        assert(lhs.size == 0)  # numrow == 0
        assert(rhs.size == 0)  # numrow == 0
        assert(A.shape[0] == 0)  # numrow == 0
        highs_wrapper(
            self.c,
            A.indptr,
            A.indices,
            A.data,
            lhs,
            rhs,
            self.lb,
            self.ub,
            self.options)

    def test_no_cols(self):
        '''Test numcol == 0 branch'''
        c = np.zeros(0)
        A = csc_matrix(np.ones((self.m, 0)))
        lb, ub = np.zeros(0), np.zeros(0)
        assert(c.size == 0)
        assert(A.shape[1] == 0)
        assert(lb.size == 0)
        assert(ub.size == 0)
        highs_wrapper(
            c,
            A.indptr,
            A.indices,
            A.data,
            self.lhs,
            self.rhs,
            lb,
            ub,
            self.options)

    def test_empty_A(self):
        '''Test nnz == 0 branch'''
        A = csc_matrix(np.zeros(self.A.shape))
        assert(A.nnz == 0)  # nnz == 0
        highs_wrapper(
            self.c,
            A.indptr,
            A.indices,
            A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            self.options)

    def test_message_level_NONE(self):
        '''Test message_level == NONE branch'''
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {'message_level': 0})

    def test_message_level_not_NONE(self):
        '''Test message_level != NONE branch'''
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {'message_level': 1})

    def test_opt_warning_bool(self):
        '''Test opt_warning bool branch'''
        opt = 'less_infeasible_DSE_check'
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            options={opt: np.pi})

    def test_opt_warning_bool2string(self):
        '''Test opt_warning bool to string branch'''
        opt = 'presolve'
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {opt: np.pi})

    def test_opt_warning_int(self):
        '''Test opt_warning int branch'''
        opt = 'allowed_simplex_cost_scale_factor'
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {opt: np.pi})

    def test_opt_warning_double(self):
        '''Test opt_warning double branch'''
        opt = 'dual_feasibility_tolerance'
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {opt: "Not a double"})

    def test_opt_warning_string(self):
        '''Test opt_warning string branch'''
        opt = 'solver'
        highs_wrapper(
            self.c,
            self.A.indptr,
            self.A.indices,
            self.A.data,
            self.lhs,
            self.rhs,
            self.lb,
            self.ub,
            {opt: np.pi})

    def test_bad_solution(self):
        '''Test bad solution branch'''
        c = np.array([-1, 1]).astype('double')
        A = csc_matrix([[-1, 2], [3, 4]]).astype('double')
        lhs = np.array([-np.inf, -np.inf]).astype('double')
        rhs = np.array([np.inf, np.inf]).astype('double')
        lb = np.array([-np.inf, -np.inf]).astype('double')
        ub = np.array([np.inf, np.inf]).astype('double')
        res = highs_wrapper(
            c,
            A.indptr,
            A.indices,
            A.data,
            lhs,
            rhs,
            lb,
            ub,
            self.options)
        assert(res['fun'] is None)  # bad solution


if __name__ == '__main__':
    unittest.main()
