from __future__ import absolute_import, print_function

from numpy.testing import TestCase, assert_equal, run_module_suite

from scipy.weave import ast_tools


class TestHarvestVariables(TestCase):
    """ Not much testing going on here, but at least it is a flame test."""
    def generic_check(self,expr,desired):
        import parser
        ast_list = parser.suite(expr).tolist()
        actual = ast_tools.harvest_variables(ast_list)
        assert_equal(actual,desired,expr)

    def test_simple_expr(self):
        # Convert simple expr to blitz
        expr = "a[:1:2] = b[:1+i+2:]"
        desired = ['a','b','i']
        self.generic_check(expr,desired)

if __name__ == "__main__":
    run_module_suite()
