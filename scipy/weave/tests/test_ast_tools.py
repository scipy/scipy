from numpy.testing import *

from scipy.weave import ast_tools
from weave_test_utils import *

class TestHarvestVariables(TestCase):
    """ Not much testing going on here, but
        at least it is a flame test.
    """
    def generic_check(self,expr,desired):
        import parser
        ast_list = parser.suite(expr).tolist()
        actual = ast_tools.harvest_variables(ast_list)
        assert_equal(actual,desired,expr)

    def test_simple_expr(self):
        """convert simple expr to blitz

           a[:1:2] = b[:1+i+2:]
        """
        expr = "a[:1:2] = b[:1+i+2:]"
        desired = ['a','b','i']
        self.generic_check(expr,desired)

if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
