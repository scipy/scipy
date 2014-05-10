from __future__ import absolute_import, print_function

from numpy.testing import TestCase, assert_, run_module_suite

from scipy.weave import inline_tools

from weave_test_utils import dec


class TestInline(TestCase):
    """These are long running tests...

    Would be useful to benchmark these things somehow.
    """
    @dec.slow
    def test_exceptions(self):
        a = 3
        code = """
               if (a < 2)
                  throw_error(PyExc_ValueError,
                              "the variable 'a' should not be less than 2");
               else
                   return_val = PyInt_FromLong(a+1);
               """
        result = inline_tools.inline(code,['a'])
        assert_(result == 4)


## Unfortunately, it is not always possible to catch distutils compiler
## errors, since SystemExit is used.  Until that is fixed, these tests
## cannot be run in the same process as the test suite.

##         try:
##             a = 1
##             result = inline_tools.inline(code,['a'])
##             assert_(1) # should've thrown a ValueError
##         except ValueError:
##             pass

##         from distutils.errors import DistutilsError, CompileError
##         try:
##             a = 'string'
##             result = inline_tools.inline(code,['a'])
##             assert_(1) # should've gotten an error
##         except:
##             # ?CompileError is the error reported, but catching it doesn't work
##             pass


if __name__ == "__main__":
    run_module_suite()
