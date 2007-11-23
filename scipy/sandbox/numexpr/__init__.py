from numexpr.info import __doc__
from numexpr.expressions import E
from numexpr.compiler import numexpr, disassemble, evaluate

def test(level=1, verbosity=1):
    from numpy.testing import NumpyTest
    NumpyTest().test(level, verbosity)
