from info import __doc__
from expressions import E
from compiler import numexpr, disassemble, evaluate

def test(level=1, verbosity=1):
    from numpy.testing import NumpyTest
    NumpyTest().test(level, verbosity)
