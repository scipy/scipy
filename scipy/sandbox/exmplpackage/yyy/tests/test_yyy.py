
import sys
from numpy.test.testing import *

set_package_path()
# make sure that all yyy symbols are imported before the del statement:
from yyy import fun
del sys.path[0]

class TestFun(NumpyTestCase):
    def check_simple(self, level=1):
        assert fun()=='Hello from yyy.fun'
    #...

if __name__ == "__main__":
    NumpyTest().run()
