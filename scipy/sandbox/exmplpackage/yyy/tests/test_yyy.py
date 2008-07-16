
from numpy.testing import *

from yyy import fun

class TestFun(TestCase):
    def test_simple(self):
        assert fun()=='Hello from yyy.fun'
    #...

if __name__ == "__main__":
    nose.run(argv=['', __file__])
