
from numpy.testing import *

from scipy.sandbox.fdfpack import periodic_finite_difference as diff


from numpy import arange, add, array,sin,cos,pi

class TestDiff(TestCase):
    def test_1(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(5,15):
                assert_array_almost_equal(diff(sin(x),m=m),cos(x))
    def test_2(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(8,15)+[6]:
                assert_array_almost_equal(diff(sin(x),k=2,m=m),-sin(x))
    def test_3(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(7,15):
                assert_array_almost_equal(diff(sin(x),k=3,m=m)/n,-cos(x)/n)
    def test_4(self):
        for n in [64,100,125,200,2000]:
            x = arange(n)*2*pi/n
            for m in range(10,15):
                assert_array_almost_equal(diff(sin(x),k=4,m=m)/n,sin(x)/n)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
