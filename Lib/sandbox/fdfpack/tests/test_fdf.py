

from numpy.testing import *
set_package_path()
from fdfpack import periodic_finite_difference as diff
restore_path()

from numpy import arange, add, array,sin,cos,pi

class test_diff(NumpyTestCase):
    def check_1(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(5,15):
                assert_array_almost_equal(diff(sin(x),m=m),cos(x))
    def check_2(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(8,15)+[6]:
                assert_array_almost_equal(diff(sin(x),k=2,m=m),-sin(x))
    def check_3(self):
        for n in [64,100,125,4000]:
            x = arange(n)*2*pi/n
            for m in range(7,15):
                assert_array_almost_equal(diff(sin(x),k=3,m=m)/n,-cos(x)/n)
    def check_4(self):
        for n in [64,100,125,200,2000]:
            x = arange(n)*2*pi/n
            for m in range(10,15):
                assert_array_almost_equal(diff(sin(x),k=4,m=m)/n,sin(x)/n)

if __name__ == "__main__":
    NumpyTest().run()
