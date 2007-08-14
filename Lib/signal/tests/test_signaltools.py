#this program corresponds to special.py

from numpy.testing import *
set_package_path()
import scipy.signal as signal
restore_path()

from numpy import array, arange

class test_convolve(NumpyTestCase):
    def check_basic(self):
        a = [3,4,5,6,5,4]
        b = [1,2,3]
        c = signal.convolve(a,b)
        assert_array_equal(c,array([3,10,22,28,32,32,23,12]))

class test_medfilt(NumpyTestCase):
    def check_basic(self):
        f = [[3,4,5],[2,3,4],[1,2,5]]
        d = signal.medfilt(f)
        assert_array_equal(d, [[0,3,0],[2,3,3],[0,2,0]])

class test_wiener(NumpyTestCase):
    def check_basic(self):
        g = array([[5,6,4,3],[3,5,6,2],[2,3,5,6],[1,6,9,7]],'d')
        correct = array([[2.16374269,3.2222222222, 2.8888888889, 1.6666666667],[2.666666667, 4.33333333333, 4.44444444444, 2.8888888888],[2.222222222, 4.4444444444, 5.4444444444, 4.801066874837],[1.33333333333, 3.92735042735, 6.0712560386, 5.0404040404]])
        h = signal.wiener(g)
        assert_array_almost_equal(h,correct,decimal=6)

class test_cspline1d_eval(NumpyTestCase):
    def check_basic(self):
        y=array([1,2,3,4,3,2,1,2,3.0])
        x=arange(len(y))
        dx=x[1]-x[0]
        cj = signal.cspline1d(y)

        x2=arange(len(y)*10.0)/10.0
        y2=signal.cspline1d_eval(cj, x2, dx=dx,x0=x[0])

        # make sure interpolated values are on knot points
        assert_array_almost_equal(y2[::10], y, decimal=5)


if __name__ == "__main__":
    NumpyTest().run()
