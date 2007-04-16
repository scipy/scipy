from numpy.testing import *
set_package_path()
import scipy.misc.pilutil as pilutil
restore_path()

import numpy as N

class test_pilutil(NumpyTestCase):
    def check_imresize(self):
        im = N.random.random((10,20))
        for T in N.sctypes['float'] + [float]:
            im1 = pilutil.imresize(im,T(1.1))
            assert_equal(im1.shape,(11,22))

if __name__ == "__main__":
    NumpyTest().run()
