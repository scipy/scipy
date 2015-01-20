from __future__ import division, print_function, absolute_import

from os.path import join as pjoin, dirname
import sys

from numpy.testing import assert_array_equal, run_module_suite

from scipy.io.matlab.mio import loadmat


test_data_path = pjoin(dirname(__file__), "data")


def test_load():
    for ver in ["6"]: # 7, 7.3 not supported
        data = loadmat(pjoin(test_data_path, "testclass_%s.mat" % ver))
        assert "obj1" in data
        assert "obj2" in data
        assert "obj3" in data


if __name__ == "__main__":
    run_module_suite(argv=sys.argv)
