from __future__ import division, print_function, absolute_import

from os.path import join as pjoin, dirname
import sys

import numpy as np
from numpy.testing import assert_array_equal, run_module_suite

from scipy.io.matlab.mio import loadmat


test_data_path = pjoin(dirname(__file__), "data")


def test_load():
    for ver in ["v6", "v7"]: # 7.3 not supported
        data = loadmat(pjoin(test_data_path, "testclass_%s.mat" % ver))

        assert_array_equal(data["obj1"]["char_field_1"], ["char_field"])
        assert_array_equal(data["obj1"]["array_field_2"].item().item(),
                           np.arange(6.)[None])
        field_item = data["obj1"]["any_field_3"].item().item()
        assert_array_equal(field_item["char_field_1"].item().item(),
                           ["char_field"])
        assert_array_equal(field_item["array_field_2"].item().item(),
                           np.arange(6.)[None])
        assert_array_equal(field_item["any_field_3"], [[0]])

        assert_array_equal(data["obj2"]["char_field_1"], ["foo"])
        assert_array_equal(data["obj2"]["array_field_2"].item().item(),
                           np.arange(6.)[None])
        field_item = data["obj2"]["any_field_3"].item().item()
        assert_array_equal(field_item["a"][0, 0], [[1.]])
        assert_array_equal(field_item[0, 0], field_item[0, 1])


if __name__ == "__main__":
    run_module_suite(argv=sys.argv)
