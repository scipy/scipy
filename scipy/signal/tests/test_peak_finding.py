
import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_almost_equal, assert_array_equal, assert_array_almost_equal, \
    assert_raises, assert_
from scipy.signal._peak_finding import argrelmax


def _gen_gaussians(center_locs, sigmas, total_length):
    num_peaks = len(sigmas)
    xdata = np.arange(0, total_length).astype(float)
    out_data = np.zeros(total_length, dtype=float)
    for ind, sigma in enumerate(sigmas):
        tmp = (xdata - center_locs[ind]) / sigma
        out_data += np.exp(-(tmp ** 2))
    return out_data


def _gen_gaussians_even(sigmas, total_length):
    num_peaks = len(sigmas)
    delta = total_length / (num_peaks + 1)
    center_locs = np.linspace(delta, total_length - delta, num=num_peaks).astype(int)
    out_data = _gen_gaussians(center_locs, sigmas, total_length)
    return out_data, center_locs


class TestArgrelmax(TestCase):

    def test_highorder(self, order=2):
        sigmas = [1.0, 2.0, 10.0, 5.0, 15.0]
        test_data, act_locs = _gen_gaussians_even(sigmas, 500)
        test_data[act_locs + order] = test_data[act_locs] * 0.99999
        test_data[act_locs - order] = test_data[act_locs] * 0.99999
        rel_max_matr = argrelmax(test_data, order=2, mode='clip')
        rel_max_locs = np.where(rel_max_matr)[0]

        assert_(len(rel_max_locs) == len(act_locs))
        assert_((rel_max_locs == act_locs).all())


    def test_2d_gaussians(self):
        sigmas = [1.0, 2.0, 10.0]
        test_data, act_locs = _gen_gaussians_even(sigmas, 100)
        rot_factor = 20
        rot_range = np.arange(0, len(test_data)) - rot_factor
        test_data_2 = np.vstack([test_data, test_data[rot_range]])
        rel_max_matr = argrelmax(test_data_2, axis=1, order=1)
        rel_max_rows, rel_max_cols = np.where(rel_max_matr)

        for rw in xrange(0, test_data_2.shape[0]):
            inds = (rel_max_rows == rw)

            assert_(len(rel_max_cols[inds]) == len(act_locs))
            assert_((act_locs == (rel_max_cols[inds] - rot_factor * rw)).all())


if __name__ == "__main__":
    run_module_suite()