from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal, assert_array_equal, assert_allclose, \
        run_module_suite

from scipy.interpolate import griddata


class TestGriddata(object):
    def test_fill_value(self):
        x = [(0,0), (0,1), (1,0)]
        y = [1, 2, 3]

        yi = griddata(x, y, [(1,1), (1,2), (0,0)], fill_value=-1)
        assert_array_equal(yi, [-1., -1, 1])

        yi = griddata(x, y, [(1,1), (1,2), (0,0)])
        assert_array_equal(yi, [np.nan, np.nan, 1])

    def test_alternative_call(self):
        x = np.array([(0,0), (-0.5,-0.5), (-0.5,0.5), (0.5, 0.5), (0.25, 0.3)],
                     dtype=np.double)
        y = (np.arange(x.shape[0], dtype=np.double)[:,None]
             + np.array([0,1])[None,:])

        for method in ('nearest', 'linear', 'cubic'):
            yi = griddata((x[:,0], x[:,1]), y, (x[:,0], x[:,1]), method=method)
            assert_allclose(y, yi, atol=1e-14, err_msg=method)

    def test_multivalue_2d(self):
        x = np.array([(0,0), (-0.5,-0.5), (-0.5,0.5), (0.5, 0.5), (0.25, 0.3)],
                     dtype=np.double)
        y = (np.arange(x.shape[0], dtype=np.double)[:,None]
             + np.array([0,1])[None,:])

        for method in ('nearest', 'linear', 'cubic'):
            yi = griddata(x, y, x, method=method)
            assert_allclose(y, yi, atol=1e-14, err_msg=method)

    def test_multipoint_2d(self):
        x = np.array([(0,0), (-0.5,-0.5), (-0.5,0.5), (0.5, 0.5), (0.25, 0.3)],
                     dtype=np.double)
        y = np.arange(x.shape[0], dtype=np.double)

        xi = x[:,None,:] + np.array([0,0,0])[None,:,None]

        for method in ('nearest', 'linear', 'cubic'):
            yi = griddata(x, y, xi, method=method)

            assert_equal(yi.shape, (5, 3), err_msg=method)
            assert_allclose(yi, np.tile(y[:,None], (1, 3)),
                            atol=1e-14, err_msg=method)

    def test_complex_2d(self):
        x = np.array([(0,0), (-0.5,-0.5), (-0.5,0.5), (0.5, 0.5), (0.25, 0.3)],
                     dtype=np.double)
        y = np.arange(x.shape[0], dtype=np.double)
        y = y - 2j*y[::-1]

        xi = x[:,None,:] + np.array([0,0,0])[None,:,None]

        for method in ('nearest', 'linear', 'cubic'):
            yi = griddata(x, y, xi, method=method)

            assert_equal(yi.shape, (5, 3), err_msg=method)
            assert_allclose(yi, np.tile(y[:,None], (1, 3)),
                            atol=1e-14, err_msg=method)

    def test_1d(self):
        x = np.array([1, 2.5, 3, 4.5, 5, 6])
        y = np.array([1, 2, 0, 3.9, 2, 1])

        for method in ('nearest', 'linear', 'cubic'):
            assert_allclose(griddata(x, y, x, method=method), y,
                            err_msg=method, atol=1e-14)
            assert_allclose(griddata(x.reshape(6, 1), y, x, method=method), y,
                            err_msg=method, atol=1e-14)
            assert_allclose(griddata((x,), y, (x,), method=method), y,
                            err_msg=method, atol=1e-14)

    def test_1d_unsorted(self):
        x = np.array([2.5, 1, 4.5, 5, 6, 3])
        y = np.array([1, 2, 0, 3.9, 2, 1])

        for method in ('nearest', 'linear', 'cubic'):
            assert_allclose(griddata(x, y, x, method=method), y,
                            err_msg=method, atol=1e-10)
            assert_allclose(griddata(x.reshape(6, 1), y, x, method=method), y,
                            err_msg=method, atol=1e-10)
            assert_allclose(griddata((x,), y, (x,), method=method), y,
                            err_msg=method, atol=1e-10)


if __name__ == "__main__":
    run_module_suite()
