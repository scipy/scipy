import numpy as np

from scipy.integrate._cubature import cub, ProductRule, GaussKronrod21
from numpy.testing import assert_allclose


def test_cub_simple():
    # Currently assumes that f is a function on arrays of shape
    #   (eval_points, dim_input) -> (eval_points, dim_output)
    n = np.arange(10)

    def f(arr):
        x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        return np.power((x + y + z)[:, np.newaxis], n)

    # For example:
    #   f(np.array([
    #       [1,1,1],
    #       [1,2,3],
    #   ]))
    # Outputs:
    #   array([[       1,        3,        9,       27,       81,      243,
    #            729,     2187,     6561,    19683],
    #      [       1,        6,       36,      216,     1296,     7776,
    #          46656,   279936,  1679616, 10077696]])

    exact = (-3 * 2**(n+3) + 3**(n+3) + 3)/(n**3 + 6*n**2 + 11*n + 6)
    rule = ProductRule([GaussKronrod21(), GaussKronrod21(), GaussKronrod21()])

    est, err = cub(
        f,
        np.array([[0, 1], [0, 1], [0, 1]]),
        rule,
        1e-10
    )

    # Outputs:
    #     [  1.           1.5          2.5          4.5          8.6
    #   17.25        36.01190476  77.75       172.73333333 393.3       ]
    #     [7.22078647e-17 1.16982873e-16 1.17539738e-16 1.14560620e-16
    #  7.53115776e-16 1.83786051e-16 1.04254382e-16 5.00106071e-15
    #  1.19304710e-14 1.97483751e-14]

    assert_allclose(est, exact, rtol=0, atol=1e-10, verbose=True)

# TODO: more test integrals
# TODO: test that product rules are calculated properly
# TODO: test that inconsistent dimensions are reported
