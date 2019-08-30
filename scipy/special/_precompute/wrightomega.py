from __future__ import division, print_function, absolute_import

import numpy as np

try:
    import mpmath
    import matplotlib.pyplot as plt
except ImportError:
    pass


def mpmath_wrightomega(x):
    return mpmath.lambertw(mpmath.exp(x), mpmath.mpf('-0.5'))


def wrightomega_series_error(x):
    series = x
    desired = mpmath_wrightomega(x)
    return abs(series - desired) / desired


def main():
    desired_error = 2 * np.finfo(float).eps
    for x in [1e5, 1e10, 1e15, 1e20]:
        error = wrightomega_series_error(x)
        print(error, error < desired_error)


if __name__ == '__main__':
    main()
