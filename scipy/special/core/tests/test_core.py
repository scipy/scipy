import os

import numpy as np
from numpy.testing import *
from scipy.special import arccosh, arcsinh, arctanh, erf, erfc, log1p, expm1, \
jn, jv, yn, yv, iv, kv, kn

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

def id(x):
    return x

TEST_MAPPING_1 = {
    'acosh_data.txt': (id, arccosh),
    'asinh_data.txt': (id, arcsinh),
    #assoc_legendre_p.txt
    'atanh_data.txt': (id, arctanh),
    'erf_data.txt': (id, erf, erfc),
    #'erf_inv_data.txt': (id, erf),
    'erf_large_data.txt': (id, erf, erfc),
    'erf_small_data.txt': (id, erf, erfc),
    'log1p_expm1_data.txt': (id, log1p, expm1),
}

TEST_MAPPING_2 = {
    'bessel_i_data.txt': (id, iv),
    'bessel_i_int_data.txt': (id, iv),
    'bessel_j_data.txt': (id, jv),
    'bessel_j_int_data.txt': (id, jn),
    'bessel_j_large_data.txt': (id, jn),
    'bessel_k_data.txt': (id, kv),
    'bessel_k_int_data.txt': (id, kn),
    'bessel_y01_data.txt': (id, yn),
    'bessel_yn_data.txt': (id, yn),
    'bessel_yv_data.txt': (id, yv),
}

def scaled_error(x, ref, eps):
    aerr = np.abs(x-ref)
    rerr = aerr * np.abs(ref)
    return rerr / eps

def _test_factory_one_arg(funcs, datafile, dtype=np.double):
    data = np.fromfile(os.path.join(DATA_DIR, datafile), sep=",", dtype=dtype)
    data = data.reshape((data.size / len(funcs), len(funcs)))
    eps = np.finfo(dtype).eps
    for i in range(1, len(funcs)):
        assert_array_almost_equal(funcs[i](dtype(data[:,0])),
                dtype(data[:,i]))

def _test_factory_two_args(func, datafile, dtype=np.double):
    data = np.fromfile(os.path.join(DATA_DIR, datafile), sep=",", dtype=dtype)
    data = data.reshape((data.size / 3, 3))
    eps = np.finfo(dtype).eps
    #print func(dtype(data[:,0]), dtype(data[:,1])), dtype(data[:, 2])

    assert_array_almost_equal(func(dtype(data[:,0]), dtype(data[:,1])),
        dtype(data[:, 2]))

def test_one_arg():
    for k, v in TEST_MAPPING_1.items():
        func = _test_factory_one_arg
        test_one_arg.__doc__ = "Boost test for %s" %  \
            ",".join([f.__name__ for f in v[1:]])
        yield func, v, k

def test_two_args():
    for k, v in TEST_MAPPING_2.items():
        func = _test_factory_two_args
        test_one_arg.__doc__ = "Boost test for %s" %  v[1].__name__
        yield func, v[1], k
