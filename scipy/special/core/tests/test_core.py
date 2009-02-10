import os

import numpy as np
from numpy.testing import *
from scipy.special import arccosh, arcsinh, arctanh, erf, erfc, log1p, expm1

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

def test_one_arg():
    for k, v in TEST_MAPPING_1.items():
        func = _test_factory_one_arg
        test_one_arg.__doc__ = "Boost test for %s" %  \
            ",".join([f.__name__ for f in v[1:]])
        yield func, v, k

