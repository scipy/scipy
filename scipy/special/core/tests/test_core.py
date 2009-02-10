import os

import numpy as np
from numpy.testing import *
from scipy.special import (
    arccosh, arcsinh, arctanh, erf, erfc, log1p, expm1, 
    jn, jv, yn, yv, iv, kv, kn, gamma, gammaln,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

def test_all():

    TESTS = [
        Data(arccosh, 'acosh_data.txt', 0, 1),

        Data(arccosh, 'acosh_data.txt', 0, 1),
        Data(arcsinh, 'asinh_data.txt', 0, 1),
        Data(arctanh, 'atanh_data.txt', 0, 1),
        Data(gamma, 'gamma_data.txt', 0, 1),
        Data(gammaln, 'gamma_data.txt', 0, 2, rtol=5e-11),

        #assoc_legendre_p.txt
        #Data(erf, 'erf_inv_data.txt', 0, 1),

        Data(erf, 'erf_data.txt', 0, 1),
        Data(erfc, 'erf_data.txt', 0, 2),
        Data(erf, 'erf_large_data.txt', 0, 1),
        Data(erfc, 'erf_large_data.txt', 0, 2),
        Data(erf, 'erf_small_data.txt', 0, 1),
        Data(erfc, 'erf_small_data.txt', 0, 2),

        Data(log1p, 'log1p_expm1_data.txt', 0, 1),
        Data(expm1, 'log1p_expm1_data.txt', 0, 2),

        Data(iv, 'bessel_i_data.txt', (0,1), 2, rtol=1e-12),
        Data(iv, 'bessel_i_int_data.txt', (0,1), 2, rtol=1e-12),

        Data(jn, 'bessel_j_int_data.txt', (0,1), 2, rtol=1e-12),
        Data(jn, 'bessel_j_large_data.txt', (0,1), 2, rtol=6e-11),
        Data(jv, 'bessel_j_data.txt', (0,1), 2, rtol=1e-12),

        Data(kn, 'bessel_k_int_data.txt', (0,1), 2, rtol=1e-12),
        Data(kv, 'bessel_k_data.txt', (0,1), 2, rtol=1e-12),

        Data(yn, 'bessel_y01_data.txt', (0,1), 2, rtol=1e-12),
        Data(yn, 'bessel_yn_data.txt', (0,1), 2, rtol=1e-12),
        Data(yv, 'bessel_yv_data.txt', (0,1), 2, rtol=1e-12),
    ]

    for test in TESTS:
        yield _test_factory, test

def _test_factory(test, dtype=np.double):
    """Boost test"""
    test.check(dtype=dtype)

class Data(object):
    def __init__(self, func, filename, param_columns, result_columns,
                 rtol=None, atol=None):
        self.func = func
        self.filename = os.path.join(DATA_DIR, filename)
        if not hasattr(param_columns, '__len__'):
            param_columns = (param_columns,)
        if not hasattr(result_columns, '__len__'):
            result_columns = (result_columns,)
        self.param_columns = tuple(param_columns)
        self.result_columns = tuple(result_columns)
        self.rtol = rtol
        self.atol = atol

    def get_tolerances(self, dtype):
        info = np.finfo(dtype)
        rtol, atol = self.rtol, self.atol
        if rtol is None:
            rtol = 5*info.eps
        if atol is None:
            atol = 5*info.tiny
        return rtol, atol

    @staticmethod
    def load_data(filename, dtype):
        f = open(filename, 'r')
        try:
            ncols = 1
            while True:
                first_line = f.readline().strip()
                if first_line:
                    ncols = len(first_line.split())
                    break
        finally:
            f.close()
        data = np.fromfile(filename, sep=' ', dtype=dtype)
        nrows = data.shape[0] / ncols
        data = data.reshape((nrows, ncols))
        return data

    def check(self, data=None, dtype=np.double):
        if data is None:
            data = Data.load_data(self.filename, dtype)

        rtol, atol = self.get_tolerances(dtype)

        params = tuple([data[:,j] for j in self.param_columns])
        wanted = tuple([data[:,j] for j in self.result_columns])
        got = self.func(*params)

        if not isinstance(got, tuple):
            got = (got,)

        assert len(got) == len(wanted)

        for output_num, (x, y) in enumerate(zip(got, wanted)):
            pinf_x = np.isinf(x) & (x > 0)
            pinf_y = np.isinf(y) & (x > 0)
            minf_x = np.isinf(x) & (x < 0)
            minf_y = np.isinf(y) & (x < 0)
            nan_x = np.isnan(x)
            nan_y = np.isnan(y)

            abs_y = np.absolute(y)
            abs_y[~np.isfinite(abs_y)] = 0
            diff = np.absolute(x - y)
            diff[~np.isfinite(diff)] = 0

            rdiff = diff / np.absolute(y)
            rdiff[~np.isfinite(rdiff)] = 0

            tol_mask = (diff < atol + rtol*abs_y)
            pinf_mask = (pinf_x == pinf_y)
            minf_mask = (minf_x == minf_y)
            nan_mask = (nan_x == nan_y)

            bad_j = ~(tol_mask & pinf_mask & minf_mask & nan_mask)

            if np.any(bad_j):
                msg = [""]
                msg.append("Max |adiff|: %g" % diff.max())
                msg.append("Max |rdiff|: %g" % rdiff.max())
                msg.append("Bad results for the following points (in output %d):"
                           % output_num)
                for j in np.where(bad_j)[0]:
                    j = int(j)
                    fmt = lambda x: "%30s" % np.array2string(x[j], precision=18)
                    a = "  ".join(map(fmt, params))
                    b = "  ".join(map(fmt, got))
                    c = "  ".join(map(fmt, wanted))
                    msg.append("%s => %s != %s" % (a, b, c))
                assert False, "\n".join(msg)

    def __repr__(self):
        return "<Boost test for %s: %s>" % (self.func.__name__,
                                            os.path.basename(self.filename))
