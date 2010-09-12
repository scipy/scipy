import os
import warnings

import numpy as np
from numpy.testing import assert_
from numpy.testing.noseclasses import KnownFailureTest

import scipy.special as sc

__all__ = ['with_special_errors', 'assert_tol_equal', 'assert_func_equal',
           'FuncData']

#------------------------------------------------------------------------------
# Enable convergence and loss of precision warnings -- turn off one by one
#------------------------------------------------------------------------------

def with_special_errors(func):
    """
    Enable special function errors (such as underflow, overflow,
    loss of precision, etc.)
    """
    def wrapper(*a, **kw):
        old_filters = list(getattr(warnings, 'filters', []))
        old_errprint = sc.errprint(1)
        warnings.filterwarnings("error", category=sc.SpecialFunctionWarning)
        try:
            return func(*a, **kw)
        finally:
            sc.errprint(old_errprint)
            setattr(warnings, 'filters', old_filters)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper

#------------------------------------------------------------------------------
# Comparing function values at many data points at once, with helpful
#------------------------------------------------------------------------------

def assert_tol_equal(a, b, rtol=1e-7, atol=0, err_msg='', verbose=True):
    """Assert that `a` and `b` are equal to tolerance ``atol + rtol*abs(b)``"""
    def compare(x, y):
        return np.allclose(x, y, rtol=rtol, atol=atol)
    a, b = np.asanyarray(a), np.asanyarray(b)
    header = 'Not equal to tolerance rtol=%g, atol=%g' % (rtol, atol)
    np.testing.utils.assert_array_compare(compare, a, b, err_msg=str(err_msg),
                                          verbose=verbose, header=header)

#------------------------------------------------------------------------------
# Comparing function values at many data points at once, with helpful
# error reports
#------------------------------------------------------------------------------

def assert_func_equal(func, results, points, rtol=None, atol=None,
                      param_filter=None, knownfailure=None,
                      vectorized=True, dtype=None):
    if hasattr(points, 'next'):
        # it's a generator
        points = list(points)

    points = np.asarray(points)
    if points.ndim == 1:
        points = points[:,None]

    if hasattr(results, '__name__'):
        # function
        if vectorized:
            results = results(*tuple(points.T))
        else:
            results = np.array([results(*tuple(p)) for p in points])
            if results.dtype == object:
                try:
                    results = results.astype(float)
                except TypeError:
                    results = results.astype(complex)
    else:
        results = np.asarray(results)

    npoints = points.shape[1]

    data = np.c_[points, results]
    fdata = FuncData(func, data, range(npoints), range(npoints, data.shape[1]),
                     rtol=rtol, atol=atol, param_filter=param_filter,
                     knownfailure=knownfailure)
    fdata.check()

class FuncData(object):
    """
    Data set for checking a special function.

    Parameters
    ----------
    func : function
        Function to test
    filename : str
        Input file name
    param_columns : int or tuple of ints
        Columns indices in which the parameters to `func` lie.
        Can be imaginary integers to indicate that the parameter
        should be cast to complex.
    result_columns : int or tuple of ints
        Column indices for expected results from `func`.
    rtol : float, optional
        Required relative tolerance. Default is 5*eps.
    atol : float, optional
        Required absolute tolerance. Default is 5*tiny.
    param_filter : function, or tuple of functions/Nones, optional
        Filter functions to exclude some parameter ranges.
        If omitted, no filtering is done.
    knownfailure : str, optional
        Known failure error message to raise when the test is run.
        If omitted, no exception is raised.

    """

    def __init__(self, func, data, param_columns, result_columns,
                 rtol=None, atol=None, param_filter=None, knownfailure=None,
                 dataname=None):
        self.func = func
        self.data = data
        self.dataname = dataname
        if not hasattr(param_columns, '__len__'):
            param_columns = (param_columns,)
        if not hasattr(result_columns, '__len__'):
            result_columns = (result_columns,)
        self.param_columns = tuple(param_columns)
        self.result_columns = tuple(result_columns)
        self.rtol = rtol
        self.atol = atol
        if not hasattr(param_filter, '__len__'):
            param_filter = (param_filter,)
        self.param_filter = param_filter
        self.knownfailure = knownfailure

    def get_tolerances(self, dtype):
        info = np.finfo(dtype)
        rtol, atol = self.rtol, self.atol
        if rtol is None:
            rtol = 5*info.eps
        if atol is None:
            atol = 5*info.tiny
        return rtol, atol

    def check(self, data=None, dtype=None):
        """Check the special function against the data."""

        if self.knownfailure:
            raise KnownFailureTest(self.knownfailure)

        if data is None:
            data = self.data

        if dtype is None:
            dtype = data.dtype
        else:
            data = data.astype(dtype)

        rtol, atol = self.get_tolerances(dtype)

        # Apply given filter functions
        if self.param_filter:
            param_mask = np.ones((data.shape[0],), np.bool_)
            for j, filter in zip(self.param_columns, self.param_filter):
                if filter:
                    param_mask &= filter(data[:,j])
            data = data[param_mask]

        # Pick parameters and results from the correct columns
        params = []
        for j in self.param_columns:
            if np.iscomplexobj(j):
                j = int(j.imag)
                params.append(data[:,j].astype(np.complex))
            else:
                params.append(data[:,j])
        wanted = tuple([data[:,j] for j in self.result_columns])

        # Evaluate
        got = self.func(*params)
        if not isinstance(got, tuple):
            got = (got,)

        # Check the validity of each output returned

        assert_(len(got) == len(wanted))

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
                # Some bad results: inform what, where, and how bad
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
                    d = fmt(rdiff)
                    msg.append("%s => %s != %s  (rdiff %s)" % (a, b, c, d))
                assert_(False, "\n".join(msg))

    def __repr__(self):
        """Pretty-printing, esp. for Nose output"""
        if np.any(map(np.iscomplexobj, self.param_columns)):
            is_complex = " (complex)"
        else:
            is_complex = ""
        if self.dataname:
            return "<Data for %s%s: %s>" % (self.func.__name__, is_complex,
                                            os.path.basename(self.dataname))
        else:
            return "<Data for %s%s>" % (self.func.__name__, is_complex)
