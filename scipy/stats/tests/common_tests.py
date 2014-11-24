from __future__ import division, print_function, absolute_import

import inspect
import warnings

import numpy as np
import numpy.testing as npt

from scipy.lib._version import NumpyVersion
from scipy import stats


NUMPY_BELOW_1_7 = NumpyVersion(np.__version__) < '1.7.0'


def check_normalization(distfn, args, distname):
    norm_moment = distfn.moment(0, *args)
    npt.assert_allclose(norm_moment, 1.0)

    # this is a temporary plug: either ncf or expect is problematic;
    # best be marked as a knownfail, but I've no clue how to do it.
    if distname == "ncf":
        atol, rtol = 1e-5, 0
    else:
        atol, rtol = 1e-7, 1e-7

    normalization_expect = distfn.expect(lambda x: 1, args=args)
    npt.assert_allclose(normalization_expect, 1.0, atol=atol, rtol=rtol,
            err_msg=distname, verbose=True)

    normalization_cdf = distfn.cdf(distfn.b, *args)
    npt.assert_allclose(normalization_cdf, 1.0)


def check_moment(distfn, arg, m, v, msg):
    m1 = distfn.moment(1, *arg)
    m2 = distfn.moment(2, *arg)
    if not np.isinf(m):
        npt.assert_almost_equal(m1, m, decimal=10, err_msg=msg +
                            ' - 1st moment')
    else:                     # or np.isnan(m1),
        npt.assert_(np.isinf(m1),
               msg + ' - 1st moment -infinite, m1=%s' % str(m1))

    if not np.isinf(v):
        npt.assert_almost_equal(m2 - m1 * m1, v, decimal=10, err_msg=msg +
                            ' - 2ndt moment')
    else:                     # or np.isnan(m2),
        npt.assert_(np.isinf(m2),
               msg + ' - 2nd moment -infinite, m2=%s' % str(m2))


def check_mean_expect(distfn, arg, m, msg):
    if np.isfinite(m):
        m1 = distfn.expect(lambda x: x, arg)
        npt.assert_almost_equal(m1, m, decimal=5, err_msg=msg +
                            ' - 1st moment (expect)')


def check_var_expect(distfn, arg, m, v, msg):
    if np.isfinite(v):
        m2 = distfn.expect(lambda x: x*x, arg)
        npt.assert_almost_equal(m2, v + m*m, decimal=5, err_msg=msg +
                            ' - 2st moment (expect)')


def check_skew_expect(distfn, arg, m, v, s, msg):
    if np.isfinite(s):
        m3e = distfn.expect(lambda x: np.power(x-m, 3), arg)
        npt.assert_almost_equal(m3e, s * np.power(v, 1.5),
                decimal=5, err_msg=msg + ' - skew')
    else:
        npt.assert_(np.isnan(s))


def check_kurt_expect(distfn, arg, m, v, k, msg):
    if np.isfinite(k):
        m4e = distfn.expect(lambda x: np.power(x-m, 4), arg)
        npt.assert_allclose(m4e, (k + 3.) * np.power(v, 2), atol=1e-5, rtol=1e-5,
                err_msg=msg + ' - kurtosis')
    else:
        npt.assert_(np.isnan(k))


def check_entropy(distfn, arg, msg):
    ent = distfn.entropy(*arg)
    npt.assert_(not np.isnan(ent), msg + 'test Entropy is nan')


def check_private_entropy(distfn, args, superclass):
    # compare a generic _entropy with the distribution-specific implementation
    npt.assert_allclose(distfn._entropy(*args),
                        superclass._entropy(distfn, *args))


def check_edge_support(distfn, args):
    # Make sure the x=self.a and self.b are handled correctly.
    x = [distfn.a, distfn.b]
    if isinstance(distfn, stats.rv_continuous):
        npt.assert_equal(distfn.cdf(x, *args), [0.0, 1.0])
        npt.assert_equal(distfn.logcdf(x, *args), [-np.inf, 0.0])

        npt.assert_equal(distfn.sf(x, *args), [1.0, 0.0])
        npt.assert_equal(distfn.logsf(x, *args), [0.0, -np.inf])

    if isinstance(distfn, stats.rv_discrete):
        x = [distfn.a - 1, distfn.b]
    npt.assert_equal(distfn.ppf([0.0, 1.0], *args), x)
    npt.assert_equal(distfn.isf([0.0, 1.0], *args), x[::-1])

    # out-of-bounds for isf & ppf
    npt.assert_(np.isnan(distfn.isf([-1, 2], *args)).all())
    npt.assert_(np.isnan(distfn.ppf([-1, 2], *args)).all())


def check_named_args(distfn, x, shape_args, defaults, meths):
    ## Check calling w/ named arguments.

    # check consistency of shapes, numargs and _parse signature
    signature = inspect.getargspec(distfn._parse_args)
    npt.assert_(signature.varargs is None)
    npt.assert_(signature.keywords is None)
    npt.assert_(signature.defaults == defaults)

    shape_argnames = signature.args[1:-len(defaults)]  # self, a, b, loc=0, scale=1
    if distfn.shapes:
        shapes_ = distfn.shapes.replace(',', ' ').split()
    else:
        shapes_ = ''
    npt.assert_(len(shapes_) == distfn.numargs)
    npt.assert_(len(shapes_) == len(shape_argnames))

    # check calling w/ named arguments
    shape_args = list(shape_args)

    vals = [meth(x, *shape_args) for meth in meths]
    npt.assert_(np.all(np.isfinite(vals)))

    names, a, k = shape_argnames[:], shape_args[:], {}
    while names:
        k.update({names.pop(): a.pop()})
        v = [meth(x, *a, **k) for meth in meths]
        npt.assert_array_equal(vals, v)
        if 'n' not in k.keys():
            # `n` is first parameter of moment(), so can't be used as named arg
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                npt.assert_equal(distfn.moment(1, *a, **k),
                                 distfn.moment(1, *shape_args))

    # unknown arguments should not go through:
    k.update({'kaboom': 42})
    npt.assert_raises(TypeError, distfn.cdf, x, **k)
