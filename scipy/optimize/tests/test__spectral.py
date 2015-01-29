from __future__ import division, absolute_import, print_function

import itertools
import functools

import numpy as np
from numpy import exp
from numpy.testing import assert_, assert_equal, assert_allclose

from scipy.optimize import root, fixpoint, minimize
from scipy.stats import poisson
from scipy.special import expit, logit


def test_performance():
    # Compare performance results to those listed in
    # [Cheng & Li, IMA J. Num. An. 29, 814 (2008)]
    # and
    # [W. La Cruz, J.M. Martinez, M. Raydan, Math. Comp. 75, 1429 (2006)].
    # and those produced by dfsane.f from M. Raydan's website.
    #
    # Where the results disagree, the largest limits are taken.

    e_a = 1e-5
    e_r = 1e-4

    table_1 = [
        dict(F=F_1, x0=x0_1, n=1000, nit=5, nfev=5),
        dict(F=F_1, x0=x0_1, n=10000, nit=2, nfev=2),
        dict(F=F_2, x0=x0_2, n=500, nit=11, nfev=11),
        dict(F=F_2, x0=x0_2, n=2000, nit=11, nfev=11),
        # dict(F=F_4, x0=x0_4, n=999, nit=243, nfev=1188),  removed: too sensitive to rounding errors
        dict(F=F_6, x0=x0_6, n=100, nit=6, nfev=6),  # Results from dfsane.f; papers list nit=3, nfev=3
        dict(F=F_7, x0=x0_7, n=99, nit=23, nfev=29),  # Must have n%3==0, typo in papers?
        dict(F=F_7, x0=x0_7, n=999, nit=23, nfev=29),  # Must have n%3==0, typo in papers?
        dict(F=F_9, x0=x0_9, n=100, nit=12, nfev=18),  # Results from dfsane.f; papers list nit=nfev=6?
        dict(F=F_9, x0=x0_9, n=1000, nit=12, nfev=18),
        dict(F=F_10, x0=x0_10, n=1000, nit=5, nfev=5),  # Results from dfsane.f; papers list nit=2, nfev=12
    ]

    # Check also scaling invariance
    for xscale, yscale, line_search in itertools.product([1.0, 1e-10, 1e10], [1.0, 1e-10, 1e10],
                                                         ['cruz', 'cheng']):
        for problem in table_1:
            n = problem['n']
            func = lambda x, n: yscale*problem['F'](x/xscale, n)
            args = (n,)
            x0 = problem['x0'](n) * xscale

            fatol = np.sqrt(n) * e_a * yscale + e_r * np.linalg.norm(func(x0, n))

            sigma_eps = 1e-10 * min(yscale/xscale, xscale/yscale)
            sigma_0 = xscale/yscale

            with np.errstate(over='ignore'):
                sol = root(func, x0, args=args,
                           options=dict(ftol=0, fatol=fatol, maxfev=problem['nfev'] + 1,
                                        sigma_0=sigma_0, sigma_eps=sigma_eps,
                                        line_search=line_search),
                           method='DF-SANE')

            err_msg = repr([xscale, yscale, line_search, problem, np.linalg.norm(func(sol.x, n)),
                            fatol, sol.success, sol.nit, sol.nfev])
            assert_(sol.success, err_msg)
            assert_(sol.nfev <= problem['nfev'] + 1, err_msg)  # nfev+1: dfsane.f doesn't count first eval
            assert_(sol.nit <= problem['nit'], err_msg)
            assert_(np.linalg.norm(func(sol.x, n)) <= fatol, err_msg)


def test_complex():
    def func(z):
        return z**2 - 1 + 2j
    x0 = 2.0j

    ftol = 1e-4
    sol = root(func, x0, tol=ftol, method='DF-SANE')

    assert_(sol.success)

    f0 = np.linalg.norm(func(x0))
    fx = np.linalg.norm(func(sol.x))
    assert_(fx <= ftol*f0)


def test_linear_definite():
    # The DF-SANE paper proves convergence for "strongly isolated"
    # solutions.
    #
    # For linear systems F(x) = A x - b = 0, with A positive or
    # negative definite, the solution is strongly isolated.

    def check_solvability(A, b, line_search='cruz'):
        func = lambda x: A.dot(x) - b
        xp = np.linalg.solve(A, b)
        eps = np.linalg.norm(func(xp)) * 1e3
        sol = root(func, b, options=dict(fatol=eps, ftol=0, maxfev=17523, line_search=line_search),
                   method='DF-SANE')
        assert_(sol.success)
        assert_(np.linalg.norm(func(sol.x)) <= eps)

    n = 90

    # Test linear pos.def. system
    np.random.seed(1234)
    A = np.arange(n*n).reshape(n, n)
    A = A + n*n * np.diag(1 + np.arange(n))
    assert_(np.linalg.eigvals(A).min() > 0)
    b = np.arange(n) * 1.0
    check_solvability(A, b, 'cruz')
    check_solvability(A, b, 'cheng')

    # Test linear neg.def. system
    check_solvability(-A, b, 'cruz')
    check_solvability(-A, b, 'cheng')


def test_shape():
    def f(x, arg):
        return x - arg

    for dt in [float, complex]:
        x = np.zeros([2,2])
        arg = np.ones([2,2], dtype=dt)

        sol = root(f, x, args=(arg,), method='DF-SANE')
        assert_(sol.success)
        assert_equal(sol.x.shape, x.shape)


def test_stagnation():
    def f(x):
        return 1 + np.random.rand()

    np.random.seed(1234)
    sol = root(f, [0], options=dict(fatol=1e-6, ftol=0, maxfev=10000, stagnation_limit=20),
               method='df-sane')
    assert_(not sol.success)
    assert_equal(sol.message, "convergence stagnated")
    assert_(sol.nfev < 1000, repr(sol.nfev))


def test_squarem():
    # fixed point should be sqrt(c)
    def func(x, c):
        return 0.5*(x + c/x)

    c = np.array([0.75, 1.0, 1.25])
    x0 = [1.1, 1.15, 0.9]
    olderr = np.seterr(all='ignore')
    try:
        sol = fixpoint(func, x0, args=(c,), method='squarem')
    finally:
        np.seterr(**olderr)
    assert_(sol.success)
    assert_allclose(sol.x, np.sqrt(c))


# Some of the test functions and initial guesses listed in
# [W. La Cruz, M. Raydan. Optimization Methods and Software, 18, 583 (2003)]

def F_1(x, n):
    g = np.zeros([n])
    i = np.arange(2, n+1)
    g[0] = exp(x[0] - 1) - 1
    g[1:] = i*(exp(x[1:] - 1) - x[1:])
    return g

def x0_1(n):
    x0 = np.empty([n])
    x0.fill(n/(n-1))
    return x0

def F_2(x, n):
    g = np.zeros([n])
    i = np.arange(2, n+1)
    g[0] = exp(x[0]) - 1
    g[1:] = 0.1*i*(exp(x[1:]) + x[:-1] - 1)
    return g

def x0_2(n):
    x0 = np.empty([n])
    x0.fill(1/n**2)
    return x0

def F_4(x, n):
    assert_equal(n % 3, 0)
    g = np.zeros([n])
    # Note: the first line is typoed in some of the references;
    # correct in original [Gasparo, Optimization Meth. 13, 79 (2000)]
    g[::3] = 0.6 * x[::3] + 1.6 * x[1::3]**3 - 7.2 * x[1::3]**2 + 9.6 * x[1::3] - 4.8
    g[1::3] = 0.48 * x[::3] - 0.72 * x[1::3]**3 + 3.24 * x[1::3]**2 - 4.32 * x[1::3] - x[2::3] + 0.2 * x[2::3]**3 + 2.16
    g[2::3] = 1.25 * x[2::3] - 0.25*x[2::3]**3
    return g

def x0_4(n):
    assert_equal(n % 3, 0)
    x0 = np.array([-1, 1/2, -1] * (n//3))
    return x0

def F_6(x, n):
    c = 0.9
    mu = (np.arange(1, n+1) - 0.5)/n
    return x - 1/(1 - c/(2*n) * (mu[:,None]*x / (mu[:,None] + mu)).sum(axis=1))

def x0_6(n):
    return np.ones([n])

def F_7(x, n):
    assert_equal(n % 3, 0)

    def phi(t):
        v = 0.5*t - 2
        v[t > -1] = ((-592*t**3 + 888*t**2 + 4551*t - 1924)/1998)[t > -1]
        v[t >= 2] = (0.5*t + 2)[t >= 2]
        return v
    g = np.zeros([n])
    g[::3] = 1e4 * x[1::3]**2 - 1
    g[1::3] = exp(-x[::3]) + exp(-x[1::3]) - 1.0001
    g[2::3] = phi(x[2::3])
    return g

def x0_7(n):
    assert_equal(n % 3, 0)
    return np.array([1e-3, 18, 1] * (n//3))

def F_9(x, n):
    g = np.zeros([n])
    i = np.arange(2, n)
    g[0] = x[0]**3/3 + x[1]**2/2
    g[1:-1] = -x[1:-1]**2/2 + i*x[1:-1]**3/3 + x[2:]**2/2
    g[-1] = -x[-1]**2/2 + n*x[-1]**3/3
    return g

def x0_9(n):
    return np.ones([n])

def F_10(x, n):
    return np.log(1 + x) - x/n

def x0_10(n):
    return np.ones([n])


class _PoissonMixtureUnboundedParameterization(object):

    @classmethod
    def encode(cls, p, mu):
        return np.concatenate((cls._encode_distribution(p), np.log(mu)))

    @classmethod
    def decode(cls, X):
        X = np.asarray(X)
        n = X.shape[0]
        encoded_p, log_mu = X[:n//2], X[n//2:]
        p = cls._decode_distribution(encoded_p)
        mu = np.exp(log_mu)
        return p, mu

    @classmethod
    def _encode_distribution(cls, X):
        Y = []
        t = 0
        for x in X[:-1]:
            Y.append(x / (1 - t))
            t += x
        return logit(Y)

    @classmethod
    def _decode_distribution(cls, Y):
        X = []
        t = 0
        for y in expit(Y):
            x = y * (1 - t)
            X.append(x)
            t += x
        X.append(1 - t)
        return np.array(X)


class _PoissonMixtureSimpleParameterization(object):

    @classmethod
    def encode(cls, p, mu):
        return np.concatenate((p[:-1], mu))

    @classmethod
    def decode(cls, X):
        X = np.asarray(X)
        n = X.shape[0]
        truncated_p, mu = X[:n//2], X[n//2:]
        p = np.concatenate((truncated_p, [1-truncated_p.sum()]))
        return p, mu


def _poisson_mixture_neg_log_likelihood(codec, data, X):
    p, mu = codec.decode(X)
    n = data.shape[0]
    counts = np.arange(n)
    likelihoods = poisson.pmf(counts[:, np.newaxis], mu).dot(p)
    return -data.dot(np.log(likelihoods))


def _poisson_mixture_em_step(codec, data, X):
    # For testing EM acceleration with the squarem algorithm.
    p, mu = codec.decode(X)
    n = data.shape[0]
    counts = np.arange(n)
    unnormalized_pi = p * poisson.pmf(counts[:, np.newaxis], mu)
    pi = unnormalized_pi / unnormalized_pi.sum(axis=1)[:, np.newaxis]
    p = data.dot(pi) / data.sum()
    mu = (data * counts).dot(pi) / data.dot(pi)
    return codec.encode(p, mu)


def test_squarem_poisson_mixtures_table_2():
    # Ravi Varadhan and Christophe Roland.
    # Simple and Globally Convergent Methods for Accelerating
    # the Convergence of Any EM Algorithm.
    # Scandinavian Journal of Statistics, Vol 35: 335--353, 2008.
    #
    # This is a poisson mixture model, with three parameters:
    # a mixture proportion and two poisson rates.
    #
    # For the given data, the maximum likelihood estimates should be:
    # mixture proportion of the first poisson rate : 0.3599
    # the first poisson rate : 1.256
    # the second poisson rate : 2.663
    # Note that this is not distinguishable from the estimates
    # (1.0 - 0.3599, 2.663, 1.256).
    #
    # The EM strategy is to first conditionally distribute the blame
    # for each of the ten counts between the two Poisson processes
    # (this is like the expectation step),
    # and then to use these condtional distributions to re-estimate the
    # overall mixing proportion and the two Poisson rates.
    #
    data = np.array([162, 267, 271, 185, 111, 61, 27, 8, 3, 1], dtype=float)
    p = np.array([0.4, 0.6])
    mu = np.array([1.0, 2.0])
    p_desired = [0.3599, 0.6401]
    mu_desired = [1.256, 2.663]

    # Test the unaccelerated EM.
    codec = _PoissonMixtureSimpleParameterization
    X0 = codec.encode(p, mu)
    iterations = 0
    prev = X0
    xtol = 1e-8
    while True:
        X = _poisson_mixture_em_step(codec, data, prev)
        iterations += 1
        if np.linalg.norm(X - prev) < xtol * np.linalg.norm(prev):
            break
        prev = X
    p_opt, mu_opt = codec.decode(X)
    assert_allclose(p_opt, p_desired, atol=1e-3)
    assert_allclose(mu_opt, mu_desired, atol=1e-3)
    # Check the ballpark number of EM iterations.
    assert_(2061-5 < iterations < 2061+5)

    # Test the accelerated EM.
    codec = _PoissonMixtureSimpleParameterization
    X0 = codec.encode(p, mu)
    f = functools.partial(_poisson_mixture_em_step, codec, data)
    sol = fixpoint(f, X0, method='squarem')
    p_opt, mu_opt = codec.decode(sol.fun)
    assert_allclose(p_opt, p_desired, atol=1e-3)
    assert_allclose(mu_opt, mu_desired, atol=1e-3)
    # Check the ballpark number of SQUAREM iterations and function evaluations.
    assert_(28-5 < sol.nit < 28+5)
    assert_(108-5 < sol.nfev < 108+5)

    # Test a direct minimization of negative log likelihood.
    codec = _PoissonMixtureUnboundedParameterization
    X0 = codec.encode(p, mu)
    obj = functools.partial(_poisson_mixture_neg_log_likelihood, codec, data)
    result = minimize(obj, X0, method='L-BFGS-B', options=dict(factr=10.0))
    p_opt, mu_opt = codec.decode(result.x)
    assert_allclose(p_opt, p_desired, atol=1e-3)
    assert_allclose(mu_opt, mu_desired, atol=1e-3)
