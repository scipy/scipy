import numpy as np
from scipy import stats, special, integrate
import matplotlib.pyplot as plt
from mpmath import mp
from scipy.integrate.tests.test_quadrature import TestTanhSinh

def compute_pairs(k):
    # Compute the abscissa-weight pairs for each level m. See [1] page 9.

    "....each level k of abscissa-weight pairs uses h = 2 **-k"
    h = 1 / (2 ** k)

    # "We find that roughly 3.6 * 2^k abscissa-weight pairs are generated at
    # "level k." The actual number per level can be generated like:
    # for i in range(10):
    #     _, xjc, wj = compute_pairs(i)
    #     # don't want want infinite weights or to evaluate f at endpoints
    #     valid = (xjc > 0) & (wj > 0) & np.isfinite(wj)
    #     print(np.sum(valid))
    # Running this code, I'm finding that the maximum index value w/ 64-bit is:
    max = int(np.ceil(6.115 * 2**k))
    # This reproduces all the integers produced by the loop above for k <= 10.
    # Note that the actual number of pairs is *half* of this (see below).

    # For iterations after the first, "....the integrand function needs to be
    # evaluated only at the odd-indexed abscissas at each level."
    j = np.arange(max) if k == 0 else np.arange(1, max, 2)
    jh = j * h

    # "In this case... the weights wj = u1/cosh(u2)^2, where..."
    pi_2 = np.pi / 2
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)
    wj = u1 / np.cosh(u2)**2

    # "We actually store 1-xj = 1/(...)."
    xjc = 1 / (np.exp(u2) * np.cosh(u2))  # complement of xj = np.tanh(u2)

    # When level k == 0, the zeroth xj corresponds with xj = 0. To simplify
    # code, the function will be evaluated there twice; each gets half weight.
    wj[0] = wj[0] / 2 if k == 0 else wj[0]
    return h, xjc, wj

def quadts(f, a, b, maxiter=10):

    if np.isinf(a) and np.isinf(b):
        def f(x, f=f):
            return f(x) + f(-x)
        a, b = 0, np.inf

    if np.isinf(a):
        def f(x, f=f):
            return f(-x)
        a, b = -b, -a

    if np.isinf(b):
        def f(x, f=f, a=a):
            return f(1/x - 1 + a)*x**-2
        a, b = 0, 1

    Sk = []

    for n in range(maxiter):
        h, xjc, wj = compute_pairs(n)

        # If we had stored xj instead of xjc, we would have
        # xj = alpha * xj + beta, where beta = (a + b)/2
        alpha = (b-a)/2
        xj = np.concatenate((-alpha * xjc + b, alpha * xjc + a))
        wj *= alpha
        wj = np.concatenate((wj, wj))

        # Temporarily comment this out so that it's easier to figure out which
        # abscissae are closest to left and right endpoints for error estimate.
        # # This is needed for functions that behave badly at the endpoints
        # i = (xj > a) & (xj < b) & (wj > 1e-100) & (wj < 1e100)
        # xj = xj[i]
        # wj = wj[i]

        # todo:
        #  store fj * wj
        #  add back endpoint protection
        #  keep track of function evaluations
        #  absolute and relative tolerance options
        #  vectorize
        fj = f(xj)
        Snm1 = 0 if not Sk else Sk[-1]
        Sn = Snm1/2 + (fj @ wj) * h  # update integral estimate

        # Check error estimate (see "5. Error Estimation, page 11")
        if n >= 2:
            Snm2, Snm1 = Sk[-2:]
            d1 = np.log10(abs(Sn - Snm1))
            d2 = np.log10(abs(Sn - Snm2))
            e1 = np.finfo(np.float64).eps
            d3 = np.log10(e1 * np.max(np.abs(wj*fj)))
            d4 = np.log10(np.max(np.reshape(np.abs(wj*fj), (2, -1))[:, -1]))
            d = np.max([d1**2/d2, 2*d1, d3, d4])
            if d < -15:
                break

        Sk.append(Sn)

    return Sn

test = TestTanhSinh()
f = test.f1
ref = f.ref
res = quadts(f, 0, f.b)

# dist = stats.alpha(a=3.57)
# a, b = -np.inf, 0.4
# res = quadts(dist.pdf, 0.4, np.inf)
# ref = dist.sf(0.4)

print(np.log10(np.abs((res-ref)/ref)))

