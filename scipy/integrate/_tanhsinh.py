import numpy as np
from scipy import stats, special, integrate
import matplotlib.pyplot as plt
from mpmath import mp
from scipy.integrate.tests.test_quadrature import TestTanhSinh

def compute_pairs(k):
    # Compute the abscissa-weight pairs for each level m
    h = 1 / (2 ** k)

    # empirically, this gives the maximum number of pairs with valid abscissae
    # (magnitude strictly less than 1) and finite weight. The list can be
    # generated like:
    # for i in range(10):
    #     _, xjc, wj = compute_pairs(i)
    #     valid = (xjc > 0) & (wj > 0) & np.isfinite(wj)
    #     print(np.sum(valid))
    # 6.115 happens to be the magic number for 64-bit and k <= 10.
    # I think this is comparable to "3.6*2**k" in the paper. If we were to
    # generate more, they would cause the function to be evaluated at the
    # endpoints of the interval or be given infinite weight. Neither of these
    # is supposed to happen.
    max = int(np.ceil(6.115 * 2**k))
    j = np.arange(max)
    jh = j * h

    pi_2 = np.pi / 2
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)

    # See [1] page 9. "We actually store 1-xj = 1/(...)."
    xjc = 1 / (np.exp(u2) * np.cosh(u2))  # complement of xj = np.tanh(u2)
    wj = u1 / np.cosh(u2)**2

    return h, xjc, wj

def quadts(f, a, b, m=5):


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

    s = 0

    for i in range(m+1):
        s /= 2
        h, xjc, wj = compute_pairs(m)

        if i == 0:
            wj[0] /= 2
        else:
            wj = wj[1::2]
            xjc = xjc[1::2]

        wj = np.concatenate((wj, wj))

        alpha = (b-a)/2
        xj = np.concatenate((-alpha * xjc + b,
                             alpha * xjc + a))
        wj *= alpha

        i = (xj > a) & (xj < b) & (wj > 1e-100) & (wj < 1e100)
        xj = xj[i]
        wj = wj[i]

        fj = f(xj)
        s += (fj @ wj) * h

    return s

m = 4

test = TestTanhSinh()
f = test.f1
ref = f.ref
res = quadts(f, 0, f.b, m=m)

# dist = stats.alpha(a=3.57)
# a, b = -np.inf, 0.4
# res = quadts(dist.pdf, 0.4, np.inf, m=m)
# ref = dist.sf(0.4)

print(np.log10(np.abs((res-ref)/ref)))

