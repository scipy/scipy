import numpy as np
from scipy import stats, special, integrate
import matplotlib.pyplot as plt
from mpmath import mp
from scipy.integrate.tests.test_quadrature import TestTanhSinh

def quadts(f, a, b, m=5):

    h = 1 / (2 ** m)
    j = np.arange(-5000, 5000)
    pi_2 = np.pi / 2

    jh = j * h
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)
    xj = np.tanh(u2)
    wj = u1 / np.cosh(u2)**2

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

    alpha = (b-a)/2
    beta = (a+b)/2

    xj = alpha*xj + beta
    wj *= alpha

    i = (xj > a) & (xj < b) & (wj > 1e-300) & (wj < 1e300)
    xj = xj[i]
    wj = wj[i]
    print(len(xj))

    fj = f(xj)
    s = (fj @ wj) * h
    return s

m = 10

# test = TestTanhSinh()
# f = test.f14
# ref = f.ref
# res = quadts(f, 0, f.b, m=m)

dist = stats.alpha(a=3.57)
a, b = -np.inf, 0.4
res = quadts(dist.pdf, 0.4, np.inf, m=m)
ref = dist.sf(0.4)

print(np.log10(np.abs((res-ref)/ref)))

