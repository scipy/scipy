# Pink noise generator

from scipy.base import *
from scipy import stats

def pink1d(n, rvs=stats.norm.rvs):
    k = min(int(floor(log(n)/log(2))), 6)
    pink = zeros((n,), Float)
    m = 1
    for i in range(k):
        p = int(ceil(float(n) / m))
        pink += repeat(rvs(size=p), m)[:n]
        m <<= 1

    return pink/k
