import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import math


class Gamma:
    def __init__(self, p):
        self.p = p

    def pdf(self, x):
        if x < 0:
            return 0.0
        return x**(self.p - 1) * math.exp(-x)

    @staticmethod
    def support():
        return 0, np.inf


def u_bound(x, p, center):
    if x < 0:
        return 0
    return (x - center) * x**((p-1)/2) * math.exp(-x/2)


# bounding rectangle for Gamma(p) shifted by center
def rectangle(p, center):
    h = (p+1+center)/2
    k = np.sqrt(h**2 - center*(p-1))
    u_min, u_max = u_bound(h-k, p, center), u_bound(h+k, p, center)
    v_max = math.sqrt((p-1)**(p-1) * math.exp(-(p-1)))
    return u_min, u_max, v_max


urng = np.random.default_rng()
p = 2.2
dist = Gamma(p)
u_min, u_max, v_max = rectangle(p, p-1)
rng = stats.NaiveRatioUniforms(dist, center=p-1, v_max=v_max,
                               u_min=u_min, u_max=u_max, random_state=urng)
rvs = rng.rvs(1000)
x = np.linspace(rvs.min()-0.1, rvs.max()+0.1, num=500)
fx = stats.gamma.pdf(x, p)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, fx, "r-", label=" Gamma({}) pdf".format(p))
ax.hist(rvs, bins=30, density=True, alpha=0.8, label="rvs")
plt.xlabel("x")
plt.title("Samples drawn using NROU method with mode shift.")
plt.legend()
plt.show()
