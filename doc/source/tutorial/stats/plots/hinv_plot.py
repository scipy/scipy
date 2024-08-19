import numpy as np
from scipy.stats.sampling import NumericalInverseHermite
from scipy.stats import norm
from scipy.special import ndtr
import matplotlib.pyplot as plt


class StandardNormal:
    def pdf(self, x):
        return 1/np.sqrt(2*np.pi) * np.exp(-x**2 / 2)

    def cdf(self, x):
        return ndtr(x)


dist = StandardNormal()
urng = np.random.default_rng()
rng = NumericalInverseHermite(dist, random_state=urng)
rvs = rng.rvs(10000)
x = np.linspace(rvs.min()-0.1, rvs.max()+0.1, 1000)
fx = norm.pdf(x)
plt.plot(x, fx, 'r-', lw=2, label='true distribution')
plt.hist(rvs, bins=20, density=True, alpha=0.8, label='random variates')
plt.xlabel('x')
plt.ylabel('PDF(x)')
plt.title('Numerical Inverse Hermite Samples')
plt.legend()
plt.show()
