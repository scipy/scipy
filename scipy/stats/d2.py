from scipy.stats import norm
from scipy.integrate import quad
import numpy as np


def d2(n=2, method="exact"):
    """
    Computes the d2 statistic used for measuring
    relative range standard deviation approximation.

    The d2 statistic is often found in statistical process
    control tables.

    That is, standard deviation can be approximated by
    R / d2 where R is the range of the data (max - min)
    and d2 is computed as below.

    See also:
    https://v8doc.sas.com/sashtml/qc/chapc/sect9.htm

    :param n:
        Number of distributions considered for computing
        d2. All are normal with mean=0 and std=1.
    :param method:
        Optional. Method to use. Default "exact".
            - "exact" computes the exact infinite integral.
            - "random" uses normal distributions generated
              on the fly to compute d2.
    :return:
        d2:
            the expectation value [that is, average here]
            of the ranges of the distributions
    """
    if method == "exact":
        def f(x, n):
            return 1 - (1 - norm.cdf(x)) ** n - (norm.cdf(x)) ** n

        d2 = quad(f, -np.inf, np.inf, args=(n))[0]
        return d2

    elif method == "random":
        x = {}
        # slots to fill with normally distributed samples
        for i in range(n):
            x[i] = norm.rvs(size=100000, loc=0, scale=1)

        x = np.vstack([x[i] for i in x])

        maxs = np.amax(x, axis=0)
        mins = np.amin(x, axis=0)

        r = maxs - mins

        d2 = np.average(r)
        return d2

