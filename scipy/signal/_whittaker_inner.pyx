import numpy as np


def _solve_WH_order2(double[:] y, double lamb):
    """Solve Whittaker-Henderson smoothing of order 2 according to Weinert."""
    n = y.shape[0]
    lamb = 1.0 / lamb

    b = np.empty(n)
    e = np.empty(n)
    f = np.empty(n)
    x = np.empty(n)

    d = 1 + lamb
    f[0] = 1 / d
    mu = 2
    e[0] = mu / d
    b[0] = lamb * y[0] / d
    mu_old = mu

    if n == 3:
        d = 4 + lamb - mu_old * e[0]
        mu = 2 - e[0]
    else:
        d = 5 + lamb - mu_old * e[0]
        mu = 4 - e[0]
    f[1] = 1 / d
    e[1] = mu / d
    b[1] = (lamb * y[1] + mu_old * b[0]) / d
    mu_old = mu

    for i in range(2, n - 2):
        d = 6 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 4 - e[i-1]
        e[i] = mu / d
        b[i] = (lamb * y[i] + mu_old * b[i-1] - b[i-2]) / d
        mu_old = mu

    if n >= 4:
        i = n - 2
        d = 5 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 2 - e[i-1]
        e[i] = mu / d
        b[i] = (lamb * y[i] + mu_old * b[i-1] - b[i-2]) / d
        mu_old = mu

    i = n - 1
    d = 1 + lamb - mu_old * e[i-1] - f[i-2]
    f[i] = 1 / d
    b[i] = (lamb * y[i] + mu_old * b[i-1] - b[i-2]) / d

    x[n-1] = b[n-1]
    x[n-2] = b[n-2] + e[n-2] * x[n-1]
    for i in range(n - 3, -1, -1):
        x[i] = b[i] + e[i] * x[i+1] - f[i] * x[i+2]
    return x
