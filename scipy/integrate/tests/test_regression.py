import numpy as np

from scipy.integrate import odeint


def test_issue_8217():
    t = np.linspace(0,0,10)

    def func(x):
        return x**2

    sol = odeint(func, [1.], t)
