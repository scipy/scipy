import numpy as np
import pyprima
from pyprima import minimize
from scipy.optimize import NonlinearConstraint
import sys
import debugpy


pyprima.common.linalg.USE_NAIVE_MATH = False

# debugpy.listen(5678)
# debugpy.wait_for_client()
np.set_printoptions(precision=19, floatmode='fixed', suppress=False)

def f(x):
    return np.sum(x**2)

# eps = np.finfo(float).eps
lb = [-1, None, 1, None, -0.5]
ub = [-0.5, -0.5, None, None, -0.5]
bounds = [(a, b) for a, b in zip(lb, ub)]

# lb_as_nlc = NonlinearConstraint(lambda x: x - lb, 0, np.inf)
# ub_as_nlc = NonlinearConstraint(lambda x: ub - x, 0, np.inf)

# res = minimize(f, x0=[1, 2, 3, 4, 5], method='cobyla', bounds=bounds, options={'rhobeg': 0.03997997997997998, 'iprint': 1})
# sys.exit()
# these are converted to Bounds internally

# debugpy.breakpoint()
failed = 0
n = 100
for i, rhobeg in enumerate(np.linspace(2e-2, 10, n)):
    res = minimize(f, x0=np.array([1, 2, 3, 4, 5]), bounds=bounds, options={'rhobeg': rhobeg})
    # print(res.nf)
    ref = [-0.5, -0.5, 1, 0, -0.5]
    ref = 1.75
    try:
        diff = (abs(res.fun - ref))
    except AttributeError:
        diff = (abs(res.f - ref))
    try:
        nf = res.nf
    except AttributeError:
        nf = res.nfev
    # print(f"diff={diff}")
    if diff > 1:
        print(f"rhobeg={rhobeg} FAILED with diff={diff} and nf={nf}")
        failed += 1
        # break
    if i % int(n/10) == 0:
        print(f"rhobeg={rhobeg} passed with diff={diff} and nf={nf}")

print("num failed: ", failed)