import numpy as np
from scipy.optimize._optimize import _minimize_cmaes, rosen

mean = np.ones(5)
sigma = 1

_minimize_cmaes(func=rosen, mean=mean, sigma=sigma)