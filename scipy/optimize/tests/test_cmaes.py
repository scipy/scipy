import numpy as np
#from scipy.optimize_optimize import _minimize_cmaes, rosen
from scipy.optimize._optimize import _minimize_cmaes, rosen
dim = 40
mean = 3 * np.ones(40)
sigma = 2

_minimize_cmaes(func=rosen, mean=mean, sigma=sigma)