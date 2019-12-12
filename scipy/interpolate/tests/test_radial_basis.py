"""Test the new RBF function"""

from __future__ import division, print_function, absolute_import


import numpy as np
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_almost_equal)

import itertools

from scipy.interpolate import RadialBasisFunction

def check_exact_interp():

    """
    Check for exact interpolation varying:
        * Dimensions
        * Kernels
        * Agumenting polynomials
    """
    kernels = [
        # Named scalable
        'multiquadric','inverse','gaussian',
        
        # Named splines
        'linear','cubic','quintic','thin_plate',

        # integer splines
        1,2,3,4,5,
        
        # callable -- odd gaussian-like kernel
        lambda r,eps: np.exp(-r**4/eps**2)]
    
    dims = [None,1,2,3,5,10]
    degrees = [0,1,2,3]
    
    N = 300 # Needs to be big enough for highest dim and highest degree poly
    
    np.random.seed(3101)
    X = np.random.uniform(size=(N,dims[-1]))
    
    for kernel,dim,degree in itertools.product(kernels,dims,degrees):
        if dim is None:
            x = X[:,0] # shape N,
            f = np.sin(2*np.pi*x)
        else:
            x = X[:,:dim]
            f = np.sin(2*np.pi*x).sum(axis=1)
        
        epsilon = None # Use the default
        
        rbf = RadialBasisFunction(x,f,
                                  degree=degree,
                                  kernel=kernel,
                                  epsilon=epsilon)
        y = rbf(x)
        assert_array_almost_equal(f,y,decimal=4) # Especially for higher orders
        #print(kernel,dim,degree,np.max(np.abs(y-f)),np.linalg.norm(y-f)/np.sqrt(N))

if __name__ == '__main__':
    check_exact_interp()
