
"""
    Parameters
    ----------
    A : ndarray or other linear operator
"""
import numpy as np
def norm(A, axis = 1):   
    if hasattr(A, "power"): 
        new_A = A.power(2).sum(axis = axis)
        return np.sqrt(new_A)
    else:
        raise "%s doesn't not support power(scalar)"%type(A)
