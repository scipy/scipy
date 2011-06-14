"""
Continuous to discrete transformations for state-space and transfer function.
"""

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# March 29, 2011

import numpy as np
import numpy.linalg
import scipy.linalg

from ltisys import tf2ss, ss2tf, zpk2ss, ss2zpk

__all__ = ['cont2discrete']


def _mrdivide(b,a):
    """Convenience function for matrix divides"""
    s = np.linalg.solve(a.transpose(), b.transpose())
    return s.transpose()

def cont2discrete(sys, dt, method="zoh"):
    """Transform a continuous to a discrete state-space system.

    Parameters
    -----------
    sys : a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:
        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)
        
    dt : float
        The discretization time step.
    method : {"bilinear", "zoh"}
        Which method to use, bilinear or zero-order hold ("zoh", the default).

    Returns
    -------
    sysd : tuple containing the discrete system
        Based on the input type, the output will be of the form
        
        (num, den, dt)   for transfer function input
        (zeros, poles, gain, dt)   for zeros-poles-gain input
        (A, B, C, D, dt) for state-space system input

    Notes
    -----
    By default, the routine uses a Zero-Order Hold (zoh) method to perform
    the transformation.  Alternatively, Tustin's bilinear approximation can
    be used.

    The Zero-Order Hold (zoh) method is based on:
    http://en.wikipedia.org/wiki/Discretization#Discretization_of_linear_state_space_models

    Tustin's bilinear approximation is based on:
    http://techteach.no/publications/discretetime_signals_systems/discrete.pdf

    """
    if len(sys) == 2:
        sysd = cont2discrete(tf2ss(sys[0], sys[1]), dt, method=method)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 3:
        sysd = cont2discrete(zpk2ss(sys[0], sys[1], sys[2]), dt, method=method)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 4:
        a, b, c, d = sys
    else:
        raise ValueError("First argument must either be a tuple of 2 (tf) "
                         "or 4 (ss) arrays.")

    if method == 'bilinear':
        # Compute the term (2/dt)*I
        itv = 2.0 / dt * np.eye(a.shape[0])

        # Solve for Ad
        ad = _mrdivide((itv + a), (itv - a))

        # Solve for Bd using a linear solver to avoid direct inversion
        iab = np.linalg.solve((itv - a), b)
        tk = 2.0 / dt
        bd = tk * iab

        # Similarly solve for the output equation matrices
        cd = 2.0 * _mrdivide(c, (itv - a))
        dd = d + np.dot(c, iab)

    elif method == 'zoh':
        # Build an exponential matrix
        em_upper = np.hstack((a, b))

        # Need to stack zeros under the a and b matrices
        em_lower = np.hstack((np.zeros((b.shape[1], a.shape[0])),
                              np.zeros((b.shape[1], b.shape[1])) ))

        em = np.vstack((em_upper, em_lower))
        ms = scipy.linalg.expm(dt * em)

        # Dispose of the lower rows
        ms = ms[:a.shape[0], :]

        ad = ms[:, 0:a.shape[1]]
        bd = ms[:, a.shape[1]:]

        cd = c
        dd = d

    else:
        raise ValueError("Unknown transformation method.")

    return ad, bd, cd, dd, dt
