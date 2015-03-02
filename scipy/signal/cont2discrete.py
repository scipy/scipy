"""
Continuous to discrete transformations for state-space and transfer function.
"""
from __future__ import division, print_function, absolute_import

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# March 29, 2011

import numpy as np
from scipy import linalg

from .ltisys import tf2ss, ss2tf, zpk2ss, ss2zpk

__all__ = ['cont2discrete']


def cont2discrete(sys, dt, method="zoh", alpha=None):
    """
    Transform a continuous to a discrete state-space system.

    Parameters
    ----------
    sys : a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

           * 2: (num, den)
           * 3: (zeros, poles, gain)
           * 4: (A, B, C, D)

    dt : float
        The discretization time step.
    method : {"gbt", "bilinear", "euler", "backward_diff", "zoh"}, optional
        Which method to use:

           * gbt: generalized bilinear transformation
           * bilinear: Tustin's approximation ("gbt" with alpha=0.5)
           * euler: Euler (or forward differencing) method ("gbt" with alpha=0)
           * backward_diff: Backwards differencing ("gbt" with alpha=1.0)
           * zoh: zero-order hold (default)

    alpha : float within [0, 1], optional
        The generalized bilinear transformation weighting parameter, which
        should only be specified with method="gbt", and is ignored otherwise

    Returns
    -------
    sysd : tuple containing the discrete system
        Based on the input type, the output will be of the form

        * (num, den, dt)   for transfer function input
        * (zeros, poles, gain, dt)   for zeros-poles-gain input
        * (A, B, C, D, dt) for state-space system input

    Notes
    -----
    By default, the routine uses a Zero-Order Hold (zoh) method to perform
    the transformation.  Alternatively, a generalized bilinear transformation
    may be used, which includes the common Tustin's bilinear approximation,
    an Euler's method technique, or a backwards differencing technique.

    The Zero-Order Hold (zoh) method is based on [1]_, the generalized bilinear
    approximation is based on [2]_ and [3]_.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Discretization#Discretization_of_linear_state_space_models

    .. [2] http://techteach.no/publications/discretetime_signals_systems/discrete.pdf

    .. [3] G. Zhang, X. Chen, and T. Chen, Digital redesign via the generalized
        bilinear transformation, Int. J. Control, vol. 82, no. 4, pp. 741-754,
        2009.
        (http://www.ece.ualberta.ca/~gfzhang/research/ZCC07_preprint.pdf)

    """
    if len(sys) == 2:
        sysd = cont2discrete(tf2ss(sys[0], sys[1]), dt, method=method,
                             alpha=alpha)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 3:
        sysd = cont2discrete(zpk2ss(sys[0], sys[1], sys[2]), dt, method=method,
                             alpha=alpha)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 4:
        a, b, c, d = sys
    else:
        raise ValueError("First argument must either be a tuple of 2 (tf), "
                         "3 (zpk), or 4 (ss) arrays.")

    if method == 'gbt':
        if alpha is None:
            raise ValueError("Alpha parameter must be specified for the "
                             "generalized bilinear transform (gbt) method")
        elif alpha < 0 or alpha > 1:
            raise ValueError("Alpha parameter must be within the interval "
                             "[0,1] for the gbt method")

    if method == 'gbt':
        # This parameter is used repeatedly - compute once here
        ima = np.eye(a.shape[0]) - alpha*dt*a
        ad = linalg.solve(ima, np.eye(a.shape[0]) + (1.0-alpha)*dt*a)
        bd = linalg.solve(ima, dt*b)

        # Similarly solve for the output equation matrices
        cd = linalg.solve(ima.transpose(), c.transpose())
        cd = cd.transpose()
        dd = d + alpha*np.dot(c, bd)

    elif method == 'bilinear' or method == 'tustin':
        return cont2discrete(sys, dt, method="gbt", alpha=0.5)

    elif method == 'euler' or method == 'forward_diff':
        return cont2discrete(sys, dt, method="gbt", alpha=0.0)

    elif method == 'backward_diff':
        return cont2discrete(sys, dt, method="gbt", alpha=1.0)

    elif method == 'zoh':
        # Build an exponential matrix
        em_upper = np.hstack((a, b))

        # Need to stack zeros under the a and b matrices
        em_lower = np.hstack((np.zeros((b.shape[1], a.shape[0])),
                              np.zeros((b.shape[1], b.shape[1]))))

        em = np.vstack((em_upper, em_lower))
        ms = linalg.expm(dt * em)

        # Dispose of the lower rows
        ms = ms[:a.shape[0], :]

        ad = ms[:, 0:a.shape[1]]
        bd = ms[:, a.shape[1]:]

        cd = c
        dd = d

    else:
        raise ValueError("Unknown transformation method '%s'" % method)

    return ad, bd, cd, dd, dt
