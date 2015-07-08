"""
dltisys - Code related to discrete linear time-invariant systems
"""

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# April 4, 2011
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.interpolate import interp1d
from .ltisys import tf2ss, zpk2ss

__all__ = ['dlsim', 'dstep', 'dimpulse']


def _system_to_statespace(system):
    """
    Return a discrete state-space system from a 3, 4, or 5-tuple input
    """
    if len(system) == 3:
        A, B, C, D = tf2ss(system[0], system[1])
        dt = system[2]
    elif len(system) == 4:
        A, B, C, D = zpk2ss(system[0], system[1], system[2])
        dt = system[3]
    elif len(system) == 5:
        A, B, C, D, dt = system
    else:
        raise ValueError("System argument should be a discrete transfer " +
                         "function, zeros-poles-gain specification, or " +
                         "state-space system")
    return A, B, C, D, dt


def dlsim(system, u, t=None, x0=None):
    """
    Simulate output of a discrete-time linear system.

    Parameters
    ----------
    system : tuple of array_like
        A tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

          - 3: (num, den, dt)
          - 4: (zeros, poles, gain, dt)
          - 5: (A, B, C, D, dt)

    u : array_like
        An input array describing the input at each time `t` (interpolation is
        assumed between given times).  If there are multiple inputs, then each
        column of the rank-2 array represents an input.
    t : array_like, optional
        The time steps at which the input is defined.  If `t` is given, it
        must be the same length as `u`, and the final value in `t` determines
        the number of steps returned in the output.
    x0 : array_like, optional
        The initial conditions on the state vector (zero by default).

    Returns
    -------
    tout : ndarray
        Time values for the output, as a 1-D array.
    yout : ndarray
        System response, as a 1-D array.
    xout : ndarray, optional
        Time-evolution of the state-vector.  Only generated if the input is a
        state-space systems.

    See Also
    --------
    lsim, dstep, dimpulse, cont2discrete

    Examples
    --------
    A simple integrator transfer function with a discrete time step of 1.0
    could be implemented as:

    >>> from scipy import signal
    >>> tf = ([1.0,], [1.0, -1.0], 1.0)
    >>> t_in = [0.0, 1.0, 2.0, 3.0]
    >>> u = np.asarray([0.0, 0.0, 1.0, 1.0])
    >>> t_out, y = signal.dlsim(tf, u, t=t_in)
    >>> y
    array([ 0.,  0.,  0.,  1.])

    """
    A, B, C, D, dt = _system_to_statespace(system)
    u = np.asarray(u)

    # Check for the statespace realization
    # Case of system y=u
    if A.shape[0] == 0:
        A = np.asarray([0])
        B = np.asarray([0])
        C = np.asarray([0])

    if u.ndim == 1:
        u = np.atleast_2d(u).T

    if t is None:
        out_samples = len(u)
        stoptime = (out_samples - 1) * dt
    else:
        stoptime = t[-1]
        out_samples = int(np.floor(stoptime / dt)) + 1

    # Pre-build output arrays
    xout = np.zeros((out_samples, A.shape[0]))
    yout = np.zeros((out_samples, C.shape[0]))
    tout = np.linspace(0.0, stoptime, num=out_samples)

    # Check initial condition
    if x0 is None:
        xout[0, :] = np.zeros((A.shape[0],))
    else:
        xout[0, :] = np.asarray(x0)

    # Pre-interpolate inputs into the desired time steps
    if t is None:
        u_dt = u
    else:
        if len(u.shape) == 1:
            u = u[:, np.newaxis]

        u_dt_interp = interp1d(t, u.transpose(), copy=False, bounds_error=True)
        u_dt = u_dt_interp(tout).transpose()

    # Simulate the system
    for i in range(0, out_samples - 1):
        xout[i+1, :] = np.dot(A, xout[i, :]) + np.dot(B, u_dt[i, :])
        yout[i, :] = np.dot(C, xout[i, :]) + np.dot(D, u_dt[i, :])

    # Last point
    yout[out_samples-1, :] = (np.dot(C, xout[out_samples-1, :]) +
                              np.dot(D, u_dt[out_samples-1, :]))

    if len(system) == 5:
        return tout, yout, xout
    else:
        return tout, yout


def dimpulse(system, x0=None, t=None, n=None):
    """Impulse response of discrete-time system.

    Parameters
    ----------
    system : tuple of array_like
        A tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

          * 3: (num, den, dt)
          * 4: (zeros, poles, gain, dt)
          * 5: (A, B, C, D, dt)

    x0 : array_like, optional
        Initial state-vector.  Defaults to zero.
    t : array_like, optional
        Time points.  Computed if not given.
    n : int, optional
        The number of time points to compute (if `t` is not given).

    Returns
    -------
    tout : ndarray
        Time values for the output, as a 1-D array.
    yout : ndarray
        Impulse response of system.  Each element of the tuple represents
        the output of the system based on an impulse in each input.

    See Also
    --------
    impulse, dstep, dlsim, cont2discrete

    """
    # Determine the system type and set number of inputs and time steps
    A, B, C, D, dt = _system_to_statespace(system)
    n_inputs = B.shape[1]

    # Default to 100 samples if unspecified
    if n is None:
        n = 100

    # If time is not specified, use the number of samples
    # and system dt
    if t is None:
        t = np.linspace(0, n * dt, n, endpoint=False)
    else:
        t = np.asarray(t)

    # For each input, implement a step change
    yout = None
    for i in range(0, n_inputs):
        u = np.zeros((t.shape[0], n_inputs))
        u[0, i] = 1.0

        one_output = dlsim((A, B, C, D, dt), u, t=t, x0=x0)

        if yout is None:
            yout = (one_output[1],)
        else:
            yout = yout + (one_output[1],)

        tout = one_output[0]

    return tout, yout


def dstep(system, x0=None, t=None, n=None):
    """Step response of discrete-time system.

    Parameters
    ----------
    system : tuple of array_like
        A tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

          * 3: (num, den, dt)
          * 4: (zeros, poles, gain, dt)
          * 5: (A, B, C, D, dt)

    x0 : array_like, optional
        Initial state-vector.  Defaults to zero.
    t : array_like, optional
        Time points.  Computed if not given.
    n : int, optional
        The number of time points to compute (if `t` is not given).

    Returns
    -------
    tout : ndarray
        Output time points, as a 1-D array.
    yout : ndarray
        Step response of system.  Each element of the tuple represents
        the output of the system based on a step response to each input.

    See Also
    --------
    step, dimpulse, dlsim, cont2discrete

    """
    # Determine the system type and set number of inputs and time steps
    A, B, C, D, dt = _system_to_statespace(system)
    n_inputs = B.shape[1]

    # Default to 100 samples if unspecified
    if n is None:
        n = 100

    # If time is not specified, use the number of samples
    # and system dt
    if t is None:
        t = np.linspace(0, n * dt, n, endpoint=False)
    else:
        t = np.asarray(t)

    # For each input, implement a step change
    yout = None
    for i in range(0, n_inputs):
        u = np.zeros((t.shape[0], n_inputs))
        u[:, i] = np.ones((t.shape[0],))

        one_output = dlsim((A, B, C, D, dt), u, t=t, x0=x0)

        if yout is None:
            yout = (one_output[1],)
        else:
            yout = yout + (one_output[1],)

        tout = one_output[0]

    return tout, yout
