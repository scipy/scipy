"""
dltisys - Code related to discrete linear time-invariant systems
"""

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# April 4, 2011
from __future__ import division, print_function, absolute_import

from .filter_design import tf2zpk, zpk2tf, normalize, freqs
import numpy as np
from scipy.interpolate import interp1d
from .ltisys import tf2ss, ss2tf, zpk2ss, ss2zpk

__all__ = ['dlsim', 'dstep', 'dimpulse', 'dlti']


class dlti(object):
    """Discrete Linear Time Invariant class which simplifies representation.

    Parameters
    ----------
    args : arguments
        The `dlti` class can be instantiated with either 2, 3 or 4 arguments.
        The following gives the number of elements in the tuple and the
        interpretation:

            * 2: (numerator, denominator)
            * 3: (zeros, poles, gain)
            * 4: (A, B, C, D)

        Each argument can be an array or sequence.

    Notes
    -----
    `dlti` instances have all types of representations available; for example
    after creating an instance z with ``(zeros, poles, gain)`` the transfer
    function representation (numerator, denominator) can be accessed as
    ``z.num`` and ``z.den``.

    """

    def __init__(self, *args, **kwords):
        """
        Initialize the DLTI system using either:

            - (numerator, denominator)
            - (zeros, poles, gain)
            - (A, B, C, D) : state-space.

        """

        N = len(args)
        if N == 2:  # Numerator denominator transfer function input
            self._num, self._den = normalize(*args)
            self._update(N)
            self.inputs = 1
            if len(self.num.shape) > 1:
                self.outputs = self.num.shape[0]
            else:
                self.outputs = 1
        elif N == 3:      # Zero-pole-gain form
            self._zeros, self._poles, self._gain = args
            self._update(N)
            # make sure we have numpy arrays
            self.zeros = numpy.asarray(self.zeros)
            self.poles = numpy.asarray(self.poles)
            self.inputs = 1
            if len(self.zeros.shape) > 1:
                self.outputs = self.zeros.shape[0]
            else:
                self.outputs = 1
        elif N == 4:       # State-space form
            self._A, self._B, self._C, self._D = abcd_normalize(*args)
            self._update(N)
            self.inputs = self.B.shape[-1]
            self.outputs = self.C.shape[0]
        else:
            raise ValueError("Needs 2, 3, or 4 arguments.")

    def __repr__(self):
        """
        Canonical representation using state-space to preserve numerical
        precision and any MIMO information
        """
        return '{0}(\n{1},\n{2},\n{3},\n{4}\n)'.format(
            self.__class__.__name__,
            repr(self.A),
            repr(self.B),
            repr(self.C),
            repr(self.D),
            )

    @property
    def num(self):
        return self._num

    @num.setter
    def num(self, value):
        self._num = value
        self._update(2)

    @property
    def den(self):
        return self._den

    @den.setter
    def den(self, value):
        self._den = value
        self._update(2)

    @property
    def zeros(self):
        return self._zeros

    @zeros.setter
    def zeros(self, value):
        self._zeros = value
        self._update(3)

    @property
    def poles(self):
        return self._poles

    @poles.setter
    def poles(self, value):
        self._poles = value
        self._update(3)

    @property
    def gain(self):
        return self._gain

    @gain.setter
    def gain(self, value):
        self._gain = value
        self._update(3)

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, value):
        self._A = value
        self._update(4)

    @property
    def B(self):
        return self._B

    @B.setter
    def B(self, value):
        self._B = value
        self._update(4)

    @property
    def C(self):
        return self._C

    @C.setter
    def C(self, value):
        self._C = value
        self._update(4)

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, value):
        self._D = value
        self._update(4)

    def _update(self, N):
        if N == 2:
            self._zeros, self._poles, self._gain = tf2zpk(self.num, self.den)
            self._A, self._B, self._C, self._D = tf2ss(self.num, self.den)
        if N == 3:
            self._num, self._den = zpk2tf(self.zeros, self.poles, self.gain)
            self._A, self._B, self._C, self._D = zpk2ss(self.zeros,
                                                        self.poles, self.gain)
        if N == 4:
            self._num, self._den = ss2tf(self.A, self.B, self.C, self.D)
            self._zeros, self._poles, self._gain = ss2zpk(self.A, self.B,
                                                          self.C, self.D)


    def dimpulse(self, x0=None, t=None, n=None):
        """
        Return the impulse response of a discrete-time system.

        """
        return impulse(self, x0=x0, t=t, n=n)

    def dstep(self, x0=None, t=None, n=None):
        """
        Return the step response of a continuous-time system.

        """
        return dstep(self, x0=x0, t=t, n=n)

def dlsim(system, u, t=None, x0=None):
    """
    Simulate output of a discrete-time linear system.

    Parameters
    ----------
    system : class instance or tuple
        An instance of the LTI class, or a tuple describing the system.
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
        The time steps at which the input is defined.  If `t` is given, the
        final value in `t` determines the number of steps returned in the
        output.
    x0 : arry_like, optional
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
    if len(system) == 3:
        a, b, c, d = tf2ss(system[0], system[1])
        dt = system[2]
    elif len(system) == 4:
        a, b, c, d = zpk2ss(system[0], system[1], system[2])
        dt = system[3]
    elif len(system) == 5:
        a, b, c, d, dt = system
    else:
        raise ValueError("System argument should be a discrete transfer " +
                         "function, zeros-poles-gain specification, or " +
                         "state-space system")

    if t is None:
        out_samples = max(u.shape)
        stoptime = (out_samples - 1) * dt
    else:
        stoptime = t[-1]
        out_samples = int(np.floor(stoptime / dt)) + 1

    # Pre-build output arrays
    xout = np.zeros((out_samples, a.shape[0]))
    yout = np.zeros((out_samples, c.shape[0]))
    tout = np.linspace(0.0, stoptime, num=out_samples)

    # Check initial condition
    if x0 is None:
        xout[0, :] = np.zeros((a.shape[1],))
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
        xout[i + 1, :] = np.dot(a, xout[i, :]) + np.dot(b, u_dt[i, :])
        yout[i, :] = np.dot(c, xout[i, :]) + np.dot(d, u_dt[i, :])

    # Last point
    yout[out_samples - 1, :] = np.dot(c, xout[out_samples - 1, :]) + \
                            np.dot(d, u_dt[out_samples - 1, :])

    if len(system) == 5:
        return tout, yout, xout
    else:
        return tout, yout


def dimpulse(system, x0=None, t=None, n=None):
    """Impulse response of discrete-time system.

    Parameters
    ----------
    system : tuple
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
    t : ndarray
        A 1-D array of time points.
    yout : tuple of array_like
        Impulse response of system.  Each element of the tuple represents
        the output of the system based on an impulse in each input.

    See Also
    --------
    impulse, dstep, dlsim, cont2discrete

    """
    # Determine the system type and set number of inputs and time steps
    if len(system) == 3:
        n_inputs = 1
        dt = system[2]
    elif len(system) == 4:
        n_inputs = 1
        dt = system[3]
    elif len(system) == 5:
        n_inputs = system[1].shape[1]
        dt = system[4]
    else:
        raise ValueError("System argument should be a discrete transfer " +
                         "function, zeros-poles-gain specification, or " +
                         "state-space system")

    # Default to 100 samples if unspecified
    if n is None:
        n = 100

    # If time is not specified, use the number of samples
    # and system dt
    if t is None:
        t = np.linspace(0, n * dt, n, endpoint=False)

    # For each input, implement a step change
    yout = None
    for i in range(0, n_inputs):
        u = np.zeros((t.shape[0], n_inputs))
        u[0, i] = 1.0

        one_output = dlsim(system, u, t=t, x0=x0)

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
    system : a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

          * 3: (num, den, dt)
          * 4: (zeros, poles, gain, dt)
          * 5: (A, B, C, D, dt)

    x0 : array_like, optional
        Initial state-vector (default is zero).
    t : array_like, optional
        Time points (computed if not given).
    n : int, optional
        Number of time points to compute if `t` is not given.

    Returns
    -------
    t : ndarray
        Output time points, as a 1-D array.
    yout : tuple of array_like
        Step response of system.  Each element of the tuple represents
        the output of the system based on a step response to each input.

    See Also
    --------
    step, dimpulse, dlsim, cont2discrete

    """
    # Determine the system type and set number of inputs and time steps
    if len(system) == 3:
        n_inputs = 1
        dt = system[2]
    elif len(system) == 4:
        n_inputs = 1
        dt = system[3]
    elif len(system) == 5:
        n_inputs = system[1].shape[1]
        dt = system[4]
    else:
        raise ValueError("System argument should be a discrete transfer " +
                         "function, zeros-poles-gain specification, or " +
                         "state-space system")

    # Default to 100 samples if unspecified
    if n is None:
        n = 100

    # If time is not specified, use the number of samples
    # and system dt
    if t is None:
        t = np.linspace(0, n * dt, n, endpoint=False)

    # For each input, implement a step change
    yout = None
    for i in range(0, n_inputs):
        u = np.zeros((t.shape[0], n_inputs))
        u[:, i] = np.ones((t.shape[0],))

        one_output = dlsim(system, u, t=t, x0=x0)

        if yout is None:
            yout = (one_output[1],)
        else:
            yout = yout + (one_output[1],)

        tout = one_output[0]

    return tout, yout
