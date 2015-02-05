"""
ltisys -- a collection of classes and functions for modeling linear
time invariant systems.
"""
from __future__ import division, print_function, absolute_import

#
# Author: Travis Oliphant 2001
#
# Feb 2010: Warren Weckesser
#   Rewrote lsim2 and added impulse2.
# Aug 2013: Juan Luis Cano
#   Rewrote abcd_normalize.
# Jan 2015: Irvin Probst irvin DOT probst AT ensta-bretagne DOT fr
#   Added pole placement
#

#these two imports are the only one used by place_poles and related functions
#a quick rewrite of the other functions in this file would allow to clean
#up the other imports
import warnings
import numpy as np

import numpy
from numpy import (r_, eye, real, atleast_1d, atleast_2d, poly,
                   squeeze, diag, asarray, product, zeros, array,
                   dot, transpose, ones, zeros_like, linspace, nan_to_num)

from scipy import integrate, interpolate, linalg
from scipy.lib.six import xrange

from .filter_design import tf2zpk, zpk2tf, normalize, freqs


__all__ = ['tf2ss', 'ss2tf', 'abcd_normalize', 'zpk2ss', 'ss2zpk', 'lti',
           'lsim', 'lsim2', 'impulse', 'impulse2', 'step', 'step2', 'bode',
           'freqresp', 'place_poles']


def tf2ss(num, den):
    """Transfer function to state-space representation.

    Parameters
    ----------
    num, den : array_like
        Sequences representing the numerator and denominator polynomials.
        The denominator needs to be at least as long as the numerator.

    Returns
    -------
    A, B, C, D : ndarray
        State space representation of the system, in controller canonical
        form.

    """
    # Controller canonical state-space representation.
    #  if M+1 = len(num) and K+1 = len(den) then we must have M <= K
    #  states are found by asserting that X(s) = U(s) / D(s)
    #  then Y(s) = N(s) * X(s)
    #
    #   A, B, C, and D follow quite naturally.
    #
    num, den = normalize(num, den)   # Strips zeros, checks arrays
    nn = len(num.shape)
    if nn == 1:
        num = asarray([num], num.dtype)
    M = num.shape[1]
    K = len(den)
    if M > K:
        msg = "Improper transfer function. `num` is longer than `den`."
        raise ValueError(msg)
    if M == 0 or K == 0:  # Null system
        return (array([], float), array([], float), array([], float),
                array([], float))

    # pad numerator to have same number of columns has denominator
    num = r_['-1', zeros((num.shape[0], K - M), num.dtype), num]

    if num.shape[-1] > 0:
        D = num[:, 0]
    else:
        D = array([], float)

    if K == 1:
        return array([], float), array([], float), array([], float), D

    frow = -array([den[1:]])
    A = r_[frow, eye(K - 2, K - 1)]
    B = eye(K - 1, 1)
    C = num[:, 1:] - num[:, 0] * den[1:]
    return A, B, C, D


def _none_to_empty_2d(arg):
    if arg is None:
        return zeros((0, 0))
    else:
        return arg


def _atleast_2d_or_none(arg):
    if arg is not None:
        return atleast_2d(arg)


def _shape_or_none(M):
    if M is not None:
        return M.shape
    else:
        return (None,) * 2


def _choice_not_none(*args):
    for arg in args:
        if arg is not None:
            return arg


def _restore(M, shape):
    if M.shape == (0, 0):
        return zeros(shape)
    else:
        if M.shape != shape:
            raise ValueError("The input arrays have incompatible shapes.")
        return M


def abcd_normalize(A=None, B=None, C=None, D=None):
    """Check state-space matrices and ensure they are two-dimensional.

    If enough information on the system is provided, that is, enough
    properly-shaped arrays are passed to the function, the missing ones
    are built from this information, ensuring the correct number of
    rows and columns. Otherwise a ValueError is raised.

    Parameters
    ----------
    A, B, C, D : array_like, optional
        State-space matrices. All of them are None (missing) by default.

    Returns
    -------
    A, B, C, D : array
        Properly shaped state-space matrices.

    Raises
    ------
    ValueError
        If not enough information on the system was provided.

    """
    A, B, C, D = map(_atleast_2d_or_none, (A, B, C, D))

    MA, NA = _shape_or_none(A)
    MB, NB = _shape_or_none(B)
    MC, NC = _shape_or_none(C)
    MD, ND = _shape_or_none(D)

    p = _choice_not_none(MA, MB, NC)
    q = _choice_not_none(NB, ND)
    r = _choice_not_none(MC, MD)
    if p is None or q is None or r is None:
        raise ValueError("Not enough information on the system.")

    A, B, C, D = map(_none_to_empty_2d, (A, B, C, D))
    A = _restore(A, (p, p))
    B = _restore(B, (p, q))
    C = _restore(C, (r, p))
    D = _restore(D, (r, q))

    return A, B, C, D


def ss2tf(A, B, C, D, input=0):
    """State-space to transfer function.

    Parameters
    ----------
    A, B, C, D : ndarray
        State-space representation of linear system.
    input : int, optional
        For multiple-input systems, the input to use.

    Returns
    -------
    num : 2-D ndarray
        Numerator(s) of the resulting transfer function(s).  `num` has one row
        for each of the system's outputs. Each row is a sequence representation
        of the numerator polynomial.
    den : 1-D ndarray
        Denominator of the resulting transfer function(s).  `den` is a sequence
        representation of the denominator polynomial.

    """
    # transfer function is C (sI - A)**(-1) B + D
    A, B, C, D = map(asarray, (A, B, C, D))
    # Check consistency and make them all rank-2 arrays
    A, B, C, D = abcd_normalize(A, B, C, D)

    nout, nin = D.shape
    if input >= nin:
        raise ValueError("System does not have the input specified.")

    # make MOSI from possibly MOMI system.
    B = B[:, input:input + 1]
    D = D[:, input:input + 1]

    try:
        den = poly(A)
    except ValueError:
        den = 1

    if (product(B.shape, axis=0) == 0) and (product(C.shape, axis=0) == 0):
        num = numpy.ravel(D)
        if (product(D.shape, axis=0) == 0) and (product(A.shape, axis=0) == 0):
            den = []
        return num, den

    num_states = A.shape[0]
    type_test = A[:, 0] + B[:, 0] + C[0, :] + D
    num = numpy.zeros((nout, num_states + 1), type_test.dtype)
    for k in range(nout):
        Ck = atleast_2d(C[k, :])
        num[k] = poly(A - dot(B, Ck)) + (D[k] - 1) * den

    return num, den


def zpk2ss(z, p, k):
    """Zero-pole-gain representation to state-space representation

    Parameters
    ----------
    z, p : sequence
        Zeros and poles.
    k : float
        System gain.

    Returns
    -------
    A, B, C, D : ndarray
        State space representation of the system, in controller canonical
        form.

    """
    return tf2ss(*zpk2tf(z, p, k))


def ss2zpk(A, B, C, D, input=0):
    """State-space representation to zero-pole-gain representation.

    Parameters
    ----------
    A, B, C, D : ndarray
        State-space representation of linear system.
    input : int, optional
        For multiple-input systems, the input to use.

    Returns
    -------
    z, p : sequence
        Zeros and poles.
    k : float
        System gain.

    """
    return tf2zpk(*ss2tf(A, B, C, D, input=input))


class lti(object):
    """Linear Time Invariant class which simplifies representation.

    Parameters
    ----------
    args : arguments
        The `lti` class can be instantiated with either 2, 3 or 4 arguments.
        The following gives the number of elements in the tuple and the
        interpretation:

            * 2: (numerator, denominator)
            * 3: (zeros, poles, gain)
            * 4: (A, B, C, D)

        Each argument can be an array or sequence.

    Notes
    -----
    `lti` instances have all types of representations available; for example
    after creating an instance s with ``(zeros, poles, gain)`` the transfer
    function representation (numerator, denominator) can be accessed as
    ``s.num`` and ``s.den``.

    """
    def __init__(self, *args, **kwords):
        """
        Initialize the LTI system using either:

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

    def impulse(self, X0=None, T=None, N=None):
        """
        Return the impulse response of a continuous-time system.
        See `scipy.signal.impulse` for details.
        """
        return impulse(self, X0=X0, T=T, N=N)

    def step(self, X0=None, T=None, N=None):
        """
        Return the step response of a continuous-time system.
        See `scipy.signal.step` for details.
        """
        return step(self, X0=X0, T=T, N=N)

    def output(self, U, T, X0=None):
        """
        Return the response of a continuous-time system to input `U`.
        See `scipy.signal.lsim` for details.
        """
        return lsim(self, U, T, X0=X0)

    def bode(self, w=None, n=100):
        """
        Calculate Bode magnitude and phase data of a continuous-time system.

        Returns a 3-tuple containing arrays of frequencies [rad/s], magnitude
        [dB] and phase [deg]. See `scipy.signal.bode` for details.

        Notes
        -----

        .. versionadded:: 0.11.0

        Examples
        --------
        >>> from scipy import signal
        >>> import matplotlib.pyplot as plt

        >>> s1 = signal.lti([1], [1, 1])
        >>> w, mag, phase = s1.bode()

        >>> plt.figure()
        >>> plt.semilogx(w, mag)    # Bode magnitude plot
        >>> plt.figure()
        >>> plt.semilogx(w, phase)  # Bode phase plot
        >>> plt.show()

        """
        return bode(self, w=w, n=n)

    def freqresp(self, w=None, n=10000):
        """
        Calculate the frequency response of a continuous-time system.

        Returns a 2-tuple containing arrays of frequencies [rad/s] and
        complex magnitude.
        See `scipy.signal.freqresp` for details.

        """
        return freqresp(self, w=w, n=n)


def lsim2(system, U=None, T=None, X0=None, **kwargs):
    """
    Simulate output of a continuous-time linear system, by using
    the ODE solver `scipy.integrate.odeint`.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)

    U : array_like (1D or 2D), optional
        An input array describing the input at each time T.  Linear
        interpolation is used between given times.  If there are
        multiple inputs, then each column of the rank-2 array
        represents an input.  If U is not given, the input is assumed
        to be zero.
    T : array_like (1D or 2D), optional
        The time steps at which the input is defined and at which the
        output is desired.  The default is 101 evenly spaced points on
        the interval [0,10.0].
    X0 : array_like (1D), optional
        The initial condition of the state vector.  If `X0` is not
        given, the initial conditions are assumed to be 0.
    kwargs : dict
        Additional keyword arguments are passed on to the function
        `odeint`.  See the notes below for more details.

    Returns
    -------
    T : 1D ndarray
        The time values for the output.
    yout : ndarray
        The response of the system.
    xout : ndarray
        The time-evolution of the state-vector.

    Notes
    -----
    This function uses `scipy.integrate.odeint` to solve the
    system's differential equations.  Additional keyword arguments
    given to `lsim2` are passed on to `odeint`.  See the documentation
    for `scipy.integrate.odeint` for the full list of arguments.

    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)

    if X0 is None:
        X0 = zeros(sys.B.shape[0], sys.A.dtype)

    if T is None:
        # XXX T should really be a required argument, but U was
        # changed from a required positional argument to a keyword,
        # and T is after U in the argument list.  So we either: change
        # the API and move T in front of U; check here for T being
        # None and raise an exception; or assign a default value to T
        # here.  This code implements the latter.
        T = linspace(0, 10.0, 101)

    T = atleast_1d(T)
    if len(T.shape) != 1:
        raise ValueError("T must be a rank-1 array.")

    if U is not None:
        U = atleast_1d(U)
        if len(U.shape) == 1:
            U = U.reshape(-1, 1)
        sU = U.shape
        if sU[0] != len(T):
            raise ValueError("U must have the same number of rows "
                             "as elements in T.")

        if sU[1] != sys.inputs:
            raise ValueError("The number of inputs in U (%d) is not "
                             "compatible with the number of system "
                             "inputs (%d)" % (sU[1], sys.inputs))
        # Create a callable that uses linear interpolation to
        # calculate the input at any time.
        ufunc = interpolate.interp1d(T, U, kind='linear',
                                     axis=0, bounds_error=False)

        def fprime(x, t, sys, ufunc):
            """The vector field of the linear system."""
            return dot(sys.A, x) + squeeze(dot(sys.B, nan_to_num(ufunc([t]))))
        xout = integrate.odeint(fprime, X0, T, args=(sys, ufunc), **kwargs)
        yout = dot(sys.C, transpose(xout)) + dot(sys.D, transpose(U))
    else:
        def fprime(x, t, sys):
            """The vector field of the linear system."""
            return dot(sys.A, x)
        xout = integrate.odeint(fprime, X0, T, args=(sys,), **kwargs)
        yout = dot(sys.C, transpose(xout))

    return T, squeeze(transpose(yout)), xout


def _cast_to_array_dtype(in1, in2):
    """Cast array to dtype of other array, while avoiding ComplexWarning.

    Those can be raised when casting complex to real.
    """
    if numpy.issubdtype(in2.dtype, numpy.float):
        # dtype to cast to is not complex, so use .real
        in1 = in1.real.astype(in2.dtype)
    else:
        in1 = in1.astype(in2.dtype)

    return in1


def lsim(system, U, T, X0=None, interp=1):
    """
    Simulate output of a continuous-time linear system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)

    U : array_like
        An input array describing the input at each time `T`
        (interpolation is assumed between given times).  If there are
        multiple inputs, then each column of the rank-2 array
        represents an input.
    T : array_like
        The time steps at which the input is defined and at which the
        output is desired.
    X0 :
        The initial conditions on the state vector (zero by default).
    interp : {1, 0}
        Whether to use linear (1) or zero-order hold (0) interpolation.

    Returns
    -------
    T : 1D ndarray
        Time values for the output.
    yout : 1D ndarray
        System response.
    xout : ndarray
        Time-evolution of the state-vector.

    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    U = atleast_1d(U)
    T = atleast_1d(T)
    if len(U.shape) == 1:
        U = U.reshape((U.shape[0], 1))
    sU = U.shape
    if len(T.shape) != 1:
        raise ValueError("T must be a rank-1 array.")
    if sU[0] != len(T):
        raise ValueError("U must have the same number of rows "
                         "as elements in T.")
    if sU[1] != sys.inputs:
        raise ValueError("System does not define that many inputs.")

    if X0 is None:
        X0 = zeros(sys.B.shape[0], sys.A.dtype)

    xout = zeros((len(T), sys.B.shape[0]), sys.A.dtype)
    xout[0] = X0
    A = sys.A
    AT, BT = transpose(sys.A), transpose(sys.B)
    dt = T[1] - T[0]
    lam, v = linalg.eig(A)
    vt = transpose(v)
    vti = linalg.inv(vt)
    GT = dot(dot(vti, diag(numpy.exp(dt * lam))), vt)
    GT = _cast_to_array_dtype(GT, xout)

    ATm1 = linalg.inv(AT)
    ATm2 = dot(ATm1, ATm1)
    I = eye(A.shape[0], dtype=A.dtype)
    GTmI = GT - I
    F1T = dot(dot(BT, GTmI), ATm1)
    if interp:
        F2T = dot(BT, dot(GTmI, ATm2) / dt - ATm1)

    for k in xrange(1, len(T)):
        dt1 = T[k] - T[k - 1]
        if dt1 != dt:
            dt = dt1
            GT = dot(dot(vti, diag(numpy.exp(dt * lam))), vt)
            GT = _cast_to_array_dtype(GT, xout)
            GTmI = GT - I
            F1T = dot(dot(BT, GTmI), ATm1)
            if interp:
                F2T = dot(BT, dot(GTmI, ATm2) / dt - ATm1)

        xout[k] = dot(xout[k - 1], GT) + dot(U[k - 1], F1T)
        if interp:
            xout[k] = xout[k] + dot((U[k] - U[k - 1]), F2T)

    yout = (squeeze(dot(U, transpose(sys.D))) +
            squeeze(dot(xout, transpose(sys.C))))
    return T, squeeze(yout), squeeze(xout)


def _default_response_times(A, n):
    """Compute a reasonable set of time samples for the response time.

    This function is used by `impulse`, `impulse2`, `step` and `step2`
    to compute the response time when the `T` argument to the function
    is None.

    Parameters
    ----------
    A : ndarray
        The system matrix, which is square.
    n : int
        The number of time samples to generate.

    Returns
    -------
    t : ndarray
        The 1-D array of length `n` of time samples at which the response
        is to be computed.
    """
    # Create a reasonable time interval.
    # TODO: This could use some more work.
    # For example, what is expected when the system is unstable?
    vals = linalg.eigvals(A)
    r = min(abs(real(vals)))
    if r == 0.0:
        r = 1.0
    tc = 1.0 / r
    t = linspace(0.0, 7 * tc, n)
    return t


def impulse(system, X0=None, T=None, N=None):
    """Impulse response of continuous-time system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple of array_like
        describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    X0 : array_like, optional
        Initial state-vector.  Defaults to zero.
    T : array_like, optional
        Time points.  Computed if not given.
    N : int, optional
        The number of time points to compute (if `T` is not given).

    Returns
    -------
    T : ndarray
        A 1-D array of time points.
    yout : ndarray
        A 1-D array containing the impulse response of the system (except for
        singularities at zero).

    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    if X0 is None:
        B = sys.B
    else:
        B = sys.B + X0
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)
    else:
        T = asarray(T)

    h = zeros(T.shape, sys.A.dtype)
    s, v = linalg.eig(sys.A)
    vi = linalg.inv(v)
    C = sys.C
    for k in range(len(h)):
        es = diag(numpy.exp(s * T[k]))
        eA = dot(dot(v, es), vi)
        eA = _cast_to_array_dtype(eA, h)
        h[k] = squeeze(dot(dot(C, eA), B))

    return T, h


def impulse2(system, X0=None, T=None, N=None, **kwargs):
    """
    Impulse response of a single-input, continuous-time linear system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple of array_like
        describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    X0 : 1-D array_like, optional
        The initial condition of the state vector.  Default: 0 (the
        zero vector).
    T : 1-D array_like, optional
        The time steps at which the input is defined and at which the
        output is desired.  If `T` is not given, the function will
        generate a set of time samples automatically.
    N : int, optional
        Number of time points to compute.  Default: 100.
    kwargs : various types
        Additional keyword arguments are passed on to the function
        `scipy.signal.lsim2`, which in turn passes them on to
        `scipy.integrate.odeint`; see the latter's documentation for
        information about these arguments.

    Returns
    -------
    T : ndarray
        The time values for the output.
    yout : ndarray
        The output response of the system.

    See Also
    --------
    impulse, lsim2, integrate.odeint

    Notes
    -----
    The solution is generated by calling `scipy.signal.lsim2`, which uses
    the differential equation solver `scipy.integrate.odeint`.

    .. versionadded:: 0.8.0

    Examples
    --------
    Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = u(t)

    >>> from scipy import signal
    >>> system = ([1.0], [1.0, 2.0, 1.0])
    >>> t, y = signal.impulse2(system)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t, y)

    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    B = sys.B
    if B.shape[-1] != 1:
        raise ValueError("impulse2() requires a single-input system.")
    B = B.squeeze()
    if X0 is None:
        X0 = zeros_like(B)
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)

    # Move the impulse in the input to the initial conditions, and then
    # solve using lsim2().
    ic = B + X0
    Tr, Yr, Xr = lsim2(sys, T=T, X0=ic, **kwargs)
    return Tr, Yr


def step(system, X0=None, T=None, N=None):
    """Step response of continuous-time system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple of array_like
        describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    X0 : array_like, optional
        Initial state-vector (default is zero).
    T : array_like, optional
        Time points (computed if not given).
    N : int
        Number of time points to compute if `T` is not given.

    Returns
    -------
    T : 1D ndarray
        Output time points.
    yout : 1D ndarray
        Step response of system.

    See also
    --------
    scipy.signal.step2

    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)
    else:
        T = asarray(T)
    U = ones(T.shape, sys.A.dtype)
    vals = lsim(sys, U, T, X0=X0)
    return vals[0], vals[1]


def step2(system, X0=None, T=None, N=None, **kwargs):
    """Step response of continuous-time system.

    This function is functionally the same as `scipy.signal.step`, but
    it uses the function `scipy.signal.lsim2` to compute the step
    response.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple of array_like
        describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    X0 : array_like, optional
        Initial state-vector (default is zero).
    T : array_like, optional
        Time points (computed if not given).
    N : int
        Number of time points to compute if `T` is not given.
    kwargs : various types
        Additional keyword arguments are passed on the function
        `scipy.signal.lsim2`, which in turn passes them on to
        `scipy.integrate.odeint`.  See the documentation for
        `scipy.integrate.odeint` for information about these arguments.

    Returns
    -------
    T : 1D ndarray
        Output time points.
    yout : 1D ndarray
        Step response of system.

    See also
    --------
    scipy.signal.step

    Notes
    -----
    .. versionadded:: 0.8.0
    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)
    else:
        T = asarray(T)
    U = ones(T.shape, sys.A.dtype)
    vals = lsim2(sys, U, T, X0=X0, **kwargs)
    return vals[0], vals[1]


def bode(system, w=None, n=100):
    """
    Calculate Bode magnitude and phase data of a continuous-time system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    w : array_like, optional
        Array of frequencies (in rad/s). Magnitude and phase data is calculated
        for every value in this array. If not given a reasonable set will be
        calculated.
    n : int, optional
        Number of frequency points to compute if `w` is not given. The `n`
        frequencies are logarithmically spaced in an interval chosen to
        include the influence of the poles and zeros of the system.

    Returns
    -------
    w : 1D ndarray
        Frequency array [rad/s]
    mag : 1D ndarray
        Magnitude array [dB]
    phase : 1D ndarray
        Phase array [deg]

    Notes
    -----

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> s1 = signal.lti([1], [1, 1])
    >>> w, mag, phase = signal.bode(s1)

    >>> plt.figure()
    >>> plt.semilogx(w, mag)    # Bode magnitude plot
    >>> plt.figure()
    >>> plt.semilogx(w, phase)  # Bode phase plot
    >>> plt.show()

    """
    w, y = freqresp(system, w=w, n=n)

    mag = 20.0 * numpy.log10(abs(y))
    phase = numpy.unwrap(numpy.arctan2(y.imag, y.real)) * 180.0 / numpy.pi

    return w, mag, phase


def freqresp(system, w=None, n=10000):
    """Calculate the frequency response of a continuous-time system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:

            * 2 (num, den)
            * 3 (zeros, poles, gain)
            * 4 (A, B, C, D)

    w : array_like, optional
        Array of frequencies (in rad/s). Magnitude and phase data is
        calculated for every value in this array. If not given a reasonable
        set will be calculated.
    n : int, optional
        Number of frequency points to compute if `w` is not given. The `n`
        frequencies are logarithmically spaced in an interval chosen to
        include the influence of the poles and zeros of the system.

    Returns
    -------
    w : 1D ndarray
        Frequency array [rad/s]
    H : 1D ndarray
        Array of complex magnitude values

    Examples
    --------
    # Generating the Nyquist plot of a transfer function

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> s1 = signal.lti([], [1, 1, 1], [5])
    # transfer function: H(s) = 5 / (s-1)^3

    >>> w, H = signal.freqresp(s1)

    >>> plt.figure()
    >>> plt.plot(H.real, H.imag, "b")
    >>> plt.plot(H.real, -H.imag, "r")
    >>> plt.show()
    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)

    if sys.inputs != 1 or sys.outputs != 1:
        raise ValueError("freqresp() requires a SISO (single input, single "
                         "output) system.")

    if w is not None:
        worN = w
    else:
        worN = n

    # In the call to freqs(), sys.num.ravel() is used because there are
    # cases where sys.num is a 2-D array with a single row.
    w, h = freqs(sys.num.ravel(), sys.den, worN=worN)

    return w, h
    
    
def _valid_inputs(A,B,poles, method):
    """
    Check the poles come in complex conjugage pairs
    Check shapes of A,B and poles are compatible.
    Check the method chosen is compatible with provided poles
    Return update method to use and ordered poles
    """
    
    #will raise ValueError if poles do not come in complex conjugates pairs   
    poles = _order_complex_poles(poles)    
    
    if len(A.shape) > 2:
        raise ValueError("A must be a 2D array/matrix.")
    if len(B.shape) > 2:
        raise ValueError("B must be a 2D array/matrix")
    if A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if len(poles) > A.shape[0]:
        raise ValueError("maximum number of poles is %d but you asked for %d" % (A.shape[0],len(poles)))
    if len(poles) < A.shape[0]:
        raise ValueError("number of poles is %d but you should provide %d" % (len(poles),A.shape[0]))
    
    r = matrix_rank(B)
    for p in poles:
        if sum(p == poles) > r:
            raise ValueError("at least one of the requested pole is repeated more than rank(B) times")
    #choose update method
    update_xj = _YT
    if method == "KNV0":
        update_xj = _KNV0
        if not all(isreal(poles)):
            raise ValueError("Complex poles are not supported by KNV0")
    return update_xj, poles


def _order_complex_poles(poles):
    """
    Check we have complex conjugates pairs and reorder P according to YT, ie
    real_poles,complex_i,conjugate complex_i, ....
    """
    
    ordered_poles = list(sort(poles[isreal(poles)]))
    for p in poles[imag(poles) < 0]:
        if conj(p) in poles:
            ordered_poles.extend((p,conj(p)))
    if poles.shape[0] != len(ordered_poles):
        raise ValueError("Complex poles must come with their conjugates")
        
    return array(ordered_poles)    
    
def _KNV0(B,ker_pole,transfer_matrix,j,poles):
    """
    Algorithm "KNV0" Kautsky et Al. Robust pole 
    assignment in linear state feedback,Int journal of Control
    1985, vol 41 p 1129->1155
    http://la.epfl.ch/files/content/sites/la/files/users/105941/public/KautskyNicholsDooren
    """
      
    #remove xj form the base
    transfer_matrix_not_j = delete(transfer_matrix,j,axis=1)
    #if we QR this matrix in full mode Q=Q0|Q1
    #then Q1 will be a single column orthogonnal to
    #Q0, that's what we are looking for !
    Q,R = linalg.qr(transfer_matrix_not_j, mode="full") 
    mat_ker_pj = dot(ker_pole[j],ker_pole[j].T)
    yj = dot(mat_ker_pj,Q[:,-1])

    #if Q[:,-1] is "almost" orthogonnal to kerP[j] its
    #projection into kerP[j] will yield a vector
    #close to 0. As we are looking for a vector in kerP[j]
    #simply stick with transfer_matrix_j (unless someone provides me with
    #a better choice ? mailto: irvin DOT probst AT ensta-bretagne DOT fr)

    if not allclose(yj,0):
        xj = yj/linalg.norm(yj)
        if isreal(poles[j]):
            transfer_matrix[:,j] = xj
#        else:
#            #KNV does not support complex poles, using YT technique
#            #the two lines below seem to work 9 ou of 10 times but
#            #it is not reliable enough.
#            #if you have any idea on how to improve this mailto:
#            #irvin DOT probst AT ensta-bretagne DOT fr
            
#            transfer_matrix[:,j]=real(xj)
#            transfer_matrix[:,j+1]=imag(xj)
# ADD THIS at the beginning of this function if you wish to test complex support            
#    removed the two lines below as complex are not supported    
#    if ~np.isreal(P[j]) and (j>=B.shape[0]-1 or P[j]!=np.conj(P[j+1])):
#        return

def _YT_real(ker_pole,Q,transfer_matrix,i,j):
    """
    Applies algorithm from YT section 6.1 page 19 related to real pairs
    """
    #The variables' name is as close as possible to their names in the original paper
    #imo => i minus one => i-1
    
    #step 1 page 19    
    u = Q[:,-2,newaxis]
    v = Q[:,-1,newaxis]
    
    #step 2 page 19
    m = dot(dot(ker_pole[i].T,dot(u,v.T)-dot(v,u.T)),ker_pole[j])

    #step 3 page 19    
    um,sm,vm = linalg.svd(m)
    #mu1, mu2 two first columns of U => 2 first lines of U.T
    mu1,mu2 = um.T[:2,:,newaxis]
    #VM is V.T with numpy we want the first two lines of V.T
    nu1,nu2 = vm[:2,:,newaxis]
    
    #what follows is a rough python translation of the formulas
    #in section 6.2 page 20 (step 4)
    transfer_matrix_j_mo_transfer_matrix_j = \
        vstack((transfer_matrix[:,i,newaxis],transfer_matrix[:,j,newaxis])) 
        
    if not allclose(sm[0],sm[1]):
        ker_pole_imo_mu1 = dot(ker_pole[i],mu1)
        ker_pole_i_nu1 = dot(ker_pole[j],nu1)
        ker_pole_mu_nu = vstack((ker_pole_imo_mu1,ker_pole_i_nu1))    
    else:     
        ker_pole_ij = vstack((hstack((ker_pole[i],zeros(ker_pole[i].shape))),
                           hstack((zeros(ker_pole[j].shape),ker_pole[j]))))
        mu_nu_matrix = vstack((hstack((mu1,mu2)),
                            hstack((nu1,nu2))))  
        ker_pole_mu_nu = dot(ker_pole_ij,mu_nu_matrix)
    transfer_matrix_ij = dot(dot(ker_pole_mu_nu,ker_pole_mu_nu.T),
                            transfer_matrix_j_mo_transfer_matrix_j)
    if not allclose(transfer_matrix_ij,0):
        transfer_matrix_ij = sqrt(2)*transfer_matrix_ij/linalg.norm(transfer_matrix_ij)        
        transfer_matrix[:,i] = transfer_matrix_ij[:transfer_matrix[:,i].shape[0],0]
        transfer_matrix[:,j] = transfer_matrix_ij[transfer_matrix[:,i].shape[0]:,0]
    else:
        #as in knv0 if transfer_matrix_j_mo_transfer_matrix_j is orthogonal to Vect{ker_pole_mu_nu}
        #assign transfer_matrixi/transfer_matrix_j to ker_pole_mu_nu and iterate. As we are looking
        #for a vector in Vect{Matker_pole_MU_NU} (see section 6.1 page 19)
        #this might help (that's a guess, not a claim !)
        transfer_matrix[:,i] = ker_pole_mu_nu[:transfer_matrix[:,i].shape[0],0]
        transfer_matrix[:,j] = ker_pole_mu_nu[transfer_matrix[:,i].shape[0]:,0]
    
def _YT_complex(ker_pole,Q,transfer_matrix,i,j):
    """
    Applies algorithm from YT section 6.2 page 20 related to complex pairs
    """
    #The variables' name is as close as possible to their names in the original paper
    #step 1 page 20
    ur = sqrt(2)*Q[:,-2,newaxis]
    ui = sqrt(2)*Q[:,-1,newaxis]
    u = ur+1j*ui
    v = conj(u)

    #step 2 page 20
    ker_pole_ij=ker_pole[i]
    m = dot(dot(conj(ker_pole_ij.T),dot(u,v.T)-dot(v,u.T)),ker_pole_ij)
   
    #step 3 page 20
    e_val,e_vec = linalg.eig(m)
    #sort eigenvalues according to their module
    e_val_idx = argsort(abs(e_val))
    mu1 = e_vec[:,e_val_idx[-1],newaxis]
    mu2 = e_vec[:,e_val_idx[-2],newaxis]
    
 
    #what follows is a rough python translation of the formulas
    #in section 6.2 page 20 (step 4)
    
    #remember transfer_matrix_i has been split as transfer_matrix[i]=real(transfer_matrix_i)
    #and transfer_matrix[j]=imag(transfer_matrix_i)
    
    transfer_matrix_j_mo_transfer_matrix_j = transfer_matrix[:,i,newaxis]+\
                                            1j*transfer_matrix[:,j,newaxis]


    if not allclose(abs(e_val[e_val_idx[-1]]),abs(e_val[e_val_idx[-2]])):
        ker_pole_mu = dot(ker_pole_ij,mu1)
        ker_pole_mu_star=dot(conj(mu1.T),conj(ker_pole_ij.T))
    else:
        mu1_mu2_matrix = hstack((mu1,mu2))     
        ker_pole_mu = dot(ker_pole_ij,mu1_mu2_matrix)
        ker_pole_mu_star=dot(conj(mu1_mu2_matrix.T),conj(ker_pole_ij.T))
    transfer_matrix_i_j = dot(dot(ker_pole_mu,ker_pole_mu_star),\
                            transfer_matrix_j_mo_transfer_matrix_j)
    
    if not allclose(transfer_matrix_i_j,0):
        transfer_matrix_i_j = transfer_matrix_i_j/linalg.norm(transfer_matrix_i_j) 
        transfer_matrix[:,i] = real(transfer_matrix_i_j[:,0])
        transfer_matrix[:,j] = imag(transfer_matrix_i_j[:,0])
    else:
        #same idea as in YT_real
        transfer_matrix[:,i] = real(ker_pole_mu[:,0])
        transfer_matrix[:,j] = imag(ker_pole_mu[:,0])


    
    
def _YT(B,ker_pole,transfer_matrix,j_main_loop,poles):
    """
    Algorithm "YT" Tits, Yang. Globally Convergent
    Algorithms for Robust Pole Assignment by State Feedback
    http://drum.lib.umd.edu/handle/1903/5598
    The poles P have to be sorted accordingly to section 6.2 page 20
    """
    #to be compatible with KNV main loop and not mess up with
    #yang tits indices from the paper
    
    j = 1+2*j_main_loop
    i = j-1

    #odd number of real poles
    if poles[isreal(poles)].shape[0] % 2:
        if i < B.shape[0] and isreal(poles[i]):
            if i == 0:
                #only one real pole (and more complex to come 
                #or we wouldn't be here anyway)
                return _KNV0(B,ker_pole,transfer_matrix,i,poles)
            elif j < B.shape[0] and ~isreal(poles[j]):
                #we are on the last real pole switch with the first one
                j = 0
        else:  # both poles are complex but we are shifted on the right by one
            i -= 1
            j -= 1

    if j >= B.shape[0]:
        #nothing to be done
        return
    #remove xi and xj form the base, same as QNV method 0 but we remove
    #two vectors instead of one    
    transfer_matrix_not_i_j = delete(transfer_matrix,(i,j),axis=1)
    Q,_ = linalg.qr(transfer_matrix_not_i_j, mode="full")
    
    if isreal(poles[i]):
        _YT_real(ker_pole,Q,transfer_matrix,i,j)
    else:
        _YT_complex(ker_pole,Q,transfer_matrix,i,j) 

def place_poles(A,B,poles, method="YT", maxtry=20, force_maxtry=False, return_poles=False):
    """
    Compute K such as eigenvalues(A-dot(B,K))=P.    
    
    K is the gain matrix such as the plant described by the linear system AX+BU
    will have its closed-loop poles, i.e the eigenvalues A-B*K, as close as possible to
    those asked for in poles.

    Parameters
    ----------
    
    A, B : ndarray
        State-space representation of linear system AX+BU.
        
    poles    : array_like
        Desired real poles and/or complex conjugates poles.
        Complex poles are only supported with method="YT" (default)
    
    method: string (default "YT") 
        Which method to choose to find K, see References and Notes.
        "KNV0": Kautsky, Nichols, Van Dooren update method 0
        "YT": Yang Tits 
    
    maxtry: integer, optionnal (default 20)
        Maximum number of iterations to compute K.
        
    force_matry: boolean, optional (default False)
        By default K will be returned as soon as an iteration does not yield
        any significant change on the determinant of the eigenvalues of A-BK.
        If force_maxtry is set to True K will be returned after maxtry steps.
        

    return_poles: boolean, optional (default False)
        Whether to return or not the actual placed poles.
        
    Returns
    -------
    K : 1D ndarray
        The closed loop matrix such as the eigenvalues of A-BK are as close
        as possible to the requested poles P.
    Pa : 1D ndarray (optional)
        The poles corresponding to A-BK.
        
    References
    ----------
    Kautsky et Al. Robust pole assignment in linear state feedback,
            Int journal of Control 1985, vol 41 p 1129->1155
    Tits, Yang. Globally Convergent Algorithms for Robust Pole Assignment by
            State Feedback http://drum.lib.umd.edu/handle/1903/5598
    
    Notes
    -----
    Tits, Yang (YT) is an update of Kautsky et Al. (KNV) original paper. KNV 
    relies on rank-1 updates to find the transfer matrix X such as 
    X*diag(poles)*inv(X)=A-B*K whereas YT uses rank-2 updates. This 
    yields faster convergence to a solution, furthermore YT algorithm
    supports complex poles whereas KNV does not in its original        
    version. Only update method 0 proposed by KNV has been implemented here
    hence the name "KNV0".
    As the solution to this problem is not unique both methods try to make
    the vectors in X "as orthogonnal as possible to the space spanned by the
    remaining vectors" (page 1145 KNV) in order to maximise the determinant of
    this matrix.
    KNV extended to complex poles is used in Matlab's place function, YT is 
    distributed under a non-free licence by Slicot under the name "robpole".
    It is unclear and undocumented how KNV0 has been extended to complex poles
    (help wanted), therefore only YT supports them in this implementation.
    Using the default method "YT" should be fine in most cases, "KNV0" is only
    provided because it is needed by "YT" in some specific cases. However as 
    the solution to the problem of pole placement is not unique "KNV0" and "YT"
    may yield different solutions with sometimes a better conditioning of the
    transfer matrix X using one method or the other.
    
    
    Example
    --------
    #Test real pole placement using KNV and YT algorithm
    #Example number 1 from section 4 of the reference KNV publication
    
    A=np.array([1.380,-0.2077,6.715,-5.676, -0.5814,-4.290,0,0.6750,
            1.067,4.273,-6.654,5.893, 0.0480,4.273,1.343,-2.104]).reshape(4,4)
    B=np.array([0,5.679,1.136,1.136,0,0,-3.146,0]).reshape(4,2)
    P=np.array((-0.2,-0.5,-5.0566,-8.6659))
        
    #compute K with KNV method 0
    K1,P1=place(A,B,P, method="KNV0",return_poles=True)
    #compute K with YT    
    K2,P2=place(A,B,P, return_poles=True)
    #note K1 and K2 may differ as the solution is not unique


    """

    #play it safe safe with whatever the users might pass as poles
    poles = sort(array(poles).flatten())[::-1]  

    #move away all the inputs checking, it only adds noise to the code
    update_xj,poles = _valid_inputs(A,B,poles,method)

    #Step A: QR decomposition of B page 1132 KNV
    #The variables' name is as close as possible to their names in the original paper
    u,z = linalg.qr(B,mode="full")
    u0 = u[:,:B.shape[1]]
    u1 = u[:,B.shape[1]:]
    z = z[:B.shape[1],:]
    
    #if the solution is unique
    if B.shape[0] == matrix_rank(B):
        #if B is square and full rank there is only one solution 
        #such as (A+BK)=diag(P) i.e BK=diag(P)-A
        #for complex poles we use the following trick
        #
        # |a -b| has for eigenvalues a+b and a-b
        # |b a|
        #
        # |a+bi 0| has the obvious eigenvalues a+bi and a-bi
        # |0 a-bi|
        #
        # e.g solving the first one in R gives the solution 
        # for the second one in C
        
        diag_poles = zeros(A.shape)
        idx = 0
        while idx < poles.shape[0]:
            p = poles[idx]
            diag_poles[idx,idx] = real(p)
            if ~isreal(p) and idx < poles.shape[0]-1:
                diag_poles[idx,idx+1] = -imag(p)
                diag_poles[idx+1,idx+1] = real(p)
                diag_poles[idx+1,idx] = imag(p)
                idx += 1  # skip next one
            idx += 1
        gain_matrix = linalg.lstsq(B,diag_poles-A)[0]

    else:        
        #step A (p1144 KNV) and beginnig of step F: decompose dot(U1.T,A-P[i]*I).T
        #and build our set of transfer_matrix vectors in the same loop
        ker_pole = []
        
        #flag to skip the conjugate of a complex pole
        skip_conjugate = False
        #select orthonormal base ker_pole for each Pole and vectors for transfer_matrix
        for j in range(B.shape[0]):
            if skip_conjugate:
                skip_conjugate = False
                continue
            pole_space_j = dot(u1.T,A-poles[j]*eye(B.shape[0])).T

            #after QR Q=Q0|Q1 
            #only Q0 is used to reconstruct  the qr'ed (dot Q,R) matrix.
            #Q1 is orthogonnal to Q0 and will be multiplied by the zeros in
            #R when using mode "full". In default mode Q1 and the zeros
            #in R are not computed
            Q,_ = linalg.qr(pole_space_j, mode="full")
            ker_pole_j = Q[:,pole_space_j.shape[1]:]
            #we want to select one vector in ker_pole_j to build the transfer
            #matrix transfer_matrix
            #however qr returns sometimes vectors with zeros on the same
            #line for each pole and this yields very long convergence times.
            #Or some other times a set of vectors, one with zero imaginary 
            #part and one (or several) with imaginary parts. After trying 
            #many ways to select the best possible one (eg ditch vectors 
            #with zero imaginary part for complex poles) I ended up summing
            #all vectors in ker_pole_j, this solves 100% of the problems and is
            #still a valid choice for transfer_matrix. Indeed for complex poles we are sure
            #to have a non zero imaginary part that way, and the problem of
            #lines full of zeros in transfer_matrix is solved too as when a vector from ker_pole_j
            #has a zero the other one(s) (when ker_pole_j.shape[1]>1) for sure won't
            #have a zero there.

            transfer_matrix_j = npsum(ker_pole_j,axis=1)[:,newaxis]
            transfer_matrix_j = transfer_matrix_j / linalg.norm(transfer_matrix_j)
            if ~isreal(poles[j]):  # complex pole
                transfer_matrix_j = hstack([real(transfer_matrix_j), imag(transfer_matrix_j)])
                ker_pole.extend([ker_pole_j, ker_pole_j])

                #skip next pole as it is the conjugate                
                skip_conjugate = True
                
            else: #real pole nothing to do     
                ker_pole.append(ker_pole_j)

            if j == 0:
                transfer_matrix = transfer_matrix_j
            else:
                transfer_matrix = hstack((transfer_matrix,transfer_matrix_j))


        if B.shape[1] > 1:  # otherwise there is nothing we can optimize from transfer_matrix computed above
            stop = False
            while maxtry > 0 and not stop:
                transfer_matrix_before = transfer_matrix.copy()
                for j in range(B.shape[0]):
                    update_xj(B,ker_pole,transfer_matrix,j,poles)
                if not force_maxtry:
                    det_transfer_matrix = max(sqrt(spacing(1)),linalg.det(transfer_matrix))
                    det_transfer_matrixb = linalg.det(transfer_matrix_before)
                    if abs((det_transfer_matrix-det_transfer_matrixb)/det_transfer_matrix) < 1e-3:  # convergence test from YT page 21                       
                        stop = True
                maxtry -= 1
        #reconstruct transfer_matrix to match complex conjugate pairs, ie transfer_matrix_j/transfer_matrix_j+1
        #are Re(Complex_pole), Im(Complex_pole) now and will be
        #Re-Im/Re+Im after

        transfer_matrix = transfer_matrix.astype(complex)
        idx = 0  
        while idx < poles.shape[0]-1:
            if ~isreal(poles[idx]):
                rel = transfer_matrix[:,idx].copy()
                img = transfer_matrix[:,idx+1].copy()
                transfer_matrix[:,idx] = rel-1j*img
                transfer_matrix[:,idx+1] = rel+1j*img
                idx += 1  # skip next one
            idx += 1
            
        try:
            m=linalg.solve(transfer_matrix.T,dot(diag(poles),transfer_matrix.T)).T
        except LinAlgError:
            raise ValueError("The poles you've chosen can't be placed")

        gain_matrix = linalg.solve(z,dot(u0.T,m-A))
        
    #Beware Kautsky solves A+BK but the usual form is A-BK
    gain_matrix = -gain_matrix
    #K still contains complex with ~=0j imaginary parts, get rid of them
    gain_matrix = real(gain_matrix) 

    return gain_matrix, linalg.eig(A-dot(B,gain_matrix))[0]
