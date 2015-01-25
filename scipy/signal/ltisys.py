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
#

from .filter_design import tf2zpk, zpk2tf, normalize, freqs
import numpy
from numpy import (product, zeros, array, dot, transpose, ones,
                   nan_to_num, zeros_like, linspace)
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy.linalg as linalg
from scipy._lib.six import xrange
from numpy import (r_, eye, real, atleast_1d, atleast_2d, poly,
                   squeeze, diag, asarray)

__all__ = ['tf2ss', 'ss2tf', 'abcd_normalize', 'zpk2ss', 'ss2zpk', 'lti',
           'lsim', 'lsim2', 'impulse', 'impulse2', 'step', 'step2', 'bode',
           'freqresp']


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
