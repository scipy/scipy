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
#imports added for place
from numpy import isreal, imag, newaxis, hstack, sort, delete, \
    allclose, sqrt, conj, vstack, spacing, argsort, sum as npsum
#for whatever reason I can't find it through scipy.linalg
from numpy.linalg import matrix_rank
from numpy.linalg.linalg import LinAlgError


__all__ = ['tf2ss', 'ss2tf', 'abcd_normalize', 'zpk2ss', 'ss2zpk', 'lti',
           'lsim', 'lsim2', 'impulse', 'impulse2', 'step', 'step2', 'bode',
           'freqresp', 'place']


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
    X0 : array_like, optional
        The initial conditions on the state vector (zero by default).
    interp : {1, 0}, optional
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
    N : int, optional
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
    N : int, optional
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
    
    
def _valid_inputs(A,B,P):
    """
    Check shapes of A,B and P are compatible.
    """
    if len(A.shape) > 2:
        raise ValueError("A must be a 2D array/matrix.")
    if len(B.shape) > 2:
        raise ValueError("B must be a 2D array/matrix")
    if A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if B.shape[1] > B.shape[0]:
        raise ValueError("B must have nb col  < = nb lines")
    if len(P) > A.shape[0]:
        raise ValueError("maximum number of poles is %d but you asked for %d" % (A.shape[0],len(P)))
    if len(P) < A.shape[0]:
        raise ValueError("number of poles is %d but you should provide %d" % (len(P),A.shape[0]))
    
    r = matrix_rank(B)
    if r != B.shape[1]:
        raise ValueError("rank of B should be equal to its number of columns")

    for p in P:
        if sum(p == P) > r:
            raise ValueError("at least one of the requested pole is repeated more than rank(B) times")

    
def _order_complex_poles(P):
    """
    Check we have complex conjugates pairs and reorder P according to YT, ie
    real_poles,complex_i,conjugate complex_i, ....
    """
    
    ordered_P = list(sort(P[isreal(P)]))
    for p in P[imag(P) < 0]:
        if conj(p) in P:
            ordered_P.extend((p,conj(p)))
    if P.shape[0] != len(ordered_P):
        raise ValueError("Complex poles must come with their conjugates")
        
    return array(ordered_P)    
    
def _KNV0(B,KerP,X,j,P):
    """
    Algorithm "KNV0" Kautsky et Al. Robust pole 
    assignment in linear state feedback,Int journal of Control
    1985, vol 41 p 1129->1155
    http://la.epfl.ch/files/content/sites/la/files/users/105941/public/KautskyNicholsDooren
    """
      
    #remove xj form the base
    X_not_j = delete(X,j,axis=1)
    #if we QR this matrix in full mode Q=Q0|Q1
    #then Q1 will be a single column orthogonnal to
    #Q0, that's what we are looking for !
    Q,R = linalg.qr(X_not_j, mode="full") 
    MatKerPj = dot(KerP[j],KerP[j].T)
    Yj = dot(MatKerPj,Q[:,-1])

    #if Q[:,-1] is "almost" orthogonnal to kerP[j] its
    #projection into kerP[j] will yield a vector
    #close to 0. As we are looking for a vector in kerP[j]
    #simply stick with Xj (unless someone provides me with
    #a better choice ? mailto: irvin DOT probst AT ensta-bretagne DOT fr)

    if not allclose(Yj,0):
        xj = Yj/linalg.norm(Yj)
        if isreal(P[j]):
            X[:,j] = xj
#        else:
#            #KNV does not support complex poles, using YT technique
#            #the two lines below seem to work 9 ou of 10 times but
#            #it is not reliable enough.
#            #if you have any idea on how to improve this mailto:
#            #irvin DOT probst AT ensta-bretagne DOT fr
            
#            X[:,j]=real(xj)
#            X[:,j+1]=imag(xj)
# ADD THIS at the beginning of this function if you wish to test complex support            
#    removed the two lines below as complex are not supported    
#    if ~np.isreal(P[j]) and (j>=B.shape[0]-1 or P[j]!=np.conj(P[j+1])):
#        return

def _YT_real(KerP,Q,X,i,j):
    """
    Applies algorithm from YT section 6.1 page 19 related to real pairs
    """
    #step 1 page 19
    u = Q[:,-2,newaxis]
    v = Q[:,-1,newaxis]
    #step 2 page 19
    M = dot(dot(KerP[i].T,dot(u,v.T)-dot(v,u.T)),KerP[j])

    #step 3 page 19    
    UM,SM,VM = linalg.svd(M)
    #mu1, mu2 two first columns of U => 2 first lines of U.T
    mu1,mu2 = UM.T[:2,:,newaxis]
    #VM is V.T with numpy we want the first two lines of V.T
    nu1,nu2 = VM[:2,:,newaxis]
    
    #what follows is a rough python translation of the formulas
    #in section 6.2 page 20 (step 4)
    VectXJmo_Xj = vstack((X[:,i,newaxis],X[:,j,newaxis]))    
    if not allclose(SM[0],SM[1]):
        KerPimo_MU1 = dot(KerP[i],mu1)
        KerPi_NU1 = dot(KerP[j],nu1)
        MatKerP_MU_NU = vstack((KerPimo_MU1,KerPi_NU1))    
    else:     
        MatKerP = vstack((hstack((KerP[i],zeros(KerP[i].shape))),
                           hstack((zeros(KerP[j].shape),KerP[j]))))
        MatMU_NU = vstack((hstack((mu1,mu2)),
                            hstack((nu1,nu2))))  
        MatKerP_MU_NU = dot(MatKerP,MatMU_NU)
    Xi = dot(dot(MatKerP_MU_NU,MatKerP_MU_NU.T),VectXJmo_Xj)
    if not allclose(Xi,0):
        Xi = sqrt(2)*Xi/linalg.norm(Xi)        
        X[:,i] = Xi[:X[:,i].shape[0],0]
        X[:,j] = Xi[X[:,i].shape[0]:,0]
    else:
        #as in knv0 if VectXJmo_Xj is orthogonal to Vect{MatKerP_MU_NU}
        #assign Xi/Xj to MatKerP_MU_NU and iterate. As we are looking
        #for a vector in Vect{MatKerP_MU_NU} (see section 6.1 page 19)
        #this might help (that's a guess, not a claim !)
        X[:,i] = MatKerP_MU_NU[:X[:,i].shape[0],0]
        X[:,j] = MatKerP_MU_NU[X[:,i].shape[0]:,0]
    
def _YT_complex(KerP,Q,X,i,j):
    """
    Applies algorithm from YT section 6.2 page 20 related to complex pairs
    """
    #step 1 page 20
    uR = sqrt(2)*Q[:,-2,newaxis]
    uI = sqrt(2)*Q[:,-1,newaxis]
    u = uR+1j*uI
    v = conj(u)

    #step 2 page 20
    M = dot(dot(KerP[i].T,dot(u,v.T)-dot(v,u.T)),KerP[i])
   
    #step 3 page 20
    e_val,e_vec = linalg.eig(M)
    SM = e_val
    mu1 = e_vec[:,0,newaxis]
    mu2 = e_vec[:,1,newaxis]
 
    #what follows is a rough python translation of the formulas
    #in section 6.2 page 20 (step 4)
    
    #remember Xi has been split as X[i]=real(xi) and X[j]=imag(xi)
    xi = X[:,i,newaxis]+1j*X[:,j,newaxis]
    VectXJ = xi

    #sort eigenvalues according to their module
    sm_idx = argsort(abs(SM))
    
    if not allclose(abs(SM[sm_idx[-1]]),abs(SM[[sm_idx[-2]]])):
        MatKerP_MU = dot(KerP[i],mu1)
    else:
        MatMU1_MU2 = hstack((mu1,mu2))     
        MatKerP_MU = dot(KerP[i],MatMU1_MU2)
    Xi = dot(dot(MatKerP_MU,MatKerP_MU.T),VectXJ)
    
    if not allclose(Xi,0):
        Xi = Xi/linalg.norm(Xi) 
        X[:,i] = real(Xi[:,0])
        X[:,j] = imag(Xi[:,0])
    else:
        #same idea as in YT_real
        X[:,i] = real(MatKerP_MU[:,0])
        X[:,j] = imag(MatKerP_MU[:,0])
    
    
def _YT(B,KerP,X,j_main_loop,P):
#def YT(B,KerP,X,i,j):
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
    if P[isreal(P)].shape[0] % 2:
        if i < B.shape[0] and isreal(P[i]):
            if i == 0:
                #only one real pole (and more complex to come 
                #or we wouldn't be here anyway)
                return _KNV0(B,KerP,X,i,P)
            elif j < B.shape[0] and ~isreal(P[j]):
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
    X_not_j = delete(X,(i,j),axis=1)
    Q,_ = linalg.qr(X_not_j, mode="full")
    
    if isreal(P[i]):
        _YT_real(KerP,Q,X,i,j)
    else:
        _YT_complex(KerP,Q,X,i,j) 

def place(A,B,P, method="YT", maxtry=20, force_maxtry=False, return_poles=False):
    """
    Compute K such as eigenvalues(A-dot(B,K))=P.

    Parameters
    ----------
    
    A, B : ndarray
        State-space representation of linear system AX+BU.
        
    P    : array_like
        Desired real poles and/or complex conjugates poles.
        Complex poles are only supported with method="YT" (default)
    
    method: string (default "YT")        
        Which method to choose to find K.
        "KNV0": Kautsky et Al. Robust pole assignment in linear state feedback,
            Int journal of Control 1985, vol 41 p 1129->1155)        
        "YT" (Tits, Yang. Globally Convergent Algorithms for Robust Pole
            Assignment by State Feedback http://drum.lib.umd.edu/handle/1903/5598). 
    
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
    #©Irvin Probst, Ensta Bretagne, 2014
    #changelog:
    #2/12/2014 
    #-initial post in scipy-dev ML
    #3/12/2014
    #-when B is square we don't need to run the optimisation
    #algorithms as there is only one solution => linalg.solve
    #-added normalization of Xj vectors 
    #-added Yang Tits optimisation of KNV method 0
    #-removed useless KNV method 1
    #4/12/2014
    #-use bette variable names (S->Ker)
    #5/12/2014
    #-finished Yang Tits
    #-code clean up
    #-improved convergence check on X to stop the main loop earlier if
    #  we have found a solution (is that really useful ?)
    #8/12/2014 
    #-added complex support for YT, deduced how to make it work for KNV0 too
    #9/12/2014
    #-fixed *again* convergence test
    #10->12/12/2014
    #-bugfixes, bugfixes, bugfixes...
    #-removed complex support from KNV0 as it is not reliable enough 
    #17/12/2014
    #-improved transfer matrix initial selection

    #play it safe safe with whatever the users might pass as P
    P = sort(array(P).flatten())[::-1]  
   
    #will return False if poles do not come in complex conjugates pairs   
    P = _order_complex_poles(P)
    
    #move away all the inputs checking, it only adds noise to the code
    _valid_inputs(A,B,P)

    #Step A: QR decomposition of B page 1132 KNV
    U,Z = linalg.qr(B,mode="full")
    U0 = U[:,:B.shape[1]]
    U1 = U[:,B.shape[1]:]
    Z = Z[:B.shape[1],:]
    
    I = eye(B.shape[0]) 
    #if B is square there is almost nothing to do
    if B.shape[0] == B.shape[1]:
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
        
        D = zeros(A.shape)
        idx = 0
        while idx < P.shape[0]:
            p = P[idx]
            D[idx,idx] = real(p)
            if ~isreal(p) and idx < P.shape[0]-1:
                D[idx,idx+1] = -imag(p)
                D[idx+1,idx+1] = real(p)
                D[idx+1,idx] = imag(p)
                idx += 1  # skip next one
            idx += 1
        K = linalg.solve(B,D-A)
    else:        
        #step A (p1144 KNV) and beginnig of step F: decompose dot(U1.T,A-P[i]*I).T
        #and build our set of X vectors in the same loop
        KerP = []
        
        #flag to skip the conjugate of a complex pole
        skip_conjugate = False
        #select orthonormal base KerP for each Pole and vectors for X
        for j in range(B.shape[0]):
            if skip_conjugate:
                skip_conjugate = False
                continue
            MatU1_A_Pi = dot(U1.T,A-P[j]*I).T

            #after QR Q=Q0|Q1 
            #only Q0 is used to reconstruct  the qr'ed (dot Q,R) matrix.
            #Q1 is orthogonnal to Q0 and will be multiplied by the zeros in
            #R when using mode "full". In default mode Q1 and the zeros
            #in R are not computed
            Q,_ = linalg.qr(MatU1_A_Pi, mode="full")
            KerPj = Q[:,MatU1_A_Pi.shape[1]:]
            
            #we want to select one vector in KerPj to build the transfer
            #matrix X
            #however qr returns sometimes vectors with zeros on the same
            #line for each pole and this yields very long convergence times.
            #Or some other times a set of vectors, one with zero imaginary 
            #part and one (or several) with imaginary parts. After trying 
            #many ways to select the best possible one (eg ditch vectors 
            #with zero imaginary part for complex poles) I ended up summing
            #all vectors in KerPj, this solves 100% of the problems and is
            #still a valid choice for X. Indeed for complex poles we are sure
            #to have a non zero imaginary part that way, and the problem of
            #lines full of zeros in X is solved too as when a vector from KerPj
            #has a zero the other one(s) (when KerPj.shape[1]>1) for sure won't
            #have a zero there.

            Xj = npsum(KerPj,axis=1)[:,newaxis]
            Xj = Xj / linalg.norm(Xj)
            if ~isreal(P[j]):  # complex pole
                Xj = hstack([real(Xj), imag(Xj)])
                KerP.extend([KerPj, KerPj])

                #skip next pole as it is the conjugate                
                skip_conjugate = True
                
            else: #real pole nothing to do     
                KerP.append(KerPj)

            if j == 0:
                X = Xj
            else:
                X = hstack((X,Xj))

        #choose update method
        update_xj = _YT
        if method == "KNV0":
            update_xj = _KNV0
            if not all(isreal(P)):
                raise ValueError("Complex poles are not supported by KNV0")  

        if B.shape[1] > 1:  # otherwise there is nothing we can optimize from X computed above
            stop = False
            while maxtry > 0 and not stop:      
                X_before = X.copy()
                for j in range(B.shape[0]):
                    update_xj(B,KerP,X,j,P)
                if not force_maxtry:
                    detX = max(sqrt(spacing(1)),linalg.det(X))
                    detXb = linalg.det(X_before)
                    if abs((detX-detXb)/detX) < 1e-3:  # convergence test from YT page 21                       
                        stop = True
                maxtry -= 1

        #reconstruct X to match complex conjugate pairs, ie Xj/Xj+1
        #are Re(Complex_pole), Im(Complex_pole) now and will be
        #Re-Im/Re+Im after

        X = X.astype(complex)
        idx = 0  
        while idx < P.shape[0]-1:
            if ~isreal(P[idx]):
                rel = X[:,idx].copy()
                img = X[:,idx+1].copy()
                X[:,idx] = rel-1j*img
                X[:,idx+1] = rel+1j*img
                idx += 1  # skip next one
            idx += 1
            
        try:
            M=linalg.solve(X.T,dot(diag(P),X.T)).T
        except LinAlgError:
            raise ValueError("The poles you've chosen can't be placed")

        K = linalg.solve(Z,dot(U0.T,M-A))
        
    #Beware Kautsky solves A+BK but the usual form is A-BK
    K = -K
    #K still contains complex with ~=0j imaginary parts, get rid of them
    K = real(K) 

    if return_poles:
        Pa,_ = linalg.eig(A-dot(B,K))
        return K, Pa
    else:
        return K
