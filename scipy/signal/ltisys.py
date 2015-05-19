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
# Mar 2015: Clancy Rowley
#   Rewrote lsim
# May 2015: Felix Berkenkamp
#   Split lti class into subclasses
#

import warnings
import numpy as np

#np.linalg.qr fails on some tests with LinAlgError: zgeqrf returns -7
#use scipy's qr until this is solved

from scipy.linalg import qr as s_qr


import numpy
from numpy import (r_, eye, real, atleast_1d, atleast_2d, poly,
                   squeeze, diag, asarray, product, zeros, array,
                   dot, transpose, ones, zeros_like, linspace, nan_to_num)
import copy

from scipy import integrate, interpolate, linalg
from scipy._lib.six import xrange

from .filter_design import tf2zpk, zpk2tf, normalize, freqs


__all__ = ['tf2ss', 'ss2tf', 'abcd_normalize', 'zpk2ss', 'ss2zpk', 'lti',
           'TransferFunction', 'ZerosPolesGain', 'StateSpace', 'lsim',
           'lsim2', 'impulse', 'impulse2', 'step', 'step2', 'bode',
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
    """
    Linear Time Invariant system base class.

    Parameters
    ----------
    *system : arguments
        The `lti` class can be instantiated with either 2, 3 or 4 arguments.
        The following gives the number of arguments and the corresponding
        subclass that is created:

            * 2: `TransferFunction`:  (numerator, denominator)
            * 3: `ZerosPolesGain`: (zeros, poles, gain)
            * 4: `StateSpace`:  (A, B, C, D)

        Each argument can be an array or a sequence.

    Notes
    -----
    `lti` instances do not exist directly. Instead, `lti` creates an instance
    of one of its subclasses: `StateSpace`, `TransferFunction` or
    `ZerosPolesGain`.

    Changing the value of properties that are not directly part of the current
    system representation (such as the `zeros` of a `StateSpace` system) is
    very inefficient and may lead to numerical inaccuracies.

    """
    def __new__(cls, *system):
        """Create an instance of the appropriate subclass."""
        if cls is lti:
            N = len(system)
            if N == 2:
                return super(lti, cls).__new__(TransferFunction)
            elif N == 3:
                return super(lti, cls).__new__(ZerosPolesGain)
            elif N == 4:
                return super(lti, cls).__new__(StateSpace)
            else:
                raise ValueError('Needs 2, 3 or 4 arguments.')
        # __new__ was called from a subclass, let it call its own functions
        return super(lti, cls).__new__(cls)

    def __init__(self, *system):
        """
        Initialize the `lti` baseclass.

        The heavy lifting is done by the subclasses.
        """
        self.inputs = None
        self.outputs = None

    @property
    def num(self):
        """Numerator of the `TransferFunction` system."""
        return self.to_tf().num

    @num.setter
    def num(self, num):
        obj = self.to_tf()
        obj.num = num
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def den(self):
        """Denominator of the `TransferFunction` system."""
        return self.to_tf().den

    @den.setter
    def den(self, den):
        obj = self.to_tf()
        obj.den = den
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def zeros(self):
        """Zeros of the `ZerosPolesGain` system."""
        return self.to_zpk().zeros

    @zeros.setter
    def zeros(self, zeros):
        obj = self.to_zpk()
        obj.zeros = zeros
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def poles(self):
        """Poles of the `ZerosPolesGain` system."""
        return self.to_zpk().poles

    @poles.setter
    def poles(self, poles):
        obj = self.to_zpk()
        obj.poles = poles
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def gain(self):
        """Gain of the `ZerosPolesGain` system."""
        return self.to_zpk().gain

    @gain.setter
    def gain(self, gain):
        obj = self.to_zpk()
        obj.gain = gain
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def A(self):
        """A matrix of the `StateSpace` system."""
        return self.to_ss().A

    @A.setter
    def A(self, A):
        obj = self.to_ss()
        obj.A = A
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def B(self):
        """B matrix of the `StateSpace` system."""
        return self.to_ss().B

    @B.setter
    def B(self, B):
        obj = self.to_ss()
        obj.B = B
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def C(self):
        """C matrix of the `StateSpace` system."""
        return self.to_ss().C

    @C.setter
    def C(self, C):
        obj = self.to_ss()
        obj.C = C
        source_class = type(self)
        self._copy(source_class(obj))

    @property
    def D(self):
        """D matrix of the `StateSpace` system."""
        return self.to_ss().D

    @D.setter
    def D(self, D):
        obj = self.to_ss()
        obj.D = D
        source_class = type(self)
        self._copy(source_class(obj))

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


class TransferFunction(lti):
    """Linear Time Invariant system class in transfer function form.

    Represents the system as the transfer function
    :math:`H(s)=\sum_i b[i] s^i / \sum_j a[j] s^i`, where :math:`a` are
    elements of the numerator `num` and :math:`b` are the elements of the
    denominator `den`.

    Parameters
    ----------
    *system : arguments
        The `TransferFunction` class can be instantiated with 1 or 2 arguments.
        The following gives the number of input arguments and their
        interpretation:

            * 1: `lti` system: (`StateSpace`, `TransferFunction` or
              `ZerosPolesGain`)
            * 2: array_like: (numerator, denominator)

    Notes
    -----
    Changing the value of properties that are not part of the
    `TransferFunction` system representation (such as the `A`, `B`, `C`, `D`
    state-space matrices) is very inefficient and may lead to numerical
    inaccuracies.

    """
    def __new__(cls, *system):
        """Handle object conversion if input is an instance of lti."""
        if len(system) == 1 and isinstance(system[0], lti):
            return system[0].to_tf()

        # No special conversion needed
        return super(TransferFunction, cls).__new__(cls)

    def __init__(self, *system):
        """Initialize the state space LTI system."""
        # Conversion of lti instances is handled in __new__
        if isinstance(system[0], lti):
            return

        super(TransferFunction, self).__init__(self, *system)

        self._num = None
        self._den = None

        self.num, self.den = normalize(*system)

    def __repr__(self):
        """Return representation of the system's transfer function"""
        return '{0}(\n{1},\n{2}\n)'.format(
            self.__class__.__name__,
            repr(self.num),
            repr(self.den),
            )

    @property
    def num(self):
        return self._num

    @num.setter
    def num(self, num):
        self._num = atleast_1d(num)

        # Update dimensions
        if len(self.num.shape) > 1:
            self.outputs, self.inputs = self.num.shape
        else:
            self.outputs = 1
            self.inputs = 1

    @property
    def den(self):
        return self._den

    @den.setter
    def den(self, den):
        self._den = atleast_1d(den)

    def _copy(self, system):
        """
        Copy the parameters of another `TransferFunction` object

        Parameters
        ----------
        system : `TransferFunction`
            The `StateSpace` system that is to be copied

        """
        self.num = system.num
        self.den = system.den

    def to_tf(self):
        """
        Return a copy of the current `TransferFunction` system.

        Returns
        -------
        sys : instance of `TransferFunction`
            The current system (copy)

        """
        return copy.deepcopy(self)

    def to_zpk(self):
        """
        Convert system representation to `ZerosPolesGain`.

        Returns
        -------
        sys : instance of `ZerosPolesGain`
            Zeros, poles, gain representation of the current system

        """
        return ZerosPolesGain(*tf2zpk(self.num, self.den))

    def to_ss(self):
        """
        Convert system representation to `StateSpace`.

        Returns
        -------
        sys : instance of `StateSpace`
            State space model of the current system

        """
        return StateSpace(*tf2ss(self.num, self.den))


class ZerosPolesGain(lti):
    """
    Linear Time Invariant system class in zeros, poles, gain form.

    Represents the system as the transfer function
    :math:`H(s)=k \prod_i (s - z[i]) / \prod_j (s - p[j])`, where :math:`k` is
    the `gain`, :math:`z` are the `zeros` and :math:`p` are the `poles`.

    Parameters
    ----------
    *system : arguments
        The `ZerosPolesGain` class can be instantiated with 1 or 3 arguments.
        The following gives the number of input arguments and their
        interpretation:

            * 1: `lti` system: (`StateSpace`, `TransferFunction` or
              `ZerosPolesGain`)
            * 3: array_like: (zeros, poles, gain)

    Notes
    -----
    Changing the value of properties that are not part of the
    `ZerosPolesGain` system representation (such as the `A`, `B`, `C`, `D`
    state-space matrices) is very inefficient and may lead to numerical
    inaccuracies.

    """
    def __new__(cls, *system):
        """Handle object conversion if input is an instance of `lti`"""
        if len(system) == 1 and isinstance(system[0], lti):
            return system[0].to_zpk()

        # No special conversion needed
        return super(ZerosPolesGain, cls).__new__(cls)

    def __init__(self, *system):
        """Initialize the zeros, poles, gain LTI system."""
        # Conversion of lti instances is handled in __new__
        if isinstance(system[0], lti):
            return

        super(ZerosPolesGain, self).__init__(self, *system)

        self._zeros = None
        self._poles = None
        self._gain = None

        self.zeros, self.poles, self.gain = system

    def __repr__(self):
        """Return representation of the `ZerosPolesGain` system"""
        return '{0}(\n{1},\n{2},\n{3}\n)'.format(
            self.__class__.__name__,
            repr(self.zeros),
            repr(self.poles),
            repr(self.gain),
            )

    @property
    def zeros(self):
        return self._zeros

    @zeros.setter
    def zeros(self, zeros):
        self._zeros = atleast_1d(zeros)

        # Update dimensions
        if len(self.zeros.shape) > 1:
            self.outputs, self.inputs = self.zeros.shape
        else:
            self.outputs = 1
            self.inputs = 1

    @property
    def poles(self):
        return self._poles

    @poles.setter
    def poles(self, poles):
        self._poles = atleast_1d(poles)

    @property
    def gain(self):
        return self._gain

    @gain.setter
    def gain(self, gain):
        self._gain = gain

    def _copy(self, system):
        """
        Copy the parameters of another `ZerosPolesGain` system.

        Parameters
        ----------
        system : instance of `ZerosPolesGain`
            The zeros, poles gain system that is to be copied

        """
        self.poles = system.poles
        self.zeros = system.zeros
        self.gain = system.gain

    def to_tf(self):
        """
        Convert system representation to `TransferFunction`.

        Returns
        -------
        sys : instance of `TransferFunction`
            Transfer function of the current system

        """
        return TransferFunction(*zpk2tf(self.zeros, self.poles, self.gain))

    def to_zpk(self):
        """
        Return a copy of the current 'ZerosPolesGain' system.

        Returns
        -------
        sys : instance of `ZerosPolesGain`
            The current system (copy)

        """
        return copy.deepcopy(self)

    def to_ss(self):
        """
        Convert system representation to `StateSpace`.

        Returns
        -------
        sys : instance of `StateSpace`
            State space model of the current system

        """
        return StateSpace(*zpk2ss(self.zeros, self.poles, self.gain))


class StateSpace(lti):
    """
    Linear Time Invariant system class in state-space form.

    Represents the system as the first order differential equation
    :math:`\dot{x} = A x + B u`.

    Parameters
    ----------
    *system : arguments
        The `StateSpace` class can be instantiated with 1 or 4 arguments.
        The following gives the number of input arguments and their
        interpretation:

            * 1: `lti` system: (`StateSpace`, `TransferFunction` or
              `ZerosPolesGain`)
            * 4: array_like: (A, B, C, D)

    Notes
    -----
    Changing the value of properties that are not part of the
    `StateSpace` system representation (such as `zeros` or `poles`) is very
    inefficient and may lead to numerical inaccuracies.

    """
    def __new__(cls, *system):
        """Handle object conversion if input is an instance of `lti`"""
        if len(system) == 1 and isinstance(system[0], lti):
            return system[0].to_ss()

        # No special conversion needed
        return super(StateSpace, cls).__new__(cls)

    def __init__(self, *system):
        """Initialize the state space LTI system."""
        # Conversion of lti instances is handled in __new__
        if isinstance(system[0], lti):
            return

        super(StateSpace, self).__init__(self, *system)

        self._A = None
        self._B = None
        self._C = None
        self._D = None

        self.A, self.B, self.C, self.D = abcd_normalize(*system)

    def __repr__(self):
        """Return representation of the `StateSpace` system."""
        return '{0}(\n{1},\n{2},\n{3},\n{4}\n)'.format(
            self.__class__.__name__,
            repr(self.A),
            repr(self.B),
            repr(self.C),
            repr(self.D),
            )

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, A):
        self._A = _atleast_2d_or_none(A)

    @property
    def B(self):
        return self._B

    @B.setter
    def B(self, B):
        self._B = _atleast_2d_or_none(B)
        self.inputs = self.B.shape[-1]

    @property
    def C(self):
        return self._C

    @C.setter
    def C(self, C):
        self._C = _atleast_2d_or_none(C)
        self.outputs = self.C.shape[0]

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, D):
        self._D = _atleast_2d_or_none(D)

    def _copy(self, system):
        """
        Copy the parameters of another `StateSpace` system.

        Parameters
        ----------
        system : instance of `StateSpace`
            The state-space system that is to be copied

        """
        self.A = system.A
        self.B = system.B
        self.C = system.C
        self.D = system.D

    def to_tf(self, **kwargs):
        """
        Convert system representation to `TransferFunction`.

        Parameters
        ----------
        kwargs : dict, optional
            Additional keywords passed to `ss2zpk`

        Returns
        -------
        sys : instance of `TransferFunction`
            Transfer function of the current system

        """
        return TransferFunction(*ss2tf(self._A, self._B, self._C, self._D,
                                       **kwargs))

    def to_zpk(self, **kwargs):
        """
        Convert system representation to `ZerosPolesGain`.

        Parameters
        ----------
        kwargs : dict, optional
            Additional keywords passed to `ss2zpk`

        Returns
        -------
        sys : instance of `ZerosPolesGain`
            Zeros, poles, gain representation of the current system

        """
        return ZerosPolesGain(*ss2zpk(self._A, self._B, self._C, self._D,
                                      **kwargs))

    def to_ss(self):
        """
        Return a copy of the current `StateSpace` system.

        Returns
        -------
        sys : instance of `StateSpace`
            The current system (copy)

        """
        return copy.deepcopy(self)


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
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()

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


def lsim(system, U, T, X0=None, interp=True):
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
        represents an input.  If U = 0 or None, a zero input is used.
    T : array_like
        The time steps at which the input is defined and at which the
        output is desired.  Must be nonnegative, increasing, and equally spaced.
    X0 : array_like, optional
        The initial conditions on the state vector (zero by default).
    interp : bool, optional
        Whether to use linear (True, the default) or zero-order-hold (False)
        interpolation for the input array.

    Returns
    -------
    T : 1D ndarray
        Time values for the output.
    yout : 1D ndarray
        System response.
    xout : ndarray
        Time evolution of the state vector.

    Examples
    --------
    Simulate a double integrator y'' = u, with a constant input u = 1

    >>> from scipy import signal
    >>> system = signal.lti([[0., 1.], [0., 0.]], [[0.], [1.]], [[1., 0.]], 0.)
    >>> t = np.linspace(0, 5)
    >>> u = np.ones_like(t)
    >>> tout, y, x = signal.lsim(system, u, t)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t, y)
    """
    if isinstance(system, lti):
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()
    T = atleast_1d(T)
    if len(T.shape) != 1:
        raise ValueError("T must be a rank-1 array.")

    A, B, C, D = map(np.asarray, (sys.A, sys.B, sys.C, sys.D))
    n_states = A.shape[0]
    n_inputs = B.shape[1]

    n_steps = T.size
    if X0 is None:
        X0 = zeros(n_states, sys.A.dtype)
    xout = zeros((n_steps, n_states), sys.A.dtype)

    if T[0] == 0:
        xout[0] = X0
    elif T[0] > 0:
        # step forward to initial time, with zero input
        xout[0] = dot(X0, linalg.expm(transpose(A) * T[0]))
    else:
        raise ValueError("Initial time must be nonnegative")

    no_input = (U is None
                or (isinstance(U, (int, float)) and U == 0.)
                or not np.any(U))

    if n_steps == 1:
        yout = squeeze(dot(xout, transpose(C)))
        if not no_input:
            yout += squeeze(dot(U, transpose(D)))
        return T, squeeze(yout), squeeze(xout)

    dt = T[1] - T[0]
    if not np.allclose((T[1:] - T[:-1]) / dt, 1.0):
        warnings.warn("Non-uniform timesteps are deprecated. Results may be "
                      "slow and/or inaccurate.", DeprecationWarning)
        return lsim2(system, U, T, X0)

    if no_input:
        # Zero input: just use matrix exponential
        # take transpose because state is a row vector
        expAT_dt = linalg.expm(transpose(A) * dt)
        for i in xrange(1, n_steps):
            xout[i] = dot(xout[i-1], expAT_dt)
        yout = squeeze(dot(xout, transpose(C)))
        return T, squeeze(yout), squeeze(xout)

    # Nonzero input
    U = atleast_1d(U)
    if U.ndim == 1:
        U = U[:, np.newaxis]

    if U.shape[0] != n_steps:
        raise ValueError("U must have the same number of rows "
                         "as elements in T.")

    if U.shape[1] != n_inputs:
        raise ValueError("System does not define that many inputs.")

    if not interp:
        # Zero-order hold
        # Algorithm: to integrate from time 0 to time dt, we solve
        #   xdot = A x + B u,  x(0) = x0
        #   udot = 0,          u(0) = u0.
        #
        # Solution is
        #   [ x(dt) ]       [ A*dt   B*dt ] [ x0 ]
        #   [ u(dt) ] = exp [  0     0    ] [ u0 ]
        M = np.vstack([np.hstack([A * dt, B * dt]),
                       np.zeros((n_inputs, n_states + n_inputs))])
        # transpose everything because the state and input are row vectors
        expMT = linalg.expm(transpose(M))
        Ad = expMT[:n_states, :n_states]
        Bd = expMT[n_states:, :n_states]
        for i in xrange(1, n_steps):
            xout[i] = dot(xout[i-1], Ad) + dot(U[i-1], Bd)
    else:
        # Linear interpolation between steps
        # Algorithm: to integrate from time 0 to time dt, with linear
        # interpolation between inputs u(0) = u0 and u(dt) = u1, we solve
        #   xdot = A x + B u,        x(0) = x0
        #   udot = (u1 - u0) / dt,   u(0) = u0.
        #
        # Solution is
        #   [ x(dt) ]       [ A*dt  B*dt  0 ] [  x0   ]
        #   [ u(dt) ] = exp [  0     0    I ] [  u0   ]
        #   [u1 - u0]       [  0     0    0 ] [u1 - u0]
        M = np.vstack([np.hstack([A * dt, B * dt,
                                  np.zeros((n_states, n_inputs))]),
                       np.hstack([np.zeros((n_inputs, n_states + n_inputs)),
                                  np.identity(n_inputs)]),
                       np.zeros((n_inputs, n_states + 2 * n_inputs))])
        expMT = linalg.expm(transpose(M))
        Ad = expMT[:n_states, :n_states]
        Bd1 = expMT[n_states+n_inputs:, :n_states]
        Bd0 = expMT[n_states:n_states + n_inputs, :n_states] - Bd1
        for i in xrange(1, n_steps):
            xout[i] = (dot(xout[i-1], Ad) + dot(U[i-1], Bd0) + dot(U[i], Bd1))

    yout = (squeeze(dot(xout, transpose(C))) + squeeze(dot(U, transpose(D))))
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
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()
    if X0 is None:
        X = squeeze(sys.B)
    else:
        X = squeeze(sys.B + X0)
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)
    else:
        T = asarray(T)

    _, h, _ = lsim(sys, 0., T, X, interp=False)
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
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()
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
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()
    if N is None:
        N = 100
    if T is None:
        T = _default_response_times(sys.A, N)
    else:
        T = asarray(T)
    U = ones(T.shape, sys.A.dtype)
    vals = lsim(sys, U, T, X0=X0, interp=False)
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
        sys = system.to_ss()
    else:
        sys = lti(*system).to_ss()
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
        sys = system.to_tf()
    else:
        sys = lti(*system).to_tf()

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


# This class will be used by place_poles to return its results
# see http://code.activestate.com/recipes/52308/
class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def _valid_inputs(A, B, poles, method, rtol, maxiter):
    """
    Check the poles come in complex conjugage pairs
    Check shapes of A, B and poles are compatible.
    Check the method chosen is compatible with provided poles
    Return update method to use and ordered poles

    """
    poles = np.asarray(poles)
    if poles.ndim > 1:
        raise ValueError("Poles must be a 1D array like.")
    # Will raise ValueError if poles do not come in complex conjugates pairs
    poles = _order_complex_poles(poles)
    if A.ndim > 2:
        raise ValueError("A must be a 2D array/matrix.")
    if B.ndim > 2:
        raise ValueError("B must be a 2D array/matrix")
    if A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if len(poles) > A.shape[0]:
        raise ValueError("maximum number of poles is %d but you asked for %d" %
                         (A.shape[0], len(poles)))
    if len(poles) < A.shape[0]:
        raise ValueError("number of poles is %d but you should provide %d" %
                         (len(poles), A.shape[0]))
    r = np.linalg.matrix_rank(B)
    for p in poles:
        if sum(p == poles) > r:
            raise ValueError("at least one of the requested pole is repeated "
                             "more than rank(B) times")
    # Choose update method
    update_loop = _YT_loop
    if method not in ('KNV0','YT'):
        raise ValueError("The method keyword must be one of 'YT' or 'KNV0'")

    if method == "KNV0":
        update_loop = _KNV0_loop
        if not all(np.isreal(poles)):
            raise ValueError("Complex poles are not supported by KNV0")

    if maxiter < 1:
        raise ValueError("maxiter must be at least equal to 1")

    # We do not check rtol <= 0 as the user can use a negative rtol to
    # force maxiter iterations
    if rtol > 1:
        raise ValueError("rtol can not be greater than 1")

    return update_loop, poles


def _order_complex_poles(poles):
    """
    Check we have complex conjugates pairs and reorder P according to YT, ie
    real_poles, complex_i, conjugate complex_i, ....
    The lexicographic sort on the complex poles is added to help the user to
    compare sets of poles.
    """
    ordered_poles = np.sort(poles[np.isreal(poles)])
    im_poles = []
    for p in np.sort(poles[np.imag(poles) < 0]):
        if np.conj(p) in poles:
            im_poles.extend((p, np.conj(p)))

    ordered_poles = np.hstack((ordered_poles, im_poles))

    if poles.shape[0] != len(ordered_poles):
        raise ValueError("Complex poles must come with their conjugates")
    return ordered_poles


def _KNV0(B, ker_pole, transfer_matrix, j, poles):
    """
    Algorithm "KNV0" Kautsky et Al. Robust pole
    assignment in linear state feedback, Int journal of Control
    1985, vol 41 p 1129->1155
    http://la.epfl.ch/files/content/sites/la/files/
        users/105941/public/KautskyNicholsDooren

    """
    # Remove xj form the base
    transfer_matrix_not_j = np.delete(transfer_matrix, j, axis=1)
    # If we QR this matrix in full mode Q=Q0|Q1
    # then Q1 will be a single column orthogonnal to
    # Q0, that's what we are looking for !

    # After merge of gh-4249 great speed improvements could be achieved
    # using QR updates instead of full QR in the line below

    # To debug with numpy qr uncomment the line below
    # Q, R = np.linalg.qr(transfer_matrix_not_j, mode="complete")
    Q, R = s_qr(transfer_matrix_not_j, mode="full")

    mat_ker_pj = np.dot(ker_pole[j], ker_pole[j].T)
    yj = np.dot(mat_ker_pj, Q[:, -1])

    # If Q[:, -1] is "almost" orthogonal to ker_pole[j] its
    # projection into ker_pole[j] will yield a vector
    # close to 0.  As we are looking for a vector in ker_pole[j]
    # simply stick with transfer_matrix[:, j] (unless someone provides me with
    # a better choice ?)

    if not np.allclose(yj, 0):
        xj = yj/np.linalg.norm(yj)
        transfer_matrix[:, j] = xj

        # KNV does not support complex poles, using YT technique the two lines
        # below seem to work 9 out of 10 times but it is not reliable enough:
        # transfer_matrix[:, j]=real(xj)
        # transfer_matrix[:, j+1]=imag(xj)

        # Add this at the beginning of this function if you wish to test
        # complex support:
        #    if ~np.isreal(P[j]) and (j>=B.shape[0]-1 or P[j]!=np.conj(P[j+1])):
        #        return
        # Problems arise when imag(xj)=>0 I have no idea on how to fix this


def _YT_real(ker_pole, Q, transfer_matrix, i, j):
    """
    Applies algorithm from YT section 6.1 page 19 related to real pairs
    """
    # step 1 page 19
    u = Q[:, -2, np.newaxis]
    v = Q[:, -1, np.newaxis]

    # step 2 page 19
    m = np.dot(np.dot(ker_pole[i].T, np.dot(u, v.T) -
        np.dot(v, u.T)), ker_pole[j])

    # step 3 page 19
    um, sm, vm = np.linalg.svd(m)
    # mu1, mu2 two first columns of U => 2 first lines of U.T
    mu1, mu2 = um.T[:2, :, np.newaxis]
    # VM is V.T with numpy we want the first two lines of V.T
    nu1, nu2 = vm[:2, :, np.newaxis]

    # what follows is a rough python translation of the formulas
    # in section 6.2 page 20 (step 4)
    transfer_matrix_j_mo_transfer_matrix_j = np.vstack((
            transfer_matrix[:, i, np.newaxis],
            transfer_matrix[:, j, np.newaxis]))

    if not np.allclose(sm[0], sm[1]):
        ker_pole_imo_mu1 = np.dot(ker_pole[i], mu1)
        ker_pole_i_nu1 = np.dot(ker_pole[j], nu1)
        ker_pole_mu_nu = np.vstack((ker_pole_imo_mu1, ker_pole_i_nu1))
    else:
        ker_pole_ij = np.vstack((
                                np.hstack((ker_pole[i],
                                           np.zeros(ker_pole[i].shape))),
                                np.hstack((np.zeros(ker_pole[j].shape),
                                                    ker_pole[j]))
                                ))
        mu_nu_matrix = np.vstack(
            (np.hstack((mu1, mu2)), np.hstack((nu1, nu2)))
            )
        ker_pole_mu_nu = np.dot(ker_pole_ij, mu_nu_matrix)
    transfer_matrix_ij = np.dot(np.dot(ker_pole_mu_nu, ker_pole_mu_nu.T),
                             transfer_matrix_j_mo_transfer_matrix_j)
    if not np.allclose(transfer_matrix_ij, 0):
        transfer_matrix_ij = (np.sqrt(2)*transfer_matrix_ij /
                              np.linalg.norm(transfer_matrix_ij))
        transfer_matrix[:, i] = transfer_matrix_ij[
            :transfer_matrix[:, i].shape[0], 0
            ]
        transfer_matrix[:, j] = transfer_matrix_ij[
            transfer_matrix[:, i].shape[0]:, 0
            ]
    else:
        # As in knv0 if transfer_matrix_j_mo_transfer_matrix_j is orthogonal to
        # Vect{ker_pole_mu_nu} assign transfer_matrixi/transfer_matrix_j to
        # ker_pole_mu_nu and iterate. As we are looking for a vector in
        # Vect{Matker_pole_MU_NU} (see section 6.1 page 19) this might help
        # (that's a guess, not a claim !)
        transfer_matrix[:, i] = ker_pole_mu_nu[
            :transfer_matrix[:, i].shape[0], 0
            ]
        transfer_matrix[:, j] = ker_pole_mu_nu[
            transfer_matrix[:, i].shape[0]:, 0
            ]


def _YT_complex(ker_pole, Q, transfer_matrix, i, j):
    """
    Applies algorithm from YT section 6.2 page 20 related to complex pairs
    """
    # step 1 page 20
    ur = np.sqrt(2)*Q[:, -2, np.newaxis]
    ui = np.sqrt(2)*Q[:, -1, np.newaxis]
    u = ur + 1j*ui

    # step 2 page 20
    ker_pole_ij = ker_pole[i]
    m = np.dot(np.dot(np.conj(ker_pole_ij.T), np.dot(u, np.conj(u).T) -
               np.dot(np.conj(u), u.T)), ker_pole_ij)

    # step 3 page 20
    e_val, e_vec = np.linalg.eig(m)
    # sort eigenvalues according to their module
    e_val_idx = np.argsort(np.abs(e_val))
    mu1 = e_vec[:, e_val_idx[-1], np.newaxis]
    mu2 = e_vec[:, e_val_idx[-2], np.newaxis]

    # what follows is a rough python translation of the formulas
    # in section 6.2 page 20 (step 4)

    # remember transfer_matrix_i has been split as
    # transfer_matrix[i]=real(transfer_matrix_i) and
    # transfer_matrix[j]=imag(transfer_matrix_i)
    transfer_matrix_j_mo_transfer_matrix_j = (
        transfer_matrix[:, i, np.newaxis] +
        1j*transfer_matrix[:, j, np.newaxis]
        )
    if not np.allclose(np.abs(e_val[e_val_idx[-1]]),
                              np.abs(e_val[e_val_idx[-2]])):
        ker_pole_mu = np.dot(ker_pole_ij, mu1)
    else:
        mu1_mu2_matrix = np.hstack((mu1, mu2))
        ker_pole_mu = np.dot(ker_pole_ij, mu1_mu2_matrix)
    transfer_matrix_i_j = np.dot(np.dot(ker_pole_mu, np.conj(ker_pole_mu.T)),
                              transfer_matrix_j_mo_transfer_matrix_j)

    if not np.allclose(transfer_matrix_i_j, 0):
        transfer_matrix_i_j = (transfer_matrix_i_j /
            np.linalg.norm(transfer_matrix_i_j))
        transfer_matrix[:, i] = np.real(transfer_matrix_i_j[:, 0])
        transfer_matrix[:, j] = np.imag(transfer_matrix_i_j[:, 0])
    else:
        # same idea as in YT_real
        transfer_matrix[:, i] = np.real(ker_pole_mu[:, 0])
        transfer_matrix[:, j] = np.imag(ker_pole_mu[:, 0])


def _YT_loop(ker_pole, transfer_matrix, poles, B, maxiter, rtol):
    """
    Algorithm "YT" Tits, Yang. Globally Convergent
    Algorithms for Robust Pole Assignment by State Feedback
    http://drum.lib.umd.edu/handle/1903/5598
    The poles P have to be sorted accordingly to section 6.2 page 20

    """
    # The IEEE edition of the YT paper gives useful information on the
    # optimal update order for the real poles in order to minimize the number
    # of times we have to loop over all poles, see page 1442
    nb_real = poles[np.isreal(poles)].shape[0]
    # hnb => Half Nb Real
    hnb = nb_real // 2

    # Stick to the indices in the paper and then remove one to get numpy array
    # index it is a bit easier to link the code to the paper this way even if it
    # is not very clean. The paper is unclear about what should be done when
    # there is only one real pole => use KNV0 on this real pole seem to work
    if nb_real > 0:
        #update the biggest real pole with the smallest one
        update_order = [[nb_real], [1]]
    else:
        update_order = [[],[]]

    r_comp = np.arange(nb_real+1, len(poles)+1, 2)
    # step 1.a
    r_p = np.arange(1, hnb+nb_real % 2)
    update_order[0].extend(2*r_p)
    update_order[1].extend(2*r_p+1)
    # step 1.b
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 1.c
    r_p = np.arange(1, hnb+1)
    update_order[0].extend(2*r_p-1)
    update_order[1].extend(2*r_p)
    # step 1.d
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 2.a
    r_j = np.arange(2, hnb+nb_real % 2)
    for j in r_j:
        for i in range(1, hnb+1):
            update_order[0].append(i)
            update_order[1].append(i+j)
    # step 2.b
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 2.c
    r_j = np.arange(2, hnb+nb_real % 2)
    for j in r_j:
        for i in range(hnb+1, nb_real+1):
            idx_1 = i+j
            if idx_1 > nb_real:
                idx_1 = i+j-nb_real
            update_order[0].append(i)
            update_order[1].append(idx_1)
    # step 2.d
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 3.a
    for i in range(1, hnb+1):
        update_order[0].append(i)
        update_order[1].append(i+hnb)
    # step 3.b
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)

    update_order = np.array(update_order).T-1
    stop = False
    nb_try = 0
    while nb_try < maxiter and not stop:
        det_transfer_matrixb = np.abs(np.linalg.det(transfer_matrix))
        for i, j in update_order:
            if i == j:
                assert i == 0, "i!=0 for KNV call in YT"
                assert np.isreal(poles[i]), "calling KNV on a complex pole"
                _KNV0(B, ker_pole, transfer_matrix, i, poles)
            else:
                transfer_matrix_not_i_j = np.delete(transfer_matrix, (i, j),
                                                    axis=1)
                # after merge of gh-4249 great speed improvements could be
                # achieved using QR updates instead of full QR in the line below

                #to debug with numpy qr uncomment the line below
                #Q, _ = np.linalg.qr(transfer_matrix_not_i_j, mode="complete")
                Q, _ = s_qr(transfer_matrix_not_i_j, mode="full")

                if np.isreal(poles[i]):
                    assert np.isreal(poles[j]), "mixing real and complex " + \
                        "in YT_real" + str(poles)
                    _YT_real(ker_pole, Q, transfer_matrix, i, j)
                else:
                    assert ~np.isreal(poles[i]), "mixing real and complex " + \
                        "in YT_real" + str(poles)
                    _YT_complex(ker_pole, Q, transfer_matrix, i, j)

        det_transfer_matrix = np.max((np.sqrt(np.spacing(1)),
                                  np.abs(np.linalg.det(transfer_matrix))))
        cur_rtol = np.abs(
            (det_transfer_matrix -
             det_transfer_matrixb) /
            det_transfer_matrix)
        if cur_rtol < rtol and det_transfer_matrix > np.sqrt(np.spacing(1)):
            # Convergence test from YT page 21
            stop = True
        nb_try += 1
    return stop, cur_rtol, nb_try


def _KNV0_loop(ker_pole, transfer_matrix, poles, B, maxiter, rtol):
    """
    Loop over all poles one by one and apply KNV method 0 algorithm
    """
    # This method is useful only because we need to be able to call
    # _KNV0 from YT without looping over all poles, otherwise it would
    # have been fine to mix _KNV0_loop and _KNV0 in a single function
    stop = False
    nb_try = 0
    while nb_try < maxiter and not stop:
        det_transfer_matrixb = np.abs(np.linalg.det(transfer_matrix))
        for j in range(B.shape[0]):
            _KNV0(B, ker_pole, transfer_matrix, j, poles)

        det_transfer_matrix = np.max((np.sqrt(np.spacing(1)),
                                  np.abs(np.linalg.det(transfer_matrix))))
        cur_rtol = np.abs((det_transfer_matrix - det_transfer_matrixb) /
                       det_transfer_matrix)
        if cur_rtol < rtol and det_transfer_matrix > np.sqrt(np.spacing(1)):
            # Convergence test from YT page 21
            stop = True

        nb_try += 1
    return stop, cur_rtol, nb_try


def place_poles(A, B, poles, method="YT", rtol=1e-3, maxiter=30):
    """
    Compute K such that eigenvalues (A - dot(B, K))=poles.

    K is the gain matrix such as the plant described by the linear system
    ``AX+BU`` will have its closed-loop poles, i.e the eigenvalues ``A - B*K``,
    as close as possible to those asked for in poles.

    SISO, MISO and MIMO systems are supported.

    Parameters
    ----------
    A, B : ndarray
        State-space representation of linear system ``AX + BU``.
    poles : array_like
        Desired real poles and/or complex conjugates poles.
        Complex poles are only supported with ``method="YT"`` (default).
    method: {'YT', 'KNV0'}, optional
        Which method to choose to find the gain matrix K. One of:

            - 'YT': Yang Tits
            - 'KNV0': Kautsky, Nichols, Van Dooren update method 0

        See References and Notes for details on the algorithms.
    rtol: float, optional
        After each iteration the determinant of the eigenvectors of
        ``A - B*K`` is compared to its previous value, when the relative
        error between these two values becomes lower than `rtol` the algorithm
        stops.  Default is 1e-3.
    maxiter: int, optional
        Maximum number of iterations to compute the gain matrix.
        Default is 30.

    Returns
    -------
    full_state_feedback : Bunch object
        full_state_feedback is composed of:
            gain_matrix : 1-D ndarray
                The closed loop matrix K such as the eigenvalues of ``A-BK``
                are as close as possible to the requested poles.
            computed_poles : 1-D ndarray
                The poles corresponding to ``A-BK`` sorted as first the real
                poles in increasing order, then the complex congugates in
                lexicographic order.
            requested_poles : 1-D ndarray
                The poles the algorithm was asked to place sorted as above,
                they may differ from what was achieved.
            X : 2D ndarray
                The transfer matrix such as ``X * diag(poles) = (A - B*K)*X``
                (see Notes)
            rtol : float
                The relative tolerance achieved on ``det(X)`` (see Notes).
                `rtol` will be NaN if the optimisation algorithms can not run,
                i.e when ``B.shape[1] == 1``, or 0 when the solution is unique.
            nb_iter : int
                The number of iterations performed before converging.
                `nb_iter` will be NaN if the optimisation algorithms can
                not run, i.e when ``B.shape[1] == 1``, or 0 when the solution
                is unique.

    Notes
    -----
    The Tits and Yang (YT), [2]_ paper is an update of the original Kautsky et
    al. (KNV) paper [1]_.  KNV relies on rank-1 updates to find the transfer
    matrix X such that ``X * diag(poles) = (A - B*K)*X``, whereas YT uses
    rank-2 updates. This yields on average more robust solutions (see [2]_
    pp 21-22), furthermore the YT algorithm supports complex poles whereas KNV
    does not in its original version.  Only update method 0 proposed by KNV has
    been implemented here, hence the name ``'KNV0'``.

    KNV extended to complex poles is used in Matlab's ``place`` function, YT is
    distributed under a non-free licence by Slicot under the name ``robpole``.
    It is unclear and undocumented how KNV0 has been extended to complex poles
    (Tits and Yang claim on page 14 of their paper that their method can not be
    used to extend KNV to complex poles), therefore only YT supports them in
    this implementation.

    As the solution to the problem of pole placement is not unique for MIMO
    systems, both methods start with a tentative transfer matrix which is
    altered in various way to increase its determinant.  Both methods have been
    proven to converge to a stable solution, however depending on the way the
    initial transfer matrix is chosen they will converge to different
    solutions and therefore there is absolutely no guarantee that using
    ``'KNV0'`` will yield results similar to Matlab's or any other
    implementation of these algorithms.

    Using the default method ``'YT'`` should be fine in most cases; ``'KNV0'``
    is only provided because it is needed by ``'YT'`` in some specific cases.
    Furthermore ``'YT'`` gives on average more robust results than ``'KNV0'``
    when ``abs(det(X))`` is used as a robustness indicator.

    [2]_ is available as a technical report on the following URL:
    http://drum.lib.umd.edu/handle/1903/5598

    References
    ----------
    .. [1] J. Kautsky, N.K. Nichols and P. van Dooren, "Robust pole assignment
           in linear state feedback", International Journal of Control, Vol. 41
           pp. 1129-1155, 1985.
    .. [2] A.L. Tits and Y. Yang, "Globally convergent algorithms for robust
           pole assignment by state feedback, IEEE Transactions on Automatic
           Control, Vol. 41, pp. 1432-1452, 1996.

    Examples
    --------
    A simple example demonstrating real pole placement using both KNV and YT
    algorithms.  This is example number 1 from section 4 of the reference KNV
    publication ([1]_):

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> A = np.array([[ 1.380,  -0.2077,  6.715, -5.676  ],
    ...               [-0.5814, -4.290,   0,      0.6750 ],
    ...               [ 1.067,   4.273,  -6.654,  5.893  ],
    ...               [ 0.0480,  4.273,   1.343, -2.104  ]])
    >>> B = np.array([[ 0,      5.679 ],
    ...               [ 1.136,  1.136 ],
    ...               [ 0,      0,    ],
    ...               [-3.146,  0     ]])
    >>> P = np.array([-0.2, -0.5, -5.0566, -8.6659])

    Now compute K with KNV method 0, with the default YT method and with the YT
    method while forcing 100 iterations of the algorithm and print some results
    after each call.

    >>> fsf1 = signal.place_poles(A, B, P, method='KNV0')
    >>> fsf1.gain_matrix
    array([[ 0.20071427, -0.96665799,  0.24066128, -0.10279785],
           [ 0.50587268,  0.57779091,  0.51795763, -0.41991442]])

    >>> fsf2 = signal.place_poles(A, B, P)  # uses YT method
    >>> fsf2.computed_poles
    array([-8.6659, -5.0566, -0.5   , -0.2   ])

    >>> fsf3 = signal.place_poles(A, B, P, rtol=-1, maxiter=100)
    >>> fsf3.X
    array([[ 0.52072442+0.j, -0.08409372+0.j, -0.56847937+0.j,  0.74823657+0.j],
           [-0.04977751+0.j, -0.80872954+0.j,  0.13566234+0.j, -0.29322906+0.j],
           [-0.82266932+0.j, -0.19168026+0.j, -0.56348322+0.j, -0.43815060+0.j],
           [ 0.22267347+0.j,  0.54967577+0.j, -0.58387806+0.j, -0.40271926+0.j]])

    The absolute value of the determinant of X is a good indicator to check the
    robustness of the results, both ``'KNV0'`` and ``'YT'`` aim at maximizing
    it.  Below a comparison of the robustness of the results above:

    >>> abs(np.linalg.det(fsf1.X)) < abs(np.linalg.det(fsf2.X))
    True
    >>> abs(np.linalg.det(fsf2.X)) < abs(np.linalg.det(fsf3.X))
    True

    Now a simple example for complex poles:

    >>> A = np.array([[ 0,  7/3.,  0,   0   ],
    ...               [ 0,   0,    0,  7/9. ],
    ...               [ 0,   0,    0,   0   ],
    ...               [ 0,   0,    0,   0   ]])
    >>> B = np.array([[ 0,  0 ],
    ...               [ 0,  0 ],
    ...               [ 1,  0 ],
    ...               [ 0,  1 ]])
    >>> P = np.array([-3, -1, -2-1j, -2+1j]) / 3.
    >>> fsf = signal.place_poles(A, B, P, method='YT')

    We can plot the desired and computed poles in the complex plane:

    >>> t = np.linspace(0, 2*np.pi, 401)
    >>> plt.plot(np.cos(t), np.sin(t), 'k--')  # unit circle
    >>> plt.plot(fsf.requested_poles.real, fsf.requested_poles.imag,
    ...          'wo', label='Desired')
    >>> plt.plot(fsf.computed_poles.real, fsf.computed_poles.imag, 'bx',
    ...          label='Placed')
    >>> plt.grid()
    >>> plt.axis('image')
    >>> plt.axis([-1.1, 1.1, -1.1, 1.1])
    >>> plt.legend(bbox_to_anchor=(1.05, 1), loc=2, numpoints=1)

    """
    # Move away all the inputs checking, it only adds noise to the code
    update_loop, poles = _valid_inputs(A, B, poles, method, rtol, maxiter)

    # The current value of the relative tolerance we achieved
    cur_rtol = np.nan
    # The number of iterations needed before converging
    nb_iter = np.nan

    # Step A: QR decomposition of B page 1132 KN
    # to debug with numpy qr uncomment the line below
    # u, z = np.linalg.qr(B, mode="complete")
    u, z = s_qr(B, mode="full")
    rankB = np.linalg.matrix_rank(B)
    u0 = u[:, :rankB]
    u1 = u[:, rankB:]
    z = z[:rankB, :]

    # If the solution is unique
    if B.shape[0] == rankB:
        # if B is square and full rank there is only one solution
        # such as (A+BK)=diag(P) i.e BK=diag(P)-A
        # if B has as many lines as its rank (but not square) the solution
        # is the same as above using least squares
        # => use lstsq in both cases
        # for complex poles we use the following trick
        #
        # |a -b| has for eigenvalues a+b and a-b
        # |b a|
        #
        # |a+bi 0| has the obvious eigenvalues a+bi and a-bi
        # |0 a-bi|
        #
        # e.g solving the first one in R gives the solution
        # for the second one in C
        diag_poles = np.zeros(A.shape)
        idx = 0
        while idx < poles.shape[0]:
            p = poles[idx]
            diag_poles[idx, idx] = np.real(p)
            if ~np.isreal(p):
                diag_poles[idx, idx+1] = -np.imag(p)
                diag_poles[idx+1, idx+1] = np.real(p)
                diag_poles[idx+1, idx] = np.imag(p)
                idx += 1  # skip next one
            idx += 1
        gain_matrix = np.linalg.lstsq(B, diag_poles-A)[0]
        transfer_matrix = np.eye(A.shape[0])
        cur_rtol = 0
        nb_iter = 0
    else:
        # step A (p1144 KNV) and begining of step F: decompose
        # dot(U1.T, A-P[i]*I).T and build our set of transfer_matrix vectors
        # in the same loop
        ker_pole = []

        # flag to skip the conjugate of a complex pole
        skip_conjugate = False
        # select orthonormal base ker_pole for each Pole and vectors for
        # transfer_matrix
        for j in range(B.shape[0]):
            if skip_conjugate:
                skip_conjugate = False
                continue
            pole_space_j = np.dot(u1.T, A-poles[j]*np.eye(B.shape[0])).T

            # after QR Q=Q0|Q1
            # only Q0 is used to reconstruct  the qr'ed (dot Q, R) matrix.
            # Q1 is orthogonnal to Q0 and will be multiplied by the zeros in
            # R when using mode "complete". In default mode Q1 and the zeros
            # in R are not computed

            # To debug with numpy qr uncomment the line below
            # Q, _ = np.linalg.qr(pole_space_j, mode="complete")
            Q, _ = s_qr(pole_space_j, mode="full")

            ker_pole_j = Q[:, pole_space_j.shape[1]:]

            # We want to select one vector in ker_pole_j to build the transfer
            # matrix, however qr returns sometimes vectors with zeros on the same
            # line for each pole and this yields very long convergence times.
            # Or some other times a set of vectors, one with zero imaginary
            # part and one (or several) with imaginary parts. After trying
            # many ways to select the best possible one (eg ditch vectors
            # with zero imaginary part for complex poles) I ended up summing
            # all vectors in ker_pole_j, this solves 100% of the problems and is
            # still a valid choice for transfer_matrix. Indeed for complex poles
            # we are sure to have a non zero imaginary part that way, and the
            # problem of lines full of zeros in transfer_matrix is solved too as
            # when a vector from ker_pole_j has a zero the other one(s)
            # (when ker_pole_j.shape[1]>1) for sure won't have a zero there.
            transfer_matrix_j = np.sum(ker_pole_j, axis=1)[:, np.newaxis]
            transfer_matrix_j = (transfer_matrix_j /
                                 np.linalg.norm(transfer_matrix_j))
            if ~np.isreal(poles[j]):  # complex pole
                transfer_matrix_j = np.hstack([np.real(transfer_matrix_j),
                                               np.imag(transfer_matrix_j)])
                ker_pole.extend([ker_pole_j, ker_pole_j])

                # Skip next pole as it is the conjugate
                skip_conjugate = True
            else:  # real pole, nothing to do
                ker_pole.append(ker_pole_j)

            if j == 0:
                transfer_matrix = transfer_matrix_j
            else:
                transfer_matrix = np.hstack((transfer_matrix, transfer_matrix_j))

        if rankB > 1:  # otherwise there is nothing we can optimize
            stop, cur_rtol, nb_iter = update_loop(ker_pole, transfer_matrix,
                                                  poles, B, maxiter, rtol)
            if not stop and rtol > 0:
                # if rtol<=0 the user has probably done that on purpose,
                # don't annoy him
                err_msg = (
                    "Convergence was not reached after maxiter iterations.\n"
                    "You asked for a relative tolerance of %f we got %f" %
                    (rtol, cur_rtol)
                    )
                warnings.warn(err_msg)

        # reconstruct transfer_matrix to match complex conjugate pairs,
        # ie transfer_matrix_j/transfer_matrix_j+1 are
        # Re(Complex_pole), Im(Complex_pole) now and will be Re-Im/Re+Im after
        transfer_matrix = transfer_matrix.astype(complex)
        idx = 0
        while idx < poles.shape[0]-1:
            if ~np.isreal(poles[idx]):
                rel = transfer_matrix[:, idx].copy()
                img = transfer_matrix[:, idx+1]
                # rel will be an array referencing a column of transfer_matrix
                # if we don't copy() it will changer after the next line and
                # and the line after will not yield the correct value
                transfer_matrix[:, idx] = rel-1j*img
                transfer_matrix[:, idx+1] = rel+1j*img
                idx += 1  # skip next one
            idx += 1

        try:
            m = np.linalg.solve(transfer_matrix.T, np.dot(np.diag(poles),
                                                          transfer_matrix.T)).T
            gain_matrix = np.linalg.solve(z, np.dot(u0.T, m-A))
        except np.linalg.LinAlgError:
            raise ValueError("The poles you've chosen can't be placed. "
                             "Check the controllability matrix and try "
                             "another set of poles")

    # Beware: Kautsky solves A+BK but the usual form is A-BK
    gain_matrix = -gain_matrix
    # K still contains complex with ~=0j imaginary parts, get rid of them
    gain_matrix = np.real(gain_matrix)

    full_state_feedback = Bunch()
    full_state_feedback.gain_matrix = gain_matrix
    full_state_feedback.computed_poles = _order_complex_poles(
        np.linalg.eig(A - np.dot(B, gain_matrix))[0]
        )
    full_state_feedback.requested_poles = poles
    full_state_feedback.X = transfer_matrix
    full_state_feedback.rtol = cur_rtol
    full_state_feedback.nb_iter = nb_iter

    return full_state_feedback
