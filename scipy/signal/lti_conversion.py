"""
ltisys -- a collection of functions to convert linear time invariant systems
from one representation to another.
"""
from __future__ import division, print_function, absolute_import

import warnings

import numpy
import numpy as np
from numpy import (r_, eye, atleast_2d, poly, dot,
                   asarray, product, zeros, array, outer)
from scipy import linalg


__all__ = ['BadCoefficients', 'abcd_normalize', 'cont2discrete', 'normalize',
           'sos2tf', 'sos2zpk', 'ss2tf', 'ss2zpk', 'tf2sos', 'tf2ss', 'tf2zpk',
           'zpk2sos', 'zpk2ss', 'zpk2tf']


class BadCoefficients(UserWarning):
    """Warning about badly conditioned filter coefficients"""
    pass


def tf2ss(num, den):
    r"""Transfer function to state-space representation.

    Parameters
    ----------
    num, den : array_like
        Sequences representing the coefficients of the numerator and
        denominator polynomials, in order of descending degree. The
        denominator needs to be at least as long as the numerator.

    Returns
    -------
    A, B, C, D : ndarray
        State space representation of the system, in controller canonical
        form.

    Examples
    --------
    Convert the transfer function:

    .. math:: H(s) = \frac{s^2 + 3s + 3}{s^2 + 2s + 1}

    >>> num = [1, 3, 3]
    >>> den = [1, 2, 1]

    to the state-space representation:

    .. math::

        \dot{\textbf{x}}(t) =
        \begin{bmatrix} -2 & -1 \\ 1 & 0 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \\ 0 \end{bmatrix} \textbf{u}(t) \\

        \textbf{y}(t) = \begin{bmatrix} 1 & 2 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \end{bmatrix} \textbf{u}(t)

    >>> from scipy.signal import tf2ss
    >>> A, B, C, D = tf2ss(num, den)
    >>> A
    array([[-2., -1.],
           [ 1.,  0.]])
    >>> B
    array([[ 1.],
           [ 0.]])
    >>> C
    array([[ 1.,  2.]])
    >>> D
    array([[ 1.]])
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
        D = atleast_2d(num[:, 0])

    else:
        # We don't assign it an empty array because this system
        # is not 'null'. It just doesn't have a non-zero D
        # matrix. Thus, it should have a non-zero shape so that
        # it can be operated on by functions like 'ss2tf'
        D = array([[0]], float)

    if K == 1:
        D = D.reshape(num.shape)

        return (zeros((1, 1)), zeros((1, D.shape[1])),
                zeros((D.shape[0], 1)), D)

    frow = -array([den[1:]])
    A = r_[frow, eye(K - 2, K - 1)]
    B = eye(K - 1, 1)
    C = num[:, 1:] - outer(num[:, 0], den[1:])
    D = D.reshape((C.shape[0], B.shape[1]))

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
        See `ss2tf` for format.

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
    r"""State-space to transfer function.

    A, B, C, D defines a linear state-space system with `p` inputs,
    `q` outputs, and `n` state variables.

    Parameters
    ----------
    A : array_like
        State (or system) matrix of shape ``(n, n)``
    B : array_like
        Input matrix of shape ``(n, p)``
    C : array_like
        Output matrix of shape ``(q, n)``
    D : array_like
        Feedthrough (or feedforward) matrix of shape ``(q, p)``
    input : int, optional
        For multiple-input systems, the index of the input to use.

    Returns
    -------
    num : 2-D ndarray
        Numerator(s) of the resulting transfer function(s).  `num` has one row
        for each of the system's outputs. Each row is a sequence representation
        of the numerator polynomial.
    den : 1-D ndarray
        Denominator of the resulting transfer function(s).  `den` is a sequence
        representation of the denominator polynomial.

    Examples
    --------
    Convert the state-space representation:

    .. math::

        \dot{\textbf{x}}(t) =
        \begin{bmatrix} -2 & -1 \\ 1 & 0 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \\ 0 \end{bmatrix} \textbf{u}(t) \\

        \textbf{y}(t) = \begin{bmatrix} 1 & 2 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \end{bmatrix} \textbf{u}(t)

    >>> A = [[-2, -1], [1, 0]]
    >>> B = [[1], [0]]  # 2-dimensional column vector
    >>> C = [[1, 2]]    # 2-dimensional row vector
    >>> D = 1

    to the transfer function:

    .. math:: H(s) = \frac{s^2 + 3s + 3}{s^2 + 2s + 1}

    >>> from scipy.signal import ss2tf
    >>> ss2tf(A, B, C, D)
    (array([[1, 3, 3]]), array([ 1.,  2.,  1.]))
    """
    # transfer function is C (sI - A)**(-1) B + D

    # Check consistency and make them all rank-2 arrays
    A, B, C, D = abcd_normalize(A, B, C, D)

    nout, nin = D.shape
    if input >= nin:
        raise ValueError("System does not have the input specified.")

    # make SIMO from possibly MIMO system.
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

    A, B, C, D defines a linear state-space system with `p` inputs,
    `q` outputs, and `n` state variables.

    Parameters
    ----------
    A : array_like
        State (or system) matrix of shape ``(n, n)``
    B : array_like
        Input matrix of shape ``(n, p)``
    C : array_like
        Output matrix of shape ``(q, n)``
    D : array_like
        Feedthrough (or feedforward) matrix of shape ``(q, p)``
    input : int, optional
        For multiple-input systems, the index of the input to use.

    Returns
    -------
    z, p : sequence
        Zeros and poles.
    k : float
        System gain.

    """
    return tf2zpk(*ss2tf(A, B, C, D, input=input))


def tf2zpk(b, a):
    r"""Return zero, pole, gain (z, p, k) representation from a numerator,
    denominator representation of a linear filter.

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.

    Returns
    -------
    z : ndarray
        Zeros of the transfer function.
    p : ndarray
        Poles of the transfer function.
    k : float
        System gain.

    Notes
    -----
    If some values of `b` are too close to 0, they are removed. In that case,
    a BadCoefficients warning is emitted.

    The `b` and `a` arrays are interpreted as coefficients for positive,
    descending powers of the transfer function variable.  So the inputs
    :math:`b = [b_0, b_1, ..., b_M]` and :math:`a =[a_0, a_1, ..., a_N]`
    can represent an analog filter of the form:

    .. math::

        H(s) = \frac
        {b_0 s^M + b_1 s^{(M-1)} + \cdots + b_M}
        {a_0 s^N + a_1 s^{(N-1)} + \cdots + a_N}

    or a discrete-time filter of the form:

    .. math::

        H(z) = \frac
        {b_0 z^M + b_1 z^{(M-1)} + \cdots + b_M}
        {a_0 z^N + a_1 z^{(N-1)} + \cdots + a_N}

    This "positive powers" form is found more commonly in controls
    engineering.  If `M` and `N` are equal (which is true for all filters
    generated by the bilinear transform), then this happens to be equivalent
    to the "negative powers" discrete-time form preferred in DSP:

    .. math::

        H(z) = \frac
        {b_0 + b_1 z^{-1} + \cdots + b_M z^{-M}}
        {a_0 + a_1 z^{-1} + \cdots + a_N z^{-N}}

    Although this is true for common filters, remember that this is not true
    in the general case.  If `M` and `N` are not equal, the discrete-time
    transfer function coefficients must first be converted to the "positive
    powers" form before finding the poles and zeros.

    """
    b, a = normalize(b, a)
    b = (b + 0.0) / a[0]
    a = (a + 0.0) / a[0]
    k = b[0]
    b /= b[0]
    z = np.roots(b)
    p = np.roots(a)
    return z, p, k


def zpk2tf(z, p, k):
    """
    Return polynomial transfer function representation from zeros and poles

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    """
    z = np.atleast_1d(z)
    k = np.atleast_1d(k)
    if len(z.shape) > 1:
        temp = poly(z[0])
        b = zeros((z.shape[0], z.shape[1] + 1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]] * z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * poly(z[i])
    else:
        b = k * poly(z)
    a = np.atleast_1d(poly(p))

    # Use real output if possible.  Copied from numpy.poly, since
    # we can't depend on a specific version of numpy.
    if issubclass(b.dtype.type, numpy.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = numpy.asarray(z, complex)
        pos_roots = numpy.compress(roots.imag > 0, roots)
        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if numpy.all(numpy.sort_complex(neg_roots) ==
                         numpy.sort_complex(pos_roots)):
                b = b.real.copy()

    if issubclass(a.dtype.type, numpy.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = numpy.asarray(p, complex)
        pos_roots = numpy.compress(roots.imag > 0, roots)
        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if numpy.all(numpy.sort_complex(neg_roots) ==
                         numpy.sort_complex(pos_roots)):
                a = a.real.copy()

    return b, a


def tf2sos(b, a, pairing='nearest'):
    """
    Return second-order sections from transfer function representation

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.
    pairing : {'nearest', 'keep_odd'}, optional
        The method to use to combine pairs of poles and zeros into sections.
        See `zpk2sos`.

    Returns
    -------
    sos : ndarray
        Array of second-order filter coefficients, with shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    See Also
    --------
    zpk2sos, sosfilt

    Notes
    -----
    It is generally discouraged to convert from TF to SOS format, since doing
    so usually will not improve numerical precision errors. Instead, consider
    designing filters in ZPK format and converting directly to SOS. TF is
    converted to SOS by first converting to ZPK format, then converting
    ZPK to SOS.

    .. versionadded:: 0.16.0
    """
    return zpk2sos(*tf2zpk(b, a), pairing=pairing)


def sos2tf(sos):
    """
    Return a single transfer function from a series of second-order sections

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    Notes
    -----
    .. versionadded:: 0.16.0
    """
    sos = np.asarray(sos)
    b = [1.]
    a = [1.]
    n_sections = sos.shape[0]
    for section in range(n_sections):
        b = np.polymul(b, sos[section, :3])
        a = np.polymul(a, sos[section, 3:])
    return b, a


def sos2zpk(sos):
    """
    Return zeros, poles, and gain of a series of second-order sections

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    z : ndarray
        Zeros of the transfer function.
    p : ndarray
        Poles of the transfer function.
    k : float
        System gain.

    Notes
    -----
    .. versionadded:: 0.16.0
    """
    sos = np.asarray(sos)
    n_sections = sos.shape[0]
    z = np.empty(n_sections*2, np.complex128)
    p = np.empty(n_sections*2, np.complex128)
    k = 1.
    for section in range(n_sections):
        zpk = tf2zpk(sos[section, :3], sos[section, 3:])
        z[2*section:2*(section+1)] = zpk[0]
        p[2*section:2*(section+1)] = zpk[1]
        k *= zpk[2]
    return z, p, k


def zpk2sos(z, p, k, pairing='nearest'):
    """
    Return second-order sections from zeros, poles, and gain of a system

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.
    pairing : {'nearest', 'keep_odd'}, optional
        The method to use to combine pairs of poles and zeros into sections.
        See Notes below.

    Returns
    -------
    sos : ndarray
        Array of second-order filter coefficients, with shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    See Also
    --------
    sosfilt

    Notes
    -----
    The algorithm used to convert ZPK to SOS format is designed to
    minimize errors due to numerical precision issues. The pairing
    algorithm attempts to minimize the peak gain of each biquadratic
    section. This is done by pairing poles with the nearest zeros, starting
    with the poles closest to the unit circle.

    *Algorithms*

    The current algorithms are designed specifically for use with digital
    filters. (The output coefficents are not correct for analog filters.)

    The steps in the ``pairing='nearest'`` and ``pairing='keep_odd'``
    algorithms are mostly shared. The ``nearest`` algorithm attempts to
    minimize the peak gain, while ``'keep_odd'`` minimizes peak gain under
    the constraint that odd-order systems should retain one section
    as first order. The algorithm steps and are as follows:

    As a pre-processing step, add poles or zeros to the origin as
    necessary to obtain the same number of poles and zeros for pairing.
    If ``pairing == 'nearest'`` and there are an odd number of poles,
    add an additional pole and a zero at the origin.

    The following steps are then iterated over until no more poles or
    zeros remain:

    1. Take the (next remaining) pole (complex or real) closest to the
       unit circle to begin a new filter section.

    2. If the pole is real and there are no other remaining real poles [#]_,
       add the closest real zero to the section and leave it as a first
       order section. Note that after this step we are guaranteed to be
       left with an even number of real poles, complex poles, real zeros,
       and complex zeros for subsequent pairing iterations.

    3. Else:

        1. If the pole is complex and the zero is the only remaining real
           zero*, then pair the pole with the *next* closest zero
           (guaranteed to be complex). This is necessary to ensure that
           there will be a real zero remaining to eventually create a
           first-order section (thus keeping the odd order).

        2. Else pair the pole with the closest remaining zero (complex or
           real).

        3. Proceed to complete the second-order section by adding another
           pole and zero to the current pole and zero in the section:

            1. If the current pole and zero are both complex, add their
               conjugates.

            2. Else if the pole is complex and the zero is real, add the
               conjugate pole and the next closest real zero.

            3. Else if the pole is real and the zero is complex, add the
               conjugate zero and the real pole closest to those zeros.

            4. Else (we must have a real pole and real zero) add the next
               real pole closest to the unit circle, and then add the real
               zero closest to that pole.

    .. [#] This conditional can only be met for specific odd-order inputs
           with the ``pairing == 'keep_odd'`` method.

    .. versionadded:: 0.16.0

    Examples
    --------

    Design a 6th order low-pass elliptic digital filter for a system with a
    sampling rate of 8000 Hz that has a pass-band corner frequency of
    1000 Hz.  The ripple in the pass-band should not exceed 0.087 dB, and
    the attenuation in the stop-band should be at least 90 dB.

    In the following call to `signal.ellip`, we could use ``output='sos'``,
    but for this example, we'll use ``output='zpk'``, and then convert to SOS
    format with `zpk2sos`:

    >>> from scipy import signal
    >>> z, p, k = signal.ellip(6, 0.087, 90, 1000/(0.5*8000), output='zpk')

    Now convert to SOS format.

    >>> sos = signal.zpk2sos(z, p, k)

    The coefficients of the numerators of the sections:

    >>> sos[:, :3]
    array([[ 0.0014154 ,  0.00248707,  0.0014154 ],
           [ 1.        ,  0.72965193,  1.        ],
           [ 1.        ,  0.17594966,  1.        ]])

    The symmetry in the coefficients occurs because all the zeros are on the
    unit circle.

    The coefficients of the denominators of the sections:

    >>> sos[:, 3:]
    array([[ 1.        , -1.32543251,  0.46989499],
           [ 1.        , -1.26117915,  0.6262586 ],
           [ 1.        , -1.25707217,  0.86199667]])

    The next example shows the effect of the `pairing` option.  We have a
    system with three poles and three zeros, so the SOS array will have
    shape (2, 6).  The means there is, in effect, an extra pole and an extra
    zero at the origin in the SOS representation.

    >>> z1 = np.array([-1, -0.5-0.5j, -0.5+0.5j])
    >>> p1 = np.array([0.75, 0.8+0.1j, 0.8-0.1j])

    With ``pairing='nearest'`` (the default), we obtain

    >>> signal.zpk2sos(z1, p1, 1)
    array([[ 1.  ,  1.  ,  0.5 ,  1.  , -0.75,  0.  ],
           [ 1.  ,  1.  ,  0.  ,  1.  , -1.6 ,  0.65]])

    The first section has the zeros {-0.5-0.05j, -0.5+0.5j} and the poles
    {0, 0.75}, and the second section has the zeros {-1, 0} and poles
    {0.8+0.1j, 0.8-0.1j}.  Note that the extra pole and zero at the origin
    have been assigned to different sections.

    With ``pairing='keep_odd'``, we obtain:

    >>> signal.zpk2sos(z1, p1, 1, pairing='keep_odd')
    array([[ 1.  ,  1.  ,  0.  ,  1.  , -0.75,  0.  ],
           [ 1.  ,  1.  ,  0.5 ,  1.  , -1.6 ,  0.65]])

    The extra pole and zero at the origin are in the same section.
    The first section is, in effect, a first-order section.

    """
    # TODO in the near future:
    # 1. Add SOS capability to `filtfilt`, `freqz`, etc. somehow (#3259).
    # 2. Make `decimate` use `sosfilt` instead of `lfilter`.
    # 3. Make sosfilt automatically simplify sections to first order
    #    when possible. Note this might make `sosfiltfilt` a bit harder (ICs).
    # 4. Further optimizations of the section ordering / pole-zero pairing.
    # See the wiki for other potential issues.

    valid_pairings = ['nearest', 'keep_odd']
    if pairing not in valid_pairings:
        raise ValueError('pairing must be one of %s, not %s'
                         % (valid_pairings, pairing))
    if len(z) == len(p) == 0:
        return array([[k, 0., 0., 1., 0., 0.]])

    # ensure we have the same number of poles and zeros, and make copies
    p = np.concatenate((p, np.zeros(max(len(z) - len(p), 0))))
    z = np.concatenate((z, np.zeros(max(len(p) - len(z), 0))))
    n_sections = (max(len(p), len(z)) + 1) // 2
    sos = zeros((n_sections, 6))

    if len(p) % 2 == 1 and pairing == 'nearest':
        p = np.concatenate((p, [0.]))
        z = np.concatenate((z, [0.]))
    assert len(p) == len(z)

    # Ensure we have complex conjugate pairs
    # (note that _cplxreal only gives us one element of each complex pair):
    z = np.concatenate(_cplxreal(z))
    p = np.concatenate(_cplxreal(p))

    p_sos = np.zeros((n_sections, 2), np.complex128)
    z_sos = np.zeros_like(p_sos)
    for si in range(n_sections):
        # Select the next "worst" pole
        p1_idx = np.argmin(np.abs(1 - np.abs(p)))
        p1 = p[p1_idx]
        p = np.delete(p, p1_idx)

        # Pair that pole with a zero

        if np.isreal(p1) and np.isreal(p).sum() == 0:
            # Special case to set a first-order section
            z1_idx = _nearest_real_complex_idx(z, p1, 'real')
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)
            p2 = z2 = 0
        else:
            if not np.isreal(p1) and np.isreal(z).sum() == 1:
                # Special case to ensure we choose a complex zero to pair
                # with so later (setting up a first-order section)
                z1_idx = _nearest_real_complex_idx(z, p1, 'complex')
                assert not np.isreal(z[z1_idx])
            else:
                # Pair the pole with the closest zero (real or complex)
                z1_idx = np.argmin(np.abs(p1 - z))
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)

            # Now that we have p1 and z1, figure out what p2 and z2 need to be
            if not np.isreal(p1):
                if not np.isreal(z1):  # complex pole, complex zero
                    p2 = p1.conj()
                    z2 = z1.conj()
                else:  # complex pole, real zero
                    p2 = p1.conj()
                    z2_idx = _nearest_real_complex_idx(z, p1, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
            else:
                if not np.isreal(z1):  # real pole, complex zero
                    z2 = z1.conj()
                    p2_idx = _nearest_real_complex_idx(p, z1, 'real')
                    p2 = p[p2_idx]
                    assert np.isreal(p2)
                else:  # real pole, real zero
                    # pick the next "worst" pole to use
                    idx = np.where(np.isreal(p))[0]
                    assert len(idx) > 0
                    p2_idx = idx[np.argmin(np.abs(np.abs(p[idx]) - 1))]
                    p2 = p[p2_idx]
                    # find a real zero to match the added pole
                    assert np.isreal(p2)
                    z2_idx = _nearest_real_complex_idx(z, p2, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
                p = np.delete(p, p2_idx)
        p_sos[si] = [p1, p2]
        z_sos[si] = [z1, z2]
    assert len(p) == len(z) == 0  # we've consumed all poles and zeros
    del p, z

    # Construct the system, reversing order so the "worst" are last
    p_sos = np.reshape(p_sos[::-1], (n_sections, 2))
    z_sos = np.reshape(z_sos[::-1], (n_sections, 2))
    gains = np.ones(n_sections)
    gains[0] = k
    for si in range(n_sections):
        x = zpk2tf(z_sos[si], p_sos[si], gains[si])
        sos[si] = np.concatenate(x)
    return sos


def _cplxreal(z, tol=None):
    """
    Split into complex and real parts, combining conjugate pairs.

    The 1D input vector `z` is split up into its complex (`zc`) and real (`zr`)
    elements.  Every complex element must be part of a complex-conjugate pair,
    which are combined into a single number (with positive imaginary part) in
    the output.  Two complex numbers are considered a conjugate pair if their
    real and imaginary parts differ in magnitude by less than ``tol * abs(z)``.

    Parameters
    ----------
    z : array_like
        Vector of complex numbers to be sorted and split
    tol : float, optional
        Relative tolerance for testing realness and conjugate equality.
        Default is ``100 * spacing(1)`` of `z`'s data type (i.e. 2e-14 for
        float64)

    Returns
    -------
    zc : ndarray
        Complex elements of `z`, with each pair represented by a single value
        having positive imaginary part, sorted first by real part, and then
        by magnitude of imaginary part.  The pairs are averaged when combined
        to reduce error.
    zr : ndarray
        Real elements of `z` (those having imaginary part less than
        `tol` times their magnitude), sorted by value.

    Raises
    ------
    ValueError
        If there are any complex numbers in `z` for which a conjugate
        cannot be found.

    See Also
    --------
    _cplxpair

    Examples
    --------
    >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
    >>> zc, zr = _cplxreal(a)
    >>> print zc
    [ 1.+1.j  2.+1.j  2.+1.j  2.+2.j]
    >>> print zr
    [ 1.  3.  4.]
    """

    z = np.atleast_1d(z)
    if z.size == 0:
        return z, z
    elif z.ndim != 1:
        raise ValueError('_cplxreal only accepts 1D input')

    if tol is None:
        # Get tolerance from dtype of input
        tol = 100 * np.finfo((1.0 * z).dtype).eps

    # Sort by real part, magnitude of imaginary part (speed up further sorting)
    z = z[np.lexsort((abs(z.imag), z.real))]

    # Split reals from conjugate pairs
    real_indices = abs(z.imag) <= tol * abs(z)
    zr = z[real_indices].real

    if len(zr) == len(z):
        # Input is entirely real
        return array([]), zr

    # Split positive and negative halves of conjugates
    z = z[~real_indices]
    zp = z[z.imag > 0]
    zn = z[z.imag < 0]

    if len(zp) != len(zn):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Find runs of (approximately) the same real part
    same_real = np.diff(zp.real) <= tol * abs(zp[:-1])
    diffs = numpy.diff(np.concatenate(([0], same_real, [0])))
    run_starts = numpy.where(diffs > 0)[0]
    run_stops = numpy.where(diffs < 0)[0]

    # Sort each run by their imaginary parts
    for i in range(len(run_starts)):
        start = run_starts[i]
        stop = run_stops[i] + 1
        for chunk in (zp[start:stop], zn[start:stop]):
            chunk[...] = chunk[np.lexsort([abs(chunk.imag)])]

    # Check that negatives match positives
    if any(abs(zp - zn.conj()) > tol * abs(zn)):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Average out numerical inaccuracy in real vs imag parts of pairs
    zc = (zp + zn.conj()) / 2

    return zc, zr


def _cplxpair(z, tol=None):
    """
    Sort into pairs of complex conjugates.

    Complex conjugates in `z` are sorted by increasing real part.  In each
    pair, the number with negative imaginary part appears first.

    If pairs have identical real parts, they are sorted by increasing
    imaginary magnitude.

    Two complex numbers are considered a conjugate pair if their real and
    imaginary parts differ in magnitude by less than ``tol * abs(z)``.  The
    pairs are forced to be exact complex conjugates by averaging the positive
    and negative values.

    Purely real numbers are also sorted, but placed after the complex
    conjugate pairs.  A number is considered real if its imaginary part is
    smaller than `tol` times the magnitude of the number.

    Parameters
    ----------
    z : array_like
        1-dimensional input array to be sorted.
    tol : float, optional
        Relative tolerance for testing realness and conjugate equality.
        Default is ``100 * spacing(1)`` of `z`'s data type (i.e. 2e-14 for
        float64)

    Returns
    -------
    y : ndarray
        Complex conjugate pairs followed by real numbers.

    Raises
    ------
    ValueError
        If there are any complex numbers in `z` for which a conjugate
        cannot be found.

    See Also
    --------
    _cplxreal

    Examples
    --------
    >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
    >>> z = _cplxpair(a)
    >>> print(z)
    [ 1.-1.j  1.+1.j  2.-1.j  2.+1.j  2.-1.j  2.+1.j  2.-2.j  2.+2.j  1.+0.j
      3.+0.j  4.+0.j]
    """

    z = np.atleast_1d(z)
    if z.size == 0 or np.isrealobj(z):
        return np.sort(z)

    if z.ndim != 1:
        raise ValueError('z must be 1-dimensional')

    zc, zr = _cplxreal(z, tol)

    # Interleave complex values and their conjugates, with negative imaginary
    # parts first in each pair
    zc = np.dstack((zc.conj(), zc)).flatten()
    z = np.append(zc, zr)
    return z


def _nearest_real_complex_idx(fro, to, which):
    """Get the next closest real or complex element based on distance"""
    assert which in ('real', 'complex')
    order = np.argsort(np.abs(fro - to))
    mask = np.isreal(fro[order])
    if which == 'complex':
        mask = ~mask
    return order[np.where(mask)[0][0]]


def cont2discrete(system, dt, method="zoh", alpha=None):
    """
    Transform a continuous to a discrete state-space system.

    Parameters
    ----------
    system : a tuple describing the system or an instance of `lti`
        The following gives the number of elements in the tuple and
        the interpretation:

            * 1: (instance of `lti`)
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
    if len(system) == 1:
        return system.to_discrete()
    if len(system) == 2:
        sysd = cont2discrete(tf2ss(system[0], system[1]), dt, method=method,
                             alpha=alpha)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 3:
        sysd = cont2discrete(zpk2ss(system[0], system[1], system[2]), dt,
                             method=method, alpha=alpha)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 4:
        a, b, c, d = system
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
        return cont2discrete(system, dt, method="gbt", alpha=0.5)

    elif method == 'euler' or method == 'forward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=0.0)

    elif method == 'backward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=1.0)

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


def normalize(b, a):
    """Normalize numerator/denominator of a continuous-time transfer function.

    If values of `b` are too close to 0, they are removed. In that case, a
    BadCoefficients warning is emitted.

    Parameters
    ----------
    b: array_like
        Numerator of the transfer function. Can be a 2d array to normalize
        multiple transfer functions.
    a: array_like
        Denominator of the transfer function. At most 1d.

    Returns
    -------
    num: array
        The numerator of the normalized transfer function. At least a 1d
        array. A 2d-array if the input `num` is a 2d array.
    den: 1d-array
        The denominator of the normalized transfer function.

    Notes
    -----
    Coefficients for both the numerator and denominator should be specified in
    descending exponent order (e.g., ``s^2 + 3s + 5`` would be represented as
    ``[1, 3, 5]``).
    """
    num, den = b, a

    den = np.atleast_1d(den)
    num = np.atleast_2d(_align_nums(num))

    if den.ndim != 1:
        raise ValueError("Denominator polynomial must be rank-1 array.")
    if num.ndim > 2:
        raise ValueError("Numerator polynomial must be rank-1 or"
                         " rank-2 array.")
    if np.all(den == 0):
        raise ValueError("Denominator must have at least on nonzero element.")

    # Trim leading zeros in denominator, leave at least one.
    den = np.trim_zeros(den, 'f')

    # Normalize transfer function
    num, den = num / den[0], den / den[0]

    # Count numerator columns that are all zero
    leading_zeros = 0
    for col in num.T:
        if np.allclose(col, 0, atol=1e-14):
            leading_zeros += 1
        else:
            break

    # Trim leading zeros of numerator
    if leading_zeros > 0:
        warnings.warn("Badly conditioned filter coefficients (numerator): the "
                      "results may be meaningless", BadCoefficients)
        # Make sure at least one column remains
        if leading_zeros == num.shape[1]:
            leading_zeros -= 1
        num = num[:, leading_zeros:]

    # Squeeze first dimension if singular
    if num.shape[0] == 1:
        num = num[0, :]

    return num, den


def _align_nums(nums):
    """Aligns the shapes of multiple numerators.

    Given an array of numerator coefficient arrays [[a_1, a_2,...,
    a_n],..., [b_1, b_2,..., b_m]], this function pads shorter numerator
    arrays with zero's so that all numerators have the same length. Such
    alignment is necessary for functions like 'tf2ss', which needs the
    alignment when dealing with SIMO transfer functions.

    Parameters
    ----------
    nums: array_like
        Numerator or list of numerators. Not necessarily with same length.

    Returns
    -------
    nums: array
        The numerator. If `nums` input was a list of numerators then a 2d
        array with padded zeros for shorter numerators is returned. Otherwise
        returns ``np.asarray(nums)``.
    """
    try:
        # The statement can throw a ValueError if one
        # of the numerators is a single digit and another
        # is array-like e.g. if nums = [5, [1, 2, 3]]
        nums = asarray(nums)

        if not np.issubdtype(nums.dtype, np.number):
            raise ValueError("dtype of numerator is non-numeric")

        return nums

    except ValueError:
        nums = [np.atleast_1d(num) for num in nums]
        max_width = max(num.size for num in nums)

        # pre-allocate
        aligned_nums = np.zeros((len(nums), max_width))

        # Create numerators with padded zeros
        for index, num in enumerate(nums):
            aligned_nums[index, -num.size:] = num

        return aligned_nums
