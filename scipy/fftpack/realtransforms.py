"""
Real spectrum tranforms (DCT, DST, MDCT)
"""
from __future__ import division, print_function, absolute_import


__all__ = ['dct', 'idct', 'dst', 'idst']

import numpy as np
from scipy.fftpack import _fftpack
from scipy.fftpack.basic import _datacopied

import atexit
atexit.register(_fftpack.destroy_ddct1_cache)
atexit.register(_fftpack.destroy_ddct2_cache)
atexit.register(_fftpack.destroy_dct1_cache)
atexit.register(_fftpack.destroy_dct2_cache)

atexit.register(_fftpack.destroy_ddst1_cache)
atexit.register(_fftpack.destroy_ddst2_cache)
atexit.register(_fftpack.destroy_dst1_cache)
atexit.register(_fftpack.destroy_dst2_cache)


def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=0):
    """
    Return the Discrete Cosine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3}, optional
        Type of the DCT (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis over which to compute the transform.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    y : ndarray of real
        The transformed input array.

    See Also
    --------
    idct

    Notes
    -----
    For a single dimension array ``x``, ``dct(x, norm='ortho')`` is equal to
    MATLAB ``dct(x)``.

    There are theoretically 8 types of the DCT, only the first 3 types are
    implemented in scipy. 'The' DCT generally refers to DCT type 2, and 'the'
    Inverse DCT generally refers to DCT type 3.

    **type I**

    There are several definitions of the DCT-I; we use the following
    (for ``norm=None``)::

                                         N-2
      y[k] = x[0] + (-1)**k x[N-1] + 2 * sum x[n]*cos(pi*k*n/(N-1))
                                         n=1

    Only None is supported as normalization mode for DCT-I. Note also that the
    DCT-I is only supported for input size > 1

    **type II**

    There are several definitions of the DCT-II; we use the following
    (for ``norm=None``)::


                N-1
      y[k] = 2* sum x[n]*cos(pi*k*(2n+1)/(2*N)), 0 <= k < N.
                n=0

    If ``norm='ortho'``, ``y[k]`` is multiplied by a scaling factor `f`::

      f = sqrt(1/(4*N)) if k = 0,
      f = sqrt(1/(2*N)) otherwise.

    Which makes the corresponding matrix of coefficients orthonormal
    (``OO' = Id``).

    **type III**

    There are several definitions, we use the following
    (for ``norm=None``)::

                        N-1
      y[k] = x[0] + 2 * sum x[n]*cos(pi*(k+0.5)*n/N), 0 <= k < N.
                        n=1

    or, for ``norm='ortho'`` and 0 <= k < N::

                                          N-1
      y[k] = x[0] / sqrt(N) + sqrt(1/N) * sum x[n]*cos(pi*(k+0.5)*n/N)
                                          n=1

    The (unnormalized) DCT-III is the inverse of the (unnormalized) DCT-II, up
    to a factor `2N`. The orthonormalized DCT-III is exactly the inverse of
    the orthonormalized DCT-II.

    References
    ----------
    .. [1] 'A Fast Cosine Transform in One and Two Dimensions', by J. 
           Makhoul, `IEEE Transactions on acoustics, speech and signal 
           processing` vol. 28(1), pp. 27-34, 
           http://dx.doi.org/10.1109/TASSP.1980.1163351 (1980).
    .. [2] Wikipedia, "Discrete cosine transform",
           http://en.wikipedia.org/wiki/Discrete_cosine_transform

    """
    if type == 1 and norm is not None:
        raise NotImplementedError(
              "Orthonormalization not yet supported for DCT-I")
    return _dct(x, type, n, axis, normalize=norm, overwrite_x=overwrite_x)


def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=0):
    """
    Return the Inverse Discrete Cosine Transform of an arbitrary type sequence.

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3}, optional
        Type of the DCT (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis over which to compute the transform.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    idct : ndarray of real
        The transformed input array.

    See Also
    --------
    dct

    Notes
    -----
    For a single dimension array `x`, ``idct(x, norm='ortho')`` is equal to
    MATLAB ``idct(x)``.

    'The' IDCT is the IDCT of type 2, which is the same as DCT of type 3.

    IDCT of type 1 is the DCT of type 1, IDCT of type 2 is the DCT of type
    3, and IDCT of type 3 is the DCT of type 2. For the definition of these
    types, see `dct`.

    """
    if type == 1 and norm is not None:
        raise NotImplementedError(
              "Orthonormalization not yet supported for IDCT-I")
    # Inverse/forward type table
    _TP = {1:1, 2:3, 3:2}
    return _dct(x, _TP[type], n, axis, normalize=norm, overwrite_x=overwrite_x)


def _dct(x, type, n=None, axis=-1, overwrite_x=0, normalize=None):
    """
    Return Discrete Cosine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        input array.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis along which the dct is computed. (default=-1)
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    z : real ndarray

    """
    tmp = np.asarray(x)
    if not np.isrealobj(tmp):
        raise TypeError("1st argument must be real sequence")

    if n is None:
        n = tmp.shape[axis]
    else:
        raise NotImplemented("Padding/truncating not yet implemented")

    if tmp.dtype == np.double:
        if type == 1:
            f = _fftpack.ddct1
        elif type == 2:
            f = _fftpack.ddct2
        elif type == 3:
            f = _fftpack.ddct3
        else:
            raise ValueError("Type %d not understood" % type)
    elif tmp.dtype == np.float32:
        if type == 1:
            f = _fftpack.dct1
        elif type == 2:
            f = _fftpack.dct2
        elif type == 3:
            f = _fftpack.dct3
        else:
            raise ValueError("Type %d not understood" % type)
    else:
        raise ValueError("dtype %s not supported" % tmp.dtype)

    if normalize:
        if normalize == "ortho":
            nm = 1
        else:
            raise ValueError("Unknown normalize mode %s" % normalize)
    else:
        nm = 0

    if type == 1 and n < 2:
        raise ValueError("DCT-I is not defined for size < 2")

    overwrite_x = overwrite_x or _datacopied(tmp, x)

    if axis == -1 or axis == len(tmp.shape) - 1:
        return f(tmp, n, nm, overwrite_x)
    #else:
    #    raise NotImplementedError("Axis arg not yet implemented")

    tmp = np.swapaxes(tmp, axis, -1)
    tmp = f(tmp, n, nm, overwrite_x)
    return np.swapaxes(tmp, axis, -1)


###########

def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=0):
    """
    Return the Discrete Sine Transform of arbitrary type sequence x.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3}, optional
        Type of the DST (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis over which to compute the transform.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    dst : ndarray of reals
        The transformed input array.

    See Also
    --------
    idst

    Notes
    -----
    For a single dimension array ``x``.

    There are theoretically 8 types of the DST for different combinations of
    even/odd boundary conditions and boundary off sets [1]_, only the first
    3 types are implemented in scipy.

    **type I**

    There are several definitions of the DST-I; we use the following
    for ``norm=None``.  DST-I assumes the input is odd around n=-1 and n=N. ::

                 N-1
      y[k] = 2 * sum x[n]*sin(pi*(k+1)*(n+1)/(N+1))
                 n=0

    Only None is supported as normalization mode for DCT-I. Note also that the
    DCT-I is only supported for input size > 1
    The (unnormalized) DCT-I is its own inverse, up to a factor `2(N+1)`.

    **type II**

    There are several definitions of the DST-II; we use the following
    for ``norm=None``.  DST-II assumes the input is odd around n=-1/2 and 
    n=N-1/2; the output is odd around k=-1 and even around k=N-1 ::

                N-1
      y[k] = 2* sum x[n]*sin(pi*(k+1)*(n+0.5)/N), 0 <= k < N.
                n=0

    if ``norm='ortho'``, ``y[k]`` is multiplied by a scaling factor `f` ::

        f = sqrt(1/(4*N)) if k == 0
        f = sqrt(1/(2*N)) otherwise.

    **type III**

    There are several definitions of the DST-III, we use the following
    (for ``norm=None``).  DST-III assumes the input is odd around n=-1
    and even around n=N-1 ::

                                 N-2
      y[k] = x[N-1]*(-1)**k + 2* sum x[n]*sin(pi*(k+0.5)*(n+1)/N), 0 <= k < N.
                                 n=0

    The (unnormalized) DCT-III is the inverse of the (unnormalized) DCT-II, up
    to a factor `2N`.  The orthonormalized DST-III is exactly the inverse of
    the orthonormalized DST-II.

    References
    ----------
    .. [1] Wikipedia, "Discrete sine transform",
           http://en.wikipedia.org/wiki/Discrete_sine_transform

    """
    if type == 1 and norm is not None:
        raise NotImplementedError(
              "Orthonormalization not yet supported for IDCT-I")
    return _dst(x, type, n, axis, normalize=norm, overwrite_x=overwrite_x)


def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=0):
    """
    Return the Inverse Discrete Sine Transform of an arbitrary type sequence.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    x : array_like
        The input array.
    type : {1, 2, 3}, optional
        Type of the DST (see Notes). Default type is 2.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis over which to compute the transform.
    norm : {None, 'ortho'}, optional
        Normalization mode (see Notes). Default is None.
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    idst : ndarray of real
        The transformed input array.

    See Also
    --------
    dst

    Notes
    -----
    'The' IDST is the IDST of type 2, which is the same as DST of type 3.

    IDST of type 1 is the DST of type 1, IDST of type 2 is the DST of type
    3, and IDST of type 3 is the DST of type 2. For the definition of these
    types, see `dst`.

    """
    if type == 1 and norm is not None:
        raise NotImplementedError(
              "Orthonormalization not yet supported for IDCT-I")
    # Inverse/forward type table
    _TP = {1:1, 2:3, 3:2}
    return _dst(x, _TP[type], n, axis, normalize=norm, overwrite_x=overwrite_x)


def _dst(x, type, n=None, axis=-1, overwrite_x=0, normalize=None):
    """
    Return Discrete Sine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        input array.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis along which the dst is computed. (default=-1)
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    z : real ndarray

    """
    tmp = np.asarray(x)
    if not np.isrealobj(tmp):
        raise TypeError("1st argument must be real sequence")

    if n is None:
        n = tmp.shape[axis]
    else:
        raise NotImplemented("Padding/truncating not yet implemented")

    if tmp.dtype == np.double:
        if type == 1:
            f = _fftpack.ddst1
        elif type == 2:
            f = _fftpack.ddst2
        elif type == 3:
            f = _fftpack.ddst3
        else:
            raise ValueError("Type %d not understood" % type)
    elif tmp.dtype == np.float32:
        if type == 1:
            f = _fftpack.dst1
        elif type == 2:
            f = _fftpack.dst2
        elif type == 3:
            f = _fftpack.dst3
        else:
            raise ValueError("Type %d not understood" % type)
    else:
        raise ValueError("dtype %s not supported" % tmp.dtype)

    if normalize:
        if normalize == "ortho":
            nm = 1
        else:
            raise ValueError("Unknown normalize mode %s" % normalize)
    else:
        nm = 0

    if type == 1 and n < 2:
        raise ValueError("DST-I is not defined for size < 2")

    overwrite_x = overwrite_x or _datacopied(tmp, x)

    if axis == -1 or axis == len(tmp.shape) - 1:
        return f(tmp, n, nm, overwrite_x)
    #else:
    #    raise NotImplementedError("Axis arg not yet implemented")

    tmp = np.swapaxes(tmp, axis, -1)
    tmp = f(tmp, n, nm, overwrite_x)
    return np.swapaxes(tmp, axis, -1)
