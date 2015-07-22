"""
=====================================================
Distance computations (:mod:`scipy.spatial.distance`)
=====================================================

.. sectionauthor:: Damian Eads

Function Reference
------------------

Distance matrix computation from a collection of raw observation vectors
stored in a rectangular array.

.. autosummary::
   :toctree: generated/

   pdist   -- pairwise distances between observation vectors.
   cdist   -- distances between two collections of observation vectors
   squareform -- convert distance matrix to a condensed one and vice versa

Predicates for checking the validity of distance matrices, both
condensed and redundant. Also contained in this module are functions
for computing the number of observations in a distance matrix.

.. autosummary::
   :toctree: generated/

   is_valid_dm -- checks for a valid distance matrix
   is_valid_y  -- checks for a valid condensed distance matrix
   num_obs_dm  -- # of observations in a distance matrix
   num_obs_y   -- # of observations in a condensed distance matrix

Distance functions between two vectors ``u`` and ``v``. Computing
distances over a large collection of vectors is inefficient for these
functions. Use ``pdist`` for this purpose.

.. autosummary::
   :toctree: generated/

   braycurtis       -- the Bray-Curtis distance.
   canberra         -- the Canberra distance.
   chebyshev        -- the Chebyshev distance.
   cityblock        -- the Manhattan distance.
   correlation      -- the Correlation distance.
   cosine           -- the Cosine distance.
   dice             -- the Dice dissimilarity (boolean).
   euclidean        -- the Euclidean distance.
   hamming          -- the Hamming distance (boolean).
   jaccard          -- the Jaccard distance (boolean).
   kulsinski        -- the Kulsinski distance (boolean).
   mahalanobis      -- the Mahalanobis distance.
   matching         -- the matching dissimilarity (boolean).
   minkowski        -- the Minkowski distance.
   rogerstanimoto   -- the Rogers-Tanimoto dissimilarity (boolean).
   russellrao       -- the Russell-Rao dissimilarity (boolean).
   seuclidean       -- the normalized Euclidean distance.
   sokalmichener    -- the Sokal-Michener dissimilarity (boolean).
   sokalsneath      -- the Sokal-Sneath dissimilarity (boolean).
   sqeuclidean      -- the squared Euclidean distance.
   wminkowski       -- the weighted Minkowski distance.
   yule             -- the Yule dissimilarity (boolean).

"""

# Copyright (C) Damian Eads, 2007-2008. New BSD License.

from __future__ import division, print_function, absolute_import

import warnings
import numpy as np

from scipy._lib.six import callable, string_types
from scipy._lib.six import xrange

from . import _distance_wrap
from ..linalg import norm
import collections


def _copy_array_if_base_present(a):
    """
    Copies the array if its base points to a parent array.
    """
    if a.base is not None:
        return a.copy()
    elif np.issubsctype(a, np.float32):
        return np.array(a, dtype=np.double)
    else:
        return a


def _copy_arrays_if_base_present(T):
    """
    Accepts a tuple of arrays T. Copies the array T[i] if its base array
    points to an actual array. Otherwise, the reference is just copied.
    This is useful if the arrays are being passed to a C function that
    does not do proper striding.
    """
    l = [_copy_array_if_base_present(a) for a in T]
    return l


def _convert_to_bool(X):
    if X.dtype != bool:
        X = X.astype(bool)
    if not X.flags.contiguous:
        X = X.copy()
    return X


def _convert_to_double(X):
    if X.dtype != np.double:
        X = X.astype(np.double)
    if not X.flags.contiguous:
        X = X.copy()
    return X


def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u


def minkowski(u, v, p):
    """
    Computes the Minkowski distance between two 1-D arrays.

    The Minkowski distance between 1-D arrays `u` and `v`,
    is defined as

    .. math::

       {||u-v||}_p = (\\sum{|u_i - v_i|^p})^{1/p}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    p : int
        The order of the norm of the difference :math:`{||u-v||}_p`.

    Returns
    -------
    d : double
        The Minkowski distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if p < 1:
        raise ValueError("p must be at least 1")
    dist = norm(u - v, ord=p)
    return dist


def wminkowski(u, v, p, w):
    """
    Computes the weighted Minkowski distance between two 1-D arrays.

    The weighted Minkowski distance between `u` and `v`, defined as

    .. math::

       \\left(\\sum{(w_i |u_i - v_i|^p)}\\right)^{1/p}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    p : int
        The order of the norm of the difference :math:`{||u-v||}_p`.
    w : (N,) array_like
        The weight vector.

    Returns
    -------
    wminkowski : double
        The weighted Minkowski distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    w = _validate_vector(w)
    if p < 1:
        raise ValueError("p must be at least 1")
    dist = norm(w * (u - v), ord=p)
    return dist


def euclidean(u, v):
    """
    Computes the Euclidean distance between two 1-D arrays.

    The Euclidean distance between 1-D arrays `u` and `v`, is defined as

    .. math::

       {||u-v||}_2

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    euclidean : double
        The Euclidean distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    dist = norm(u - v)
    return dist


def sqeuclidean(u, v):
    """
    Computes the squared Euclidean distance between two 1-D arrays.

    The squared Euclidean distance between `u` and `v` is defined as

    .. math::

       {||u-v||}_2^2.


    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    sqeuclidean : double
        The squared Euclidean distance between vectors `u` and `v`.

    """
    # Preserve float dtypes, but convert everything else to np.float64
    # for stability.
    utype, vtype = None, None
    if not (hasattr(u, "dtype") and np.issubdtype(u.dtype, np.inexact)):
        utype = np.float64
    if not (hasattr(v, "dtype") and np.issubdtype(v.dtype, np.inexact)):
        vtype = np.float64

    u = _validate_vector(u, dtype=utype)
    v = _validate_vector(v, dtype=vtype)
    u_v = u - v

    return np.dot(u_v, u_v)


def cosine(u, v):
    """
    Computes the Cosine distance between 1-D arrays.

    The Cosine distance between `u` and `v`, is defined as

    .. math::

       1 - \\frac{u \\cdot v}
                {||u||_2 ||v||_2}.

    where :math:`u \\cdot v` is the dot product of :math:`u` and
    :math:`v`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    cosine : double
        The Cosine distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    dist = 1.0 - np.dot(u, v) / (norm(u) * norm(v))
    return dist


def correlation(u, v):
    """
    Computes the correlation distance between two 1-D arrays.

    The correlation distance between `u` and `v`, is
    defined as

    .. math::

       1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
               {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

    where :math:`\\bar{u}` is the mean of the elements of `u`
    and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    correlation : double
        The correlation distance between 1-D array `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    umu = u.mean()
    vmu = v.mean()
    um = u - umu
    vm = v - vmu
    dist = 1.0 - np.dot(um, vm) / (norm(um) * norm(vm))
    return dist


def hamming(u, v):
    """
    Computes the Hamming distance between two 1-D arrays.

    The Hamming distance between 1-D arrays `u` and `v`, is simply the
    proportion of disagreeing components in `u` and `v`. If `u` and `v` are
    boolean vectors, the Hamming distance is

    .. math::

       \\frac{c_{01} + c_{10}}{n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    hamming : double
        The Hamming distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.shape != v.shape:
        raise ValueError('The 1d arrays must have equal lengths.')
    return (u != v).mean()


def jaccard(u, v):
    """
    Computes the Jaccard-Needham dissimilarity between two boolean 1-D arrays.

    The Jaccard-Needham dissimilarity between 1-D boolean arrays `u` and `v`,
    is defined as

    .. math::

       \\frac{c_{TF} + c_{FT}}
            {c_{TT} + c_{FT} + c_{TF}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    jaccard : double
        The Jaccard distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    dist = (np.double(np.bitwise_and((u != v),
                                     np.bitwise_or(u != 0, v != 0)).sum())
            / np.double(np.bitwise_or(u != 0, v != 0).sum()))
    return dist


def kulsinski(u, v):
    """
    Computes the Kulsinski dissimilarity between two boolean 1-D arrays.

    The Kulsinski dissimilarity between two boolean 1-D arrays `u` and `v`,
    is defined as

    .. math::

         \\frac{c_{TF} + c_{FT} - c_{TT} + n}
              {c_{FT} + c_{TF} + n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    kulsinski : double
        The Kulsinski distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    n = float(len(u))
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)

    return (ntf + nft - ntt + n) / (ntf + nft + n)


def seuclidean(u, v, V):
    """
    Returns the standardized Euclidean distance between two 1-D arrays.

    The standardized Euclidean distance between `u` and `v`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    V : (N,) array_like
        `V` is an 1-D array of component variances. It is usually computed
        among a larger collection vectors.

    Returns
    -------
    seuclidean : double
        The standardized Euclidean distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    V = _validate_vector(V, dtype=np.float64)
    if V.shape[0] != u.shape[0] or u.shape[0] != v.shape[0]:
        raise TypeError('V must be a 1-D array of the same dimension '
                        'as u and v.')
    return np.sqrt(((u - v) ** 2 / V).sum())


def cityblock(u, v):
    """
    Computes the City Block (Manhattan) distance.

    Computes the Manhattan distance between two 1-D arrays `u` and `v`,
    which is defined as

    .. math::

       \\sum_i {\\left| u_i - v_i \\right|}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    cityblock : double
        The City Block (Manhattan) distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    return abs(u - v).sum()


def mahalanobis(u, v, VI):
    """
    Computes the Mahalanobis distance between two 1-D arrays.

    The Mahalanobis distance between 1-D arrays `u` and `v`, is defined as

    .. math::

       \\sqrt{ (u-v) V^{-1} (u-v)^T }

    where ``V`` is the covariance matrix.  Note that the argument `VI`
    is the inverse of ``V``.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    VI : ndarray
        The inverse of the covariance matrix.

    Returns
    -------
    mahalanobis : double
        The Mahalanobis distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    VI = np.atleast_2d(VI)
    delta = u - v
    m = np.dot(np.dot(delta, VI), delta)
    return np.sqrt(m)


def chebyshev(u, v):
    """
    Computes the Chebyshev distance.

    Computes the Chebyshev distance between two 1-D arrays `u` and `v`,
    which is defined as

    .. math::

       \\max_i {|u_i-v_i|}.

    Parameters
    ----------
    u : (N,) array_like
        Input vector.
    v : (N,) array_like
        Input vector.

    Returns
    -------
    chebyshev : double
        The Chebyshev distance between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    return max(abs(u - v))


def braycurtis(u, v):
    """
    Computes the Bray-Curtis distance between two 1-D arrays.

    Bray-Curtis distance is defined as

    .. math::

       \\sum{|u_i-v_i|} / \\sum{|u_i+v_i|}

    The Bray-Curtis distance is in the range [0, 1] if all coordinates are
    positive, and is undefined if the inputs are of length zero.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    braycurtis : double
        The Bray-Curtis distance between 1-D arrays `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v, dtype=np.float64)
    return abs(u - v).sum() / abs(u + v).sum()


def canberra(u, v):
    """
    Computes the Canberra distance between two 1-D arrays.

    The Canberra distance is defined as

    .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    canberra : double
        The Canberra distance between vectors `u` and `v`.

    Notes
    -----
    When `u[i]` and `v[i]` are 0 for given i, then the fraction 0/0 = 0 is
    used in the calculation.

    """
    u = _validate_vector(u)
    v = _validate_vector(v, dtype=np.float64)
    olderr = np.seterr(invalid='ignore')
    try:
        d = np.nansum(abs(u - v) / (abs(u) + abs(v)))
    finally:
        np.seterr(**olderr)
    return d


def _nbool_correspond_all(u, v):
    if u.dtype != v.dtype:
        raise TypeError("Arrays being compared must be of the same data type.")

    if u.dtype == int or u.dtype == np.float_ or u.dtype == np.double:
        not_u = 1.0 - u
        not_v = 1.0 - v
        nff = (not_u * not_v).sum()
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
        ntt = (u * v).sum()
    elif u.dtype == bool:
        not_u = ~u
        not_v = ~v
        nff = (not_u & not_v).sum()
        nft = (not_u & v).sum()
        ntf = (u & not_v).sum()
        ntt = (u & v).sum()
    else:
        raise TypeError("Arrays being compared have unknown type.")

    return (nff, nft, ntf, ntt)


def _nbool_correspond_ft_tf(u, v):
    if u.dtype == int or u.dtype == np.float_ or u.dtype == np.double:
        not_u = 1.0 - u
        not_v = 1.0 - v
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
    else:
        not_u = ~u
        not_v = ~v
        nft = (not_u & v).sum()
        ntf = (u & not_v).sum()
    return (nft, ntf)


def yule(u, v):
    """
    Computes the Yule dissimilarity between two boolean 1-D arrays.

    The Yule dissimilarity is defined as

    .. math::

         \\frac{R}{c_{TT} * c_{FF} + \\frac{R}{2}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2.0 * c_{TF} * c_{FT}`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    yule : double
        The Yule dissimilarity between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
    return float(2.0 * ntf * nft) / float(ntt * nff + ntf * nft)


def matching(u, v):
    """
    Computes the Matching dissimilarity between two boolean 1-D arrays.

    The Matching dissimilarity between two boolean 1-D arrays
    `u` and `v`, is defined as

    .. math::

       \\frac{c_{TF} + c_{FT}}{n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    matching : double
        The Matching dissimilarity between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(nft + ntf) / float(len(u))


def dice(u, v):
    """
    Computes the Dice dissimilarity between two boolean 1-D arrays.

    The Dice dissimilarity between `u` and `v`, is

    .. math::

         \\frac{c_{TF} + c_{FT}}
              {2c_{TT} + c_{FT} + c_{TF}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) ndarray, bool
        Input 1-D array.
    v : (N,) ndarray, bool
        Input 1-D array.

    Returns
    -------
    dice : double
        The Dice dissimilarity between 1-D arrays `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(ntf + nft) / float(2.0 * ntt + ntf + nft)


def rogerstanimoto(u, v):
    """
    Computes the Rogers-Tanimoto dissimilarity between two boolean 1-D arrays.

    The Rogers-Tanimoto dissimilarity between two boolean 1-D arrays
    `u` and `v`, is defined as

    .. math::
       \\frac{R}
            {c_{TT} + c_{FF} + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2(c_{TF} + c_{FT})`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    rogerstanimoto : double
        The Rogers-Tanimoto dissimilarity between vectors
        `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
    return float(2.0 * (ntf + nft)) / float(ntt + nff + (2.0 * (ntf + nft)))


def russellrao(u, v):
    """
    Computes the Russell-Rao dissimilarity between two boolean 1-D arrays.

    The Russell-Rao dissimilarity between two boolean 1-D arrays, `u` and
    `v`, is defined as

    .. math::

      \\frac{n - c_{TT}}
           {n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    russellrao : double
        The Russell-Rao dissimilarity between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    return float(len(u) - ntt) / float(len(u))


def sokalmichener(u, v):
    """
    Computes the Sokal-Michener dissimilarity between two boolean 1-D arrays.

    The Sokal-Michener dissimilarity between boolean 1-D arrays `u` and `v`,
    is defined as

    .. math::

       \\frac{R}
            {S + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`, :math:`R = 2 * (c_{TF} + c_{FT})` and
    :math:`S = c_{FF} + c_{TT}`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    sokalmichener : double
        The Sokal-Michener dissimilarity between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == bool:
        ntt = (u & v).sum()
        nff = (~u & ~v).sum()
    else:
        ntt = (u * v).sum()
        nff = ((1.0 - u) * (1.0 - v)).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(2.0 * (ntf + nft)) / float(ntt + nff + 2.0 * (ntf + nft))


def sokalsneath(u, v):
    """
    Computes the Sokal-Sneath dissimilarity between two boolean 1-D arrays.

    The Sokal-Sneath dissimilarity between `u` and `v`,

    .. math::

       \\frac{R}
            {c_{TT} + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2(c_{TF} + c_{FT})`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.

    Returns
    -------
    sokalsneath : double
        The Sokal-Sneath dissimilarity between vectors `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    denom = ntt + 2.0 * (ntf + nft)
    if denom == 0:
        raise ValueError('Sokal-Sneath dissimilarity is not defined for '
                            'vectors that are entirely false.')
    return float(2.0 * (ntf + nft)) / denom


def pdist(X, metric='euclidean', p=2, w=None, V=None, VI=None):
    """
    Pairwise distances between observations in n-dimensional space.

    The following are common calling conventions.

    1. ``Y = pdist(X, 'euclidean')``

       Computes the distance between m points using Euclidean distance
       (2-norm) as the distance metric between the points. The points
       are arranged as m n-dimensional row vectors in the matrix X.

    2. ``Y = pdist(X, 'minkowski', p)``

       Computes the distances using the Minkowski distance
       :math:`||u-v||_p` (p-norm) where :math:`p \\geq 1`.

    3. ``Y = pdist(X, 'cityblock')``

       Computes the city block or Manhattan distance between the
       points.

    4. ``Y = pdist(X, 'seuclidean', V=None)``

       Computes the standardized Euclidean distance. The standardized
       Euclidean distance between two n-vectors ``u`` and ``v`` is

       .. math::

          \\sqrt{\\sum {(u_i-v_i)^2 / V[x_i]}}


       V is the variance vector; V[i] is the variance computed over all
       the i'th components of the points.  If not passed, it is
       automatically computed.

    5. ``Y = pdist(X, 'sqeuclidean')``

       Computes the squared Euclidean distance :math:`||u-v||_2^2` between
       the vectors.

    6. ``Y = pdist(X, 'cosine')``

       Computes the cosine distance between vectors u and v,

       .. math::

          1 - \\frac{u \\cdot v}
                   {{||u||}_2 {||v||}_2}

       where :math:`||*||_2` is the 2-norm of its argument ``*``, and
       :math:`u \\cdot v` is the dot product of ``u`` and ``v``.

    7. ``Y = pdist(X, 'correlation')``

       Computes the correlation distance between vectors u and v. This is

       .. math::

          1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
                   {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

       where :math:`\\bar{v}` is the mean of the elements of vector v,
       and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.

    8. ``Y = pdist(X, 'hamming')``

       Computes the normalized Hamming distance, or the proportion of
       those vector elements between two n-vectors ``u`` and ``v``
       which disagree. To save memory, the matrix ``X`` can be of type
       boolean.

    9. ``Y = pdist(X, 'jaccard')``

       Computes the Jaccard distance between the points. Given two
       vectors, ``u`` and ``v``, the Jaccard distance is the
       proportion of those elements ``u[i]`` and ``v[i]`` that
       disagree where at least one of them is non-zero.

    10. ``Y = pdist(X, 'chebyshev')``

       Computes the Chebyshev distance between the points. The
       Chebyshev distance between two n-vectors ``u`` and ``v`` is the
       maximum norm-1 distance between their respective elements. More
       precisely, the distance is given by

       .. math::

          d(u,v) = \\max_i {|u_i-v_i|}

    11. ``Y = pdist(X, 'canberra')``

       Computes the Canberra distance between the points. The
       Canberra distance between two points ``u`` and ``v`` is

       .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}


    12. ``Y = pdist(X, 'braycurtis')``

       Computes the Bray-Curtis distance between the points. The
       Bray-Curtis distance between two points ``u`` and ``v`` is


       .. math::

            d(u,v) = \\frac{\\sum_i {u_i-v_i}}
                          {\\sum_i {u_i+v_i}}

    13. ``Y = pdist(X, 'mahalanobis', VI=None)``

       Computes the Mahalanobis distance between the points. The
       Mahalanobis distance between two points ``u`` and ``v`` is
       :math:`(u-v)(1/V)(u-v)^T` where :math:`(1/V)` (the ``VI``
       variable) is the inverse covariance. If ``VI`` is not None,
       ``VI`` will be used as the inverse covariance matrix.

    14. ``Y = pdist(X, 'yule')``

       Computes the Yule distance between each pair of boolean
       vectors. (see yule function documentation)

    15. ``Y = pdist(X, 'matching')``

       Computes the matching distance between each pair of boolean
       vectors. (see matching function documentation)

    16. ``Y = pdist(X, 'dice')``

       Computes the Dice distance between each pair of boolean
       vectors. (see dice function documentation)

    17. ``Y = pdist(X, 'kulsinski')``

       Computes the Kulsinski distance between each pair of
       boolean vectors. (see kulsinski function documentation)

    18. ``Y = pdist(X, 'rogerstanimoto')``

       Computes the Rogers-Tanimoto distance between each pair of
       boolean vectors. (see rogerstanimoto function documentation)

    19. ``Y = pdist(X, 'russellrao')``

       Computes the Russell-Rao distance between each pair of
       boolean vectors. (see russellrao function documentation)

    20. ``Y = pdist(X, 'sokalmichener')``

       Computes the Sokal-Michener distance between each pair of
       boolean vectors. (see sokalmichener function documentation)

    21. ``Y = pdist(X, 'sokalsneath')``

       Computes the Sokal-Sneath distance between each pair of
       boolean vectors. (see sokalsneath function documentation)

    22. ``Y = pdist(X, 'wminkowski')``

       Computes the weighted Minkowski distance between each pair of
       vectors. (see wminkowski function documentation)

    23. ``Y = pdist(X, f)``

       Computes the distance between all pairs of vectors in X
       using the user supplied 2-arity function f. For example,
       Euclidean distance between the vectors could be computed
       as follows::

         dm = pdist(X, lambda u, v: np.sqrt(((u-v)**2).sum()))

       Note that you should avoid passing a reference to one of
       the distance functions defined in this library. For example,::

         dm = pdist(X, sokalsneath)

       would calculate the pair-wise distances between the vectors in
       X using the Python function sokalsneath. This would result in
       sokalsneath being called :math:`{n \\choose 2}` times, which
       is inefficient. Instead, the optimized C version is more
       efficient, and we call it using the following syntax.::

         dm = pdist(X, 'sokalsneath')

    Parameters
    ----------
    X : ndarray
        An m by n array of m original observations in an
        n-dimensional space.
    metric : str or function, optional
        The distance metric to use. The distance function can
        be 'braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'kulsinski', 'mahalanobis', 'matching',
        'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
        'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'.
    w : ndarray, optional
        The weight vector (for weighted Minkowski).
    p : double, optional
        The p-norm to apply (for Minkowski, weighted and unweighted)
    V : ndarray, optional
        The variance vector (for standardized Euclidean).
    VI : ndarray, optional
        The inverse of the covariance matrix (for Mahalanobis).

    Returns
    -------
    Y : ndarray
        Returns a condensed distance matrix Y.  For
        each :math:`i` and :math:`j` (where :math:`i<j<n`), the
        metric ``dist(u=X[i], v=X[j])`` is computed and stored in entry ``ij``.

    See Also
    --------
    squareform : converts between condensed distance matrices and
                 square distance matrices.

    Notes
    -----
    See ``squareform`` for information on how to calculate the index of
    this entry or to convert the condensed distance matrix to a
    redundant square matrix.

    """

#         21. Y = pdist(X, 'test_Y')
#
#           Computes the distance between all pairs of vectors in X
#           using the distance metric Y but with a more succinct,
#           verifiable, but less efficient implementation.

    X = np.asarray(X, order='c')

    # The C code doesn't do striding.
    [X] = _copy_arrays_if_base_present([_convert_to_double(X)])

    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1)) // 2, dtype=np.double)

    wmink_names = ['wminkowski', 'wmi', 'wm', 'wpnorm']
    if w is None and (metric == wminkowski or metric in wmink_names):
        raise ValueError('weighted minkowski requires a weight '
                            'vector `w` to be given.')

    if callable(metric):
        if metric == minkowski:
            def dfun(u, v):
                return minkowski(u, v, p)
        elif metric == wminkowski:
            def dfun(u, v):
                return wminkowski(u, v, p, w)
        elif metric == seuclidean:
            def dfun(u, v):
                return seuclidean(u, v, V)
        elif metric == mahalanobis:
            def dfun(u, v):
                return mahalanobis(u, v, V)
        else:
            dfun = metric

        k = 0
        for i in xrange(0, m - 1):
            for j in xrange(i + 1, m):
                dm[k] = dfun(X[i], X[j])
                k = k + 1

    elif isinstance(metric, string_types):
        mstr = metric.lower()

        #if X.dtype != np.double and \
        #       (mstr != 'hamming' and mstr != 'jaccard'):
        #    TypeError('A double array must be passed.')
        if mstr in set(['euclidean', 'euclid', 'eu', 'e']):
            _distance_wrap.pdist_euclidean_wrap(_convert_to_double(X), dm)
        elif mstr in set(['sqeuclidean', 'sqe', 'sqeuclid']):
            _distance_wrap.pdist_sqeuclidean_wrap(_convert_to_double(X), dm)
        elif mstr in set(['cityblock', 'cblock', 'cb', 'c']):
            _distance_wrap.pdist_city_block_wrap(X, dm)
        elif mstr in set(['hamming', 'hamm', 'ha', 'h']):
            if X.dtype == bool:
                _distance_wrap.pdist_hamming_bool_wrap(_convert_to_bool(X), dm)
            else:
                _distance_wrap.pdist_hamming_wrap(_convert_to_double(X), dm)
        elif mstr in set(['jaccard', 'jacc', 'ja', 'j']):
            if X.dtype == bool:
                _distance_wrap.pdist_jaccard_bool_wrap(_convert_to_bool(X), dm)
            else:
                _distance_wrap.pdist_jaccard_wrap(_convert_to_double(X), dm)
        elif mstr in set(['chebychev', 'chebyshev', 'cheby', 'cheb', 'ch']):
            _distance_wrap.pdist_chebyshev_wrap(_convert_to_double(X), dm)
        elif mstr in set(['minkowski', 'mi', 'm']):
            _distance_wrap.pdist_minkowski_wrap(_convert_to_double(X), dm, p)
        elif mstr in wmink_names:
            w = _convert_to_double(np.asarray(w))
            _distance_wrap.pdist_weighted_minkowski_wrap(_convert_to_double(X),
                                                         dm, p, w)
        elif mstr in set(['seuclidean', 'se', 's']):
            if V is not None:
                V = np.asarray(V, order='c')
                if type(V) != np.ndarray:
                    raise TypeError('Variance vector V must be a numpy array')
                if V.dtype != np.double:
                    raise TypeError('Variance vector V must contain doubles.')
                if len(V.shape) != 1:
                    raise ValueError('Variance vector V must '
                                     'be one-dimensional.')
                if V.shape[0] != n:
                    raise ValueError('Variance vector V must be of the same '
                            'dimension as the vectors on which the distances '
                            'are computed.')
                # The C code doesn't do striding.
                [VV] = _copy_arrays_if_base_present([_convert_to_double(V)])
            else:
                VV = np.var(X, axis=0, ddof=1)
            _distance_wrap.pdist_seuclidean_wrap(_convert_to_double(X), VV, dm)
        elif mstr in set(['cosine', 'cos']):
            norms = _row_norms(X)
            _distance_wrap.pdist_cosine_wrap(_convert_to_double(X), dm, norms)
        elif mstr in set(['old_cosine', 'old_cos']):
            norms = _row_norms(X)
            nV = norms.reshape(m, 1)
            # The numerator u * v
            nm = np.dot(X, X.T)
            # The denom. ||u||*||v||
            de = np.dot(nV, nV.T)
            dm = 1.0 - (nm / de)
            dm[xrange(0, m), xrange(0, m)] = 0.0
            dm = squareform(dm)
        elif mstr in set(['correlation', 'co']):
            X2 = X - X.mean(1)[:, np.newaxis]
            norms = _row_norms(X2)
            _distance_wrap.pdist_cosine_wrap(_convert_to_double(X2),
                                             _convert_to_double(dm),
                                             _convert_to_double(norms))
        elif mstr in set(['mahalanobis', 'mahal', 'mah']):
            if VI is not None:
                VI = _convert_to_double(np.asarray(VI, order='c'))
                if type(VI) != np.ndarray:
                    raise TypeError('VI must be a numpy array.')
                if VI.dtype != np.double:
                    raise TypeError('The array must contain 64-bit floats.')
                [VI] = _copy_arrays_if_base_present([VI])
            else:
                if m <= n:
                    # There are fewer observations than the dimension of
                    # the observations.
                    raise ValueError("The number of observations (%d) is too "
                                     "small; the covariance matrix is "
                                     "singular. For observations with %d "
                                     "dimensions, at least %d observations "
                                     "are required." % (m, n, n + 1))
                V = np.atleast_2d(np.cov(X.T))
                VI = _convert_to_double(np.linalg.inv(V).T.copy())
            # (u-v)V^(-1)(u-v)^T
            _distance_wrap.pdist_mahalanobis_wrap(_convert_to_double(X),
                                                  VI, dm)
        elif mstr == 'canberra':
            _distance_wrap.pdist_canberra_wrap(_convert_to_double(X), dm)
        elif mstr == 'braycurtis':
            _distance_wrap.pdist_bray_curtis_wrap(_convert_to_double(X), dm)
        elif mstr == 'yule':
            _distance_wrap.pdist_yule_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'matching':
            _distance_wrap.pdist_matching_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'kulsinski':
            _distance_wrap.pdist_kulsinski_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'dice':
            _distance_wrap.pdist_dice_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'rogerstanimoto':
            _distance_wrap.pdist_rogerstanimoto_bool_wrap(_convert_to_bool(X),
                                                          dm)
        elif mstr == 'russellrao':
            _distance_wrap.pdist_russellrao_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'sokalmichener':
            _distance_wrap.pdist_sokalmichener_bool_wrap(_convert_to_bool(X),
                                                         dm)
        elif mstr == 'sokalsneath':
            _distance_wrap.pdist_sokalsneath_bool_wrap(_convert_to_bool(X), dm)
        elif metric == 'test_euclidean':
            dm = pdist(X, euclidean)
        elif metric == 'test_sqeuclidean':
            if V is None:
                V = np.var(X, axis=0, ddof=1)
            else:
                V = np.asarray(V, order='c')
            dm = pdist(X, lambda u, v: seuclidean(u, v, V))
        elif metric == 'test_braycurtis':
            dm = pdist(X, braycurtis)
        elif metric == 'test_mahalanobis':
            if VI is None:
                V = np.cov(X.T)
                VI = np.linalg.inv(V)
            else:
                VI = np.asarray(VI, order='c')
            [VI] = _copy_arrays_if_base_present([VI])
            # (u-v)V^(-1)(u-v)^T
            dm = pdist(X, (lambda u, v: mahalanobis(u, v, VI)))
        elif metric == 'test_canberra':
            dm = pdist(X, canberra)
        elif metric == 'test_cityblock':
            dm = pdist(X, cityblock)
        elif metric == 'test_minkowski':
            dm = pdist(X, minkowski, p=p)
        elif metric == 'test_wminkowski':
            dm = pdist(X, wminkowski, p=p, w=w)
        elif metric == 'test_cosine':
            dm = pdist(X, cosine)
        elif metric == 'test_correlation':
            dm = pdist(X, correlation)
        elif metric == 'test_hamming':
            dm = pdist(X, hamming)
        elif metric == 'test_jaccard':
            dm = pdist(X, jaccard)
        elif metric == 'test_chebyshev' or metric == 'test_chebychev':
            dm = pdist(X, chebyshev)
        elif metric == 'test_yule':
            dm = pdist(X, yule)
        elif metric == 'test_matching':
            dm = pdist(X, matching)
        elif metric == 'test_dice':
            dm = pdist(X, dice)
        elif metric == 'test_kulsinski':
            dm = pdist(X, kulsinski)
        elif metric == 'test_rogerstanimoto':
            dm = pdist(X, rogerstanimoto)
        elif metric == 'test_russellrao':
            dm = pdist(X, russellrao)
        elif metric == 'test_sokalsneath':
            dm = pdist(X, sokalsneath)
        elif metric == 'test_sokalmichener':
            dm = pdist(X, sokalmichener)
        else:
            raise ValueError('Unknown Distance Metric: %s' % mstr)
    else:
        raise TypeError('2nd argument metric must be a string identifier '
                        'or a function.')
    return dm


def squareform(X, force="no", checks=True):
    """
    Converts a vector-form distance vector to a square-form distance
    matrix, and vice-versa.

    Parameters
    ----------
    X : ndarray
        Either a condensed or redundant distance matrix.
    force : str, optional
        As with MATLAB(TM), if force is equal to 'tovector' or 'tomatrix',
        the input will be treated as a distance matrix or distance vector
        respectively.
    checks : bool, optional
        If `checks` is set to False, no checks will be made for matrix
        symmetry nor zero diagonals. This is useful if it is known that
        ``X - X.T1`` is small and ``diag(X)`` is close to zero.
        These values are ignored any way so they do not disrupt the
        squareform transformation.

    Returns
    -------
    Y : ndarray
        If a condensed distance matrix is passed, a redundant one is
        returned, or if a redundant one is passed, a condensed distance
        matrix is returned.

    Notes
    -----

    1. v = squareform(X)

       Given a square d-by-d symmetric distance matrix X,
       ``v=squareform(X)`` returns a ``d * (d-1) / 2`` (or
       `${n \\choose 2}$`) sized vector v.

      v[{n \\choose 2}-{n-i \\choose 2} + (j-i-1)] is the distance
      between points i and j. If X is non-square or asymmetric, an error
      is returned.

    2. X = squareform(v)

      Given a d*d(-1)/2 sized v for some integer d>=2 encoding distances
      as described, X=squareform(v) returns a d by d distance matrix X. The
      X[i, j] and X[j, i] values are set to
      v[{n \\choose 2}-{n-i \\choose 2} + (j-u-1)] and all
      diagonal elements are zero.

    """

    X = _convert_to_double(np.asarray(X, order='c'))

    if not np.issubsctype(X, np.double):
        raise TypeError('A double array must be passed.')

    s = X.shape

    if force.lower() == 'tomatrix':
        if len(s) != 1:
            raise ValueError("Forcing 'tomatrix' but input X is not a "
                             "distance vector.")
    elif force.lower() == 'tovector':
        if len(s) != 2:
            raise ValueError("Forcing 'tovector' but input X is not a "
                             "distance matrix.")

    # X = squareform(v)
    if len(s) == 1:
        if X.shape[0] == 0:
            return np.zeros((1, 1), dtype=np.double)

        # Grab the closest value to the square root of the number
        # of elements times 2 to see if the number of elements
        # is indeed a binomial coefficient.
        d = int(np.ceil(np.sqrt(X.shape[0] * 2)))

        # Check that v is of valid dimensions.
        if d * (d - 1) / 2 != int(s[0]):
            raise ValueError('Incompatible vector size. It must be a binomial '
                             'coefficient n choose 2 for some integer n >= 2.')

        # Allocate memory for the distance matrix.
        M = np.zeros((d, d), dtype=np.double)

        # Since the C code does not support striding using strides.
        # The dimensions are used instead.
        [X] = _copy_arrays_if_base_present([X])

        # Fill in the values of the distance matrix.
        _distance_wrap.to_squareform_from_vector_wrap(M, X)

        # Return the distance matrix.
        return M
    elif len(s) == 2:
        if s[0] != s[1]:
            raise ValueError('The matrix argument must be square.')
        if checks:
            is_valid_dm(X, throw=True, name='X')

        # One-side of the dimensions is set here.
        d = s[0]

        if d <= 1:
            return np.array([], dtype=np.double)

        # Create a vector.
        v = np.zeros((d * (d - 1)) // 2, dtype=np.double)

        # Since the C code does not support striding using strides.
        # The dimensions are used instead.
        [X] = _copy_arrays_if_base_present([X])

        # Convert the vector to squareform.
        _distance_wrap.to_vector_from_squareform_wrap(X, v)
        return v
    else:
        raise ValueError(('The first argument must be one or two dimensional '
                         'array. A %d-dimensional array is not '
                         'permitted') % len(s))


def is_valid_dm(D, tol=0.0, throw=False, name="D", warning=False):
    """
    Returns True if input array is a valid distance matrix.

    Distance matrices must be 2-dimensional numpy arrays containing
    doubles. They must have a zero-diagonal, and they must be symmetric.

    Parameters
    ----------
    D : ndarray
        The candidate object to test for validity.
    tol : float, optional
        The distance matrix should be symmetric. `tol` is the maximum
        difference between entries ``ij`` and ``ji`` for the distance
        metric to be considered symmetric.
    throw : bool, optional
        An exception is thrown if the distance matrix passed is not valid.
    name : str, optional
        The name of the variable to checked. This is useful if
        throw is set to True so the offending variable can be identified
        in the exception message when an exception is thrown.
    warning : bool, optional
        Instead of throwing an exception, a warning message is
        raised.

    Returns
    -------
    valid : bool
        True if the variable `D` passed is a valid distance matrix.

    Notes
    -----
    Small numerical differences in `D` and `D.T` and non-zeroness of
    the diagonal are ignored if they are within the tolerance specified
    by `tol`.

    """
    D = np.asarray(D, order='c')
    valid = True
    try:
        s = D.shape
        if D.dtype != np.double:
            if name:
                raise TypeError(('Distance matrix \'%s\' must contain doubles '
                                 '(double).') % name)
            else:
                raise TypeError('Distance matrix must contain doubles '
                                '(double).')
        if len(D.shape) != 2:
            if name:
                raise ValueError(('Distance matrix \'%s\' must have shape=2 '
                                 '(i.e. be two-dimensional).') % name)
            else:
                raise ValueError('Distance matrix must have shape=2 (i.e. '
                                 'be two-dimensional).')
        if tol == 0.0:
            if not (D == D.T).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' must be '
                                     'symmetric.') % name)
                else:
                    raise ValueError('Distance matrix must be symmetric.')
            if not (D[xrange(0, s[0]), xrange(0, s[0])] == 0).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' diagonal must '
                                     'be zero.') % name)
                else:
                    raise ValueError('Distance matrix diagonal must be zero.')
        else:
            if not (D - D.T <= tol).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' must be '
                                      'symmetric within tolerance %5.5f.')
                                     % (name, tol))
                else:
                    raise ValueError('Distance matrix must be symmetric within'
                                     ' tolerance %5.5f.' % tol)
            if not (D[xrange(0, s[0]), xrange(0, s[0])] <= tol).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' diagonal must be'
                                      ' close to zero within tolerance %5.5f.')
                                     % (name, tol))
                else:
                    raise ValueError(('Distance matrix \'%s\' diagonal must be'
                                      ' close to zero within tolerance %5.5f.')
                                     % tol)
    except Exception as e:
        if throw:
            raise
        if warning:
            warnings.warn(str(e))
        valid = False
    return valid


def is_valid_y(y, warning=False, throw=False, name=None):
    """
    Returns True if the input array is a valid condensed distance matrix.

    Condensed distance matrices must be 1-dimensional
    numpy arrays containing doubles. Their length must be a binomial
    coefficient :math:`{n \\choose 2}` for some positive integer n.

    Parameters
    ----------
    y : ndarray
        The condensed distance matrix.
    warning : bool, optional
        Invokes a warning if the variable passed is not a valid
        condensed distance matrix. The warning message explains why
        the distance matrix is not valid.  `name` is used when
        referencing the offending variable.
    throw : bool, optional
        Throws an exception if the variable passed is not a valid
        condensed distance matrix.
    name : bool, optional
        Used when referencing the offending variable in the
        warning or exception message.

    """
    y = np.asarray(y, order='c')
    valid = True
    try:
        if type(y) != np.ndarray:
            if name:
                raise TypeError(('\'%s\' passed as a condensed distance '
                                 'matrix is not a numpy array.') % name)
            else:
                raise TypeError('Variable is not a numpy array.')
        if y.dtype != np.double:
            if name:
                raise TypeError(('Condensed distance matrix \'%s\' must '
                                 'contain doubles (double).') % name)
            else:
                raise TypeError('Condensed distance matrix must contain '
                                'doubles (double).')
        if len(y.shape) != 1:
            if name:
                raise ValueError(('Condensed distance matrix \'%s\' must '
                                  'have shape=1 (i.e. be one-dimensional).')
                                 % name)
            else:
                raise ValueError('Condensed distance matrix must have shape=1 '
                                 '(i.e. be one-dimensional).')
        n = y.shape[0]
        d = int(np.ceil(np.sqrt(n * 2)))
        if (d * (d - 1) / 2) != n:
            if name:
                raise ValueError(('Length n of condensed distance matrix '
                                  '\'%s\' must be a binomial coefficient, i.e.'
                                  'there must be a k such that '
                                  '(k \choose 2)=n)!') % name)
            else:
                raise ValueError('Length n of condensed distance matrix must '
                                 'be a binomial coefficient, i.e. there must '
                                 'be a k such that (k \choose 2)=n)!')
    except Exception as e:
        if throw:
            raise
        if warning:
            warnings.warn(str(e))
        valid = False
    return valid


def num_obs_dm(d):
    """
    Returns the number of original observations that correspond to a
    square, redundant distance matrix.

    Parameters
    ----------
    d : ndarray
        The target distance matrix.

    Returns
    -------
    num_obs_dm : int
        The number of observations in the redundant distance matrix.

    """
    d = np.asarray(d, order='c')
    is_valid_dm(d, tol=np.inf, throw=True, name='d')
    return d.shape[0]


def num_obs_y(Y):
    """
    Returns the number of original observations that correspond to a
    condensed distance matrix.

    Parameters
    ----------
    Y : ndarray
        Condensed distance matrix.

    Returns
    -------
    n : int
        The number of observations in the condensed distance matrix `Y`.

    """
    Y = np.asarray(Y, order='c')
    is_valid_y(Y, throw=True, name='Y')
    k = Y.shape[0]
    if k == 0:
        raise ValueError("The number of observations cannot be determined on "
                         "an empty distance matrix.")
    d = int(np.ceil(np.sqrt(k * 2)))
    if (d * (d - 1) / 2) != k:
        raise ValueError("Invalid condensed distance matrix passed. Must be "
                         "some k where k=(n choose 2) for some n >= 2.")
    return d


def _row_norms(X):
    norms = np.einsum('ij,ij->i', X, X)
    return np.sqrt(norms, out=norms)


def _cosine_cdist(XA, XB, dm):
    XA = _convert_to_double(XA)
    XB = _convert_to_double(XB)

    normsA = _row_norms(XA)
    normsB = _row_norms(XB)

    np.dot(XA, XB.T, out=dm)

    dm /= normsA.reshape(-1, 1)
    dm /= normsB
    dm *= -1
    dm += 1


def cdist(XA, XB, metric='euclidean', p=2, V=None, VI=None, w=None):
    """
    Computes distance between each pair of the two collections of inputs.

    The following are common calling conventions:

    1. ``Y = cdist(XA, XB, 'euclidean')``

       Computes the distance between :math:`m` points using
       Euclidean distance (2-norm) as the distance metric between the
       points. The points are arranged as :math:`m`
       :math:`n`-dimensional row vectors in the matrix X.

    2. ``Y = cdist(XA, XB, 'minkowski', p)``

       Computes the distances using the Minkowski distance
       :math:`||u-v||_p` (:math:`p`-norm) where :math:`p \\geq 1`.

    3. ``Y = cdist(XA, XB, 'cityblock')``

       Computes the city block or Manhattan distance between the
       points.

    4. ``Y = cdist(XA, XB, 'seuclidean', V=None)``

       Computes the standardized Euclidean distance. The standardized
       Euclidean distance between two n-vectors ``u`` and ``v`` is

       .. math::

          \\sqrt{\\sum {(u_i-v_i)^2 / V[x_i]}}.

       V is the variance vector; V[i] is the variance computed over all
       the i'th components of the points. If not passed, it is
       automatically computed.

    5. ``Y = cdist(XA, XB, 'sqeuclidean')``

       Computes the squared Euclidean distance :math:`||u-v||_2^2` between
       the vectors.

    6. ``Y = cdist(XA, XB, 'cosine')``

       Computes the cosine distance between vectors u and v,

       .. math::

          1 - \\frac{u \\cdot v}
                   {{||u||}_2 {||v||}_2}

       where :math:`||*||_2` is the 2-norm of its argument ``*``, and
       :math:`u \\cdot v` is the dot product of :math:`u` and :math:`v`.

    7. ``Y = cdist(XA, XB, 'correlation')``

       Computes the correlation distance between vectors u and v. This is

       .. math::

          1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
                   {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

       where :math:`\\bar{v}` is the mean of the elements of vector v,
       and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.


    8. ``Y = cdist(XA, XB, 'hamming')``

       Computes the normalized Hamming distance, or the proportion of
       those vector elements between two n-vectors ``u`` and ``v``
       which disagree. To save memory, the matrix ``X`` can be of type
       boolean.

    9. ``Y = cdist(XA, XB, 'jaccard')``

       Computes the Jaccard distance between the points. Given two
       vectors, ``u`` and ``v``, the Jaccard distance is the
       proportion of those elements ``u[i]`` and ``v[i]`` that
       disagree where at least one of them is non-zero.

    10. ``Y = cdist(XA, XB, 'chebyshev')``

       Computes the Chebyshev distance between the points. The
       Chebyshev distance between two n-vectors ``u`` and ``v`` is the
       maximum norm-1 distance between their respective elements. More
       precisely, the distance is given by

       .. math::

          d(u,v) = \\max_i {|u_i-v_i|}.

    11. ``Y = cdist(XA, XB, 'canberra')``

       Computes the Canberra distance between the points. The
       Canberra distance between two points ``u`` and ``v`` is

       .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}.

    12. ``Y = cdist(XA, XB, 'braycurtis')``

       Computes the Bray-Curtis distance between the points. The
       Bray-Curtis distance between two points ``u`` and ``v`` is


       .. math::

            d(u,v) = \\frac{\\sum_i (u_i-v_i)}
                          {\\sum_i (u_i+v_i)}

    13. ``Y = cdist(XA, XB, 'mahalanobis', VI=None)``

       Computes the Mahalanobis distance between the points. The
       Mahalanobis distance between two points ``u`` and ``v`` is
       :math:`(u-v)(1/V)(u-v)^T` where :math:`(1/V)` (the ``VI``
       variable) is the inverse covariance. If ``VI`` is not None,
       ``VI`` will be used as the inverse covariance matrix.

    14. ``Y = cdist(XA, XB, 'yule')``

       Computes the Yule distance between the boolean
       vectors. (see `yule` function documentation)

    15. ``Y = cdist(XA, XB, 'matching')``

       Computes the matching distance between the boolean
       vectors. (see `matching` function documentation)

    16. ``Y = cdist(XA, XB, 'dice')``

       Computes the Dice distance between the boolean vectors. (see
       `dice` function documentation)

    17. ``Y = cdist(XA, XB, 'kulsinski')``

       Computes the Kulsinski distance between the boolean
       vectors. (see `kulsinski` function documentation)

    18. ``Y = cdist(XA, XB, 'rogerstanimoto')``

       Computes the Rogers-Tanimoto distance between the boolean
       vectors. (see `rogerstanimoto` function documentation)

    19. ``Y = cdist(XA, XB, 'russellrao')``

       Computes the Russell-Rao distance between the boolean
       vectors. (see `russellrao` function documentation)

    20. ``Y = cdist(XA, XB, 'sokalmichener')``

       Computes the Sokal-Michener distance between the boolean
       vectors. (see `sokalmichener` function documentation)

    21. ``Y = cdist(XA, XB, 'sokalsneath')``

       Computes the Sokal-Sneath distance between the vectors. (see
       `sokalsneath` function documentation)


    22. ``Y = cdist(XA, XB, 'wminkowski')``

       Computes the weighted Minkowski distance between the
       vectors. (see `wminkowski` function documentation)

    23. ``Y = cdist(XA, XB, f)``

       Computes the distance between all pairs of vectors in X
       using the user supplied 2-arity function f. For example,
       Euclidean distance between the vectors could be computed
       as follows::

         dm = cdist(XA, XB, lambda u, v: np.sqrt(((u-v)**2).sum()))

       Note that you should avoid passing a reference to one of
       the distance functions defined in this library. For example,::

         dm = cdist(XA, XB, sokalsneath)

       would calculate the pair-wise distances between the vectors in
       X using the Python function `sokalsneath`. This would result in
       sokalsneath being called :math:`{n \\choose 2}` times, which
       is inefficient. Instead, the optimized C version is more
       efficient, and we call it using the following syntax::

         dm = cdist(XA, XB, 'sokalsneath')

    Parameters
    ----------
    XA : ndarray
        An :math:`m_A` by :math:`n` array of :math:`m_A`
        original observations in an :math:`n`-dimensional space.
        Inputs are converted to float type.
    XB : ndarray
        An :math:`m_B` by :math:`n` array of :math:`m_B`
        original observations in an :math:`n`-dimensional space.
        Inputs are converted to float type.
    metric : str or callable, optional
        The distance metric to use.  If a string, the distance function can be
        'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
        'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
        'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
        'wminkowski', 'yule'.
    w : ndarray, optional
        The weight vector (for weighted Minkowski).
    p : scalar, optional
        The p-norm to apply (for Minkowski, weighted and unweighted)
    V : ndarray, optional
        The variance vector (for standardized Euclidean).
    VI : ndarray, optional
        The inverse of the covariance matrix (for Mahalanobis).

    Returns
    -------
    Y : ndarray
        A :math:`m_A` by :math:`m_B` distance matrix is returned.
        For each :math:`i` and :math:`j`, the metric
        ``dist(u=XA[i], v=XB[j])`` is computed and stored in the
        :math:`ij` th entry.

    Raises
    ------
    ValueError
        An exception is thrown if `XA` and `XB` do not have
        the same number of columns.

    Examples
    --------
    Find the Euclidean distances between four 2-D coordinates:

    >>> from scipy.spatial import distance
    >>> coords = [(35.0456, -85.2672),
    ...           (35.1174, -89.9711),
    ...           (35.9728, -83.9422),
    ...           (36.1667, -86.7833)]
    >>> distance.cdist(coords, coords, 'euclidean')
    array([[ 0.    ,  4.7044,  1.6172,  1.8856],
           [ 4.7044,  0.    ,  6.0893,  3.3561],
           [ 1.6172,  6.0893,  0.    ,  2.8477],
           [ 1.8856,  3.3561,  2.8477,  0.    ]])


    Find the Manhattan distance from a 3-D point to the corners of the unit
    cube:

    >>> a = np.array([[0, 0, 0],
                      [0, 0, 1],
                      [0, 1, 0],
                      [0, 1, 1],
                      [1, 0, 0],
                      [1, 0, 1],
                      [1, 1, 0],
                      [1, 1, 1]])
    >>> b = np.array([[ 0.1,  0.2,  0.4]])
    >>> distance.cdist(a, b, 'cityblock')
    array([[ 0.7],
           [ 0.9],
           [ 1.3],
           [ 1.5],
           [ 1.5],
           [ 1.7],
           [ 2.1],
           [ 2.3]])

    """

#         21. Y = cdist(XA, XB, 'test_Y')
#
#           Computes the distance between all pairs of vectors in X
#           using the distance metric Y but with a more succint,
#           verifiable, but less efficient implementation.

    XA = np.asarray(XA, order='c')
    XB = np.asarray(XB, order='c')

    #if np.issubsctype(X, np.floating) and not np.issubsctype(X, np.double):
    #    raise TypeError('Floating point arrays must be 64-bit (got %r).' %
    #    (X.dtype.type,))

    # The C code doesn't do striding.
    [XA] = _copy_arrays_if_base_present([_convert_to_double(XA)])
    [XB] = _copy_arrays_if_base_present([_convert_to_double(XB)])

    s = XA.shape
    sB = XB.shape

    if len(s) != 2:
        raise ValueError('XA must be a 2-dimensional array.')
    if len(sB) != 2:
        raise ValueError('XB must be a 2-dimensional array.')
    if s[1] != sB[1]:
        raise ValueError('XA and XB must have the same number of columns '
                         '(i.e. feature dimension.)')

    mA = s[0]
    mB = sB[0]
    n = s[1]
    dm = np.zeros((mA, mB), dtype=np.double)

    if callable(metric):
        if metric == minkowski:
            for i in xrange(0, mA):
                for j in xrange(0, mB):
                    dm[i, j] = minkowski(XA[i, :], XB[j, :], p)
        elif metric == wminkowski:
            for i in xrange(0, mA):
                for j in xrange(0, mB):
                    dm[i, j] = wminkowski(XA[i, :], XB[j, :], p, w)
        elif metric == seuclidean:
            for i in xrange(0, mA):
                for j in xrange(0, mB):
                    dm[i, j] = seuclidean(XA[i, :], XB[j, :], V)
        elif metric == mahalanobis:
            for i in xrange(0, mA):
                for j in xrange(0, mB):
                    dm[i, j] = mahalanobis(XA[i, :], XB[j, :], V)
        else:
            for i in xrange(0, mA):
                for j in xrange(0, mB):
                    dm[i, j] = metric(XA[i, :], XB[j, :])
    elif isinstance(metric, string_types):
        mstr = metric.lower()

        #if XA.dtype != np.double and \
        #       (mstr != 'hamming' and mstr != 'jaccard'):
        #    TypeError('A double array must be passed.')
        if mstr in set(['euclidean', 'euclid', 'eu', 'e']):
            _distance_wrap.cdist_euclidean_wrap(_convert_to_double(XA),
                                                _convert_to_double(XB), dm)
        elif mstr in set(['sqeuclidean', 'sqe', 'sqeuclid']):
            _distance_wrap.cdist_sqeuclidean_wrap(_convert_to_double(XA),
                                                _convert_to_double(XB), dm)
        elif mstr in set(['cityblock', 'cblock', 'cb', 'c']):
            _distance_wrap.cdist_city_block_wrap(_convert_to_double(XA),
                                                 _convert_to_double(XB), dm)
        elif mstr in set(['hamming', 'hamm', 'ha', 'h']):
            if XA.dtype == bool:
                _distance_wrap.cdist_hamming_bool_wrap(_convert_to_bool(XA),
                                                       _convert_to_bool(XB),
                                                       dm)
            else:
                _distance_wrap.cdist_hamming_wrap(_convert_to_double(XA),
                                                  _convert_to_double(XB), dm)
        elif mstr in set(['jaccard', 'jacc', 'ja', 'j']):
            if XA.dtype == bool:
                _distance_wrap.cdist_jaccard_bool_wrap(_convert_to_bool(XA),
                                                       _convert_to_bool(XB),
                                                       dm)
            else:
                _distance_wrap.cdist_jaccard_wrap(_convert_to_double(XA),
                                                  _convert_to_double(XB), dm)
        elif mstr in set(['chebychev', 'chebyshev', 'cheby', 'cheb', 'ch']):
            _distance_wrap.cdist_chebyshev_wrap(_convert_to_double(XA),
                                                _convert_to_double(XB), dm)
        elif mstr in set(['minkowski', 'mi', 'm', 'pnorm']):
            _distance_wrap.cdist_minkowski_wrap(_convert_to_double(XA),
                                                _convert_to_double(XB), dm, p)
        elif mstr in set(['wminkowski', 'wmi', 'wm', 'wpnorm']):
            _distance_wrap.cdist_weighted_minkowski_wrap(_convert_to_double(XA),
                                                         _convert_to_double(XB),
                                                         dm, p,
                                                         _convert_to_double(w))
        elif mstr in set(['seuclidean', 'se', 's']):
            if V is not None:
                V = np.asarray(V, order='c')
                if type(V) != np.ndarray:
                    raise TypeError('Variance vector V must be a numpy array')
                if V.dtype != np.double:
                    raise TypeError('Variance vector V must contain doubles.')
                if len(V.shape) != 1:
                    raise ValueError('Variance vector V must be '
                                     'one-dimensional.')
                if V.shape[0] != n:
                    raise ValueError('Variance vector V must be of the same '
                                     'dimension as the vectors on which the '
                                     'distances are computed.')
                # The C code doesn't do striding.
                [VV] = _copy_arrays_if_base_present([_convert_to_double(V)])
            else:
                X = np.vstack([XA, XB])
                VV = np.var(X, axis=0, ddof=1)
                X = None
                del X
            _distance_wrap.cdist_seuclidean_wrap(_convert_to_double(XA),
                                                 _convert_to_double(XB), VV, dm)
        elif mstr in set(['cosine', 'cos']):
            _cosine_cdist(XA, XB, dm)
        elif mstr in set(['correlation', 'co']):
            XA = np.array(XA, dtype=np.double, copy=True)
            XB = np.array(XB, dtype=np.double, copy=True)
            XA -= XA.mean(axis=1)[:, np.newaxis]
            XB -= XB.mean(axis=1)[:, np.newaxis]
            _cosine_cdist(XA, XB, dm)
        elif mstr in set(['mahalanobis', 'mahal', 'mah']):
            if VI is not None:
                VI = _convert_to_double(np.asarray(VI, order='c'))
                if type(VI) != np.ndarray:
                    raise TypeError('VI must be a numpy array.')
                if VI.dtype != np.double:
                    raise TypeError('The array must contain 64-bit floats.')
                [VI] = _copy_arrays_if_base_present([VI])
            else:
                m = mA + mB
                if m <= n:
                    # There are fewer observations than the dimension of
                    # the observations.
                    raise ValueError("The number of observations (%d) is too "
                                     "small; the covariance matrix is "
                                     "singular. For observations with %d "
                                     "dimensions, at least %d observations "
                                     "are required." % (m, n, n + 1))
                X = np.vstack([XA, XB])
                V = np.atleast_2d(np.cov(X.T))
                X = None
                del X
                VI = _convert_to_double(np.linalg.inv(V).T.copy())
            # (u-v)V^(-1)(u-v)^T
            _distance_wrap.cdist_mahalanobis_wrap(_convert_to_double(XA),
                                                  _convert_to_double(XB),
                                                  VI, dm)
        elif mstr == 'canberra':
            _distance_wrap.cdist_canberra_wrap(_convert_to_double(XA),
                                               _convert_to_double(XB), dm)
        elif mstr == 'braycurtis':
            _distance_wrap.cdist_bray_curtis_wrap(_convert_to_double(XA),
                                                  _convert_to_double(XB), dm)
        elif mstr == 'yule':
            _distance_wrap.cdist_yule_bool_wrap(_convert_to_bool(XA),
                                                _convert_to_bool(XB), dm)
        elif mstr == 'matching':
            _distance_wrap.cdist_matching_bool_wrap(_convert_to_bool(XA),
                                                    _convert_to_bool(XB), dm)
        elif mstr == 'kulsinski':
            _distance_wrap.cdist_kulsinski_bool_wrap(_convert_to_bool(XA),
                                                     _convert_to_bool(XB), dm)
        elif mstr == 'dice':
            _distance_wrap.cdist_dice_bool_wrap(_convert_to_bool(XA),
                                                _convert_to_bool(XB), dm)
        elif mstr == 'rogerstanimoto':
            _distance_wrap.cdist_rogerstanimoto_bool_wrap(_convert_to_bool(XA),
                                                          _convert_to_bool(XB),
                                                          dm)
        elif mstr == 'russellrao':
            _distance_wrap.cdist_russellrao_bool_wrap(_convert_to_bool(XA),
                                                      _convert_to_bool(XB), dm)
        elif mstr == 'sokalmichener':
            _distance_wrap.cdist_sokalmichener_bool_wrap(_convert_to_bool(XA),
                                                         _convert_to_bool(XB),
                                                         dm)
        elif mstr == 'sokalsneath':
            _distance_wrap.cdist_sokalsneath_bool_wrap(_convert_to_bool(XA),
                                                       _convert_to_bool(XB),
                                                       dm)
        elif metric == 'test_euclidean':
            dm = cdist(XA, XB, euclidean)
        elif metric == 'test_seuclidean':
            if V is None:
                V = np.var(np.vstack([XA, XB]), axis=0, ddof=1)
            else:
                V = np.asarray(V, order='c')
            dm = cdist(XA, XB, lambda u, v: seuclidean(u, v, V))
        elif metric == 'test_sqeuclidean':
            dm = cdist(XA, XB, lambda u, v: sqeuclidean(u, v))
        elif metric == 'test_braycurtis':
            dm = cdist(XA, XB, braycurtis)
        elif metric == 'test_mahalanobis':
            if VI is None:
                X = np.vstack([XA, XB])
                V = np.cov(X.T)
                VI = np.linalg.inv(V)
                X = None
                del X
            else:
                VI = np.asarray(VI, order='c')
            [VI] = _copy_arrays_if_base_present([VI])
            # (u-v)V^(-1)(u-v)^T
            dm = cdist(XA, XB, (lambda u, v: mahalanobis(u, v, VI)))
        elif metric == 'test_canberra':
            dm = cdist(XA, XB, canberra)
        elif metric == 'test_cityblock':
            dm = cdist(XA, XB, cityblock)
        elif metric == 'test_minkowski':
            dm = cdist(XA, XB, minkowski, p=p)
        elif metric == 'test_wminkowski':
            dm = cdist(XA, XB, wminkowski, p=p, w=w)
        elif metric == 'test_correlation':
            dm = cdist(XA, XB, correlation)
        elif metric == 'test_hamming':
            dm = cdist(XA, XB, hamming)
        elif metric == 'test_jaccard':
            dm = cdist(XA, XB, jaccard)
        elif metric == 'test_chebyshev' or metric == 'test_chebychev':
            dm = cdist(XA, XB, chebyshev)
        elif metric == 'test_yule':
            dm = cdist(XA, XB, yule)
        elif metric == 'test_matching':
            dm = cdist(XA, XB, matching)
        elif metric == 'test_dice':
            dm = cdist(XA, XB, dice)
        elif metric == 'test_kulsinski':
            dm = cdist(XA, XB, kulsinski)
        elif metric == 'test_rogerstanimoto':
            dm = cdist(XA, XB, rogerstanimoto)
        elif metric == 'test_russellrao':
            dm = cdist(XA, XB, russellrao)
        elif metric == 'test_sokalsneath':
            dm = cdist(XA, XB, sokalsneath)
        elif metric == 'test_sokalmichener':
            dm = cdist(XA, XB, sokalmichener)
        else:
            raise ValueError('Unknown Distance Metric: %s' % mstr)
    else:
        raise TypeError('2nd argument metric must be a string identifier '
                        'or a function.')
    return dm
