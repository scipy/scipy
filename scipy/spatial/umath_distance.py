from __future__ import division, absolute_import, print_function

__all__ = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'cosine',
           'correlation', 'dice', 'euclidean', 'hamming', 'jaccard',
           'kulsinski', 'matching', 'rogerstanimoto', 'russellrao',
           'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule',
           'greatcircle']
           
import numpy as np
from scipy.spatial import _umath_distance

def _broadcast_has_zero_strides(u, v):
    """
    Broadcasts the shapes of its two input arrays and returns `True` if
    there are any non matching dimensions of size larger than 1. This
    will translate to a zero stride when the arrays are broadcasted.
    """
    u_shape = (1,) * (v.ndim - u.ndim) + u.shape
    v_shape = (1,) * (u.ndim - v.ndim) + v.shape
    return any(ud != vd and (ud > 1 or vd > 1)
               for ud, vd in zip(u_shape, v_shape))

def _imerge(a, b):
    """
    Merges two iterables into a single generator alternating elements
    from each.
    """
    for i, j in zip(a, b):
        yield i
        yield j

def _cache_array(arr, cache_gufunc):
    """
    Assembles a tuple with an array and the returns of running that array
    through the passed gufunc.
    """
    args = cache_gufunc(arr)
    return (arr,) + ((args,) if cache_gufunc.nout < 2 else args)

def _generic_uv_distance(u, v=None, out=None, cdist=False, dtype=None,
                         dist_gufunc=None, pdist_gufunc=None,
                         cached=False, lowmem=False, cache_gufunc=None,
                         cached_dist_gufunc=None, cached_pdist_gufunc= None):
    """
    Wrapper to provide a common interface for all gufunc distance
    calculations of two input vectors with no additional parameters.

    Parameters
    ----------
    u : (..., N) array_like
        The first input array.
    v : (..., N) array_like, optional
        The second input array. If None, then all pairwise distances
        between vectors in the last two dimensions of `u` will be computed.
    out: (...) ndarray, optional
        Optional output array to store the calculated distance. Must be of
        the right shape
    cdist : bool, optional
        If True, computes the distance between each pair of vectors in the
        last two dimensions of `u` and `v`. Ignored if `v` is None.
    dtype: optional
        Data type to convert non-array inputs to.
    dist_gufunc: gufunc with `(d),(d)->()` signature
        Gufunc that computes the distance on the last dimension of its
        inputs, and broadcasts on the rest.
    pdist_gufunc: gufunc with `(n,d)->(p)` signature
        Gufunc that computes all pairwise distances over the last two
        dimensions of its input, and broadcasts on the rest.
    cached: bool, optional
        If True, then a second set of gufuncs must be provided to precompute
        some data from the input vectors, and to reuse the precomputed data
        when individual vectors are used multiple times, i.e. in cdist and
        pdist type calculations. The cached set of gufuncs are only used if
        there are vectors being reused in the actual calculations, determined
        by analyzing the shape of the input arrays.
    lowmem: bool, optional
        If caching functions are defined, this flag allows the user to
        disable their use. For properly defined gufuncs it should result
        in a slower calculation with a more efficient memory use.
    cache_gufunc: gufunc with `(d)->...` signature, optional
        Gufunc to compute the values to cache. Its (possibly multiple)
        return values will be passed as additional arguments to
        `cached_pdist_gufunc` and `cached_dist_gufunc`, in the latter case
        alternating values cached from `u` and `v`.
    cached_dist_gufunc: gufunc with `(d),(d),...->()` signature, optional
        Gufunc that computes the distance on the last dimension of its first
        two inputs, and broadcasts on the rest, using the cached data
        passed in the additional arguments, which alternates values for the
        first and second argument, e.g. `f(u, v, u_cache1, v_cache1, ...)`.
    cached_pdist_gufunc: gufunc with `(n,d),...->(p)` signature, optional
        Gufunc that computes all pairwise distances over the last two
        dimensions of its first input, and broadcasts on the rest, using
        the cached data passed in the additional arguments.

    Returns
    -------
    distance : (...) ndarray
        The distance between vectors `u` and `v`.
    """
    # Convert `u` to array
    dt = None if isinstance(u, np.ndarray) else dtype
    u = np.asarray(u, dtype=dt)
    scalar_output = u.ndim < 2
    u.shape = (1,) * (2 - u.ndim) + u.shape
    
    if v is None:
        if u.shape[-2] == 1:
            raise ValueError('not enough vectors to do pairwise distances')
        # The pairwise distance gufunc relies on an array of the proper
        # shape to be passed in as a return array with the `out` argument
        p_len = u.shape[-2] * (u.shape[-2] - 1) // 2
        out_shape = u.shape[:-2] + (p_len,)
        if out is None:
            out = np.empty(out_shape, dtype=np.double)
        elif not isinstance(out, np.ndarray):
            raise TypeError('return arrays must be of ArrayType')
        elif out.ndim < 1 or out.shape[-1] != p_len:
            raise ValueError('wrong return array shape for pairwise distance')
        func = pdist_gufunc
        args = (u,)
        if cached and not lowmem:
            args = _cache_array(u, cache_gufunc)
            func = cached_pdist_gufunc
        # numpy 1.5 only supports out as a positional (not keyword) argument
        args = args + (out,)
        return func(*args)
    else:
        dt = None if isinstance(v, np.ndarray) else dtype
        v = np.asarray(v, dtype=dt)
        scalar_output &= v.ndim < 2
        v.shape = (1,) * (2 - v.ndim) + v.shape
        if cdist:
            u = u[..., np.newaxis, :]
            v = v[..., np.newaxis, :, :]
            scalar_output = False
        func = dist_gufunc
        args = (u, v)
        kwargs = {} if out is None else {'out' : out}
        if cached and not lowmem and _broadcast_has_zero_strides(*args):
            args = (_cache_array(j, cache_gufunc) for j in (u, v))
            args = tuple(_imerge(*args))
            func = cached_dist_gufunc
        # numpy 1.5 only supports out as a positional (not keyword) argument
        if out is not None:
            args = args + (out,)
        ret = func(*args)
        if scalar_output:
            ret = float(ret.squeeze())
        return ret
        
def braycurtis(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._braycurtis,
                                _umath_distance._pairwise_braycurtis)

def canberra(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._canberra,
                                _umath_distance._pairwise_canberra)

def chebyshev(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._chebyshev,
                                _umath_distance._pairwise_chebyshev)

def cityblock(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._cityblock,
                                _umath_distance._pairwise_cityblock)

def euclidean(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._euclidean,
                                _umath_distance._pairwise_euclidean)

def sqeuclidean(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._sqeuclidean,
                                _umath_distance._pairwise_sqeuclidean)

def dice(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._dice,
                                _umath_distance._pairwise_dice)

def hamming(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                None,
                                _umath_distance._hamming,
                                _umath_distance._pairwise_hamming)

def jaccard(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._jaccard,
                                _umath_distance._pairwise_jaccard)

def kulsinski(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._kulsinski,
                                _umath_distance._pairwise_kulsinski)

# matching is the same as hamming for bools
def matching(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._hamming,
                                _umath_distance._pairwise_hamming)

def rogerstanimoto(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._rogerstanimoto,
                                _umath_distance._pairwise_rogerstanimoto)

def russellrao(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._russellrao,
                                _umath_distance._pairwise_russellrao)

# sokalmichener is the same as hamming for bools
def sokalmichener(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._hamming,
                                _umath_distance._pairwise_hamming)
                                
def sokalsneath(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._sokalsneath,
                                _umath_distance._pairwise_sokalsneath)

def yule(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.bool,
                                _umath_distance._yule,
                                _umath_distance._pairwise_yule)

def greatcircle(u, v=None, out=None, cdist=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._greatcircle,
                                _umath_distance._pairwise_greatcircle)

def cosine(u, v=None, out=None, cdist=False, lowmem=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._cosine,
                                _umath_distance._pairwise_cosine,
                                True, lowmem,
                                _umath_distance._norm,
                                _umath_distance._cached_cosine,
                                _umath_distance._cached_pairwise_cosine)

def correlation(u, v=None, out=None, cdist=False, lowmem=False):
    return _generic_uv_distance(u, v, out, cdist,
                                np.double,
                                _umath_distance._correlation,
                                _umath_distance._pairwise_correlation,
                                True, lowmem,
                                _umath_distance._mean_ssqdm,
                                _umath_distance._cached_correlation,
                                _umath_distance._cached_pairwise_correlation)


_doc_header = """
    Computes the {0} distance between two arrays of vectors
    using gufuncs.
    """

_doc_parameters = """

    Parameters
    ----------
    u : (..., d) array_like
        Input array.
    v : (..., d) array_like, optional
        Input array. If `None`, then all pairwise distances between vectors
        in the last two dimensions of `u` will be computed.
    out: (...) ndarray, optional
        Optional output array to store the calculated distance. Must be of
        the right shape.
    cdist : bool, optional
        If True, computes the distance between each pair of vectors in the
        last two dimensiond of `u` and `v`. Ignored if `v` is None.
    """

_doc_correlation_lowmem = """
    lowmem: bool, optional
        If True, no attempt will be made to speed up the calculation by
        cacheing the means and sums of squared differences from the mean
        of the vectors. The default False will typically result in faster
        calculations, with a possibly larger memory footprint.
    """

_doc_returns = """    

    Returns
    -------
    {0} : (...) ndarray
        The {1} distance between vectors `u` and `v`.
    """

_doc_cosine_lowmem = """
    lowmem: bool, optional
        If True, no attempt will be made to speed up the calculation by
        cacheing the norms of the vectors. The default False will typically
        result in faster calculations, with a possibly larger memory
        footprint.
    """
                                
_docs = [(braycurtis, 'Bray-Curtis', 'braycurtis', ''),
         (canberra, 'Canberra', 'canberra', ''),
         (chebyshev, 'Chebyshev', 'chebyshev', ''),
         (cityblock, 'city block (Manhattan)', 'cityblock', ''),
         (correlation, 'correlation', 'correlation', _doc_correlation_lowmem),
         (cosine, 'cosine', 'cosine', _doc_cosine_lowmem),
         (dice, 'Dice', 'dice', ''),
         (euclidean, 'Euclidean', 'euclidean', ''),
         (hamming, 'Hamming', 'hamming', ''),
         (jaccard, 'Jaccard', 'jaccard', ''),
         (kulsinski, 'Kulsinski', 'kulsinski', ''),
         (matching, 'matching', 'matching', ''),
         (rogerstanimoto, 'Rogers-Tanimoto', 'rogerstanimoto', ''),
         (russellrao, 'Russell-Rao', 'russellrao', ''),
         (sokalmichener, 'Sokal-Michener', 'sokalmichener', ''),
         (sokalsneath, 'Sokal-Sneath', 'sokalsneath', ''),
         (yule, 'Yule', 'yule', ''),
         (greatcircle, 'great-circle', 'greatcircle', '')]

for _func, _Name, _name, _lowmem in _docs:
    _func.__doc__ = ''.join((_doc_header.format(_Name),
                             _doc_parameters, _lowmem,
                             _doc_returns.format(_name, _Name)))