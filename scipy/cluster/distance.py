"""
Distance matrix computation from a collection of raw observation vectors

 pdist              computes distances between each observation pair.

Distance functions between two vectors u and v

 braycurtis         the Bray-Curtis distance.
 canberra           the Canberra distance.
 chebyshev          the Chebyshev distance.
 cityblock          the Manhattan distance.
 correlation        the Correlation distance.
 cosine             the Cosine distance.
 dice               the Dice dissimilarity (boolean).
 euclidean          the Euclidean distance.
 hamming            the Hamming distance (boolean).
 jaccard            the Jaccard distance (boolean).
 kulsinski          the Kulsinski distance (boolean).
 mahalanobis        the Mahalanobis distance.
 matching           the matching dissimilarity (boolean).
 minkowski          the Minkowski distance.
 rogerstanimoto     the Rogers-Tanimoto dissimilarity (boolean).
 russellrao         the Russell-Rao dissimilarity (boolean).
 seuclidean         the normalized Euclidean distance.
 sokalmichener      the Sokal-Michener dissimilarity (boolean).
 sokalsneath        the Sokal-Sneath dissimilarity (boolean).
 sqeuclidean        the squared Euclidean distance.
 yule               the Yule dissimilarity (boolean).

Copyright (C) Damian Eads, 2007-2008. New BSD License.

"""

import numpy as np
import _distance_wrap
import types

def _copy_array_if_base_present(a):
    """
    Copies the array if its base points to a parent array.
    """
    if a.base is not None:
        return a.copy()
    elif np.issubsctype(a, np.float32):
        return array(a, dtype=np.double)
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
    if X.dtype != np.bool:
        X = np.bool_(X)
    if not X.flags.contiguous:
        X = X.copy()
    return X

def _convert_to_double(X):
    if X.dtype != np.double:
        X = np.double(X)
    if not X.flags.contiguous:
        X = X.copy()
    return X

def minkowski(u, v, p):
    """
    d = minkowski(u, v, p)

      Returns the Minkowski distance between two vectors u and v,

        ||u-v||_p = (\sum {|u_i - v_i|^p})^(1/p).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if p < 1:
        raise ValueError("p must be at least 1")
    return (abs(u-v)**p).sum() ** (1.0 / p)

def euclidean(u, v):
    """
    d = euclidean(u, v)

      Computes the Euclidean distance between two n-vectors u and v, ||u-v||_2
    """
    u = np.asarray(u)
    v = np.asarray(v)
    q=np.matrix(u-v)
    return np.sqrt((q*q.T).sum())

def sqeuclidean(u, v):
    """
    d = sqeuclidean(u, v)

      Computes the squared Euclidean distance between two n-vectors u and v,
        (||u-v||_2)^2.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return ((u-v)*(u-v).T).sum()

def cosine(u, v):
    """
    d = cosine(u, v)

      Computes the Cosine distance between two n-vectors u and v,
        (1-uv^T)/(||u||_2 * ||v||_2).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return (1.0 - (np.dot(u, v.T) / \
                   (np.sqrt(np.dot(u, u.T)) * np.sqrt(np.dot(v, v.T)))))

def correlation(u, v):
    """
    d = correlation(u, v)

      Computes the correlation distance between two n-vectors u and v,

            1 - (u - n|u|_1)(v - n|v|_1)^T
            --------------------------------- ,
            |(u - n|u|_1)|_2 |(v - n|v|_1)|^T

      where |*|_1 is the Manhattan norm and n is the common dimensionality
      of the vectors.
    """
    umu = u.mean()
    vmu = v.mean()
    um = u - umu
    vm = v - vmu
    return 1.0 - (np.dot(um, vm) /
                  (np.sqrt(np.dot(um, um)) \
                   * np.sqrt(np.dot(vm, vm))))

def hamming(u, v):
    """
    d = hamming(u, v)

      Computes the Hamming distance between two n-vectors u and v,
      which is simply the proportion of disagreeing components in u
      and v. If u and v are boolean vectors, the hamming distance is

         (c_{01} + c_{10}) / n

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return (u != v).mean()

def jaccard(u, v):
    """
    d = jaccard(u, v)

      Computes the Jaccard-Needham dissimilarity between two boolean
      n-vectors u and v, which is

              c_{TF} + c_{FT}
         ------------------------
         c_{TT} + c_{FT} + c_{TF}

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return (np.double(np.bitwise_and((u != v),
                     np.bitwise_or(u != 0, v != 0)).sum())
            /  np.double(np.bitwise_or(u != 0, v != 0).sum()))

def kulsinski(u, v):
    """
    d = kulsinski(u, v)

      Computes the Kulsinski dissimilarity between two boolean n-vectors
      u and v, which is

         c_{TF} + c_{FT} - c_{TT} + n
         ----------------------------
              c_{FT} + c_{TF} + n

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    n = len(u)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)

    return (ntf + nft - ntt + n) / (ntf + nft + n)

def seuclidean(u, v, V):
    """
    d = seuclidean(u, v, V)

      Returns the standardized Euclidean distance between two
      n-vectors u and v. V is a m-dimensional vector of component
      variances. It is usually computed among a larger collection vectors.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    V = np.asarray(V)
    if len(V.shape) != 1 or V.shape[0] != u.shape[0] or u.shape[0] != v.shape[0]:
        raise TypeError('V must be a 1-D array of the same dimension as u and v.')
    return np.sqrt(((u-v)**2 / V).sum())

def cityblock(u, v):
    """
    d = cityblock(u, v)

      Computes the Manhattan distance between two n-vectors u and v,
         \sum {u_i-v_i}.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return abs(u-v).sum()

def mahalanobis(u, v, VI):
    """
    d = mahalanobis(u, v, VI)

      Computes the Mahalanobis distance between two n-vectors u and v,
        (u-v)VI(u-v)^T
      where VI is the inverse covariance matrix.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    VI = np.asarray(VI)
    return np.sqrt(np.dot(np.dot((u-v),VI),(u-v).T).sum())

def chebyshev(u, v):
    """
    d = chebyshev(u, v)

      Computes the Chebyshev distance between two n-vectors u and v,
        \max {|u_i-v_i|}.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return max(abs(u-v))

def braycurtis(u, v):
    """
    d = braycurtis(u, v)

      Computes the Bray-Curtis distance between two n-vectors u and v,
        \sum{|u_i-v_i|} / \sum{|u_i+v_i|}.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return abs(u-v).sum() / abs(u+v).sum()

def canberra(u, v):
    """
    d = canberra(u, v)

      Computes the Canberra distance between two n-vectors u and v,
        \sum{|u_i-v_i|} / \sum{|u_i|+|v_i}.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    return abs(u-v).sum() / (abs(u).sum() + abs(v).sum())

def _nbool_correspond_all(u, v):
    if u.dtype != v.dtype:
        raise TypeError("Arrays being compared must be of the same data type.")

    if u.dtype == np.int or u.dtype == np.float_ or u.dtype == np.double:
        not_u = 1.0 - u
        not_v = 1.0 - v
        nff = (not_u * not_v).sum()
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
        ntt = (u * v).sum()
    elif u.dtype == np.bool:
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
    if u.dtype == np.int or u.dtype == np.float_ or u.dtype == np.double:
        not_u = 1.0 - u
        not_v = 1.0 - v
        nff = (not_u * not_v).sum()
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
        ntt = (u * v).sum()
    else:
        not_u = ~u
        not_v = ~v
        nft = (not_u & v).sum()
        ntf = (u & not_v).sum()
    return (nft, ntf)

def yule(u, v):
    """
    d = yule(u, v)
      Computes the Yule dissimilarity between two boolean n-vectors u and v,

                  R
         ---------------------
         c_{TT} + c_{FF} + R/2

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n, and

         R = 2.0 * (c_{TF} + c_{FT}).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
    print nff, nft, ntf, ntt
    return float(2.0 * ntf * nft) / float(ntt * nff + ntf * nft)

def matching(u, v):
    """
    d = matching(u, v)

      Computes the Matching dissimilarity between two boolean n-vectors
      u and v, which is

         (c_{TF} + c_{FT}) / n

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(nft + ntf) / float(len(u))

def dice(u, v):
    """
    d = dice(u, v)

      Computes the Dice dissimilarity between two boolean n-vectors
      u and v, which is

                c_{TF} + c_{FT}
         ----------------------------
         2 * c_{TT} + c_{FT} + c_{TF}

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if u.dtype == np.bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(ntf + nft) / float(2.0 * ntt + ntf + nft)

def rogerstanimoto(u, v):
    """
    d = rogerstanimoto(u, v)

      Computes the Rogers-Tanimoto dissimilarity between two boolean
      n-vectors u and v,

                  R
         -------------------
         c_{TT} + c_{FF} + R

      where c_{ij} is the number of occurrences of

         u[k] == i and v[k] == j

      for k < n, and

         R = 2.0 * (c_{TF} + c_{FT}).

    """
    u = np.asarray(u)
    v = np.asarray(v)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
    return float(2.0 * (ntf + nft)) / float(ntt + nff + (2.0 * (ntf + nft)))

def russellrao(u, v):
    """
    d = russellrao(u, v)

      Computes the Russell-Rao dissimilarity between two boolean n-vectors
      u and v, (n - c_{TT}) / n where c_{ij} is the number of occurrences
      of u[k] == i and v[k] == j for k < n.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if u.dtype == np.bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    return float(len(u) - ntt) / float(len(u))

def sokalmichener(u, v):
    """
    d = sokalmichener(u, v)

      Computes the Sokal-Michener dissimilarity between two boolean vectors
      u and v, 2R / (S + 2R) where c_{ij} is the number of occurrences of
      u[k] == i and v[k] == j for k < n and R = 2 * (c_{TF} + c{FT}) and
      S = c_{FF} + c_{TT}.
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if u.dtype == np.bool:
        ntt = (u & v).sum()
        nff = (~u & ~v).sum()
    else:
        ntt = (u * v).sum()
        nff = ((1.0 - u) * (1.0 - v)).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(2.0 * (ntf + nft))/float(ntt + nff + 2.0 * (ntf + nft))

def sokalsneath(u, v):
    """
    d = sokalsneath(u, v)

      Computes the Sokal-Sneath dissimilarity between two boolean vectors
      u and v, 2R / (c_{TT} + 2R) where c_{ij} is the number of occurrences
      of u[k] == i and v[k] == j for k < n and R = 2 * (c_{TF} + c{FT}).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if u.dtype == np.bool:
        ntt = (u & v).sum()
    else:
        ntt = (u * v).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(2.0 * (ntf + nft))/float(ntt + 2.0 * (ntf + nft))


def pdist(X, metric='euclidean', p=2, V=None, VI=None):
    """ Y = pdist(X, method='euclidean', p=2)

           Computes the distance between m original observations in
           n-dimensional space. Returns a condensed distance matrix Y.
           For each i and j (i<j), the metric dist(u=X[i], v=X[j]) is
           computed and stored in the ij'th entry. See squareform
           to learn how to retrieve this entry.

        1. Y = pdist(X)

          Computes the distance between m points using Euclidean distance
          (2-norm) as the distance metric between the points. The points
          are arranged as m n-dimensional row vectors in the matrix X.

        2. Y = pdist(X, 'minkowski', p)

          Computes the distances using the Minkowski distance ||u-v||_p
          (p-norm) where p>=1.

        3. Y = pdist(X, 'cityblock')

          Computes the city block or Manhattan distance between the
          points.

        4. Y = pdist(X, 'seuclidean', V=None)

          Computes the standardized Euclidean distance. The standardized
          Euclidean distance between two n-vectors u and v is

            sqrt(\sum {(u_i-v_i)^2 / V[x_i]}).

          V is the variance vector; V[i] is the variance computed over all
          the i'th components of the points. If not passed, it is
          automatically computed.

        5. Y = pdist(X, 'sqeuclidean')

          Computes the squared Euclidean distance ||u-v||_2^2 between
          the vectors.

        6. Y = pdist(X, 'cosine')

          Computes the cosine distance between vectors u and v,

               1 - uv^T
             -----------
             |u|_2 |v|_2

          where |*|_2 is the 2 norm of its argument *.

        7. Y = pdist(X, 'correlation')

          Computes the correlation distance between vectors u and v. This is

            1 - (u - n|u|_1)(v - n|v|_1)^T
            --------------------------------- ,
            |(u - n|u|_1)|_2 |(v - n|v|_1)|^T

          where |*|_1 is the Manhattan (or 1-norm) of its argument *,
          and n is the common dimensionality of the vectors.

        8. Y = pdist(X, 'hamming')

          Computes the normalized Hamming distance, or the proportion
          of those vector elements between two n-vectors u and v which
          disagree. To save memory, the matrix X can be of type boolean.

        9. Y = pdist(X, 'jaccard')

          Computes the Jaccard distance between the points. Given two
          vectors, u and v, the Jaccard distance is the proportion of
          those elements u_i and v_i that disagree where at least one
          of them is non-zero.

        10. Y = pdist(X, 'chebyshev')

          Computes the Chebyshev distance between the points. The
          Chebyshev distance between two n-vectors u and v is the maximum
          norm-1 distance between their respective elements. More
          precisely, the distance is given by

            d(u,v) = max {|u_i-v_i|}.

        11. Y = pdist(X, 'canberra')

          Computes the Canberra distance between the points. The
          Canberra distance between two points u and v is

                      |u_1-v_1|     |u_2-v_2|           |u_n-v_n|
            d(u,v) = ----------- + ----------- + ... + -----------
                     |u_1|+|v_1|   |u_2|+|v_2|         |u_n|+|v_n|

        12. Y = pdist(X, 'braycurtis')

          Computes the Bray-Curtis distance between the points. The
          Bray-Curtis distance between two points u and v is

                     |u_1-v_1| + |u_2-v_2| + ... + |u_n-v_n|
            d(u,v) = ---------------------------------------
                     |u_1+v_1| + |u_2+v_2| + ... + |u_n+v_n|

        13. Y = pdist(X, 'mahalanobis', VI=None)

          Computes the Mahalanobis distance between the points. The
          Mahalanobis distance between two points u and v is
                (u-v)(1/V)(u-v)^T
          where (1/V) is the inverse covariance. If VI is not None,
          VI will be used as the inverse covariance matrix.

        14. Y = pdist(X, 'yule')

          Computes the Yule distance between each pair of boolean
          vectors. (see yule function documentation)

        15. Y = pdist(X, 'matching')

          Computes the matching distance between each pair of boolean
          vectors. (see matching function documentation)

        16. Y = pdist(X, 'dice')

          Computes the Dice distance between each pair of boolean
          vectors. (see dice function documentation)

        17. Y = pdist(X, 'kulsinski')

          Computes the Kulsinski distance between each pair of
          boolean vectors. (see kulsinski function documentation)

        17. Y = pdist(X, 'rogerstanimoto')

          Computes the Rogers-Tanimoto distance between each pair of
          boolean vectors. (see rogerstanimoto function documentation)

        18. Y = pdist(X, 'russellrao')

          Computes the Russell-Rao distance between each pair of
          boolean vectors. (see russellrao function documentation)

        19. Y = pdist(X, 'sokalmichener')

          Computes the Sokal-Michener distance between each pair of
          boolean vectors. (see sokalmichener function documentation)

        20. Y = pdist(X, 'sokalsneath')

          Computes the Sokal-Sneath distance between each pair of
          boolean vectors. (see sokalsneath function documentation)

        21. Y = pdist(X, f)

          Computes the distance between all pairs of vectors in X
          using the user supplied 2-arity function f. For example,
          Euclidean distance between the vectors could be computed
          as follows,

            dm = pdist(X, (lambda u, v: np.sqrt(((u-v)*(u-v).T).sum())))

          Note that you should avoid passing a reference to one of
          the distance functions defined in this library. For example,

            dm = pdist(X, sokalsneath)

          would calculate the pair-wise distances between the vectors
          in X using the Python function sokalsneath. This would result
          in sokalsneath being called {n \choose 2} times, which is
          inefficient. Instead, the optimized C version is more
          efficient, and we call it using the following syntax.

            dm = pdist(X, 'sokalsneath')
       """
#         21. Y = pdist(X, 'test_Y')
#
#           Computes the distance between all pairs of vectors in X
#           using the distance metric Y but with a more succint,
#           verifiable, but less efficient implementation.


    X = np.asarray(X)

    #if np.issubsctype(X, np.floating) and not np.issubsctype(X, np.double):
    #    raise TypeError('Floating point arrays must be 64-bit (got %r).' %
    #    (X.dtype.type,))

    # The C code doesn't do striding.
    [X] = _copy_arrays_if_base_present([_convert_to_double(X)])

    s = X.shape

    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.');

    m = s[0]
    n = s[1]
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)

    mtype = type(metric)
    if mtype is types.FunctionType:
        k = 0
        if metric == minkowski:
            for i in xrange(0, m - 1):
                for j in xrange(i+1, m):
                    dm[k] = minkowski(X[i, :], X[j, :], p)
                    k = k + 1
        elif metric == seuclidean:
            for i in xrange(0, m - 1):
                for j in xrange(i+1, m):
                    dm[k] = seuclidean(X[i, :], X[j, :], V)
                    k = k + 1
        elif metric == mahalanobis:
            for i in xrange(0, m - 1):
                for j in xrange(i+1, m):
                    dm[k] = mahalanobis(X[i, :], X[j, :], V)
                    k = k + 1
        else:
            for i in xrange(0, m - 1):
                for j in xrange(i+1, m):
                    dm[k] = metric(X[i, :], X[j, :])
                    k = k + 1

    elif mtype is types.StringType:
        mstr = metric.lower()

        #if X.dtype != np.double and \
        #       (mstr != 'hamming' and mstr != 'jaccard'):
        #    TypeError('A double array must be passed.')
        if mstr in set(['euclidean', 'euclid', 'eu', 'e']):
            _distance_wrap.pdist_euclidean_wrap(_convert_to_double(X), dm)
        elif mstr in set(['sqeuclidean', 'sqe', 'sqeuclid']):
            _distance_wrap.pdist_euclidean_wrap(_convert_to_double(X), dm)
            dm = dm ** 2.0
        elif mstr in set(['cityblock', 'cblock', 'cb', 'c']):
            _distance_wrap.pdist_city_block_wrap(X, dm)
        elif mstr in set(['hamming', 'hamm', 'ha', 'h']):
            if X.dtype == np.bool:
                _distance_wrap.pdist_hamming_bool_wrap(_convert_to_bool(X), dm)
            else:
                _distance_wrap.pdist_hamming_wrap(_convert_to_double(X), dm)
        elif mstr in set(['jaccard', 'jacc', 'ja', 'j']):
            if X.dtype == np.bool:
                _distance_wrap.pdist_jaccard_bool_wrap(_convert_to_bool(X), dm)
            else:
                _distance_wrap.pdist_jaccard_wrap(_convert_to_double(X), dm)
        elif mstr in set(['chebychev', 'chebyshev', 'cheby', 'cheb', 'ch']):
            _distance_wrap.pdist_chebyshev_wrap(_convert_to_double(X), dm)
        elif mstr in set(['minkowski', 'mi', 'm']):
            _distance_wrap.pdist_minkowski_wrap(_convert_to_double(X), dm, p)
        elif mstr in set(['seuclidean', 'se', 's']):
            if V is not None:
                V = np.asarray(V)
                if type(V) != np.ndarray:
                    raise TypeError('Variance vector V must be a numpy array')
                if V.dtype != np.double:
                    raise TypeError('Variance vector V must contain doubles.')
                if len(V.shape) != 1:
                    raise ValueError('Variance vector V must be one-dimensional.')
                if V.shape[0] != n:
                    raise ValueError('Variance vector V must be of the same dimension as the vectors on which the distances are computed.')
                # The C code doesn't do striding.
                [VV] = _copy_arrays_if_base_present([_convert_to_double(V)])
            else:
                VV = np.var(X, axis=0, ddof=1)
            _distance_wrap.pdist_seuclidean_wrap(_convert_to_double(X), VV, dm)
        # Need to test whether vectorized cosine works better.
        # Find out: Is there a dot subtraction operator so I can
        # subtract matrices in a similar way to multiplying them?
        # Need to get rid of as much unnecessary C code as possible.
        elif mstr in set(['cosine', 'cos']):
            norms = np.sqrt(np.sum(X * X, axis=1))
            _distance_wrap.pdist_cosine_wrap(_convert_to_double(X), dm, norms)
        elif mstr in set(['old_cosine', 'old_cos']):
            norms = np.sqrt(np.sum(X * X, axis=1))
            nV = norms.reshape(m, 1)
            # The numerator u * v
            nm = np.dot(X, X.T)
            # The denom. ||u||*||v||
            de = np.dot(nV, nV.T);
            dm = 1.0 - (nm / de)
            dm[xrange(0,m),xrange(0,m)] = 0.0
            dm = squareform(dm)
        elif mstr in set(['correlation', 'co']):
            X2 = X - X.mean(1)[:,np.newaxis]
            #X2 = X - np.matlib.repmat(np.mean(X, axis=1).reshape(m, 1), 1, n)
            norms = np.sqrt(np.sum(X2 * X2, axis=1))
            _distance_wrap.pdist_cosine_wrap(_convert_to_double(X2), _convert_to_double(dm), _convert_to_double(norms))
        elif mstr in set(['mahalanobis', 'mahal', 'mah']):
            if VI is not None:
                VI = _convert_to_double(np.asarray(VI))
                if type(VI) != np.ndarray:
                    raise TypeError('VI must be a numpy array.')
                if VI.dtype != np.double:
                    raise TypeError('The array must contain 64-bit floats.')
                [VI] = _copy_arrays_if_base_present([VI])
            else:
                V = np.cov(X.T)
                VI = _convert_to_double(np.linalg.inv(V).T.copy())
            # (u-v)V^(-1)(u-v)^T
            _distance_wrap.pdist_mahalanobis_wrap(_convert_to_double(X), VI, dm)
        elif mstr == 'canberra':
            _distance_wrap.pdist_canberra_wrap(_convert_to_double(X), dm)
        elif mstr == 'braycurtis':
            _distance_wrap.pdist_bray_curtis_wrap(_convert_to_bool(X), dm)
        elif mstr == 'yule':
            _distance_wrap.pdist_yule_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'matching':
            _distance_wrap.pdist_matching_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'kulsinski':
            _distance_wrap.pdist_kulsinski_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'dice':
            _distance_wrap.pdist_dice_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'rogerstanimoto':
            _distance_wrap.pdist_rogerstanimoto_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'russellrao':
            _distance_wrap.pdist_russellrao_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'sokalmichener':
            _distance_wrap.pdist_sokalmichener_bool_wrap(_convert_to_bool(X), dm)
        elif mstr == 'sokalsneath':
            _distance_wrap.pdist_sokalsneath_bool_wrap(_convert_to_bool(X), dm)
        elif metric == 'test_euclidean':
            dm = pdist(X, euclidean)
        elif metric == 'test_sqeuclidean':
            if V is None:
                V = np.var(X, axis=0, ddof=1)
            else:
                V = np.asarray(V)
            dm = pdist(X, lambda u, v: seuclidean(u, v, V))
        elif metric == 'test_braycurtis':
            dm = pdist(X, braycurtis)
        elif metric == 'test_mahalanobis':
            if VI is None:
                V = np.cov(X.T)
                VI = np.linalg.inv(V)
            else:
                VI = np.asarray(VI)
            [VI] = _copy_arrays_if_base_present([VI])
            # (u-v)V^(-1)(u-v)^T
            dm = pdist(X, (lambda u, v: mahalanobis(u, v, VI)))
        elif metric == 'test_canberra':
            dm = pdist(X, canberra)
        elif metric == 'test_cityblock':
            dm = pdist(X, cityblock)
        elif metric == 'test_minkowski':
            dm = pdist(X, minkowski, p)
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
        raise TypeError('2nd argument metric must be a string identifier or a function.')
    return dm
