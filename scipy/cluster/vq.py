""" Vector Quantization Module

    Provides several routines used in creating a code book from a set of
    observations and comparing a set of observations to a code book.

    All routines expect an "observation vector" to be stored in each
    row of the obs matrix.  Similarly the codes are stored row wise
    in the code book matrix.

    whiten(obs) --
        Normalize a group of observations on a per feature basis
    vq(obs,code_book) --
        Calculate code book membership of obs
    kmeans(obs,k_or_guess,iter=20,thresh=1e-5) --
        Train a codebook for mimimum distortion using the kmeans algorithm
    kmeans2
        Similar to kmeans, but with several initialization methods.

"""
__docformat__ = 'restructuredtext'

__all__ = ['whiten', 'vq', 'kmeans', 'kmeans2']

# TODO:
#   - implements high level method for running several times kmeans with
#   different initialialization
#   - warning: what happens if different number of clusters ? For now, emit a
#   warning, but it is not great, because I am not sure it really make sense to
#   succeed in this case (maybe an exception is better ?)
import warnings

from numpy.random import randint
from numpy import shape, zeros, sqrt, argmin, minimum, array, \
     newaxis, arange, compress, equal, common_type, single, double, take, \
     std, mean
import numpy as N

class ClusterError(Exception):
    pass

def whiten(obs):
    """ Normalize a group of observations on a per feature basis.

    Before running kmeans algorithms, it is beneficial to "whiten", or
    scale, the observation data on a per feature basis.  This is done
    by dividing each feature by its standard deviation across all
    observations.

    :Parameters:
        obs : ndarray
            Each row of the array is an observation.  The
            columns are the "features" seen during each observation
            ::

                      #   f0    f1    f2
                obs = [[  1.,   1.,   1.],  #o0
                       [  2.,   2.,   2.],  #o1
                       [  3.,   3.,   3.],  #o2
                       [  4.,   4.,   4.]]) #o3

            XXX perhaps should have an axis variable here.

    :Returns:
        result : ndarray
            Contains the values in obs scaled by the standard devation
            of each column.

    Examples
    --------

    >>> from numpy import array
    >>> from scipy.cluster.vq import whiten
    >>> features  = array([[  1.9,2.3,1.7],
    ...                    [  1.5,2.5,2.2],
    ...                    [  0.8,0.6,1.7,]])
    >>> whiten(features)
    array([[ 3.41250074,  2.20300046,  5.88897275],
           [ 2.69407953,  2.39456571,  7.62102355],
           [ 1.43684242,  0.57469577,  5.88897275]])

    """
    std_dev = std(obs, axis=0)
    return obs / std_dev

def vq(obs, code_book):
    """ Vector Quantization: assign features sets to codes in a code book.

    Vector quantization determines which code in the code book best represents
    an observation of a target.  The features of each observation are compared
    to each code in the book, and assigned the one closest to it.  The
    observations are contained in the obs array. These features should be
    "whitened," or nomalized by the standard deviation of all the features
    before being quantized.  The code book can be created using the kmeans
    algorithm or something similar.

    :Parameters:
        obs : ndarray
            Each row of the array is an observation.  The columns are the
            "features" seen during each observation The features must be
            whitened first using the whiten function or something equivalent.
        code_book : ndarray.
            The code book is usually generated using the kmeans algorithm.
            Each row of the array holds a different code, and the columns are
            the features of the code.

            ::

                            #   f0    f1    f2   f3
                code_book = [[  1.,   2.,   3.,   4.],  #c0
                             [  1.,   2.,   3.,   4.],  #c1
                             [  1.,   2.,   3.,   4.]]) #c2

    :Returns:
        code : ndarray
            If obs is a NxM array, then a length N array is returned that holds
            the selected code book index for each observation.
        dist : ndarray
            The distortion (distance) between the observation and its nearest
            code

    Notes
    -----
    This currently forces 32 bit math precision for speed.  Anyone know
    of a situation where this undermines the accuracy of the algorithm?

    Examples
    --------
    >>> from numpy import array
    >>> from scipy.cluster.vq import vq
    >>> code_book = array([[1.,1.,1.],
    ...                    [2.,2.,2.]])
    >>> features  = array([[  1.9,2.3,1.7],
    ...                    [  1.5,2.5,2.2],
    ...                    [  0.8,0.6,1.7]])
    >>> vq(features,code_book)
    (array([1, 1, 0],'i'), array([ 0.43588989,  0.73484692,  0.83066239]))

    """
    try:
        import _vq
        ct = common_type(obs, code_book)
        c_obs = obs.astype(ct)
        c_code_book = code_book.astype(ct)
        if ct is single:
            results = _vq.vq(c_obs, c_code_book)
        elif ct is double:
            results = _vq.vq(c_obs, c_code_book)
        else:
            results = py_vq(obs, code_book)
    except ImportError:
        results = py_vq(obs, code_book)
    return results

def py_vq(obs, code_book):
    """ Python version of vq algorithm.

    The algorithm simply computes the euclidian distance between each
    observation and every frame in the code_book.

    :Parameters:
        obs : ndarray
            Expects a rank 2 array. Each row is one observation.
        code_book : ndarray
            Code book to use. Same format than obs. Should have same number of
            features (eg columns) than obs.

    :Note:
        This function is slower than the C versions, but it works for
        all input types.  If the inputs have the wrong types for the
        C versions of the function, this one is called as a last resort.

        Its about 20 times slower than the C versions.

    :Returns:
        code : ndarray
            code[i] gives the label of the ith obversation, that its code is
            code_book[code[i]].
        mind_dist : ndarray
            min_dist[i] gives the distance between the ith observation and its
            corresponding code.

    """
    # n = number of observations
    # d = number of features
    if N.ndim(obs) == 1:
        if not N.ndim(obs) == N.ndim(code_book):
            raise ValueError(
                    "Observation and code_book should have the same rank")
        else:
            return _py_vq_1d(obs, code_book)
    else:
        (n, d) = shape(obs)

    # code books and observations should have same number of features and same
    # shape
    if not N.ndim(obs) == N.ndim(code_book):
        raise ValueError("Observation and code_book should have the same rank")
    elif not d == code_book.shape[1]:
        raise ValueError("Code book(%d) and obs(%d) should have the same " \
                         "number of features (eg columns)""" %
                         (code_book.shape[1], d))

    code = zeros(n, dtype=int)
    min_dist = zeros(n)
    for i in range(n):
        dist = N.sum((obs[i] - code_book) ** 2, 1)
        code[i] = argmin(dist)
        min_dist[i] = dist[code[i]]

    return code, sqrt(min_dist)

def _py_vq_1d(obs, code_book):
    """ Python version of vq algorithm for rank 1 only.

    :Parameters:
        obs : ndarray
            Expects a rank 1 array. Each item is one observation.
        code_book : ndarray
            Code book to use. Same format than obs. Should rank 1 too.

    :Returns:
        code : ndarray
            code[i] gives the label of the ith obversation, that its code is
            code_book[code[i]].
        mind_dist : ndarray
            min_dist[i] gives the distance between the ith observation and its
            corresponding code.

    """
    raise RuntimeError("_py_vq_1d buggy, do not use rank 1 arrays for now")
    n = obs.size
    nc = code_book.size
    dist = N.zeros((n, nc))
    for i in range(nc):
        dist[:, i] = N.sum(obs - code_book[i])
    print dist
    code = argmin(dist)
    min_dist = dist[code]

    return code, sqrt(min_dist)

def py_vq2(obs, code_book):
    """2nd Python version of vq algorithm.

    The algorithm simply computes the euclidian distance between each
    observation and every frame in the code_book/

    :Parameters:
        obs : ndarray
            Expect a rank 2 array. Each row is one observation.
        code_book : ndarray
            Code book to use. Same format than obs. Should have same number of
            features (eg columns) than obs.

    :Note:
        This could be faster when number of codebooks is small, but it becomes
        a real memory hog when codebook is large.  It requires NxMxO storage
        where N=number of obs, M = number of features, and O = number of codes.

    :Returns:
        code : ndarray
            code[i] gives the label of the ith obversation, that its code is
            code_book[code[i]].
        mind_dist : ndarray
            min_dist[i] gives the distance between the ith observation and its
            corresponding code.

    """
    d = shape(obs)[1]

    # code books and observations should have same number of features
    if not d == code_book.shape[1]:
        raise ValueError("""
            code book(%d) and obs(%d) should have the same
            number of features (eg columns)""" % (code_book.shape[1], d))

    diff = obs[newaxis, :, :] - code_book[:,newaxis,:]
    dist = sqrt(N.sum(diff * diff, -1))
    code = argmin(dist, 0)
    min_dist = minimum.reduce(dist, 0) #the next line I think is equivalent
                                      #  - and should be faster
    #min_dist = choose(code,dist) # but in practice, didn't seem to make
                                  # much difference.
    return code, min_dist

def _kmeans(obs, guess, thresh=1e-5):
    """ "raw" version of kmeans.

    :Returns:
        code_book :
            the lowest distortion codebook found.
        avg_dist :
            the average distance a observation is from a code in the book.
            Lower means the code_book matches the data better.

    :SeeAlso:
        - kmeans : wrapper around kmeans

    XXX should have an axis variable here.

    Examples
    --------

    Note: not whitened in this example.

    >>> from numpy import array
    >>> from scipy.cluster.vq import _kmeans
    >>> features  = array([[ 1.9,2.3],
    ...                    [ 1.5,2.5],
    ...                    [ 0.8,0.6],
    ...                    [ 0.4,1.8],
    ...                    [ 1.0,1.0]])
    >>> book = array((features[0],features[2]))
    >>> _kmeans(features,book)
    (array([[ 1.7       ,  2.4       ],
           [ 0.73333333,  1.13333333]]), 0.40563916697728591)

    """

    code_book = array(guess, copy = True)
    nc = code_book.shape[0]
    avg_dist = []
    diff = thresh+1.
    while diff > thresh:
        #compute membership and distances between obs and code_book
        obs_code, distort = vq(obs, code_book)
        avg_dist.append(mean(distort, axis=-1))
        #recalc code_book as centroids of associated obs
        if(diff > thresh):
            has_members = []
            for i in arange(nc):
                cell_members = compress(equal(obs_code, i), obs, 0)
                if cell_members.shape[0] > 0:
                    code_book[i] = mean(cell_members, 0)
                    has_members.append(i)
            #remove code_books that didn't have any members
            code_book = take(code_book, has_members, 0)
        if len(avg_dist) > 1:
            diff = avg_dist[-2] - avg_dist[-1]
    #print avg_dist
    return code_book, avg_dist[-1]

def kmeans(obs, k_or_guess, iter=20, thresh=1e-5):
    """Generate a code book with minimum distortion.

    :Parameters:
        obs : ndarray
            Each row of the array is an observation.  The columns are the
            "features" seen during each observation The features must be
            whitened first using the whiten function or something equivalent.
        k_or_guess : int or ndarray
            If integer, it is the number of code book elements.  If a 2D array,
            the array is used as the intial guess for the code book.  The array
            should have k rows, and the same number of columns (features) as
            the obs array.
        iter : int
            The number of times to restart the kmeans algorithm with a new
            initial guess.  If k_or_guess is a 2D array (codebook), this
            argument is ignored and only 1 iteration is run.
        thresh : float
            Terminate each kmeans run when the distortion change from one
            iteration to the next is less than this value.
    :Returns:
        codesbook : ndarray
            The codes that best fit the observation
        distortion : float
            The distortion between the observations and the codes.

    :SeeAlso:
        - kmeans2: similar function, but with more options for initialization,
          and returns label of each observation

    Examples
    --------

    >>> from numpy import array
    >>> from scipy.cluster.vq import vq, kmeans, whiten
    >>> features  = array([[ 1.9,2.3],
    ...                    [ 1.5,2.5],
    ...                    [ 0.8,0.6],
    ...                    [ 0.4,1.8],
    ...                    [ 0.1,0.1],
    ...                    [ 0.2,1.8],
    ...                    [ 2.0,0.5],
    ...                    [ 0.3,1.5],
    ...                    [ 1.0,1.0]])
    >>> whitened = whiten(features)
    >>> book = array((whitened[0],whitened[2]))
    >>> kmeans(whitened,book)
    (array([[ 2.3110306 ,  2.86287398],
           [ 0.93218041,  1.24398691]]), 0.85684700941625547)

    >>> from numpy import random
    >>> random.seed((1000,2000))
    >>> codes = 3
    >>> kmeans(whitened,codes)
    (array([[ 2.3110306 ,  2.86287398],
           [ 1.32544402,  0.65607529],
           [ 0.40782893,  2.02786907]]), 0.5196582527686241)

    """
    if int(iter) < 1:
        raise ValueError, 'iter must be >= to 1.'
    if type(k_or_guess) == type(array([])):
        guess = k_or_guess
        result = _kmeans(obs, guess, thresh = thresh)
    else:
        #initialize best distance value to a large value
        best_dist = 100000
        No = obs.shape[0]
        k = k_or_guess
        #print 'kmeans iter: ',
        for i in range(iter):
            #the intial code book is randomly selected from observations
            guess = take(obs, randint(0, No, k), 0)
            book, dist = _kmeans(obs, guess, thresh = thresh)
            if dist < best_dist:
                best_book = book
                best_dist = dist
        result = best_book, best_dist
    return result

def _kpoints(data, k):
    """Pick k points at random in data (one row = one observation).

    This is done by taking the k first values of a random permutation of 1..N
    where N is the number of observation.

    :Parameters:
        data : ndarray
            Expect a rank 1 or 2 array. Rank 1 are assumed to describe one
            dimensional data, rank 2 multidimensional data, in which case one
            row is one observation.
        k : int
            Number of samples to generate.

    """
    if data.ndim > 1:
        n = data.shape[0]
    else:
        n = data.size

    p = N.random.permutation(n)
    x = data[p[:k], :].copy()

    return x

def _krandinit(data, k):
    """Returns k samples of a random variable which parameters depend on data.

    More precisely, it returns k observations sampled from a Gaussian random
    variable which mean and covariances are the one estimated from data.

    :Parameters:
        data : ndarray
            Expect a rank 1 or 2 array. Rank 1 are assumed to describe one
            dimensional data, rank 2 multidimensional data, in which case one
            row is one observation.
        k : int
            Number of samples to generate.

    """
    mu  = N.mean(data, 0)
    cov = N.atleast_2d(N.cov(data, rowvar = 0))

    # k rows, d cols (one row = one obs)
    # Generate k sample of a random variable ~ Gaussian(mu, cov)
    x = N.random.randn(k, mu.size)
    x = N.dot(x, N.linalg.cholesky(cov).T) + mu

    return x

_valid_init_meth = {'random': _krandinit, 'points': _kpoints}

def _missing_warn():
    """Print a warning when called."""
    warnings.warn("One of the clusters is empty. "
                 "Re-run kmean with a different initialization.")

def _missing_raise():
    """raise a ClusterError when called."""
    raise ClusterError, "One of the clusters is empty. "\
                        "Re-run kmean with a different initialization."

_valid_miss_meth = {'warn': _missing_warn, 'raise': _missing_raise}

def kmeans2(data, k, iter = 10, thresh = 1e-5, minit = 'random',
        missing = 'warn'):
    """Classify a set of points into k clusters using kmean algorithm.

    The algorithm works by minimizing the euclidian distance between data points
    of cluster means. This version is more complete than kmean (has several
    initialisation methods).

    :Parameters:
        data : ndarray
            Expect a rank 1 or 2 array. Rank 1 are assumed to describe one
            dimensional data, rank 2 multidimensional data, in which case one
            row is one observation.
        k : int or ndarray
            Number of clusters. If minit arg is 'matrix', or if a ndarray is
            given instead, it is interpreted as initial cluster to use instead.
        niter : int
            Number of iterations to run.
        thresh : float
            (not used yet).
        minit : string
            Method for initialization. Available methods are random, points and
            uniform:

            random uses k points drawn from a Gaussian random generator which
            mean and variances are estimated from the data.

            points choses k points at random from the points in data.

            uniform choses k points from the data such are they form a uniform
            grid od the dataset (not supported yet).

            matrix means that k has to be interpreted as initial clusters
            (format is the same than data).

    :Returns:
        clusters : ndarray
            the found clusters (one cluster per row).
        label : ndarray
            label[i] gives the label of the ith obversation, that its centroid is
            cluster[label[i]].

    """
    if missing not in _valid_miss_meth.keys():
        raise ValueError("Unkown missing method: %s" % str(missing))
    # If data is rank 1, then we have 1 dimension problem.
    nd  = N.ndim(data)
    if nd == 1:
        d = 1
        #raise ValueError("Input of rank 1 not supported yet")
    elif nd == 2:
        d = data.shape[1]
    else:
        raise ValueError("Input of rank > 2 not supported")

    # If k is not a single value, then it should be compatible with data's
    # shape
    if N.size(k) > 1 or minit == 'matrix':
        if not nd == N.ndim(k):
            raise ValueError("k is not an int and has not same rank than data")
        if d == 1:
            nc = len(k)
        else:
            (nc, dc) = k.shape
            if not dc == d:
                raise ValueError("k is not an int and has not same rank than\
                        data")
        clusters = k.copy()
    else:
        nc = int(k)
        if not nc == k:
            warnings.warn("k was not an integer, was converted.")
        try:
            init = _valid_init_meth[minit]
        except KeyError:
            raise ValueError("unknown init method %s" % str(minit))
        clusters = init(data, k)

    assert not iter == 0
    return _kmeans2(data, clusters, iter, nc, _valid_miss_meth[missing])

def _kmeans2(data, code, niter, nc, missing):
    """ "raw" version of kmeans2. Do not use directly.

    Run kmeans with a given initial codebook.  """
    for i in range(niter):
        # Compute the nearest neighbour for each obs
        # using the current code book
        label = vq(data, code)[0]
        # Update the code by computing centroids using the new code book
        for j in range(nc):
            mbs = N.where(label==j)
            if mbs[0].size > 0:
                code[j] = N.mean(data[mbs], axis=0)
            else:
                missing()

    return code, label

if __name__  == '__main__':
    pass
    #import _vq
    #a = N.random.randn(4, 2)
    #b = N.random.randn(2, 2)

    #print _vq.vq(a, b)
    #print _vq.vq(N.array([[1], [2], [3], [4], [5], [6.]]),
    #        N.array([[2.], [5.]]))
    #print _vq.vq(N.array([1, 2, 3, 4, 5, 6.]), N.array([2., 5.]))
    #_vq.vq(a.astype(N.float32), b.astype(N.float32))
    #_vq.vq(a, b.astype(N.float32))
    #_vq.vq([0], b)
