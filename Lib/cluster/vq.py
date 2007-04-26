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

"""
__docformat__ = 'restructuredtext'

__all__ = ['whiten', 'vq', 'kmeans']


from numpy.random import randint
from numpy import shape, zeros, sqrt, argmin, minimum, array, \
     newaxis, arange, compress, equal, common_type, single, double, take, \
     std, mean
import numpy as N

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
            results = _vq.float_vq(c_obs, c_code_book)
        elif ct is double:
            results = _vq.double_vq(c_obs, c_code_book)
        else:
            results = py_vq(obs, code_book)
    except ImportError:
        results = py_vq(obs, code_book)
    return results

def py_vq(obs, code_book):
    """ Python version of vq algorithm.

    The algorithm simply computes the euclidian distance between each
    observation and every frame in the code_book/

    :Parameters:
        obs : ndarray
            Expect a rank 2 array. Each row is one observation.
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
    (n, d)  = shape(obs)

    # code books and observations should have same number of features
    if not d == code_book.shape[1]:
        raise ValueError("""
            code book(%d) and obs(%d) should have the same 
            number of features (eg columns)""" % (code_book.shape[1], d))
    
    code        = zeros(n, dtype = int)
    min_dist    = zeros(n)
    for i in range(n):
        dist        = N.sum((obs[i] - code_book) ** 2, 1)
        code[i]     = argmin(dist)
        min_dist[i] = dist[code[i]]

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
    
    diff = obs[newaxis, :, :] - code_book[:, newaxis, :]
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
    Nc = code_book.shape[0]
    avg_dist = []
    diff = thresh+1.
    while diff > thresh:
        #compute membership and distances between obs and code_book
        obs_code, distort = vq(obs, code_book)
        avg_dist.append(mean(distort, axis=-1))
        #recalc code_book as centroids of associated obs
        if(diff > thresh):
            has_members = []
            for i in arange(Nc):
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
    """ Generate a code book with minimum distortion.

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

    ("Not checked carefully for accuracy..." he said sheepishly)

    >>> from numpy import array
    >>> from scipy.cluster.vq import vq, kmeans
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

    >>> import RandomArray
    >>> RandomArray.seed(1000,2000)
    >>> codes = 3
    >>> kmeans(whitened,codes)
    (array([[ 2.3110306 ,  2.86287398],
           [ 1.32544402,  0.65607529],
           [ 0.40782893,  2.02786907]]), 0.5196582527686241)

    """
    if int(iter) < 1:
        raise ValueError, 'iter must be >= to 1.'
    if type(k_or_guess) == type(array([])):
        guess   = k_or_guess
        result  = _kmeans(obs, guess, thresh = thresh)
    else:
        #initialize best distance value to a large value
        best_dist = 100000
        No = obs.shape[0]
        k = k_or_guess
        #print 'kmeans iter: ',
        for i in range(iter):
            #the intial code book is randomly selected from observations
            guess       = take(obs, randint(0, No, k), 0)
            book, dist  = _kmeans(obs, guess, thresh = thresh)
            if dist < best_dist:
                best_book = book
                best_dist = dist
        result = best_book, best_dist
    return result
