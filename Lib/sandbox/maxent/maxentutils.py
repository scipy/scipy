"""maxentutils.py: Utility routines for the maximum entropy module.  Most
of them are either Python replacements for the corresponding Fortran
routines or wrappers around matrices to allow the maxent module to
manipulate ndarrays, scipy sparse matrices, and PySparse matrices a
common interface.

Perhaps the logsumexp() function belongs under the utils/ branch where other
modules can access it more easily.

Copyright: Ed Schofield, 2003-2006
License: BSD-style (see LICENSE.txt in main source directory)
"""

__author__ = "Ed Schofield"
__version__ = '2.0-alpha3'

from __future__ import division
import random, math, bisect, cmath
import numpy
from numpy import log, exp, asarray
from scipy import sparse

def logsumexp(a):
    """Compute the log of the sum of exponentials log(e^{a_1}+...e^{a_n})
    of the components of the array a, avoiding numerical overflow.
    """
    a = asarray(a)
    a_max = a.max()
    return a_max + log((exp(a-a_max)).sum())


def _logsumexpcomplex(values):
    """A version of logsumexp that should work if the values passed are
    complex-numbered, such as the output of robustarraylog().  So we
    expect:

    cmath.exp(logsumexpcomplex(robustarraylog(values))) ~= sum(values)

    except for a small rounding error in both real and imag components.
    The output is complex.  (To recover just the real component, use
    A.real, where A is the complex return value.)
    """
    if len(values) == 0:
        return 0.0
    iterator = iter(values)
    # Get the first element
    while True:
        # Loop until we have a value greater than -inf
        try:
            b_i = iterator.next() + 0j
        except StopIteration:
            # empty
            return float('-inf')
        if b_i.real != float('-inf'):
            break
    
    # Now the rest
    for a_i in iterator:
        a_i += 0j
        if b_i.real > a_i.real:
            increment = robustlog(1.+cmath.exp(a_i - b_i))
            # print "Increment is " + str(increment)
            b_i = b_i + increment
        else:
            increment = robustlog(1.+cmath.exp(b_i - a_i))
            # print "Increment is " + str(increment)
            b_i = a_i + increment
            
    return b_i
    

def logsumexp_naive(values):
    """For testing logsumexp().  Subject to numerical overflow for large
    values (e.g. 720).
    """

    s = 0.0
    for x in values:
        s += math.exp(x)
    return math.log(s)


def robustlog(x):
    """Returns log(x) if x > 0, the complex log cmath.log(x) if x < 0,
    or float('-inf') if x == 0.
    """
    if x == 0.:
        return float('-inf')
    elif type(x) is complex or (type(x) is float and x < 0):
        return cmath.log(x)
    else:
        return math.log(x)


def _robustarraylog(x):
    """ An array version of robustlog.  Operates on a real array x.
    """
    arraylog = emptyarray(len(x), numpy.Complex64)
    for i in range(len(x)):
        xi = x[i]
        if xi > 0:
            arraylog[i] = math.log(xi)
        elif xi == 0.:
            arraylog[i] = float('-inf')
        else:
            arraylog[i] = cmath.log(xi)
    return arraylog

#try:
#    from logsumexp import logsumexp, logsumexpcomplex, robustarraylog
#except:
#    print "Warning: could not load the fast FORTRAN library for logsumexp()."
#    logsumexp = _logsumexp
#    logsumexpcomplex = _logsumexpcomplex
#    robustarraylog = _robustarraylog
#    pass


def arrayexp(x):
    """Returns the elementwise antilog of the real array x.  We try to
    exponentiate with numpy.exp() and, if that fails, with python's
    math.exp().  numpy.exp() is about 10 times faster but throws an
    OverflowError exception for numerical underflow (e.g. exp(-800),
    whereas python's math.exp() just returns zero, which is much more
    helpful.
    """
    try:
        ex = numpy.exp(x)
    except OverflowError:
        print "Warning: OverflowError using numpy.exp(). Using slower Python"\
              " routines instead!"
        ex = numpy.empty(len(x), float)
        for j in range(len(x)):
            ex[j] = math.exp(x[j])
    return ex
 
def arrayexpcomplex(x):
    """Returns the elementwise antilog of the vector x.  We try to
    exponentiate with numpy.exp() and, if that fails, with python's
    math.exp().  numpy.exp() is about 10 times faster but throws an
    OverflowError exception for numerical underflow (e.g. exp(-800),
    whereas python's math.exp() just returns zero, which is much more
    helpful.
    """
    try:
        ex = numpy.exp(x).real
    except OverflowError:
        ex = numpy.empty(len(x),numpy.Float64)
        try:
            for j in range(len(x)):
                ex[j] = math.exp(x[j])
        except TypeError:
            # Perhaps x[j] is complex.  If so, try using the complex
            # exponential and returning the real part.
            for j in range(len(x)):
                ex[j] = cmath.exp(x[j]).real
    return ex
 

def sample_wr(population, k):
    """Chooses k random elements (with replacement) from a population.
    (From the Python Cookbook).
    """
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [population[_int(_random() * n)] for i in xrange(k)]


def densefeatures(f, x):
    """Returns a 1xN dense matrix of 64-bit floats of non-zero
    evaluations of the functions fi in the list f at the point x.
    """
    
    return numpy.array([fi(x) for fi in f])

def densefeaturematrix(f, sample):
    """Returns an (m x n) dense matrix of non-zero evaluations of the
    scalar functions fi in the list f at the points x_1,...,x_n in the
    list sample.
    """
    
    # Was: return numpy.array([[fi(x) for fi in f] for x in sample])

    m = len(f)
    n = len(sample)
    
    F = numpy.empty((m, n), numpy.Float64)
    for i in xrange(m):
        f_i = f[i]
        for i in xrange(m):
            x = sample[j]
            F[i,j] = f_i(x)

     #for j in xrange(n):
     #   x = sample[j]
     #   for i in xrange(m):
     #       F[j,i] = f[i](x)
            
    return F


def sparsefeaturematrix(f, sample, format='ll_mat'):
    """Returns an MxN sparse matrix of non-zero evaluations of the scalar
    or vector functions fi in the list f at the points x_1,...,x_N in the
    sequence 'sample'.
    """

    M = len(f)
    N = len(sample)
    if format == 'll_mat':
        import spmatrix
        sparseF = spmatrix.ll_mat(M, N)
    elif format in ('dok_matrix', 'csc_matrix', 'csr_matrix'):
        sparseF = sparse.dok_matrix((M,N))
        sparseF._validate = False               # speed hack
    else:
        raise ValueError, "sparse matrix format not recognized"

    for i in xrange(M):
        f_i = f[i]
        for j in xrange(N):
            x = sample[j]
            f_i_x = f_i(x)
            if f_i_x != 0:
                sparseF[i,j] = f_i_x

     #for j in xrange(N):
     #   x = sample[j]
     #   for i in xrange(M):
     #       f_i_x = f[i](x)
     #       if f_i_x != 0:
     #           sparseF[j,i] = f_i_x
                    
    #except TypeError:
    #    raise
    #    # The error was (probably) because some fi are not scalar functions
    #    # but vectors.
    #    
    #    for i in xrange(M):
    #        sparseF_i = f[i](sample)
    #        for j in range(len(sparseF_i)):
    #            if sparseF_i[j]:
    #                sparseF[j, i] = sparseF_i[j]
        
    if format == 'csc_matrix':
        return sparseF.tocsc()
    elif format == 'csr_matrix':
        return sparseF.tocsr()
    else:
        return sparseF



def dotprod(u,v):
    """This is a wrapper around general dense or sparse dot products.
    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.
    
    Returns the dot product of the (1xM) sparse array u with the (Mx1)
    (dense) numpy array v.
    """
    #print "Taking the dot product u.v, where"
    #print "u has shape " + str(u.shape)
    #print "v = " + str(v)
    
    try:
        dotprod = numpy.array([0.0])  # a 1x1 array.  Required by spmatrix.
        u.matvec(v, dotprod)
        return dotprod[0]		# extract the scalar
    except AttributeError:
        # Assume u is a dense array.
        return numpy.dot(u,v)



def innerprod(u,v):
    """This is a wrapper around general dense or sparse dot products.
    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.
    
    Returns the inner product of the (MxN) sparse matrix u with the
    (Nx1) (sparse or dense) vector v.  This is a wrapper for u.dot(v) for
    dense arrays and spmatrix objects, and for u.matvec(v, result) for
    PySparse matrices.
    """
    
    # We assume u is sparse.
    (M, N) = u.shape
    vshape = v.shape
    try:
        (P,) = vshape
    except ValueError:
        (P,Q) = vshape
        # if Q != 1:
        #     raise TypeError, "second argument must have shape (N,) rather than (N,P)"
    if N != P:
        raise TypeError, "matrix dimensions are incompatible"
    try:
        v.matvec
    except AttributeError:
        # It looks like v is dense
        try:
            # See if u is sparse
            u.matvec
        except AttributeError:
            # It looks like u is dense
            return numpy.dot(u,v)
        else:
            # Assume u is sparse
            if sparse.isspmatrix(u):
                innerprod = u.matvec(v)   # This returns a float32 type. Why???
                return innerprod
            else:
                # Assume PySparse format
                innerprod = numpy.empty(M, numpy.Float64)
                u.matvec(v, innerprod)
                return innerprod
    else:
        # v looks like it's sparse, so we take the inner prod with
        # spmatrix.matrixmultiply().
        try:
            return spmatrix.matrixmultiply(u,v)
        except:
            raise TypeError, "can't discern correct types for the inner product.  If you are trying to multiply a dense matrix by a sparse vector, this is unsupported."


def innerprodtranspose(U,V):
    """This is a wrapper around general dense or sparse dot products.
    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.
    
    Computes U^T V, where U is a dense or sparse matrix and V is a numpy
    array.  If U is sparse, V must be a rank-1 array, not a matrix.  This
    function is efficient for large matrices U.  This is a wrapper for
    u.T.dot(v) for dense arrays and spmatrix objects, and for
    u.matvec_transp(v, result) for pysparse matrices.
    """

    (m, n) = U.shape
    #pdb.set_trace()
    try:
        # See if U is a PySparse matrix
        U.matvec_transp
    except AttributeError:
        # See if U is a scipy.sparse.spmatrix
        if sparse.isspmatrix(U):
            innerprodtranspose = U.rmatvec(V)
            return innerprodtranspose
        else:
            # Assume U is a dense matrix
            if isinstance(V, numpy.ndarray):
                # V is also dense
                if len(V.shape) == 1:
                    # We can't transpose a rank-1 matrix into a row vector, so
                    # we 'reshape' it.
                    Vm = V.shape[0]
                    Vcolumn = numpy.reshape(V, (1, Vm))
                    x = numpy.dot(Vcolumn, U)
                    return numpy.reshape(x, (n,))
                else:
                    #(Vm, Vn) = V.shape
                    # Assume Vm == m
                    x = numpy.dot(numpy.transpose(V), U)
                    return numpy.transpose(x)
            else:
                raise TypeError, "V must be a dense array"
    else:
        # U looks like a PySparse matrix
        if len(V.shape) == 1:
            innerprod = numpy.empty(n, numpy.Float64)
            U.matvec_transp(V, innerprod)
        else:
            raise TypeError, "innerprodtranspose(U,V) requires that V be " \
                    "a vector (rank-1 dense array) if U is sparse."
        return innerprod



def columnmeans(A):
    """This is a wrapper for general dense or sparse dot products.
    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.

    Returns a dense (1xN) vector with the column averages of A, which can
    be an MxN sparse or dense matrix or a list of M (1xN) sparse
    matrices.
    
    >>> a = numpy.array([[1,2],[3,4]],'d')
    >>> columnmeans(a)
    array([ 2.,  3.])
    """
    if type(A) is numpy.ndarray:
        return numpy.average(A)
    else:
        # Assume it's sparse
        try:
            m = A.shape[0]
        except AttributeError:
            raise TypeError, \
                    "columnmeans() only works with sparse and dense arrays"
        columnsum = innerprodtranspose(A, numpy.ones(m))
        return columnsum / float(m)

def columnvariances(A):
    """This is a wrapper for general dense or sparse dot products.  It
    is not necessary except as a common interface for supporting ndarray,
    scipy spmatrix, and PySparse arrays.

    Returns a dense (1xN) vector with unbiased estimators for the column
    variances for each column of the MxN sparse or dense matrix A.  (The
    normalization is by (M-1).)
    
    >>> a = numpy.array([[1,2],[3,4]],'d')
    >>> columnvariances(a)
    array([ 2.,  2.])
    """
    if type(A) is numpy.ndarray:
        return numpy.std(A,0)**2
    else:
        try:
            m = A.shape[0]
        except AttributeError:
            raise TypeError, \
                    "columnvariances() only works with sparse and dense arrays"
        means = columnmeans(A)
        return columnmeans((A-means)**2) * (m/(m-1.0))

def var(A):
    """Returns an unbiased estimator of the variance of a (1xN) vector or
    list.

    Examples:
    >>> A = numpy.array([1.,1.])
    >>> var(A) == 0
    True
    >>> B = numpy.array([1.])
    >>> var(B) == 0
    True
    >>> C = numpy.array([1.,2.,3.])
    >>> var(C) == 1
    True
    >>> D = numpy.array(range(100))
    >>> abs(var(D) - (841+2./3)) < 1e-10
    True
    >>> E = range(5)   # test with a list
    >>> var(E)
    2.5
    """
    try:
        mean = numpy.average(A)
        return sum((A-mean)**2) / (A.shape[0] - 1.0)
    except ZeroDivisionError:
        # A has only one element
        return 0.0
    except AttributeError:
        # Perhaps A is a list, not an array
        return sum((numpy.array(A)-mean)**2) / (len(A) - 1.0)



class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class DivergenceError(Error):
    """Exception raised if the entropy dual has no finite minimum.
    """
    def __init__(self, message):
        self.message = message
        Error.__init__(self)
        
    def __str__(self):
        return repr(self.message)
    
def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

