"""
Utility routines for the maximum entropy module.

Most of them are either Python replacements for the corresponding Fortran
routines or wrappers around matrices to allow the maxent module to
manipulate ndarrays, scipy sparse matrices, and PySparse matrices a
common interface.

Perhaps the logsumexp() function belongs under the utils/ branch where other
modules can access it more easily.

Copyright: Ed Schofield, 2003-2006
License: BSD-style (see LICENSE.txt in main source directory)

"""

# Future imports must come before any code in 2.5
from __future__ import division

__author__ = "Ed Schofield"
__version__ = '2.0'

import random
import math
import cmath
import numpy
from numpy import log, exp, asarray, ndarray, empty
from scipy import sparse
from scipy.misc import logsumexp


def _logsumexpcomplex(values):
    """A version of logsumexp that should work if the values passed are
    complex-numbered, such as the output of robustarraylog().  So we
    expect:

    cmath.exp(logsumexpcomplex(robustarraylog(values))) ~= sum(values,axis=0)

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
    arraylog = empty(len(x), numpy.complex64)
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
    """
    Returns the elementwise antilog of the real array x.

    We try to exponentiate with numpy.exp() and, if that fails, with
    python's math.exp().  numpy.exp() is about 10 times faster but throws
    an OverflowError exception for numerical underflow (e.g. exp(-800),
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
    """
    Returns the elementwise antilog of the vector x.

    We try to exponentiate with numpy.exp() and, if that fails, with python's
    math.exp().  numpy.exp() is about 10 times faster but throws an
    OverflowError exception for numerical underflow (e.g. exp(-800),
    whereas python's math.exp() just returns zero, which is much more
    helpful.

    """
    try:
        ex = numpy.exp(x).real
    except OverflowError:
        ex = numpy.empty(len(x), float)
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
    """Returns a dense array of non-zero evaluations of the functions fi
    in the list f at the point x.
    """

    return numpy.array([fi(x) for fi in f])

def densefeaturematrix(f, sample):
    """Returns an (m x n) dense array of non-zero evaluations of the
    scalar functions fi in the list f at the points x_1,...,x_n in the
    list sample.
    """

    # Was: return numpy.array([[fi(x) for fi in f] for x in sample])

    m = len(f)
    n = len(sample)

    F = numpy.empty((m, n), float)
    for i in xrange(m):
        f_i = f[i]
        for j in xrange(n):
            x = sample[j]
            F[i,j] = f_i(x)

     #for j in xrange(n):
     #   x = sample[j]
     #   for i in xrange(m):
     #       F[j,i] = f[i](x)

    return F


def sparsefeatures(f, x, format='csc_matrix'):
    """ Returns an Mx1 sparse matrix of non-zero evaluations of the
    scalar functions f_1,...,f_m in the list f at the point x.

    If format='ll_mat', the PySparse module (or a symlink to it) must be
    available in the Python site-packages/ directory.  A trimmed-down
    version, patched for NumPy compatibility, is available in the SciPy
    sandbox/pysparse directory.
    """
    m = len(f)
    if format == 'll_mat':
        import spmatrix
        sparsef = spmatrix.ll_mat(m, 1)
    elif format in ('dok_matrix', 'csc_matrix', 'csr_matrix'):
        sparsef = sparse.dok_matrix((m, 1))

    for i in xrange(m):
        f_i_x = f[i](x)
        if f_i_x != 0:
            sparsef[i, 0] = f_i_x

    if format == 'csc_matrix':
        print "Converting to CSC matrix ..."
        return sparsef.tocsc()
    elif format == 'csr_matrix':
        print "Converting to CSR matrix ..."
        return sparsef.tocsr()
    else:
        return sparsef

def sparsefeaturematrix(f, sample, format='csc_matrix'):
    """Returns an (m x n) sparse matrix of non-zero evaluations of the scalar
    or vector functions f_1,...,f_m in the list f at the points
    x_1,...,x_n in the sequence 'sample'.

    If format='ll_mat', the PySparse module (or a symlink to it) must be
    available in the Python site-packages/ directory.  A trimmed-down
    version, patched for NumPy compatibility, is available in the SciPy
    sandbox/pysparse directory.
    """

    m = len(f)
    n = len(sample)
    if format == 'll_mat':
        import spmatrix
        sparseF = spmatrix.ll_mat(m, n)
    elif format in ('dok_matrix', 'csc_matrix', 'csr_matrix'):
        sparseF = sparse.dok_matrix((m, n))
    else:
        raise ValueError("sparse matrix format not recognized")

    for i in xrange(m):
        f_i = f[i]
        for j in xrange(n):
            x = sample[j]
            f_i_x = f_i(x)
            if f_i_x != 0:
                sparseF[i,j] = f_i_x

    if format == 'csc_matrix':
        return sparseF.tocsc()
    elif format == 'csr_matrix':
        return sparseF.tocsr()
    else:
        return sparseF



def dotprod(u,v):
    """
    This is a wrapper around general dense or sparse dot products.

    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.

    Returns the dot product of the (1 x m) sparse array u with the
    (m x 1) (dense) numpy array v.

    """
    #print "Taking the dot product u.v, where"
    #print "u has shape " + str(u.shape)
    #print "v = " + str(v)

    try:
        dotprod = numpy.array([0.0])  # a 1x1 array.  Required by spmatrix.
        u.matvec(v, dotprod)
        return dotprod[0]               # extract the scalar
    except AttributeError:
        # Assume u is a dense array.
        return numpy.dot(u,v)



def innerprod(A,v):
    """
    This is a wrapper around general dense or sparse dot products.

    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.

    Returns the inner product of the (m x n) dense or sparse matrix A
    with the n-element dense array v.  This is a wrapper for A.dot(v) for
    dense arrays and spmatrix objects, and for A.matvec(v, result) for
    PySparse matrices.

    """

    # We assume A is sparse.
    (m, n) = A.shape
    vshape = v.shape
    try:
        (p,) = vshape
    except ValueError:
        (p, q) = vshape
    if n != p:
        raise TypeError("matrix dimensions are incompatible")
    if isinstance(v, ndarray):
        try:
            # See if A is sparse
            A.matvec
        except AttributeError:
            # It looks like A is dense
            return numpy.dot(A, v)
        else:
            # Assume A is sparse
            if sparse.isspmatrix(A):
                innerprod = A.matvec(v)   # This returns a float32 type. Why???
                return innerprod
            else:
                # Assume PySparse format
                innerprod = numpy.empty(m, float)
                A.matvec(v, innerprod)
                return innerprod
    elif sparse.isspmatrix(v):
        return A * v
    else:
        raise TypeError("unsupported types for inner product")


def innerprodtranspose(A,v):
    """
    This is a wrapper around general dense or sparse dot products.

    It is not necessary except as a common interface for supporting
    ndarray, scipy spmatrix, and PySparse arrays.

    Computes A^T V, where A is a dense or sparse matrix and V is a numpy
    array.  If A is sparse, V must be a rank-1 array, not a matrix.  This
    function is efficient for large matrices A.  This is a wrapper for
    u.T.dot(v) for dense arrays and spmatrix objects, and for
    u.matvec_transp(v, result) for pysparse matrices.

    """

    (m, n) = A.shape
    #pdb.set_trace()
    if hasattr(A, 'matvec_transp'):
        # A looks like a PySparse matrix
        if len(v.shape) == 1:
            innerprod = numpy.empty(n, float)
            A.matvec_transp(v, innerprod)
        else:
            raise TypeError("innerprodtranspose(A,v) requires that v be "
                    "a vector (rank-1 dense array) if A is sparse.")
        return innerprod
    elif sparse.isspmatrix(A):
        return (A.conj().transpose() * v).transpose()
    else:
        # Assume A is dense
        if isinstance(v, numpy.ndarray):
            # v is also dense
            if len(v.shape) == 1:
                # We can't transpose a rank-1 matrix into a row vector, so
                # we reshape it.
                vm = v.shape[0]
                vcolumn = numpy.reshape(v, (1, vm))
                x = numpy.dot(vcolumn, A)
                return numpy.reshape(x, (n,))
            else:
                #(vm, vn) = v.shape
                # Assume vm == m
                x = numpy.dot(numpy.transpose(v), A)
                return numpy.transpose(x)
        else:
            raise TypeError("unsupported types for inner product")


def rowmeans(A):
    """
    This is a wrapper for general dense or sparse dot products.

    It is only necessary as a common interface for supporting ndarray,
    scipy spmatrix, and PySparse arrays.

    Returns a dense (m x 1) vector representing the mean of the rows of A,
    which be an (m x n) sparse or dense matrix.

    >>> a = numpy.array([[1,2],[3,4]], float)
    >>> rowmeans(a)
    array([ 1.5,  3.5])

    """
    if type(A) is numpy.ndarray:
        return A.mean(1)
    else:
        # Assume it's sparse
        try:
            n = A.shape[1]
        except AttributeError:
            raise TypeError("rowmeans() only works with sparse and dense "
                            "arrays")
        rowsum = innerprod(A, numpy.ones(n, float))
        return rowsum / float(n)

def columnmeans(A):
    """
    This is a wrapper for general dense or sparse dot products.

    It is only necessary as a common interface for supporting ndarray,
    scipy spmatrix, and PySparse arrays.

    Returns a dense (1 x n) vector with the column averages of A, which can
    be an (m x n) sparse or dense matrix.

    >>> a = numpy.array([[1,2],[3,4]],'d')
    >>> columnmeans(a)
    array([ 2.,  3.])

    """
    if type(A) is numpy.ndarray:
        return A.mean(0)
    else:
        # Assume it's sparse
        try:
            m = A.shape[0]
        except AttributeError:
            raise TypeError("columnmeans() only works with sparse and dense "
                            "arrays")
        columnsum = innerprodtranspose(A, numpy.ones(m, float))
        return columnsum / float(m)

def columnvariances(A):
    """
    This is a wrapper for general dense or sparse dot products.

    It is not necessary except as a common interface for supporting ndarray,
    scipy spmatrix, and PySparse arrays.

    Returns a dense (1 x n) vector with unbiased estimators for the column
    variances for each column of the (m x n) sparse or dense matrix A.  (The
    normalization is by (m - 1).)

    >>> a = numpy.array([[1,2], [3,4]], 'd')
    >>> columnvariances(a)
    array([ 2.,  2.])

    """
    if type(A) is numpy.ndarray:
        return numpy.std(A,0)**2
    else:
        try:
            m = A.shape[0]
        except AttributeError:
            raise TypeError("columnvariances() only works with sparse "
                            "and dense arrays")
        means = columnmeans(A)
        return columnmeans((A-means)**2) * (m/(m-1.0))

def flatten(a):
    """Flattens the sparse matrix or dense array/matrix 'a' into a
    1-dimensional array
    """
    if sparse.isspmatrix(a):
        return a.A.flatten()
    else:
        return numpy.asarray(a).flatten()

class DivergenceError(Exception):
    """Exception raised if the entropy dual has no finite minimum.
    """
    def __init__(self, message):
        self.message = message
        Exception.__init__(self)

    def __str__(self):
        return repr(self.message)

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
