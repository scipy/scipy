#!/usr/bin/env python 

import numpy as np
np.seterr( all='raise')
from numpy import random 

from expokit import dgpadm, zgpadm, dspadm, zhpadm, \
        dgexpv, zgexpv
from scipy.linalg import expm as pyexpm
from scipy.sparse.base import isspmatrix

minarrsize  = 5
maxarrsize = 50
stride = 5

def expmwrap(A, t=None):
    # dgpadm does an internal check for square matrices
    #  so assume ldh == m
    # constant input variables
    if hasattr( A, 'shape' ):
        m   = A.shape[0]
        ldh = A.shape[1] 
        dtype = A.dtype
    else:
        m   = len(A)
        #ldh = len(A[0])
        dtype = type(A[0][0])
    #t = np.array([1.])[0]
    #ideg = 6
    #lwsp = 4 * m*m + ideg + 1
    # output integers needed to get result.
    iflag = np.array([0])
    # Workspace / Output array
    #wsp = np.zeros( lwsp, dtype=dtype )
    if isspmatrix(A):
        def matvec( v ):
            return A.dot( v )
        #matvec = lambda v : A.dot(v)
        itrace = np.array([0])
        tol = np.array([1e-7])
        # Should use the general Krylov methods here: dgexpv & zgexpv
        if dtype in ('float64', 'float32', float,):
            # w,wsp = dgexpv(v,tol,anorm,matvec,itrace,iflag,[m,t,matvec_extra_args])
            v = np.ones(ldh, dtype=dtype) 
            anorm = np.linalg.norm(A.todense(),np.inf)
            w, wsp = dgexpv(v, tol, anorm, matvec, itrace, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zgexpv( A.todense(), iexph, ns, iflag, t=t )
        print 'done dgexpv'
    else:
        iexph = np.array([0])
        ns    = np.array([0])
        if dtype in ('float64', 'float32', float,):
            wsp = dgpadm( A, iexph, ns, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zgpadm( A, iexph, ns, iflag, t=t )

    if iflag[0] != 0:
        raise IOError( "DGPADM returned error code: {0}".format( iflag[0] ) )
    return np.reshape( wsp[ iexph[0]-1 : iexph[0] + m*m -1 ], (m,m), order='F' )

def sanity_test(array):
    padm = expmwrap( array )
    pypd = pyexpm( array )
    test = np.all(padm == pypd)
    if test is False:
        print "Matrices differ!\n"
        print "EXPOKIT:\n{0}\n\nexpm:\n{1}\n\n{2}".format( padm, pypd, test )
        assert padm.shape == pypd.shape
        for xi,x in enumerate(pypd):
            for yi,y in enumerate(x):
                assert '{0:.10f}'.format(padm[xi][yi]) == '{0:.10f}'.format(y), \
                        'expokit:{0}, expm: {1} at ({2},{3})'.format(padm[xi][yi], y, xi, yi)
    print 'padm result:\n',padm
    print 'pyex result:\n',pypd
    print 'arrays equal to 10 decimal places :)'

def time_dgpadm():
    return expmwrap( array )

def time_pyexpm():
    return pyexpm( array )

def make_arrays(arrsize):
    """Return a real and complex array, using host's default data-types."""
    array = random.rand(arrsize,arrsize)
    carray = np.dot( array, np.eye(arrsize)*np.complex(1+1j) )
    return array, carray

def make_small_arrays(arrsize):
    """Make a real and complex array, using data-types used by default 
    32bit machines."""
    array = np.array( random.rand(arrsize,arrsize), dtype = np.float32 )
    carray = np.array( np.dot( array, np.eye(arrsize)*np.complex(1+1j) ), dtype=np.complex64 )
    return array, carray

def make_sparse_arrays(arrsize=50):
    from scipy.sparse import diags
    print 'diags:-'
    ddiags = [random.rand( arrsize )]
    #ddiags  = np.array([i[:len(i)+1-n] for n,i in enumerate(random.rand(arrsize, arrsize))])
    #ddiags = [random.randint(arrsize)/2 for i in xrange(arrsize) ]
    #ddiags  = random.rand(arrsize, arrsize)
    print ddiags
    #print len( ddiags[0] )
    print 'offsets'
    #offsets = np.array( range( arrsize ), dtype=int )  / 2 
    offsets = np.array([0])
    print offsets
    #print abs( int( offsets[0] ) )
    #data  = np.array( [[1,2,3,4],[1,2,3,4],[1,2,3,4]], dtype=float )
    #diags = np.array( [0, -1, 2], dtype=float )
    array = diags( ddiags, offsets)

    cdata  = np.dot( ddiags,  np.eye(arrsize) * np.complex(1+1j) )
    #cdiags = np.dot( ddiags, np.eye(arrsize-2) * np.complex(1+1j) )
    carray = diags( cdata, offsets )
    return array, carray

if __name__ == '__main__':
    array, carray = make_small_arrays(minarrsize)
    #print "dgpadm small sanity test on dtype:", array.dtype
    #sanity_test(array)
    #print "zgpadm small sanity test on dtype:", carray.dtype
    #sanity_test(carray)

    array, carray = make_arrays(minarrsize)
    #print "dgpadm sanity test on dtype:", array.dtype
    #sanity_test(array)
    #print "zgpadm sanity test on dtype:", carray.dtype
    #sanity_test(carray)

    array, carray = make_sparse_arrays()
    print "dgpadm sparse sanity test on dtype:", array.dtype
    sanity_test(array)
    print "zgpadm sparse sanity test on dtype:", carray.dtype
    try:
        sanity_test(carray)
    except TypeError:
        print carray
        raise

    exit()

    import timeit
    print '{0} {1}'.format( ''.ljust(25), '----- Mean time +- S.D. -----'.rjust(37) )
    print '{0} {1} {2}'.format( 'Array size'.ljust(25), \
            'Expokit'.ljust(25), 'scipy..expm'.ljust(25) )
    for arrsize in xrange( minarrsize, maxarrsize+1, stride ):
        array, carray = make_arrays(arrsize)
        print '{0}'.format(arrsize).ljust(25),

        timer = timeit.Timer("time_dgpadm()", setup="from __main__ import time_dgpadm, np, array")
        times = timer.repeat(number=1000)
        print '{0:.03f} += {1:.03f}'.format( np.mean(times), np.std(times) ).ljust(25),

        timer = timeit.Timer("time_pyexpm()", setup="from __main__ import time_pyexpm, np, array")
        times = timer.repeat(number = 1000)
        print '{0:.03f} += {1:.03f}'.format( np.mean(times), np.std(times) ).ljust(25)

