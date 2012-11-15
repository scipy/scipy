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
            #return w
        itrace = np.array([0])
        tol = np.array([1e-7])
        v = np.ones(m, dtype=dtype) 
        anorm = np.linalg.norm(A.todense(), np.inf)
        # Should use the general Krylov methods here: dgexpv & zgexpv
        if dtype in ('float64', 'float32', float,):
            # w,wsp = dgexpv(v,tol,anorm,matvec,itrace,iflag,[m,t,matvec_extra_args])
            w, wsp = dgexpv(v, tol, anorm, matvec, itrace, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            w, wsp = zgexpv(v * 1+1j, tol, anorm, matvec, itrace, iflag, t=t )
            #wsp = zgexpv( A.todense(), iexph, ns, iflag, t=t )
        #retmat = wsp[m*(m+1)+m+(m+2)**2:]
        retmat = w
    else:
        iexph = np.array([0])
        ns    = np.array([0])
        if dtype in ('float64', 'float32', float,):
            wsp = dgpadm( A, iexph, ns, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zgpadm( A, iexph, ns, iflag, t=t )
        retmat = np.reshape( wsp[ iexph[0]-1 : iexph[0] + m*m -1 ], (m,m), order='F' )

    if iflag[0] != 0:
        raise IOError( "DGPADM returned error code: {0}".format( iflag[0] ) )
    return retmat

def sanity_test(array, precision=10):
    """Call scipy.linalg.expm and expmwrap, above.
    Compare results to specified precision, and raise
    and AssertionError if there are differences"""
    padm = expmwrap( array )
    pypd = pyexpm( array )
    test = np.all(padm == pypd)
    if test == False:
        print "np.all says matrices differ!"
        print 'EXPOKIT result:\n', padm
        print 'SciPy   result:\n', pypd
        assert padm.shape == pypd.shape, 'Shapes differ: EXPOKIT: {0}, SciPy: {1}'\
                .format( padm.shape, pypd.shape )
        for xi,x in enumerate(pypd):
            for yi,xy in enumerate(x):
                #assert '{0:.10f}'.format(padm[xi][yi]) == '{0:.10f}'.format(xy), \
                #        'expokit:{0}, expm: {1} at ({2},{3})'.format(padm[xi][yi], xy, xi, yi)
                assert '{0:.{1}f}'.format(padm[xi][yi], precision) == \
                       '{0:.{1}f}'.format(pypd[xi][yi], precision), \
                        'Matrices differ to {0} decimal places.'.format(precision) + \
                        'expokit:{0}, expm: {1} at ({2},{3})'.format( \
                            padm[xi][yi], pypd[xi][yi], xi, yi )
    print 'arrays equal to {0} decimal places :)'.format( precision )

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
    on 32bit machines."""
    array = np.array( random.rand(arrsize,arrsize), dtype = np.float32 )
    carray = np.array( np.dot( array, np.eye(arrsize)*np.complex(1+1j) ), dtype=np.complex64 )
    return array, carray

def make_sparse_arrays(arrsize=50):
    from scipy.sparse import diags
    # real sparse matrix
    ddiags = [ random.rand( arrsize ),
               random.rand(arrsize-1),
               random.rand(arrsize-1) ]
    offsets = np.array([0, 1, -1])
    array = diags( ddiags, offsets)

    # complex sparse matrix
    cdata  = [ np.dot( ddiags[0],  np.eye(arrsize)   * np.complex(1+1j) ),
               np.dot( ddiags[1],  np.eye(arrsize-1) * np.complex(1+1j) ),
               np.dot( ddiags[1],  np.eye(arrsize-1) * np.complex(1+1j) ) ]
    carray = diags( cdata, offsets )
    return array, carray

if __name__ == '__main__':
    prec = 4
    #array, carray = make_small_arrays(minarrsize)
    #print "\nDGPADM small sanity test on dtype:", array.dtype
    #sanity_test(array, precision=prec)
    #print "\nZGPADM small sanity test on dtype:", carray.dtype
    #sanity_test(carray, precision=prec)

    #array, carray = make_arrays(minarrsize)
    #print "\nDGPADM sanity test on dtype:", array.dtype
    #sanity_test(array, precision=prec)
    #print "\nZGPADM sanity test on dtype:", carray.dtype
    #sanity_test(carray, precision=prec)

    array, carray = make_sparse_arrays()
    print "\nDGEXPV sparse sanity test on dtype:", array.dtype
    sanity_test(array, precision=prec)
    print "\nZGEXPV sparse sanity test on dtype:", carray.dtype
    sanity_test(carray, precision=prec)

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

