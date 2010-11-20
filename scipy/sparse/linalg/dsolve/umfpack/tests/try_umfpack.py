#!/usr/bin/env python
# Created by: Robert Cimrman, 05.12.2005

"""Benchamrks for umfpack module"""

from optparse import OptionParser
import time
import urllib
import gzip

import numpy as np

import scipy.sparse as sp
import scipy.sparse.linalg.dsolve.umfpack as um
import scipy.linalg as nla

defaultURL = 'http://www.cise.ufl.edu/research/sparse/HBformat/'

usage = """%%prog [options] <matrix file name> [<matrix file name>, ...]

<matrix file name> can be a local or distant (gzipped) file

default url is:
        %s

supported formats are:
        triplet .. [nRow, nCol, nItem] followed by 'nItem' * [ir, ic, value]
        hb      .. Harwell-Boeing format N/A
""" % defaultURL


##
# 05.12.2005, c
def read_triplet( fd ):
    nRow, nCol = map( int, fd.readline().split() )
    nItem = int( fd.readline() )

    ij = np.zeros( (nItem,2), np.int32 )
    val = np.zeros( (nItem,), np.float64 )
    for ii, row in enumerate( fd.readlines() ):
        aux = row.split()
        ij[ii] = int( aux[0] ), int( aux[1] )
        val[ii] = float( aux[2] )

    mtx = sp.csc_matrix( (val, ij), dims = (nRow, nCol), nzmax = nItem )

    return mtx

##
# 06.12.2005, c
def read_triplet2( fd ):
    nRow, nCol = map( int, fd.readline().split() )
    nItem = int( fd.readline() )

    ij, val = io.read_array( fd,
                             columns = [(0,1), (2,)],
                             atype = (np.int32, np.float64),
                             rowsize = nItem )

    mtx = sp.csc_matrix( (val, ij), dims = (nRow, nCol), nzmax = nItem )

    return mtx


formatMap = {'triplet' : read_triplet}
##
# 05.12.2005, c
def readMatrix( matrixName, options ):

    if options.default_url:
        matrixName = defaultURL + matrixName

    print 'url:', matrixName

    if matrixName[:7] == 'http://':
        fileName, status = urllib.urlretrieve( matrixName )
##        print status
    else:
        fileName = matrixName

    print 'file:', fileName

    try:
        readMatrix = formatMap[options.format]
    except:
        raise ValueError('unsupported format: %s' % options.format)

    print 'format:', options.format

    print 'reading...'
    if fileName.endswith('.gz'):
        fd = gzip.open( fileName )
    else:
        fd = open( fileName )

    mtx = readMatrix( fd )

    fd.close()

    print 'ok'

    return mtx

##
# 05.12.2005, c
def main():
    parser = OptionParser( usage = usage )
    parser.add_option( "-c", "--compare",
                       action = "store_true", dest = "compare",
                       default = False,
                       help = "compare with default scipy.sparse solver [default: %default]" )
    parser.add_option( "-p", "--plot",
                       action = "store_true", dest = "plot",
                       default = False,
                       help = "plot time statistics [default: %default]" )
    parser.add_option( "-d", "--default-url",
                       action = "store_true", dest = "default_url",
                       default = False,
                       help = "use default url [default: %default]" )
    parser.add_option( "-f", "--format", type = type( '' ),
                       dest = "format", default = 'triplet',
                       help = "matrix format [default: %default]" )
    (options, args) = parser.parse_args()

    if (len( args ) >= 1):
        matrixNames = args;
    else:
        parser.print_help(),
        return

    sizes, nnzs, times, errors = [], [], [], []
    legends = ['umfpack', 'sparse.solve']
    for ii, matrixName in enumerate( matrixNames ):

        print '*' * 50
        mtx = readMatrix( matrixName, options )

        sizes.append( mtx.shape )
        nnzs.append( mtx.nnz )
        tts = np.zeros( (2,), dtype = np.double )
        times.append( tts )
        err = np.zeros( (2,2), dtype = np.double )
        errors.append( err )

        print 'size              : %s (%d nnz)' % (mtx.shape, mtx.nnz)

        sol0 = np.ones( (mtx.shape[0],), dtype = np.double )
        rhs = mtx * sol0

        umfpack = um.UmfpackContext()

        tt = time.clock()
        sol = umfpack( um.UMFPACK_A, mtx, rhs, autoTranspose = True )
        tts[0] = time.clock() - tt
        print "umfpack           : %.2f s" % tts[0]

        error = mtx * sol - rhs
        err[0,0] = nla.norm( error )
        print '||Ax-b||          :', err[0,0]

        error = sol0 - sol
        err[0,1] = nla.norm( error )
        print '||x - x_{exact}|| :', err[0,1]

        if options.compare:
            tt = time.clock()
            sol = sp.solve( mtx, rhs )
            tts[1] = time.clock() - tt
            print "sparse.solve      : %.2f s" % tts[1]

            error = mtx * sol - rhs
            err[1,0] = nla.norm( error )
            print '||Ax-b||          :', err[1,0]

            error = sol0 - sol
            err[1,1] = nla.norm( error )
            print '||x - x_{exact}|| :', err[1,1]

    if options.plot:
        try:
            import pylab
        except ImportError:
            raise ImportError("could not import pylab")
        times = np.array( times )
        print times
        pylab.plot( times[:,0], 'b-o' )
        if options.compare:
            pylab.plot( times[:,1], 'r-s' )
        else:
            del legends[1]

        print legends

        ax = pylab.axis()
        y2 = 0.5 * (ax[3] - ax[2])
        xrng = range( len( nnzs ) )
        for ii in xrng:
            yy = y2 + 0.4 * (ax[3] - ax[2])\
                 * np.sin( ii * 2 * np.pi / (len( xrng ) - 1) )

            if options.compare:
                pylab.text( ii+0.02, yy,
                            '%s\n%.2e err_umf\n%.2e err_sp'
                            % (sizes[ii], np.sum( errors[ii][0,:] ),
                               np.sum( errors[ii][1,:] )) )
            else:
                pylab.text( ii+0.02, yy,
                            '%s\n%.2e err_umf'
                            % (sizes[ii], np.sum( errors[ii][0,:] )) )
            pylab.plot( [ii, ii], [ax[2], ax[3]], 'k:' )

        pylab.xticks( xrng, ['%d' % (nnzs[ii] ) for ii in xrng] )
        pylab.xlabel( 'nnz' )
        pylab.ylabel( 'time [s]' )
        pylab.legend( legends )
        pylab.axis( [ax[0] - 0.05, ax[1] + 1, ax[2], ax[3]] )
        pylab.show()

if __name__ == '__main__':
    main()
