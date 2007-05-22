"""
Pure SciPy implementation of Locally Optimal Block Preconditioned Conjugate
Gradient Method (LOBPCG), see
http://www-math.cudenver.edu/~aknyazev/software/BLOPEX/

License: BSD

(c) Robert Cimrman, Andrew Knyazev
"""

import numpy as nm
import scipy as sc
import scipy.sparse as sp
import scipy.linalg as la
import scipy.io as io
import types
from symeig import symeig

def pause():
    raw_input()

def save( ar, fileName ):
    io.write_array( fileName, ar, precision = 8 )

##
# 21.05.2007, c
def as2d( ar ):
    if ar.ndim == 2:
        return ar
    else: # Assume 1!
        aux = nm.array( ar, copy = False )
        aux.shape = (ar.shape[0], 1)
        return aux

##
# 05.04.2007, c
# 10.04.2007
def makeOperator( operatorInput, expectedShape ):
    class Operator( object ):
        def __call__( self, vec ):
            return self.call( vec )
    
    operator = Operator()
    operator.obj = operatorInput
    
    if hasattr( operatorInput, 'shape' ):
        operator.shape = operatorInput.shape
        operator.dtype = operatorInput.dtype
        if operator.shape != expectedShape:
            raise ValueError, 'bad operator shape %s != %s' \
                  % (expectedShape, operator.shape)
        if sp.issparse( operatorInput ):
            def call( vec ):
                out = operator.obj * vec
                if sp.issparse( out ):
                    out = out.toarray()
                return as2d( out )
        else:
            def call( vec ):
                return as2d( nm.asarray( sc.dot( operator.obj, vec ) ) )
        operator.call = call

    elif isinstance( operatorInput, types.FunctionType ) or \
         isinstance( operatorInput, types.BuiltinFunctionType ):
        operator.shape = expectedShape
        operator.dtype = nm.float64
        operator.call = operatorInput

    return operator

##
# 05.04.2007, c
def applyConstraints( blockVectorV, factYBY, blockVectorBY, blockVectorY ):
    """Changes blockVectorV in place."""
    gramYBV = sc.dot( blockVectorBY.T, blockVectorV )
    tmp = la.cho_solve( factYBY, gramYBV )
    blockVectorV -= sc.dot( blockVectorY, tmp )

##
# 05.04.2007, c
def b_orthonormalize( operatorB, blockVectorV,
                      blockVectorBV = None, retInvR = False ):

    if blockVectorBV is None:
        if operatorB is not None:
            blockVectorBV = operatorB( blockVectorV )
        else:
            blockVectorBV = blockVectorV # Shared data!!!
    gramVBV = sc.dot( blockVectorV.T, blockVectorBV )
    gramVBV = la.cholesky( gramVBV )
    la.inv( gramVBV, overwrite_a = True )
    # gramVBV is now R^{-1}.
    blockVectorV = sc.dot( blockVectorV, gramVBV )
    if operatorB is not None:
        blockVectorBV = sc.dot( blockVectorBV, gramVBV )

    if retInvR:
        return blockVectorV, blockVectorBV, gramVBV
    else:
        return blockVectorV, blockVectorBV

##
# 04.04.2007, c
# 05.04.2007
# 06.04.2007
# 10.04.2007
def lobpcg( blockVectorX, operatorA,
            operatorB = None, operatorT = None, blockVectorY = None,
            residualTolerance = None, maxIterations = 20,
            largest = True, verbosityLevel = 0,
            retLambdaHistory = False, retResidualNormsHistory = False ):

    exitFlag = 0

    if blockVectorY is not None:
        sizeY = blockVectorY.shape[1]
    else:
        sizeY = 0

    # Block size.
    n, sizeX = blockVectorX.shape
    if sizeX > n:
        raise ValueError,\
              'the first input argument blockVectorX must be tall, not fat' +\
              ' (%d, %d)' % blockVectorX.shape

    if n < 1:
        raise ValueError,\
              'the matrix size is wrong (%d)' % n
        
    operatorA = makeOperator( operatorA, (n, n) )

    if operatorB is not None:
        operatorB = makeOperator( operatorB, (n, n) )

    if operatorT is not None:
        operatorT = makeOperator( operatorT, (n, n) )
##     if n != operatorA.shape[0]:
##         aux = 'The size (%d, %d) of operatorA is not the same as\n'+\
##               '%d - the number of rows of blockVectorX' % operatorA.shape + (n,)
##         raise ValueError, aux

##     if operatorA.shape[0] != operatorA.shape[1]:
##         raise ValueError, 'operatorA must be a square matrix (%d, %d)' %\
##               operatorA.shape

    if residualTolerance is None:
        residualTolerance = nm.sqrt( 1e-15 ) * n

    maxIterations = min( n, maxIterations )

    if verbosityLevel:
        aux = "Solving "
        if operatorB is None:
            aux += "standard"
        else:
            aux += "generalized"
        aux += " eigenvalue problem with"
        if operatorT is None:
            aux += "out"
        aux += " preconditioning\n\n"
        aux += "matrix size %d\n" % n
        aux += "block size %d\n\n" % sizeX
        if blockVectorY is None:
            aux += "No constraints\n\n"
        else:
            if sizeY > 1:
                aux += "%d constraints\n\n" % sizeY
            else:
                aux += "%d constraint\n\n" % sizeY
        print aux

    ##
    # Apply constraints to X.
    if blockVectorY is not None:

        if operatorB is not None:
            blockVectorBY = operatorB( blockVectorY )
        else:
            blockVectorBY = blockVectorY
    
        # gramYBY is a dense array.
        gramYBY = sc.dot( blockVectorY.T, blockVectorBY )
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = la.cho_factor( gramYBY )
        except:
            print 'cannot handle linear dependent constraints'
            raise

        applyConstraints( blockVectorX, gramYBY, blockVectorBY, blockVectorY )

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX = b_orthonormalize( operatorB, blockVectorX )

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = operatorA( blockVectorX )
    gramXAX = sc.dot( blockVectorX.T, blockVectorAX )
    # gramXBX is X^T * X.
    gramXBX = sc.dot( blockVectorX.T, blockVectorX )
    _lambda, eigBlockVector = symeig( gramXAX )
    ii = nm.argsort( _lambda )[:sizeX]
    if largest:
        ii = ii[::-1]
    _lambda = _lambda[ii]
    eigBlockVector = nm.asarray( eigBlockVector[:,ii] )
#    pause()
    blockVectorX = sc.dot( blockVectorX, eigBlockVector )
    blockVectorAX = sc.dot( blockVectorAX, eigBlockVector )
    if operatorB is not None:
        blockVectorBX = sc.dot( blockVectorBX, eigBlockVector )
    
    ##
    # Active index set.
    activeMask = nm.ones( (sizeX,), dtype = nm.bool )

    lambdaHistory = [_lambda]
    residualNormsHistory = []

    previousBlockSize = sizeX
    ident = nm.eye( sizeX, dtype = operatorA.dtype )
    ident0 = nm.eye( sizeX, dtype = operatorA.dtype )
    
    ##
    # Main iteration loop.
    for iterationNumber in xrange( maxIterations ):
        if verbosityLevel > 0:
            print 'iteration %d' %  iterationNumber

        aux = blockVectorBX * _lambda[nm.newaxis,:]
        blockVectorR = blockVectorAX - aux

        aux = nm.sum( blockVectorR.conjugate() * blockVectorR, 0 )
        residualNorms = nm.sqrt( aux )

        
##         if iterationNumber == 2:
##             print blockVectorAX
##             print blockVectorBX
##             print blockVectorR
##             pause()

        residualNormsHistory.append( residualNorms )

        ii = nm.where( residualNorms > residualTolerance, True, False )
        activeMask = activeMask & ii
        if verbosityLevel > 2:
            print activeMask

        currentBlockSize = activeMask.sum()
        if currentBlockSize != previousBlockSize:
            previousBlockSize = currentBlockSize
            ident = nm.eye( currentBlockSize, dtype = operatorA.dtype )


        if currentBlockSize == 0:
            failureFlag = 0 # All eigenpairs converged.
            break

        if verbosityLevel > 0:
            print 'current block size:', currentBlockSize
            print 'eigenvalue:', _lambda
            print 'residual norms:', residualNorms
        if verbosityLevel > 10:
            print eigBlockVector

        activeBlockVectorR = as2d( blockVectorR[:,activeMask] )
        
        if iterationNumber > 0:
            activeBlockVectorP = as2d( blockVectorP[:,activeMask] )
            activeBlockVectorAP = as2d( blockVectorAP[:,activeMask] )
            activeBlockVectorBP = as2d( blockVectorBP[:,activeMask] )

#        print activeBlockVectorR
        if operatorT is not None:
            ##
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = operatorT( activeBlockVectorR )

#        assert nm.all( blockVectorR == activeBlockVectorR )

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            applyConstraints( activeBlockVectorR,
                              gramYBY, blockVectorBY, blockVectorY )

#        assert nm.all( blockVectorR == activeBlockVectorR )

        ##
        # B-orthonormalize the preconditioned residuals.
#        print activeBlockVectorR

        aux = b_orthonormalize( operatorB, activeBlockVectorR )
        activeBlockVectorR, activeBlockVectorBR = aux
#        print activeBlockVectorR

        activeBlockVectorAR = operatorA( activeBlockVectorR )

        if iterationNumber > 0:
            aux = b_orthonormalize( operatorB, activeBlockVectorP,
                                    activeBlockVectorBP, retInvR = True )
            activeBlockVectorP, activeBlockVectorBP, invR = aux
            activeBlockVectorAP = sc.dot( activeBlockVectorAP, invR )

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        xaw = sc.dot( blockVectorX.T, activeBlockVectorAR )
        waw = sc.dot( activeBlockVectorR.T, activeBlockVectorAR )
        xbw = sc.dot( blockVectorX.T, activeBlockVectorBR )
        
        if iterationNumber > 0:
            xap = sc.dot( blockVectorX.T, activeBlockVectorAP )
            wap = sc.dot( activeBlockVectorR.T, activeBlockVectorAP )
            pap = sc.dot( activeBlockVectorP.T, activeBlockVectorAP )
            xbp = sc.dot( blockVectorX.T, activeBlockVectorBP )
            wbp = sc.dot( activeBlockVectorR.T, activeBlockVectorBP )
            
            gramA = nm.bmat( [[nm.diag( _lambda ), xaw, xap],
                              [xaw.T, waw, wap],
                              [xap.T, wap.T, pap]] )
            try:
                gramB = nm.bmat( [[ident0, xbw, xbp],
                                  [xbw.T, ident, wbp],
                                  [xbp.T, wbp.T, ident]] )
            except:
                print ident
                print xbw
                raise
        else:
            gramA = nm.bmat( [[nm.diag( _lambda ), xaw],
                              [xaw.T, waw]] )
            gramB = nm.bmat( [[ident0, xbw],
                              [xbw.T, ident0]] )
        try:
            assert nm.allclose( gramA.T, gramA )
        except:
            print gramA.T - gramA
            raise

        try:
            assert nm.allclose( gramB.T, gramB )
        except:
            print gramB.T - gramB
            raise

##         print nm.diag( _lambda )
##         print xaw
##         print waw
##         print xbw
##         try:
##             print xap
##             print wap
##             print pap
##             print xbp
##             print wbp
##         except:
##             pass
##         pause()

        if verbosityLevel > 10:
            save( gramA, 'gramA' )
            save( gramB, 'gramB' )
        ##
        # Solve the generalized eigenvalue problem.
#        _lambda, eigBlockVector = la.eig( gramA, gramB )
        _lambda, eigBlockVector = symeig( gramA, gramB )
        ii = nm.argsort( _lambda )[:sizeX]
        if largest:
            ii = ii[::-1]
        if verbosityLevel > 10:
            print ii
        
        _lambda = _lambda[ii].astype( nm.float64 )
        eigBlockVector = nm.asarray( eigBlockVector[:,ii].astype( nm.float64 ) )
        if verbosityLevel > 10:
            print 'lambda:', _lambda
##         # Normalize eigenvectors!
##         aux = nm.sum( eigBlockVector.conjugate() * eigBlockVector, 0 )
##         eigVecNorms = nm.sqrt( aux )
##         eigBlockVector = eigBlockVector / eigVecNorms[nm.newaxis,:]
#        eigBlockVector, aux = b_orthonormalize( operatorB, eigBlockVector )

        if verbosityLevel > 10:
            print eigBlockVector
            pause()
        ##
        # Compute Ritz vectors.
        if iterationNumber > 0:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:sizeX+currentBlockSize]
            eigBlockVectorP = eigBlockVector[sizeX+currentBlockSize:]

            pp = sc.dot( activeBlockVectorR, eigBlockVectorR )\
                 + sc.dot( activeBlockVectorP, eigBlockVectorP )

            app = sc.dot( activeBlockVectorAR, eigBlockVectorR )\
                  + sc.dot( activeBlockVectorAP, eigBlockVectorP )

            bpp = sc.dot( activeBlockVectorBR, eigBlockVectorR )\
                  + sc.dot( activeBlockVectorBP, eigBlockVectorP )
        else:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:]

            pp = sc.dot( activeBlockVectorR, eigBlockVectorR )

            app = sc.dot( activeBlockVectorAR, eigBlockVectorR )

            bpp = sc.dot( activeBlockVectorBR, eigBlockVectorR )

        if verbosityLevel > 10:
            print pp
            print app
            print bpp
            pause()
#        print pp.shape, app.shape, bpp.shape

        blockVectorX = sc.dot( blockVectorX, eigBlockVectorX ) + pp
        blockVectorAX = sc.dot( blockVectorAX, eigBlockVectorX ) + app
        blockVectorBX = sc.dot( blockVectorBX, eigBlockVectorX ) + bpp

        blockVectorP, blockVectorAP, blockVectorBP = pp, app, bpp

    aux = blockVectorBX * _lambda[nm.newaxis,:]
    blockVectorR = blockVectorAX - aux

    aux = nm.sum( blockVectorR.conjugate() * blockVectorR, 0 )
    residualNorms = nm.sqrt( aux )


    if verbosityLevel > 0:
        print 'final eigenvalue:', _lambda
        print 'final residual norms:', residualNorms


    return _lambda, eigBlockVectorX

###########################################################################
if __name__ == '__main__':
    from scipy.sparse import spdiags, speye
    import time

##     def operatorB( vec ):
##         return vec

    n = 100
    vals = [nm.arange( n, dtype = nm.float64 ) + 1]
    operatorA = spdiags( vals, 0, n, n )
    operatorB = speye( n, n )
#    operatorB[0,0] = 0
    operatorB = nm.eye( n, n )
    Y = nm.eye( n, 3 )


##    X = sc.rand( n, 3 )
    xfile = {100 : 'X.txt', 1000 : 'X2.txt', 10000 : 'X3.txt'}
    X = nm.fromfile( xfile[n], dtype = nm.float64, sep = ' ' )
    X.shape = (n, 3)

    ivals = [1./vals[0]]
    def precond( x ):
        invA = spdiags( ivals, 0, n, n )
        y = invA  * x
        if sp.issparse( y ):
            y = y.toarray()

        return as2d( y )

    tt = time.clock()
    eigs, vecs = lobpcg( X, operatorA, operatorB, blockVectorY = Y,
                         operatorT = precond,
                         residualTolerance = 1e-4, maxIterations = 40,
                         largest = False, verbosityLevel = 1 )
    print 'solution time:', time.clock() - tt
    print eigs
    
