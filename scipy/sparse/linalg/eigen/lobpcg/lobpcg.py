"""
Pure SciPy implementation of Locally Optimal Block Preconditioned Conjugate
Gradient Method (LOBPCG), see
http://www-math.cudenver.edu/~aknyazev/software/BLOPEX/

License: BSD

Depends upon symeig (http://mdp-toolkit.sourceforge.net/symeig.html) for the
moment, as the symmetric eigenvalue solvers were not available in scipy.

(c) Robert Cimrman, Andrew Knyazev

Examples in tests directory contributed by Nils Wagner.
"""

import types
from warnings import warn

import numpy as nm
import scipy as sc
import scipy.sparse as sp
import scipy.linalg as la
import scipy.io as io
from scipy.sparse.linalg import aslinearoperator, LinearOperator

try:
    from symeig import symeig
except:
    raise ImportError('lobpcg requires symeig')

def pause():
    raw_input()

def save( ar, fileName ):
    io.write_array( fileName, ar, precision = 8 )

##
# 21.05.2007, c
def as2d( ar ):
    """
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else: # Assume 1!
        aux = nm.array( ar, copy = False )
        aux.shape = (ar.shape[0], 1)
        return aux

def makeOperator( operatorInput, expectedShape ):
    """Internal. Takes a dense numpy array or a sparse matrix or 
    a function and makes an operator performing matrix * blockvector 
    products.

    Example
    -------

    A = makeOperator( arrayA, (n, n) )
    vectorB = A( vectorX )
    """
    if operatorInput is None:
        def ident(x):
            return x
        operator = LinearOperator(expectedShape, ident, matmat=ident)
    else:
        operator = aslinearoperator(operatorInput)

    if operator.shape != expectedShape:
        raise ValueError('operator has invalid shape')
    
    operator.__call__ = operator.matmat

    return operator



def applyConstraints( blockVectorV, factYBY, blockVectorBY, blockVectorY ):
    """Internal. Changes blockVectorV in place."""
    gramYBV = sc.dot( blockVectorBY.T, blockVectorV )
    tmp = la.cho_solve( factYBY, gramYBV )
    blockVectorV -= sc.dot( blockVectorY, tmp )


def b_orthonormalize( B, blockVectorV,
                      blockVectorBV = None, retInvR = False ):
    """Internal."""
    if blockVectorBV is None:
        if B is not None:
            blockVectorBV = B( blockVectorV )
        else:
            blockVectorBV = blockVectorV # Shared data!!!
    gramVBV = sc.dot( blockVectorV.T, blockVectorBV )
    gramVBV = la.cholesky( gramVBV )
    la.inv( gramVBV, overwrite_a = True )
    # gramVBV is now R^{-1}.
    blockVectorV = sc.dot( blockVectorV, gramVBV )
    if B is not None:
        blockVectorBV = sc.dot( blockVectorBV, gramVBV )

    if retInvR:
        return blockVectorV, blockVectorBV, gramVBV
    else:
        return blockVectorV, blockVectorBV

def lobpcg( blockVectorX, A,
            B = None, M = None, blockVectorY = None,
            residualTolerance = None, maxIterations = 20,
            largest = True, verbosityLevel = 0,
            retLambdaHistory = False, retResidualNormsHistory = False ):
    """Solve symmetric partial eigenproblems with optional preconditioning

    This function implements the Locally Optimal Block Preconditioned 
    Conjugate Gradient Method (LOBPCG).

    TODO write in terms of Ax=lambda B x

    Parameters
    ----------
    blockVectorX : array_like 
        initial approximation to eigenvectors shape=(n,blockSize)
    A : {dense matrix, sparse matrix, LinearOperator}
        the linear operator of the problem, usually a sparse matrix
        often called the "stiffness matrix"

    Returns
    -------
    (lambda,blockVectorV) : tuple of arrays
        blockVectorX and lambda are computed blockSize eigenpairs, where
        blockSize=size(blockVectorX,2) for the initial guess blockVectorX 
        if it is full rank.

    Optional Parameters
    -------------------
    B : {dense matrix, sparse matrix, LinearOperator}
        the right hand side operator in a generalized eigenproblem.
        by default, B = Identity
        often called the "mass matrix"
    M : {dense matrix, sparse matrix, LinearOperator}
        preconditioner to A; by default M = Identity
        M should approximate the inverse of A
    blockVectorY : array_like
        n-by-sizeY matrix of constraints, sizeY < n
        The iterations will be performed in the B-orthogonal complement 
        of the column-space of blockVectorY. blockVectorY must be full rank.

    Other Parameters
    ----------------
    residualTolerance : scalar
        solver tolerance. default: residualTolerance=n*sqrt(eps)
    maxIterations: integer
        maximum number of iterations
        by default: maxIterations=min(n,20)
    largest : boolean
        when True, solve for the largest eigenvalues, otherwise the smallest
    verbosityLevel : integer
        controls solver output.  default: verbosityLevel = 0.
    retLambdaHistory : boolean
        whether to return eigenvalue history
    retResidualNormsHistory : boolean
        whether to return history of residual norms


    Notes
    -----
    If both retLambdaHistory and retResidualNormsHistory are True, the
    return tuple has the following format:
        (lambda, blockVectorV, lambda history, residual norms history)

    """
    failureFlag = True

    if blockVectorY is not None:
        sizeY = blockVectorY.shape[1]
    else:
        sizeY = 0

    # Block size.
    if len(blockVectorX.shape) != 2:
        raise ValueError('expected rank-2 array for argument blockVectorX')

    n, sizeX = blockVectorX.shape
    if sizeX > n:
        raise ValueError('blockVectorX column dimension exceeds the row dimension')

    A = makeOperator(A, (n,n))
    B = makeOperator(B, (n,n))
    M = makeOperator(M, (n,n))

    if (n - sizeY) < (5 * sizeX):
        #warn('The problem size is small compared to the block size.' \
        #        ' Using dense eigensolver instead of LOBPCG.')

        if blockVectorY is not None:
            raise NotImplementedError('symeig does not support constraints')

        if largest:
            lohi = (n - sizeX, n)
        else:
            lohi = (1, sizeX)

        A_dense = A(nm.eye(n))

        if B is not None:
            B_dense = B(nm.eye(n))
            _lambda, eigBlockVector = symeig(A_dense, B_dense, range=lohi )
        else:
            _lambda, eigBlockVector = symeig(A_dense, range=lohi )

        return _lambda, eigBlockVector


    if residualTolerance is None:
        residualTolerance = nm.sqrt( 1e-15 ) * n

    maxIterations = min( n, maxIterations )

    if verbosityLevel:
        aux = "Solving "
        if B is None:
            aux += "standard"
        else:
            aux += "generalized"
        aux += " eigenvalue problem with"
        if M is None:
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

        if B is not None:
            blockVectorBY = B( blockVectorY )
        else:
            blockVectorBY = blockVectorY

        # gramYBY is a dense array.
        gramYBY = sc.dot( blockVectorY.T, blockVectorBY )
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = la.cho_factor( gramYBY )
        except:
            raise ValueError('cannot handle linear dependent constraints')

        applyConstraints( blockVectorX, gramYBY, blockVectorBY, blockVectorY )

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX = b_orthonormalize( B, blockVectorX )

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A( blockVectorX )
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
    if B is not None:
        blockVectorBX = sc.dot( blockVectorBX, eigBlockVector )

    ##
    # Active index set.
    activeMask = nm.ones( (sizeX,), dtype = nm.bool )

    lambdaHistory = [_lambda]
    residualNormsHistory = []

    previousBlockSize = sizeX
    ident = nm.eye( sizeX, dtype = A.dtype )
    ident0 = nm.eye( sizeX, dtype = A.dtype )

    ##
    # Main iteration loop.
    for iterationNumber in xrange( maxIterations ):
        if verbosityLevel > 0:
            print 'iteration %d' %  iterationNumber

        aux = blockVectorBX * _lambda[nm.newaxis,:]
        blockVectorR = blockVectorAX - aux

        aux = nm.sum( blockVectorR.conjugate() * blockVectorR, 0 )
        residualNorms = nm.sqrt( aux )

        residualNormsHistory.append( residualNorms )

        ii = nm.where( residualNorms > residualTolerance, True, False )
        activeMask = activeMask & ii
        if verbosityLevel > 2:
            print activeMask

        currentBlockSize = activeMask.sum()
        if currentBlockSize != previousBlockSize:
            previousBlockSize = currentBlockSize
            ident = nm.eye( currentBlockSize, dtype = A.dtype )

        if currentBlockSize == 0:
            failureFlag = False # All eigenpairs converged.
            break

        if verbosityLevel > 0:
            print 'current block size:', currentBlockSize
            print 'eigenvalue:', _lambda
            print 'residual norms:', residualNorms
        if verbosityLevel > 10:
            print eigBlockVector

        activeBlockVectorR = as2d( blockVectorR[:,activeMask] )

        if iterationNumber > 0:
            activeBlockVectorP  = as2d( blockVectorP [:,activeMask] )
            activeBlockVectorAP = as2d( blockVectorAP[:,activeMask] )
            activeBlockVectorBP = as2d( blockVectorBP[:,activeMask] )

        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M( activeBlockVectorR )

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            applyConstraints( activeBlockVectorR,
                              gramYBY, blockVectorBY, blockVectorY )

        ##
        # B-orthonormalize the preconditioned residuals.

        aux = b_orthonormalize( B, activeBlockVectorR )
        activeBlockVectorR, activeBlockVectorBR = aux

        activeBlockVectorAR = A( activeBlockVectorR )

        if iterationNumber > 0:
            aux = b_orthonormalize( B, activeBlockVectorP,
                                    activeBlockVectorBP, retInvR = True )
            activeBlockVectorP, activeBlockVectorBP, invR = aux
            activeBlockVectorAP = sc.dot( activeBlockVectorAP, invR )

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        xaw = sc.dot( blockVectorX.T,       activeBlockVectorAR )
        waw = sc.dot( activeBlockVectorR.T, activeBlockVectorAR )
        xbw = sc.dot( blockVectorX.T,       activeBlockVectorBR )

        if iterationNumber > 0:
            xap = sc.dot( blockVectorX.T,       activeBlockVectorAP )
            wap = sc.dot( activeBlockVectorR.T, activeBlockVectorAP )
            pap = sc.dot( activeBlockVectorP.T, activeBlockVectorAP )
            xbp = sc.dot( blockVectorX.T,       activeBlockVectorBP )
            wbp = sc.dot( activeBlockVectorR.T, activeBlockVectorBP )

            gramA = nm.bmat( [[nm.diag( _lambda ), xaw, xap],
                              [xaw.T, waw, wap],
                              [xap.T, wap.T, pap]] )
            try:
                gramB = nm.bmat( [[ident0,   xbw,   xbp],
                                  [ xbw.T, ident,   wbp],
                                  [ xbp.T, wbp.T, ident]] )
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

        lambdaHistory.append( _lambda )

        if verbosityLevel > 10:
            print 'lambda:', _lambda
##         # Normalize eigenvectors!
##         aux = nm.sum( eigBlockVector.conjugate() * eigBlockVector, 0 )
##         eigVecNorms = nm.sqrt( aux )
##         eigBlockVector = eigBlockVector / eigVecNorms[nm.newaxis,:]
#        eigBlockVector, aux = b_orthonormalize( B, eigBlockVector )

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

    if retLambdaHistory:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, lambdaHistory, residualNormsHistory
        else:
            return _lambda, blockVectorX, lambdaHistory
    else:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, residualNormsHistory
        else:
            return _lambda, blockVectorX

###########################################################################
if __name__ == '__main__':
    from scipy.sparse import spdiags, speye
    import time

##     def B( vec ):
##         return vec

    n = 100
    vals = [nm.arange( n, dtype = nm.float64 ) + 1]
    A = spdiags( vals, 0, n, n )
    B = speye( n, n )
#    B[0,0] = 0
    B = nm.eye( n, n )
    Y = nm.eye( n, 3 )


#    X = sc.rand( n, 3 )
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

#    precond = spdiags( ivals, 0, n, n )

    tt = time.clock()
    eigs, vecs = lobpcg( X, A, B, blockVectorY = Y,
                         M = precond,
                         residualTolerance = 1e-4, maxIterations = 40,
                         largest = False, verbosityLevel = 1 )
    print 'solution time:', time.clock() - tt
    print eigs

    print vecs
