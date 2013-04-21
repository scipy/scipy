"""
Pure SciPy implementation of Locally Optimal Block Preconditioned Conjugate
Gradient Method (LOBPCG), see
http://www-math.cudenver.edu/~aknyazev/software/BLOPEX/

License: BSD

Authors: Robert Cimrman, Andrew Knyazev

Examples in tests directory contributed by Nils Wagner.
"""

from __future__ import division, print_function, absolute_import

import sys
import numpy as np
import scipy as sp

from scipy.lib.six.moves import xrange

from scipy.sparse.linalg import aslinearoperator, LinearOperator

__all__ = ['lobpcg']

## try:
##     from symeig import symeig
## except:
##     raise ImportError('lobpcg requires symeig')


def symeig(mtxA, mtxB=None, eigenvectors=True, select=None):
    import scipy.linalg as sla
    if select is None:
        if np.iscomplexobj(mtxA):
            if mtxB is None:
                fun = sla.get_lapack_funcs('heev', arrays=(mtxA,))
            else:
                fun = sla.get_lapack_funcs('hegv', arrays=(mtxA,))
        else:
            if mtxB is None:
                fun = sla.get_lapack_funcs('syev', arrays=(mtxA,))
            else:
                fun = sla.get_lapack_funcs('sygv', arrays=(mtxA,))
##         print fun
        if mtxB is None:
            out = fun(mtxA)
        else:
            out = fun(mtxA, mtxB)
            out = out[1], out[0], out[2]
##         print w
##         print v
##         print info
##         from symeig import symeig
##         print symeig( mtxA, mtxB )
    else:
        out = sla.eig(mtxA, mtxB, right=eigenvectors)
        w = out[0]
        ii = np.argsort(w)
        w = w[slice(*select)]
        if eigenvectors:
            v = out[1][:,ii]
            v = v[:,slice(*select)]
            out = w, v, 0
        else:
            out = w, 0

    return out[:-1]


def pause():
    input()


def save(ar, fileName):
    from numpy import savetxt
    savetxt(fileName, ar, precision=8)

##
# 21.05.2007, c


def as2d(ar):
    """
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


class CallableLinearOperator(LinearOperator):
    def __call__(self, x):
        return self.matmat(x)


def makeOperator(operatorInput, expectedShape):
    """Internal. Takes a dense numpy array or a sparse matrix or
    a function and makes an operator performing matrix * blockvector
    products.

    Examples
    --------
    >>> A = makeOperator( arrayA, (n, n) )
    >>> vectorB = A( vectorX )

    """
    if operatorInput is None:
        def ident(x):
            return x
        operator = LinearOperator(expectedShape, ident, matmat=ident)
    else:
        operator = aslinearoperator(operatorInput)

    if operator.shape != expectedShape:
        raise ValueError('operator has invalid shape')

    if sys.version_info[0] >= 3:
        # special methods are looked up on the class -- so make a new one
        operator.__class__ = CallableLinearOperator
    else:
        operator.__call__ = operator.matmat

    return operator


def applyConstraints(blockVectorV, factYBY, blockVectorBY, blockVectorY):
    """Internal. Changes blockVectorV in place."""
    gramYBV = sp.dot(blockVectorBY.T, blockVectorV)
    import scipy.linalg as sla
    tmp = sla.cho_solve(factYBY, gramYBV)
    blockVectorV -= sp.dot(blockVectorY, tmp)


def b_orthonormalize(B, blockVectorV,
                      blockVectorBV=None, retInvR=False):
    """Internal."""
    import scipy.linalg as sla
    if blockVectorBV is None:
        if B is not None:
            blockVectorBV = B(blockVectorV)
        else:
            blockVectorBV = blockVectorV  # Shared data!!!
    gramVBV = sp.dot(blockVectorV.T, blockVectorBV)
    gramVBV = sla.cholesky(gramVBV)
    gramVBV = sla.inv(gramVBV, overwrite_a=True)
    # gramVBV is now R^{-1}.
    blockVectorV = sp.dot(blockVectorV, gramVBV)
    if B is not None:
        blockVectorBV = sp.dot(blockVectorBV, gramVBV)

    if retInvR:
        return blockVectorV, blockVectorBV, gramVBV
    else:
        return blockVectorV, blockVectorBV


def lobpcg(A, X,
            B=None, M=None, Y=None,
            tol=None, maxiter=20,
            largest=True, verbosityLevel=0,
            retLambdaHistory=False, retResidualNormsHistory=False):
    """Solve symmetric partial eigenproblems with optional preconditioning

    This function implements the Locally Optimal Block Preconditioned
    Conjugate Gradient Method (LOBPCG).

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The symmetric linear operator of the problem, usually a
        sparse matrix.  Often called the "stiffness matrix".
    X : array_like
        Initial approximation to the k eigenvectors. If A has
        shape=(n,n) then X should have shape shape=(n,k).
    B : {dense matrix, sparse matrix, LinearOperator}, optional
        the right hand side operator in a generalized eigenproblem.
        by default, B = Identity
        often called the "mass matrix"
    M : {dense matrix, sparse matrix, LinearOperator}, optional
        preconditioner to A; by default M = Identity
        M should approximate the inverse of A
    Y : array_like, optional
        n-by-sizeY matrix of constraints, sizeY < n
        The iterations will be performed in the B-orthogonal complement
        of the column-space of Y. Y must be full rank.

    Returns
    -------
    w : array
        Array of k eigenvalues
    v : array
        An array of k eigenvectors.  V has the same shape as X.

    Other Parameters
    ----------------
    tol : scalar, optional
        Solver tolerance (stopping criterion)
        by default: tol=n*sqrt(eps)
    maxiter : integer, optional
        maximum number of iterations
        by default: maxiter=min(n,20)
    largest : boolean, optional
        when True, solve for the largest eigenvalues, otherwise the smallest
    verbosityLevel : integer, optional
        controls solver output.  default: verbosityLevel = 0.
    retLambdaHistory : boolean, optional
        whether to return eigenvalue history
    retResidualNormsHistory : boolean, optional
        whether to return history of residual norms


    Notes
    -----
    If both retLambdaHistory and retResidualNormsHistory are True, the
    return tuple has the following format
    (lambda, V, lambda history, residual norms history)

    """
    failureFlag = True
    import scipy.linalg as sla

    blockVectorX = X
    blockVectorY = Y
    residualTolerance = tol
    maxIterations = maxiter

    if blockVectorY is not None:
        sizeY = blockVectorY.shape[1]
    else:
        sizeY = 0

    # Block size.
    if len(blockVectorX.shape) != 2:
        raise ValueError('expected rank-2 array for argument X')

    n, sizeX = blockVectorX.shape
    if sizeX > n:
        raise ValueError('X column dimension exceeds the row dimension')

    A = makeOperator(A, (n,n))
    B = makeOperator(B, (n,n))
    M = makeOperator(M, (n,n))

    if (n - sizeY) < (5 * sizeX):
        # warn('The problem size is small compared to the block size.' \
        #        ' Using dense eigensolver instead of LOBPCG.')

        if blockVectorY is not None:
            raise NotImplementedError('symeig does not support constraints')

        if largest:
            lohi = (n - sizeX, n)
        else:
            lohi = (1, sizeX)

        A_dense = A(np.eye(n))

        if B is not None:
            B_dense = B(np.eye(n))
            _lambda, eigBlockVector = symeig(A_dense, B_dense, select=lohi)
        else:
            _lambda, eigBlockVector = symeig(A_dense, select=lohi)

        return _lambda, eigBlockVector

    if residualTolerance is None:
        residualTolerance = np.sqrt(1e-15) * n

    maxIterations = min(n, maxIterations)

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
        print(aux)

    ##
    # Apply constraints to X.
    if blockVectorY is not None:

        if B is not None:
            blockVectorBY = B(blockVectorY)
        else:
            blockVectorBY = blockVectorY

        # gramYBY is a dense array.
        gramYBY = sp.dot(blockVectorY.T, blockVectorBY)
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = sla.cho_factor(gramYBY)
        except:
            raise ValueError('cannot handle linearly dependent constraints')

        applyConstraints(blockVectorX, gramYBY, blockVectorBY, blockVectorY)

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX = b_orthonormalize(B, blockVectorX)

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A(blockVectorX)
    gramXAX = sp.dot(blockVectorX.T, blockVectorAX)
    # gramXBX is X^T * X.
    gramXBX = sp.dot(blockVectorX.T, blockVectorX)

    _lambda, eigBlockVector = symeig(gramXAX)
    ii = np.argsort(_lambda)[:sizeX]
    if largest:
        ii = ii[::-1]
    _lambda = _lambda[ii]

    eigBlockVector = np.asarray(eigBlockVector[:,ii])
    blockVectorX = sp.dot(blockVectorX, eigBlockVector)
    blockVectorAX = sp.dot(blockVectorAX, eigBlockVector)
    if B is not None:
        blockVectorBX = sp.dot(blockVectorBX, eigBlockVector)

    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=np.bool)

    lambdaHistory = [_lambda]
    residualNormsHistory = []

    previousBlockSize = sizeX
    ident = np.eye(sizeX, dtype=A.dtype)
    ident0 = np.eye(sizeX, dtype=A.dtype)

    ##
    # Main iteration loop.
    for iterationNumber in xrange(maxIterations):
        if verbosityLevel > 0:
            print('iteration %d' % iterationNumber)

        aux = blockVectorBX * _lambda[np.newaxis,:]
        blockVectorR = blockVectorAX - aux

        aux = np.sum(blockVectorR.conjugate() * blockVectorR, 0)
        residualNorms = np.sqrt(aux)

        residualNormsHistory.append(residualNorms)

        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        if verbosityLevel > 2:
            print(activeMask)

        currentBlockSize = activeMask.sum()
        if currentBlockSize != previousBlockSize:
            previousBlockSize = currentBlockSize
            ident = np.eye(currentBlockSize, dtype=A.dtype)

        if currentBlockSize == 0:
            failureFlag = False  # All eigenpairs converged.
            break

        if verbosityLevel > 0:
            print('current block size:', currentBlockSize)
            print('eigenvalue:', _lambda)
            print('residual norms:', residualNorms)
        if verbosityLevel > 10:
            print(eigBlockVector)

        activeBlockVectorR = as2d(blockVectorR[:,activeMask])

        if iterationNumber > 0:
            activeBlockVectorP = as2d(blockVectorP[:,activeMask])
            activeBlockVectorAP = as2d(blockVectorAP[:,activeMask])
            activeBlockVectorBP = as2d(blockVectorBP[:,activeMask])

        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR)

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            applyConstraints(activeBlockVectorR,
                              gramYBY, blockVectorBY, blockVectorY)

        ##
        # B-orthonormalize the preconditioned residuals.

        aux = b_orthonormalize(B, activeBlockVectorR)
        activeBlockVectorR, activeBlockVectorBR = aux

        activeBlockVectorAR = A(activeBlockVectorR)

        if iterationNumber > 0:
            aux = b_orthonormalize(B, activeBlockVectorP,
                                    activeBlockVectorBP, retInvR=True)
            activeBlockVectorP, activeBlockVectorBP, invR = aux
            activeBlockVectorAP = sp.dot(activeBlockVectorAP, invR)

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        xaw = sp.dot(blockVectorX.T, activeBlockVectorAR)
        waw = sp.dot(activeBlockVectorR.T, activeBlockVectorAR)
        xbw = sp.dot(blockVectorX.T, activeBlockVectorBR)

        if iterationNumber > 0:
            xap = sp.dot(blockVectorX.T, activeBlockVectorAP)
            wap = sp.dot(activeBlockVectorR.T, activeBlockVectorAP)
            pap = sp.dot(activeBlockVectorP.T, activeBlockVectorAP)
            xbp = sp.dot(blockVectorX.T, activeBlockVectorBP)
            wbp = sp.dot(activeBlockVectorR.T, activeBlockVectorBP)

            gramA = np.bmat([[np.diag(_lambda), xaw, xap],
                              [xaw.T, waw, wap],
                              [xap.T, wap.T, pap]])

            gramB = np.bmat([[ident0, xbw, xbp],
                              [xbw.T, ident, wbp],
                              [xbp.T, wbp.T, ident]])
        else:
            gramA = np.bmat([[np.diag(_lambda), xaw],
                              [xaw.T, waw]])
            gramB = np.bmat([[ident0, xbw],
                              [xbw.T, ident]])

        try:
            assert np.allclose(gramA.T, gramA)
        except:
            print(gramA.T - gramA)
            raise

        try:
            assert np.allclose(gramB.T, gramB)
        except:
            print(gramB.T - gramB)
            raise

        if verbosityLevel > 10:
            save(gramA, 'gramA')
            save(gramB, 'gramB')

        ##
        # Solve the generalized eigenvalue problem.
#        _lambda, eigBlockVector = la.eig( gramA, gramB )
        _lambda, eigBlockVector = symeig(gramA, gramB)
        ii = np.argsort(_lambda)[:sizeX]
        if largest:
            ii = ii[::-1]
        if verbosityLevel > 10:
            print(ii)

        _lambda = _lambda[ii].astype(np.float64)
        eigBlockVector = np.asarray(eigBlockVector[:,ii].astype(np.float64))

        lambdaHistory.append(_lambda)

        if verbosityLevel > 10:
            print('lambda:', _lambda)
##         # Normalize eigenvectors!
##         aux = np.sum( eigBlockVector.conjugate() * eigBlockVector, 0 )
##         eigVecNorms = np.sqrt( aux )
##         eigBlockVector = eigBlockVector / eigVecNorms[np.newaxis,:]
#        eigBlockVector, aux = b_orthonormalize( B, eigBlockVector )

        if verbosityLevel > 10:
            print(eigBlockVector)
            pause()

        ##
        # Compute Ritz vectors.
        if iterationNumber > 0:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:sizeX+currentBlockSize]
            eigBlockVectorP = eigBlockVector[sizeX+currentBlockSize:]

            pp = sp.dot(activeBlockVectorR, eigBlockVectorR)
            pp += sp.dot(activeBlockVectorP, eigBlockVectorP)

            app = sp.dot(activeBlockVectorAR, eigBlockVectorR)
            app += sp.dot(activeBlockVectorAP, eigBlockVectorP)

            bpp = sp.dot(activeBlockVectorBR, eigBlockVectorR)
            bpp += sp.dot(activeBlockVectorBP, eigBlockVectorP)
        else:
            eigBlockVectorX = eigBlockVector[:sizeX]
            eigBlockVectorR = eigBlockVector[sizeX:]

            pp = sp.dot(activeBlockVectorR, eigBlockVectorR)
            app = sp.dot(activeBlockVectorAR, eigBlockVectorR)
            bpp = sp.dot(activeBlockVectorBR, eigBlockVectorR)

        if verbosityLevel > 10:
            print(pp)
            print(app)
            print(bpp)
            pause()

        blockVectorX = sp.dot(blockVectorX, eigBlockVectorX) + pp
        blockVectorAX = sp.dot(blockVectorAX, eigBlockVectorX) + app
        blockVectorBX = sp.dot(blockVectorBX, eigBlockVectorX) + bpp

        blockVectorP, blockVectorAP, blockVectorBP = pp, app, bpp

    aux = blockVectorBX * _lambda[np.newaxis,:]
    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conjugate() * blockVectorR, 0)
    residualNorms = np.sqrt(aux)

    if verbosityLevel > 0:
        print('final eigenvalue:', _lambda)
        print('final residual norms:', residualNorms)

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
    from scipy.sparse import spdiags, speye, issparse
    import time

##     def B( vec ):
##         return vec

    n = 100
    vals = [np.arange(n, dtype=np.float64) + 1]
    A = spdiags(vals, 0, n, n)
    B = speye(n, n)
#    B[0,0] = 0
    B = np.eye(n, n)
    Y = np.eye(n, 3)

#    X = sp.rand( n, 3 )
    xfile = {100: 'X.txt', 1000: 'X2.txt', 10000: 'X3.txt'}
    X = np.fromfile(xfile[n], dtype=np.float64, sep=' ')
    X.shape = (n, 3)

    ivals = [1./vals[0]]

    def precond(x):
        invA = spdiags(ivals, 0, n, n)
        y = invA * x
        if issparse(y):
            y = y.toarray()

        return as2d(y)

    precond = spdiags(ivals, 0, n, n)
#    precond = None
    tt = time.clock()
#    B = None
    eigs, vecs = lobpcg(X, A, B, blockVectorY=Y,
                         M=precond,
                         residualTolerance=1e-4, maxIterations=40,
                         largest=False, verbosityLevel=1)
    print('solution time:', time.clock() - tt)

    print(vecs)
    print(eigs)
