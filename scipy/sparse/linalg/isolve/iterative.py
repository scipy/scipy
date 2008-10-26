## Automatically adapted for scipy Oct 18, 2005 by


# Iterative methods using reverse-communication raw material
#   These methods solve
#   Ax = b  for x

#   where A must have A.matvec(x,*args) defined
#    or be a numeric array


__all__ = ['bicg','bicgstab','cg','cgs','gmres','qmr']

import _iterative
import numpy as np

from scipy.sparse.linalg.interface import LinearOperator
from utils import make_system

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}

def bicg(A, b, x0=None, tol=1e-5, maxiter=None, xtype=None, M=None, callback=None):
    """Use BIConjugate Gradient iteration to solve A x = b

    Inputs:

    A --   An array or an object with matvec(x) and rmatvec(x) methods
           to represent A * x and A^H * x respectively.  May also have
           psolve(b) and rpsolve(b) methods for representing solutions
           to the preconditioning equations M * x = b and
           M^H * x = b respectively.
    b --   An n-length vector

    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(xk), where xk is the
                current parameter vector.
    """
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec, rmatvec = A.matvec, A.rmatvec
    psolve, rpsolve = M.matvec, M.rmatvec
    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'bicgrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = np.zeros(6*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice2] *= sclr2
            work[slice2] += sclr1*rmatvec(work[slice1])
        elif (ijob == 3):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 4):
            work[slice1] = rpsolve(work[slice2])
        elif (ijob == 5):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 6):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return postprocess(x), info


def bicgstab(A, b, x0=None, tol=1e-5, maxiter=None, xtype=None, M=None, callback=None):
    """Use BIConjugate Gradient STABilized iteration to solve A x = b

    Inputs:

    A --   An array or an object with a matvec(x) method
           to represent A * x.  May also have a psolve(b) methods for
           representing solution to the preconditioning equation
           M * x = b.
    b --   An n-length vector

    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(xk), where xk is the
                current parameter vector.
    """
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'bicgstabrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = np.zeros(7*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            if psolve is None:
                psolve = get_psolve(A)
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return postprocess(x), info


def cg(A, b, x0=None, tol=1e-5, maxiter=None, xtype=None, M=None, callback=None):
    """Use Conjugate Gradient iteration to solve A x = b (A^H = A)

    Inputs:

    A --   An array or an object with a matvec(x) method
           to represent A * x.  May also have a psolve(b) methods for
           representing solution to the preconditioning equation
           M * x = b.
    b --   An n-length vector


    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(xk), where xk is the
                current parameter vector.
    """
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'cgrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = np.zeros(4*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return postprocess(x), info


def cgs(A, b, x0=None, tol=1e-5, maxiter=None, xtype=None, M=None, callback=None):
    """Use Conjugate Gradient Squared iteration to solve A x = b

    Inputs:

    A --   An array or an object with a matvec(x) method
           to represent A * x.  May also have a psolve(b) methods for
           representing solution to the preconditioning equation
           M * x = b.
    b --   An n-length vector


    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(xk), where xk is the
                current parameter vector.
    """
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'cgsrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = np.zeros(7*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return postprocess(x), info


def gmres(A, b, x0=None, tol=1e-5, restrt=None, maxiter=None, xtype=None, M=None, callback=None):
    """Use Generalized Minimal RESidual iteration to solve A x = b

    Inputs:

    A --   An array or an object with a matvec(x) method
           to represent A * x.  May also have a psolve(b) methods for
           representing solution to the preconditioning equation
           M * x = b.
    b --   An n-length vector


    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    restrt -- (n) When to restart (change this to get faster performance -- but
                   may not converge).
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(rk), where rk is the
                the current relative residual
    """
    A,M,x,b,postprocess = make_system(A,M,x0,b,xtype)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'gmresrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    if restrt is None:
        restrt = n
    resid = tol
    ndx1 = 1
    ndx2 = -1
    work  = np.zeros((6+restrt)*n,dtype=x.dtype)
    work2 = np.zeros((restrt+1)*(2*restrt+2),dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    old_ijob = ijob
    first_pass = True
    resid_ready = False
    iter_num = 1
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, restrt, work, work2, iter_, resid, info, ndx1, ndx2, ijob)
        #if callback is not None and iter_ > olditer:
        #    callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1): # gmres success, update last residual
            if resid_ready and callback is not None:
                callback(resid)
                resid_ready = False

            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
            if not first_pass and old_ijob==3:
                resid_ready = True

            first_pass = False
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
            if resid_ready and callback is not None:
                callback(resid)
                resid_ready = False
                iter_num = iter_num+1

        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)

        old_ijob = ijob
        ijob = 2

        if iter_num > maxiter:
            break

    return postprocess(x), info


def qmr(A, b, x0=None, tol=1e-5, maxiter=None, xtype=None, M1=None, M2=None, callback=None):
    """Use Quasi-Minimal Residual iteration to solve A x = b

    Inputs:

    A --   An array or an object with matvec(x) and rmatvec(x) methods
           to represent A * x and A^H * x respectively.  May also have
           psolve(b,<which>) and rpsolve(b,<which>) methods for
           representing solutions to the preconditioning equations
           M * x = b and M^H * x = b respectively.   The <which> argument
           may be given to specify 'left' or 'right' preconditioning.
    b --   An n-length vector

    Outputs:

    x  --  The converged solution
    info -- output result
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Optional Inputs:

    x0  -- (0) default starting guess.
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations
    xtype  --  The type of the result.  If None, then it will be
               determined from A.dtype.char and b.  If A does not have a
               typecode method then it will compute A.matvec(x0) to get a
               typecode.   To save the extra computation when A does not
               have a typecode attribute use xtype=0 for the same type as
               b or use xtype='f','d','F',or 'D'
    callback -- an optional user-supplied function to call after each
                iteration.  It is called as callback(xk), where xk is the
                current parameter vector.
    """
    A_ = A
    A,M,x,b,postprocess = make_system(A,None,x0,b,xtype)

    if M1 is None and M2 is None:
        if hasattr(A_,'psolve'):
            def left_psolve(b):
                return A_.psolve(b,'left')
            def right_psolve(b):
                return A_.psolve(b,'right')
            def left_rpsolve(b):
                return A_.rpsolve(b,'left')
            def right_rpsolve(b):
                return A_.rpsolve(b,'right')
            M1 = LinearOperator(A.shape, matvec=left_psolve, rmatvec=left_rpsolve)
            M2 = LinearOperator(A.shape, matvec=right_psolve, rmatvec=right_rpsolve)
        else:
            def id(b):
                return b
            M1 = LinearOperator(A.shape, matvec=id, rmatvec=id)
            M2 = LinearOperator(A.shape, matvec=id, rmatvec=id)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    ltr = _type_conv[x.dtype.char]
    revcom   = getattr(_iterative, ltr + 'qmrrevcom')
    stoptest = getattr(_iterative, ltr + 'stoptest2')

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = np.zeros(11*n,x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.matvec(work[slice1])
        elif (ijob == 2):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.rmatvec(work[slice1])
        elif (ijob == 3):
            work[slice1] = M1.matvec(work[slice2])
        elif (ijob == 4):
            work[slice1] = M2.matvec(work[slice2])
        elif (ijob == 5):
            work[slice1] = M1.rmatvec(work[slice2])
        elif (ijob == 6):
            work[slice1] = M2.rmatvec(work[slice2])
        elif (ijob == 7):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.matvec(x)
        elif (ijob == 8):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return postprocess(x), info
