
# Iterative methods using reverse-communication raw material
#   These methods solve
#   Ax = b  for x

#   where A must have A.matvec(x,*args) defined
#    or be a numeric array


__all__ = ['bicg','bicgstab','cg','cgs','gmres','qmr'] 
import _iterative
import scipy_base as sb

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}

class get_matvec:
    methname = 'matvec'
    def __init__(self, obj, *args):
        self.obj = obj
        self.args = args
        if isinstance(obj, sb.ArrayType):
            self.callfunc = self.type1
            return
        meth = getattr(obj,self.methname,None)
        if not callable(meth):
            raise ValueError, "Object must be an array "\
                  "or have a callable %s attribute." % (self.methname,)

        self.obj = meth
        self.callfunc = self.type2

    def __call__(self, x):
        return self.callfunc(x)

    def type1(self, x):
        return sb.dot(self.obj, x)

    def type2(self, x):
        return self.obj(x,*self.args)

class get_rmatvec(get_matvec):
    methname = 'rmatvec'
    def type1(self, x):
        return sb.dot(x, self.obj)

class get_psolve:
    methname = 'psolve'
    def __init__(self, obj, *args):
        self.obj = obj
        self.args = args
        meth = getattr(obj,self.methname,None)   
        if meth is None:  # no preconditiong available
            self.callfunc = self.type1
            return

        if not callable(meth):
            raise ValueError, "Preconditioning method %s "\
                  "must be callable." % (self.methname,)

        self.obj = meth
        self.callfunc = self.type2

    def __call__(self, x):
        return self.callfunc(x)
    
    def type1(self, x):
        return x

    def type2(self, x):
        return self.obj(x,*self.args)

class get_rpsolve(get_psolve):
    methname = 'rpsolve'

class get_psolveq(get_psolve):

    def __call__(self, x, which):
        return self.callfunc(x, which)

    def type1(self, x, which):
        return x

    def type2(self, x, which):
        return self.obj(x,which,*self.args)

class get_rpsolveq(get_psolveq):
    methname = 'rpsolve'
    
def bicg(A,b,x0=None,tol=1e-5,maxiter=None):
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

    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve, rmatvec, rpsolve = (None,)*4
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'bicgrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros(6*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
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
            if rmatvec is None:
                rmatvec = get_rmatvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*rmatvec(work[slice1])
        elif (ijob == 3):
            if psolve is None:
                psolve = get_psolve(A)
            work[slice1] = psolve(work[slice2])
        elif (ijob == 4):
            if rpsolve is None:
                rpsolve = get_rpsolve(A)
            work[slice1] = rpsolve(work[slice2])
        elif (ijob == 5):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 6):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return x, info

def bicgstab(A,b,x0=None,tol=1e-5,maxiter=None):
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

    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve = (None,)*2
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'bicgstabrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros(7*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
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

    return x, info

def cg(A,b,x0=None,tol=1e-5,maxiter=None):
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

    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve = (None,)*2
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'cgrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros(4*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
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

    return x, info


def cgs(A,b,x0=None,tol=1e-5,maxiter=None):
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

    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve = (None,)*2
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'cgsrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros(7*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
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

    return x, info

# not working yet.
def gmres(A,b,restrt=None,x0=None,tol=1e-5,maxiter=None):
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

    restrt -- (n) When to restart (change this to get faster performance -- but
                   may not converge). 
    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve = (None,)*2
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'gmresrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    if restrt is None:
        restrt = n
    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros((6+restrt)*n,typ)
    work2 = sb.zeros((restrt+1)*(2*restrt+2),typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, restrt, work, work2, iter_, resid, info, ndx1, ndx2, ijob)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            break
        elif (ijob == 1):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)            
        elif (ijob == 2):
            if psolve is None:
                psolve = get_psolve(A)
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return x, info


def qmr(A,b,x0=None,tol=1e-5,maxiter=None):
    """Use Quasi-Minimal Residucal iteration to solve A x = b

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

    x0  -- (0) default starting guess
    tol -- (1e-5) relative tolerance to achieve
    maxiter -- (10*n) maximum number of iterations

    """
    b = sb.asarray(b)
    n = len(b)
    typ = b.typecode()
    if maxiter is None:
        maxiter = n*10

    x = x0
    if x is None:
        x = sb.zeros(n,typ)
        
    matvec, psolve, rmatvec, rpsolve = (None,)*4
    ltr = _type_conv[typ]
    revcom = _iterative.__dict__[ltr+'qmrrevcom']
    stoptest = _iterative.__dict__[ltr+'stoptest2']

    resid = tol
    ndx1 = 1
    ndx2 = -1
    work = sb.zeros(11*n,typ)
    ijob = 1
    info = 0
    ftflag = True
    bnrm2 = -1.0
    iter_ = maxiter
    while 1:
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
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
            if rmatvec is None:
                rmatvec = get_rmatvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*rmatvec(work[slice1])
        elif (ijob == 3):
            if psolve is None:
                psolve = get_psolveq(A)
            work[slice1] = psolve(work[slice2],'left')
        elif (ijob == 4):
            if psolve is None:
                psolve = get_psolveq(A)
            work[slice1] = psolve(work[slice2],'right')            
        elif (ijob == 5):
            if rpsolve is None:
                rpsolve = get_rpsolveq(A)
            work[slice1] = rpsolve(work[slice2],'left')
        elif (ijob == 6):
            if rpsolve is None:
                rpsolve = get_rpsolveq(A)
            work[slice1] = rpsolve(work[slice2],'right')
        elif (ijob == 7):
            if matvec is None:
                matvec = get_matvec(A)
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 8):
            if ftflag:
                info = -1
                ftflag = False
            bnrm2, resid, info = stoptest(work[slice1], b, bnrm2, tol, info)
        ijob = 2

    return x, info


