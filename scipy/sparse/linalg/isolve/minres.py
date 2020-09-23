from numpy import sqrt, inner, zeros, inf, finfo, count_nonzero
from numpy.linalg import norm
from scipy.sparse import csc_matrix

from .utils import make_system

__all__ = ['minres']


def minres(A, b, x0=None, shift=0.0, tol=1e-5, maxiter=None,
           M=None, callback=None, show=False, check=False, tolMPD = 0.0, precon=0):
    """
    Use MINimum RESidual iteration to solve Ax=b

    MINRES minimizes norm(A*x - b) for a real symmetric matrix A.  Unlike
    the Conjugate Gradient method, A can be indefinite or singular.

    If shift != 0 then the method solves (A - shift*I)x = b

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The real symmetric N-by-N matrix of the linear system
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : {array, matrix}
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : {array, matrix}
        The converged solution.
    info : integer
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Other Parameters
    ----------------
    x0  : {array, matrix}
        Starting guess for the solution.
    tol : float
        Tolerance to achieve. The algorithm terminates when the relative
        residual is below `tol`.
    maxiter : integer
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse matrix, dense matrix, LinearOperator}
        Preconditioner for A.  The preconditioner should approximate the
        inverse of A.  Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import minres
    >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> A = A + A.T
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = minres(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True

    References
    ----------
    Solution of sparse indefinite systems of linear equations,
        C. C. Paige and M. A. Saunders (1975),
        SIAM J. Numer. Anal. 12(4), pp. 617-629.
        https://web.stanford.edu/group/SOL/software/minres/

    This file is a translation of the following MATLAB implementation:
        https://web.stanford.edu/group/SOL/software/minres/minres-matlab.zip

    """
    A, M, x, b, postprocess = make_system(A, M, x0, b)

    matvec = A.matvec
    psolve = M.matvec

    first = 'Enter minres.   '
    last = 'Exit  minres.   '

    n = A.shape[0]

    if maxiter is None:
        maxiter = 5 * n

    msg = [ ' beta2 = 0.  If M = I, b and x are eigenvectors    ',   # -1
            ' beta1 = 0.  The exact solution is x0              ',   # 0
            ' A solution to Ax = b was found, given tol        ',   # 1
            ' A least-squares solution was found, given tol    ',   # 2
            ' Reasonable accuracy achieved, given eps           ',   # 3
            ' x has converged to an eigenvector                 ',   # 4
            ' Acond has exceeded 0.1/eps                        ',   # 5
            ' The iteration limit was reached                   ',   # 6
            ' A  does not define a symmetric matrix             ',   # 7
            ' M  does not define a symmetric matrix             ',   # 8
            ' M  does not define a pos-def preconditioner       ',   # 9
            ' tol reduced max times on preconditioned system   ']   # 10

    r0 = b
    bnorm = norm(b)
    r0norm = bnorm
    x0norm = 0

    if len(x0) == n:
        y = matvec(x0) - shift*x0
        r0 = b - y
        r0norm = norm(r0)
        x0norm = norm(x0)
    
    if M.shape == 0 :
        precon = 0

    nnzr0 = count_nonzero(r0)

    sparser0 = nnzr0/n <= 0.2

    if show:
        print(first + 'Solution of symmetric Ax = b or (A-shift*I)x = b')
        print(first + 'n        =  %9g     shift  =  %22.14e' % (n,shift))
        print(first + 'itnlim   =  %9g     tol   =  %9.1e' % (maxiter,tol))
        print(first + 'precon   =  %9g     tolMPD   =  %9.1e' % (precon,tolMPD))
        print(first + 'norm(x0) =  %9.1e   norm(b) =%9.1e' % (x0norm,bnorm))
        print(first + 'count_nonzero(r0) =  %9i   norm(r0) =%9.1e' % (nnzr0,r0norm))
        print()

    istop = 0
    itn = 0
    Anorm = 0
    Acond = 0
    rnorm = 0
    xnorm = 0
    betacheck = 0
    dxnorm = 0
    done = False 
    # ynorm = 0

    xtype = x.dtype

    eps = finfo(xtype).eps

    if sparser0 :
        x = csc_matrix((n,1), dtype = float)
    else :
        x = zeros((n,1), dtype = float)

    # Set up y and v for the first Lanczos vector v1.
    # y  =  beta1 P' v1,  where  P = C**(-1).
    # v is really P' v1.
    
    y = r0
    r1 = r0

    if precon > 0 :
        y = psolve(r0)
    
    beta1 = inner(r0,y)
    
    if inner(r0,r0) == 0 :
        betacheck = inf
    else : 
        betacheck = beta1/inner(r0,r0)
    
    if beta1 == 0 :
        istop = 0
        show = True
        done = True
    elif precon > 0 and betacheck < tolMPD :
        istop = 9
        show = True
        done = True
    
    if not done :
        beta1 = sqrt(beta1)
        # see if M is symmetric
        if check and precon :
            r2 = psolve(y)
            s = inner(y,y)
            t = inner(r1,r2)
            z = abs(s-t)
            if z > 0 :
                if s == 0 and t == 0 :
                    z = inf
                else  : 
                    z = z/(abs(s) + abs(t))
                tol = eps**(1.0/3.0)
                if z > tol :
                    istop = 8
                    show = True
                    done = True
                    
        # see if A is symmetric
        if check :
            w = matvec(y)
            r2 = matvec(w)
            s = inner(w,w)
            t = inner(y,r2)
            z = abs(s - t)
            if z > 0 :
                if s == 0 and t == 0 :
                    z = inf
                else  : 
                    z = z/(abs(s) + abs(t))
                tol = eps**(1.0/3.0)
                if z > tol :
                    istop = 7
                    show = True
                    done = True

    # Initialize other quantities
    oldb = 0
    beta = beta1
    dbar = 0
    epsln = 0
    # qrnorm = beta1
    phibar = beta1
    rhs1 = beta1
    rhs2 = 0
    tnorm2 = 0
    gmax = 0
    gmin = finfo(xtype).max
    cs = -1
    sn = 0
    numtol = 1
    tol10 = tol
    r2 = r0
    if sparser0 :
        w = csc_matrix((n,1), dtype = xtype)
        w2 = csc_matrix((n,1), dtype = xtype)
    else :
        w = zeros(n, dtype = xtype)
        w2 = zeros(n, dtype = xtype)

    if show:
        print('   Itn     x(1)     Compatible    LS       norm(A)  cond(A) gbar/|A|')
        x1 = 0
        if x0norm > 0 :
            x1 = x0[0]
        print('%6g %12.5e' % (itn, x1))
        print()
        
    if not done :
        while itn < maxiter:
            itn += 1

            s = 1.0/beta
            v = s*y

            y = matvec(v)
            y = y - shift * v

            if itn >= 2:
                if oldb == 0 :
                    y = inf
                else : 
                    y = y - (beta/oldb)*r1

            alfa = inner(v,y)
            y = y - (alfa/beta)*r2
            r1 = r2
            r2 = y
            if precon > 0 :
                y = psolve(r2)
            oldb = beta
            beta = inner(r2,y)
            if inner(r2,r2) == 0 :
                betacheck = inf
            else :
                betacheck = beta/inner(r2,r2)

            if precon > 0 and betacheck < tolMPD :
                istop = 9
                show = True
                done = True
                x1 = x[0]
                if x0norm > 0 :
                    x1 = x0[0] + x1
                print('%6g %12.5e  betacheck=%10.3e ' % (itn, x1, betacheck))
                raise ValueError('non-symmetric matrix')

            beta = sqrt(beta)
            tnorm2 += alfa**2 + oldb**2 + beta**2

            if itn == 1:
                if beta/beta1 <= 10*eps:
                    istop = -1  # Terminate later

            # Apply previous rotation Qk-1 to get
            #   [deltak epslnk+1] = [cs  sn][dbark    0   ]
            #   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

            oldeps = epsln
            delta = cs * dbar + sn * alfa   # delta1 = 0         deltak
            gbar = sn * dbar - cs * alfa   # gbar 1 = alfa1     gbar k
            epsln = sn * beta     # epsln2 = 0         epslnk+1
            dbar = - cs * beta   # dbar 2 = beta2     dbar k+1
            root = norm([gbar, dbar])
            Arnorm = phibar * root

            # Compute the next plane rotation Qk

            gamma = norm([gbar, beta])       # gammak
            gamma = max(gamma, eps)
            cs = gbar / gamma             # ck
            sn = beta / gamma             # sk
            phi = cs * phibar              # phik
            phibar = sn * phibar              # phibark+1

            # Update  x.

            denom = 1.0/gamma
            w1 = w2
            w2 = w
            w = (v - oldeps*w1 - delta*w2) * denom
            x = x + phi*w

            # Go round again.

            gmax = max(gmax, gamma)
            gmin = min(gmin, gamma)
            z = rhs1 / gamma
            rhs1 = rhs2 - delta*z
            rhs2 = - epsln*z

            # Estimate various norms and test for convergence.

            Anorm = sqrt(tnorm2)
            dxnorm = norm(x)
            rnorm = phibar
            
            Acond = gmax/gmin

            if istop != 0 :
                break 

            # Estimate  cond(A).
            # In this version we look at the diagonals of  R  in the
            # factorization of the lower Hessenberg matrix,  Q * H = R,
            # where H is the tridiagonal matrix from Lanczos with one
            # extra row, beta(k+1) e_k^T.

            # See if any of the stopping criteria are satisfied.
            # In rare cases, istop is already -1 from above (Abar = const*I).

            epsx = (Anorm*dxnorm + beta1)*eps
            epsr = (Anorm*dxnorm + beta1)*tol
            if (Anorm == 0 or dxnorm and bnorm == 0) or (Anorm*dxnorm == -bnorm) :
                test1 = inf
            else :
                test1 = rnorm/(Anorm*dxnorm + bnorm)    # ||r||/(||A|| ||x|| + ||b||)
            if Anorm == 0 or rnorm == -eps :
                test2 = inf
            else :
                test2 = rnorm/(Anorm*(rnorm + eps))     # ||r||/(||A||(||r|| + ||eps||))

            t1 = 1 + test1
            t2 = 1 + test2
            if t2 <= 1 :
                istop = 2
            if t1 <= 1 :
                istop = 1
            if itn >= maxiter:
                istop = 6
            if Acond >= 0.1/eps:
                istop = 5
            if epsx >= beta1:
                istop = 3
            # if rnorm <= epsx   : istop = 2
            # if rnorm <= epsr   : istop = 1
            if test2 <= tol:
                istop = 2
            if test1 <= tol:
                istop = 1  

            if precon > 0 and istop > 0 and istop <= 5:
                if x0norm == 0 :
                    xnorm = dxnorm
                    rnormk = norm(b - matvec(x) - shift*x)
                else :
                    xt = x0 + x
                    xnorm = norm(xt)
                    rnormk = norm(b - matvec(xt) - shift*xt)

                epsr = (Anorm*xnorm + bnorm)*tol10 
                if rnormk <= epsr:
                    istop = 1
                    print('rnormk small enough: %8.1e    epsr =%8.1e' % (rnormk, epsr))
                elif numtol < 5 and tol > eps:
                    numtol = numtol + 1
                    tol = tol/10.0
                    print('\n%7g %2i: tol reduced to%8.1e' %  (itn,numtol, tol))
                    print('   rnormk =%8.1e    epsr =%8.1e' %  (rnormk, epsr))
                    istop = 0
                else :
                    istop = 10

            # See if it is time to print something.

            prnt = False
            if n <= 40:
                prnt = True
            if itn <= 10:
                prnt = True
            if itn >= maxiter-10:
                prnt = True
            if itn % 1000 == 0:
                prnt = True
            if rnorm <= 1.001*epsx:
                prnt = True
            if rnorm <= 1.001*epsr:
                prnt = True
            if Acond <= 1e-2/eps:
                prnt = True
            if istop != 0:
                prnt = True

            if show and prnt:
                x1 = x[0]
                if x0norm > 0:
                    x1 = x0[0] + x1
                str1 = '%6g %12.5e %10.3e' % (itn, x1, test1)
                str2 = ' %10.3e' % (test2,)
                str3 = ' %8.1e %8.1e %8.1e' % (Anorm, Acond, gbar/Anorm)

                print(str1 + str2 + str3)

                if itn % 10 == 0:
                    print()

            if callback is not None:
                callback(x)

            if istop != 0:
                break  # TODO check this

        if show:
            print()
            print(last + ' istop   =  %3g               itn   =%5g' % (istop,itn))
            if precon > 0 :
                print(last + ' Abarnorm   =  %12.4e      Abarcond =  %12.4e' % (Anorm,Acond))
            else : 
                print(last + ' Anorm   =  %12.4e      Acond =  %12.4e' % (Anorm,Acond))
            print(last + ' rnorm   =  %12.4e      xnorm =  %12.4e' % (rnorm,xnorm))
            print(last + ' Arnorm  =  %12.4e' % (Arnorm,))
            print(last + msg[istop+1])

        if istop == 6:
            info = maxiter
        else:
            info = 0

        return (postprocess(x),info)

if __name__ == '__main__':
    from numpy import arange
    from scipy.sparse import spdiags

    n = 10

    residuals = []

    def cb(x):
        residuals.append(norm(b - A*x))

    # A = poisson((10,),format='csr')
    A = spdiags([arange(1,n+1,dtype=float)], [0], n, n, format='csr')
    M = spdiags([1.0/arange(1,n+1,dtype=float)], [0], n, n, format='csr')
    A.psolve = M.matvec
    b = zeros(A.shape[0])
    x = minres(A,b,tol=1e-12,maxiter=None,callback=cb)
    # x = cg(A,b,x0=b,tol=1e-12,maxiter=None,callback=cb)[0]
