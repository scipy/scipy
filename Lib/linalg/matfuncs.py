#
# Author: Travis Oliphant, March 2002
#

__all__ = ['expm','expm2','expm3','cosm','sinm','tanm','coshm','sinhm',
           'tanhm','logm','funm','signm','sqrtm']

from scipy_base import asarray, Inf, dot, floor, log2, eye, diag, exp, \
     product, logical_not, ravel, transpose, conjugate, \
     cast, log, ogrid, isfinite, imag, real, absolute, amax, sign, \
     isfinite, sqrt
from Matrix import Matrix as mat
import scipy_base as sb
from basic import solve, LinAlgError, inv, norm, triu, all_mat
from decomp import eig, schur, rsf2csf, orth, eigvals, svd

eps = sb.limits.double_epsilon
feps = sb.limits.float_epsilon

def expm(A,q=7):
    """Compute the matrix exponential using Pade approximation of order q.
    """
    A = asarray(A)
    ss = A.spacesaver()
    if A.typecode() in ['f', 'F']:
        A.savespace(1)
    else:
        A.savespace(0)
 
    # Scale A so that norm is < 1/2
    val = log2(norm(A,Inf))
    e = int(floor(val))
    j = max(0,e+1)
    A = A / 2.0**j
 
    # Pade Approximation for exp(A)
    X = A
    c = 1.0/2
    N = eye(*A.shape) + c*A
    D = eye(*A.shape) - c*A
    for k in range(2,q+1):
        c = c * (q-k+1) / (k*(2*q-k+1))
        X = dot(A,X)
        cX = c*X
        N = N + cX
        if not k % 2:
            D = D + cX;
        else:
            D = D - cX;
    F = solve(D,N)
    for k in range(1,j+1):
        F = dot(F,F)
    A.savespace(ss)
    return F

def expm2(A):
    """Compute the matrix exponential using eigenvalue decomposition.
    """
    A = asarray(A)
    t = A.typecode()
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'    
    s,vr = eig(A)
    vri = inv(vr)
    return dot(dot(vr,diag(exp(s))),vri).astype(t)

def expm3(A,q=20):
    """Compute the matrix exponential using a Taylor series.of order q.
    """
    A = asarray(A)
    t = A.typecode()
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'
    A = mat(A)
    eA = eye(*A.shape,**{'typecode':t})
    trm = mat(eA)
    castfunc = cast[t]
    for k in range(1,q):
        trm *= A / castfunc(k)
        eA += trm
    return eA

_array_precision = {'i': 1, 'l': 1, 'f': 0, 'd': 1, 'F': 0, 'D': 1}
def toreal(arr,tol=None):
    """Return as real array if imaginary part is small.
    """
    if tol is None:
        tol = {0:feps*1e3, 1:eps*1e6}[_array_precision[arr.typecode()]]
    if (arr.typecode() in ['F', 'D']) and \
       sb.allclose(arr.imag, 0.0, atol=tol):
        arr = arr.real
    return arr
 
def cosm(A):
    """matrix cosine.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(0.5*(expm(1j*A) + expm(-1j*A)))
    else:
        return 0.5*(expm(1j*A) + expm(-1j*A))
        
            
def sinm(A):
    """matrix sine.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(-0.5j*(expm(1j*A) - expm(-1j*A)))
    else:
        return -0.5j*(expm(1j*A) - expm(-1j*A))
    
def tanm(A): 
    """matrix tangent.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(solve(cosm(A), sinm(A)))
    else:        
        return solve(cosm(A), sinm(A))
 
def coshm(A): 
    """matrix hyperbolic cosine.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(0.5*(expm(A) + expm(-A)))
    else:
        return 0.5*(expm(A) + expm(-A))
 
def sinhm(A):
    """matrix hyperbolic sine.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(0.5*(expm(A) - expm(-A)))
    else:
        return 0.5*(expm(A) - expm(-A))
    
def tanhm(A):
    """matrix hyperbolic tangent.
    """
    A = asarray(A)
    if A.typecode() not in ['F','D']:
        return toreal(solve(coshm(A), sinhm(A)))
    else:
        return solve(coshm(A), sinhm(A))
 
def funm(A,func,disp=1):
    """matrix function for arbitrary callable object func.
    """
    # func should take a vector of arguments (see vectorize if
    #  it needs wrapping.
 
    # Perform Shur decomposition (lapack ?gees)
    A = asarray(A)
    if len(A.shape)!=2:
        raise ValueError, "Non-matrix input to matrix function."    
    if A.typecode() in ['F', 'D']:
        cmplx_type = 1
    else:
        cmplx_type = 0
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    n,n = T.shape
    F = diag(func(diag(T)))  # apply function to diagonal elements
    F = F.astype(T.typecode()) # e.g. when F is real but T is complex

    minden = abs(T[0,0])

    # implement Algorithm 11.1.1 from Golub and Van Loan
    #                 "matrix Computations." 
    for p in range(1,n):
        for i in range(1,n-p+1):
            j = i + p
            s = T[i-1,j-1] * (F[j-1,j-1] - F[i-1,i-1])
            ksl = slice(i,j-1)
            s = s + dot(T[i-1,ksl],F[ksl,j-1]) - dot(F[i-1,ksl],T[ksl,j-1])
            den = T[j-1,j-1] - T[i-1,i-1]
            if den != 0.0:
                s = s / den
            F[i-1,j-1] = s
            minden = min(minden,abs(den))

    F = dot(dot(Z, F),transpose(conjugate(Z)))
    if not cmplx_type:
        F = toreal(F)

    tol = {0:feps, 1:eps}[_array_precision[F.typecode()]]
    if minden == 0.0:
        minden = tol
    err = min(1, max(tol,(tol/minden)*norm(triu(T,1),1)))
    if product(ravel(logical_not(isfinite(F)))):
        err = Inf
    if disp:
        if err > 1000*tol:
            print "Result may be inaccurate, approximate err =", err
        return F
    else:
        return F, err

def logm(A,disp=1):
    """Matrix logarithm, inverse of expm."""
    # Compute using general funm but then use better error estimator and
    #   make one step in improving estimate using a rotation matrix.
    A = mat(asarray(A))
    F, errest = funm(A,log,disp=0)
    errtol = 1000*eps
    # Only iterate if estimate of error is too large.
    if errest >= errtol:
        # Use better approximation of error
        errest = norm(expm(F)-A,1) / norm(A,1)
        if not isfinite(errest) or errest >= errtol:
            N,N = A.shape
            X,Y = ogrid[1:N+1,1:N+1]
            R = mat(orth(eye(N,typecode='d')+X+Y))
            F, dontcare = funm(R*A*R.H,log,disp=0)
            F = R.H*F*R
            if (norm(imag(F),1)<=1000*errtol*norm(F,1)):
                F = mat(real(F))
            E = mat(expm(F))
            temp = mat(solve(E.T,(E-A).T))
            F = F - temp.T
            errest = norm(expm(F)-A,1) / norm(A,1)
    if disp:
        if not isfinite(errest) or errest >= errtol:
            print "Result may be inaccurate, approximate err =", errest
        return F
    else:
        return F, errest

def signm(a,disp=1):
    """matrix sign"""
    def rounded_sign(x):
        rx = real(x)
        if rx.typecode()=='f':
            c =  1e3*feps*amax(x)
        else:
            c =  1e3*eps*amax(x)
        return sign( (absolute(rx) > c) * rx )
    result,errest = funm(a, rounded_sign, disp=0)
    errtol = {0:1e3*feps, 1:1e3*eps}[_array_precision[result.typecode()]]
    if errest < errtol:
        return result

    # Handle signm of defective matrices:

    # See "E.D.Denman and J.Leyva-Ramos, Appl.Math.Comp.,
    # 8:237-250,1981" for how to improve the following (currently a
    # rather naive) iteration process:

    a = asarray(a)
    #a = result # sometimes iteration converges faster but where??

    # Shifting to avoid zero eigenvalues. How to ensure that shifting does
    # not change the spectrum too much?
    vals = svd(a,compute_uv=0)
    max_sv = sb.amax(vals)
    #min_nonzero_sv = vals[(vals>max_sv*errtol).tolist().count(1)-1]
    #c = 0.5/min_nonzero_sv
    c = 0.5/max_sv
    S0 = a + c*sb.identity(a.shape[0])
    prev_errest = errest
    for i in range(100):
        iS0 = inv(S0)
        S0 = 0.5*(S0 + iS0)
        Pp=0.5*(dot(S0,S0)+S0)
        errest = norm(dot(Pp,Pp)-Pp,1)
        if errest < errtol or prev_errest==errest:
            break
        prev_errest = errest
    if disp:
        if not isfinite(errest) or errest >= errtol:
            print "Result may be inaccurate, approximate err =", errest
        return S0
    else:
        return S0, errest

def sqrtm(A,disp=1):
    """Matrix square root

    If disp is non-zero display warning if singular matrix.
    If disp is zero then return residual ||A-X*X||_F / ||A||_F

    Uses algorithm by Nicholas J. Higham
    """
    A = asarray(A)
    if len(A.shape)!=2:
        raise ValueError, "Non-matrix input to matrix function."    
    if A.typecode() in ['F', 'D']:
        cmplx_type = 1
    else:
        cmplx_type = 0
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    n,n = T.shape

    R = sb.zeros((n,n),T.typecode())
    for j in range(n):
        R[j,j] = sqrt(T[j,j])
        for i in range(j-1,-1,-1):
            s = 0
            for k in range(i+1,j):
                s = s + R[i,k]*R[k,j]
            R[i,j] = (T[i,j] - s)/(R[i,i] + R[j,j])

    R, Z = all_mat(R,Z)
    X = (Z * R * Z.H)

    if disp:
        nzeig = sb.any(sb.diag(T)==0)
        if nzeig:
            print "Matrix is singular and may not have a square root."
        return X.A
    else:
        arg2 = norm(X*X - A,'fro')**2 / norm(A,'fro')
        return X.A, arg2

    





