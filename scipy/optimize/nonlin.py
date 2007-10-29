"""
Nonlinear solvers
=================

These solvers find x for which F(x)=0. Both x and F is multidimensional.

They accept the user defined function F, which accepts a python tuple x and it
should return F(x), which can be either a tuple, or numpy array.

Example:

def F(x):
    "Should converge to x=[0,0,0,0,0]"
    import numpy
    d = numpy.array([3,2,1.5,1,0.5])
    c = 0.01
    return -d*numpy.array(x)-c*numpy.array(x)**3

from scipy import optimize
x = optimize.broyden2(F,[1,1,1,1,1])

All solvers have the parameter iter (the number of iterations to compute), some
of them have other parameters of the solver, see the particular solver for
details.

 A collection of general-purpose nonlinear multidimensional solvers.

   broyden1            --  Broyden's first method - is a quasi-Newton-Raphson
                           method for updating an approximate Jacobian and then
                           inverting it
   broyden2            --  Broyden's second method - the same as broyden1, but
                           updates the inverse Jacobian directly
   broyden3            --  Broyden's second method - the same as broyden2, but
                           instead of directly computing the inverse Jacobian,
                           it remembers how to construct it using vectors, and
                           when computing inv(J)*F, it uses those vectors to
                           compute this product, thus avoding the expensive NxN
                           matrix multiplication.
   broyden_generalized --  Generalized Broyden's method, the same as broyden2,
                           but instead of approximating the full NxN Jacobian,
                           it construct it at every iteration in a way that
                           avoids the NxN matrix multiplication.  This is not
                           as precise as broyden3.
   anderson            --  extended Anderson method, the same as the
                           broyden_generalized, but added w_0^2*I to before
                           taking inversion to improve the stability
   anderson2           --  the Anderson method, the same as anderson, but
                           formulated differently

 The broyden2 is the best. For large systems, use broyden3. excitingmixing is
 also very effective. There are some more solvers implemented (see their
 docstrings), however, those are of mediocre quality.


 Utility Functions

   norm --  Returns an L2 norm of the vector

"""

import math

import numpy

def mlog(x):
    if x==0.:
        return 13
    else:
        return math.log(x)

def norm(v):
    """Returns an L2 norm of the vector."""
    return math.sqrt(numpy.sum((numpy.array(v)**2).flat))

def myF(F,xm):
    return numpy.matrix(F(tuple(xm.flat))).T

def difference(a,b):
    m=0.
    for x,y in zip(a,b):
        m+=(x-y)**2
    return math.sqrt(m)

def sum(a,b):
    return [ai+bi for ai,bi in zip(a,b)]

def mul(C,b):
    return [C*bi for bi in b]

def solve(A,b):
    """Solve Ax=b, returns x"""
    try:
        from scipy import linalg
        return linalg.solve(A,b)
    except:
        return A.I*b

def broyden2(F, xin, iter=10, alpha=0.4, verbose = False):
    """Broyden's second method.

    Updates inverse Jacobian by an optimal formula.
    There is NxN matrix multiplication in every iteration.

    The best norm |F(x)|=0.003 achieved in ~20 iterations.

    Recommended.
    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    Gm=-alpha*numpy.matrix(numpy.identity(len(xin)))
    for n in range(iter):
        deltaxm=-Gm*Fxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
        Gm=Gm+(deltaxm-Gm*deltaFxm)*deltaFxm.T/norm(deltaFxm)**2
        if verbose:
            print "%d:  |F(x)|=%.3f"%(n+1, norm(Fxm))
    return xm.flat

def broyden3(F, xin, iter=10, alpha=0.4, verbose = False):
    """Broyden's second method.

    Updates inverse Jacobian by an optimal formula.
    The NxN matrix multiplication is avoided.

    The best norm |F(x)|=0.003 achieved in ~20 iterations.

    Recommended.
    """
    zy=[]
    def updateG(z,y):
        "G:=G+z*y.T"
        zy.append((z,y))
    def Gmul(f):
        "G=-alpha*1+z*y.T+z*y.T ..."
        s=-alpha*f
        for z,y in zy:
            s=s+z*(y.T*f)
        return s
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
#    Gm=-alpha*numpy.matrix(numpy.identity(len(xin)))
    for n in range(iter):
        #deltaxm=-Gm*Fxm
        deltaxm=Gmul(-Fxm)
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
        #Gm=Gm+(deltaxm-Gm*deltaFxm)*deltaFxm.T/norm(deltaFxm)**2
        updateG(deltaxm-Gmul(deltaFxm),deltaFxm/norm(deltaFxm)**2)
        if verbose:
            print "%d:  |F(x)|=%.3f"%(n+1, norm(Fxm))
    return xm.flat

def broyden_generalized(F, xin, iter=10, alpha=0.1, M=5, verbose = False):
    """Generalized Broyden's method.

    Computes an approximation to the inverse Jacobian from the last M
    interations. Avoids NxN matrix multiplication, it only has MxM matrix
    multiplication and inversion.

    M=0 .... linear mixing
    M=1 .... Anderson mixing with 2 iterations
    M=2 .... Anderson mixing with 3 iterations
    etc.
    optimal is M=5

    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    G0=-alpha
    dxm=[]
    dFxm=[]
    for n in range(iter):
        deltaxm=-G0*Fxm
        if M>0:
            MM=min(M,n)
            for m in range(n-MM,n):
                deltaxm=deltaxm-(float(gamma[m-(n-MM)])*dxm[m]-G0*dFxm[m])
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1

        if M>0:
            dxm.append(deltaxm)
            dFxm.append(deltaFxm)
            MM=min(M,n+1)
            a=numpy.matrix(numpy.empty((MM,MM)))
            for i in range(n+1-MM,n+1):
                for j in range(n+1-MM,n+1):
                    a[i-(n+1-MM),j-(n+1-MM)]=dFxm[i].T*dFxm[j]

            dFF=numpy.matrix(numpy.empty(MM)).T
            for k in range(n+1-MM,n+1):
                dFF[k-(n+1-MM)]=dFxm[k].T*Fxm
            gamma=a.I*dFF

        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm.flat

def anderson(F, xin, iter=10, alpha=0.1, M=5, w0=0.01, verbose = False):
    """Extended Anderson method.

    Computes an approximation to the inverse Jacobian from the last M
    interations. Avoids NxN matrix multiplication, it only has MxM matrix
    multiplication and inversion.

    M=0 .... linear mixing
    M=1 .... Anderson mixing with 2 iterations
    M=2 .... Anderson mixing with 3 iterations
    etc.
    optimal is M=5

    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    dxm=[]
    dFxm=[]
    for n in range(iter):
        deltaxm=alpha*Fxm
        if M>0:
            MM=min(M,n)
            for m in range(n-MM,n):
                deltaxm=deltaxm-(float(gamma[m-(n-MM)])*dxm[m]+alpha*dFxm[m])
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1

        if M>0:
            dxm.append(deltaxm)
            dFxm.append(deltaFxm)
            MM=min(M,n+1)
            a=numpy.matrix(numpy.empty((MM,MM)))
            for i in range(n+1-MM,n+1):
                for j in range(n+1-MM,n+1):
                    if i==j: wd=w0**2
                    else: wd=0
                    a[i-(n+1-MM),j-(n+1-MM)]=(1+wd)*dFxm[i].T*dFxm[j]

            dFF=numpy.matrix(numpy.empty(MM)).T
            for k in range(n+1-MM,n+1):
                dFF[k-(n+1-MM)]=dFxm[k].T*Fxm
            gamma=solve(a,dFF)
#            print gamma

        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm.flat

def anderson2(F, xin, iter=10, alpha=0.1, M=5, w0=0.01, verbose = False):
    """Anderson method.

    M=0 .... linear mixing
    M=1 .... Anderson mixing with 2 iterations
    M=2 .... Anderson mixing with 3 iterations
    etc.
    optimal is M=5

    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    dFxm=[]
    for n in range(iter):
        deltaxm=Fxm
        if M>0:
            MM=min(M,n)
            for m in range(n-MM,n):
                deltaxm=deltaxm+float(theta[m-(n-MM)])*(dFxm[m]-Fxm)
        deltaxm=deltaxm*alpha
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1

        if M>0:
            dFxm.append(Fxm-deltaFxm)
            MM=min(M,n+1)
            a=numpy.matrix(numpy.empty((MM,MM)))
            for i in range(n+1-MM,n+1):
                for j in range(n+1-MM,n+1):
                    if i==j: wd=w0**2
                    else: wd=0
                    a[i-(n+1-MM),j-(n+1-MM)]= \
                        (1+wd)*(Fxm-dFxm[i]).T*(Fxm-dFxm[j])

            dFF=numpy.matrix(numpy.empty(MM)).T
            for k in range(n+1-MM,n+1):
                dFF[k-(n+1-MM)]=(Fxm-dFxm[k]).T*Fxm
            theta=solve(a,dFF)
#            print gamma

        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm.flat

def broyden_modified(F, xin, iter=10, alpha=0.35, w0=0.01, wl=5, verbose = False):
    """Modified Broyden's method.

    Updates inverse Jacobian using information from all the iterations and
    avoiding the NxN matrix multiplication. The problem is with the weights,
    it converges the same or worse than broyden2 or broyden_generalized

    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    G0=alpha
    w=[]
    u=[]
    dFxm=[]
    for n in range(iter):
        deltaxm=G0*Fxm
        for i in range(n):
            for j in range(n):
                deltaxm-=w[i]*w[j]*betta[i,j]*u[j]*(dFxm[i].T*Fxm)
        xm+=deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1

        w.append(wl/norm(Fxm))

        u.append((G0*deltaFxm+deltaxm)/norm(deltaFxm))
        dFxm.append(deltaFxm/norm(deltaFxm))
        a=numpy.matrix(numpy.empty((n+1,n+1)))
        for i in range(n+1):
            for j in range(n+1):
                a[i,j]=w[i]*w[j]*dFxm[j].T*dFxm[i]
        betta=(w0**2*numpy.matrix(numpy.identity(n+1))+a).I

        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm.flat

def broyden1(F, xin, iter=10, alpha=0.1, verbose = False):
    """Broyden's first method.

    Updates Jacobian and computes inv(J) by a matrix inversion at every
    iteration. It's very slow.

    The best norm |F(x)|=0.005 achieved in ~45 iterations.
    """
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    Jm=-1/alpha*numpy.matrix(numpy.identity(len(xin)))

    for n in range(iter):
        deltaxm=solve(-Jm,Fxm)
        #!!!! What the fuck?!
        #xm+=deltaxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
        Jm=Jm+(deltaFxm-Jm*deltaxm)*deltaxm.T/norm(deltaxm)**2
        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm.flat

def broyden1_modified(F, xin, iter=10, alpha=0.1, verbose = False):
    """Broyden's first method, modified by O. Certik.

    Updates inverse Jacobian using some matrix identities at every iteration,
    its faster then newton_slow, but still not optimal.

    The best norm |F(x)|=0.005 achieved in ~45 iterations.
    """
    def inv(A,u,v):

        #interesting is that this
        #return (A.I+u*v.T).I
        #is more stable than
        #return A-A*u*v.T*A/float(1+v.T*A*u)
        Au=A*u
        return A-Au*(v.T*A)/float(1+v.T*Au)
    xm=numpy.matrix(xin).T
    Fxm=myF(F,xm)
    Jm=alpha*numpy.matrix(numpy.identity(len(xin)))
    for n in range(iter):
        deltaxm=Jm*Fxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
#        print "-------------",norm(deltaFxm),norm(deltaxm)
        deltaFxm/=norm(deltaxm)
        deltaxm/=norm(deltaxm)
        Jm=inv(Jm+deltaxm*deltaxm.T*Jm,-deltaFxm,deltaxm)

        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm

def vackar(F, xin, iter=10, alpha=0.1, verbose = False):
    """J=diag(d1,d2,...,dN)

    The best norm |F(x)|=0.005 achieved in ~110 iterations.
    """
    def myF(F,xm):
        return numpy.array(F(tuple(xm.flat))).T
    xm=numpy.array(xin)
    Fxm=myF(F,xm)
    d=1/alpha*numpy.ones(len(xin))
    Jm=numpy.matrix(numpy.diag(d))

    for n in range(iter):
        deltaxm=1/d*Fxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
        d=d-(deltaFxm+d*deltaxm)*deltaxm/norm(deltaxm)**2
        if verbose:
            print "%d:  |F(x)|=%.3f"%(n, norm(Fxm))
    return xm

def linearmixing(F,xin, iter=10, alpha=0.1, verbose = False):
    """J=-1/alpha

    The best norm |F(x)|=0.005 achieved in ~140 iterations.
    """
    def myF(F,xm):
        return numpy.array(F(tuple(xm.flat))).T
    xm=numpy.array(xin)
    Fxm=myF(F,xm)
    for n in range(iter):
        deltaxm=alpha*Fxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        Fxm=Fxm1
        if verbose:
            print "%d: |F(x)|=%.3f" %(n,norm(Fxm))

    return xm

def excitingmixing(F,xin,iter=10,alpha=0.1,alphamax=1.0, verbose = False):
    """J=-1/alpha

    The best norm |F(x)|=0.005 achieved in ~140 iterations.
    """
    def myF(F,xm):
        return numpy.array(F(tuple(xm.flat))).T
    xm=numpy.array(xin)
    beta=numpy.array([alpha]*len(xm))
    Fxm=myF(F,xm)
    for n in range(iter):
        deltaxm=beta*Fxm
        xm=xm+deltaxm
        Fxm1=myF(F,xm)
        deltaFxm=Fxm1-Fxm
        for i in range(len(xm)):
            if Fxm1[i]*Fxm[i] > 0:
                beta[i]=beta[i]+alpha
                if beta[i] > alphamax:
                    beta[i] = alphamax
            else:
                beta[i]=alpha
        Fxm=Fxm1
        if verbose:
            print "%d: |F(x)|=%.3f" %(n,norm(Fxm))

    return xm
