from numpy.testing import assert_equal, assert_almost_equal, assert_array_equal, \
        assert_array_almost_equal, assert_allclose, assert_, TestCase
from numpy import array, diff, shape, asarray, pi, sin, cos, arange, dot, \
     ravel, sqrt, inf, round
from scipy.interpolate.fitpack import splrep, splev, bisplrep, bisplev, \
     sproot, splprep, splint, spalde


def norm2(x):
    return sqrt(dot(x.T,x))

def f1(x,d=0):
    if d is None: return "sin"
    if x is None: return "sin(x)"
    if d%4 == 0: return sin(x)
    if d%4 == 1: return cos(x)
    if d%4 == 2: return -sin(x)
    if d%4 == 3: return -cos(x)

def f2(x,y=0,dx=0,dy=0):
    if x is None: return "sin(x+y)"
    d=dx+dy
    if d%4 == 0: return sin(x+y)
    if d%4 == 1: return cos(x+y)
    if d%4 == 2: return -sin(x+y)
    if d%4 == 3: return -cos(x+y)

def makepairs(x, y):
    """Helper function to create an array of pairs of x and y."""
    # Or itertools.product (>= python 2.6)
    xy = array([[a, b] for a in asarray(x) for b in asarray(y)])
    return xy.T

def put(*a):
    """Produce some output if file run directly"""
    import sys
    if hasattr(sys.modules['__main__'], '__put_prints'):
        sys.stderr.write("".join(map(str, a)) + "\n")

class TestSmokeTests(TestCase):
    """
    Smoke tests (with a few asserts) for fitpack routines -- mostly
    check that they are runnable
    """

    def check_1(self,f=f1,per=0,s=0,a=0,b=2*pi,N=20,at=0,xb=None,xe=None):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)    # nodes
        x1=a+(b-a)*arange(1,N,dtype=float)/float(N-1) # middle points of the nodes
        v,v1=f(x),f(x1)
        nk=[]

        def err_est(k, d):
            # Assume f has all derivatives < 1
            h = 1.0/float(N)
            tol = 5 * h**(.75*(k-d))
            if s > 0:
                tol += 1e5*s
            return tol

        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,xe=xe)
            if at:t=tck[0][k:-k]
            else: t=x1
            nd=[]
            for d in range(k+1):
                tol = err_est(k, d)
                err = norm2(f(t,d)-splev(t,tck,d)) / norm2(f(t,d))
                assert_(err < tol, (k, d, err, tol))
                nd.append((err, tol))
            nk.append(nd)
        put("\nf = %s  s=S_k(x;t,c)  x in [%s, %s] > [%s, %s]"%(f(None),
                                                        `round(xb,3)`,`round(xe,3)`,
                                                          `round(a,3)`,`round(b,3)`))
        if at:
            str="at knots"
        else:
            str="at the middle of nodes"
        put(" per=%d s=%s Evaluation %s"%(per,`s`,str))
        put(" k :  |f-s|^2  |f'-s'| |f''-.. |f'''-. |f''''- |f'''''")
        k=1
        for l in nk:
            put(' %d : '%k)
            for r in l:
                put(' %.1e  %.1e'%r)
            put('\n')
            k=k+1

    def check_2(self,f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)    # nodes
        v=f(x)

        def err_est(k, d):
            # Assume f has all derivatives < 1
            h = 1.0/float(N)
            tol = 5 * h**(.75*(k-d))
            if s > 0:
                tol += 1e5*s
            return tol

        nk=[]
        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,xe=xe)
            nk.append([splint(ia,ib,tck),spalde(dx,tck)])
        put("\nf = %s  s=S_k(x;t,c)  x in [%s, %s] > [%s, %s]"%(f(None),
                                                   `round(xb,3)`,`round(xe,3)`,
                                                    `round(a,3)`,`round(b,3)`))
        put(" per=%d s=%s N=%d [a, b] = [%s, %s]  dx=%s"%(per,`s`,N,`round(ia,3)`,`round(ib,3)`,`round(dx,3)`))
        put(" k :  int(s,[a,b]) Int.Error   Rel. error of s^(d)(dx) d = 0, .., k")
        k=1
        for r in nk:
            if r[0]<0: sr='-'
            else: sr=' '
            put(" %d   %s%.8f   %.1e "%(k,sr,abs(r[0]),
                                         abs(r[0]-(f(ib,-1)-f(ia,-1)))))
            d=0
            for dr in r[1]:
                err = abs(1-dr/f(dx,d))
                tol = err_est(k, d)
                assert_(err < tol, (k, d))
                put(" %.1e %.1e"%(err, tol))
                d=d+1
            put("\n")
            k=k+1

    def check_3(self,f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)    # nodes
        v=f(x)
        nk=[]
        put("  k  :     Roots of s(x) approx %s  x in [%s,%s]:"%\
              (f(None),`round(a,3)`,`round(b,3)`))
        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,xe=xe)
            roots = sproot(tck)
            if k == 3:
                assert_allclose(roots, pi*array([1, 2, 3, 4]),
                                rtol=1e-3)
            put('  %d  : %s'%(k,`roots.tolist()`))

    def check_4(self,f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)    # nodes
        x1=a+(b-a)*arange(1,N,dtype=float)/float(N-1) # middle points of the nodes
        v,v1=f(x),f(x1)
        nk=[]
        put(" u = %s   N = %d"%(`round(dx,3)`,N))
        put("  k  :  [x(u), %s(x(u))]  Error of splprep  Error of splrep "%(f(0,None)))
        for k in range(1,6):
            tckp,u=splprep([x,v],s=s,per=per,k=k,nest=-1)
            tck=splrep(x,v,s=s,per=per,k=k)
            uv=splev(dx,tckp)
            err1 = abs(uv[1]-f(uv[0]))
            err2 = abs(splev(uv[0],tck)-f(uv[0]))
            assert_(err1 < 1e-2)
            assert_(err2 < 1e-2)
            put("  %d  :  %s    %.1e           %.1e"%\
                  (k,`map(lambda x:round(x,3),uv)`,
                   err1,
                   err2))
        put("Derivatives of parametric cubic spline at u (first function):")
        k=3
        tckp,u=splprep([x,v],s=s,per=per,k=k,nest=-1)
        for d in range(1,k+1):
            uv=splev(dx,tckp,d)
            put(" %s "%(`uv[0]`))

    def check_5(self,f=f2,kx=3,ky=3,xb=0,xe=2*pi,yb=0,ye=2*pi,Nx=20,Ny=20,s=0):
        x=xb+(xe-xb)*arange(Nx+1,dtype=float)/float(Nx)
        y=yb+(ye-yb)*arange(Ny+1,dtype=float)/float(Ny)
        xy=makepairs(x,y)
        tck=bisplrep(xy[0],xy[1],f(xy[0],xy[1]),s=s,kx=kx,ky=ky)
        tt=[tck[0][kx:-kx],tck[1][ky:-ky]]
        t2=makepairs(tt[0],tt[1])
        v1=bisplev(tt[0],tt[1],tck)
        v2=f2(t2[0],t2[1])
        v2.shape=len(tt[0]),len(tt[1])
        err = norm2(ravel(v1-v2))
        assert_(err < 1e-2, err)
        put(err)

    def test_smoke_splrep_splev(self):
        put("***************** splrep/splev")
        self.check_1(s=1e-6)
        self.check_1()
        self.check_1(at=1)
        self.check_1(per=1)
        self.check_1(per=1,at=1)
        self.check_1(b=1.5*pi)
        self.check_1(b=1.5*pi,xe=2*pi,per=1,s=1e-1)

    def test_smoke_splint_spalde(self):
        put("***************** splint/spalde")
        self.check_2()
        self.check_2(per=1)
        self.check_2(ia=0.2*pi,ib=pi)
        self.check_2(ia=0.2*pi,ib=pi,N=50)

    def test_smoke_sproot(self):
        put("***************** sproot")
        self.check_3(a=0.1,b=15)

    def test_smoke_splprep_splrep_splev(self):
        put("***************** splprep/splrep/splev")
        self.check_4()
        self.check_4(N=50)

    def test_smoke_bisplrep_bisplev(self):
        put("***************** bisplev")
        self.check_5()

class TestSplev(TestCase):
    def test_1d_shape(self):
        x = [1,2,3,4,5]
        y = [4,5,6,7,8]
        tck = splrep(x, y)
        z = splev([1], tck)
        assert_equal(z.shape, (1,))
        z = splev(1, tck)
        assert_equal(z.shape, ())


if __name__ == "__main__":
    __put_prints = True
    import nose
    nose.runmodule()
