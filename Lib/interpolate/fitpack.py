#!/usr/bin/env python
"""
fitpack (dierckx in netlib) --- A Python-C wrapper to FITPACK (by P. Dierckx).
        FITPACK is a collection of FORTRAN programs for CURVE and SURFACE
        FITTING with SPLINES and TENSOR PRODUCT SPLINES. 

See
   http://www.cs.kuleuven.ac.be/cwis/research/nalag/research/topics/fitpack.html
or
   http://www.netlib.org/dierckx/index.html

Copyright 1999 Pearu Peterson all rights reserved,
Pearu Peterson <pearu@ioc.ee>          
Permission to use, modify, and distribute this software is given under the
terms of the LGPL.  See http://www.fsf.org

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Pearu Peterson

Running test programs:
    $ python fitpack.py 1 3    # run test programs 1, and 3
    $ python fitpack.py        # run all available test programs

TODO: Make interfaces to the following fitpack functions:
    For univariate splines: cocosp, concon, fourco, insert
    For bivariate splines: profil, regrid, parsur, surev
"""

__version__ = "$Revision$"[10:-1]
import _fitpack
from common_routines import *

_iermess = {0:["""\
    The spline has a residual sum of squares fp such that abs(fp-s)/s<=0.001""",None],
               -1:["""\
    The spline is an interpolating spline (fp=0)""",None],
               -2:["""\
    The spline is weighted least-squares polynomial of degree k.
    fp gives the upper bound fp0 for the smoothing factor s""",None],
               1:["""\
    The required storage space exceeds the available strorage space.
    Probably causes: nest to small or s is too small. (fp>s)""",ValueError],
               2:["""\
    A theoretically impossible results when finding a smoothin spline
    with fp = s. Probably causes: s too small. (abs(fp-s)/s>0.001)""",ValueError],
               3:["""\
    The maximal number of iterations (20) allowed for finding smoothing
    spline with fp=s has been reached. Probably causes: s too small.
    (abs(fp-s)/s>0.001)""",ValueError],
               10:["""\
    Error on input data""",ValueError],
               'unknown':["""\
    An error occured""",TypeError]}
_iermess2 = {0:["""\
    The spline has a residual sum of squares fp such that abs(fp-s)/s<=0.001""",None],
            -1:["""\
    The spline is an interpolating spline (fp=0)""",None],
            -2:["""\
    The spline is weighted least-squares polynomial of degree kx and ky.
    fp gives the upper bound fp0 for the smoothing factor s""",None],
            -3:["""\
    Warning. The coefficients of the spline have been computed as the minimal
    norm least-squares solution of a rank deficient system.""",None],
            1:["""\
    The required storage space exceeds the available strorage space.
    Probably causes: nxest or nyest to small or s is too small. (fp>s)""",ValueError],
            2:["""\
    A theoretically impossible results when finding a smoothin spline
    with fp = s. Probably causes: s too small or badly chosen eps.
    (abs(fp-s)/s>0.001)""",ValueError],
            3:["""\
    The maximal number of iterations (20) allowed for finding smoothing
    spline with fp=s has been reached. Probably causes: s too small.
    (abs(fp-s)/s>0.001)""",ValueError],
            4:["""\
    No more knots can be added because the number of B-spline coefficients
    already exceeds the number of data points m. Probably causes: either
    s or m too small. (fp>s)""",ValueError],
            5:["""\
    No more knots can be added because the additional knot would coincide
    with an old one. Probably cause: s too small or too large a weight
    to an inaccurate data point. (fp>s)""",ValueError],
            10:["""\
    Error on input data""",ValueError],
            11:["""\
    rwrk2 to small, i.e. there is not enough workspace for computing
    the minimal least-squares solution of a rank deficient system of linear
    equations.""",ValueError],
            'unknown':["""\
    An error occured""",TypeError]}
_parcur_cache = {'t': array([],'d'), 'wrk': array([],'d'), 'iwrk':array([],'i'),
                 'u': array([],'d'),'ub':0,'ue':1}
def _parcur(x,w=None,u=None,ub=None,ue=None,k=3,task=0,s=None,t=None,
           full_output=0,nest=None,per=0,quite=1):
    if task<=0:
        _parcur_cache = {'t': array([],'d'), 'wrk': array([],'d'),
                         'iwrk':array([],'i'),'u': array([],'d'),'ub':0,'ue':1}
    x=myasarray(x)
    idim,m=x.shape
    if per:
        for i in range(idim):
            if x[i][0]!=x[i][-1]:
                if quite<2:print 'Warning: Setting x[%d][%d]=x[%d][0]'%(i,m,i)
                x[i][-1]=x[i][0]
    if not 0<idim<11: raise TypeError,'0<idim<11 must hold'
    if w is None: w=ones(m,'d')
    else: w=myasarray(w)
    ipar=(u is not None)
    if ipar:
        _parcur_cache['u']=u
        if ub is None: _parcur_cache['ub']=u[0]
        else: _parcur_cache['ub']=ub
        if ue is None: _parcur_cache['ue']=u[-1]
        else: _parcur_cache['ue']=ue
    else: _parcur_cache['u']=zeros(m,'d')
    if not (1<=k<=5): raise TypeError, '1<=k=%d<=5 must hold'%(k)
    if not (-1<=task<=1): raise TypeError, 'task must be either -1,0, or 1'
    if (not len(w)==m) or (ipar==1 and (not len(u)==m)):
        raise TypeError,'Mismatch of input dimensions'
    if s is None: s=m-sqrt(2*m)
    if t is None and task==-1: raise TypeError, 'Knots must be given for task=-1'
    if t is not None:
        _parcur_cache['t']=myasarray(t)
    n=len(_parcur_cache['t'])
    if task==-1 and n<2*k+2:
        raise TypeError, 'There must be at least 2*k+2 knots for task=-1'
    if m<=k: raise TypeError, 'm>k must hold'
    if nest is None: nest=m/2

    if (task>=0 and s==0) or (nest<0):
        if per: nest=m+2*k
        else: nest=m+k+1
    nest=max(nest,2*k+3)
    u=_parcur_cache['u']
    ub=_parcur_cache['ub']
    ue=_parcur_cache['ue']
    t=_parcur_cache['t']
    wrk=_parcur_cache['wrk']
    iwrk=_parcur_cache['iwrk']
    t,c,o=_fitpack._parcur(ravel(transpose(x)),w,u,ub,ue,k,task,ipar,s,t,
                             nest,wrk,iwrk,per)
    _parcur_cache['u']=o['u']
    _parcur_cache['ub']=o['ub']
    _parcur_cache['ue']=o['ue']
    _parcur_cache['t']=t
    _parcur_cache['wrk']=o['wrk']
    _parcur_cache['iwrk']=o['iwrk']
    ier,fp,n=o['ier'],o['fp'],len(t)
    u=o['u']
    c.shape=idim,n-k-1
    tcku = [t,list(c),k],u
    if ier<=0 and not quite:
        print _iermess[ier][0]
        print "\tk=%d n=%d m=%d fp=%f s=%f"%(k,len(t),m,fp,s)
    if ier>0 and not full_output:
        if ier in [1,2,3]:
            print "Warning: "+_iermess[ier][0]
        else:
            try:
                raise _iermess[ier][1],_iermess[ier][0]
            except KeyError:
                raise _iermess['unknown'][1],_iermess['unknown'][0]
    if full_output:
        try:
            return tcku,fp,ier,_iermess[ier][0]
        except KeyError:
            return tcku,fp,ier,_iermess['unknown'][0]
    else:
        return tcku
_curfit_cache = {'t': array([],'d'), 'wrk': array([],'d'), 'iwrk':array([],'i')}
def _curfit(x,y,w=None,xb=None,xe=None,k=3,task=0,s=None,t=None,
           full_output=0,nest=None,per=0,quite=1):
    """
    splrep --- See documentation.
    """
    if task<=0:
        _curfit_cache = {'t': array([],'d'), 'wrk': array([],'d'),
                         'iwrk':array([],'i')}
    x,y=map(myasarray,[x,y])
    m=len(x)
    if w is None: w=ones(m,'d')
    else: w=myasarray(w)
    if not len(w) == m: raise TypeError,' len(w)=%d is not equal to m=%d'%(len(w),m)
    if xb is None: xb=x[0]
    if xe is None: xe=x[-1]
    if not (-1<=task<=1): raise TypeError, 'task must be either -1,0, or 1'
    if s is None: s=m-sqrt(2*m)
    if t is None and task==-1: raise TypeError, 'Knots must be given for task=-1'
    if t is not None: _curfit_cache['t']=myasarray(t)
    n=len(_curfit_cache['t'])
    if task==-1 and n<2*k+2:
        raise TypeError, 'There must be at least 2*k+2 knots for task=-1'
    if (m != len(y)) or (m != len(w)):
        raise TypeError, 'Lengths of the first three arguments (x,y,w) must be equal'
    if not (1<=k<=5):
        raise TypeError, 'Given degree of the spline (k=%d) is not supported. (1<=k<=5)'%(k)
    if m<=k: raise TypeError, 'm>k must hold'
    if nest is None: nest=m/2
    if nest<0:
        if per: nest=m+2*k
        else: nest=m+k+1
    nest=max(nest,2*k+3)
    if task>=0 and s==0:
        if per: nest=m+2*k
        else: nest=m+k+1
    if task==-1:
        _curfit_cache['t']=myasarray(t)
        if not (2*k+2<=len(t)<=min(nest,m+k+1)):
            raise TypeError, 'Number of knots n is not acceptable (2*k+2<=n<=min(nest,m+l+1))'
    t=_curfit_cache['t']
    wrk=_curfit_cache['wrk']
    iwrk=_curfit_cache['iwrk']
    t,c,o = _fitpack._curfit(x,y,w,xb,xe,k,task,s,t,nest,wrk,iwrk,per)
    _curfit_cache['t']=t
    _curfit_cache['wrk']=o['wrk']
    _curfit_cache['iwrk']=o['iwrk']
    ier,fp=o['ier'],o['fp']
    tck = [t,c,k]
    if ier<=0 and not quite:
        print _iermess[ier][0]
        print "\tk=%d n=%d m=%d fp=%f s=%f"%(k,len(t),m,fp,s)
    if ier>0 and not full_output:
        if ier in [1,2,3]:
            print "Warning: "+_iermess[ier][0]
        else:
            try:
                raise _iermess[ier][1],_iermess[ier][0]
            except KeyError:
                raise _iermess['unknown'][1],_iermess['unknown'][0]
    if full_output:
        try:
            return tck,fp,ier,_iermess[ier][0]
        except KeyError:
            return tck,fp,ier,_iermess['unknown'][0]
    else:
        return tck
splrep=_curfit
splprep=_parcur
def _ntlist(l): # return non-trivial list
    return l
    #if len(l)>1: return l
    #return l[0]
def splev(x,tck,der=0):
    """
    splev --- see documentation.
    """
    t,c,k=tck
    try:
        c[0][0]
        return map(lambda c,x=x,t=t,k=k,der=der:splev(x,[t,c,k],der),c)
    except: pass
    if not (0<=der<=k):
        raise ValueError,"0<=der=%d<=k=%d must hold"%(der,k)
    x=myasarray(x)
    y,ier=_fitpack._spl_(x,der,t,c,k)
    if ier==10: raise ValueError,"Invalid input data"
    if ier: raise TypeError,"An error occurred"
    if len(y)>1: return y
    return y[0]
def splint(a,b,tck,full_output=0):
    t,c,k=tck
    try: c[0][0];return _ntlist(map(lambda c,a=a,b=b,t=t,k=k:splint(a,b,[t,c,k]),c))
    except: pass
    aint,wrk=_fitpack._splint(t,c,k,a,b)
    if full_output: return aint,wrk
    else: return aint
def sproot(tck,mest=10):
    t,c,k=tck
    if k==4: t=t[1:-1]
    if k==5: t=t[2:-2]
    try: c[0][0];return _ntlist(map(lambda c,t=t,k=k,mest=mest:sproot([t,c,k],mest),c))
    except: pass
    if len(t)<8:
        raise TypeError,"The number of knots %d>=8"%(len(t))
    z,ier=_fitpack._sproot(t,c,k,mest)
    if ier==10:
        raise TypeError,"Invalid input data. t1<=..<=t4<t5<..<tn-3<=..<=tn must hold."
    if ier==0: return z
    if ier==1:
        print "Warning: the number of zeros exceeds mest"
        return z
    raise TypeError,"Unknown error"
def spalde(x,tck):
    t,c,k=tck
    try:
        c[0][0]
        return _ntlist(map(lambda c,x=x,t=t,k=k:spalde(x,[t,c,k]),c))
    except: pass
    try: x=x.tolist()
    except:
        try: x=list(x)
        except: x=[x]
    if len(x)>1:
        return map(lambda x,tck=tck:spalde(x,tck),x)
    d,ier=_fitpack._spalde(t,c,k,x[0])
    if ier==0: return d
    if ier==10:
        raise TypeError,"Invalid input data. t(k)<=x<=t(n-k+1) must hold."
    raise TypeError,"Unknown error"
#def _curfit(x,y,w=None,xb=None,xe=None,k=3,task=0,s=None,t=None,
#           full_output=0,nest=None,per=0,quite=1):
_surfit_cache = {'tx': array([],'d'),'ty': array([],'d'),
                 'wrk': array([],'d'), 'iwrk':array([],'i')}
def _surfit(x,y,z,w=None,xb=None,xe=None,yb=None,ye=None,kx=3,ky=3,task=0,s=None,
           eps=1e-16,tx=None,ty=None,full_output=0,nxest=None,nyest=None,quite=1):
    x,y,z=map(myasarray,[x,y,z])
    m=len(x)
    if not (m==len(y)==len(z)): raise TypeError, 'len(x)==len(y)==len(z) must hold.'
    if w is None: w=ones(m,'d')
    else: w=myasarray(w)
    if not len(w) == m: raise TypeError,' len(w)=%d is not equal to m=%d'%(len(w),m)
    if xb is None: xb=x[0]
    if xe is None: xe=x[-1]
    if yb is None: yb=y[0]
    if ye is None: ye=y[-1]
    if not (-1<=task<=1): raise TypeError, 'task must be either -1,0, or 1'
    if s is None: s=m-sqrt(2*m)
    if tx is None and task==-1: raise TypeError, 'Knots_x must be given for task=-1'
    if tx is not None: _curfit_cache['tx']=myasarray(tx)
    nx=len(_surfit_cache['tx'])
    if ty is None and task==-1: raise TypeError, 'Knots_y must be given for task=-1'
    if ty is not None: _curfit_cache['ty']=myasarray(ty)
    ny=len(_surfit_cache['ty'])
    if task==-1 and nx<2*kx+2:
        raise TypeError, 'There must be at least 2*kx+2 knots_x for task=-1'
    if task==-1 and ny<2*ky+2:
        raise TypeError, 'There must be at least 2*ky+2 knots_x for task=-1'
    if not ((1<=kx<=5) and (1<=ky<=5)): 
        raise TypeError, 'Given degree of the spline (kx,ky=%d,%d) is not supported. (1<=k<=5)'%(kx,ky)
    if m<=(kx+1)*(ky+1): raise TypeError, 'm>(kx+1)(ky+1) must hold'
    if nxest is None: nxest=kx+sqrt(m/2)
    if nyest is None: nyest=ky+sqrt(m/2)
    nxest,nyest=max(nxest,2*kx+3),max(nyest,2*ky+3)
    if task>=0 and s==0:
        nxest=int(kx+sqrt(3*m))
        nyest=int(ky+sqrt(3*m))
    if task==-1:
        _surfit_cache['tx']=myasarray(tx)
        _surfit_cache['ty']=myasarray(ty)
    tx,ty=_surfit_cache['tx'],_surfit_cache['ty']
    wrk=_surfit_cache['wrk']
    iwrk=_surfit_cache['iwrk']
    u,v,km,ne=nxest-kx-1,nyest-ky-1,max(kx,ky)+1,max(nxest,nyest)
    bx,by=kx*v+ky+1,ky*u+kx+1
    b1,b2=bx,bx+v-ky
    if bx>by: b1,b2=by,by+u-kx
    lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    lwrk2=1
    tx,ty,c,o = _fitpack._surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,task,s,eps,
                                   tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)
    _curfit_cache['tx']=tx
    _curfit_cache['ty']=ty
    _curfit_cache['wrk']=o['wrk']
    ier,fp=o['ier'],o['fp']
    tck=[tx,ty,c,kx,ky]
    if ier<=0 and not quite:
        print _iermess2[ier][0]
        print "\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s)
    ierm=min(11,max(-3,ier))
    if ierm>0 and not full_output:
        if ier in [1,2,3,4,5]:
            print "Warning: "+_iermess2[ierm][0]
            print "\tkx,ky=%d,%d nx,ny=%d,%d m=%d fp=%f s=%f"%(kx,ky,len(tx),
                                                           len(ty),m,fp,s)
        else:
            try:
                raise _iermess2[ierm][1],_iermess2[ierm][0]
            except KeyError:
                raise _iermess2['unknown'][1],_iermess2['unknown'][0]
    if full_output:
        try:
            return tck,fp,ier,_iermess2[ierm][0]
        except KeyError:
            return tck,fp,ier,_iermess2['unknown'][0]
    else:
        return tck
bisplrep=_surfit

def bisplev(x,y,tck,dx=0,dy=0):
    tx,ty,c,kx,ky=tck
    if not (0<=dx<kx): raise ValueError,"0<=dx=%d<kx=%d must hold"%(dx,kx)
    if not (0<=dy<ky): raise ValueError,"0<=dy=%d<ky=%d must hold"%(dy,ky)
    x,y=map(myasarray,[x,y])
    z,ier=_fitpack._bispev(tx,ty,c,kx,ky,x,y,dx,dy)
    if ier==10: raise ValueError,"Invalid input data"
    if ier: raise TypeError,"An error occurred"
    z.shape=len(x),len(y)
    if len(z)>1: return z
    if len(z[0])>1: return z[0]
    return z[0][0]


if __name__ == "__main__":
    import sys,string
    runtest=range(10)
    if len(sys.argv[1:])>0:
        runtest=map(string.atoi,sys.argv[1:])
    put=sys.stdout.write
    def norm2(x):
        return matrixmultiply(transpose(x),x)
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
    def test1(f=f1,per=0,s=0,a=0,b=2*pi,N=20,at=0,xb=None,xe=None):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,typecode='d')/float(N)    # nodes
        x1=a+(b-a)*arange(1,N,typecode='d')/float(N-1) # middle points of the nodes
        v,v1=f(x),f(x1)
        nk=[]
        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,nest=-1,xe=xe)
            if at:t=tck[0][k:-k]
            else: t=x1
            nd=[]
            for d in range(k+1):
                nd.append(norm2(f(t,d)-splev(t,tck,d)))
            nk.append(nd)
        print "\nf = %s  s=S_k(x;t,c)  x in [%s, %s] > [%s, %s]"%(f(None),
                                                        `round(xb,3)`,`round(xe,3)`,
                                                          `round(a,3)`,`round(b,3)`)
        if at: str="at knots"
        else: str="at the middle of nodes"
        print " per=%d s=%s Evaluation %s"%(per,`s`,str)
        print " k :  |f-s|^2  |f'-s'| |f''-.. |f'''-. |f''''- |f'''''"
        k=1
        for l in nk:
            put(' %d : '%k)
            for r in l:
                put(' %.1e'%r)
            put('\n')
            k=k+1
    def test2(f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,typecode='d')/float(N)    # nodes
        v=f(x)
        nk=[]
        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,nest=-1,xe=xe)
            nk.append([splint(ia,ib,tck),spalde(dx,tck)])
        print "\nf = %s  s=S_k(x;t,c)  x in [%s, %s] > [%s, %s]"%(f(None),
                                                   `round(xb,3)`,`round(xe,3)`,
                                                    `round(a,3)`,`round(b,3)`)
        print " per=%d s=%s N=%d [a, b] = [%s, %s]  dx=%s"%(per,`s`,N,`round(ia,3)`,`round(ib,3)`,`round(dx,3)`)
        print " k :  int(s,[a,b]) Int.Error   Rel. error of s^(d)(dx) d = 0, .., k"
        k=1
        for r in nk:
            if r[0]<0: sr='-'
            else: sr=' '
            put(" %d   %s%.8f   %.1e "%(k,sr,abs(r[0]),
                                         abs(r[0]-(f(ib,-1)-f(ia,-1)))))
            d=0
            for dr in r[1]:
                put(" %.1e "%(abs(1-dr/f(dx,d))))
                d=d+1
            put("\n")
            k=k+1
    def test3(f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,typecode='d')/float(N)    # nodes
        v=f(x)
        nk=[]
        print "  k  :     Roots of s(x) approx %s  x in [%s,%s]:"%\
              (f(None),`round(a,3)`,`round(b,3)`)
        for k in range(1,6):
            tck=splrep(x,v,s=s,per=per,k=k,nest=-1,xe=xe)
            print '  %d  : %s'%(k,`sproot(tck).tolist()`)
    def test4(f=f1,per=0,s=0,a=0,b=2*pi,N=20,xb=None,xe=None,
              ia=0,ib=2*pi,dx=0.2*pi):
        if xb is None: xb=a
        if xe is None: xe=b
        x=a+(b-a)*arange(N+1,typecode='d')/float(N)    # nodes
        x1=a+(b-a)*arange(1,N,typecode='d')/float(N-1) # middle points of the nodes
        v,v1=f(x),f(x1)
        nk=[]
        print " u = %s   N = %d"%(`round(dx,3)`,N)
        print "  k  :  [x(u), %s(x(u))]  Error of splprep  Error of splrep "%(f(0,None))
        for k in range(1,6):
            tckp,u=splprep([x,v],s=s,per=per,k=k,nest=-1)
            tck=splrep(x,v,s=s,per=per,k=k,nest=-1)
            uv=splev(dx,tckp)
            print "  %d  :  %s    %.1e           %.1e"%\
                  (k,`map(lambda x:round(x,3),uv)`,
                   abs(uv[1]-f(uv[0])),
                   abs(splev(uv[0],tck)-f(uv[0])))
        print "Derivatives of parametric cubic spline at u (first function):"
        k=3
        tckp,u=splprep([x,v],s=s,per=per,k=k,nest=-1)
        for d in range(1,k+1):
            uv=splev(dx,tckp,d)
            put(" %s "%(`uv[0]`))
        print
    def makepairs(x,y):
        x,y=map(myasarray,[x,y])
        xy=array(map(lambda x,y:map(None,len(y)*[x],y),x,len(x)*[y]))
        sh=xy.shape
        xy.shape=sh[0]*sh[1],sh[2]
        return transpose(xy)
    def test5(f=f2,kx=3,ky=3,xb=0,xe=2*pi,yb=0,ye=2*pi,Nx=20,Ny=20,s=0):
        x=xb+(xe-xb)*arange(Nx+1,typecode='d')/float(Nx)
        y=yb+(ye-yb)*arange(Ny+1,typecode='d')/float(Ny)
        xy=makepairs(x,y)
        tck=bisplrep(xy[0],xy[1],f(xy[0],xy[1]),s=s,kx=kx,ky=ky)
        tt=[tck[0][kx:-kx],tck[1][ky:-ky]]
        t2=makepairs(tt[0],tt[1])
        v1=bisplev(tt[0],tt[1],tck)
        v2=f2(t2[0],t2[1])
        v2.shape=len(tt[0]),len(tt[1])
        print norm2(ravel(v1-v2))
    if 1 in runtest:
        print """\
******************************************
\tTests of splrep and splev
******************************************"""
        test1(s=1e-6)
        test1()
        test1(at=1)
        test1(per=1)
        test1(per=1,at=1)
        test1(b=1.5*pi)
        test1(b=1.5*pi,xe=2*pi,per=1,s=1e-1)
    if 2 in runtest:
        print """\
******************************************
\tTests of splint and spalde
******************************************"""
        test2()
        test2(per=1)
        test2(ia=0.2*pi,ib=pi)
        test2(ia=0.2*pi,ib=pi,N=50)
    if 3 in runtest:
        print """\
******************************************
\tTests of sproot
******************************************"""
        test3(a=0,b=15)
        print "Note that if k is not 3, some roots are missed or incorrect"
    if 4 in runtest:
        print """\
******************************************
\tTests of splprep, splrep, and splev
******************************************"""
        test4()
        test4(N=50)
    if 5 in runtest:
        print """\
******************************************
\tTests of bisplrep, bisplev
******************************************"""
        test5()







