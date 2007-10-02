#!/usr/bin/env python
# Created by John Travers, January 2007
""" Test functions for spline.fitpack module
"""
__usage__ = """
Build spline:
  python setup_spline.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.spline.test(<level>)'
Run tests if spline is not installed:
  python tests/test_fitpack.py [<level>]
"""

import sys
from numpy.testing import *
from numpy import array, arange, around, pi, sin, ravel, zeros, asarray

set_package_path()
from spline.fitpack import splprep, splrep, splev, sproot, splint, spalde
from spline.fitpack import bisplev, bisplrep, splprep
restore_path()

set_local_path()
from dierckx_test_data import *
restore_path()

class TestSplrepSplev(NumpyTestCase):
    def check_curfit_against_dierckx_smth(self):
        x,y = curfit_test['x'],curfit_test['y']
        k,s = curfit_test_smth['k'],curfit_test_smth['s']
        iopt = curfit_test_smth['iopt']
        for i in range(len(k)):
            tck = splrep(x,y,k=k[i],task=iopt[i],s=s[i])
            out = splrep(x,y,k=k[i],task=iopt[i],s=s[i],full_output=1)
            tck,fp=out[0],out[1]
            sp = splev(x,tck)
            assert_almost_equal(fp, curfit_test_smth['fp'][i], decimal=2)
            assert_array_almost_equal(around(tck[0],1),
                                      curfit_test_smth['t'][i])
            assert_array_almost_equal(around(tck[1],4),
                                      curfit_test_smth['c'][i], decimal=3)
            assert_array_almost_equal(around(sp,1),
                                      curfit_test_smth['sp'][i])

    def check_curfit_against_dierckx_lsq(self):
        """ Test against results obtined from the pure fortran routines.
            
            Here we check simple spline creation and evaluation.
        """
        x,y = curfit_test['x'],curfit_test['y']
        k = curfit_test_lsq['k']
        for i in range(len(k)):
            t = curfit_test_lsq['t'][i]
            out = splrep(x,y,t=t,k=k[i],full_output=1)
            tck,fp=out[0],out[1]
            sp = splev(x,tck)
            assert_almost_equal(fp, curfit_test_lsq['fp'][i], decimal=2)
            assert_array_almost_equal(around(tck[1],4),
                                      curfit_test_lsq['c'][i], decimal=3)
            assert_array_almost_equal(around(sp,1),
                                      curfit_test_lsq['sp'][i])

    def check_percur_against_dierckx(self):
            x,y = percur_test['x'], percur_test['y']
            k,s = percur_test['k'], percur_test['s']
            iopt, res = percur_test['iopt'], percur_test['res']
            err = percur_test['err']
            coef, knots = percur_test['coef'], percur_test['knots']
            sp = percur_test['sp']
            for i in range(len(k)):
                if iopt[i] != -1:
                    tck,fp,ier,msg = splrep(x,y,k=k[i],task=iopt[i],s=s[i],
                                                    per=True,full_output=True)
                else:
                    tck,fp,ier,msg = splrep(x,y,t=knots[i],k=k[i],task=iopt[i],
                                                     per=True,full_output=True)
                tt,cc,kk = tck
                tt,cc = asarray(tt), asarray(cc)
                assert_almost_equal(ier,err[i])
                assert_almost_equal(fp,res[i],decimal=1)
                assert_array_almost_equal(tt,knots[i], decimal=3)
                assert_array_almost_equal(cc,coef[i], decimal=3)
                yy = asarray(splev(x,tck))
                assert_array_almost_equal(yy,sp[i], decimal=3)
            
class TestSplprepSplev(NumpyTestCase):
    def check_parcur_against_dierckx(self):
        xa,xo = parcur_test['xa'], parcur_test['xo']
        k,s = parcur_test['k'], parcur_test['s']
        u = parcur_test['u']
        ub,ue = parcur_test['ub'], parcur_test['ue']
        iopt, res = parcur_test['iopt'], parcur_test['res']
        err, ipar = parcur_test['err'], parcur_test['ipar']
        knots = parcur_test['knots']
        sx, sy = parcur_test['sx'], parcur_test['sy']
        sp = parcur_test['sp']
        x = array([xa, xo])
        for i in range(len(k)):
            if iopt[i] != -1:
                if ipar[i] == 1:
                    tcku,fp,ier,msg = splprep(x,u=u,ub=ub,ue=ue,k=k[i],
                                   task=iopt[i],s=s[i],full_output=True)
                else:
                    tcku,fp,ier,msg = splprep(x,ub=ub,ue=ue,k=k[i],
                                   task=iopt[i],s=s[i],full_output=True)
            else:
                tcku,fp,ier,msg = splprep(x,ub=ub,ue=ue,t=knots[i],
                                   k=k[i],task=iopt[i],full_output=True)
            tck,u = tcku
            tt,cc,kk = tck
            tt,cc = asarray(tt), asarray(cc)
            assert_almost_equal(ier,err[i])
            assert_almost_equal(fp,res[i],decimal=3)
            assert_array_almost_equal(tt,knots[i], decimal=3)
            assert_array_almost_equal(cc[0],sx[i], decimal=3)
            assert_array_almost_equal(cc[1],sy[i], decimal=3)
            y = asarray(splev(u,tck))
            yy = zeros(64, 'float')
            yy[0:-1:2] = y[0]
            yy[1::2] = y[1]
            assert_array_almost_equal(yy,sp[i], decimal=3)
    
    def check_clocur_against_dierckx(self):
        xa,xo = clocur_test['xa'], clocur_test['xo']
        k,s = clocur_test['k'], clocur_test['s']
        u = clocur_test['u']
        iopt, res = clocur_test['iopt'], clocur_test['res']
        err, ipar = clocur_test['err'], clocur_test['ipar']
        knots = clocur_test['knots']
        sx, sy = clocur_test['sx'], clocur_test['sy']
        sp = clocur_test['sp']
        x = array([xa, xo]) 
        for i in range(len(k)):
            if iopt[i] != -1:
                if ipar[i] == 1:
                    tcku,fp,ier,msg = splprep(x,u=u,k=k[i],task=iopt[i],
                                        s=s[i],per=True,full_output=True)
                else:
                    tcku,fp,ier,msg = splprep(x,k=k[i],task=iopt[i],
                                        s=s[i],per=True,full_output=True)
            else:
                tcku,fp,ier,msg = splprep(x,t=knots[i],k=k[i],task=iopt[i],
                                                 per=True,full_output=True)
            tck,u = tcku
            tt,cc,kk = tck
            tt,cc = asarray(tt), asarray(cc)
            assert_almost_equal(ier,err[i])
            assert_almost_equal(fp,res[i],decimal=3)
            assert_array_almost_equal(tt,knots[i], decimal=3)
            assert_array_almost_equal(cc[0],sx[i], decimal=3)
            assert_array_almost_equal(cc[1],sy[i], decimal=3)
            y = asarray(splev(u,tck))
            yy = zeros(36, 'float')
            yy[0:-1:2] = y[0,:-1]
            yy[1::2] = y[1,:-1]
            assert_array_almost_equal(yy,sp[i], decimal=3)

class TestSplintSpalde(NumpyTestCase):
    def check_splint_spalde(self):
        per = [0, 1, 0]
        N = [20, 20, 50]
        ia = [0, 0, 0.2*pi]
        ib = [0, 0, pi]
        a,b = 0,2*pi
        dx = 0.2*pi
        k = range(1,6)
        for i in range(len(per)):
            x=a+(b-a)*arange(N[i]+1,dtype=float)/float(N[i])
            v=f1(x)
            for j in range(len(k)):
                tck = splrep(x,v,k=k[j],s=0,per=per[i])
                ir = splint(ia[i],ib[i],tck)
                dr = spalde(dx,tck)
                assert_almost_equal(ir, f1(ib[i],-1)-f1(ia[i],-1), decimal=2)
                d=0
                for ddr in dr:
                    if d<k[j]-1:
                        assert_almost_equal(1, ddr/f1(dx,d), decimal=2)
                    d=d+1

class TestSplder(NumpyTestCase):
    def check_splder(self):
        N = 50
        a,b = 0,2*pi
        dx = 0.2*pi
        k = range(3,6)
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)
        v=f1(x)
        for j in range(len(k)):
            tck = splrep(x,v,k=k[j],s=0)
            dr = spalde(dx,tck)
            d2 = splev(dx,tck,der=1)
            assert_almost_equal(1, dr[1]/f1(dx,1.0), decimal=2)

class TestSproot(NumpyTestCase):
    def check_sproot(self):
        a=0
        b=15
        N=20
        x=a+(b-a)*arange(N+1,dtype=float)/float(N)
        v=f1(x)
        k=3
        tck = splrep(x,v,k=k,s=0)
        ex = array([0.0, pi, 2.0*pi, 3.0*pi, 4.0*pi])
        assert_array_almost_equal(sproot(tck),ex, decimal=3)

class TestBisplevBisplrep(NumpyTestCase):
    def test_bisplev_bisplrep(self):
        f=f2; kx=3; ky=3; xb=0; xe=2*pi
        yb=0; ye=2*pi; Nx=20; Ny=20; s=0
        x=xb+(xe-xb)*arange(Nx+1,dtype=float)/float(Nx)
        y=yb+(ye-yb)*arange(Ny+1,dtype=float)/float(Ny)
        xy=makepairs(x,y)
        tck=bisplrep(xy[0],xy[1],f(xy[0],xy[1]),s=s,kx=kx,ky=ky)
        tt=[tck[0][kx:-kx],tck[1][ky:-ky]]
        t2=makepairs(tt[0],tt[1])
        v1=bisplev(tt[0],tt[1],tck)
        v2=f2(t2[0],t2[1])
        v2.shape=len(tt[0]),len(tt[1])
        assert_almost_equal(0.0, norm2(ravel(v1-v2)), decimal=5)

class TestParcur(NumpyTestCase):
    def check_parcur(self):
        f=f1; per=0; s=0; a=0; b=2*pi;
        N=[20,50]
        ia=0; ib=2*pi; dx=0.2*pi
        xb,xe=a,b
        for i in range(len(N)):
            x=a+(b-a)*arange(N[i]+1,dtype=float)/float(N[i])
            x1=a+(b-a)*arange(1,N[i],dtype=float)/float(N[i]-1) # mid points
            v,v1=f(x),f(x1)
            for k in range(1,6):
                tckp,u=splprep([x,v],s=s,per=per,k=k)
                tck=splrep(x,v,s=s,per=per,k=k)
                uv=splev(dx,tckp)
                assert_almost_equal(0.0, around(abs(uv[1]-f(uv[0])),2), 
                                                                    decimal=1)
                assert_almost_equal(0.0, 
                        around(abs(splev(uv[0],tck)-f(uv[0])),2),decimal=1)

if __name__ == "__main__":
    NumpyTest().run()
