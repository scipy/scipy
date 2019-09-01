from __future__ import division, print_function, absolute_import

from iterative import *
import numpy as np


def test_fun(alpha, x, beta, y, A, n):
    # compute z = alpha*A*x + beta*y
    xx = x[:n]
    yy = y[:n]
    w = np.dot(A,xx)
    z = alpha*w+beta*yy
    y[:n] = z
    return


def test_fun_t(alpha, x, beta, y, A, n):
    # compute z = alpha*A*x + beta*y
    xx = x[:n]
    yy = y[:n]
    AA = np.conj(np.transpose(A))
    w = np.dot(AA,xx)
    z = alpha*w+beta*yy
    y[:n] = z
    return


def test_psolve(x,b,n):
    x[:n] = b[:n]
    return


def test_psolve_t(x,b,n):
    x[:n] = b[:n]
    return


def test_psolveq(x,b,which,n):
    x[:n] = b[:n]
    return


def test_psolveq_t(x,b,which,n):
    x[:n] = b[:n]
    return


n = 5
dA = 1.0*np.array([[2, -1, 0, 0, 0],
                   [-1, 2, -1, 0, 0],
                   [0, -1, 2, -1, 0],
                   [0, 0, -1, 2, -1],
                   [0, 0, 0, 1, 2]])
db = 1.0*np.array([0,1,1,0,0])

##zA = (1.0+0j)*np.array([[      2, -1+0.1j,       0,       0,       0],
##                        [-1+0.1j,       2, -1-0.1j,       0,       0],
##                        [      0, -1-0.1j,       2, -1+0.1j,       0],
##                        [      0,       0, -1+0.1j,       2, -1-0.1j],
##                        [      0,       0,       0,      -1,  2-0.1j]])
zA = (1.0+0j)*np.array([[2, -1 + 1j, 0, 0, 0],
                        [-1+0.1j, 2, -1-0.1j, 0, 0],
                        [0, -1 - 1j, 2, -1+0.1j, 0],
                        [0, 0, -1+0.1j, 2, -1-0.1j],
                        [0, 0, 0, -1, 2-0.1j]])
zb = (1.0+0j)*np.array([0,1,1,0,0])

dx = 0*db.copy()
zx = 0*zb.copy()

diter = 1000
dresid = 1e-6
ziter = 1000
zresid = 1e-6
drestrt = n
zrestrt = n

############### BiCG #######################
dx,diter,dresid,dinfor = dbicg(db,dx,diter,dresid,test_fun,test_fun_t,test_psolve,test_psolve_t,(dA,n),(dA,n),(n,),(n,))

zx,ziter,zresid,zinfor = zbicg(zb,zx,ziter,zresid,test_fun,test_fun_t,test_psolve,test_psolve_t,(zA,n),(zA,n),(n,),(n,))

############### BiCGSTAB ###################
#dx,diter,dresid,dinfor = dbicgstab(db,dx,diter,dresid,test_fun,test_psolve,(dA,n),(n,))
#zx,ziter,zresid,zinfor = zbicgstab(zb,zx,ziter,zresid,test_fun,test_psolve,(zA,n),(n,))

############### CG #########################
##dA = 1.0*array([[ 2,  -1,   0,   0,   0],
##                [-1,   2,  -1,   0,   0],
##                [ 0,  -1,   2,  -1,   0],
##                [ 0,   0,  -1,   2,  -1],
##                [ 0,   0,   0,  -1,   2]])
##dx = db.copy()
##zA = (1.0+0j)*array([[      2, -1+0.1j,       0,       0,       0],
##                     [-1+0.1j,       2, -1-0.1j,       0,       0],
##                     [      0, -1-0.1j,       2, -1+0.1j,       0],
##                     [      0,       0, -1+0.1j,       2, -1-0.1j],
##                     [      0,       0,       0,      -1,  2-0.1j]])
##zx = zb.copy()
##dx,diter,dresid,dinfor = dcg(db,dx,diter,dresid,test_fun,test_psolve,(dA,n),(n,))
##zx,ziter,zresid,zinfor = zcg(zb,zx,ziter,zresid,test_fun,test_psolve,(zA,n),(n,))

############### CGS ########################
#dx,diter,dresid,dinfor = dcgs(db,dx,diter,dresid,test_fun,test_psolve,(dA,n),(n,))
#zx,ziter,zresid,zinfor = zcgs(zb,zx,ziter,zresid,test_fun,test_psolve,(zA,n),(n,))

############### GMRES ######################
#dx,diter,dresid,dinfor = dgmres(db,dx,drestrt,diter,dresid,test_fun,test_psolve,(dA,n),(n,))
#zx,ziter,zresid,zinfor = zgmres(zb,zx,zrestrt,ziter,zresid,test_fun,test_psolve,(zA,n),(n,))

############### QMR ########################
#dx,diter,dresid,dinfor = dqmr(db,dx,diter,dresid,test_fun,test_fun_t,test_psolveq,test_psolveq_t,(dA,n),(dA,n),(n,),(n,))

#zx,ziter,zresid,zinfor = zqmr(zb,zx,ziter,zresid,test_fun,test_fun_t,test_psolveq,test_psolveq_t,(zA,n),(zA,n),(n,),(n,))

print()
print('**************** double *****************')
print('iter:',diter, 'resid:', dresid, 'info:',dinfor)
print('x=',dx)
print('*****************************************')
print()
print()
print('**************** complex ****************')
print('iter:',ziter, 'resid:',zresid, 'info:',zinfor)
print('x=',zx)
print('*****************************************')
print()
