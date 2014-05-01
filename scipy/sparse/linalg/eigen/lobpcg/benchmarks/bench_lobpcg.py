from __future__ import division, print_function, absolute_import

from scipy import *
from scipy.sparse.linalg import lobpcg
from symeig import symeig
from pylab import plot, show, legend, xlabel, ylabel
set_printoptions(precision=3,linewidth=90)
import time


def test(n):
    x = arange(1,n+1)
    B = diag(1./x)
    y = arange(n-1,0,-1)
    z = arange(2*n-1,0,-2)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A,B


def as2d(ar):
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = nm.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


def precond(x):
    y = linalg.cho_solve((LorU, lower),x)
    return as2d(y)

m = 10  # Blocksize
N = array(([128,256,512,1024,2048]))  # Increasing matrix size

data1 = []
data2 = []

for n in N:
    print('******', n)
    A,B = test(n)  # Mikota pair
    X = rand(n,m)
    X = linalg.orth(X)

    tt = time.clock()
    (LorU, lower) = linalg.cho_factor(A, lower=0, overwrite_a=0)
    eigs,vecs = lobpcg.lobpcg(X,A,B,operatorT=precond,
                              residualTolerance=1e-4, maxIterations=40)
    data1.append(time.clock()-tt)
    eigs = sort(eigs)
    print()
    print('Results by LOBPCG')
    print()
    print(n,eigs)

    tt = time.clock()
    w,v = symeig(A,B,range=(1,m))
    data2.append(time.clock()-tt)
    print()
    print('Results by symeig')
    print()
    print(n, w)

xlabel(r'Size $n$')
ylabel(r'Elapsed time $t$')
plot(N,data1,label='LOBPCG')
plot(N,data2,label='SYMEIG')
legend()
show()
