from __future__ import division, print_function, absolute_import

from scipy import *
from scipy.linalg import eigh, orth, cho_factor, cho_solve
from scipy.sparse.linalg import lobpcg
#from pylab import plot, show, legend, xlabel, ylabel
set_printoptions(precision=3,linewidth=90)
import time


def _mikota_pair(n):
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


#def precond(x):
    #y = cho_solve((LorU, lower),x)
    #return as2d(y)

m = 10  # Blocksize
N = array(([128,256,512,1024,2048]))  # Increasing matrix size

data1 = []
data2 = []

for n in N:
    print('******', n)
    A, B = _mikota_pair(n)
    X = rand(n,m)
    X = orth(X)

    tt = time.clock()
    # For now skip the preconditioner until the general linear function
    # interface is imporoved.
    #(LorU, lower) = cho_factor(A, lower=0, overwrite_a=0)
    precond = None
    eigs, vecs = lobpcg(A, X, B, M=precond, tol=1e-4, maxiter=40)
    elapsed = time.clock() - tt
    data1.append(elapsed)
    eigs = sort(eigs)
    print()
    print('LOBPCG results:')
    print(n, eigs, elapsed)

    tt = time.clock()
    w = eigh(A, B, eigvals_only=True, eigvals=(0, m-1))
    elapsed = time.clock() - tt
    data2.append(elapsed)
    print()
    print('eigh results:')
    print(n, w, elapsed)

#xlabel(r'Size $n$')
#ylabel(r'Elapsed time $t$')
#plot(N,data1,label='LOBPCG')
#plot(N,data2,label='SYMEIG')
#legend()
#show()
