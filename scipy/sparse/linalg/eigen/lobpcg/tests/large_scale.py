from scipy import array, arange, ones, sort, cos, pi, rand, \
     set_printoptions, r_
from scipy.sparse.linalg import lobpcg
from scipy import sparse
from pylab import loglog, show, xlabel, ylabel, title
set_printoptions(precision=8,linewidth=90)
import time

def sakurai(n):
    """ Example taken from
        T. Sakurai, H. Tadano, Y. Inadomi and U. Nagashima
        A moment-based method for large-scale generalized eigenvalue problems
        Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004) """

    A = sparse.eye( n, n )
    d0 = array(r_[5,6*ones(n-2),5])
    d1 = -4*ones(n)
    d2 =  ones(n)
    B = sparse.spdiags([d2,d1,d0,d1,d2],[-2,-1,0,1,2],n,n)

    k = arange(1,n+1)
    w_ex = sort(1./(16.*pow(cos(0.5*k*pi/(n+1)),4))) # exact eigenvalues

    return A,B, w_ex

m = 3  # Blocksize

#
# Large scale
#
n = 2500
A,B, w_ex = sakurai(n) # Mikota pair
X = rand(n,m)
data=[]
tt = time.clock()
eigs,vecs, resnh = lobpcg(X,A,B, residualTolerance = 1e-6, maxIterations =500, retResidualNormsHistory=1)
data.append(time.clock()-tt)
print 'Results by LOBPCG for n='+str(n)
print
print eigs
print
print 'Exact eigenvalues'
print
print w_ex[:m]
print
print 'Elapsed time',data[0]
loglog(arange(1,n+1),w_ex,'b.')
xlabel(r'Number $i$')
ylabel(r'$\lambda_i$')
title('Eigenvalue distribution')
show()
