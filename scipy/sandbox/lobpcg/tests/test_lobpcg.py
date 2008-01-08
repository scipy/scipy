from scipy import array, arange, ones, sort, cos, pi, rand, \
     set_printoptions, r_, diag, linalg
from scipy.sandbox import lobpcg
from symeig import symeig
from pylab import plot, show, legend, xlabel, ylabel
set_printoptions(precision=3,linewidth=90)

def check1(n):
    L = 1.0
    le=L/n
    rho = 7.85e3
    S = 1.e-4
    E = 2.1e11
    mass = rho*S*le/6.
    k = E*S/le
    A = k*(diag(r_[2.*ones(n-1),1])-diag(ones(n-1),1)-diag(ones(n-1),-1))
    B = mass*(diag(r_[4.*ones(n-1),2])+diag(ones(n-1),1)+diag(ones(n-1),-1))
    return A,B

def check2(n):
    x = arange(1,n+1)
    B = diag(1./x)
    y = arange(n-1,0,-1)
    z = arange(2*n-1,0,-2)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A,B

n = 100 # Dimension

def test_1and2():
    A,B = check1(n) # Fixed-free elastic rod
    A,B = check2(n) # Mikota pair acts as a nice test since the eigenvalues are the squares of the integers n, n=1,2,...
    
    m = 20
    V = rand(n,m)
    X = linalg.orth(V)
    
    eigs,vecs = lobpcg.lobpcg(X,A,B)
    eigs = sort(eigs)
    
    w,v=symeig(A,B)

    plot(arange(0,len(w[:m])),w[:m],'bx',label='Results by symeig')
    plot(arange(0,len(eigs)),eigs,'r+',label='Results by lobpcg')
    legend()
    xlabel(r'Eigenvalue $i$')
    ylabel(r'$\lambda_i$')
    show()
