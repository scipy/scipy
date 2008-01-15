
from scipy import *
from scipy.sandbox.multigrid.utils import diag_sparse
from scipy.sandbox.multigrid.gallery import poisson, linear_elasticity


#A = poisson( (200,200), spacing=(1,1e-2)
#aggregation = [ sa_constant_interpolation(A*A*A,epsilon=0.0) ]

#A = io.mmread("tests/sample_data/laplacian_41_3dcube.mtx").tocsr()
#A = io.mmread("laplacian_40_3dcube.mtx").tocsr()
#A = io.mmread("/home/nathan/Desktop/9pt/9pt-100x100.mtx").tocsr()
#A = io.mmread("/home/nathan/Desktop/BasisShift_W_EnergyMin_Luke/9pt-5x5.mtx").tocsr()


#D = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
#A = D * A * D

#A = io.mmread("tests/sample_data/elas30_A.mtx").tocsr()

A,B = linear_elasticity( (100,100) )

from time import clock; start = clock()
asa = adaptive_sa_solver(A,max_candidates=3,mu=5,blocks=blocks,aggregation=aggregation)
print "Adaptive Solver Construction: %s seconds" % (clock() - start); del start

scipy.random.seed(0)  #make tests repeatable
x = randn(A.shape[0])
b = A*randn(A.shape[0])
#b = zeros(A.shape[0])


print "solving"
if False:
    x_sol,residuals = asa.solver.solve(b,x0=x,maxiter=20,tol=1e-12,return_residuals=True)
else:
    residuals = []
    def add_resid(x):
        residuals.append(linalg.norm(b - A*x))
    A.psolve = asa.solver.psolve
    x_sol = linalg.cg(A,b,x0=x,maxiter=30,tol=1e-12,callback=add_resid)[0]

residuals = array(residuals)/residuals[0]

print "residuals ",residuals
print "mean convergence factor",(residuals[-1]/residuals[0])**(1.0/len(residuals))
print "last convergence factor",residuals[-1]/residuals[-2]

print
print asa.solver

print "constant Rayleigh quotient",dot(ones(A.shape[0]),A*ones(A.shape[0]))/float(A.shape[0])

def plot2d_arrows(x):
    from pylab import figure,quiver,show
    x = x.reshape(-1)
    N = (len(x)/2)**0.5
    assert(2 * N * N == len(x))
    X = linspace(-1,1,N).reshape(1,N).repeat(N,0).reshape(-1)
    Y = linspace(-1,1,N).reshape(1,N).repeat(N,0).T.reshape(-1)

    dX = x[0::2]
    dY = x[1::2]

    figure()
    quiver(X,Y,dX,dY)
    show()

def plot2d(x):
    from pylab import pcolor,figure,show
    figure()
    pcolor(x.reshape(sqrt(len(x)),sqrt(len(x))))
    show()


for c in asa.Bs[0].T:
    #plot2d(c)
    plot2d_arrows(c)
    print "candidate Rayleigh quotient",dot(c,A*c)/dot(c,c)


##W = asa.AggOps[0]*asa.AggOps[1]
##pcolor((W * rand(W.shape[1])).reshape((200,200)))

