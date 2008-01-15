from scipy import *
from scipy.sandbox.multigrid.sa import *
from scipy.sandbox.multigrid import *
from scipy.sandbox.multigrid.utils import *
from time import clock

mats = io.loadmat('/home/nathan/Work/ogroup/matrices/elasticity/simple_grid_2d/elasticity_500x500.mat')
A = mats['A'].tobsr(blocksize=(2,2))
B = mats['B']

#A = poisson_problem2D(50)
#B = None

start = clock()
sa = smoothed_aggregation_solver(A,B=B)
print "constructed solver in %s seconds" % (clock() - start)

x0 = rand(A.shape[0])
b = zeros_like(x0)
start = clock()
x,residuals = sa.solve(b,x0=x0,return_residuals=True)
print "solved in %s seconds" % (clock() - start)
    
residuals = array(residuals)/residuals[0]
avg_convergence_ratio = residuals[-1]**(1.0/len(residuals))
print "average convergence ratio",avg_convergence_ratio
print "last convergence ratio",residuals[-1]/residuals[-2]
print 'residuals',residuals
