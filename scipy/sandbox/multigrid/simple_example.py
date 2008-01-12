from scipy import rand
from scipy.sandbox.multigrid.sa import *
from scipy.sandbox.multigrid import *
from scipy.sandbox.multigrid.utils import *
from time import clock

mats = io.loadmat('/home/nathan/Work/ogroup/matrices/elasticity/simple_grid_2d/elasticity_50x50.mat')
A = mats['A'].tobsr(blocksize=(2,2))
B = mats['B']

#A = poisson_problem2D(500)
#B = None

start = clock()
sa = smoothed_aggregation_solver(A,B=B)
print "constructed solver in %s seconds" % (clock() - start)

b = rand(A.shape[0])
start = clock()
x,residuals = sa.solve(b,return_residuals=True)
print "solved in %s seconds" % (clock() - start)
print 'residuals',residuals
