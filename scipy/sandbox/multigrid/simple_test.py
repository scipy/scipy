from scipy import *
from scipy.sandbox.multigrid.sa import *
from scipy.sandbox.multigrid import *
from scipy.sandbox.multigrid.utils import *
from time import clock

A = poisson_problem2D(500)

start = clock()
sa = smoothed_aggregation_solver(A)
print "constructed solver in %s seconds" % (clock() - start)

b = rand(A.shape[0])
start = clock()
x,residuals = sa.solve(b,return_residuals=True)
print "solved in %s seconds" % (clock() - start)
print 'residuals',residuals
