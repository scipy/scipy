from multilevel import *
from scipy import *

A = poisson_problem2D(200)
rs_solver = ruge_stuben_solver(A)
b = rand(A.shape[0])
x,residuals = rs_solver.solve(b,return_residuals=True)
print 'residuals',residuals



