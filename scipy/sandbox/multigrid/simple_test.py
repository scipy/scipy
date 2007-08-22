from multilevel import *
from multigrid import *
from scipy import *

A = poisson_problem2D(200)
rs_solver = ruge_stuben_solver(A)
b = rand(A.shape[0])
x,res = rs_solver.solve(b,return_residuals=True)
print res



