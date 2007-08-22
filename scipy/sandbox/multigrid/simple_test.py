from multilevel import *
from multigrid import *
from scipy import *

A = poisson_problem(100).T
s = scalar_solver(A)
b = rand(A.shape[0])
x,res = s.solve(b,return_residuals=True)
r = (b - A*x)
print abs(r).max()




