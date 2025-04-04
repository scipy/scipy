import pycutest
from pyprima import minimize as pyprima_minimize
from prima import minimize, LinearConstraint, Bounds, NonlinearConstraint
from optiprofiler.problems import load_cutest_problem
import numpy as np
import debugpy
import sys

np.set_printoptions(precision=53, floatmode='fixed', suppress=False)

debug = input("Type Y to debug, ENTER to continue: ")
if debug == 'Y':
    print("Waiting for debugger to attach...")
    debugpy.listen(5678)
    debugpy.wait_for_client()

sys.modules['pyprima'].common.linalg.COMPARING = False
problem_name = 'POLAK3'
problem = load_cutest_problem(problem_name)

if debug == 'Y':
    debugpy.breakpoint()
constraints = []
if problem.m_linear_ub > 0:
    print("Adding linear inequality constraints")
    constraints.append(LinearConstraint(problem.a_ub, -np.inf, problem.b_ub))
if problem.m_linear_eq > 0:
    print("Adding linear equality constraints")
    constraints.append(LinearConstraint(problem.a_eq, problem.b_eq, problem.b_eq))
if problem.m_nonlinear_ub > 0:
    print("Adding nonlinear inequality constraints")
    constraints.append(NonlinearConstraint(problem.c_ub, -np.inf, np.zeros(problem.m_nonlinear_ub)))
if problem.m_nonlinear_eq > 0:
    print("Adding nonlinear equality constraints")
    constraints.append(NonlinearConstraint(problem.c_eq, np.zeros(problem.m_nonlinear_eq), np.zeros(problem.m_nonlinear_eq)))
bounds = Bounds(problem.lb, problem.ub)

x0 = problem.x0
f0 = problem.fun(x0)

nondefault_options = lambda n, f0: {
    'ftarget' : f0 - 314, # if this is doable, 3.14 otherwise
    'maxfev' : 271*n,
    # 'npt' : int(min(3.14*n, n**1.23)),
    'rhobeg' : 2.71828,
    # 'ctol': 2e-4,
    'rhoend' : 3.14e-4,
    'iprint' : 1
}


python_options = nondefault_options(len(x0), f0)
# del python_options['npt']


fortran_options = nondefault_options(len(x0), f0)

# def fun(x):
#         return x[0]**2 + abs(x[1])**3

# def con1(x):
#     return x[0]**2 + x[1]**2 - 25

# def con2(x):
#     return -con1(x)

# x0 = [np.sqrt(25 - (2.0/3)**2), 2.0/3 + 1e-4]

# bounds=None

# constraints = [NonlinearConstraint(con1, -np.inf, 0),
#                 NonlinearConstraint(con2, -np.inf, 0)]

# print(con1(x0))

fun = problem.fun

if sys.argv[1] == 'p':
    result = pyprima_minimize(fun, x0, method='cobyla', options=python_options, constraints=constraints, bounds=bounds)
    print(result.cstrv)
    from pyprima.common.linalg import matprod
    print(matprod.counter)
else:
    result = minimize(fun, x0, method='cobyla', options=fortran_options, constraints=constraints, bounds=bounds)
    print(result.maxcv)
