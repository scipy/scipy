from scipy.optimize import minimize as scipy_minimize
from pyprima import minimize, LinearConstraint, Bounds, NonlinearConstraint
from optiprofiler.problems import load_cutest_problem
from optiprofiler.utils import ProblemError
import numpy as np
import debugpy
from time import time



debugpy.breakpoint()
def get_constraints(problem):
    constraints = []
    if problem.m_linear_ub > 0:
        constraints.append(LinearConstraint(problem.a_ub, -np.inf, problem.b_ub))
    if problem.m_linear_eq > 0:
        constraints.append(LinearConstraint(problem.a_eq, -np.inf, problem.b_eq))
        constraints.append(LinearConstraint(-problem.a_eq, -np.inf, problem.b_eq))
    if problem.m_nonlinear_ub > 0:
        constraints.append(NonlinearConstraint(problem.c_ub, -np.inf, np.zeros(problem.m_nonlinear_ub)))
    if problem.m_nonlinear_eq > 0:
        constraints.append(NonlinearConstraint(problem.c_eq, -np.inf, np.zeros(problem.m_nonlinear_eq)))
        constraints.append(NonlinearConstraint(lambda x: -problem.c_eq(x), -np.inf, np.zeros(problem.m_nonlinear_eq)))
    return constraints


nondefault_options = lambda n, f0: {
    'ftarget' : f0 - 314, # if this is doable, 3.14 otherwise
    'maxfev' : 271*n,
    # 'npt' : int(min(3.14*n, n**1.23)),
    'rhobeg' : 2.71828,
    # 'ctol': 2e-4,
    'rhoend' : 3.14e-4,
    # 'iprint' : 1
}

with open('cobyla.txt') as f:
    problems = f.read().splitlines()


f = open('results.csv', 'w')
f.write(", ".join(["Problem", "PyPRIMA result", "SciPy result", "Error", "Error != 0", "PyPRIMA cstrv", "SciPy cstrv", "Error", "PyPRIMA nfev", "SciPy nfev", "PyPRIMA s", "SciPy s, Speedup (<1)/Slowdown (>1)"]))
f.write("\n")


for problem_name in problems:
    try:
        problem = load_cutest_problem(problem_name)
    except ProblemError:
        continue
    cons = get_constraints(problem)
    bounds = Bounds(problem.lb, problem.ub)
    options = {}  #nondefault_options(problem.n, problem.fun(problem.x0))
    print("Solving with bindings")
    time1 = time()
    result_1 = minimize(problem.fun, problem.x0, method='cobyla', options=options, constraints=cons, bounds=bounds)
    time1 = (time() - time1)/result_1.nf
    print("Solving with scipy")
    options = {}  #nondefault_options(problem.n, problem.fun(problem.x0))
    time2 = time()
    result_2 = scipy_minimize(problem.fun, problem.x0, method='cobyla', options=options, constraints=cons, bounds=bounds)
    time2 = (time() - time2)/result_2.nfev

    # Do some math. The math is designed such that a negative value for error means the
    # first algorithm tested is better.
    funerror = (result_1.f - result_2.fun)/abs(result_1.f) if result_1.f != 0 else 0
    funerror0 = funerror != 0
    cstrerror = (result_1.cstrv - result_2.maxcv)/abs(result_1.cstrv) if result_1.cstrv != 0 else 0
    speed_change = time1/time2  # If first one is faster, this is < 1
    f.write(f'{problem_name: <11}, {result_1.f: <23}, {result_2.fun: <23}, {funerror*100: <8.2f}, {funerror0}, {result_1.cstrv: <23}, {result_2.maxcv: <23}, {cstrerror*100: <8.2f}, {result_1.nf: <5}, {result_2.nfev: <5}, {time1: <5}, {time2: <5}, {speed_change: <5}\n')
    f.flush()


    # for result in results:
    #     problem_name, result_1_fun, presult_f, funerror, result_1_maxcv, presult_cstrv, cstrerror, result_1_nfev, presult_nf = result
    #     f.write(f'{problem_name: <11}, {result_1_fun: <23}, {presult_f: <23}, {funerror*100: <8.2f}, {result_1_maxcv: <23}, {presult_cstrv: <23}, {cstrerror*100: <8.2f}, {result_1_nfev: <5}, {presult_nf: <5}\n')
