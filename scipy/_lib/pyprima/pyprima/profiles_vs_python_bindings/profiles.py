import numpy as np
import sys
import os
from excludelist import excludelist
from optiprofiler import set_cutest_problem_options, find_cutest_problems, run_benchmark
import argparse
from time import time


nondefault_options = lambda n, f0: {
    'ftarget' : f0 - 314, # if this is doable, 3.14 otherwise
    'maxfev' : 271*n,
    'npt' : int(min(3.14*n, n**1.23)),
    'rhobeg' : 2.71828,
    'rhoend' : 3.14159*1.0e-4,
}

def get_pyprima_options(n, f0):
    if os.environ.get('NONDEFAULT_PYPRIMA') == 'True':
        options = nondefault_options(n, f0)
        # Change the option name
        options['maxfun'] = options.pop('maxfev')
        return options
    return {}


def get_python_options(n, f0):
    if os.environ.get('NONDEFAULT_PYTHON') == 'True':
        return nondefault_options(n, f0)
    return {}


def pyprima_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from pyprima import minimize, Bounds, LinearConstraint, NonlinearConstraint

    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    options = get_pyprima_options(len(x0), f0)
    if 'npt' in options:
        del options['npt']
    result = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=options)
    return result.x


def python_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint

    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    res = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=get_python_options(len(x0), f0))
    return res.x


def get_problems(description):
    cutest_problem_names = find_cutest_problems(description)
    return list(filter(lambda x: x not in excludelist(description), cutest_problem_names))


if __name__ == '__main__':
    # If we run this script from a directory other than the one that contains it, pycutest's call to importlib will fail,
    # unless we insert the current working directory into the path.
    sys.path.insert(0, os.getcwd())
    os.environ['PYCUTEST_CACHE'] = os.getcwd()
    
    parser = argparse.ArgumentParser(description='Generate performance profiles comparing PyPRIMA to PRIMA Python (bindings).')
    parser.add_argument('-j', '--n_jobs', type=int, default=None, help='Number of jobs to run in parallel')
    parser.add_argument('--default_only', action='store_true', help='Run only the default options for both PyPRIMA and PRIMA')
    args = parser.parse_args()
    

    def run_three_benchmarks(pyprima_fun, python_fun, algorithm, cutest_problem_names, default_only, n_jobs):
        '''
        Proper validation of both default and nondefault options requires 3 runs: Both default, both nondefault, and
        one default one nondefault. The first two should look identical, and so the third run confirms that our
        experiment setup is valid (i.e. it rules out a scenario where even though options are provided, both algorithms
        end up using default options anyway).
        '''
        # Sharing state with multiprocessing is hard when we can't control the function signature,
        # so we resort to using the environment to pass options.
        algorithm = algorithm.lower()
        ALGORITHM = algorithm.upper()
        project_x0 = algorithm == 'lincoa'
        os.environ['NONDEFAULT_PYPRIMA'] = "False"
        os.environ['NONDEFAULT_PYTHON'] = "False"
        run_benchmark([pyprima_fun, python_fun], [f'PyPRIMA-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_default_options', n_jobs=n_jobs, project_x0=project_x0)
        if not default_only:
            os.environ['NONDEFAULT_PYPRIMA'] = "True"
            os.environ['NONDEFAULT_PYTHON'] = "True"
            run_benchmark([pyprima_fun, python_fun], [f'PyPRIMA-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_nondefault_options', n_jobs=n_jobs, project_x0=project_x0)
            os.environ['NONDEFAULT_PYPRIMA'] = "True"
            os.environ['NONDEFAULT_PYTHON'] = "False"
            run_benchmark([pyprima_fun, python_fun], [f'PyPRIMA-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_different_options', n_jobs=n_jobs, project_x0=project_x0)

    start = time()
    print("Running profiles for COBYLA")
    with open('cobyla.txt') as f:
        cutest_problem_names = f.read().splitlines()
    cutest_problem_names = list(filter(lambda x: x not in excludelist('cobyla'), cutest_problem_names))
    run_three_benchmarks(pyprima_cobyla, python_cobyla, 'cobyla', cutest_problem_names, args.default_only, args.n_jobs)
    print(f'Completed COBYLA profile in {time() - start:.2f} seconds')
