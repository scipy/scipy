from scipy.optimize import rosen, Bounds, NonlinearConstraint, LinearConstraint
from pyprima import minimize
import numpy as np
import concurrent.futures


def test_threading():
    # This is to validate that python cobyla does not require locks.
    # We will solve a set of problems in a single threaded fashion and get their
    # solutions, then we will solve the same problems in a thread pool, and validate
    # that we get the same results.

    problem1 = {
        'fun': rosen,
        'x0': np.linspace(3.4, 17.8, 10),
        'method': 'cobyla',
        'bounds': None,
        'constraints': None
    }
    problem2 = {
        'fun': rosen,
        'x0': np.linspace(-178, -23, 7),
        'method': 'cobyla',
        'bounds': None,
        'constraints': None
    }
    nlc = NonlinearConstraint(lambda x: x[0]**2 - 5.0, 0.0, np.inf)
    problem3 = {
        'fun': rosen,
        'x0': np.linspace(-178, -23, 7),
        'method': 'cobyla',
        'bounds': None,
        'constraints': nlc
    }
    bounds = Bounds(np.array([0, 0]), np.array([10, 10]))
    lc = LinearConstraint([1]*7, 50, 50)  # Sum of optimum point is 50
    problem4 = {
        'fun': rosen,
        'x0': np.linspace(-178, -23, 7),
        'method': 'cobyla',
        'bounds': bounds,
        'constraints': lc
    }
    problems = [problem1, problem2, problem3, problem4]

    def run_problem(problem):
        return minimize(**problem)

    single_threaded_results = [run_problem(problem) for problem in problems]

    with concurrent.futures.ThreadPoolExecutor() as executor:
        multithreaded_results = executor.map(run_problem, problems)

    for single_threaded_result, multithreaded_result in zip(single_threaded_results, multithreaded_results):
        assert np.allclose(single_threaded_result.x, multithreaded_result.x)
        assert np.allclose(single_threaded_result.constr, multithreaded_result.constr)
        assert single_threaded_result.nf == multithreaded_result.nf
        assert single_threaded_result.f == multithreaded_result.f