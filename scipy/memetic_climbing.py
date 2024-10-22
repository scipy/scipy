"""
Memetic Climbing Algorithm
==========================

This module provides an implementation of the memetic climbing algorithm, which combines hill climbing with 
restarts and local search using L-BFGS-B.

Functions
---------
memetic_climbing(objective, bounds, nx, cycles=1000, local_search=True)
    Memetic algorithm that incorporates local search with hill climbing.

hollow_distribution(dimension)
    Generate initial random values for each call.

hill_climbing_with_restarts(objective, bounds, nx, cycles=1000, max_iterations=20000, step_size=20, decay_rate=0.99977, local_search=True)
    Hill climbing optimization with restarts.

Authors
-------
Dawit Gulta (dawit.lambebo@gmail.com)
Stephen Chen (sychen@yorku.ca)
Copyright (c) 2024 Dawit Gulta (dawit.lambebo@gmail.com), Stephen Chen (sychen@yorku.ca)
"""

import numpy as np
from scipy.optimize import minimize

__all__ = ['memetic_climbing']

def hollow_distribution(dimension):
    """
    Generate initial random values for each call
    Parameters
    ----------
    dimension : int
        The dimensionality of the distribution.

    Returns
    -------
    adjusted_numbers : ndarray
        An array of adjusted random values.

    Notes
    -----
    The adjustment is designed to improve exploration by ensuring that values close to
    the center of the distribution are pushed outward. To enhance exploration compared
    to a simple uniform random distribution, the 'next step' was made not just ±10
    from the current location. Instead, it is a uniform value from the range [5, 10],
    applied in either the positive or negative direction.

    """
    
    
    initial_randoms = np.random.uniform(-1, 1, size=dimension)
    adjusted_numbers = np.zeros(dimension)  # Initialize an array of zeros of the correct size

    for i in range(len(initial_randoms)):
        num = initial_randoms[i]
        if -0.5 < num < 0:
            adjusted_numbers[i] = num - 0.5
        elif 0 <= num < 0.5:
            adjusted_numbers[i] = num + 0.5
        else:
            adjusted_numbers[i] = num  # Keep the number as it is if it's already in the desired range

    return adjusted_numbers

def hill_climbing_with_restarts(objective, bounds, nx, cycles=1000, max_iterations=20000, step_size=20, decay_rate=0.99977, local_search=True):
    """
    Hill climbing optimization with restarts.

    This method repeatedly applies hill climbing with random restarts to find the
    minimum of an objective function. Optionally, it performs a local search
    using L-BFGS-B to refine the solution.

    Parameters
    ----------
    objective : callable
        The objective function to be minimized. It should accept different arrays of dimensions.
    bounds : sequence of tuples
        Bounds for the variables. Each tuple represents the (min, max) bounds
        for a single dimension.
    nx : int
        Number of dimensions for the input space.
    cycles : int, optional
        Number of cycles to perform (default is 1000).
    max_iterations : int, optional
        Maximum number of iterations (default is 20000).
    step_size : float, optional
        Step size for hill climbing (default is 20).
    decay_rate : float, optional
        Decay rate for step size (default is 0.99977).
    local_search : bool, optional
        Whether to perform local search using L-BFGS-B (default is True).

    Returns
    -------
    iteration_details_global : list
        List of details for each iteration.
    best_result : dict
        Dictionary containing the best found value, with the key:
        - 'best_value': The minimum value found for the objective function.


    Notes
    -----
    The hill climbing process starts with random points and adjusts them according
    to the specified step size. If `local_search` is enabled, a final optimization
    using L-BFGS-B is performed on the best found solution.
    """
    best_global_value = float('inf')
    iteration_details_global = []
    max_iterations = max_iterations // cycles

    lower_bounds, upper_bounds = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])

    for cycle in range(cycles):

        x0 = np.random.uniform(lower_bounds, upper_bounds, size=nx)
        best_x = np.array(x0)

        best_value = objective(best_x)

        iterations = 0

        while iterations < max_iterations:
            next_x = best_x + hollow_distribution(nx) * step_size
            next_x = np.clip(next_x, a_min=lower_bounds, a_max=upper_bounds)
            next_value = objective(next_x)

            if next_value <= best_value:
                best_x, best_value = next_x, next_value

            iterations += 1
            step_size *= decay_rate

        if local_search:
            result_bfgs = minimize(objective, best_x, bounds=bounds, method='L-BFGS-B')
            best_x_bfgs, best_value_bfgs = result_bfgs.x, result_bfgs.fun

            if best_value_bfgs <= best_value:
                best_x, best_value = best_x_bfgs, best_value_bfgs

        if best_value < best_global_value:
            best_global_value = best_value

    return iteration_details_global, {
        # 'best_x': best_global_x,
        'best_value': best_global_value,
    }

def memetic_climbing(objective, bounds, nx, cycles=1000, local_search=True):

    """Memetic algorithm that incorporates local search with hill climbing.

    Memetic algorithm that incorporates local search with hill climbing.

    This algorithm first applies a hill climbing method with multiple restarts
    and then refines the best solution found using the L-BFGS-B local search method.

    Parameters
    ----------
    objective : callable
        The objective function to be minimized. This should accept a 1D array
        of inputs and return a scalar value.
    bounds : sequence of tuple
        Bounds for the variables, with each tuple representing the (min, max)
        bounds for a single dimension.
    nx : int
        Number of dimensions for the input space.
    cycles : int, optional
        Number of cycles to perform. Default is 1000.
    local_search : bool, optional
        Whether to perform local search using L-BFGS-B after hill climbing.
        Default is True.

    Returns
    -------
    best_result : dict
        A dictionary containing the best found value and corresponding input.
        Keys include:
        - 'best_value': The minimum value found for the objective function.

    Notes
    -----
    This implementation is based on combining simple hill climbing with
    restarts and a final local optimization step using L-BFGS-B. The method
    is particularly effective for problems where the global landscape has
    multiple local minima.

    References
    ----------
    .. [1] SciPy optimization documentation: https://docs.scipy.org/doc/scipy/reference/optimize.html

    """
    if not callable(objective):
        raise TypeError("The objective must be a callable function.")

    if not isinstance(bounds, (list, tuple)) or len(bounds) != nx:
        raise ValueError(f"bounds must be a sequence of {nx} (min, max) tuples.")

    for bound in bounds:
        if not isinstance(bound, (list, tuple)) or len(bound) != 2 or bound[0] > bound[1]:
            raise ValueError("Each bound must be a (min, max) tuple with min <= max.")

    if not isinstance(nx, int) or nx <= 0:
        raise ValueError("nx must be a positive integer.")

    if not isinstance(cycles, int) or cycles <= 0:
        raise ValueError("cycles must be a positive integer.")

    if not isinstance(local_search, bool):
        raise TypeError("local_search must be a boolean.")

    try:
        iteration_details, result = hill_climbing_with_restarts(
            objective, bounds, nx, cycles=cycles, local_search=local_search)
        return result
    except Exception as e:
        raise RuntimeError(f"An error occurred during optimization: {str(e)}")

def test_memetic_climbing():
    """
    Test the memetic climbing optimization algorithm.

    This function tests the `memetic_climbing` function by applying it to a simple
    quadratic objective function. It asserts that the algorithm converges to a
    solution close to the global minimum.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If the optimization does not converge to the expected minimum.

    Examples
    --------
    >>> test_memetic_climbing()
    """
    def objective_function(x):
        return np.sum(x**2)

    bounds = [(-10, 10)] * 10

    result = memetic_climbing(objective_function, bounds, nx=10, cycles=1000, local_search=True)
    expected_value = 7.747023613798194e-10
    tolerance = 1e-8
    assert abs(result['best_value'] <= expected_value), \
        f"Optimization result is {result['best_value']} which is not within the expected tolerance of {expected_value}."