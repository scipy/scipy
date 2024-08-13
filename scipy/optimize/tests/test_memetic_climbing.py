# test_memetic_climbing.py

import numpy as np
from scipy.optimizer.memetic_climbing import memetic_climbing

def test_memetic_climbing():
    """
    Test case for the memetic_climbing algorithm.

    This test verifies that the memetic_climbing algorithm converges to
    a value close to the expected minimum for a simple objective function.

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
    assert abs(result['best_value'] - expected_value) <= tolerance, \
        f"Optimization result is {result['best_value']} which is not within the expected tolerance of {expected_value}."

# Run the test function if the script is executed directly
if __name__ == "__main__":
    test_memetic_climbing()
