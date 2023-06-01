import numpy as np
import scipy.special as sc


def test_softplus():
    """
    Test cases for the softplus function.
    """
    # Test case 1: Test softplus with a single value
    result = sc.softplus(0)
    expected = 0.6931471805599453
    assert np.allclose(result, expected)

    # Test case 2: Test softplus with an array of values
    result = sc.softplus([-1, 0, 1])
    expected = np.array([0.31326169, 0.69314718, 1.31326169])
    assert np.allclose(result, expected)

    # Test case 3: Test softplus with a large positive value
    result = sc.softplus(100)
    expected = 100.0  # The softplus of a large positive value approaches the value itself
    assert np.allclose(result, expected)

    # Test case 4: Test softplus with a negative value
    result = sc.softplus(-5)
    expected = 0.00671534848911849  # The softplus of a negative value is a small positive value
    assert np.allclose(result, expected)

    # Test case 5: Test softplus with a large negative value
    result = sc.softplus(-100)
    expected = 0.0  # The softplus of a large negative value approaches 0 due to underflow
    assert np.allclose(result, expected)

    # Test case 6: Test softplus with a very large positive value causing overflow
    result = sc.softplus(10000)
    expected = 10000.0  # The softplus of a very large positive value approaches the value itself
    assert np.allclose(result, expected)

