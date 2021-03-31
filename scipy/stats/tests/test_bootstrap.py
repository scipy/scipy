import numpy as np
import pytest
from scipy.stats import bootstrap_ci


def test_bootstrap_ci_iv():
    message = "`data` must be a sequence of samples."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(1, np.mean)

    message = "`data` must contain at least one sample."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(tuple(), np.mean)

    message = "each sample in `data` must contain two or more observations..."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3], [1]), np.mean)

    message = "`axis` must be an integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, axis=1.5)

    message = "could not convert string to float"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, confidence_level='ni')

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, n_resamples=-1000)

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, n_resamples=1000.5)

    message = "`method` must be in"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, method='ekki')

    message = "'herring' cannot be used to seed a"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci(([1, 2, 3],), np.mean, random_state='herring')
