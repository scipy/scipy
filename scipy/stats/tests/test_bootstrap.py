import numpy as np
import pytest
from scipy.stats import bootstrap_ci_1samp

def test_bootstrap_ci_1samp_iv():
    message = "`data` must contain two or more observations along `axis`."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([], np.mean)

    message = "`axis` must be an integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, axis=1.5)

    message = "could not convert string to float"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, confidence_level='ni')

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, n_resamples=-1000)

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, n_resamples=1000.5)

    message = "`method` must be in"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, method='ekki')

    message = "'herring' cannot be used to seed a"
    with pytest.raises(ValueError, match=message):
        bootstrap_ci_1samp([1, 2, 3], np.mean, random_state='herring')
