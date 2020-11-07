
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import relative_risk


# Test just the calculation of the relative risk, including edge
# cases that result in a relative risk of 0, inf or nan.
@pytest.mark.parametrize(
    'exposed_cases, exposed_total, control_cases, control_total, expected_rr',
    [(1, 4, 3, 8, 0.25 / 0.375),
     (0, 10, 5, 20, 0),
     (0, 10, 0, 20, np.nan),
     (5, 15, 0, 20, np.inf)]
)
def test_relative_risk(exposed_cases, exposed_total,
                       control_cases, control_total, expected_rr):
    result = relative_risk(exposed_cases, exposed_total,
                           control_cases, control_total)
    assert_allclose(result.relative_risk, expected_rr, rtol=1e-13)


def test_relative_risk_confidence_interval():
    exposed_cases = 16
    exposed_total = 128
    control_cases = 24
    control_total = 256
    result = relative_risk(exposed_cases, exposed_total,
                           control_cases, control_total)
    rr = result.relative_risk
    ci = result.confidence_interval(confidence_level=0.95)
    # The corresponding calculation in R using the epitools package.
    #
    # > library(epitools)
    # > c <- matrix(c(232, 112, 24, 16), nrow=2)
    # > result <- riskratio(c)
    # > result$measure
    #               risk ratio with 95% C.I.
    # Predictor  estimate     lower    upper
    #   Exposed1 1.000000        NA       NA
    #   Exposed2 1.333333 0.7347317 2.419628
    #
    # The last line is the result that we want.
    assert_allclose(rr, 4/3)
    assert_allclose(ci, (0.7347317, 2.419628), rtol=5e-7)
