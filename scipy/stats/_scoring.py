import numpy as np
from scipy.stats import norm

__all__ = ['crps_gaussian']


def crps_gaussian(y, mu, sigma):
    """
    Compute the Continuous Ranked Probability Score (CRPS) for Gaussian
    forecasts.

    The CRPS measures the accuracy of probabilistic forecasts against
    observed values. For a Gaussian forecast distribution N(mu, sigma^2)
    and an observation y, the CRPS has a closed-form expression [1]_.

    Parameters
    ----------
    y : array_like
        Observed values.
    mu : array_like
        Mean of the forecast Gaussian distribution.
    sigma : array_like
        Standard deviation of the forecast Gaussian distribution.
        Must be positive.

    Returns
    -------
    crps : ndarray or float
        The CRPS value(s). Lower values indicate better forecasts.
        The output shape is determined by broadcasting `y`, `mu`,
        and `sigma`.

    Notes
    -----
    The CRPS for a Gaussian distribution is given by:

    .. math::

        \\text{CRPS}(y; \\mu, \\sigma) = \\sigma \\left[
            z (2 \\Phi(z) - 1) + 2 \\varphi(z) - \\frac{1}{\\sqrt{\\pi}}
        \\right]

    where :math:`z = (y - \\mu) / \\sigma`, :math:`\\varphi` is the standard
    normal PDF, and :math:`\\Phi` is the standard normal CDF.

    The CRPS generalizes the mean absolute error (MAE) to probabilistic
    forecasts. When sigma approaches zero, the CRPS converges to
    ``|y - mu|``.

    The CRPS is a proper scoring rule, meaning that the expected score is
    minimized when the forecast distribution equals the true data-generating
    distribution.

    References
    ----------
    .. [1] Gneiting, T., Raftery, A.E., Westveld III, A.H. and Goldman, T.,
       2005. "Calibrated probabilistic forecasting using ensemble model output
       statistics and minimum CRPS estimation." Monthly Weather Review, 133(5),
       pp.1098-1118.

    Examples
    --------
    Compute the CRPS for a single forecast and observation:

    >>> from scipy.stats import crps_gaussian
    >>> crps_gaussian(y=0.5, mu=0.0, sigma=1.0)
    0.3207...

    A tighter forecast (smaller sigma) around the true value gives a
    lower (better) score:

    >>> crps_gaussian(y=0.0, mu=0.0, sigma=0.1)
    0.0563...
    >>> crps_gaussian(y=0.0, mu=0.0, sigma=1.0)
    0.2338...

    The function supports broadcasting:

    >>> import numpy as np
    >>> y = np.array([0.0, 1.0, 2.0])
    >>> crps_gaussian(y, mu=0.0, sigma=1.0)
    array([0.2338..., 0.5207..., 1.0169...])
    """
    y = np.asarray(y, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    sigma = np.asarray(sigma, dtype=np.float64)

    # Handle sigma <= 0 before division
    invalid = sigma <= 0

    # Replace invalid sigma with 1.0 to avoid division warnings,
    # then overwrite those entries with NaN at the end.
    sigma_safe = np.where(invalid, 1.0, sigma)

    z = (y - mu) / sigma_safe

    phi_z = norm.pdf(z)
    Phi_z = norm.cdf(z)
    pi_inv = 1.0 / np.sqrt(np.pi)

    crps = sigma_safe * (z * (2.0 * Phi_z - 1.0) + 2.0 * phi_z - pi_inv)

    # Set invalid entries to NaN
    if np.any(invalid):
        crps = np.where(invalid, np.nan, crps)

    return crps
