from __future__ import division
import numpy as np
from scipy.special import logsumexp


def softmax(X, theta=1.0, axis=None):
    r"""
    Softmax function

    Compute the softmax of each element along an axis of X.

    .. math:: \sigma(X)_j = \frac{e^{\theta x_j}}{\sum_{k=1}^{K}e^{\theta x_k}} \textrm{ for } j = 1, \ldots, K

    Parameters
    ----------
    X : array_like
        Input array.
    theta : float, optional
            parameter, used as a multiplier prior to exponentiation. Default = 1.0
    axis : int, optional
           axis to compute values along. Default is the first non-singleton axis.

    Returns
    -------
    s : ndarray
        an array the same size as X. The result will sum to 1 along the specified axis.
    """
    X = np.asarray(X)

    # make X at least 2d
    y = np.atleast_2d(X)

    # find axis
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # multiply y against the theta parameter,
    y = y * float(theta)

    # compute in log space for numerical stability
    ax_sum = logsumexp(y, axis=axis)
    sigma = np.exp(y - ax_sum)

    # flatten if X was 1D
    if X.ndim == 1:
        sigma = sigma.flatten()

    return sigma
