# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:06:05 2021

@author: matth
"""

from __future__ import annotations
import math
import numpy as np
from scipy import special
from typing import Optional, Union

__all__ = ['entropy', 'differential_entropy']


def entropy(pk, qk=None, base=None, axis=0):
    """Calculate the entropy of a distribution for given probability values.

    If only probabilities `pk` are given, the entropy is calculated as
    ``S = -sum(pk * log(pk), axis=axis)``.

    If `qk` is not None, then compute the Kullback-Leibler divergence
    ``S = sum(pk * log(pk / qk), axis=axis)``.

    This routine will normalize `pk` and `qk` if they don't sum to 1.

    Parameters
    ----------
    pk : sequence
        Defines the (discrete) distribution. ``pk[i]`` is the (possibly
        unnormalized) probability of event ``i``.
    qk : sequence, optional
        Sequence against which the relative entropy is computed. Should be in
        the same format as `pk`.
    base : float, optional
        The logarithmic base to use, defaults to ``e`` (natural logarithm).
    axis: int, optional
        The axis along which the entropy is calculated. Default is 0.

    Returns
    -------
    S : float
        The calculated entropy.

    Examples
    --------

    >>> from scipy.stats import entropy

    Bernoulli trial with different p.
    The outcome of a fair coin is the most uncertain:

    >>> entropy([1/2, 1/2], base=2)
    1.0

    The outcome of a biased coin is less uncertain:

    >>> entropy([9/10, 1/10], base=2)
    0.46899559358928117

    Relative entropy:

    >>> entropy([1/2, 1/2], qk=[9/10, 1/10])
    0.5108256237659907

    """
    pk = np.asarray(pk)
    pk = 1.0*pk / np.sum(pk, axis=axis, keepdims=True)
    if qk is None:
        vec = special.entr(pk)
    else:
        qk = np.asarray(qk)
        pk, qk = np.broadcast_arrays(pk, qk)
        qk = 1.0*qk / np.sum(qk, axis=axis, keepdims=True)
        vec = special.rel_entr(pk, qk)
    S = np.sum(vec, axis=axis)
    if base is not None:
        S /= np.log(base)
    return S


def differential_entropy(
    values: np.typing.ArrayLike,
    *,
    window_length: Optional[int] = None,
    base: Optional[float] = None,
    axis: int = 0,
) -> Union[np.number, np.ndarray]:
    r"""Given a sample of a distribution, calculate the differential entropy.

    This routine uses the Vasicek estimator of the differential entropy. Given
    a sorted random sample :math:`X_1, \ldots X_n` this is defined as:

    .. math::
        \frac{1}{n}\sum_1^n \log\left[ \frac{n}{2m} (X_{i+m} - X_{i-m}) \right]

    where :math:`m` is the window length parameter, :math:`X_{i} = X_1` if
    :math:`i < 1` and :math:`X_{i} = X_n` if :math:`i > n`.

    Parameters
    ----------
    values : sequence
        Samples of the (continuous) distribution.
    window_length : int, optional
        Window length for computing Vasicek estimate. Must be an integer
        between 1 and half of the sample size. If ``None`` (the default) it
        uses the heuristic value

        .. math::
            \left \lfloor \sqrt{n} + 0.5 \right \rfloor

        where :math:`n` is the sample size. This heuristic was originally
        proposed in [2]_ and has become common in the literature.
    base : float, optional
        The logarithmic base to use, defaults to ``e`` (natural logarithm).
    axis : int, optional
        The axis along which the differential entropy is calculated.
        Default is 0.

    Returns
    -------
    entropy : float
        The calculated differential entropy.

    Notes
    -----
    This function will converge to the true differential entropy in the limit
    when

    .. math::
        n \to \infty, \quad m \to \infty, \quad \frac{m}{n} \to 0

    The optimal choice of ``window_length`` for a given sample size, depends on
    the (unknown) distribution. In general, the smoother the density of the
    distribution, the larger is such optimal value of ``window_length`` [1]_.

    References
    ----------
    .. [1] Vasicek, O. (1976). A test for normality based on sample entropy.
           Journal of the Royal Statistical Society:
           Series B (Methodological), 38(1), 54-59.
    .. [2] Crzcgorzewski, P., & Wirczorkowski, R. (1999). Entropy-based
           goodness-of-fit test for exponentiality. Communications in
           Statistics-Theory and Methods, 28(5), 1183-1202.

    Examples
    --------
    >>> from scipy.stats import differential_entropy, norm

    Entropy of a standard normal distribution:

    >>> # We use high numbered seeds in order to have good entropy in the RNG
    >>> rng = np.random.default_rng(seed=148632)
    >>> values = rng.standard_normal(100)
    >>> differential_entropy(values)
    1.401904073487716

    Compare with the true entropy:

    >>> float(norm.entropy())
    1.4189385332046727

    """
    values = np.asarray(values)
    values = np.moveaxis(values, axis, 0)

    if window_length is None:
        window_length = math.floor(math.sqrt(len(values)) + 0.5)

    if not 2 <= 2 * window_length < len(values):
        raise ValueError(
            f"Window length ({window_length}) must be positive and less "
            f"than half the sample size ({len(values)}).",
        )

    sorted_data = np.sort(values, axis=0)

    repeats = np.array(
        (window_length + 1,)
        + ((1,) * (len(sorted_data) - 2))
        + (window_length + 1,),
    )

    padded_data = np.repeat(
        sorted_data,
        repeats=repeats,
        axis=0,
    )

    differences = (
        padded_data[2 * window_length:] -
        padded_data[:-2 * window_length]
    )

    logs = np.log(
        len(differences) * differences / (2 * window_length),
    )
    if base is not None:
        logs /= np.log(base)

    return np.mean(logs, axis=0)
