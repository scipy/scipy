"""
Hypoexponential distribution: the sum of independent exponential random
variables with (possibly different, possibly repeated) scale parameters.
"""
import numpy as np
from scipy._lib import doccer
from scipy.special import gammaln, xlogy

from ._multivariate import (multi_rv_generic, multi_rv_frozen, _squeeze_output,
                            _doc_random_state)

__all__ = ['hypoexpon']

_EPS = np.finfo(float).eps

# Relative rounding-error bound above which the floating point closed form is
# not trusted and the point is recomputed with the phase-type method.
_CLOSED_FORM_RTOL = 1e-12

# Beyond this value of ``x * max(rate)`` the uniformization series is not
# used; its rounding error grows roughly like ``x * max(rate) * eps``.
_SERIES_MAX_LAMBDA = 1e5
_SERIES_MAX_BLOCK_ENTRIES = 2 ** 22
_SQUARING_MAX_BLOCK_ENTRIES = 2 ** 22

# Rough cost model (seconds) used only to split the evaluation points between
# the two phase-type variants, see `_phase_type_plan`. Only ratios matter.
_SERIES_FIXED_PER_K = 6.5e-6   # one step of the c_k recursion
_SERIES_PER_TERM = 4e-9        # one Poisson term for one point
_SERIES_PER_BLOCK = 5e-5       # fixed overhead of one block of points
_SQUARING_PER_POINT = 2e-6     # one point costs this ...
_SQUARING_PER_POINT_N3 = 1e-8  # ... plus this times (n + 1)**3


def _hypoexpon_check_parameters(scales):
    scales = np.asarray(scales, dtype=float)
    if scales.ndim != 1:
        raise ValueError("Parameter vector `scales` must be one dimensional, "
                         f"but scales.shape = {scales.shape}.")
    if scales.size == 0:
        raise ValueError("Parameter vector `scales` must not be empty.")
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("All entries of `scales` must be positive and finite.")
    return scales


def _closed_form(x, scales, kind):
    r"""Closed form for distinct scales, vectorized over scales and points.

    With :math:`w_i = \prod_{j \ne i} \omega_i / (\omega_i - \omega_j)`,

    .. math::

        f(x) = \sum_i w_i \omega_i^{-1} e^{-x/\omega_i}, \quad
        F(x) = \sum_i w_i (1 - e^{-x/\omega_i}), \quad
        S(x) = \sum_i w_i e^{-x/\omega_i}.

    When scales are close, the weights are huge with alternating signs and
    the sum suffers catastrophic cancellation. For each point an upper bound
    on the relative rounding error of the result is therefore also returned,

        4 n eps * sum_i |term_i| / |sum_i term_i|,

    (the differences of close floats are exact, so each weight is accurate
    to ~2 n eps, and summing n rounded terms loses at most eps per term
    relative to sum |term_i|). The caller recomputes points whose bound is
    too large with the phase-type method. Overflowed weights give non-finite
    values, which the caller must treat as failures as well.
    """
    n = scales.size
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        # ratio[i, j] = omega_i / (omega_i - omega_j); the diagonal is set to 1
        # so that the row product is w_i = prod_{j != i} ratio[i, j]
        ratio = scales[:, None] / (scales[:, None] - scales[None, :])
        np.fill_diagonal(ratio, 1.0)
        weights = np.prod(ratio, axis=1)
        if kind == 'pdf':
            basis = np.exp(-x[None, :] / scales[:, None]) / scales[:, None]
        elif kind == 'cdf':
            # expm1 keeps full relative accuracy for small x
            basis = -np.expm1(-x[None, :] / scales[:, None])
        else:
            basis = np.exp(-x[None, :] / scales[:, None])
        terms = weights[:, None] * basis
        value = np.sum(terms, axis=0)
        rel_err = 4 * n * _EPS * np.sum(np.abs(terms), axis=0) / np.abs(value)
    return value, rel_err


def _series_blocks(lam_sorted, n):
    """Cut sorted lambdas into blocks [start, stop) whose range is at most
    ``10 sqrt(lambda)``, so that one Poisson window around the block's mode
    covers all of them. Returns a list of ``(start, stop, half_width)``."""
    blocks = []
    start = 0
    m = lam_sorted.size
    while start < m:
        lo = lam_sorted[start]
        stop = int(np.searchsorted(lam_sorted, lo + 10 * np.sqrt(lo + 1),
                                   side='right'))
        hi = lam_sorted[stop - 1]
        half = int(12 * np.sqrt(hi + 1) + n + 30)
        max_points = max(1, _SERIES_MAX_BLOCK_ENTRIES // (2 * half + 1))
        stop = min(stop, start + max_points)
        blocks.append((start, stop, half))
        start = stop
    return blocks


def _state_probability(v, kind):
    """Probability of the event described by `kind` given the distribution
    `v` over the n + 1 states of the chain (axis -1): being in the last
    transient phase (pdf, to be multiplied by its rate), having been absorbed
    (cdf) or not having been absorbed yet (sf)."""
    if kind == 'pdf':
        return v[..., -2]
    elif kind == 'cdf':
        return v[..., -1]
    return v[..., :-1].sum(axis=-1)


def _phase_type_series(x, scales, kind):
    r"""Uniformization series, see `_phase_type`.

    .. math::

        g(x) = \sum_k \mathrm{Poisson}(k; \lambda) c_k, \quad
        \lambda = x \mu, \quad c_k = (\alpha P^k)_{\mathrm{col}}, \quad
        P = I + Q / \mu.

    The :math:`c_k` do not depend on x and are obtained by K vector-bidiagonal
    products, O(K n) in total with K ~ lambda_max. The points are sorted by
    lambda and cut into blocks (`_series_blocks`); within a block the Poisson
    weights are generated from the weight at the block's mode k0 by the
    recursions w_{k+1} = w_k lambda / (k + 1) and w_{k-1} = w_k k / lambda,
    as cumulative products over windows of indices, vectorized over both
    points and indices. Each sweep is extended, for the points that need it,
    until a rigorous bound on the remaining tail (0 <= c_k <= 1, geometric
    decay of the weights) is below eps times the value accumulated so far.
    Finally the result is divided by the sum of the weights used (= 1 up to
    the negligible tails), which removes the rounding error of the starting
    weight. Everything is nonnegative, so there is no cancellation.
    """
    # States 0..n-1 are "currently in exponential phase i", state n is
    # absorbing. Uniformization with rate mu = max(rate) turns the
    # continuous-time chain into a discrete one watched at the ticks of a
    # Poisson(mu) clock: at each tick the chain in phase i moves to phase
    # i + 1 with probability move[i] = rate_i / mu and stays otherwise.
    rates = 1.0 / scales
    n = rates.size
    mu = rates.max()
    stay = 1.0 - rates / mu
    move = rates / mu
    lam = x * mu
    out = np.empty(x.size)

    order = np.argsort(lam)
    lam_sorted = lam[order]
    lam_max = lam_sorted[-1]
    # Poisson weights beyond the mode + 40 sigma underflow to exactly 0, so
    # the upward sweep always terminates before K.
    K = int(np.ceil(lam_max + 40 * np.sqrt(lam_max + 1) + n + 50))

    # Step 1: c_k = probability of the event after k ticks, k = 0..K.
    # v is the distribution over the n + 1 states after k ticks; one tick is
    # v <- v P with P bidiagonal, i.e. O(n) per step.
    c = np.empty(K + 1)
    v = np.zeros(n + 1)
    v[0] = 1.0
    for k in range(K + 1):
        c[k] = _state_probability(v, kind)
        nxt = np.empty(n + 1)
        nxt[:n] = v[:n] * stay
        nxt[1:n] += v[:n - 1] * move[:n - 1]
        nxt[n] = v[n] + v[n - 1] * move[n - 1]
        v = nxt

    # Step 2: for each block of points, sum_k Poisson(k; lambda) c_k.
    for start, stop, half in _series_blocks(lam_sorted, n):
        idx = order[start:stop]
        lb = lam_sorted[start:stop]
        k0 = int(lb[0])
        w0 = np.exp(xlogy(k0, lb) - lb - gammaln(k0 + 1))
        acc = w0 * c[k0]
        wsum = w0.copy()

        # upward sweep from k0; `active` are the points whose tail bound is
        # not yet negligible
        active = np.arange(lb.size)
        w_last = w0.copy()
        k_from, k_to = k0 + 1, min(K, k0 + half)
        while k_from <= k_to and active.size:
            ks = np.arange(k_from, k_to + 1)
            W = w_last[active, None] * np.cumprod(lb[active, None] / ks[None, :],
                                                  axis=1)
            acc[active] += W @ c[ks]
            wsum[active] += W.sum(axis=1)
            w_last[active] = W[:, -1]
            # k_to > lambda here, so the remaining weights decay at least
            # geometrically with ratio lambda / (k_to + 1) < 1; with c_k <= 1
            # the remaining sum is at most w_last (k_to + 1) / (k_to + 1 - lam)
            tail = w_last[active] * (k_to + 1) / (k_to + 1 - lb[active])
            active = active[tail > _EPS * acc[active]]
            k_from, k_to = k_to + 1, min(K, k_to + half)

        # downward sweep from k0
        active = np.arange(lb.size)
        w_last = w0.copy()
        k_from, k_to = k0 - 1, max(0, k0 - half)
        while k_from >= k_to and active.size:
            ks = np.arange(k_from, k_to - 1, -1)
            W = w_last[active, None] * np.cumprod((ks[None, :] + 1)
                                                  / lb[active, None], axis=1)
            acc[active] += W @ c[ks]
            wsum[active] += W.sum(axis=1)
            w_last[active] = W[:, -1]
            if k_to == 0:
                break
            # k_to < lambda here: the weights decay with ratio k / lambda < 1
            tail = w_last[active] * lb[active] / (lb[active] - k_to)
            active = active[tail > _EPS * acc[active]]
            k_from, k_to = k_to - 1, max(0, k_to - half)

        # Step 3: normalize by the weights used (sum to 1 up to the tails),
        # cancelling the rounding error of w0.
        out[idx] = acc / wsum
    return out


def _phase_type_squaring(x, scales, kind):
    """Scaling-and-squaring variant of the uniformization, see `_phase_type`.
    O(len(x) (n + log2(x mu)) n^3)."""
    rates = 1.0 / scales
    n = rates.size
    mu = rates.max()
    # B = Q + mu I is nonnegative (bidiagonal: mu - rate_i on the diagonal,
    # rate_i above it, mu for the absorbing state) with row sums mu, so
    # expm(h Q) = exp(-h mu) expm(h B) and the Taylor series of expm(h B) has
    # only nonnegative terms.
    idx = np.arange(n)
    B = np.zeros((n + 1, n + 1))
    B[idx, idx] = mu - rates
    B[idx, idx + 1] = rates
    B[n, n] = mu
    out = np.empty(x.size)

    # Choose s so that h = x / 2^s has h mu <= 1; then the Taylor series of
    # expm(h B) converges in ~n + 20 terms and expm(x Q) = expm(h Q)^(2^s).
    squarings = np.zeros(x.size, dtype=int)
    big = x * mu > 1
    squarings[big] = np.ceil(np.log2(x[big] * mu)).astype(int)

    block = max(1, _SQUARING_MAX_BLOCK_ENTRIES // (n + 1) ** 2)
    for s in np.unique(squarings):
        where = np.flatnonzero(squarings == s)
        for start in range(0, where.size, block):
            sel = where[start:start + block]
            h = x[sel] / 2.0 ** s
            A = h[:, None, None] * B
            # Taylor series, all terms nonnegative; stop once the newest term
            # no longer changes any entry, but not before k = n because the
            # first n powers of a bidiagonal matrix fill in the corner entries
            E = np.broadcast_to(np.eye(n + 1), A.shape).copy()
            term = E.copy()
            for k in range(1, n + 100):
                term = term @ A / k
                E += term
                if k >= n and np.all(term <= _EPS * E):
                    break
            E *= np.exp(-h * mu)[:, None, None]
            for _ in range(s):
                E = E @ E
                # E is stochastic; renormalizing the rows stops the row-sum
                # rounding error from doubling at every squaring
                E /= E.sum(axis=2, keepdims=True)
            out[sel] = _state_probability(E[:, 0, :], kind)
    return out


def _phase_type_plan(lam, n):
    """Choose the threshold on lambda = x mu below which points are handled by
    the series variant (the others by scaling and squaring) by minimizing the
    modelled cost over a few candidate thresholds."""
    lam_sorted = np.sort(lam)
    m = lam_sorted.size
    sq_point = _SQUARING_PER_POINT + _SQUARING_PER_POINT_N3 * (n + 1) ** 3
    per_point = np.cumsum(_SERIES_PER_TERM
                          * (34 * np.sqrt(lam_sorted) + 2 * n + 60))
    best_cost, threshold = m * sq_point, -1.0
    candidates = np.unique(np.minimum(m, np.geomspace(1, m, 24).astype(int)))
    for j in candidates:
        L = lam_sorted[j - 1]
        if L > _SERIES_MAX_LAMBDA:
            break
        n_blocks = len(_series_blocks(lam_sorted[:j], n))
        cost = (_SERIES_FIXED_PER_K * (L + 40 * np.sqrt(L + 1) + n + 50)
                + per_point[j - 1] + _SERIES_PER_BLOCK * n_blocks
                + (m - j) * sq_point)
        if cost < best_cost:
            best_cost, threshold = cost, L
    return threshold


def _phase_type(x, scales, kind):
    r"""Phase-type representation, valid for any scales including repeated
    ones.

    The distribution is the absorption time of the Markov chain
    0 -> 1 -> ... -> n (state i leaves at rate 1/scale_i, state n absorbing)
    with generator Q and initial distribution alpha = (1, 0, ..., 0), so

        pdf(x) = rate_{n-1} expm(x Q)[0, n-1],   cdf(x) = expm(x Q)[0, n].

    expm(x Q) is evaluated by uniformization rather than `scipy.linalg.expm`:
    with mu = max(rate) and P = I + Q / mu (a nonnegative, bidiagonal
    stochastic matrix),

        expm(x Q) = exp(-x mu) sum_k (x mu)^k / k! P^k,

    a sum of nonnegative terms, so no cancellation occurs and every entry is
    accurate to about machine precision, including tiny values for very small
    or very large x. Two variants, with lambda = x mu:

    * `_phase_type_series` sums the series directly: a fixed O(lambda_max n)
      part plus O(sqrt(lambda) + n) per point, i.e. essentially linear in
      both the number of points and the number of scales. Its rounding error
      is ~lambda eps, so it is only used for lambda <= 1e5.
    * `_phase_type_squaring` uses scaling and squaring of expm(h Q) with
      h mu <= 1: O((n + log2 lambda) n^3) per point, accurate for any lambda.

    The points are split between the two by a small cost model.
    """
    n = scales.size
    mu = (1.0 / scales).max()
    lam = x * mu
    out = np.empty(x.size)
    threshold = _phase_type_plan(lam, n)
    use_series = lam <= threshold
    if np.any(use_series):
        out[use_series] = _phase_type_series(x[use_series], scales, kind)
    if not np.all(use_series):
        out[~use_series] = _phase_type_squaring(x[~use_series], scales, kind)
    if kind == 'pdf':
        out /= scales[-1]
    else:
        np.clip(out, 0.0, 1.0, out=out)
    return out


def _hypoexpon_eval(x, scales, kind):
    """pdf/cdf/sf for arbitrary x (any shape): closed form where its error
    bound allows, phase-type otherwise; repeated scales always use the
    phase-type method. Handles x <= 0, x = inf and nan."""
    x = np.asarray(x, dtype=float)
    shape = x.shape
    x = x.ravel()
    out = np.full(x.size, np.nan)
    if kind == 'pdf':
        out[x < 0] = 0.0
        out[x == np.inf] = 0.0
        out[x == 0] = 1.0 / scales[0] if scales.size == 1 else 0.0
    elif kind == 'cdf':
        out[x <= 0] = 0.0
        out[x == np.inf] = 1.0
    else:
        out[x <= 0] = 1.0
        out[x == np.inf] = 0.0
    todo = np.flatnonzero(np.isfinite(x) & (x > 0))
    if todo.size:
        xt = x[todo]
        distinct = np.unique(scales).size == scales.size
        if distinct:
            value, rel_err = _closed_form(xt, scales, kind)
            good = np.isfinite(value) & (value > 0) & (rel_err <= _CLOSED_FORM_RTOL)
            bad = np.flatnonzero(~good)
            if bad.size:
                value[bad] = _phase_type(xt[bad], scales, kind)
        else:
            value = _phase_type(xt, scales, kind)
        out[todo] = value
    return out.reshape(shape)


_hypoexpon_doc_default_callparams = """\
scales : array_like
    Scale parameters :math:`\\omega_1, \\ldots, \\omega_n` (the reciprocals of
    the rates) of the exponential summands. Must be positive and finite; they
    need not be distinct.
"""

_hypoexpon_doc_callparams_note = ""

_hypoexpon_doc_frozen_callparams = ""

_hypoexpon_doc_frozen_callparams_note = """\
See class definition for a detailed description of parameters."""

hypoexpon_docdict_params = {
    '_hypoexpon_doc_default_callparams': _hypoexpon_doc_default_callparams,
    '_hypoexpon_doc_callparams_note': _hypoexpon_doc_callparams_note,
    '_doc_random_state': _doc_random_state
}

hypoexpon_docdict_noparams = {
    '_hypoexpon_doc_default_callparams': _hypoexpon_doc_frozen_callparams,
    '_hypoexpon_doc_callparams_note': _hypoexpon_doc_frozen_callparams_note,
    '_doc_random_state': _doc_random_state
}


class hypoexpon_gen(multi_rv_generic):
    r"""A hypoexponential random variable.

    The hypoexponential (or generalized Erlang) distribution is the
    distribution of the sum of :math:`n` independent exponential random
    variables with scale parameters :math:`\omega_1, \ldots, \omega_n`.

    .. versionadded:: 2.0.0

    Parameters
    ----------
    %(_hypoexpon_doc_default_callparams)s
    %(_doc_random_state)s

    Methods
    -------
    pdf(x, scales)
        Probability density function.
    logpdf(x, scales)
        Log of the probability density function.
    cdf(x, scales)
        Cumulative distribution function.
    sf(x, scales)
        Survival function.
    rvs(scales, size=1, random_state=None)
        Draw random samples.
    mean(scales)
        Mean of the distribution.
    var(scales)
        Variance of the distribution.

    Notes
    -----
    The random variable :math:`X = X_1 + \cdots + X_n`, where the
    :math:`X_i` are independent and exponentially distributed with scale
    :math:`\omega_i` (rate :math:`\lambda_i = 1/\omega_i`), has support
    :math:`x \ge 0`, mean :math:`\sum_i \omega_i` and variance
    :math:`\sum_i \omega_i^2`. For :math:`n = 1` it is the exponential
    distribution and for equal scales it is the Erlang distribution
    (`scipy.stats.erlang`).

    When the scales are distinct the probability density function is [1]_

    .. math::

        f(x) = \sum_{i=1}^n \frac{e^{-x/\omega_i}}{\omega_i}
               \prod_{j \ne i} \frac{\omega_i}{\omega_i - \omega_j},
               \qquad x \ge 0.

    The products in this formula are large and of alternating sign when the
    scales are close, so the sum suffers from catastrophic cancellation in
    floating point arithmetic; for repeated scales the formula does not
    apply at all. `hypoexpon` evaluates the closed form only where a bound
    on its rounding error shows the result to be accurate, and otherwise
    uses the phase-type representation [1]_ [2]_: :math:`X` is the absorption
    time of a Markov chain that visits the states :math:`1, \ldots, n` in
    order, leaving state :math:`i` at rate :math:`\lambda_i`. The density and
    distribution function are entries of the matrix exponential
    :math:`e^{xQ}` of the generator, which is evaluated by uniformization
    [3]_ as a sum of nonnegative terms, so that the result is accurate to
    nearly machine precision for any scales, including repeated and nearly
    equal ones.

    References
    ----------
    .. [1] "Hypoexponential distribution", Wikipedia,
           https://en.wikipedia.org/wiki/Hypoexponential_distribution
    .. [2] M. Bladt and B. F. Nielsen, "Matrix-Exponential Distributions in
           Applied Probability", Springer, 2017.
    .. [3] W. J. Stewart, "Probability, Markov Chains, Queues, and
           Simulation", Princeton University Press, 2009, Section 10.4
           (uniformization).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import hypoexpon
    >>> scales = [1., 2., 4.]

    Evaluate the density and the distribution function:

    >>> hypoexpon.pdf(3., scales)
    0.10837656446820126
    >>> hypoexpon.cdf(3., scales)
    0.17002049019819898

    The mean and variance are the sums of the scales and of their squares:

    >>> hypoexpon.mean(scales), hypoexpon.var(scales)
    (7.0, 21.0)

    Repeated or nearly equal scales are supported; with equal scales the
    distribution is the Erlang distribution:

    >>> from scipy.stats import erlang
    >>> hypoexpon.pdf(2.5, [1., 1., 1.])
    0.25651562069968376
    >>> erlang.pdf(2.5, 3)
    0.25651562069968376
    >>> hypoexpon.pdf(1., [1., 1.00000001, 1.00000002])
    0.18393971690692681

    Draw random samples:

    >>> rng = np.random.default_rng()
    >>> x = hypoexpon.rvs(scales, size=5, random_state=rng)
    >>> x.shape
    (5,)

    Alternatively, the object may be called (as a function) to fix the
    scales, returning a "frozen" hypoexponential random variable:

    >>> rv = hypoexpon(scales)
    >>> rv.cdf(3.)
    0.17002049019819898

    """

    def __init__(self, seed=None):
        super().__init__(seed)
        self.__doc__ = doccer.docformat(self.__doc__, hypoexpon_docdict_params)

    def __call__(self, scales, seed=None):
        return hypoexpon_frozen(scales, seed=seed)

    def pdf(self, x, scales):
        """Probability density function of the hypoexponential distribution.

        Parameters
        ----------
        x : array_like
            Quantiles.
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        pdf : ndarray or scalar
            Probability density function evaluated at `x`.

        """
        scales = _hypoexpon_check_parameters(scales)
        return _squeeze_output(_hypoexpon_eval(x, scales, 'pdf'))

    def logpdf(self, x, scales):
        """Log of the probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles.
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        logpdf : ndarray or scalar
            Log of the probability density function evaluated at `x`.

        """
        scales = _hypoexpon_check_parameters(scales)
        with np.errstate(divide='ignore'):
            out = np.log(_hypoexpon_eval(x, scales, 'pdf'))
        return _squeeze_output(out)

    def cdf(self, x, scales):
        """Cumulative distribution function.

        Parameters
        ----------
        x : array_like
            Quantiles.
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        cdf : ndarray or scalar
            Cumulative distribution function evaluated at `x`.

        """
        scales = _hypoexpon_check_parameters(scales)
        return _squeeze_output(_hypoexpon_eval(x, scales, 'cdf'))

    def sf(self, x, scales):
        """Survival function (complement of the cdf).

        Parameters
        ----------
        x : array_like
            Quantiles.
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        sf : ndarray or scalar
            Survival function evaluated at `x`.

        """
        scales = _hypoexpon_check_parameters(scales)
        return _squeeze_output(_hypoexpon_eval(x, scales, 'sf'))

    def mean(self, scales):
        """Mean of the hypoexponential distribution.

        Parameters
        ----------
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        mean : float
            The mean, ``sum(scales)``.

        """
        scales = _hypoexpon_check_parameters(scales)
        return float(np.sum(scales))

    def var(self, scales):
        """Variance of the hypoexponential distribution.

        Parameters
        ----------
        %(_hypoexpon_doc_default_callparams)s

        Returns
        -------
        var : float
            The variance, ``sum(scales**2)``.

        """
        scales = _hypoexpon_check_parameters(scales)
        return float(np.sum(scales ** 2))

    def rvs(self, scales, size=1, random_state=None):
        """Draw random samples from a hypoexponential distribution.

        Parameters
        ----------
        %(_hypoexpon_doc_default_callparams)s
        size : int or tuple of ints, optional
            Shape of the output. Default is 1.
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray
            Random variates of shape `size`.

        """
        scales = _hypoexpon_check_parameters(scales)
        random_state = self._get_random_state(random_state)
        size = (size,) if np.isscalar(size) else tuple(size)
        # a sum of independent exponentials with the given scales
        e = random_state.standard_exponential(size=size + (scales.size,))
        return e @ scales


hypoexpon = hypoexpon_gen()


class hypoexpon_frozen(multi_rv_frozen):
    __class_getitem__ = None  # pyrefly:ignore[bad-assignment]

    def __init__(self, scales, seed=None):
        self.scales = _hypoexpon_check_parameters(scales)
        self._dist = hypoexpon_gen(seed)

    def pdf(self, x):
        return self._dist.pdf(x, self.scales)

    def logpdf(self, x):
        return self._dist.logpdf(x, self.scales)

    def cdf(self, x):
        return self._dist.cdf(x, self.scales)

    def sf(self, x):
        return self._dist.sf(x, self.scales)

    def mean(self):
        return self._dist.mean(self.scales)

    def var(self):
        return self._dist.var(self.scales)

    def rvs(self, size=1, random_state=None):
        return self._dist.rvs(self.scales, size, random_state)


# Set frozen generator docstrings from corresponding docstrings in
# hypoexpon_gen and fill in default strings in class docstrings
for name in ['pdf', 'logpdf', 'cdf', 'sf', 'mean', 'var', 'rvs']:
    method = hypoexpon_gen.__dict__[name]
    method_frozen = hypoexpon_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, hypoexpon_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, hypoexpon_docdict_params)
