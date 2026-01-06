import numpy as np
from scipy.optimize import isotonic_regression

class grenander:
    """Grenander estimator of a nonincreasing density.

    The Grenander estimator is the nonparametric maximum likelihood estimator 
    of a nonincreasing density with a known lower bound on the support.

    Parameters
    ----------
    dataset : array_like
        1D sample of observations. Must lie strictly above `support_min`.
        Duplicates are not allowed (ties create zero-width intervals).
    support_min : float, default=0.0
        Known lower bound of the support.
    assume_sorted : bool, default=False
        If True, assume `dataset` is already sorted ascending.

    Attributes
    ----------
    fitted : ndarray, shape (n,)
        Estimated density values corresponding to observations in `dataset`.
    knots : ndarray, shape (B+1,)
        Bin edges / knot locations of the least concave majorant (LCM) CDF.
    heights : ndarray, shape (B+1,)
        CDF values at `knots` (piecewise-linear between knots).
    slopes : ndarray, shape (B,)
        Density levels on intervals [knots[j], knots[j+1]).
    support_min : float
        Lower support bound.

    See Also
    --------
    scipy.optimize.isotonic_regression : Isotonic regression 
    scipy.stats.ecdf : Empirical cumulative distribution function
    scipy.stats.gaussian_kde : Gaussian kernel density estimation

    Notes
    -----
    The Grenander estimator [2]_ is the nonparametric maximum likelihood 
    estimator of a nonincreasing density on a known, lower-bounded support. 
    Equivalently, it solves

    .. math::

        \\text{maximize} \\frac{1}{n}\\sum_{i=1}^n \\log f(x_i)

    subject to

    .. math::
        
        f : [0,\\infty) \\to [0,\\infty) \\text{ nonincreasing, and } 
        \\int_0^\\infty f(x)\\,dx = 1.

    The estimator may be characterized geometrically as the left derivative of 
    the least concave majorant (LCM) of the empirical cumulative distribution 
    function [3]_. As a result, the fitted density is a piecewise-constant 
    function with data-adaptive bin widths, as well as a decreasing histogram 
    whose breakpoints are determined by the pool-adjacent-violators algorithm.
    Under mild regularity conditions, the Grenander estimator is consistent and
    converges at the nonstandard rate :math:`n^{-1/3}` [1]_ [4]_. 
    
    References 
    ---------- 
    .. [1] L. Birgé. "Estimating a density under order restrictions: 
           Nonasymptotic minimax risk", The Annals of Statistics, Vol. 15, pp. 
           995-1012, 1987.

    .. [2] U. Grenander, "On the theory of mortality measurement: part II", 
           Scandinavian Actuarial Journal, Vol. 39, pp. 125-153, 1956. 
           
    .. [3] P. Groeneboom and G. Jongbloed, "Nonparametric estimation under 
           shape constraints", Cambridge University Press, New York, 2014. 
           
    .. [4] B.L.S.P. Rao, "Estimation of a unimodal density", Sankhyā: The 
           Indian Journal of Statistics, Series A, Vol. 31, pp. 23-36, 1969.

    Examples
    --------
    Estimate a decreasing density from exponential data.

    >>> import numpy as np
    >>> from scipy.stats import expon, grenander
    >>> rng = np.random.default_rng(1234567891)
    >>> x = rng.exponential(size=100)
    >>> g = grenander(x)

    The estimator provides the fitted density values at the observations:

    >>> g.fitted.shape
    (100,)

    Evaluate the PDF and CDF at specific points:

    >>> g.pdf([0.5, 1.0, 2.0])
    array([0.52828842, 0.23170026, 0.21893189])
    >>> g.cdf([0.5, 1.0, 2.0])
    array([0.37589365, 0.60242482, 0.83153771])

    The estimator is piecewise constant with data-adaptive breakpoints:

    >>> len(g.slopes)  # number of constant pieces
    15
    >>> g.knots  # breakpoints of the histogram
    array([0.        , 0.00903572, 0.02683709, 0.05754583, 0.46991483,
           0.73492159, 0.77473536, 0.98953465, 1.11901228, 1.9016237 ,
           2.08432892, 2.29638573, 2.68628217, 3.38105853, 4.67484783,
           5.93619219])

    Convert to a distribution object for additional functionality:

    >>> dist = stats.make_distribution(g)()
    >>> samples = dist.sample(500)
    >>> median = dist.icdf(0.5)

    Visualize the estimated density alongside the true exponential density:

    >>> import matplotlib.pyplot as plt
    >>> xs = np.linspace(0, 6, 500)
    >>> plt.plot(xs, stats.expon.pdf(xs), 'k-', label='True density')
    >>> plt.plot(xs, g.pdf(xs), drawstyle='steps-post', label='Grenander')
    >>> plt.xlabel('x')
    >>> plt.ylabel('Density')
    >>> plt.legend()
    >>> plt.show()
    """
    __make_distribution_version__ = "1.16.0"
    parameters = tuple()

    def __init__(self,
                 dataset,
                 support_min=0.0,
                 assume_sorted=False):
        self.support_min = support_min
        params = self._fit(dataset, assume_sorted)
        self.__dict__.update(params)

    @property
    def support(self):
        """Support dict for make_distribution compatibility."""
        return {
            'endpoints': (self.support_min, self.knots[-1]),
            'inclusive': (True, True)
        }

    def pdf(self, x):
        """Evaluate the estimated density at `x`."""
        x = np.asarray(x, dtype=float)
        indices = np.searchsorted(self.knots, x, side='right') - 1
        indices = np.clip(indices, 0, len(self.slopes) - 1)
        
        in_support = (self.support_min <= x) & (x <= self.knots[-1])
        
        return np.where(in_support, self.slopes[indices], 0.0)
    
    __call__ = pdf

    def cdf(self, x):
        """Evaluate the estimated CDF at `x`."""
        x = np.asarray(x, dtype=float)
        return np.interp(x, self.knots, self.heights, 
                         left=0.0, right=1.0)
    
    def icdf(self, q):
        """Percent point function (inverse CDF)."""
        q = np.asarray(q, dtype=float)
        return np.interp(q, self.heights, self.knots,
                         left=self.support_min,
                         right=self.knots[-1])

    def _fit(self, dataset, assume_sorted=False):
        """Internal fitting method."""
        x = np.asarray(dataset, dtype=float)
        n = x.size
        if n == 0:
            raise ValueError("`dataset` has no observations.")
        if np.any(x <= self.support_min):
            raise ValueError(
                f"All observations in `dataset` must be strictly greater than "
                f"support_min={self.support_min}"
            )
        if not assume_sorted:
            I = x.argsort()
            x = x[I]
        if np.any(x[1:] == x[:-1]):
            raise ValueError(
                "`dataset` contains duplicate values. The Grenander estimator "
                "requires strictly distinct observations."
            )

        # Collect gaps in order statistics
        W = np.concat([[self.support_min], x])
        w = np.diff(W)

        # Raw ECDF slopes on each interval are 1/(n*w); LCM slopes are the
        # isotonic (nonincreasing) regression of these with weights w
        raw = np.ones(n)/(n*w)
        res = isotonic_regression(raw, 
                                weights=w.copy(), 
                                increasing=False)
        fitted = res.x
        blocks  = res.blocks

        params = dict()
        params['knots']   = W[blocks]
        params['heights'] = blocks/n
        params['slopes']  = fitted[blocks[:-1]]

        # Return fitted values in original order
        if not assume_sorted:
            fitted = fitted[I.argsort()]

        params['fitted'] = fitted

        return params