import warnings
from scipy.stats.sampling import RatioUniforms

def rvs_ratio_uniforms(pdf, umax, vmin, vmax, size=1, c=0, random_state=None):
    """
    Generate random samples from a probability density function using the
    ratio-of-uniforms method.

    Parameters
    ----------
    pdf : callable
        A function with signature `pdf(x)` that is proportional to the
        probability density function of the distribution.
    umax : float
        The upper bound of the bounding rectangle in the u-direction.
    vmin : float
        The lower bound of the bounding rectangle in the v-direction.
    vmax : float
        The upper bound of the bounding rectangle in the v-direction.
    size : int or tuple of ints, optional
        Defining number of random variates (default is 1).
    c : float, optional.
        Shift parameter of ratio-of-uniforms method, see Notes. Default is 0.
    random_state : {None, int, `numpy.random.Generator`,
                    `numpy.random.RandomState`}, optional

        If `seed` is None (or `np.random`), the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``Generator`` or ``RandomState`` instance then
        that instance is used.

    Returns
    -------
    rvs : ndarray
        The random variates distributed according to the probability
        distribution defined by the pdf.

    Notes
    -----
    Please refer to `scipy.stats.sampling.RatioUniforms` for the documentation
    and use ``RatioUniforms`` instead of this function which is deprecated.
    """
    warnings.warn("Please use `RatioUniforms` from the "
                  "`scipy.stats.sampling` namespace. The "
                  "`scipy.stats.rvs_ratio_uniforms` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)
    gen = RatioUniforms(pdf, umax, vmin, vmax, c, random_state)
    return gen.rvs(size)
