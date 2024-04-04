# Docstrings for generated ufuncs
#
# The syntax is designed to look like the function add_newdoc is being
# called from numpy.lib, but in this file add_newdoc puts the
# docstrings in a dictionary. This dictionary is used in
# _generate_pyx.py to generate the docstrings for the ufuncs in
# scipy.special at the C level when the ufuncs are created at compile
# time.

docdict: dict[str, str] = {}


def get(name):
    return docdict.get(name)


def add_newdoc(name, doc):
    docdict[name] = doc

add_newdoc(
    "_beta_pdf",
    """
    _beta_pdf(x, a, b)

    Probability density function of beta distribution.

    Parameters
    ----------
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration
    a, b : array_like
           Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_invgauss_ppf",
    """
    _invgauss_ppf(x, mu)

    Percent point function of inverse gaussian distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    mu : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_invgauss_isf",
    """
    _invgauss_isf(x, mu, s)

    Inverse survival function of inverse gaussian distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    mu : array_like
        Positive, real-valued parameters
    s : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_pdf",
    """
    _ncx2_pdf(x, k, l)

    Probability density function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_cdf",
    """
    _ncx2_cdf(x, k, l)

    Cumulative density function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_ppf",
    """
    _ncx2_ppf(x, k, l)

    Percent point function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_sf",
    """
    _ncx2_sf(x, k, l)

    Survival function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_isf",
    """
    _ncx2_isf(x, k, l)

    Inverse survival function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_pdf",
    """
    _ncf_pdf(x, v1, v2, l)

    Probability density function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_cdf",
    """
    _ncf_cdf(x, v1, v2, l)

    Cumulative density function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_ppf",
    """
    _ncf_ppf(x, v1, v2, l)

    Percent point function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_sf",
    """
    _ncf_sf(x, v1, v2, l)

    Survival function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_isf",
    """
    _ncf_isf(x, v1, v2, l)

    Inverse surivial function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_mean",
    """
    _ncf_mean(v1, v2, l)

    Mean of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_variance",
    """
    _ncf_variance(v1, v2, l)

    Variance of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_skewness",
    """
    _ncf_skewness(v1, v2, l)

    Skewness of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_kurtosis_excess",
    """
    _ncf_kurtosis_excess(v1, v2, l)

    Kurtosis excess of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_cdf",
    """
    _nct_cdf(x, v, l)

    Cumulative density function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_ppf",
    """
    _nct_ppf(x, v, l)

    Percent point function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_sf",
    """
    _nct_sf(x, v, l)

    Survival function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_isf",
    """
    _nct_isf(x, v, l)

    Inverse surivial function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_mean",
    """
    _nct_mean(v, l)

    Mean of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_variance",
    """
    _nct_variance(v, l)

    Variance of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_skewness",
    """
    _nct_skewness(v, l)

    Skewness of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_kurtosis_excess",
    """
    _nct_kurtosis_excess(v, l)

    Kurtosis excess of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_cdf",
    """
    _skewnorm_cdf(x, l, sc, sh)

    Cumulative density function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_ppf",
    """
    _skewnorm_ppf(x, l, sc, sh)

    Percent point function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_isf",
    """
    _skewnorm_isf(x, l, sc, sh)

    Inverse surivial function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_pmf",
    """
    _binom_pmf(x, n, p)

    Probability mass function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_liks
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_cdf",
    """
    _binom_cdf(x, n, p)

    Cumulative density function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_liks
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_ppf",
    """
    _binom_ppf(x, n, p)

    Percent point function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_liks
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_sf",
    """
    _binom_sf(x, n, p)

    Survival function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_liks
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_isf",
    """
    _binom_isf(x, n, p)

    Inverse survival function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_liks
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)
