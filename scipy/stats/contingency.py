"""Some functions for working with contingency tables (i.e. cross tabulations).
"""


from __future__ import division, print_function, absolute_import

from functools import reduce
import numpy as np
from .stats import power_divergence, chisquare
import math
from warnings import warn

__all__ = ['margins', 'expected_freq', 'chi2_contingency', 'association_test']


def margins(a):
    """Return a list of the marginal sums of the array `a`.

    Parameters
    ----------
    a : ndarray
        The array for which to compute the marginal sums.

    Returns
    -------
    margsums : list of ndarrays
        A list of length `a.ndim`.  `margsums[k]` is the result
        of summing `a` over all axes except `k`; it has the same
        number of dimensions as `a`, but the length of each axis
        except axis `k` will be 1.

    Examples
    --------
    >>> a = np.arange(12).reshape(2, 6)
    >>> a
    array([[ 0,  1,  2,  3,  4,  5],
           [ 6,  7,  8,  9, 10, 11]])
    >>> m0, m1 = margins(a)
    >>> m0
    array([[15],
           [51]])
    >>> m1
    array([[ 6,  8, 10, 12, 14, 16]])

    >>> b = np.arange(24).reshape(2,3,4)
    >>> m0, m1, m2 = margins(b)
    >>> m0
    array([[[ 66]],
           [[210]]])
    >>> m1
    array([[[ 60],
            [ 92],
            [124]]])
    >>> m2
    array([[[60, 66, 72, 78]]])
    """
    margsums = []
    ranged = list(range(a.ndim))
    for k in ranged:
        marg = np.apply_over_axes(np.sum, a, [j for j in ranged if j != k])
        margsums.append(marg)
    return margsums


def expected_freq(observed):
    """
    Compute the expected frequencies from a contingency table.

    Given an n-dimensional contingency table of observed frequencies,
    compute the expected frequencies for the table based on the marginal
    sums under the assumption that the groups associated with each
    dimension are independent.

    Parameters
    ----------
    observed : array_like
        The table of observed frequencies.  (While this function can handle
        a 1-D array, that case is trivial.  Generally `observed` is at
        least 2-D.)

    Returns
    -------
    expected : ndarray of float64
        The expected frequencies, based on the marginal sums of the table.
        Same shape as `observed`.

    Examples
    --------
    >>> observed = np.array([[10, 10, 20],[20, 20, 20]])
    >>> from scipy.stats import expected_freq
    >>> expected_freq(observed)
    array([[ 12.,  12.,  16.],
           [ 18.,  18.,  24.]])

    """
    # Typically `observed` is an integer array. If `observed` has a large
    # number of dimensions or holds large values, some of the following
    # computations may overflow, so we first switch to floating point.
    observed = np.asarray(observed, dtype=np.float64)

    # Create a list of the marginal sums.
    margsums = margins(observed)

    # Create the array of expected frequencies.  The shapes of the
    # marginal sums returned by apply_over_axes() are just what we
    # need for broadcasting in the following product.
    d = observed.ndim
    expected = reduce(np.multiply, margsums) / observed.sum() ** (d - 1)
    return expected


def chi2_contingency(observed, correction=True, lambda_=None):
    """Chi-square test of independence of variables in a contingency table.

    This function computes the chi-square statistic and p-value for the
    hypothesis test of independence of the observed frequencies in the
    contingency table [1]_ `observed`.  The expected frequencies are computed
    based on the marginal sums under the assumption of independence; see
    `scipy.stats.contingency.expected_freq`.  The number of degrees of
    freedom is (expressed using numpy functions and attributes)::

        dof = observed.size - sum(observed.shape) + observed.ndim - 1


    Parameters
    ----------
    observed : array_like
        The contingency table. The table contains the observed frequencies
        (i.e. number of occurrences) in each category.  In the two-dimensional
        case, the table is often described as an "R x C table".
    correction : bool, optional
        If True, *and* the degrees of freedom is 1, apply Yates' correction
        for continuity.  The effect of the correction is to adjust each
        observed value by 0.5 towards the corresponding expected value.
    lambda_ : float or str, optional.
        By default, the statistic computed in this test is Pearson's
        chi-squared statistic [2]_.  `lambda_` allows a statistic from the
        Cressie-Read power divergence family [3]_ to be used instead.  See
        `power_divergence` for details.

    Returns
    -------
    chi2 : float
        The test statistic.
    p : float
        The p-value of the test
    dof : int
        Degrees of freedom
    expected : ndarray, same shape as `observed`
        The expected frequencies, based on the marginal sums of the table.

    See Also
    --------
    contingency.expected_freq
    fisher_exact
    chisquare
    power_divergence

    Notes
    -----
    An often quoted guideline for the validity of this calculation is that
    the test should be used only if the observed and expected frequencies
    in each cell are at least 5.

    This is a test for the independence of different categories of a
    population. The test is only meaningful when the dimension of
    `observed` is two or more.  Applying the test to a one-dimensional
    table will always result in `expected` equal to `observed` and a
    chi-square statistic equal to 0.

    This function does not handle masked arrays, because the calculation
    does not make sense with missing values.

    Like stats.chisquare, this function computes a chi-square statistic;
    the convenience this function provides is to figure out the expected
    frequencies and degrees of freedom from the given contingency table.
    If these were already known, and if the Yates' correction was not
    required, one could use stats.chisquare.  That is, if one calls::

        chi2, p, dof, ex = chi2_contingency(obs, correction=False)

    then the following is true::

        (chi2, p) == stats.chisquare(obs.ravel(), f_exp=ex.ravel(),
                                     ddof=obs.size - 1 - dof)

    The `lambda_` argument was added in version 0.13.0 of scipy.

    References
    ----------
    .. [1] "Contingency table", http://en.wikipedia.org/wiki/Contingency_table
    .. [2] "Pearson's chi-squared test",
           http://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
    .. [3] Cressie, N. and Read, T. R. C., "Multinomial Goodness-of-Fit
           Tests", J. Royal Stat. Soc. Series B, Vol. 46, No. 3 (1984),
           pp. 440-464.

    Examples
    --------
    A two-way example (2 x 3):

    >>> from scipy.stats import chi2_contingency
    >>> obs = np.array([[10, 10, 20], [20, 20, 20]])
    >>> chi2_contingency(obs)
    (2.7777777777777777,
     0.24935220877729619,
     2,
     array([[ 12.,  12.,  16.],
            [ 18.,  18.,  24.]]))

    Perform the test using the log-likelihood ratio (i.e. the "G-test")
    instead of Pearson's chi-squared statistic.

    >>> g, p, dof, expctd = chi2_contingency(obs, lambda_="log-likelihood")
    >>> g, p
    (2.7688587616781319, 0.25046668010954165)

    A four-way example (2 x 2 x 2 x 2):

    >>> obs = np.array(
    ...     [[[[12, 17],
    ...        [11, 16]],
    ...       [[11, 12],
    ...        [15, 16]]],
    ...      [[[23, 15],
    ...        [30, 22]],
    ...       [[14, 17],
    ...        [15, 16]]]])
    >>> chi2_contingency(obs)
    (8.7584514426741897,
     0.64417725029295503,
     11,
     array([[[[ 14.15462386,  14.15462386],
              [ 16.49423111,  16.49423111]],
             [[ 11.2461395 ,  11.2461395 ],
              [ 13.10500554,  13.10500554]]],
            [[[ 19.5591166 ,  19.5591166 ],
              [ 22.79202844,  22.79202844]],
             [[ 15.54012004,  15.54012004],
              [ 18.10873492,  18.10873492]]]]))
    """
    observed = np.asarray(observed)
    if np.any(observed < 0):
        raise ValueError("All values in `observed` must be nonnegative.")
    if observed.size == 0:
        raise ValueError("No data; `observed` has size 0.")

    expected = expected_freq(observed)
    if np.any(expected == 0):
        # Include one of the positions where expected is zero in
        # the exception message.
        zeropos = list(zip(*np.where(expected == 0)))[0]
        raise ValueError("The internally computed table of expected "
                         "frequencies has a zero element at %s." % (zeropos,))

    # The degrees of freedom
    dof = expected.size - sum(expected.shape) + expected.ndim - 1

    if dof == 0:
        # Degenerate case; this occurs when `observed` is 1D (or, more
        # generally, when it has only one nontrivial dimension).  In this
        # case, we also have observed == expected, so chi2 is 0.
        chi2 = 0.0
        p = 1.0
    else:
        if dof == 1 and correction:
            # Adjust `observed` according to Yates' correction for continuity.
            observed = observed + 0.5 * np.sign(expected - observed)

        chi2, p = power_divergence(observed, expected,
                                   ddof=observed.size - 1 - dof, axis=None,
                                   lambda_=lambda_)

    return chi2, p, dof, expected


def _reshape_array(arr):
    """Converts Nd array to 2d array.

    Parameters
    ----------

    arr : array_like
        Observed values in n-way array

    Returns
    -------
    arrx : array_like 
        Observed values in reshaped to n by n array.
    nrows : int
        Number of rows in reshaped array.
    ncols : int
        Number of columns in reshaped array.

    Notes
    -----
    Required if the array has more that 2 dimensions in order to 
    calculate the number of rows and columns
    """
    arr = np.array(arr)
    shape = arr.shape
    if len(shape) > 2:
        if shape[0] == 1:
            nrows = shape[1]
            ncols = np.prod([i for i in shape[2:]])
        elif shape[0] == 0:
            raise IndexError("Invalid Array Shape: Empty Array Found")
        else:
            nrows = shape[0]
            ncols = np.prod([i for i in shape[1:]])
        arrx = arr.reshape(nrows, ncols)
        return arrx, nrows, ncols

    elif len(shape) == 2:
        nrows = shape[0]
        ncols = shape[1]
        arrx = arr
        return arrx, nrows, ncols

    else:
        raise IndexError("Invalid Array Shape: Vectors or 1D arrays are not supported")


def _association_bias_correction(phi_squared, n_rows, n_cols, n_obs):
    """Corrects bias in estimated value of phi squared derived from nxn contingency 
    tables.

    Parameters
    ----------
    phi_squared : float
        Empirical Phi Squared Value
    n_rows : int
        Unadjusted number of rows
    n_cols : int
        Unadjusted number of columns
    n_obs : int
        Total number of observations
    
    Returns
    -------
    adj_phisq : float
        Unbiased Phi Squared value
    nrows_hat : float
        Unbiased number of rows
    ncols_hat : float
        Unbiased number of columns


    References
    ----------
    .. [1] Bergsma, Wicher, "A bias-correction for Cramer's V and Tschuprow's T", 
           London School of Econ. and Pol. Sci., pp. 2-4. http://stats.lse.ac.uk/bergsma/pdf/cramerV3.pdf
    .. [2] https://en.wikipedia.org/wiki/Cramer's_V
    

    Notes
    ------
    Improves accuracy of estimators with tables > 2x2 and smaller sample sizes [1]
    
    Cramer's V can be a heavily biased estimator of its population counterpart 
    and will tend to overestimate the strength of association. The adjusted statistic
    estimates the same population quantity as the original statistic but with typically 
    much smaller mean squared error. [2]
    """
    phi_squared = phi_squared - (((n_rows - 1) * (n_cols - 1)) / (n_obs - 1))
    adj_phisq = max(0., phi_squared)
    nrows_hat = float(n_rows - (((n_rows - 1) ** 2) / (n_obs - 1)))
    ncols_hat = float(n_cols - (((n_cols - 1) ** 2) / (n_obs - 1)))
    return adj_phisq, nrows_hat, ncols_hat


def association_test(stat, observed, chi2_stat=None, correct_bias=True):
    """Calculates degree of association between variables that are nominal or greater. Allows for 
    specification of one of three related methods, Tschuprow's T, Cramer's V, and phi.

    Parameters
    ----------
    stat : {"V", "T", "phi"}
        The association test statistic. 
    observed : iterable object
        The contingency table. 
    chi2_stat : float or None, optional (default = None)
        The chi squared statistic. If equal to None, chi squared value will be automatically calculated using
        scipy.contingency.chi2_contingency() method.
    correct_bias : boolean, optional (default = True)
        If True, bias correction will be applied to phi as per Bergsma (2013). 
    
    Returns
    -------
    value : float
        Value of the test statistic
        


    References
    ----------
    .. [1] "Tschuprow's T",
           https://en.wikipedia.org/wiki/Tschuprow's_T
    .. [2] Bergsma, Wicher, "A bias-correction for Cramer's V and Tschuprow's T", 
           London School of Econ. and Pol. Sci., pp. 5-7. 
           http://stats.lse.ac.uk/bergsma/pdf/cramerV3.pdf
    .. [3] Tschuprow, A. A. (1939) Principles of the Mathematical Theory of Correlation; 
           translated by M. Kantorowitsch. W. Hodge & Co.
    .. [4] "Cramer's V", https://en.wikipedia.org/wiki/Cramer's_V
    .. [5] "Nominal Association: Phi and Cramer's V",
           http://www.people.vcu.edu/~pdattalo/702SuppRead/MeasAssoc/NominalAssoc.html


    Examples
    --------
    
    2-way Example
    
    >>> obs = [[100, 150], [203, 322], [42, 7], [32, 21]]
    
    Cramer's V with bias correction
    >>> association_test(stat="V", observed=obs)
    0.46927187061981274
    
    Cramer's V without bias correction
    >>> association_test(stat="V", observed=obs, correct_bias=False)
    0.4726408338900912
    
    Tschuprow's T with bias correction
    >>> association_test(stat="T", observed=obs)
    0.35677355915423664
    
    Tschuprow's T without bias correction
    >>> association_test(stat="T", observed=obs, correct_bias=False)
    0.3591293720858179
    
    Phi with bias correction
    >>> association_test(stat="phi", observed=obs)
    0.46900394489393604
    
    Phi without bias correction
    >>> association_test(stat="phi", observed=obs, correct_bias=False)
    0.4726408338900912
    
    4-way Example
    
    >>> obs = [[[[56, 23],
    ...          [21, 45]],
    ...         [[13, 16],
    ...          [76, 99]],
    ...         [[21, 22],
    ...          [41, 44]]]]
    
    Cramer's V with bias correction
    >>> association_test(stat="V", observed=obs)
    0.17427841264890478
    
    Cramer's V without bias correction
    >>> association_test(stat="V", observed=obs, correct_bias=False)
    0.1911749686107999
    
    Tschuprow's T with bias correction
    >>> association_test(stat="T", observed=obs)
    0.15756149073389059
    
    Tschuprow's T without bias correction
    >>> association_test(stat="T", observed=obs, correct_bias=False)
    0.172746084676768
    
    Phi with bias correction
    >>> association_test(stat="phi", observed=obs)
    0.24594856216264802
    
    Phi without bias correction
    >>> association_test(stat="phi", observed=obs, correct_bias=False)
    0.27036223339564397

    Notes
    ------
    Cramer's V and Tschuprow's T measure degree to which two variables are related, or the level of 
    their association. This differs from correlation, although many often mistakenly consider them equivalent.
    Correlation measures in what way two variables are related, whereas, association measures
    how related the variables are. As such, association does not subsume independent variables, and is 
    rather a test of independence. Where a value of 1.0 = perfect association or dependent variables, and 
    0.0 = no association or entirely independent variables.
    
    Both the Cramer's V and Tschuprow's T are extensions of the phi coefficient. Moreover, due 
    to the close relationship between the Cramer's V and Tschuprow's T the returned values can often
    be similar or even equivalent. They are likely to diverge more as the array shape diverges from a 2x2.
    As is seen in the examples above. 
    
    """
    arr, n_rows, n_cols = _reshape_array(arr=observed)
    flat_array = arr.flatten()

    if len([i for i in flat_array if type(i) == np.float64]) > 0:
        raise TypeError("Array must be counts not frequencies.")
    elif len([i for i in flat_array if i < 5]) > 0:
        warn("The Chi-Squared Statistic is invalid if any counts are < 5")
    else:
        pass

    n_obs = flat_array.sum(0, dtype=np.int64)

    if chi2_stat is None:
        phi2 = float(tuple(chisquare(f_obs=arr))[0][0]) / n_obs
    elif type(chi2_stat) is float or type(chi2_stat) is int:
        phi2 = float(chi2_stat) / n_obs
    else:
        raise TypeError("Invalid chi2_stat value")

    if correct_bias is True:
        phi_sq, nrows, ncols = _association_bias_correction(phi_squared=phi2, n_rows=n_rows, n_cols=n_cols, n_obs=n_obs)
    elif correct_bias is False:
        phi_sq = phi2
        nrows = n_rows
        ncols = n_cols
    else:
        raise TypeError("invalid argument type: 'correct_bias' must be boolean")

    if stat.lower() == "v":
        value = math.sqrt(phi_sq / min(ncols - 1, nrows - 1))
    elif stat.lower() == "t":
        value = math.sqrt(phi_sq / math.sqrt((nrows - 1) * (ncols - 1)))
    elif stat.lower() == "phi":
        value = math.sqrt(phi_sq)
    else:
        raise ValueError("Invalid argument value: 'stat' must be t, v or phi")

    return value
