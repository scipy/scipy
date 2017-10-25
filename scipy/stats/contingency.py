"""Some functions for working with contingency tables (i.e. cross tabulations).
"""


from __future__ import division, print_function, absolute_import

from functools import reduce
import numpy as np
from .stats import power_divergence
import math


__all__ = ['margins', 'expected_freq', 'chi2_contingency', 'associationTests']


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


class associationTests(object):
    """
    This class contains 2 related measures of association for nominal data, 
    Cramer's V and Tschuprow's T. Each measure calculates 
    """

    def __init__(self, observed, chi2_stat=None, **kwargs):
        """
        Both the Cramer's V and Tschuprow's T are empirically determined using the chi-squared 
        statistic 

        :param chi2_stat: float, required (depending on contingency table dims)
        :param observed: np.array or iterable, required 
            NDArray of observed frequencies in integer format
        :param kwargs: 
            :kwargs n_cols: int, optional
                The number of columns in your observed contingency table
                   if left unset will try to derive from len(observed.T)
                Important if you have a complex array
                (i.e. more than [[*,*,...,*], [*,*,...*]] )

            :kwargs n_rows: int, optional
                The number of rows in your observed contingency table
                  if left unset will try and derive from len(observed)
                Important if you have a complex array 
                (i.e. more than [[*,*,...,*], [*,*,...*]] )
        """
        if isinstance(self._check_array_structure(array=observed), np.array):
            self.observed = observed
        else:
            self._check_array_structure(array=observed)
        invalid_kwargs = associationTests._check_kwargs(kwargs=dict(**kwargs))
        if len(invalid_kwargs) == 0:
            pass
        else:
            raise KeyError("Invalid keyword arguments: %s" % invalid_kwargs)

        self._n_cols = kwargs.get("n_cols", None)
        self._n_rows = kwargs.get("n_rows", None)
        self._n_obs = kwargs.get("n_obs", None)
        if self._n_cols is None or self._n_rows is None:
            try:
                status = self._check_array_structure(array=observed)
            except Exception as err:
                print(err)
            else:
                if type(status) == Exception:
                    raise status
                else:
                    pass
        self.phi_sq = self._calculate_phi_sq(chi2_stat=chi2_stat)

    @property
    def n_rows(self):
        return self._n_rows

    @property
    def n_cols(self):
        return self._n_cols

    @property
    def n_obs(self):
        return self._n_obs

    @n_obs.setter
    def n_obs(self, n):
        self._n_obs = n

    @n_rows.setter
    def n_rows(self, n):
        self._n_rows = n

    @n_cols.setter
    def n_cols(self, n):
        self._n_cols = n

    @staticmethod
    def _check_kwargs(kwargs):
        invalid = []
        for k in kwargs.keys():
            if k not in ["n_rows", "n_cols", "n_obs"]:
                invalid.append(k)
            else:
                pass
        return invalid

    def _check_array_structure(self, array):
        array = np.asarray(array, dtype=np.int64)
        if np.any([_ for _ in array if type(_) != np.int]):
            raise TypeError("Invalid datatype in array")
        if np.any([_ for _ in array if _ < 0]):
            return ValueError("All values in `observed` must be nonnegative.")
        if array.size == 0:
            return ValueError("No data; `observed` has size 0.")
        shape = array.shape
        if len(shape) > 2:
            return ValueError("Currently only two dimensional arrays are supported\n"
                              "Please specify the n_rows, n_cols and n_obs explicitly in the"
                              "class init method")
        elif np.any([i for i in shape if i <= 1]):
            raise ValueError("observed array must have at least 2 rows and cols")
        else:
            self.n_rows = shape[0]
            self.n_cols = shape[1]
            self.n_obs = sum(array.flatten())

            return "Valid"

    def _calculate_phi_sq(self, chi2_stat):
        """
        calculates phi squared from 2x2 matrix

        :return: 
        """
        if self._n_rows == 2 and self._n_cols == 2:
            if len([_ for _ in self.observed.flatten if _ == 0]) == 0:

                a = self.observed[0][0]
                b = self.observed[0][1]
                c = self.observed[1][0]
                d = self.observed[1][1]

                x = (((a * d) - (c * b)) ** 2) / ((a + b) * (c + d) * (a + c) * (b + d))
                return x

            else:
                raise ValueError("Invalid array value. All counts must be > 0")
        else:
            if chi2_stat is None:
                chi2 = chi2_contingency(observed=self.observed)[0]
                return chi2 / self.n_obs

    def _bias_correction(self):
        """
        Corrects bias in estimated value of phi squared derived from nxn contingency 
        tables.

        :return: 

        References
        ----------
        Bergsma, Wicher, "A bias-correction for Cramer's V and Tschuprow's T", 
        London School of Econ. and Pol. Sci., pp. 2-4. 
        http://stats.lse.ac.uk/bergsma/pdf/cramerV3.pdf

        Notes
        ------
        Improves accuracy of estimators with tables > 2x2 and 
        smaller sample sizes, Bergsma (2013)
        """
        phi_sq = self.phi_sq - (1 / (self.n_obs - 1)) * (self.n_rows - 1) * (self.n_cols - 1)
        phi_sq = max(0., phi_sq)
        rows = self.n_rows - (1 / (self.n_obs - 1)) * ((self.n_rows - 1) ** 2)
        cols = self.n_cols - (1 / (self.n_obs - 1)) * ((self.n_cols - 1) ** 2)
        return phi_sq, rows, cols

    def cramers_v(self, correct_bias=True):
        """
        This function computes the Cramer's V, which measures the association
        between sets of variables that are nominal or greater. The returned value 
        is a float in the range of 0.0 to 1.0. Where a value of 0.0 implies complete independence
        among the observed frequencies and a value of 1.0 indicates complete dependence.
        The value of Cramer's V can be heavily biased when sample values are meaningfully
        different from one another, the sample size is large or the number of observed categories
        is high. The function will default to using a bias correction measure however, this can be 
        set to false explicitly.


        :param correct_bias: boolean, required, default: True
        :return: 

        References
        ----------
        .. [1] "Cramér's V", https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V
        .. [2] "Nominal Association: Phi and Cramer's V",
               http://www.people.vcu.edu/~pdattalo/702SuppRead/MeasAssoc/NominalAssoc.html
        .. [3] Bergsma, Wicher, "A bias-correction for Cramer's V and Tschuprow's T", 
               London School of Econ. and Pol. Sci., pp. 2-4. 
               http://stats.lse.ac.uk/bergsma/pdf/cramerV3.pdf


        Examples
        --------
        3x4 Contingency Table - Multinomial
        >>> from scipy.stats import chi2_contingency, associationTests
        >>> obs = np.array([[1, 1, 4, 4], [2, 5, 5, 4], [6, 4, 4, 1]])
        >>> chi2 = chi2_contingency(obs)[0]
        >>> assnTest = associationTests(observed=obs, chi2_stat=chi2)

            with bias correction
        >>> assnTest.cramers_v()
        0.29059448273788646
            without bias correction
        >>> assnTest.cramers_v(correct_bias=False)
        0.5846526848230882


        2x2 Contingency Table
        >>> from scipy.stats import chi2_contingency, associationTests
        >>> obs = np.array([[2, 10], [13, 3]])
        >>> chi2 = chi2_contingency(obs)[0]
        >>> assnTest = associationTests(observed=obs, chi2_stat=chi2)
        chi2_stat could be set to none (as it is ignored in 2x2 matrices) 
        however, it is safer to include it if you are unsure of the matrix shape

            with bias correction
        >>> assnTest.cramers_v()
        0.34062536756037565

            without bias correction
        >>> assnTest.cramers_v(correct_bias=False)
        0.6408461287109103


        Notes (From: [2] "Nominal Association: Phi and Cramer's V")
        ------

        Interpretation: 

        V may be viewed as the association between two variables as a 
        percentage of their maximum possible variation.V2 is the mean square 
        canonical correlation between the variables. For 2-by-2 tables, 
        V = phi (hence some packages like Systat print V only for larger tables).

        Meaning of association: 

        V defines a perfect relationship as one 
        which is predictive or ordered monotonic, and defines a null 
        relationship as statistical independence, as discussed in the 
        section on association. However, the more unequal the marginals, 
        the more V will be less than 1.0.

        Symmetry: 

        V is a symmetrical measure. 
        It does not matter which is the independent (column) variable.

        Data level: 

        V may be used with nominal data or higher.

        Other features: 

        V can reach 1.0 only when the two variables have equal marginals. 

        """

        if correct_bias is True:
            phi_sq, nrows, ncols = self._bias_correction()
        elif correct_bias is False:
            nrows = self.n_rows
            ncols = self.n_cols
            phi_sq = self.phi_sq
        else:
            raise TypeError("correct_bias must be boolean")

        v = math.sqrt(phi_sq / min(ncols - 1, nrows - 1))
        return v

    def tschuprows_t(self, correct_bias=True):
        """
        This function returns the value of Tschuprow's T. The value ranges from 0
        to 1. A value of 0.0 occurs when there is complete independence among the values
        in the contingency table, and a value of 1.0 is returned when there is complete 
        dependence.

        :param correct_bias: boolean, required, default=True

        Because, as in the function for Cramer's V, the chi-squared statistic 
        is used to empirically estimate the value of phi, in this case for contingency tables 
        with greater than 2x2 dimensions. Normally, this produces consistent measurements 
        however, on occasion it can result in a negative phi value (especially as the dimesions of
        the contingency table and sample sizes increase), which would otherwise be impossible. 
        In order to account for this, and other biases introduced in the empirical estimate,
        this argument allows for the measure can be corrected for any such biases.

        :return: 

        References
        ----------
        .. [1] "Tschuprow's T",
               https://en.wikipedia.org/wiki/Tschuprow%27s_T
        .. [2] Bergsma, Wicher, "A bias-correction for Cramer's V and Tschuprow's T", 
               London School of Econ. and Pol. Sci., pp. 5-7. 
               http://stats.lse.ac.uk/bergsma/pdf/cramerV3.pdf
        .. [3] Tschuprow, A. A. (1939) Principles of the Mathematical Theory of Correlation; 
               translated by M. Kantorowitsch. W. Hodge & Co.

        Examples
        --------

        3x4 Contingency Table - Multinomial (with and without bias adjustment)
        >>> from scipy.stats import chi2_contingency, associationTests

        >>> obs = np.array([[1, 1, 4, 4], [2, 5, 5, 4], [6, 4, 4, 1]])
        >>> chi2 = chi2_contingency(obs)[0]
        >>> assnTest = associationTests(observed=obs, chi2_stat=chi2)

            with bias correction
        >>> assnTest.tschuprows_t(correct_bias=True)
        0.2704286376641383
            without bias correction
        >>> assnTest.tschuprows_t(correct_bias=False)
        0.5282933374220176


        2x2 Contingency Table (with and without bias adjustment)
        >>> from scipy.stats import chi2_contingency, association_tests
        >>> obs = np.array([[2, 10], [13, 3]])
        >>> chi2 = chi2_contingency(obs)[0]
        >>> assnTest = association_tests(observed=obs, chi2_stat=chi2)

            with bias correction
        >>> assnTest.tschuprows_t(correct_bias=True)
        0.34062536756037565
            without bias correction
        >>> assnTest.tschuprows_t(correct_bias=False)
        0.6408461287109103


        Notes
        ------

        If the observed set is multinomial this function returns 
        the empircal (or estimated) value of Tshuprow's T which uses 
        Pearsons Chi-Squared Statistic to calculat phi. 
        """

        if correct_bias is True:
            phi_sq, nrows, ncols = self._bias_correction()
        elif correct_bias is False:
            nrows = self.n_rows
            ncols = self.n_cols
            phi_sq = self.phi_sq
        else:
            raise TypeError("correct_bias must be boolean")

        t = math.sqrt(phi_sq / math.sqrt((nrows - 1) * (ncols - 1)))
        return t