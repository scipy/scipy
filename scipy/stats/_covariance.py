from functools import cached_property

import numpy as np
from scipy import linalg


__all__ = ["Covariance"]


class Covariance:
    """
    Representation of a covariance matrix

    Calculations involving covariance matrices (e.g. data whitening,
    multivariate normal function evaluation) are often performed more
    efficiently using a decomposition of the covariance matrix instead of the
    covariance metrix itself. This class allows the user to construct an
    object representing a covariance matrix using any of several
    decompositions and perform calculations using a common interface.

    Examples
    --------
    The most common use of the `Covariance` class is to call one of the
    factory methods to create a `Covariance` object, then pass that
    representation of the `Covariance` matrix as a shape parameter of a
    multivariate distribution.

    For instance, the multivariate normal distribution can accept an array
    representing a covariance matrix:

    >>> from scipy import stats
    >>> d = [1, 2, 3]
    >>> A = np.diag(d)  # a diagonal covariance matrix
    >>> x = [4, -2, 5]  # a point of interest
    >>> dist = stats.multivariate_normal(mean=[0, 0, 0], cov=A)
    >>> dist.pdf(x)
    4.9595685102808205e-08

    but the calculations are performed in a very generic way that does not
    take advantage of any special properties of the covariance matrix. Because
    our covariance matrix is diagonal, we can use ``Covariance.from_diagonal``
    to create an object representing the covariance matrix, and
    `multivariate_normal` can use this to compute the probability density
    function more efficiently.

    >>> cov = stats.Covariance.from_diagonal(d)
    >>> dist = stats.multivariate_normal(mean=[0, 0, 0], cov=cov)
    >>> dist.pdf(x)
    4.9595685102808205e-08

    """

    @staticmethod
    def from_diagonal(diagonal):
        r"""
        Return a representation of a covariance matrix from its diagonal.

        Parameters
        ----------
        diagonal : array_like
            The diagonal elements of a diagonal matrix.

        Notes
        -----
        Let the diagonal elements of a diagonal covariance matrix :math:`D` be
        stored in the vector :math:`d`.

        When all elements of :math:`d` are strictly positive, whitening of a
        data point :math:`x` is performed by computing
        :math:`x \cdot d^{-1/2}`, where the inverse square root can be taken
        element-wise.
        :math:`\log\det{D}` is calculated as :math:`-2 \sum(\log{d})`,
        where the :math:`\log` operation is performed element-wise.

        This `Covariance` class supports singular covariance matrices. When
        computing ``_log_pdet``, non-positive elements of :math:`d` are
        ignored. Whitening is not well defined when the point to be whitened
        does not lie in the span of the columns of the covariance matrix. The
        convention taken here is to treat the inverse square root of
        non-positive elements of :math:`d` as zeros.
        """
        return CovViaDiagonal(diagonal)

    @staticmethod
    def from_precision(precision, covariance=None):
        r"""
        Return a representation of a covariance from its precision matrix.

        Parameters
        ----------
        precision : array_like
            The precision matrix; that is, the inverse of a square, symmetric,
            positive definite covariance matrix.
        covariance : array_like, optional
            The square, symmetric, positive definite covariance matrix. If not
            provided, this may need to be calculated (e.g. to evaluate the
            cumulative distribution function of
            `scipy.stats.multivariate_normal`) by inverting `precision`.

        Notes
        -----
        Let the covariance matrix be :math:`A`, its precision matrix be
        :math:`P = A^{-1}`, and :math:`L` be the lower Cholesky factor such
        that :math:`L L^T = P`.
        Whitening of a data point :math:`x` is performed by computing
        :math:`x^T L`. :math:`\log\det{A}` is calculated as
        :math:`-2tr(\log{P})`, where the :math:`\log` operation is performed
        element-wise.

        This `Covariance` class does not support singular covariance matrices
        because the precision matrix does not exist for a singular covariance
        matrix.
        """
        return CovViaPrecision(precision, covariance)

    def whiten(self, x):
        """
        Perform a whitening transformation on data.

        "Whitening" ("white" as in "white noise", in which each frequency has
        equal magnitude) transforms a set of random variables into a new set of
        random variables with unit-diagonal covariance. When a whitening
        transform is applied to a sample of points distributed according to
        a multivariate normal distribution with zero mean, the covariance of
        the transformed sample is approximately the identity matrix.

        Parameters
        ----------
        x : array_like
            An array of points. The last dimension must correspond with the
            dimensionality of the space, i.e., the number of columns in the
            covariance matrix.

        Returns
        -------
        x_ : array_like
            The transformed array of points.

        References
        ----------
        .. [1] "Whitening Transformation". Wikipedia.
               https://en.wikipedia.org/wiki/Whitening_transformation

        Examples
        --------
        >>> import numpy as np
        >>> from scipy import stats
        >>> rng = np.random.default_rng()
        >>> n = 3
        >>> A = rng.random(size=(n, n))
        >>> cov_array = A @ A.T  # make matrix symmetric positive definite
        >>> cov_object = stats.CovViaPrecision(np.linalg.inv(cov_array))
        >>> x = rng.multivariate_normal(np.zeros(n), cov_array, size=(10000))
        >>> x_ = cov_object.whiten(x)
        >>> np.cov(x_, rowvar=False)  # near-identity covariance
        array([[0.97862122, 0.00893147, 0.02430451],
               [0.00893147, 0.96719062, 0.02201312],
               [0.02430451, 0.02201312, 0.99206881]])

        """
        return self._whiten(np.asarray(x))

    @property
    def log_pdet(self):
        """
        Log of the pseudo-determinant of the covariance matrix
        """
        return np.array(self._log_pdet, dtype=float)[()]

    @property
    def rank(self):
        """
        Rank of the covariance matrix
        """
        return np.array(self._rank, dtype=int)[()]

    @property
    def covariance(self):
        """
        Explicit representation of the covariance matrix
        """
        return self._covariance

    @property
    def dimensionality(self):
        """
        Dimensionality of the vector space
        """
        return np.array(self._dimensionality, dtype=int)[()]

    @property
    def shape(self):
        """
        Shape of the covariance array
        """
        return self._shape

    def _validate_matrix(self, A, name):
        A = np.atleast_2d(A)
        m, n = A.shape[-2:]
        if m != n or A.ndim != 2 or not (np.issubdtype(A.dtype, np.integer) or
                                         np.issubdtype(A.dtype, np.floating)):
            message = (f"The input `{name}` must be a square, "
                       "two-dimensional array of real numbers.")
            raise ValueError(message)
        return A

    def _validate_vector(self, A, name):
        A = np.atleast_1d(A)
        if A.ndim != 1 or not (np.issubdtype(A.dtype, np.integer) or
                               np.issubdtype(A.dtype, np.floating)):
            message = (f"The input `{name}` must be a one-dimensional array "
                       "of real numbers.")
            raise ValueError(message)
        return A


class CovViaPrecision(Covariance):

    def __init__(self, precision, covariance=None):
        precision = self._validate_matrix(precision, 'precision')
        if covariance is not None:
            covariance = self._validate_matrix(covariance, 'covariance')
            message = "`precision.shape` must equal `covariance.shape`."
            if precision.shape != covariance.shape:
                raise ValueError(message)

        self._chol_P = np.linalg.cholesky(precision)
        self._log_pdet = -2*np.log(np.diag(self._chol_P)).sum(axis=-1)
        self._rank = precision.shape[-1]  # must be full rank if invertible
        self._precision = precision
        self._cov_matrix = covariance
        self._dimensionality = self._rank
        self._shape = precision.shape
        self._allow_singular = False

    def _whiten(self, x):
        return x @ self._chol_P

    @cached_property
    def _covariance(self):
        n = self._dimensionality
        return (linalg.cho_solve((self._chol_P, True), np.eye(n))
                if self._cov_matrix is None else self._cov_matrix)


def _dot_diag(x, d):
    # If d were a full diagonal matrix, x @ d would always do what we want
    # This is for when `d` is compressed to include only the diagonal elements
    return x * d if x.ndim < 2 else x * np.expand_dims(d, -2)


class CovViaDiagonal(Covariance):

    def __init__(self, diagonal):
        diagonal = self._validate_vector(diagonal, 'diagonal')

        i_zero = diagonal <= 0
        positive_diagonal = np.array(diagonal, dtype=np.float64)

        positive_diagonal[i_zero] = 1  # ones don't affect determinant
        self._log_pdet = np.sum(np.log(positive_diagonal), axis=-1)

        psuedo_reciprocals = 1 / np.sqrt(positive_diagonal)
        psuedo_reciprocals[i_zero] = 0

        self._LP = psuedo_reciprocals
        self._rank = positive_diagonal.shape[-1] - i_zero.sum(axis=-1)
        self._covariance = np.apply_along_axis(np.diag, -1, diagonal)
        self._dimensionality = diagonal.shape[-1]
        self._i_zero = i_zero
        self._shape = self._covariance.shape
        self._allow_singular = True

    def _whiten(self, x):
        return _dot_diag(x, self._LP)

    def _support_mask(self, x):
        """
        Check whether x lies in the support of the distribution.
        """
        return ~np.any(_dot_diag(x, self._i_zero), axis=-1)


class CovViaPSD(Covariance):
    """
    Representation of a covariance provided via an instance of _PSD
    """

    def __init__(self, psd):
        self._LP = psd.U
        self._log_pdet = psd.log_pdet
        self._rank = psd.rank
        self._covariance = psd._M
        self._dimensionality = psd._M.shape[-1]
        self._shape = psd._M.shape
        self._psd = psd
        self._allow_singular = False  # by default

    def _whiten(self, x):
        return x @ self._LP

    def _support_mask(self, x):
        return self._psd._support_mask(x)
