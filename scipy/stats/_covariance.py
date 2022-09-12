from functools import cached_property

import numpy as np
from scipy import linalg


__all__ = ["CovViaPrecision"]


class Covariance():
    """
    Representation of a covariance matrix
    """

    def whiten(self, x):
        """
        Perform a whitening transformation on data.

        "Whitening" ("white" as in "white noise", in which each frequency has
        equal magnitude) transforms a set of random variables into a new set of
        random variable with unit-diagonal covariance. When a whitening
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
        >>> cov_array = A @ A.T
        >>> cov_object = stats.CovViaPrecision(np.linalg.inv(cov_array))
        >>> x = rng.multivariate_normal(np.zeros(n), cov_array, size=(10000))
        >>> x_ = cov_object.whiten(x)
        >>> np.cov(x_, rowvar=False)  # near-identity covariance
        array([[0.97862122, 0.00893147, 0.02430451],
               [0.00893147, 0.96719062, 0.02201312],
               [0.02430451, 0.02201312, 0.99206881]])

        """
        return self._whiten(np.asarray(x))

    @cached_property
    def log_pdet(self):
        """
        Log of the pseudo-determinant of the covariance matrix
        """
        return np.array(self._log_pdet, dtype=float)[()]

    @cached_property
    def rank(self):
        """
        Rank of the covariance matrix
        """
        return np.array(self._rank, dtype=int)[()]

    @cached_property
    def covariance(self):
        """
        Explicit representation of the covariance matrix
        """
        return self._covariance

    @cached_property
    def dimensionality(self):
        """
        Dimensionality of the vector space
        """
        return np.array(self._dimensionality, dtype=int)[()]

    @cached_property
    def shape(self):
        """
        Shape of the covariance array
        """
        return self._shape

    def _validate_matrix(self, A, name):
        A = np.atleast_2d(A)
        m, n = A.shape[-2:]
        if m != n or A.ndim != 2 or not np.issubdtype(A.dtype, np.number):
            message = (f"The input `{name}` must be a square, "
                       "two-dimensional array of numbers.")
            raise ValueError(message)
        return A

    def _validate_vector(self, A, name):
        A = np.atleast_1d(A)
        if A.ndim != 1 or not np.issubdtype(A.dtype, np.number):
            message = (f"The input `{name}` must be a one-dimensional array "
                       "of numbers.")
            raise ValueError(message)
        return A


class CovViaPrecision(Covariance):
    r"""
    Representation of a covariance provided via its precision matrix.

    Notes
    -----
    Let the covariance matrix be :math:`A`, its precision matrix be
    :math:`P = A^{-1}`, and :math:`L` be the lower Cholesky factor such that
    :math:`L L^T = P`.
    Whitening of a data point :math:`x` is performed by computing
    :math:`x^T L`. :math:`\log\det{A}` is calculated as :math:`-2tr(\log{P})`,
    where the :math:`\log` operation is performed element-wise.

    This `Covariance` class does not support singular covariance matrices
    because the precision matrix does not exist for a singular covariance
    matrix.

    """

    def __init__(self, precision, covariance=None):
        """
        Parameters
        ----------
        precision : array_like
            The precision matrix; that is, the inverse of a square, symmetric,
            positive definite covariance matrix.
        covariance : array_like, optional
            The square, symmetric, positive definite covariance matrix. If not
            provided, this may need to be calculated (e.g. to evaluate the
            cumulative distribution function of
            `scipy.stats.multivariate_normal`)
        """
        precision = self._validate_matrix(precision, 'precision')
        if covariance is not None:
            covariance = self._validate_matrix(covariance, 'covariance')
            message = "`precision.shape` must equal `covariance.shape`."
            if precision.shape != covariance.shape:
                raise ValueError(message)

        self._chol_P = np.linalg.cholesky(precision)
        self._log_pdet = -2*np.log(np.diag(self._chol_P)).sum(axis=-1)
        self._rank = precision.shape[-1]  # must be full rank in invertible
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
