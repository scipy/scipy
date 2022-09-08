from functools import cached_property

import numpy as np


__all__ = ["Covariance", "CovViaPrecision", "CovViaPSD"]


def _T(A):
    if A.ndim < 2:
        return A
    else:
        return np.swapaxes(A, -1, -2)


def _dot_diag(x, d):
    # If d were a full diagonal matrix, x @ d would always do what we want
    # This is for when `d` is compressed to include only the diagonal elements
    return x * d if x.ndim < 2 else x * np.expand_dims(d, -2)


class Covariance():
    """
    Representation of a covariance matrix as needed by multivariate_normal
    """
    # The last two axes of matrix-like input represent the dimensionality
    # In the case of diagonal elements or eigenvalues, in which a diagonal
    # matrix has been reduced to one dimension, the last dimension
    # of the input represents the dimensionality. Internally, it
    # will be expanded in the second to last axis as needed.
    # Matrix math works, but instead of the fundamental matrix-vector
    # operation being A@x, think of x are row vectors that pre-multiply A.

    def whiten(self, x):
        """
        Right multiplication by the left square root of the precision matrix.
        """
        return self._whiten(x)

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
    def A(self):
        """
        Explicit representation of the covariance matrix
        """
        return self._A

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
    """
    Representation of a covariance provided via the precision matrix
    """

    def __init__(self, precision, covariance=None):
        precision = self._validate_matrix(precision, 'precision')
        if covariance is not None:
            covariance = self._validate_matrix(covariance, 'covariance')
            message = "`precision.shape` must equal `covariance.shape`."
            if precision.shape != covariance.shape:
                raise ValueError(message)

        self._LP = np.linalg.cholesky(precision)
        self._log_pdet = -2*np.log(np.diag(self._LP)).sum(axis=-1)
        self._rank = precision.shape[-1]  # must be full rank in invertible
        self._precision = precision
        self._covariance = covariance
        self._dimensionality = self._rank
        self._shape = precision.shape
        self._allow_singular = False

    def _whiten(self, x):
        return x @ self._LP

    @cached_property
    def _A(self):
        return (np.linalg.inv(self._precision) if self._covariance is None
                else self._covariance)


class CovViaPSD(Covariance):
    """
    Representation of a covariance provided via an instance of _PSD
    """

    def __init__(self, psd):
        self._LP = psd.U
        self._log_pdet = psd.log_pdet
        self._rank = psd.rank
        self._A = psd._M
        self._dimensionality = psd._M.shape[-1]
        self._shape = psd._M.shape
        self._psd = psd
        self._allow_singular = False  # by default

    def _whiten(self, x):
        return x @ self._LP

    def _support_mask(self, x):
        return self._psd._support_mask(x)
