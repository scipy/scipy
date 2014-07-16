#
# Author: Joris Vankerschaver 2013
#
from __future__ import division, print_function, absolute_import

from scipy.misc import doccer
from functools import wraps
import numpy as np
import scipy.linalg

__all__ = ['multivariate_normal']


_LOG_2PI = np.log(2 * np.pi)


def _process_parameters(dim, mean, cov):
    """
    Infer dimensionality from mean or covariance matrix, ensure that
    mean and covariance are full vector resp. matrix.

    """

    # Try to infer dimensionality
    if dim is None:
        if mean is None:
            if cov is None:
                dim = 1
            else:
                cov = np.asarray(cov, dtype=float)
                if cov.ndim < 2:
                    dim = 1
                else:
                    dim = cov.shape[0]
        else:
            mean = np.asarray(mean, dtype=float)
            dim = mean.size
    else:
        if not np.isscalar(dim):
            raise ValueError("Dimension of random variable must be a scalar.")

    # Check input sizes and return full arrays for mean and cov if necessary
    if mean is None:
        mean = np.zeros(dim)
    mean = np.asarray(mean, dtype=float)

    if cov is None:
        cov = 1.0
    cov = np.asarray(cov, dtype=float)

    if dim == 1:
        mean.shape = (1,)
        cov.shape = (1, 1)

    if mean.ndim != 1 or mean.shape[0] != dim:
        raise ValueError("Array 'mean' must be vector of length %d." % dim)
    if cov.ndim == 0:
        cov = cov * np.eye(dim)
    elif cov.ndim == 1:
        cov = np.diag(cov)
    elif cov.ndim == 2 and cov.shape != (dim, dim):
        raise ValueError("Array 'cov' must be square if it is two dimensional, "
                "but cov.shape = %s" % str(cov.shape))
    elif cov.ndim > 2:
        raise ValueError("Array 'cov' must be at most two-dimensional, "
                "but cov.ndim = %d" % cov.ndim)

    return dim, mean, cov


def _process_quantiles(x, dim):
    """
    Adjust quantiles array so that last axis labels the components of
    each data point.

    """
    x = np.asarray(x, dtype=float)

    if x.ndim == 0:
        x = x[np.newaxis]
    elif x.ndim == 1:
        if dim == 1:
            x = x[:, np.newaxis]
        else:
            x = x[np.newaxis, :]

    return x


def _squeeze_output(out):
    """
    Remove single-dimensional entries from array and convert to scalar,
    if necessary.

    """
    out = out.squeeze()
    if out.ndim == 0:
        out = out[()]
    return out


def _pinv_1d(v, eps=1e-5):
    """
    A helper function for computing the pseudoinverse.

    Parameters
    ----------
    v : iterable of numbers
        This may be thought of as a vector of eigenvalues or singular values.
    eps : float
        Elements of v smaller than eps are considered negligible.

    Returns
    -------
    v_pinv : 1d float ndarray
        A vector of pseudo-inverted numbers.

    """
    return np.array([0 if abs(x) < eps else 1/x for x in v], dtype=float)


class _PSD(object):
    """
    Compute coordinated functions of a symmetric positive semidefinite matrix.

    This class addresses two issues.  Firstly it allows the pseudoinverse,
    the logarithm of the pseudo-determinant, and the rank of the matrix
    to be computed using one call to eigh instead of three.
    Secondly it allows these functions to be computed in a way
    that gives mutually compatible results.
    All of the functions are computed with a common understanding as to
    which of the eigenvalues are to be considered negligibly small.
    The functions are designed to coordinate with scipy.linalg.pinvh()
    but not necessarily with np.linalg.det() or with np.linalg.matrix_rank().

    Parameters
    ----------
    M : 2d array-like
        Symmetric positive semidefinite matrix.
    cond, rcond : float, optional
        Cutoff for small eigenvalues.
        Singular values smaller than rcond * largest_eigenvalue are
        considered zero.
        If None or -1, suitable machine precision is used.
    lower : bool, optional
        Whether the pertinent array data is taken from the lower
        or upper triangle of M. (Default: lower)
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite
        numbers. Disabling may give a performance gain, but may result
        in problems (crashes, non-termination) if the inputs do contain
        infinities or NaNs.
    allow_singular : bool, optional
        Whether to allow a singular matrix.  (Default: True)

    Notes
    -----
    The arguments are similar to those of scipy.linalg.pinvh().

    """
    def __init__(self, M, cond=None, rcond=None, lower=True,
            check_finite=True, allow_singular=True):
        # Compute the symmetric eigendecomposition.
        # Note that eigh takes care of array conversion, chkfinite,
        # and assertion that the matrix is square.
        s, u = scipy.linalg.eigh(M, lower=lower, check_finite=check_finite)

        # This looks redundant but its purpose is to keep the code
        # compatible with the scipy.linalg.pinvh function
        # so that it may be more easily reorganized later if desired.
        if rcond is not None:
            cond = rcond
        if cond in [None, -1]:
            t = u.dtype.char.lower()
            factor = {'f': 1E3, 'd': 1E6}
            cond = factor[t] * np.finfo(t).eps

        eps = cond * np.max(abs(s))
        if np.min(s) < -eps:
            raise ValueError('the input matrix must be positive semidefinite')
        d = s[s > eps]
        if len(d) < len(s) and not allow_singular:
            raise np.linalg.LinAlgError('singular matrix')
        s_pinv = _pinv_1d(s, eps)
        U = np.multiply(u, np.sqrt(s_pinv))

        # Initialize the eagerly precomputed attributes.
        self.rank = len(d)
        self.U = U
        self.log_pdet = np.sum(np.log(d))

        # Initialize an attribute to be lazily computed.
        self._pinv = None

    @property
    def pinv(self):
        if self._pinv is None:
            self._pinv = np.dot(self.U, self.U.T)
        return self._pinv


_doc_default_callparams = """\
mean : array_like, optional
    Mean of the distribution (default zero)
cov : array_like, optional
    Covariance matrix of the distribution (default one)
allow_singular : bool, optional
    Whether to allow a singular covariance matrix.  (Default: False)
"""

_doc_callparams_note = \
"""Setting the parameter `mean` to `None` is equivalent to having `mean`
be the zero-vector. The parameter `cov` can be a scalar, in which case
the covariance matrix is the identity times that value, a vector of
diagonal entries for the covariance matrix, or a two-dimensional
array_like.
"""

_doc_frozen_callparams = ""

_doc_frozen_callparams_note = \
"""See class definition for a detailed description of parameters."""

docdict_params = {
    '_doc_default_callparams': _doc_default_callparams,
    '_doc_callparams_note': _doc_callparams_note
}

docdict_noparams = {
    '_doc_default_callparams': _doc_frozen_callparams,
    '_doc_callparams_note': _doc_frozen_callparams_note
}


class multivariate_normal_gen(object):
    r"""
    A multivariate normal random variable.

    The `mean` keyword specifies the mean. The `cov` keyword specifies the
    covariance matrix.

    Methods
    -------
    pdf(x, mean=None, cov=1, allow_singular=False)
        Probability density function.
    logpdf(x, mean=None, cov=1, allow_singular=False)
        Log of the probability density function.
    rvs(mean=None, cov=1, allow_singular=False, size=1)
        Draw random samples from a multivariate normal distribution.
    entropy()
        Compute the differential entropy of the multivariate normal.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s

    Alternatively, the object may be called (as a function) to fix the mean
    and covariance parameters, returning a "frozen" multivariate normal
    random variable:

    rv = multivariate_normal(mean=None, cov=1, allow_singular=False)
        - Frozen object with the same methods but holding the given
          mean and covariance fixed.

    Notes
    -----
    %(_doc_callparams_note)s

    The covariance matrix `cov` must be a (symmetric) positive
    semi-definite matrix. The determinant and inverse of `cov` are computed
    as the pseudo-determinant and pseudo-inverse, respectively, so
    that `cov` does not need to have full rank.

    The probability density function for `multivariate_normal` is

    .. math::

        f(x) = \frac{1}{\sqrt{(2 \pi)^k \det \Sigma}} \exp\left( -\frac{1}{2} (x - \mu)^T \Sigma^{-1} (x - \mu) \right),

    where :math:`\mu` is the mean, :math:`\Sigma` the covariance matrix,
    and :math:`k` is the dimension of the space where :math:`x` takes values.

    .. versionadded:: 0.14.0

    Examples
    --------
    >>> from scipy.stats import multivariate_normal
    >>> x = np.linspace(0, 5, 10, endpoint=False)
    >>> y = multivariate_normal.pdf(x, mean=2.5, cov=0.5); y
    array([ 0.00108914,  0.01033349,  0.05946514,  0.20755375,  0.43939129,
            0.56418958,  0.43939129,  0.20755375,  0.05946514,  0.01033349])
    >>> plt.plot(x, y)

    The input quantiles can be any shape of array, as long as the last
    axis labels the components.  This allows us for instance to
    display the frozen pdf for a non-isotropic random variable in 2D as
    follows:

    >>> x, y = np.mgrid[-1:1:.01, -1:1:.01]
    >>> pos = np.empty(x.shape + (2,))
    >>> pos[:, :, 0] = x; pos[:, :, 1] = y
    >>> rv = multivariate_normal([0.5, -0.2], [[2.0, 0.3], [0.3, 0.5]])
    >>> plt.contourf(x, y, rv.pdf(pos))

    """

    def __init__(self):
        self.__doc__ = doccer.docformat(self.__doc__, docdict_params)

    def __call__(self, mean=None, cov=1, allow_singular=False):
        """
        Create a frozen multivariate normal distribution.

        See `multivariate_normal_frozen` for more information.

        """
        return multivariate_normal_frozen(mean, cov,
                allow_singular=allow_singular)

    def _logpdf(self, x, mean, prec_U, log_det_cov, rank):
        """
        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability
            density function
        mean : ndarray
            Mean of the distribution
        prec_U : ndarray
            A decomposition such that np.dot(prec_U, prec_U.T)
            is the precision matrix, i.e. inverse of the covariance matrix.
        log_det_cov : float
            Logarithm of the determinant of the covariance matrix
        rank : int
            Rank of the covariance matrix.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        dev = x - mean
        maha = np.sum(np.square(np.dot(dev, prec_U)), axis=-1)
        return -0.5 * (rank * _LOG_2PI + log_det_cov + maha)

    def logpdf(self, x, mean, cov, allow_singular=False):
        """
        Log of the multivariate normal probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        x = _process_quantiles(x, dim)
        psd = _PSD(cov, allow_singular=allow_singular)
        out = self._logpdf(x, mean, psd.U, psd.log_pdet, psd.rank)
        return _squeeze_output(out)

    def pdf(self, x, mean, cov, allow_singular=False):
        """
        Multivariate normal probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        x = _process_quantiles(x, dim)
        psd = _PSD(cov, allow_singular=allow_singular)
        out = np.exp(self._logpdf(x, mean, psd.U, psd.log_pdet, psd.rank))
        return _squeeze_output(out)

    def rvs(self, mean=None, cov=1, size=1):
        """
        Draw random samples from a multivariate normal distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer, optional
            Number of samples to draw (default 1).

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `N`), where `N` is the
            dimension of the random variable.

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        out = np.random.multivariate_normal(mean, cov, size)
        return _squeeze_output(out)

    def entropy(self, mean=None, cov=1):
        """
        Compute the differential entropy of the multivariate normal.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        h : scalar
            Entropy of the multivariate normal distribution

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        return 0.5 * np.log(np.linalg.det(2 * np.pi * np.e * cov))

multivariate_normal = multivariate_normal_gen()


class multivariate_normal_frozen(object):
    def __init__(self, mean=None, cov=1, allow_singular=False):
        """
        Create a frozen multivariate normal distribution.

        Parameters
        ----------
        mean : array_like, optional
            Mean of the distribution (default zero)
        cov : array_like, optional
            Covariance matrix of the distribution (default one)
        allow_singular : bool, optional
            If this flag is True then tolerate a singular
            covariance matrix (default False).

        Examples
        --------
        When called with the default parameters, this will create a 1D random
        variable with mean 0 and covariance 1:

        >>> from scipy.stats import multivariate_normal
        >>> r = multivariate_normal()
        >>> r.mean
        array([ 0.])
        >>> r.cov
        array([[1.]])

        """
        self.dim, self.mean, self.cov = _process_parameters(None, mean, cov)
        self.cov_info = _PSD(self.cov, allow_singular=allow_singular)
        self._mnorm = multivariate_normal_gen()

    def logpdf(self, x):
        x = _process_quantiles(x, self.dim)
        out = self._mnorm._logpdf(x, self.mean,
                self.cov_info.U, self.cov_info.log_pdet, self.cov_info.rank)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def rvs(self, size=1):
        return self._mnorm.rvs(self.mean, self.cov, size)

    def entropy(self):
        """
        Computes the differential entropy of the multivariate normal.

        Returns
        -------
        h : scalar
            Entropy of the multivariate normal distribution

        """
        log_pdet = self.cov_info.log_pdet
        rank = self.cov_info.rank
        return 0.5 * (rank * (_LOG_2PI + 1) + log_pdet)


# Set frozen generator docstrings from corresponding docstrings in
# multivariate_normal_gen and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'rvs']:
    method = multivariate_normal_gen.__dict__[name]
    method_frozen = multivariate_normal_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(method.__doc__, docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, docdict_params)


def _dirichlet_check_parameters(a):
    if min(a) <= 0:
        raise ValueError("All parameters must be greater than 0")
    elif a.ndim != 1:
        raise ValueError("Parameter vector 'a' must be one dimensional, " +
                         "but a.shape = %s." % str(a.shape))
    return a


def _dirichlet_check_input(a, x):

    if x.shape != a.shape:
        raise ValueError("Vector 'x' must have the same shape as" +
                         " parameter vector 'a', but a.shape = %s and " +
                         "x.shape = %s." % (a.shape, x.shape))
    if min(x) < 0:
        raise ValueError("Each entry in 'x' must be greater or equal zero.")

    if max(x) > 1:
        raise ValueError("Each entry in 'x' must be smaller or equal one.")

    if np.sum(x) > 1:
        raise ValueError("The input vector 'x' must lie within the normal simplex.")

    return x


class dirichlet_gen(object):
    def __init__(self):
        pass

    def __call__(self, a, x=None):
        return dirichlet_frozen(a, x)


    def _B(self, a):
        return np.sum(np.gammaln(a)) - np.gammaln(np.sum(a))

    def _logpdf(self, x, a):
        a = _dirichlet_check_parameters(a)
        x = _dirichlet_check_input(a, x)

        B = self._B(a)
        out = - B + np.sum((a - 1) * np.log(x))
        return _squeeze_output(out)

    def _pdf(self, x, a):
        a = _dirichlet_check_parameters(a)
        x = _dirichlet_check_input(a, x)

        return np.exp(self._logpdf(a, x))

    def _mean(self, a):
        a = _dirichlet_check_parameters(a)

        out = a/(np.sum(a))
        return _squeeze_output(out)

    def _var(self, a):
        a = _dirichlet_check_parameters(a)

        a0 = np.sum(a)
        out =  a * (a0 - 1)/((a0 * a0) * (a0 + 1))
        return _squeeze_output(out)

    def _entropy(self, a):
        a = _dirichlet_check_parameters(a)

        a0 = np.sum(a)
        B = self._B(a)
        K = a.shape[0]

        out =  np.log(B) + (a0 - K) * np.psi(a0) - np.sum((a - 1) * np.psi(a))
        return _squeeze_output(out)

    def _rvs(self, a):
        return np.random.dirichlet(a)


dirichlet = dirichlet_gen()


class dirichlet_frozen(object):
    def __init__(self, a):
        self.a = _dirichlet_check_parameters(a)
        self._dirichlet = dirichlet_gen()

    def logpdf(self, x):
        return self._dirichlet._logpdf(x, self.a)

    def pdf(self, x):
        return self._dirichlet._pdf(x, self.a)

    def mean(self):
        return self._dirichlet._mean(self.a)

    def var(self):
        return self._dirichlet._var(self.a)

    def entropy(self):
        return self._dirichlet._entropy(self.a)

    def rvs(self):
        return self._dirichlet._rvs(self.a)






