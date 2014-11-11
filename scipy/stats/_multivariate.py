#
# Author: Joris Vankerschaver 2013
#
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
from scipy.misc import doccer
from scipy.special import gammaln, psi, multigammaln


__all__ = ['multivariate_normal', 'dirichlet', 'wishart', 'invwishart']

_LOG_2PI = np.log(2 * np.pi)
_LOG_2 = np.log(2)
_LOG_PI = np.log(np.pi)


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
        raise ValueError("Array 'mean' must be a vector of length %d." % dim)
    if cov.ndim == 0:
        cov = cov * np.eye(dim)
    elif cov.ndim == 1:
        cov = np.diag(cov)
    elif cov.ndim == 2 and cov.shape != (dim, dim):
        rows, cols = cov.shape
        if rows != cols:
            msg = ("Array 'cov' must be square if it is two dimensional,"
                   " but cov.shape = %s." % str(cov.shape))
        else:
            msg = ("Dimension mismatch: array 'cov' is of shape %s,"
                   " but 'mean' is a vector of length %d.")
            msg = msg % (str(cov.shape), len(mean))
        raise ValueError(msg)
    elif cov.ndim > 2:
        raise ValueError("Array 'cov' must be at most two-dimensional,"
                         " but cov.ndim = %d" % cov.ndim)

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


def _eigvalsh_to_eps(spectrum, cond=None, rcond=None):
    """
    Determine which eigenvalues are "small" given the spectrum.

    This is for compatibility across various linear algebra functions
    that should agree about whether or not a Hermitian matrix is numerically
    singular and what is its numerical matrix rank.
    This is designed to be compatible with scipy.linalg.pinvh.

    Parameters
    ----------
    spectrum : 1d ndarray
        Array of eigenvalues of a Hermitian matrix.
    cond, rcond : float, optional
        Cutoff for small eigenvalues.
        Singular values smaller than rcond * largest_eigenvalue are
        considered zero.
        If None or -1, suitable machine precision is used.

    Returns
    -------
    eps : float
        Magnitude cutoff for numerical negligibility.

    """
    if rcond is not None:
        cond = rcond
    if cond in [None, -1]:
        t = spectrum.dtype.char.lower()
        factor = {'f': 1E3, 'd': 1E6}
        cond = factor[t] * np.finfo(t).eps
    eps = cond * np.max(abs(spectrum))
    return eps


def _pinv_1d(v, eps=1e-5):
    """
    A helper function for computing the pseudoinverse.

    Parameters
    ----------
    v : iterable of numbers
        This may be thought of as a vector of eigenvalues or singular values.
    eps : float
        Values with magnitude no greater than eps are considered negligible.

    Returns
    -------
    v_pinv : 1d float ndarray
        A vector of pseudo-inverted numbers.

    """
    return np.array([0 if abs(x) <= eps else 1/x for x in v], dtype=float)


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

        eps = _eigvalsh_to_eps(s, cond, rcond)
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
    >>> import matplotlib.pyplot as plt
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
        # TODO replace with np.linalg.slogdet with Numpy 1.5.x not necessary
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
        out = self._mnorm._logpdf(x, self.mean, self.cov_info.U,
                                  self.cov_info.log_pdet, self.cov_info.rank)
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

_dirichlet_doc_default_callparams = """\
alpha : array_like
    The concentration parameters. The number of entries determines the
    dimensionality of the distribution.
"""
_dirichlet_doc_frozen_callparams = ""

_dirichlet_doc_frozen_callparams_note = \
    """See class definition for a detailed description of parameters."""

dirichlet_docdict_params = {
    '_dirichlet_doc_default_callparams': _dirichlet_doc_default_callparams,
}

dirichlet_docdict_noparams = {
    '_dirichlet_doc_default_callparams': _dirichlet_doc_frozen_callparams,
}


def _dirichlet_check_parameters(alpha):
    alpha = np.asarray(alpha)
    if np.min(alpha) <= 0:
        raise ValueError("All parameters must be greater than 0")
    elif alpha.ndim != 1:
        raise ValueError("Parameter vector 'a' must be one dimensional, " +
                         "but a.shape = %s." % str(alpha.shape))
    return alpha


def _dirichlet_check_input(alpha, x):
    x = np.asarray(x)

    if x.shape[0] + 1 != alpha.shape[0] and x.shape[0] != alpha.shape[0]:
        raise ValueError("Vector 'x' must have one entry less then the" +
                         " parameter vector 'a', but alpha.shape = " +
                         "%s and " % alpha.shape +
                         "x.shape = %s." % x.shape)

    if x.shape[0] != alpha.shape[0]:
        xk = np.array([1 - np.sum(x, 0)])
        if xk.ndim == 1:
            x = np.append(x, xk)
        elif xk.ndim == 2:
            x = np.vstack((x, xk))
        else:
            raise ValueError("The input must be one dimensional or a two "
                             "dimensional matrix containing the entries.")

    if np.min(x) < 0:
        raise ValueError("Each entry in 'x' must be greater or equal zero.")

    if np.max(x) > 1:
        raise ValueError("Each entry in 'x' must be smaller or equal one.")

    if (np.abs(np.sum(x, 0) - 1.0) > 10e-10).any():
        raise ValueError("The input vector 'x' must lie within the normal " +
                         "simplex. but sum(x)=%f." % np.sum(x, 0))

    return x


def _lnB(alpha):
    r"""
    Internal helper function to compute the log of the useful quotient

    .. math::
        B(\alpha) = \frac{\prod_{i=1}{K}\Gamma(\alpha_i)}{\Gamma\left(\sum_{i=1}^{K}\alpha_i\right)}

    Parameters
    ----------
    %(_dirichlet_doc_default_callparams)s

    Returns
    -------
    B : scalar
        Helper quotient, internal use only

    """
    return np.sum(gammaln(alpha)) - gammaln(np.sum(alpha))


class dirichlet_gen(object):
    r"""
    A Dirichlet random variable.

    The `alpha` keyword specifies the concentration parameters of the
    distribution.

    .. versionadded:: 0.15.0

    Methods
    -------
    pdf(x, alpha)
        Probability density function.
    logpdf(x, alpha)
        Log of the probability density function.
    rvs(alpha, size=1)
        Draw random samples from a Dirichlet distribution.
    mean(alpha)
        The mean of the Dirichlet distribution
    var(alpha)
        The variance of the Dirichlet distribution
    entropy(alpha)
        Compute the differential entropy of the multivariate normal.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_dirichlet_doc_default_callparams)s

    Alternatively, the object may be called (as a function) to fix
    concentration parameters, returning a "frozen" Dirichlet
    random variable:

    rv = dirichlet(alpha)
        - Frozen object with the same methods but holding the given
          concentration parameters fixed.

    Notes
    -----
    Each :math:`\alpha` entry must be positive. The distribution has only
    support on the simplex defined by

    .. math::
        \sum_{i=1}^{K} x_i \le 1


    The probability density function for `dirichlet` is

    .. math::

        f(x) = \frac{1}{\mathrm{B}(\boldsymbol\alpha)} \prod_{i=1}^K x_i^{\alpha_i - 1}

    where

    .. math::
        \mathrm{B}(\boldsymbol\alpha) = \frac{\prod_{i=1}^K \Gamma(\alpha_i)}{\Gamma\bigl(\sum_{i=1}^K \alpha_i\bigr)}

    and :math:`\boldsymbol\alpha=(\alpha_1,\ldots,\alpha_K)`, the
    concentration parameters and :math:`K` is the dimension of the space
    where :math:`x` takes values.

    """

    def __init__(self):
        self.__doc__ = doccer.docformat(self.__doc__, dirichlet_docdict_params)

    def __call__(self, alpha):
        return dirichlet_frozen(alpha)

    def _logpdf(self, x, alpha):
        """
        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability
            density function
        %(_dirichlet_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        lnB = _lnB(alpha)
        return - lnB + np.sum((np.log(x.T) * (alpha - 1)).T, 0)

    def logpdf(self, x, alpha):
        """
        Log of the Dirichlet probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
        %(_dirichlet_doc_default_callparams)s

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`
        """
        alpha = _dirichlet_check_parameters(alpha)
        x = _dirichlet_check_input(alpha, x)

        out = self._logpdf(x, alpha)
        return _squeeze_output(out)

    def pdf(self, x, alpha):
        """
        The Dirichlet probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
        %(_dirichlet_doc_default_callparams)s

        Returns
        -------
        pdf : ndarray
            The probability density function evaluated at `x`
        """
        alpha = _dirichlet_check_parameters(alpha)
        x = _dirichlet_check_input(alpha, x)

        out = np.exp(self._logpdf(x, alpha))
        return _squeeze_output(out)

    def mean(self, alpha):
        """
        Compute the mean of the dirichlet distribution.

        Parameters
        ----------
        %(_dirichlet_doc_default_callparams)s

        Returns
        -------
        mu : scalar
            Mean of the Dirichlet distribution

        """
        alpha = _dirichlet_check_parameters(alpha)

        out = alpha / (np.sum(alpha))
        return _squeeze_output(out)

    def var(self, alpha):
        """
        Compute the variance of the dirichlet distribution.

        Parameters
        ----------
        %(_dirichlet_doc_default_callparams)s

        Returns
        -------
        v : scalar
            Variance of the Dirichlet distribution

        """

        alpha = _dirichlet_check_parameters(alpha)

        alpha0 = np.sum(alpha)
        out = (alpha * (alpha0 - alpha)) / ((alpha0 * alpha0) * (alpha0 + 1))
        return out

    def entropy(self, alpha):
        """
        Compute the differential entropy of the dirichlet distribution.

        Parameters
        ----------
        %(_dirichlet_doc_default_callparams)s

        Returns
        -------
        h : scalar
            Entropy of the Dirichlet distribution

        """

        alpha = _dirichlet_check_parameters(alpha)

        alpha0 = np.sum(alpha)
        lnB = _lnB(alpha)
        K = alpha.shape[0]

        out = lnB + (alpha0 - K) * scipy.special.psi(alpha0) - np.sum(
            (alpha - 1) * scipy.special.psi(alpha))
        return _squeeze_output(out)

    def rvs(self, alpha, size=1):
        """
        Draw random samples from a Dirichlet distribution.

        Parameters
        ----------
        %(_dirichlet_doc_default_callparams)s
        size : integer, optional
            Number of samples to draw (default 1).


        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `N`), where `N` is the
            dimension of the random variable.

        """
        alpha = _dirichlet_check_parameters(alpha)
        return np.random.dirichlet(alpha, size=size)


dirichlet = dirichlet_gen()


class dirichlet_frozen(object):
    def __init__(self, alpha):
        self.alpha = _dirichlet_check_parameters(alpha)
        self._dirichlet = dirichlet_gen()

    def logpdf(self, x):
        return self._dirichlet.logpdf(x, self.alpha)

    def pdf(self, x):
        return self._dirichlet.pdf(x, self.alpha)

    def mean(self):
        return self._dirichlet.mean(self.alpha)

    def var(self):
        return self._dirichlet.var(self.alpha)

    def entropy(self):
        return self._dirichlet.entropy(self.alpha)

    def rvs(self, size=1):
        return self._dirichlet.rvs(self.alpha, size)


# Set frozen generator docstrings from corresponding docstrings in
# multivariate_normal_gen and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'rvs', 'mean', 'var', 'entropy']:
    method = dirichlet_gen.__dict__[name]
    method_frozen = dirichlet_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, dirichlet_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, dirichlet_docdict_params)


_wishart_doc_default_callparams = """\
df : integer
    Degrees of freedom, must be greater than or equal to dimension of the
    scale matrix
scale : array_like
    Symmetric positive definite scale matrix of the distribution
"""

_wishart_doc_callparams_note = ""

_wishart_doc_frozen_callparams = ""

_wishart_doc_frozen_callparams_note = \
    """See class definition for a detailed description of parameters."""

wishart_docdict_params = {
    '_doc_default_callparams': _wishart_doc_default_callparams,
    '_doc_callparams_note': _wishart_doc_callparams_note
}

wishart_docdict_noparams = {
    '_doc_default_callparams': _wishart_doc_frozen_callparams,
    '_doc_callparams_note': _wishart_doc_frozen_callparams_note
}


class wishart_gen(object):
    r"""
    A Wishart random variable.

    The `df` keyword specifies the degrees of freedom. The `scale` keyword
    specifies the scale matrix, which must be symmetric and positive definite.
    In this context, the scale matrix is often interpreted in terms of a
    multivariate normal precision matrix (the inverse of the covariance
    matrix).

    Methods
    -------
    pdf(x, df, scale)
        Probability density function.
    logpdf(x, df, scale)
        Log of the probability density function.
    rvs(df, scale, size=1)
        Draw random samples from a Wishart distribution.
    entropy()
        Compute the differential entropy of the Wishart distribution.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s

    Alternatively, the object may be called (as a function) to fix the degrees
    of freedom and scale parameters, returning a "frozen" Wishart random
    variable:

    rv = wishart(df=1, scale=1)
        - Frozen object with the same methods but holding the given
          degrees of freedom and scale fixed.

    Notes
    -----
    %(_doc_callparams_note)s

    The scale matrix `scale` must be a symmetric positive definite
    matrix. Singular matrices, including the symmetric positive semi-definite
    case, are not supported.

    The Wishart distribution is often denoted

    .. math::

        W_p(\nu, \Sigma)

    where :math:`\nu` is the degrees of freedom and :math:`\Sigma` is the
    :math:`p \times p` scale matrix.

    The probability density function for `wishart` has support over positive
    definite matrices :math:`S`; if :math:`S \sim W_p(\nu, \Sigma)`, then
    its PDF is given by:

    .. math::

        f(S) = \frac{|S|^{\frac{\nu - p - 1}{2}}}{2^{ \frac{\nu p}{2} } |\Sigma|^\frac{\nu}{2} \Gamma_p \left ( \frac{\nu}{2} \right )}
        e^{-tr(\Sigma^{-1} S) / 2}

    If :math:`S \sim W_p(\nu, \Sigma)` (Wishart) then
    :math:`S^{-1} \sim W_p^{-1}(\nu, \Sigma^{-1})` (inverse Wishart).

    If the scale matrix is 1-dimensional and equal to one, then the Wishart
    distribution :math:`W_1(\nu, 1)` collapses to the :math:`\chi^2(\nu)`
    distribution.

    .. versionadded:: 0.15.0

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import wishart, chi2
    >>> x = np.linspace(1e-5, 8, 100)
    >>> w = wishart.pdf(x, df=3, scale=1); w[:5]
    array([ 0.00126156,  0.10892176,  0.14793434,  0.17400548,  0.1929669 ])
    >>> c = chi2.pdf(x, 3); c[:5]
    array([ 0.00126156,  0.10892176,  0.14793434,  0.17400548,  0.1929669 ])
    >>> plt.plot(x, w)

    The input quantiles can be any shape of array, as long as the last
    axis labels the components.

    References
    ----------
    .. [1] Eaton, Morris L. 2007.
       "Multivariate Statistics: A Vector Space Approach."
       Lecture Notes-Monograph Series 53 (January): i - 512.
    .. [2] Smith, W. B., and R. R. Hocking. 1972.
       "Algorithm AS 53: Wishart Variate Generator."
       Journal of the Royal Statistical Society.
       Series C (Applied Statistics) 21 (3): 341-45. doi:10.2307/2346290.

    """

    def __init__(self):
        self.__doc__ = doccer.docformat(self.__doc__, wishart_docdict_params)

    def __call__(self, df=None, scale=None):
        """
        Create a frozen Wishart distribution.

        See `wishart_frozen` for more information.

        """
        return wishart_frozen(df, scale)

    def _process_parameters(self, df, scale):
        if scale is None:
            scale = 1.0
        scale = np.asarray(scale, dtype=float)

        if scale.ndim == 0:
            scale = scale[np.newaxis,np.newaxis]
        elif scale.ndim == 1:
            scale = np.diag(scale)
        elif scale.ndim == 2 and not scale.shape[0] == scale.shape[1]:
            raise ValueError("Array 'scale' must be square if it is two"
                             " dimensional, but scale.scale = %s."
                             % str(scale.shape))
        elif scale.ndim > 2:
            raise ValueError("Array 'scale' must be at most two-dimensional,"
                             " but scale.ndim = %d" % scale.ndim)

        dim = scale.shape[0]

        if df is None:
            df = dim
        elif not np.isscalar(df):
            raise ValueError("Degrees of freedom must be a scalar.")
        elif df < dim:
            raise ValueError("Degrees of freedom cannot be less than dimension"
                             " of scale matrix, but df = %d" % df)

        return dim, df, scale

    def _process_quantiles(self, x, dim):
        """
        Adjust quantiles array so that last axis labels the components of
        each data point.
        """
        x = np.asarray(x, dtype=float)

        if x.ndim == 0:
            x = x * np.eye(dim)[:, :, np.newaxis]
        if x.ndim == 1:
            if dim == 1:
                x = x[np.newaxis, np.newaxis, :]
            else:
                x = np.diag(x)[:, :, np.newaxis]
        elif x.ndim == 2:
            if not x.shape[0] == x.shape[1]:
                raise ValueError("Quantiles must be square if they are two"
                                 " dimensional, but x.shape = %s."
                                 % str(x.shape))
            x = x[:, :, np.newaxis]
        elif x.ndim == 3:
            if not x.shape[0] == x.shape[1]:
                raise ValueError("Quantiles must be square in the first two"
                                 " dimensions if they are three dimensional"
                                 ", but x.shape = %s." % str(x.shape))
        elif x.ndim > 3:
            raise ValueError("Quantiles must be at most two-dimensional with"
                             " an additional dimension for multiple"
                             "components, but x.ndim = %d" % x.ndim)

        # Now we have 3-dim array; should have shape [dim, dim, *]
        if not x.shape[0:2] == (dim, dim):
            raise ValueError('Quantiles have incompatible dimensions: should'
                             ' be %s, got %s.' % ((dim, dim), x.shape[0:2]))

        return x

    def _process_size(self, size):
        size = np.array(size, dtype=float)

        if size.ndim == 0:
            size = size[np.newaxis]
        elif size.ndim > 1:
            raise ValueError('Size must be an integer or tuple of integers;'
                             ' thus must have dimension <= 1.'
                             ' Got size.ndim = %s' % str(tuple(size)))
        n = size.prod()
        shape = tuple(size)
        
        return n, shape

    def _logpdf(self, x, dim, df, scale, log_det_scale, C):
        """
        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability
            density function
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        scale : ndarray
            Scale matrix
        log_det_scale : float
            Logarithm of the determinant of the scale matrix
        C : ndarray
            Cholesky factorization of the scale matrix, lower triagular.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        # log determinant of x
        # Note: x has components along the last axis, so that x.T has
        # components alone the 0-th axis. Then since det(A) = det(A'), this
        # gives us a 1-dim vector of determinants
        # TODO slogdet is unavailable as long as Numpy 1.5.x is supported
        # s, log_det_x = np.linalg.slogdet(x.T)

        # Retrieve tr(scale^{-1} x)
        log_det_x = np.zeros(x.shape[-1])
        scale_inv_x = np.zeros(x.shape)
        tr_scale_inv_x = np.zeros(x.shape[-1])
        for i in range(x.shape[-1]):
            log_det_x[i] = np.log(np.linalg.det(x[:,:,i]))
            scale_inv_x[:,:,i] = scipy.linalg.cho_solve((C, True), x[:,:,i])
            tr_scale_inv_x[i] = scale_inv_x[:,:,i].diagonal().sum()

        # Log PDF
        out = (
            (0.5 * (df - dim - 1) * log_det_x - 0.5 * tr_scale_inv_x) -
            (0.5 * df * dim * _LOG_2 + 0.5 * df * log_det_scale + multigammaln(0.5*df, dim))
        )

        return out

    def logpdf(self, x, df, scale):
        """
        Log of the Wishart probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
            Each quantile must be a symmetric positive definite matrix.
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        """
        dim, df, scale = self._process_parameters(df, scale)
        x = self._process_quantiles(x, dim)

        # Cholesky decomposition of V, log(det(V)), solve V^{-1} x = `v_inv_x`
        C, lower = scipy.linalg.cho_factor(scale, lower=True)
        log_det_scale = 2 * np.sum(np.log(C.diagonal()))

        out = self._logpdf(x, dim, df, scale, log_det_scale, C)
        return _squeeze_output(out)

    def pdf(self, x, df, scale):
        """
        Wishart probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
            Each quantile must be a symmetric positive definite matrix.
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        """
        return np.exp(self.logpdf(x, df, scale))

    def _mean(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'mean' instead.

        """
        return df * scale

    def mean(self, df, scale):
        """
        Mean of the Wishart distribution

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        mean : float
            The mean of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._mean(dim, df, scale)
        return _squeeze_output(out)

    def _mode(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'mode' instead.

        """
        if df >= dim + 1:
            out = (df-dim-1) * scale
        else:
            out = None
        return out

    def mode(self, df, scale):
        """
        Mode of the Wishart distribution

        Only valid if the degrees of freedom are greater than the dimension of
        the scale matrix.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        mode : float or None
            The Mode of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._mode(dim, df, scale)
        return _squeeze_output(out) if out is not None else out

    def _var(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'var' instead.

        """
        var = scale**2
        diag = scale.diagonal()[np.newaxis, :]  # 1 x dim array
        var += np.outer(diag, diag)
        var *= df
        return var

    def var(self, df, scale):
        """
        Variance of the Wishart distribution

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        var : float
            The variance of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._var(dim, df, scale)
        return _squeeze_output(out)

    def _rvs(self, n, shape, dim, df, scale, C):
        """
        Parameters
        ----------
        n : integer
            Number of variates to generate
        shape : iterable
            Shape of the variates to generate
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        scale : ndarray
            Scale matrix
        C : ndarray
            Cholesky factorization of the scale matrix, lower triagular.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        # Random normal variates for off-diagonal elements
        n_tril = (dim*(dim-1)/2)
        covariances = np.random.normal(size=n*n_tril).reshape(shape+(n_tril,))

        # Random chi-square variates for diagonal elements
        variances = np.r_[
            [np.random.chisquare(df-(i+1)+1, size=n)**0.5 for i in range(dim)]
        ].reshape((dim,) + shape[::-1]).T

        # Create the A matri(ces)
        A = np.zeros(shape + scale.shape)

        # Input the covariances
        size_idx = tuple([slice(None,None,None)]*len(shape))
        tril_idx = np.tril_indices(dim, k=-1)
        A[size_idx + tril_idx] = covariances

        # Input the variances
        diag_idx = np.diag_indices(dim)
        A[size_idx + diag_idx] = variances

        # Calculate C A A' C'
        for index in np.ndindex(shape):
            CA = np.dot(C, A[index])
            A[index] = np.dot(CA, CA.T)

        return A

    def rvs(self, df, scale, size=1):
        """
        Draw random samples from a Wishart distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer or iterable of integers, optional
            Number of samples to draw (default 1).

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        rvs : ndarray
            Random variates of shape (`size`) + (`dim`, `dim), where `dim` is
            the dimension of the scale matrix.

        """
        n, shape = self._process_size(size)
        dim, df, scale = self._process_parameters(df, scale)

        # Cholesky decomposition of scale
        C, lower = scipy.linalg.cho_factor(scale, lower=True)

        out = self._rvs(n, shape, dim, df, scale, C)
        return _squeeze_output(out)

    def _entropy(self, dim, df, log_det_scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        log_det_scale : float
            Logarithm of the determinant of the scale matrix

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'entropy' instead.

        """
        return (
            0.5 * (dim+1) * log_det_scale +
            0.5 * dim * (dim+1) * _LOG_2 +
            multigammaln(0.5*df, dim) -
            0.5 * (df - dim - 1) * np.sum(
                [psi(0.5*(df + 1 - (i+1))) for i in range(dim)]
            ) + 
            0.5 * df * dim
        )

    def entropy(self, df, scale):
        """
        Compute the differential entropy of the Wishart.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        h : scalar
            Entropy of the Wishart distribution

        """
        dim, df, scale = self._process_parameters(df, scale)

        # TODO replace with np.linalg.slogdet when Numpy 1.5.x not necessary
        log_det_scale = np.log(np.linalg.det(scale))
        
        return self._entropy(dim, df, log_det_scale)
wishart = wishart_gen()


class wishart_frozen(object):
    """
    Create a frozen Wishart distribution.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution
    scale : array_like
        Scale matrix of the distribution
    """
    def __init__(self, df, scale):
        self._wishart = wishart_gen()
        self.dim, self.df, self.scale = self._wishart._process_parameters(
            df, scale
        )
        self.C, lower = scipy.linalg.cho_factor(self.scale, lower=True)
        self.log_det_scale = 2 * np.sum(np.log(self.C.diagonal()))

    def logpdf(self, x):
        x = self._wishart._process_quantiles(x, self.dim)

        out = self._wishart._logpdf(x, self.dim, self.df, self.scale,
                                    self.log_det_scale, self.C)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def mean(self):
        out = self._wishart._mean(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def mode(self):
        out = self._wishart._mode(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def var(self):
        out = self._wishart._var(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def rvs(self, size=1):
        n, shape = self._wishart._process_size(size)
        out = self._wishart._rvs(n, shape,
                                 self.dim, self.df, self.scale, self.C)
        return _squeeze_output(out)

    def entropy(self):
        return self._wishart._entropy(self.dim, self.df, self.log_det_scale)

# Set frozen generator docstrings from corresponding docstrings in
# Wishart and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'mean', 'mode', 'var', 'rvs', 'entropy']:
    method = wishart_gen.__dict__[name]
    method_frozen = wishart_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, wishart_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, wishart_docdict_params)


class invwishart_gen(wishart_gen):
    r"""
    An inverse Wishart random variable.

    The `df` keyword specifies the degrees of freedom. The `scale` keyword
    specifies the scale matrix, which must be symmetric and positive definite.
    In this context, the scale matrix is often interpreted in terms of a
    multivariate normal covariance matrix.

    Methods
    -------
    pdf(x, df, scale)
        Probability density function.
    logpdf(x, df, scale)
        Log of the probability density function.
    rvs(df, scale, size=1)
        Draw random samples from an inverse Wishart distribution.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s

    Alternatively, the object may be called (as a function) to fix the degrees
    of freedom and scale parameters, returning a "frozen" inverse Wishart
    random variable:

    rv = invwishart(df=1, scale=1)
        - Frozen object with the same methods but holding the given
          degrees of freedom and scale fixed.

    Notes
    -----
    %(_doc_callparams_note)s

    The scale matrix `scale` must be a symmetric positive definite
    matrix. Singular matrices, including the symmetric positive semi-definite
    case, are not supported.

    The inverse Wishart distribution is often denoted

    .. math::

        W_p^{-1}(\nu, \Psi)

    where :math:`\nu` is the degrees of freedom and :math:`\Psi` is the
    :math:`p \times p` scale matrix.

    The probability density function for `invwishart` has support over positive
    definite matrices :math:`S`; if :math:`S \sim W^{-1}_p(\nu, \Sigma)`,
    then its PDF is given by:

    .. math::

        f(S) = \frac{|\Sigma|^\frac{\nu}{2}}{2^{ \frac{\nu p}{2} } |S|^{\frac{\nu + p + 1}{2}} \Gamma_p \left ( \frac{\nu}{2} \right )}
        e^{-tr(\Sigma S^{-1}) / 2}

    If :math:`S \sim W_p^{-1}(\nu, \Psi)` (inverse Wishart) then
    :math:`S^{-1} \sim W_p(\nu, \Psi^{-1})` (Wishart).

    If the scale matrix is 1-dimensional and equal to one, then the inverse
    Wishart distribution :math:`W_1(\nu, 1)` collapses to the
    inverse Gamma distribution with parameters :math:`shape = \frac{\nu}{2}`
    and :math:`scale=\frac{1}{2}`.

    .. versionadded:: 0.15.0

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import invwishart, invgamma
    >>> x = np.linspace(0.01, 1, 100)
    >>> iw = invwishart.pdf(x, df=6, scale=1); iw[:3]
    array([  1.20546865e-15,   5.42497807e-06,   4.45813929e-03])
    >>> ig = invgamma.pdf(x, 6/2., scale=1./2); c[:3]
    array([  1.20546865e-15,   5.42497807e-06,   4.45813929e-03])
    >>> plt.plot(x, iw)

    The input quantiles can be any shape of array, as long as the last
    axis labels the components.

    References
    ----------
    .. [1] Eaton, Morris L. 2007.
       "Multivariate Statistics: A Vector Space Approach."
       Lecture Notes-Monograph Series 53 (January): i - 512.

    """

    def __init__(self):
        self.__doc__ = doccer.docformat(self.__doc__, wishart_docdict_params)

    def __call__(self, df=None, scale=None):
        """
        Create a frozen inverse Wishart distribution.

        See `invwishart_frozen` for more information.

        """
        return invwishart_frozen(df, scale)

    def _logpdf(self, x, dim, df, scale, log_det_scale):
        """
        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability
            density function.
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        scale : ndarray
            Scale matrix
        log_det_scale : float
            Logarithm of the determinant of the scale matrix

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        log_det_x = np.zeros(x.shape[-1])
        scale_x_inv = np.zeros(x.shape)
        tr_scale_x_inv = np.zeros(x.shape[-1])

        for i in range(x.shape[-1]):
            C, lower = scipy.linalg.cho_factor(x[:,:,i], lower=True)

            log_det_x[i] = 2 * np.sum(np.log(C.diagonal()))

            scale_x_inv[:,:,i] = scipy.linalg.cho_solve((C, True), scale).T
            tr_scale_x_inv[i] = scale_x_inv[:,:,i].diagonal().sum()

        # Log PDF
        out = (
            (0.5 * df * log_det_scale - 0.5 * tr_scale_x_inv) -
            (0.5 * df * dim * _LOG_2 + 0.5 * (df + dim + 1) * log_det_x) -
            multigammaln(0.5*df, dim)
        )

        return out

    def logpdf(self, x, df, scale):
        """
        Log of the inverse Wishart probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
            Each quantile must be a symmetric positive definite matrix.
        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        """
        dim, df, scale = self._process_parameters(df, scale)
        x = self._process_quantiles(x, dim)

        # TODO replace with np.linalg.slogdet when Numpy 1.5.x not necessary
        log_det_scale = np.log(np.linalg.det(scale))

        out = self._logpdf(x, dim, df, scale, log_det_scale)
        return _squeeze_output(out)

    def pdf(self, x, df, scale):
        """
        Inverse Wishart probability density function.

        Parameters
        ----------
        x : array_like
            Quantiles, with the last axis of `x` denoting the components.
            Each quantile must be a symmetric positive definite matrix.

        %(_doc_default_callparams)s

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        """
        return np.exp(self.logpdf(x, df, scale))

    def _mean(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'mean' instead.

        """
        if df > dim + 1:
            out = scale / (df - dim - 1)
        else:
            out = None
        return out

    def mean(self, df, scale):
        """
        Mean of the inverse Wishart distribution

        Only valid if the degrees of freedom are greater than the dimension of
        the scale matrix plus one.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        mean : float or None
            The mean of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._mean(dim, df, scale)
        return _squeeze_output(out) if out is not None else out

    def _mode(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'mode' instead.

        """
        return scale / (df + dim + 1)

    def mode(self, df, scale):
        """
        Mode of the inverse Wishart distribution

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        mode : float
            The Mode of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._mode(dim, df, scale)
        return _squeeze_output(out)

    def _var(self, dim, df, scale):
        """
        Parameters
        ----------
        dim : int
            Dimension of the scale matrix
        %(_doc_default_callparams)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'var' instead.
        """
        if df > dim + 3:
            var = (df - dim + 1) * scale**2
            diag = scale.diagonal()[np.newaxis, :]  # 1 x dim array
            var += (df - dim - 1) * np.outer(diag, diag)
            var /= (df - dim) * (df - dim - 1)**2 * (df - dim - 3)
        else:
            var = None
        return var

    def var(self, df, scale):
        """
        Variance of the inverse Wishart distribution

        Only valid if the degrees of freedom are greater than the dimension of
        the scale matrix plus three.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        var : float
            The variance of the distribution
        """
        dim, df, scale = self._process_parameters(df, scale)
        out = self._var(dim, df, scale)
        return _squeeze_output(out) if out is not None else out

    def _rvs(self, n, shape, dim, df, inv_scale, C, eye):
        """
        Parameters
        ----------
        n : integer
            Number of variates to generate
        shape : iterable
            Shape of the variates to generate
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        inv_scale : ndarray
            Inverse of the scale matrix
        C : ndarray
            Cholesky factorization of the inverse scale matrix, lower
            triagular.
        eye : ndarray
            A `dim`-shape identity matrix.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        # Get random draws from a Wishart with parameter `inv_scale`
        out = super(invwishart_gen, self)._rvs(n, shape, dim, df, inv_scale, C)

        # Invert the draws to get draws from inverse Wishart
        for index in np.ndindex(shape):
            out[index] = scipy.linalg.cho_factor(out[index], lower=True)
            out[index] = scipy.linalg.cho_solve((out[index], True), eye)

        return out

    def rvs(self, df, scale, size=1):
        """
        Draw random samples from an inverse Wishart distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer or iterable of integers, optional
            Number of samples to draw (default 1).

        Notes
        -----
        %(_doc_callparams_note)s

        Returns
        -------
        rvs : ndarray
            Random variates of shape (`size`) + (`dim`, `dim), where `dim` is
            the dimension of the scale matrix.

        """
        n, shape = self._process_size(size)
        dim, df, scale = self._process_parameters(df, scale)

        # Invert the scale
        eye = np.eye(dim)
        L, lower = scipy.linalg.cho_factor(scale, lower=True)
        inv_scale = scipy.linalg.cho_solve((L, lower), eye)
        # Cholesky decomposition of inverted scale
        C, lower = scipy.linalg.cho_factor(inv_scale, lower=True)

        out = self._rvs(n, shape, dim, df, inv_scale, C, eye)

        return _squeeze_output(out)

    def entropy(self):
        # Need to find reference for inverse Wishart entropy
        raise AttributeError

invwishart = invwishart_gen()

class invwishart_frozen(object):
    def __init__(self, df, scale):
        """
        Create a frozen inverse Wishart distribution.

        Parameters
        ----------
        df : array_like
            Degrees of freedom of the distribution
        scale : array_like
            Scale matrix of the distribution
        """
        self._invwishart = invwishart_gen()
        self.dim, self.df, self.scale = self._invwishart._process_parameters(
            df, scale
        )
        self.eye = np.eye(self.dim)

        # Get the determinant via Cholesky factorization
        C, lower = scipy.linalg.cho_factor(self.scale, lower=True)
        self.log_det_scale = 2 * np.sum(np.log(C.diagonal()))

        # Get the inverse using the Cholesky factorization
        self.inv_scale = scipy.linalg.cho_solve((C, lower), self.eye)

        # Get the Cholesky factorization of the inverse scale
        self.C, lower = scipy.linalg.cho_factor(self.inv_scale, lower=True)

    def logpdf(self, x):
        x = self._invwishart._process_quantiles(x, self.dim)
        out = self._invwishart._logpdf(x, self.dim, self.df, self.scale,
                                       self.log_det_scale)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def mean(self):
        out = self._invwishart._mean(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def mode(self):
        out = self._invwishart._mode(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def var(self):
        out = self._invwishart._var(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def rvs(self, size=1):
        n, shape = self._invwishart._process_size(size)

        out = self._invwishart._rvs(n, shape, self.dim, self.df,
                                    self.inv_scale, self.C, self.eye)

        return _squeeze_output(out)

    def entropy(self):
        # Need to find reference for inverse Wishart entropy
        raise AttributeError

# Set frozen generator docstrings from corresponding docstrings in
# inverse Wishart and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'mean', 'mode', 'var', 'rvs']:
    method = invwishart_gen.__dict__[name]
    method_frozen = wishart_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, wishart_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, wishart_docdict_params)
