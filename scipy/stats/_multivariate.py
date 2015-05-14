#
# Author: Joris Vankerschaver 2013
#
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
from scipy.misc import doccer
from scipy.special import gammaln, psi, multigammaln
from scipy._lib._util import check_random_state, _asarray_validated
from scipy.optimize import minimize
from scipy.optimize import fmin_l_bfgs_b


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
    M : array_like
        Symmetric positive semidefinite matrix (2-D).
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

_doc_random_state = """\
random_state : None or int or np.random.RandomState instance, optional
    If int or RandomState, use it for drawing the random variates.
    If None (or np.random), the global np.random state is used.
    Default is None.
"""

_doc_frozen_callparams = ""

_doc_frozen_callparams_note = \
    """See class definition for a detailed description of parameters."""

docdict_params = {
    '_doc_default_callparams': _doc_default_callparams,
    '_doc_callparams_note': _doc_callparams_note,
    '_doc_random_state': _doc_random_state
}

docdict_noparams = {
    '_doc_default_callparams': _doc_frozen_callparams,
    '_doc_callparams_note': _doc_frozen_callparams_note,
    '_doc_random_state': _doc_random_state
}


class multi_rv_generic(object):
    """
    Class which encapsulates common functionality between all multivariate
    distributions.

    """
    def __init__(self, seed=None):
        super(multi_rv_generic, self).__init__()
        self._random_state = check_random_state(seed)

    @property
    def random_state(self):
        """ Get or set the RandomState object for generating random variates.

        This can be either None or an existing RandomState object.

        If None (or np.random), use the RandomState singleton used by np.random.
        If already a RandomState instance, use it.
        If an int, use a new RandomState instance seeded with seed.

        """
        return self._random_state

    @random_state.setter
    def random_state(self, seed):
        self._random_state = check_random_state(seed)

    def _get_random_state(self, random_state):
        if random_state is not None:
            return check_random_state(random_state)
        else:
            return self._random_state


class multi_rv_frozen(object):
    """
    Class which encapsulates common functionality between all frozen
    multivariate distributions.
    """
    @property
    def random_state(self):
        return self._dist._random_state

    @random_state.setter
    def random_state(self, seed):
        self._dist._random_state = check_random_state(seed)


class multivariate_normal_gen(multi_rv_generic):
    r"""
    A multivariate normal random variable.

    The `mean` keyword specifies the mean. The `cov` keyword specifies the
    covariance matrix.

    Methods
    -------
    ``pdf(x, mean=None, cov=1, allow_singular=False)``
        Probability density function.
    ``logpdf(x, mean=None, cov=1, allow_singular=False)``
        Log of the probability density function.
    ``rvs(mean=None, cov=1, size=1, random_state=None)``
        Draw random samples from a multivariate normal distribution.
    ``entropy()``
        Compute the differential entropy of the multivariate normal.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s
    %(_doc_random_state)s

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

        f(x) = \frac{1}{\sqrt{(2 \pi)^k \det \Sigma}}
               \exp\left( -\frac{1}{2} (x - \mu)^T \Sigma^{-1} (x - \mu) \right),

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
    >>> fig1 = plt.figure()
    >>> ax = fig1.add_subplot(111)
    >>> ax.plot(x, y)

    The input quantiles can be any shape of array, as long as the last
    axis labels the components.  This allows us for instance to
    display the frozen pdf for a non-isotropic random variable in 2D as
    follows:

    >>> x, y = np.mgrid[-1:1:.01, -1:1:.01]
    >>> pos = np.empty(x.shape + (2,))
    >>> pos[:, :, 0] = x; pos[:, :, 1] = y
    >>> rv = multivariate_normal([0.5, -0.2], [[2.0, 0.3], [0.3, 0.5]])
    >>> fig2 = plt.figure()
    >>> ax2 = fig2.add_subplot(111)
    >>> ax2.contourf(x, y, rv.pdf(pos))

    """

    def __init__(self, seed=None):
        super(multivariate_normal_gen, self).__init__(seed)
        self.__doc__ = doccer.docformat(self.__doc__, docdict_params)

    def __call__(self, mean=None, cov=1, allow_singular=False, seed=None):
        """
        Create a frozen multivariate normal distribution.

        See `multivariate_normal_frozen` for more information.

        """
        return multivariate_normal_frozen(mean, cov,
                                          allow_singular=allow_singular,
                                          seed=seed)

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

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

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

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        x = _process_quantiles(x, dim)
        psd = _PSD(cov, allow_singular=allow_singular)
        out = np.exp(self._logpdf(x, mean, psd.U, psd.log_pdet, psd.rank))
        return _squeeze_output(out)

    def rvs(self, mean=None, cov=1, size=1, random_state=None):
        """
        Draw random samples from a multivariate normal distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `N`), where `N` is the
            dimension of the random variable.

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, mean, cov = _process_parameters(None, mean, cov)

        random_state = self._get_random_state(random_state)
        out = random_state.multivariate_normal(mean, cov, size)
        return _squeeze_output(out)

    def entropy(self, mean=None, cov=1):
        """
        Compute the differential entropy of the multivariate normal.

        Parameters
        ----------
        %(_doc_default_callparams)s

        Returns
        -------
        h : scalar
            Entropy of the multivariate normal distribution

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, mean, cov = _process_parameters(None, mean, cov)
        _, logdet = np.linalg.slogdet(2 * np.pi * np.e * cov)
        return 0.5 * logdet


multivariate_normal = multivariate_normal_gen()


class multivariate_normal_frozen(multi_rv_frozen):
    def __init__(self, mean=None, cov=1, allow_singular=False, seed=None):
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
        seed : None or int or np.random.RandomState instance, optional
            This parameter defines the RandomState object to use for drawing
            random variates.
            If None (or np.random), the global np.random state is used.
            If integer, it is used to seed the local RandomState instance
            Default is None.

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
        self._dist = multivariate_normal_gen(seed)

    def logpdf(self, x):
        x = _process_quantiles(x, self.dim)
        out = self._dist._logpdf(x, self.mean, self.cov_info.U,
                                 self.cov_info.log_pdet, self.cov_info.rank)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def rvs(self, size=1, random_state=None):
        return self._dist.rvs(self.mean, self.cov, size, random_state)

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
    '_doc_random_state': _doc_random_state
}

dirichlet_docdict_noparams = {
    '_dirichlet_doc_default_callparams': _dirichlet_doc_frozen_callparams,
    '_doc_random_state': _doc_random_state
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


class dirichlet_gen(multi_rv_generic):
    r"""
    A Dirichlet random variable.

    The `alpha` keyword specifies the concentration parameters of the
    distribution.

    .. versionadded:: 0.15.0

    Methods
    -------
    ``pdf(x, alpha)``
        Probability density function.
    ``logpdf(x, alpha)``
        Log of the probability density function.
    ``rvs(alpha, size=1, random_state=None)``
        Draw random samples from a Dirichlet distribution.
    ``mean(alpha)``
        The mean of the Dirichlet distribution
    ``var(alpha)``
        The variance of the Dirichlet distribution
    ``entropy(alpha)``
        Compute the differential entropy of the multivariate normal.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_dirichlet_doc_default_callparams)s
    %(_doc_random_state)s

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

        \mathrm{B}(\boldsymbol\alpha) = \frac{\prod_{i=1}^K \Gamma(\alpha_i)}
                                     {\Gamma\bigl(\sum_{i=1}^K \alpha_i\bigr)}

    and :math:`\boldsymbol\alpha=(\alpha_1,\ldots,\alpha_K)`, the
    concentration parameters and :math:`K` is the dimension of the space
    where :math:`x` takes values.

    """

    def __init__(self, seed=None):
        super(dirichlet_gen, self).__init__(seed)
        self.__doc__ = doccer.docformat(self.__doc__, dirichlet_docdict_params)

    def __call__(self, alpha, seed=None):
        return dirichlet_frozen(alpha, seed=seed)

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
            Log of the probability density function evaluated at `x`.

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
            The probability density function evaluated at `x`.

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

        out = lnB + (alpha0 - K) * psi(alpha0) - np.sum((alpha - 1) * psi(alpha))
        return _squeeze_output(out)

    def rvs(self, alpha, size=1, random_state=None):
        """
        Draw random samples from a Dirichlet distribution.

        Parameters
        ----------
        %(_dirichlet_doc_default_callparams)s
        size : int, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `N`), where `N` is the
            dimension of the random variable.

        """
        alpha = _dirichlet_check_parameters(alpha)
        random_state = self._get_random_state(random_state)
        return random_state.dirichlet(alpha, size=size)

    def fit(self, data, x0=None):
        """
        Estimate the Dirichlet distribution parameter values.

        Parameters
        ----------
        data : 2d array-like
            A sequence of multivariate observations.
        x0 : 1d array-like, optional
            Initial guess of the parameter vector alpha.

        Returns
        -------
        shape, loc, scale : tuple of 1d arrays
            The shape is the alpha vector maximum likelihood estimate.
            The loc and scale are hard-coded as the zero-vector
            and the one-vector respectively.

        References
        ----------
        .. [1] Jonathan Huang, "Maximum Likelihood Estimation of Dirichlet
           Distribution Parameters"
           http://jonathan-huang.org/research/dirichlet/dirichlet.pdf

        Notes
        -----
        This fit is computed by maximizing a log likelihood function.

        .. versionadded:: 0.17.0

        """
        X = _asarray_validated(data)
        if X.ndim != 2:
            raise NotImplementedError('expected a sequence of observations')
        n, k = X.shape
        if x0 is None:
            x0 = np.ones(k)
        else:
            x0 = _asarray_validated(x0)
            if x0.shape[0] != k:
                raise ValueError('the shape of the initial parameter vector '
                                 'guess is incompatible with the shape of '
                                 'the data')
            if np.any(np.less(x0, 0)):
                raise ValueError('the initial guess must be positive')
        # Compute the sufficient statistics
        log_p_hat = np.log(X).mean(axis=0)

        # Use log of alpha for max likelihood estimation.
        x0 = np.log(x0)

        def func(log_alpha):
            a = np.exp(log_alpha)
            scaled_ll = np.expm1(log_alpha).dot(log_p_hat) - _lnB(a)
            scaled_score = a * (psi(np.sum(a)) - psi(a) + log_p_hat)
            return -scaled_ll, -scaled_score

        log_a, f, d = fmin_l_bfgs_b(func, x0)
        a = np.exp(log_a)

        return a, np.zeros_like(a), np.ones_like(a)


dirichlet = dirichlet_gen()


class dirichlet_frozen(multi_rv_frozen):
    def __init__(self, alpha, seed=None):
        self.alpha = _dirichlet_check_parameters(alpha)
        self._dist = dirichlet_gen(seed)

    def logpdf(self, x):
        return self._dist.logpdf(x, self.alpha)

    def pdf(self, x):
        return self._dist.pdf(x, self.alpha)

    def mean(self):
        return self._dist.mean(self.alpha)

    def var(self):
        return self._dist.var(self.alpha)

    def entropy(self):
        return self._dist.entropy(self.alpha)

    def rvs(self, size=1, random_state=None):
        return self._dist.rvs(self.alpha, size, random_state)


# Set frozen generator docstrings from corresponding docstrings in
# multivariate_normal_gen and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'rvs', 'mean', 'var', 'entropy']:
    method = dirichlet_gen.__dict__[name]
    method_frozen = dirichlet_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, dirichlet_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, dirichlet_docdict_params)


_wishart_doc_default_callparams = """\
df : int
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
    '_doc_callparams_note': _wishart_doc_callparams_note,
    '_doc_random_state': _doc_random_state
}

wishart_docdict_noparams = {
    '_doc_default_callparams': _wishart_doc_frozen_callparams,
    '_doc_callparams_note': _wishart_doc_frozen_callparams_note,
    '_doc_random_state': _doc_random_state
}


class wishart_gen(multi_rv_generic):
    r"""
    A Wishart random variable.

    The `df` keyword specifies the degrees of freedom. The `scale` keyword
    specifies the scale matrix, which must be symmetric and positive definite.
    In this context, the scale matrix is often interpreted in terms of a
    multivariate normal precision matrix (the inverse of the covariance
    matrix).

    Methods
    -------
    ``pdf(x, df, scale)``
        Probability density function.
    ``logpdf(x, df, scale)``
        Log of the probability density function.
    ``rvs(df, scale, size=1, random_state=None)``
        Draw random samples from a Wishart distribution.
    ``entropy()``
        Compute the differential entropy of the Wishart distribution.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s
    %(_doc_random_state)s

    Alternatively, the object may be called (as a function) to fix the degrees
    of freedom and scale parameters, returning a "frozen" Wishart random
    variable:

    rv = wishart(df=1, scale=1)
        - Frozen object with the same methods but holding the given
          degrees of freedom and scale fixed.

    See Also
    --------
    invwishart, chi2

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

        f(S) = \frac{|S|^{\frac{\nu - p - 1}{2}}}{2^{ \frac{\nu p}{2} }
               |\Sigma|^\frac{\nu}{2} \Gamma_p \left ( \frac{\nu}{2} \right )}
               \exp\left( -tr(\Sigma^{-1} S) / 2 \right)

    If :math:`S \sim W_p(\nu, \Sigma)` (Wishart) then
    :math:`S^{-1} \sim W_p^{-1}(\nu, \Sigma^{-1})` (inverse Wishart).

    If the scale matrix is 1-dimensional and equal to one, then the Wishart
    distribution :math:`W_1(\nu, 1)` collapses to the :math:`\chi^2(\nu)`
    distribution.

    .. versionadded:: 0.16.0

    References
    ----------
    .. [1] M.L. Eaton, "Multivariate Statistics: A Vector Space Approach",
           Wiley, 1983.
    .. [2] W.B. Smith and R.R. Hocking, "Algorithm AS 53: Wishart Variate
           Generator", Applied Statistics, vol. 21, pp. 341-345, 1972.

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

    """

    def __init__(self, seed=None):
        super(wishart_gen, self).__init__(seed)
        self.__doc__ = doccer.docformat(self.__doc__, wishart_docdict_params)

    def __call__(self, df=None, scale=None, seed=None):
        """
        Create a frozen Wishart distribution.

        See `wishart_frozen` for more information.

        """
        return wishart_frozen(df, scale, seed)

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
        size = np.asarray(size)

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

        # Retrieve tr(scale^{-1} x)
        log_det_x = np.zeros(x.shape[-1])
        scale_inv_x = np.zeros(x.shape)
        tr_scale_inv_x = np.zeros(x.shape[-1])
        for i in range(x.shape[-1]):
            _, log_det_x[i] = self._cholesky_logdet(x[:,:,i])
            scale_inv_x[:,:,i] = scipy.linalg.cho_solve((C, True), x[:,:,i])
            tr_scale_inv_x[i] = scale_inv_x[:,:,i].trace()

        # Log PDF
        out = ((0.5 * (df - dim - 1) * log_det_x - 0.5 * tr_scale_inv_x) -
               (0.5 * df * dim * _LOG_2 + 0.5 * df * log_det_scale +
                multigammaln(0.5*df, dim)))

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

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, df, scale = self._process_parameters(df, scale)
        x = self._process_quantiles(x, dim)

        # Cholesky decomposition of scale, get log(det(scale))
        C, log_det_scale = self._cholesky_logdet(scale)

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

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

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
        diag = scale.diagonal()  # 1 x dim array
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

    def _standard_rvs(self, n, shape, dim, df, random_state):
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
        random_state : np.random.RandomState instance
            RandomState used for drawing the random variates.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        # Random normal variates for off-diagonal elements
        n_tril = dim * (dim-1) // 2
        covariances = random_state.normal(
            size=n*n_tril).reshape(shape+(n_tril,))

        # Random chi-square variates for diagonal elements
        variances = np.r_[[random_state.chisquare(df-(i+1)+1, size=n)**0.5
             for i in range(dim)]].reshape((dim,) + shape[::-1]).T

        # Create the A matri(ces) - lower triangular
        A = np.zeros(shape + (dim, dim))

        # Input the covariances
        size_idx = tuple([slice(None,None,None)]*len(shape))
        tril_idx = np.tril_indices(dim, k=-1)
        A[size_idx + tril_idx] = covariances

        # Input the variances
        diag_idx = np.diag_indices(dim)
        A[size_idx + diag_idx] = variances

        return A

    def _rvs(self, n, shape, dim, df, C, random_state):
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
            Cholesky factorization of the scale matrix, lower triangular.
        %(_doc_random_state)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        random_state = self._get_random_state(random_state)
        # Calculate the matrices A, which are actually lower triangular
        # Cholesky factorizations of a matrix B such that B ~ W(df, I)
        A = self._standard_rvs(n, shape, dim, df, random_state)

        # Calculate SA = C A A' C', where SA ~ W(df, scale)
        # Note: this is the product of a (lower) (lower) (lower)' (lower)'
        #       or, denoting B = AA', it is C B C' where C is the lower
        #       triangular Cholesky factorization of the scale matrix.
        #       this appears to conflict with the instructions in [1]_, which
        #       suggest that it should be D' B D where D is the lower
        #       triangular factorization of the scale matrix. However, it is
        #       meant to refer to the Bartlett (1933) representation of a
        #       Wishart random variate as L A A' L' where L is lower triangular
        #       so it appears that understanding D' to be upper triangular
        #       is either a typo in or misreading of [1]_.
        for index in np.ndindex(shape):
            CA = np.dot(C, A[index])
            A[index] = np.dot(CA, CA.T)

        return A

    def rvs(self, df, scale, size=1, random_state=None):
        """
        Draw random samples from a Wishart distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer or iterable of integers, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray
            Random variates of shape (`size`) + (`dim`, `dim), where `dim` is
            the dimension of the scale matrix.

        Notes
        -----
        %(_doc_callparams_note)s

        """
        n, shape = self._process_size(size)
        dim, df, scale = self._process_parameters(df, scale)

        # Cholesky decomposition of scale
        C = scipy.linalg.cholesky(scale, lower=True)

        out = self._rvs(n, shape, dim, df, C, random_state)

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

        Returns
        -------
        h : scalar
            Entropy of the Wishart distribution

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, df, scale = self._process_parameters(df, scale)
        _, log_det_scale = self._cholesky_logdet(scale)
        return self._entropy(dim, df, log_det_scale)

    def _cholesky_logdet(self, scale):
        """
        Compute Cholesky decomposition and determine (log(det(scale)).

        Parameters
        ----------
        scale : ndarray
            Scale matrix.

        Returns
        -------
        c_decomp : ndarray
            The Cholesky decomposition of `scale`.
        logdet : scalar
            The log of the determinant of `scale`.

        Notes
        -----
        This computation of ``logdet`` is equivalent to
        ``np.linalg.slogdet(scale)``.  It is ~2x faster though.

        """
        c_decomp = scipy.linalg.cholesky(scale, lower=True)
        logdet = 2 * np.sum(np.log(c_decomp.diagonal()))
        return c_decomp, logdet
wishart = wishart_gen()


class wishart_frozen(multi_rv_frozen):
    """
    Create a frozen Wishart distribution.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution
    scale : array_like
        Scale matrix of the distribution
    seed : None or int or np.random.RandomState instance, optional
        This parameter defines the RandomState object to use for drawing
        random variates.
        If None (or np.random), the global np.random state is used.
        If integer, it is used to seed the local RandomState instance
        Default is None.

    """
    def __init__(self, df, scale, seed=None):
        self._dist = wishart_gen(seed)
        self.dim, self.df, self.scale = self._dist._process_parameters(
            df, scale)
        self.C, self.log_det_scale = self._dist._cholesky_logdet(self.scale)

    def logpdf(self, x):
        x = self._dist._process_quantiles(x, self.dim)

        out = self._dist._logpdf(x, self.dim, self.df, self.scale,
                                 self.log_det_scale, self.C)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def mean(self):
        out = self._dist._mean(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def mode(self):
        out = self._dist._mode(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def var(self):
        out = self._dist._var(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def rvs(self, size=1, random_state=None):
        n, shape = self._dist._process_size(size)
        out = self._dist._rvs(n, shape, self.dim, self.df,
                              self.C, random_state)
        return _squeeze_output(out)

    def entropy(self):
        return self._dist._entropy(self.dim, self.df, self.log_det_scale)

# Set frozen generator docstrings from corresponding docstrings in
# Wishart and fill in default strings in class docstrings
for name in ['logpdf', 'pdf', 'mean', 'mode', 'var', 'rvs', 'entropy']:
    method = wishart_gen.__dict__[name]
    method_frozen = wishart_frozen.__dict__[name]
    method_frozen.__doc__ = doccer.docformat(
        method.__doc__, wishart_docdict_noparams)
    method.__doc__ = doccer.docformat(method.__doc__, wishart_docdict_params)


from numpy import asarray_chkfinite, asarray
from scipy.linalg.misc import LinAlgError
from scipy.linalg.lapack import get_lapack_funcs
def _cho_inv_batch(a, check_finite=True):
    """
    Invert the matrices a_i, using a Cholesky factorization of A, where
    a_i resides in the last two dimensions of a and the other indices describe
    the index i.

    Overwrites the data in a.

    Parameters
    ----------
    a : array
        Array of matrices to invert, where the matrices themselves are stored
        in the last two dimensions.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    x : array
        Array of inverses of the matrices ``a_i``.

    See also
    --------
    scipy.linalg.cholesky : Cholesky factorization of a matrix

    """
    if check_finite:
        a1 = asarray_chkfinite(a)
    else:
        a1 = asarray(a)
    if len(a1.shape) < 2 or a1.shape[-2] != a1.shape[-1]:
        raise ValueError('expected square matrix in last two dimensions')

    potrf, potri = get_lapack_funcs(('potrf','potri'), (a1,))

    tril_idx = np.tril_indices(a.shape[-2], k=-1)
    triu_idx = np.triu_indices(a.shape[-2], k=1)
    for index in np.ndindex(a1.shape[:-2]):

        # Cholesky decomposition
        a1[index], info = potrf(a1[index], lower=True, overwrite_a=False,
                                clean=False)
        if info > 0:
            raise LinAlgError("%d-th leading minor not positive definite"
                              % info)
        if info < 0:
            raise ValueError('illegal value in %d-th argument of internal'
                             ' potrf' % -info)
        # Inversion
        a1[index], info = potri(a1[index], lower=True, overwrite_c=False)
        if info > 0:
            raise LinAlgError("the inverse could not be computed")
        if info < 0:
            raise ValueError('illegal value in %d-th argument of internal'
                             ' potrf' % -info)

        # Make symmetric (dpotri only fills in the lower triangle)
        a1[index][triu_idx] = a1[index][tril_idx]

    return a1


class invwishart_gen(wishart_gen):
    r"""
    An inverse Wishart random variable.

    The `df` keyword specifies the degrees of freedom. The `scale` keyword
    specifies the scale matrix, which must be symmetric and positive definite.
    In this context, the scale matrix is often interpreted in terms of a
    multivariate normal covariance matrix.

    Methods
    -------
    ``pdf(x, df, scale)``
        Probability density function.
    ``logpdf(x, df, scale)``
        Log of the probability density function.
    ``rvs(df, scale, size=1, random_state=None)``
        Draw random samples from an inverse Wishart distribution.

    Parameters
    ----------
    x : array_like
        Quantiles, with the last axis of `x` denoting the components.
    %(_doc_default_callparams)s
    %(_doc_random_state)s

    Alternatively, the object may be called (as a function) to fix the degrees
    of freedom and scale parameters, returning a "frozen" inverse Wishart
    random variable:

    rv = invwishart(df=1, scale=1)
        - Frozen object with the same methods but holding the given
          degrees of freedom and scale fixed.

    See Also
    --------
    wishart

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

        f(S) = \frac{|\Sigma|^\frac{\nu}{2}}{2^{ \frac{\nu p}{2} }
               |S|^{\frac{\nu + p + 1}{2}} \Gamma_p \left(\frac{\nu}{2} \right)}
               \exp\left( -tr(\Sigma S^{-1}) / 2 \right)

    If :math:`S \sim W_p^{-1}(\nu, \Psi)` (inverse Wishart) then
    :math:`S^{-1} \sim W_p(\nu, \Psi^{-1})` (Wishart).

    If the scale matrix is 1-dimensional and equal to one, then the inverse
    Wishart distribution :math:`W_1(\nu, 1)` collapses to the
    inverse Gamma distribution with parameters shape = :math:`\frac{\nu}{2}`
    and scale = :math:`\frac{1}{2}`.

    .. versionadded:: 0.16.0

    References
    ----------
    .. [1] M.L. Eaton, "Multivariate Statistics: A Vector Space Approach",
           Wiley, 1983.
    .. [2] M.C. Jones, "Generating Inverse Wishart Matrices", Communications in
           Statistics - Simulation and Computation, vol. 14.2, pp.511-514, 1985.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import invwishart, invgamma
    >>> x = np.linspace(0.01, 1, 100)
    >>> iw = invwishart.pdf(x, df=6, scale=1)
    >>> iw[:3]
    array([  1.20546865e-15,   5.42497807e-06,   4.45813929e-03])
    >>> ig = invgamma.pdf(x, 6/2., scale=1./2)
    >>> ig[:3]
    array([  1.20546865e-15,   5.42497807e-06,   4.45813929e-03])
    >>> plt.plot(x, iw)

    The input quantiles can be any shape of array, as long as the last
    axis labels the components.

    """

    def __init__(self, seed=None):
        super(invwishart_gen, self).__init__(seed)
        self.__doc__ = doccer.docformat(self.__doc__, wishart_docdict_params)

    def __call__(self, df=None, scale=None, seed=None):
        """
        Create a frozen inverse Wishart distribution.

        See `invwishart_frozen` for more information.

        """
        return invwishart_frozen(df, scale, seed)

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
        #scale_x_inv = np.zeros(x.shape)
        x_inv = np.copy(x).T
        if dim > 1:
            _cho_inv_batch(x_inv)  # works in-place
        else:
            x_inv = 1./x_inv
        tr_scale_x_inv = np.zeros(x.shape[-1])

        for i in range(x.shape[-1]):
            C, lower = scipy.linalg.cho_factor(x[:,:,i], lower=True)

            log_det_x[i] = 2 * np.sum(np.log(C.diagonal()))

            #scale_x_inv[:,:,i] = scipy.linalg.cho_solve((C, True), scale).T
            tr_scale_x_inv[i] = np.dot(scale, x_inv[i]).trace()

        # Log PDF
        out = ((0.5 * df * log_det_scale - 0.5 * tr_scale_x_inv) -
               (0.5 * df * dim * _LOG_2 + 0.5 * (df + dim + 1) * log_det_x) -
               multigammaln(0.5*df, dim))

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

        Returns
        -------
        pdf : ndarray
            Log of the probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

        """
        dim, df, scale = self._process_parameters(df, scale)
        x = self._process_quantiles(x, dim)
        _, log_det_scale = self._cholesky_logdet(scale)
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

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at `x`

        Notes
        -----
        %(_doc_callparams_note)s

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
            diag = scale.diagonal()  # 1 x dim array
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

    def _rvs(self, n, shape, dim, df, C, random_state):
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
        C : ndarray
            Cholesky factorization of the scale matrix, lower triagular.
        %(_doc_random_state)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        random_state = self._get_random_state(random_state)
        # Get random draws A such that A ~ W(df, I)
        A = super(invwishart_gen, self)._standard_rvs(n, shape, dim,
                                                      df, random_state)

        # Calculate SA = (CA)'^{-1} (CA)^{-1} ~ iW(df, scale)
        eye = np.eye(dim)
        trtrs = get_lapack_funcs(('trtrs'), (A,))

        for index in np.ndindex(A.shape[:-2]):
            # Calculate CA
            CA = np.dot(C, A[index])
            # Get (C A)^{-1} via triangular solver
            if dim > 1:
                CA, info = trtrs(CA, eye, lower=True)
                if info > 0:
                    raise LinAlgError("Singular matrix.")
                if info < 0:
                    raise ValueError('Illegal value in %d-th argument of'
                                     ' internal trtrs' % -info)
            else:
                CA = 1. / CA
            # Get SA
            A[index] = np.dot(CA.T, CA)

        return A

    def rvs(self, df, scale, size=1, random_state=None):
        """
        Draw random samples from an inverse Wishart distribution.

        Parameters
        ----------
        %(_doc_default_callparams)s
        size : integer or iterable of integers, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray
            Random variates of shape (`size`) + (`dim`, `dim), where `dim` is
            the dimension of the scale matrix.

        Notes
        -----
        %(_doc_callparams_note)s

        """
        n, shape = self._process_size(size)
        dim, df, scale = self._process_parameters(df, scale)

        # Invert the scale
        eye = np.eye(dim)
        L, lower = scipy.linalg.cho_factor(scale, lower=True)
        inv_scale = scipy.linalg.cho_solve((L, lower), eye)
        # Cholesky decomposition of inverted scale
        C = scipy.linalg.cholesky(inv_scale, lower=True)

        out = self._rvs(n, shape, dim, df, C, random_state)

        return _squeeze_output(out)

    def entropy(self):
        # Need to find reference for inverse Wishart entropy
        raise AttributeError

invwishart = invwishart_gen()

class invwishart_frozen(multi_rv_frozen):
    def __init__(self, df, scale, seed=None):
        """
        Create a frozen inverse Wishart distribution.

        Parameters
        ----------
        df : array_like
            Degrees of freedom of the distribution
        scale : array_like
            Scale matrix of the distribution
        seed : None or int or np.random.RandomState instance, optional
            This parameter defines the RandomState object to use for drawing
            random variates.
            If None (or np.random), the global np.random state is used.
            If integer, it is used to seed the local RandomState instance
            Default is None.

        """
        self._dist = invwishart_gen(seed)
        self.dim, self.df, self.scale = self._dist._process_parameters(
            df, scale
        )

        # Get the determinant via Cholesky factorization
        C, lower = scipy.linalg.cho_factor(self.scale, lower=True)
        self.log_det_scale = 2 * np.sum(np.log(C.diagonal()))

        # Get the inverse using the Cholesky factorization
        eye = np.eye(self.dim)
        self.inv_scale = scipy.linalg.cho_solve((C, lower), eye)

        # Get the Cholesky factorization of the inverse scale
        self.C = scipy.linalg.cholesky(self.inv_scale, lower=True)

    def logpdf(self, x):
        x = self._dist._process_quantiles(x, self.dim)
        out = self._dist._logpdf(x, self.dim, self.df, self.scale,
                                 self.log_det_scale)
        return _squeeze_output(out)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def mean(self):
        out = self._dist._mean(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def mode(self):
        out = self._dist._mode(self.dim, self.df, self.scale)
        return _squeeze_output(out)

    def var(self):
        out = self._dist._var(self.dim, self.df, self.scale)
        return _squeeze_output(out) if out is not None else out

    def rvs(self, size=1, random_state=None):
        n, shape = self._dist._process_size(size)

        out = self._dist._rvs(n, shape, self.dim, self.df,
                              self.C, random_state)

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
