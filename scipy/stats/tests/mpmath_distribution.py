import numpy as np
from mpmath import mp

mp.dps = 100  # default in case the user forgets to set it

class Distribution:
    # Minimalist distribution infrastructure for generating reference values
    # of distribution functions. No edge case handling.

    @staticmethod
    def _make_mpf_array(x):
        shape = np.shape(x)
        x = np.asarray(x, dtype=np.float64).ravel()
        return np.asarray([mp.mpf(xi) for xi in x]).reshape(shape)[()]

    def __init__(self, **kwargs):
        self._params = {key:self._make_mpf_array(val)
                        for key, val in kwargs.items()}

    def _pdf(self, x):
        raise NotImplementedError("_pdf must be overridden.")

    def _cdf(self, x, **kwargs):
        a, _ = self._support(**kwargs)
        return mp.quad(lambda x: self._pdf(x, **kwargs), (a, x))

    def _sf(self, x, **kwargs):
        _, b = self._support(**kwargs)
        return mp.quad(lambda x: self._pdf(x, **kwargs), (x, b))

    def _logpdf(self, x, **kwargs):
        return mp.log(self._pdf(x, **kwargs))

    def _logcdf(self, x, **kwargs):
        return mp.log(self._cdf(x, **kwargs))

    def _logsf(self, x, **kwargs):
        return mp.log(self._sf(x, **kwargs))

    def _support(self, **kwargs):
        return -mp.inf, mp.inf

    def _entropy(self, **kwargs):
        def integrand(x):
            logpdf = self._logpdf(x, **kwargs)
            pdf = mp.exp(logpdf)
            return  -pdf*logpdf

        a, b = self._support(**kwargs)
        return mp.quad(integrand, (a, b))

    def _mean(self, **kwargs):
        return self._moment(order=1, center=0, **kwargs)

    def _var(self, **kwargs):
        mu = self._mean(**kwargs)
        return self._moment(order=2, center=mu, **kwargs)

    def _skew(self, **kwargs):
        mu = self._mean(**kwargs)
        u2 = self._moment(order=2, center=mu, **kwargs)
        sigma = mp.sqrt(u2)
        u3 = self._moment(order=3, center=mu, **kwargs)
        return u3 / sigma**3

    def _kurtosis(self, **kwargs):
        mu = self._mean(**kwargs)
        u2 = self._moment(order=2, center=mu, **kwargs)
        u4 = self._moment(order=4, center=mu, **kwargs)
        return u4 / u2**2 - 3

    def _moment(self, order, center=None, **kwargs):
        def integrand(x):
            return self._pdf(x, **kwargs)*(x - center)**order

        if center is None:
            center = self._mean(**kwargs)

        a, b = self._support(**kwargs)
        return mp.quad(integrand, (a, b))

    def pdf(self, x):
        fun = np.vectorize(self._pdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def cdf(self, x):
        fun = np.vectorize(self._cdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def sf(self, x):
        fun = np.vectorize(self._sf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def logpdf(self, x):
        fun = np.vectorize(self._logpdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def logcdf(self, x):
        fun = np.vectorize(self._logcdf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def logsf(self, x):
        fun = np.vectorize(self._logsf)
        x = self._make_mpf_array(x)
        res = fun(x, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def support(self):
        fun = np.vectorize(self._support)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def entropy(self):
        fun = np.vectorize(self._entropy)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def mean(self):
        fun = np.vectorize(self._mean)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def var(self):
        fun = np.vectorize(self._var)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def skew(self):
        fun = np.vectorize(self._skew)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def kurtosis(self):
        fun = np.vectorize(self._kurtosis)
        res = fun(**self._params)
        return np.asarray(res, dtype=np.float64)[()]

    def moment(self, order):
        fun = np.vectorize(self._moment)
        order = self._make_mpf_array(order)
        res = fun(order, **self._params)
        return np.asarray(res, dtype=np.float64)[()]

class Normal(Distribution):

    # always override __init__ so IDEs hint at shape parameter names
    def __init__(self, *, mu, sigma):
        super().__init__(mu=mu, sigma=sigma)

    def _pdf(self, x, mu, sigma):
        return mp.npdf(x, mu, sigma)

class SkewNormal(Distribution):

    def __init__(self, *, a):
        super().__init__(a=a)

    def _pdf(self, x, a):
        return 2 * mp.npdf(x) * mp.ncdf(a * x)
