---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(distribution_infrastructure)=
```{eval-rst}
.. jupyterlite:: ../../_contents/hypothesis_bartlett.ipynb
   :new_tab: True
```

## Distribution Infrastructure (New)
_SciPy has two distribution infrastructures. This tutorial is for the (much) newer one, which has many structural
advantages, but fewer pre-defined distributions. For the old infrastructure, see :doc:`probability_distributions`._ 
### Basics

+++

Distributions (or, more accurately, distribution families) are classes named according to `CamelCase` conventions. They must be instantiated before use, with parameters passed as keyword-only arguments.
*Instances* of the distribution classes can be thought of as random variables, which are commonly denoted in mathematics using capital letters.

```{code-cell} ipython3
from scipy import stats
X = stats.LogUniform(a=1, b=2)
X
```

```{code-cell} ipython3
X.support()
```

Distributions can support multiple parameterizations, resolving requests like [gh-4538](https://github.com/scipy/scipy/issues/4538). For instance, it is also natural to parameterize the log-uniform distribution using the logarithms of the support endpoints. (If a log-uniform random variable is supported on $[a, b]$, its logarithm follows a uniform distribution with support $[\log(a), \log(b)$].)

```{code-cell} ipython3
import numpy as np
Y = stats.LogUniform(log_a=np.log(1), log_b=np.log(2))
Y
```

After being defined, these two random variables are essentially equivalent. As a weak example:

```{code-cell} ipython3
X.support() == Y.support()
```

All parameters of the distribution underlying the random variable are available as attributes.

```{code-cell} ipython3
X.a, X.b, X.log_a, X.log_b
```

Currently, distribution parameters are not intended to be changed, since the additional overhead of instantiating a new random variable is small. Nonetheless, support for modification of parameters is one of many planned enhancements.

+++

### Defining a distribution

+++

Minimal information is needed to fully define a distribution class. For example, a class representing a uniform distribution parameterized by the lower and upper ends of the support might look like this.

```{code-cell} ipython3
from scipy.stats._distribution_infrastructure import (ContinuousDistribution, _RealDomain,
                                                      _RealParameter, _Parameterization, oo)

class UniformDistribution(ContinuousDistribution):
    _a_param = _RealParameter('a', domain=_RealDomain(endpoints=(-oo, oo)))
    _b_param = _RealParameter('b', domain=_RealDomain(endpoints=('a', oo)))
    _x_param = _RealParameter('x', domain=_RealDomain(endpoints=('a', 'b'), inclusive=(True, True)))

    _parameterizations = [_Parameterization(_a_param, _b_param)]
    _variable = _x_param

    def _pdf_formula(self, x, *, a, b, **kwargs):
        return np.ones_like(x)/(b-a)
```

The infrastructure automatically validates numerical distribution parameters and method arguments based on their abstract definitions.

```{code-cell} ipython3
a, b = 1, 3
X = UniformDistribution(a=a, b=b)
```

```{code-cell} ipython3
X.support()
```

```{code-cell} ipython3
x = np.arange(a - 0.5, b + 0.51, 0.5)
x, X.pdf(x)
```

```{code-cell} ipython3
X.cdf(x)
```

```{code-cell} ipython3
X.icdf([-0.5, 0.5, 1.5])  # there are no numbers for which the CDF is negative or greater than 1
```

Above, note that the domain of the argument `x` was set to be inclusive, so the PDF *at* the limits of the support was `0.5` rather than `0.0`. On the other hand, the domains of parameters `a` and `b` are exclusive (by default). Rather than raising errors, out-of-domain shapes and NaNs result in methods returning NaNs. This allows for valid calculations to proceed normally.

```{code-cell} ipython3
X = UniformDistribution(a=[b, a, np.nan],
                        b=[a, b, np.nan])  # recall that the domain of b is (a, oo)
X.pdf(x[:, np.newaxis])
```

Besides input validation, the parameter information is used to draw numerical values of parameters for property-based tests:

```{code-cell} ipython3
UniformDistribution._a_param.typical = _RealDomain(endpoints=(-2, 0))
UniformDistribution._b_param.typical = _RealDomain(endpoints=(0, 2))
X = UniformDistribution._draw(sizes=[(3,), (2, 1)])
X.a, X.b
```

and to generate documentation:

```{code-cell} ipython3
from scipy.stats._distribution_infrastructure import _combine_docs
UniformDistribution.__doc__ = "The pdf is :math:`1 / (b-a)`..."
print(_combine_docs(UniformDistribution))
```

<a id='Transformations'></a>
### Transformations

Transformations can be applied to random variables. For instance, shifted and scaled versions can be created using `ShiftedScaledDistribution`.

```{code-cell} ipython3
from scipy.stats._distribution_infrastructure import ShiftedScaledDistribution
x = 1.
loc = np.asarray([1, 2, 3])
scale = np.asarray([2, 3])[:, np.newaxis]
X = stats.Normal()
Y = stats.ShiftedScaledDistribution(X, loc=loc, scale=scale)
np.testing.assert_equal(Y.cdf(x), X.cdf((x-loc)/scale))
```

```{code-cell} ipython3
Y.loc, Y.scale
```

For convenience, a `ShiftedScaledDistribution` can be created simply by performing an elementary arithmetic operation between a random variable and a numerical array.

```{code-cell} ipython3
Y = X*scale + loc
np.testing.assert_equal(Y.cdf(x), X.cdf((x-loc)/scale))
```

 There are several advantages of this architecture compared to building transformations directly into the `ContinuousDistribution` class:
- It allows distributions to use common parameterizations. By contrast, `rv_continuous` requires parameterizations to consider `loc` and `scale` or risk overparameterization (e.g. [gh-14716](https://github.com/scipy/scipy/issues/14716)). For example,
  - `stats.uniform` does not allow parameterization with the left and right support endpoints; it only accepts `loc` and `scale`.
  - `stats.loguniform` accepts the left and right support endpoints as shape parameters `a` and `b`; consequently, `a`, `b`, and `scale` are not independent parameters.
- Any overhead associated with a transformation is avoided unless the transformation is intentionally applied. (Although this is possible to achieve even if the transformation capabilities are built into the class, it may require special care.)
- It is highly extensible. For instance, transformations can also be used to generically define:
  - truncated distributions
  - half/double distributions
  - wrapped distributions
  - order statistic distributions
  - $\log$/$\exp$ transformed distributions

  and these transformations can be applied in any order.
- It avoids common pitfalls when fitting distributions to data. For instance, in the current infrastructure:
  - Users often forget to fix the location of distributions which almost always have fixed locations. This often results in poor fits or unexpected values of fit parameters.
  - It is impossible to fix the truncation points of truncated distributions because the loc-scale transformation is applied *after* the shape parameters truncate the support. It is more naturable to use the distribution if the these transformations are applied in the opposite order.

+++

Negative scale (multiplication by negative values) is supported. This eliminates the need to have separate left and right versions of some distributions (e.g. Weibull, Gumbel).

```{code-cell} ipython3
X = stats.LogUniform(a=1, b=2)
Y = stats.ShiftedScaledDistribution(X, loc=0, scale=-1)
X.support(), Y.support()
```

### Performance
#### Overhead
I've been careful to reduce overhead where possible.

```{code-cell} ipython3
x = 1.
X = stats.Normal()
%timeit X.pdf(x)
```

```{code-cell} ipython3
dist = stats.norm()  # old infrastructure
%timeit dist.pdf(x)
```

Even though these are meant to be instantiated once and used many times, instantiation followed by use is still tends to be faster than in the old infrastructure.

```{code-cell} ipython3
%timeit stats.Normal().pdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.norm.pdf(x)  # old infrastructure
```

If there's still too much overhead, the user can disable input validation.

```{code-cell} ipython3
X = stats.Normal(iv_policy='skip_all')
%timeit X.pdf(x)
```

```{code-cell} ipython3
%timeit stats.Normal(iv_policy='skip_all').pdf(x)
```

Overhead increases when shape parameters are invalid, need to be broadcast, or need to be converted to a floating point type for calculations. In these cases, there has been substantial effort to keep the overhead low and provide performance comparable to or better than `rv_continuous`.

+++

#### Numerical calculations
Another important aspect of performance is that of methods for which analytical formulas are not available. For example, the Gauss hypergeometric distribution can be defined as follows.

```{code-cell} ipython3
from scipy.stats._distribution_infrastructure import (ContinuousDistribution, _RealDomain,
                                                      _RealParameter, _Parameterization, oo)
from scipy import special

class GaussHyper(ContinuousDistribution):
    """Gauss hypergeometric distribution"""

    _a_param = _RealParameter('a', domain=_RealDomain(endpoints=(0, oo)))
    _b_param = _RealParameter('b', domain=_RealDomain(endpoints=(0, oo)))
    _c_param = _RealParameter('c', domain=_RealDomain(endpoints=(-oo, oo)))
    _z_param = _RealParameter('z', domain=_RealDomain(endpoints=(-1, oo)))
    _x_param = _RealParameter('x', domain=_RealDomain(endpoints=(0, 1), inclusive=(True, True)))

    _parameterizations = [_Parameterization(_a_param, _b_param, _c_param, _z_param)]
    _variable = _x_param

    def _pdf_formula(self, x, *, a, b, c, z, **kwargs):
        Cinv = special.gamma(a) * special.gamma(b) / special.gamma(a + b) * special.hyp2f1(c, a, a + b, -z)
        return 1.0 / Cinv * x ** (a - 1.0) * (1.0 - x) ** (b - 1.0) / (1.0 + z * x) ** c

a, b, c, z = 1.5, 2.5, 2, 0
X = GaussHyper(a=a, b=b, c=c, z=z)
x = 0.5
```

For scalar shapes and argument, performance of the new and old infrastructures are comparable.

```{code-cell} ipython3
%timeit X.cdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.gausshyper.cdf(x, a, b, c, z)  # old infrastructure
```

```{code-cell} ipython3
%timeit X.icdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.gausshyper.ppf(x, a, b, c, z)  # old infrastructure
```

```{code-cell} ipython3
np.testing.assert_allclose(X.cdf(x), stats.gausshyper.cdf(x, a, b, c, z))
np.testing.assert_allclose(X.icdf(x), stats.gausshyper.ppf(x, a, b, c, z))
```

But the quadrature and rootfinding code of the new infrastructure is vectorized (and eventually will be Array-API compatible), so it is much faster when arrays are involved.

```{code-cell} ipython3
x = np.linspace(0, 1, 1000)
```

```{code-cell} ipython3
%timeit X.cdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.gausshyper.cdf(x, a, b, c, z)  # old infrastructure
```

```{code-cell} ipython3
%timeit X.icdf(x)  # new infrastructure
```

```{code-cell} ipython3
# Warning: takes a long time
%timeit -r 1 -n 1 stats.gausshyper.ppf(x, a, b, c, z)  # old infrastructure
```

There are plans for the new infrastructure to use interpolation for additional performance gains with very large arrays.

+++

### Distribution properties
The new infrastructure has the distribution "properties" one would expect. `mode`, `skewness`, `kurtosis`, and `logentropy` are new.

```{code-cell} ipython3
X = stats.Normal()
X.mean(), X.median(), X.mode()
```

```{code-cell} ipython3
X.standard_deviation(), X.variance()
```

```{code-cell} ipython3
X.skewness(), X.kurtosis()  # *Pearson* kurtosis
```

```{code-cell} ipython3
X.entropy(), X.logentropy()
```

Note that the `logentropy` method returns a complex value because the entropy can be negative. The logarithm of a negative number is the logarithm of the number's magnitude plus an odd multiple of $\pi i$.

```{code-cell} ipython3
Y = stats.LogUniform(a=1, b=2)
Y.entropy(), Y.logentropy()
```

These are implemented as methods rather than `@property`s because they accept arguments. For instance, the entropy can be computed using the analytical formula, by exponentiating the log-entropy, or by quadrature.

```{code-cell} ipython3
X.entropy(), X.entropy(method='logexp'), X.entropy(method='quadrature')
```

### Distribution functions
Functions of the distributions underlying the random variables follow a consistent naming scheme.
- prefix `i` is for "inverse"
- prefix `c` is for "complementary"
- prefix `log` is for "logarithm of"

```{code-cell} ipython3
x = 1.
np.testing.assert_allclose(X.icdf(X.cdf(x)), x)
np.testing.assert_allclose(X.iccdf(X.ccdf(x)), x)
np.testing.assert_allclose(X.ilogcdf(X.logcdf(x)), x)
np.testing.assert_allclose(X.ilogccdf(X.logccdf(x)), x)
```

Note the addition of new methods for the inverse of the logarithm of distribution functions. These are useful when the argument of `icdf` would be too small or too close to `1.0` to represent accurately using floating point numbers.

```{code-cell} ipython3
np.testing.assert_allclose(X.ilogcdf(X.logcdf(-1000.)), -1000)
np.testing.assert_allclose(X.ilogccdf(X.logccdf(1000.)), 1000)
```

The distribution methods also have two-argument versions.

```{code-cell} ipython3
x1, x2 = 1., 2.
np.testing.assert_allclose(X.cdf(x1, x2),
                           X.cdf(x2) - X.cdf(x1))
np.testing.assert_allclose(X.ccdf(x1, x2),
                           1 - X.cdf(x1, x2))
np.testing.assert_allclose(X.logcdf(x1, x2),
                           np.log(X.cdf(x1, x2)))
np.testing.assert_allclose(X.logccdf(x1, x2),
                           np.log(X.ccdf(x1, x2)))
```

Besides convenience, this avoids catastropic cancellation where possible.

```{code-cell} ipython3
x1, x2 = 20., 20.5
X.cdf(20, 20.5), X.cdf(x2) - X.cdf(x1)
```

For numerically challenging cases, there are alternative `method` options available.

```{code-cell} ipython3
eps = 1e-100
res = X.logcdf(0., eps, method='quadrature')
ref = X.logpdf(0.) + np.log(eps)
np.testing.assert_equal(res, ref)
```

All distribution functions from the old distribution infrastructure are available in the new infrastructure (albeit under different names) with the following exceptions.
- `interval` is not available as a separate method, but the same values can be calculated using `iccdf` and `icdf`. However, the probability interval is in some sense an inverse of the two-argument `cdf`, so we could consider adding the capabilities to `icdf`.
- `expect` will not be supported. In the old infrastructure, this was little more than a light wrapper around an integrator, and we cannot do much better in general cases. The bug report to convenience ratio was too unfavorable to justify inclusion in the new infrastructure.

+++

### Random Sampling
Technically, "observe" might be a better name for this method, since instances like `X` represent a random variable. In any case, `sample` is easier to interpret than `rvs`:

```{code-cell} ipython3
X.sample()
```

Currently, a Generator can be passed either during construction or when calling the `sample` method.

```{code-cell} ipython3
rng = np.random.default_rng(872438745698345)
X = stats.Normal(rng=rng)
sample1 = X.sample()

rng2 = np.random.default_rng(872438745698345)
sample2 = X.sample(rng=rng2)

np.testing.assert_equal(sample1, sample2)
```

The parameter that controls the shape of the sample is called `shape`.

```{code-cell} ipython3
X.sample(shape=(2, 3))
```

`QMCEngine`s can also be used. Each slice along the last axis is generated from an independent low-discrepancy sequence. (*Note: currently, this is not the way it works, but that is what is slated to happen.*)

```{code-cell} ipython3
qrng = stats.qmc.Halton
n_observations = 10000
sample1 = X.sample(shape=(n_observations,), qmc_engine=qrng)
# Verify a property we would expect to hold exactly
np.testing.assert_equal((sample1 > 0).sum(), n_observations/2)
```

An important change is that the user does not need to consider the shape of the distribution parameters when specifying the `shape` of the sample. Instead, the shape of the output array is the specified `shape` concatenated with the distribution shape.

```{code-cell} ipython3
n_observations = 4
X_temp = stats.LogUniform(a=[0.5, 0.9],
                          b=[[1], [2], [3]])
sample = X_temp.sample(shape=n_observations)
sample.shape == (n_observations,) + X_temp._shape
```

### Moments

The `moment` method can compute raw, central, and standard moments of any order.

```{code-cell} ipython3
np.testing.assert_equal(X.moment(order=1, kind='raw'),
                        X.mean())
np.testing.assert_equal(X.moment(order=2, kind='central'),
                        X.variance())
np.testing.assert_equal(X.moment(order=3, kind='standardized'),
                        X.skewness())

X.moment(order=10, kind='standardized')
```

### Fitting
There is a draft of a generalized `fit` method. The method would unify techniques like maximum likelihood estimation with other needs, such as inverting distribution functions with respect to distribution parameters. We begin by initializing a normal distribution.

```{code-cell} ipython3
from scipy.stats._new_distributions import Normal
X = Normal(mu=-1, sigma=0.5)
```

Suppose we know the desired mean and standard deviation and wish to fit the `mu` and `sigma` parameters of the distribution to achieve them.

```{code-cell} ipython3
parameters = ['mu', 'sigma']
objective = {'f': lambda: [X.mean(), X.standard_deviation()],
             'output': [0.5, 1.5]}
X.fit(parameters, objective)
print(X.mean(), X.standard_deviation())
```

Or if we know the desired values of the `pdf` and `cdf` when the argument is `0`:

```{code-cell} ipython3
objective = dict(f=lambda x: [X.pdf(x), X.cdf(x)],
                 input=[0.],
                 output=[0.5, 0.35])
X.fit(parameters, objective)
X.pdf(0), X.cdf(0)
```

Of course, we can still perform maximum likelihood optimization.

```{code-cell} ipython3
data = X.sample(1000, rng=rng)
objective = dict(f=X.llf, input=(data,))
X.fit(parameters, objective)
X.mu, X.sigma
```

Currently, `fit` relies entirely on generic optimization procedures. In future work, the behavior can be overridden depending on the distribution, parameters, and objectives..

+++

### Visualization

We can visualize the results of the fit above using the convenience method, `plot`.

```{code-cell} ipython3
import matplotlib.pyplot as plt
ax = X.plot()
ax = plt.hist(data, density=True, alpha=0.5, bins=50)
```

```{code-cell} ipython3
X = stats.Normal(mu=[1, 2, 3], sigma=[1, 2, 3])
X.plot(y='cdf')
```

The `plot` method is relatively flexible, with a signature inspired by grammar of graphics. For instance, with the argument `x` on [-10, 10], plot the `pdf` against the `cdf`.

```{code-cell} ipython3
X.plot('cdf', 'pdf', t=('x', -10, 10))
```

### Order statistics distributions
There is draft support for distributions of [order statistics](https://en.wikipedia.org/wiki/Order_statistic) of distributions, partially to demonstrate the flexibility of distribution transformations. For example, we can plot the probability density functions of the order statistics of a normal distribution with sample size 4.

```{code-cell} ipython3
from scipy.stats._new_distributions import OrderStatisticDistribution
n = 4
r = np.arange(1, n+1)
X = stats.Normal()
Y = OrderStatisticDistribution(X, r=r, n=n)
Y.plot()
```

Compute the expected values of these order statistics.

```{code-cell} ipython3
Y.mean()
```

The `OrderStatisticDistribution` can be shifted and scaled, or we can generate an `OrderStatsticDistribution` from a shifted and scaled distribution. (In this case, the order of operations doesn't matter, but that is not the case for all transformations.)

```{code-cell} ipython3
loc, scale= 1, 2
Y1 = stats.ShiftedScaledDistribution(OrderStatisticDistribution(stats.Normal(), r=r, n=n), loc=loc, scale=scale)
Y2 = OrderStatisticDistribution(stats.ShiftedScaledDistribution(stats.Normal(), loc=loc, scale=scale), r=r, n=n)
np.testing.assert_allclose(Y1.mean(), Y.mean()*scale+loc)
np.testing.assert_allclose(Y2.mean(), Y.mean()*scale+loc)
```
