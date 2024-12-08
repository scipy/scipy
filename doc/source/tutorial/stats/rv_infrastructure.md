---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(rv_infrastructure)=
```{eval-rst}
.. jupyterlite:: ../../_contents/hypothesis_bartlett.ipynb
   :new_tab: True
```

## Random Variable Transition Guide

+++

### Background

+++

Prior to SciPy 1.15, all of SciPy's continuous probability distributions (e.g. `scipy.stats.norm`) have been instances of subclasses of `scipy.stats.rv_continuous`.

```{code-cell} ipython3
from scipy import stats
dist = stats.norm
type(dist)
```

```{code-cell} ipython3
isinstance(dist, stats.rv_continuous)
```

There were two obvious ways to work these objects.

According to the more common way, both arguments (e.g. `x`) and distribution parameters (e.g. `loc`, `scale`) were provided as inputs to methods of the object.

```{code-cell} ipython3
x, loc, scale = 1., 0., 1.
dist.pdf(x, loc, scale)
```

The less common approach was to invoke the `__call__` method of the distribution option, which returned an instance of `rv_continuous_frozen`, regardless of the original class.

```{code-cell} ipython3
frozen = stats.norm()
type(frozen)
```

Methods of this new object accept only arguments, but not the distribution parameters.

```{code-cell} ipython3
frozen.pdf(x)
```

In a sense, the instances of `rv_continuous` like `norm` represented "distribution families", which require parameters to identify a particular probability distribution, and an instance of `rv_continuous_frozen` was akin to a "random variable" - a mathematical object that follows a particular probability distribution.

Both approaches are valid and have advantages in certain situations. For instance, `stats.norm.pdf(x)` may appear more natural than `stats.norm().pdf(x)` for simple invocations. However, the former approach has a few inherent disadvantages; e.g., all of SciPy's 125 continuous distributions have to be instantiated at import time, and distribution parameters must be validated every time a method is called.

To address these and other shortcomings, gh-15928 proposed a new, separate infrastructure based on the latter (random variable) approach. This transition guide documents how users of `rv_continuous` and `rv_continuous_frozen` can migrate to the new infrastructure.

+++

### Basics

+++

In the new infrastructure, distributions families are classes named according to `CamelCase` conventions. They must be instantiated before use, with parameters passed as keyword-only arguments.
*Instances* of the distribution family classes can be thought of as random variables, which are commonly denoted in mathematics using capital letters.

```{code-cell} ipython3
from scipy import stats
X = stats.Normal()
X
```

```{code-cell} ipython3
X.pdf(x)
```

For simple calls like this (e.g. the argument is a valid float), call to methods of the new random variables will typically be faster than comparable calls to the old distribution methods.

```{code-cell} ipython3
%timeit dist.pdf(x)
```

```{code-cell} ipython3
%timeit frozen.pdf(x)
```

```{code-cell} ipython3
%timeit X.pdf(x)
```

Besides those above, many other calls to methods of the old distribution `norm` (represented as `dist`), a "frozen" distribution created by `norm` (represented as `frozen`), and an instance of the new `Normal` class (represented as `X`) appear similar. Calls to the following methods of `dist` and `norm` can be applied directly to `X`.

+++

- `pdf` (probability density function)
- `logpdf` (logarithm of probability density function)
- `cdf` (cumulative distribution function)
- `logcdf` (logarithm of cumulative distribution function)
- `entropy` (differential entropy)
- `median`
- `mean`
- `support`

+++

Others methods have new names, but are otherwise very similar.
- `sf` (survival function) $\rightarrow$ `ccdf` (complementary cumulative distribution function)
- `logsf` $\rightarrow$ `logccdf`
- `ppf` (percent point function) $\rightarrow$ `icdf` (inverse cumulative distribution function)
- `isf` (inverse survive function) $\rightarrow$ `iccdf` (inverse complementary cumulative distribution function)
- `std` $\rightarrow$ `standard_deviation`
- `var` $\rightarrow$ `variance`

+++

The new infrastructure has several new methods in the same vein as those above.

- `ilogcdf` (inverse of the logarithm of the cumulative distribution function)
- `ilogccdf` (inverse of the logarithm of the complementary cumulative distribution function)
- `logentropy` (logarithm of the entropy)
- `mode` (mode of the distribution)
- `skewness`
- `kurtosis` (*non-excess* kurtosis; see "Standardized Moments" below)

+++

And it has a new `plot` method for convenience

```{code-cell} ipython3
X.plot();
```

Most of the remaining methods of the old infrastructure (`rvs`, `moment`, `stats`, `interval`, `fit`, `nnlf`, `fit_loc_scale`, and `expect`) can be replaced, but some care is required. Before describing the replacements, we briefly mention how to work with random variables that are not distributed according to the standard normal (zero mean, unit standard deviation): almost all old distribution objects can be converted into a new distribution class with `scipy.stats.make_distribution`, and the new distribution class can be instantiated by passing the shape parameters as keyword arguments. For instance, consider the [Weibull distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weibull_min.html#scipy.stats.weibull_min). We can create a new class that is an abstraction of the distribution family like:

```{code-cell} ipython3
dist = stats.weibull_min
Weibull = stats.make_distribution(dist)
```

According to the documentation, the shape parameter of `weibull_min` was denoted `c`, so we can instantiate a random variable by passing `c` as a keyword argument.

```{code-cell} ipython3
c = 2.
X = Weibull(c=c)
X.plot();
```

Previously, all distributions inherited `loc` and `scale` parameters; now, random variables can be shifted and scaled with arithmetic operators.

```{code-cell} ipython3
Y = 2*X + 1
Y.plot();  # note the change in the abscissae
```

A separate distribution, `weibull_max`, was provided as the reflection of `weibull_min` about the origin. Now, this is just `-X`.

```{code-cell} ipython3
Y = -X
Y.plot();
```

#### Moments

+++

The previous infrastructure offered a `moments` method for raw moments. When only the order is specified, the new `moment` method is a drop-in replacement.

```{code-cell} ipython3
dist.moment(1, c), X.moment(1)  # first raw moment of the Weibull distribution with shape c
```

However, the previous infrastructure also had a `stats` method, which provided various statistics of the distribution. The following statistics were associated with the indicated characters.

+++

- Mean (first raw moment about the origin, `'m'`)
- Variance (second central moment, `'v'`)
- Skewness (third standarized moment, `'s'`)
- Excess kurtosis (fourh standarized moment minus 3, `'k'`)

For example:

```{code-cell} ipython3
dist.stats(c, moments='mvsk')
```

Now, moments of any `order` and `kind` (raw, central, and standardized) can be computed by passing the appropriate arguments to the new `moment` method.

```{code-cell} ipython3
X.moment(order=1, kind='raw'), X.moment(order=2, kind='central'), X.moment(order=3, kind='standardized'), X.moment(order=4, kind='standardized') 
```

Note the difference in definition of [kurtosis](https://en.wikipedia.org/wiki/Kurtosis). Previously, the "excess" (AKA "Fisher's") kurtosis was provided. As a matter of convention (rather than as a mathematical necessity), this is the standardized fourth moment shifted by a constant value (`3`) to give a value of `0.0` for the normal distribution.

```{code-cell} ipython3
stats.norm.stats(moments='k')
```

The new `moment` and `kurtosis` methods do not observe this convention by default; they report the standardized fourth moment according to the standard mathematical definition. This is also known as the "non-excess" (or "Pearson's") kurtosis, and the standard normal has a value of `3.0`.

```{code-cell} ipython3
stats.Normal().moment(4, kind='standardized'), stats.Normal().kurtosis()
```

For convenience, the `kurtosis` method offers a `convention` argument to select between the two.

```{code-cell} ipython3
stats.Normal().kurtosis(convention='non-excess'), stats.Normal().kurtosis(convention='excess') 
```

#### Random Variates

+++

The old `rvs` method has been replaced with `sample`, but there are two differences that should be noted.

First, `random_state` is replaced by `rng` per [SPEC 7](https://scientific-python.org/specs/spec-0007/). A pattern to control the random state in the past has been the use of `numpy.random.seed` or passing integers or instances of `numpy.RandomState`.

```{code-cell} ipython3
import numpy as np
np.random.seed(1)
dist = stats.norm
dist.rvs(), dist.rvs(random_state=1)
```

Now, use the `rng` argument with integers or instances of `numpy.Generator`.

```{code-cell} ipython3
X = stats.Normal()
X.sample(rng=1), X.sample(rng=np.random.default_rng(1))
```

Second, the argument `shape` replaces argument `size`.

```{code-cell} ipython3
dist.rvs(size=(3, 4))
```

```{code-cell} ipython3
X.sample(shape=(3, 4))
```

Besides the difference in name, there is a difference in behavior when array shape parameters are involved. Previously, the shape of distribution parameter arrays had to be included in `size`.

```{code-cell} ipython3
dist.rvs(size=(3, 4, 2), loc=[0, 1]).shape  # `loc` has shape (2,)
```

Now, the shape of the parameter arrays is considered to be a property of the random variable object itself. Specifying the shape of array shape parameters would be redundant, so it is not included when specifying the `shape` of the sample.

```{code-cell} ipython3
Y = stats.Normal(mu = [0, 1])
Y.sample(shape=(3, 4)).shape  # the sample has shape (3, 4); each element is of shape (2,)
```

#### Probability Mass (AKA "Confidence") Intervals

+++

The old infrastructure offered an `interval` method which, by default, would return a symmetric (about the median) interval containing a specified percentage of the distribution's probability mass.

```{code-cell} ipython3
a = 0.95
dist.interval(confidence=a)
```

Now, call the inverse CDF and complementary CDF methods with the desired probability mass.

```{code-cell} ipython3
p = 1 - a
X.icdf(p/2), X.iccdf(p/2)
```

#### Likelihood Function

+++

The old infrastructure offered a function to compute the *n*egative *l*og-*l*ikelihood *f*unction, erroneously named `nnlf` (instead of `nllf`). It accepts the parameters of the distribution as a tuple and observations as an array.

```{code-cell} ipython3
mu = 0
sigma = 1
data = stats.norm.rvs(size=100, loc=mu, scale=sigma)
stats.norm.nnlf((mu, sigma), data)
```

Now, simply compute the negative [log-likehood function](https://en.wikipedia.org/wiki/Likelihood_function) according to its mathematical definition.

```{code-cell} ipython3
X = stats.Normal(mu=mu, sigma=sigma)
-X.logpdf(data).sum()
```

#### Expected Values

+++

The `expect` method of the old infrastructure estimates a definite integral of an arbitrary function weighted by the PDF of the distribution. For instance, the fourth moment of the distribution is given by:

```{code-cell} ipython3
def f(x): return x**4
stats.norm.expect(f, lb=-np.inf, ub=np.inf)
```

This provides little added convencience over what the source code does, which is to use `scipy.integrate.quad` to perform the integration numerically.

```{code-cell} ipython3
from scipy import integrate
def f(x): return x**4 * stats.norm.pdf(x)
integrate.quad(f, a=-np.inf, b=np.inf)  # integral estimate, estimate of the error
```

Of course, this can just as easily be done with the new infrastructure.

```{code-cell} ipython3
def f(x): return x**4 * X.pdf(x)
integrate.quad(f, a=-np.inf, b=np.inf)  # integral estimate, estimate of the error
```

The `conditional` argument simply scales the result by inverse of the probability mass contained within the interval.

```{code-cell} ipython3
a, b = -1, 3
def f(x): return x**4
stats.norm.expect(f, lb=a, ub=b, conditional=True)
```

Using `quad` directly, that is:

```{code-cell} ipython3
def f(x): return x**4 * stats.norm.pdf(x)
prob = stats.norm.cdf(b) - stats.norm.cdf(a)
integrate.quad(f, a=a, b=b)[0]/prob
```

And with the new random variables:

```{code-cell} ipython3
integrate.quad(f, a=a, b=b)[0] / X.cdf(a, b)
```

Note that this is actually simplified because the `cdf` method of the new random variables accepts two arguments to compute the probability mass within the specified interval.

+++

#### Fitting

+++

The old infrastructure offered a function to estimate location and scale parameters of the distribution from data.

```{code-cell} ipython3
dist = stats.weibull_min
c, loc, scale = 0.5, 3., 4.
data = dist.rvs(size=100, c=c, loc=loc, scale=scale)
dist.fit_loc_scale(data, c)
```

Based on the source code, it is easy to replicate the method.

```{code-cell} ipython3
X = Weibull(c=c)
def fit_loc_scale(X, data):
    m, v = X.mean(), X.variance()
    m_, v_ = data.mean(), data.var()
    scale = np.sqrt(v_ / v)
    loc = m_ - scale*m
    return loc, scale

fit_loc_scale(X, data)
```

Note that the estimates are quite poor in this example, and poor performance of the heuristic factored into the decision not to provide a replacement.

The last method of the old infrastructure is `fit`, which estimates the location, scale, and shape parameters of an underlying distribution given a sample from the distribution.

```{code-cell} ipython3
params = stats.weibull_min.fit(data)
params
```

The convenience of this method is a blessing a curse. When it works, the simplicity is much appreciated. For some distributions, analytical expressions for the [maximum likelihood estimates](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) (MLE) of the parameters have been programmed manually, and for these distributions, the `fit` method is quite reliable. For most distributions, however, fitting is performed by numerical minimization of the negative log-likelihood function. This is not guaranteed to be successful - both because of inherent limitations of MLE and the state of the art of numerical optimization - and in these cases, the user of `fit` is stuck. However, a modicum of understanding goes a long way toward ensuring success, so here we present a mini-tutorial of fitting using the new infrastructure. 

First, note that MLE is one of several potential strategies for fitting distributions to data. It does nothing more than find the values of the distribution parameters for which the negative log-likelihood (NLLF) is minimized. Note that for the original data, the NLLF is:

```{code-cell} ipython3
stats.weibull_min.nnlf((c, loc, scale), data)
```

The NLLF of the parameters estimated using `fit` are lower:

```{code-cell} ipython3
stats.weibull_min.nnlf(params, data)
```

Therefore, `fit` considers its job done, whether or not this satisfies the user's notions of a "good fit" or whether the parameters are within reasonable bounds.

At its simplest, then, it is not so hard to replicate what `fit` does using tools available in SciPy. As discussed above, the NLLF is given by:

```{code-cell} ipython3
def nllf(params):
    c, loc, scale = params
    X = Weibull(c=c) * scale + loc
    return -X.logpdf(data).sum()

nllf(params)
```

To perform the minimization, we can use `scipy.optimize.minimize`.

```{code-cell} ipython3
from scipy import optimize
eps = 1e-10  # numerical tolerance to avoid invalid parameters
lb = [eps, -10, eps]  # lower bound on `c`, `location`, and `scale`
ub = [10, np.min(data)-eps, 10]  # upper bound on `c`, `location`, and `scale`
x0 = [1, 1, 1]  # guess to get optimization started
bounds = optimize.Bounds(lb, ub)
res = optimize.minimize(nllf, x0, bounds=bounds)
res
```

Although the value of `fun` (the NLLF) is not quite as low as with the solution provided by `weibull_min.fit`, the PDFs are not so different.

```{code-cell} ipython3
import matplotlib.pyplot as plt
x = np.linspace(2, 10, 300)

c, loc, scale = res.x
X = Weibull(c=c)*scale + loc
plt.plot(x, X.pdf(x), '-', label='numerical optimization')

c, loc, scale = params
Y = Weibull(c=c)*scale + loc
plt.plot(x, Y.pdf(x), '--', label='`scipy.stats.fit`')

plt.hist(data, bins=np.linspace(2, 10, 30), density=True, alpha=0.1)
plt.xlabel('x')
plt.ylabel('pdf(x)')
plt.legend()
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
