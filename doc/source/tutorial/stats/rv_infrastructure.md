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

The less common approach was to invoke the `__call__` method of the distribution object, which returned an instance of `rv_continuous_frozen`, regardless of the original class.

```{code-cell} ipython3
frozen = stats.norm()
type(frozen)
```

Methods of this new object accept only arguments, not the distribution parameters.

```{code-cell} ipython3
frozen.pdf(x)
```

In a sense, the instances of `rv_continuous` like `norm` represented "distribution families", which require parameters to identify a particular probability distribution, and an instance of `rv_continuous_frozen` was akin to a "random variable" - a mathematical object that follows a particular probability distribution.

Both approaches are valid and have advantages in certain situations. For instance, `stats.norm.pdf(x)` may appear more natural than `stats.norm().pdf(x)` for simple invocations. However, the former approach has a few inherent disadvantages; e.g., all of SciPy's 125 continuous distributions have to be instantiated at import time, distribution parameters must be validated every time a method is called, and documentation of methods must either a) be generated separately for every method of every distribution or b) omit the shape parameters that are unique for each distribution.

To address these and other shortcomings, gh-15928 proposed a new, separate infrastructure based on the latter (random variable) approach. This transition guide documents how users of `rv_continuous` and `rv_continuous_frozen` can migrate to the new infrastructure.

+++

### Basics

+++

In the new infrastructure, distributions families are classes named according to `CamelCase` conventions. They must be instantiated before use, with parameters passed as keyword-only arguments.
*Instances* of the distribution family classes can be thought of as random variables, which are commonly denoted in mathematics using capital letters.

```{code-cell} ipython3
from scipy import stats
X = stats.Normal(mu=0, sigma=1)
X
```

Once instantiated, shape parameters can be read (but not written) as attributes.

```{code-cell} ipython3
X.mu, X.sigma
```

Note that the documentation of [`scipy.stats.Normal`](https://scipy.github.io/devdocs/reference/generated/scipy.stats.Normal.html#scipy.stats.Normal) contains links to more detailed documentation of each of its methods. (Compare, for instance, against the documentation of [`scipy.stats.norm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html).)

```{code-cell} ipython3
X.pdf(x)
```

For simple calls like this (e.g. the argument is a valid float), calls to methods of the new random variables will typically be faster than comparable calls to the old distribution methods.

```{code-cell} ipython3
%timeit dist.pdf(x, loc=mu, scale=sigma)
```

```{code-cell} ipython3
%timeit frozen.pdf(x)
```

```{code-cell} ipython3
%timeit X.pdf(x)
```

Note that the calls to `frozen.pdf` and `X.pdf` are identical, and the call to `dist.pdf` is very similar to the call to `X.pdf` - the only difference is that the call to `dist.pdf` includes the shape parameters, whereas in the new infrastructure, shape parameters are only provided when the random variable is instantiated.

Besides `pdf`, several other methods of the new infrastructure are essentially the same as the old methods.

+++

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

The new infrastructure has several *new* methods in the same vein as those above.

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

Most of the remaining methods of the old infrastructure (`rvs`, `moment`, `stats`, `interval`, `fit`, `nnlf`, `fit_loc_scale`, and `expect`) can be replaced, but some care is required. Before describing the replacements, we briefly mention how to work with random variables that are not normally distributed: almost all old distribution objects can be converted into a new distribution class with `scipy.stats.make_distribution`, and the new distribution class can be instantiated by passing the shape parameters as keyword arguments. For instance, consider the [Weibull distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weibull_min.html#scipy.stats.weibull_min). We can create a new class that is an abstraction of the distribution family like:

```{code-cell} ipython3
Weibull = stats.make_distribution(stats.weibull_min)
```

According to the documentation of [`weibull_min`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weibull_min.html#scipy.stats.weibull_min), the shape parameter is denoted `c`, so we can instantiate a random variable by passing `c` as a keyword argument to the new `Weibull` class.

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
stats.weibull_min.moment(1, c), X.moment(1)  # first raw moment of the Weibull distribution with shape c
```

However, the previous infrastructure also had a `stats` method, which provided various statistics of the distribution. The following statistics were associated with the indicated characters.

+++

- Mean (first raw moment about the origin, `'m'`)
- Variance (second central moment, `'v'`)
- Skewness (third standarized moment, `'s'`)
- Excess kurtosis (fourh standarized moment minus 3, `'k'`)

For example:

```{code-cell} ipython3
stats.weibull_min.stats(c, moments='mvsk')
```

Now, moments of any `order` and `kind` (raw, central, and standardized) can be computed by passing the appropriate arguments to the new `moment` method.

```{code-cell} ipython3
X.moment(order=1, kind='raw'), X.moment(order=2, kind='central'), X.moment(order=3, kind='standardized'), X.moment(order=4, kind='standardized') 
```

Note the difference in definition of [kurtosis](https://en.wikipedia.org/wiki/Kurtosis). Previously, the "excess" (AKA "Fisher's") kurtosis was provided. As a matter of convention (rather than the natural mathematical definition), this is the standardized fourth moment shifted by a constant value (`3`) to give a value of `0.0` for the normal distribution.

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

The old `rvs` method has been replaced with `sample`, but there are a two differences that should be noted.

First, `random_state` is replaced by `rng` per [SPEC 7](https://scientific-python.org/specs/spec-0007/). A pattern to control the random state in the past has been to use `numpy.random.seed` or to pass integer seeds to the `random_state` argument. The integer seeds were converted to instances of `numpy.random.RandomState`, so behavior for a given integer seed would be identical in these two cases:

```{code-cell} ipython3
import numpy as np
np.random.seed(1)
dist = stats.norm
dist.rvs(), dist.rvs(random_state=1)
```

In the new infrastructure, pass instances of `numpy.Generator` to the `rng` argument of `sample` to control the random state. Note that integer seeds are now converted to instances of `numpy.random.Generator`, not instances of `numpy.random.RandomState`.

```{code-cell} ipython3
X = stats.Normal()
rng = np.random.default_rng(1)  # instantiate a numpy.random.Generator
X.sample(rng=rng), X.sample(rng=1)
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

```{code-cell} ipython3

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

This provides little added convencience over what the source code does: use `scipy.integrate.quad` to perform the integration numerically.

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

Note that the `cdf` method of the new random variables accepts two arguments to compute the probability mass within the specified interval. In many cases, this can avoid subtractive error ("catastrophic cancellation") when the probability mass is very small.

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

## New Features

+++

<a id='Transformations'></a>
### Transformations

Mathematical transformations can be applied to random variables. For instance, many elementary arithmetic operations (`+`, `-`, `*`, `/`, `**`) between random variables and scalars work.

For instance, we can see that the reciprocal of a Weibull-distributed RV is distributed according to `scipy.stats.invweibull`. (An "inverse" distribution is typically the distribution underlying the reciprocal of a random variable.)

```{code-cell} ipython3
c = 10.6

X = Weibull(c=10.6)  
Y = 1 / X  # compare to `invweibull`
Y.plot();

x = np.linspace(0.8, 2.05, 300)
plt.plot(x, stats.invweibull(c=c).pdf(x), '--')
```

`scipy.stats.chis2` describes a sum of the squares of normally-distributed random variables.

```{code-cell} ipython3
X = stats.Normal()
Y = X**2  # compare to chi2
Y.plot(t=('x', 0, 5));

x = np.linspace(0, 5, 300)
plt.plot(x, stats.chi2(df=1).pdf(x), '--')
```

`scipy.stats.foldcauchy` describes the absolute value of a Cauchy-distributed random variable. (A "folded" distribution is the distribution underlying the absolute value of a random variable.)

```{code-cell} ipython3
Cauchy = stats.make_distribution(stats.cauchy)
c = 4.72

X = Cauchy() + c  
Y = abs(X)  # compare to `foldcauchy`
Y.plot(t=('x', 0, 60));

x = np.linspace(0, 60, 300)
plt.plot(x, stats.foldcauchy(c=c).pdf(x), '--')
```

`scipy.stats.lognormal` describes the exponential of a normally distributed random variable. It is so-named because the logarithm of the resulting random variable is normally distributed. (In general, a "log" distribution is typically the distribution underlying the *exponential* of a random variable.)

```{code-cell} ipython3
u, s = 1, 0.5

X = stats.Normal()*s + u
Y = stats.exp(X)  # compare to `lognorm`
Y.plot(t=('x', 0, 9));

x = np.linspace(0, 9, 300)
plt.plot(x, stats.lognorm(s=s, scale=np.exp(u)).pdf(x), '--')
```

`scipy.stats.loggamma` describes the logarithm of of a Gamma-distributed random variable. Note that the more common name of this distribution would be [exp-gamma](https://reference.wolfram.com/language/ref/ExpGammaDistribution.html.en#:~:text=The%20term%20exp%2Dgamma%20is,of%20qualitatively%20similar%20probability%20distributions.) because the exponential of the RV is gamma-distributed.

```{code-cell} ipython3
Gamma = stats.make_distribution(stats.gamma)
a = 0.414

X = Gamma(a=a)  
Y = stats.log(X)  # compare to `loggamma`
Y.plot();

x = np.linspace(-17.5, 2, 300)
plt.plot(x, stats.loggamma(c=a).pdf(x), '--')
```

`scipy.stats.truncnorm` is the distribution underlying a truncated normal random variable. Note that the truncation transformation can be applied either before or after shifting and scaling the normally-distributed random variable, which can make it much easier to achieve the desired result than `scipy.stats.truncnorm` (which inherently shifts and scales *after* truncation).

```{code-cell} ipython3
a, b = 0.1, 2

X = stats.Normal()  
Y = stats.truncate(X, a, b)  # compare to `truncnorm`
Y.plot();

x = np.linspace(a, b, 300)
plt.plot(x, stats.truncnorm(a, b).pdf(x), '--')
```

`scipy.stats.dgamma` is a mixture of two gamma-distributed RVs, one of which is reflected about the origin. (Typically, a "double" distribution is the distribution underlying a mixture of RVs, one of which is reflected.)

```{code-cell} ipython3
a = 1.1
X = Gamma(a=a)
Y = stats.Mixture((X, -X), weights=[0.5, 0.5])
# plot method not available for mixtures

x = np.linspace(-4, 4, 300)
plt.plot(x, Y.pdf(x))
plt.plot(x, stats.dgamma(a=a).pdf(x), '--')
```

*Limitations*:

While most arithmetic transformations between random variables and Python scalars or NumPy arrays are supported, there are a few restrictions.

- Raising a random variable to the power of an argument is only implemented when the argument is a positive integer.
- Raising an argument to the power of a random variable is only implemented when the argument is a positive scalar other than 1.
- Division by a random variable is only implemented when the support is either non-negative or non-positive.
- The logarithm of a random variable is only implemented when the support is non-negative. (The logarithm of negative values is imaginary.)

Arithmetic operations between two random variables are *not* yet supported. Note that such operations are much more complex mathematically; for instance, the PDF of the sum of the two random variables involves convolution of the two PDFs.

Also, while transformations are composable, a) truncation and b) shifting/scaling can each be done only once. For instance, `Y = 2 * stats.exp(X + 1)` will raise an error because this would require shifting before exponentiation and scaling after exponentiation; this is treated as "shifting/scaling" twice. However, a mathematical simplification is possible here (and in many cases) to avoid the problem: `Y = (2*math.exp(2)) * stats.exp(X)` is equivalent and requires only one scaling operation.

Although the transformations are fairly robust, they all reply on generic implementations which may cause numerical difficulties. If you are concerned about the accuracy of the results of a transformation, consider comparing the resulting PDF against a histogram of a random sample.

+++

X = stats.Normal()
Y = X**3
x = np.linspace(-5, 5, 300)
plt.plot(x, Y.pdf(x), label='pdf')
plt.hist(X.sample(100000)**3, density=True, bins=np.linspace(-5, 5, 100), alpha=0.5);
plt.ylim(0, 2)

+++

### Quasi-Monte Carlo Sampling

+++

Random variables enable generation of quasi-random, low-discrepancy samples from statistical distributions.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
X = stats.Normal()

rng = np.random.default_rng(7824387278234)
qrng = stats.qmc.Sobol(1, rng=rng)  # instantiate a QMCEngine

bins = np.linspace(-3.5, 3.5, 31)
plt.hist(X.sample(512, rng=qrng), bins, alpha=0.5, label='quasi-random')
plt.hist(X.sample(512, rng=rng), bins, alpha=0.5, label='pseudo-random')
plt.title('Histogram of normally-distributed sample')
plt.legend()
plt.show()
```

Note that when generating multiple samples (e.g. `len(shape) > 1`),  each slice along the zeroth axis is an independent, low-discrepancy sequency. The following generates two independent QMC samples, each of length 512.

```{code-cell} ipython3
samples = X.sample((512, 2), rng=qrng)
plt.hist(samples[:, 0], bins, alpha=0.5, label='sample 0')
plt.hist(samples[:, 1], bins, alpha=0.5, label='sample 1')
plt.title('Histograms of normally-distributed samples')
plt.legend()
plt.show()
```

The result is quite different if we generate 512 independent QMC sequences, each of length 2.

```{code-cell} ipython3
samples = X.sample((2, 512), rng=qrng)
plt.hist(samples[0], bins, alpha=0.5, label='sample 0')
plt.hist(samples[1], bins, alpha=0.5, label='sample 1')
plt.title('Histograms of normally-distributed samples')
plt.legend()
plt.show()
```

### Accuracy

+++

For some distributions, like `scipy.stats.norm`, almost all methods of the distribution are customized to ensure accurate computation. For others, like `scipy.stats.gausshyper`, little more than the PDF is defined, and the other methods must be computed numerically based on the PDF. For example, the survival function (complementary CDF) of the Gauss hypergeometric distribution is calculated by numerically integrating the PDF from `0` (the left end of the support) to `x` to get the CDF, then the complement is taken.

```{code-cell} ipython3
a, b, c, z = 1.5, 2.5, 2, 0
x = 0.5
frozen = stats.gausshyper(a=a, b=b, c=c, z=z)
frozen.sf(x)
```

```{code-cell} ipython3
frozen.sf(x) == 1 - integrate.quad(frozen.pdf, 0, x)[0]
```

However, another aproach would be to numerically integrate the PDF from `x` to `1` (the right end of the support).

```{code-cell} ipython3
integrate.quad(frozen.pdf, x, 1)[0]
```

These are distince but equally valid approaches, so assuming the PDF is accurate, it is unlikely that the two results are inaccurate in the same way. Therefore, we can estimate the accuracy of the results by comparing them.

```{code-cell} ipython3
res1 = frozen.sf(x)
res2 = integrate.quad(frozen.pdf, x, 1)[0]
abs((res1 - res2) / res1)
```

The new infrastructure is aware of several different ways for computing most quantities. For example, this diagram illustrates the relationships between various distribution functions.

![image.png](attachment:c7742b7b-e07f-4bc4-9e05-3566bc7bde2b.png)

It follows a decision tree to choose what is expected to be the most accurate way of estimating the quantity. These decision trees are subject to change, but an example for computing the complementary CDF might look something like:

![image.png](attachment:ee923abe-2ba7-43ec-815a-04e69e6f2b73.png)

and for calculating a moment:

![image.png](attachment:68d09366-2cca-477c-ad70-72541b64b29a.png)

However, you can override the method it uses with the `method` argument.

```{code-cell} ipython3
GaussHyper = stats.make_distribution(stats.gausshyper)
X = GaussHyper(a=a, b=b, c=c, z=z)
```

```{code-cell} ipython3
X.ccdf(x, method='quadrature') == integrate.tanhsinh(X.pdf, x, 1).integral
```

```{code-cell} ipython3
X.ccdf(x, method='complement') == (1 - integrate.tanhsinh(X.pdf, 0, x).integral)
```

```{code-cell} ipython3
X.ccdf(x, method='logexp') == np.exp(integrate.tanhsinh(lambda x: X.logpdf(x), x, 1, log=True).integral)
```

When in doubt, consider trying different `method`s to compute a given quantity.

+++

### Performance

+++

Consider the performance of calculations involving the Gauss hypergeometric distribution. For scalar shapes and argument, performance of the new and old infrastructures are comparable. The new infrastructure is not expected to be much faster; although it reduces the "overhead" of parameter validation, it is not significantly faster at numerical integration of scalar quantities, and the latter dominates the execution time here.

```{code-cell} ipython3
%timeit X.cdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.gausshyper.cdf(x, a, b, c, z)  # old infrastructure
```

```{code-cell} ipython3
np.isclose(X.cdf(x), stats.gausshyper.cdf(x, a, b, c, z))
```

But the new infrastructure is much faster when arrays are involved. This is because the underlying integrator ([`scipy.integrate.tanhsinh`](https://scipy.github.io/devdocs/reference/generated/scipy.integrate.tanhsinh.html#scipy.integrate.tanhsinh)) and root finder ([`scipy.optimize.elementwise.find_root`](https://scipy.github.io/devdocs/reference/generated/scipy.optimize.elementwise.find_root.html)) of the new infrastructure are natively vectorized, whereas the old routines (`scipy.integrate.quad` and `scipy.optimize.brentq`) are not.

```{code-cell} ipython3
x = np.linspace(0, 1, 1000)
```

```{code-cell} ipython3
%timeit X.cdf(x)  # new infrastructure
```

```{code-cell} ipython3
%timeit stats.gausshyper.cdf(x, a, b, c, z)  # old infrastructure
```

Comparable performance improvements can be expected for array calculations whenever the underlying functions are computed by either numerical quadrature or inversion (root-finding).

+++

### Visualization

We can readily visualize functions of the distribution underlying random variables using the convenience method, `plot`.

```{code-cell} ipython3
import matplotlib.pyplot as plt
ax = X.plot()
```

```{code-cell} ipython3
X.plot(y='cdf')
```

The `plot` method is quite flexible, with a signature inspired by grammar of graphics. For instance, with the argument `x` on [-10, 10], plot the `pdf` against the `cdf`.

```{code-cell} ipython3
X.plot('cdf', 'pdf', t=('x', -10, 10))
```

### Order statistics distributions
There is support for distributions of [order statistics](https://en.wikipedia.org/wiki/Order_statistic) of random samples from distribution. For example, we can plot the probability density functions of the order statistics of a normal distribution with sample size 4.

```{code-cell} ipython3
n = 4
r = np.arange(1, n+1)
X = stats.Normal()
Y = stats.order_statistic(X, r=r, n=n)
Y.plot()
```

Compute the expected values of these order statistics:

```{code-cell} ipython3
Y.mean()
```
