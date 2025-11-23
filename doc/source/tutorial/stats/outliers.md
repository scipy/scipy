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

+++ {"tags": ["jupyterlite_sphinx_strip"]}

```{eval-rst}
.. notebooklite:: outliers.md
   :new_tab: True
```

(outliers)=

+++

# Trimming and Winsorization Transition Guide

An outlier is an observation that differs substantially from other observations in a sample. There are many references available about identifying outliers, whether to make adjustments for outliers, and if so, what adjustments to make. This guide does not make any recommendations on these topics. Instead, the focus is how to perform two of the most common procedures:

- trimming: removing outliers from the dataset, and
- winsorization: replacing outliers with a more central value,

using tools available in SciPy, with an emphasis on transitioning from use of legacy functions to the preferred approach.

## Motivation
In the past, SciPy has offered several "convenience functions" that combine trimming with computation of a statistic. Consider, for instance, `scipy.stats.trim_mean`.

```{code-cell} ipython3
import numpy as np
from scipy import stats
np.set_printoptions(linewidth=120)
x = np.arange(20.)
stats.trim_mean(x, proportiontocut=0.1)
```

Without carefully reading the documentation, it is not immediately obvious how the operation was performed. What would be your guess: does `proportiontocut` refer to the proportion to cut from each tail, or is it the total proportion of the data to remove?

In this case, we would have gotten the same result either way:  < Albert - shall I just remove all the `np.mean` calculations here? I'm thinking they don't add anything after all.>

```{code-cell} ipython3
# Remove 2/20 observations total before taking the `mean`, one in each tail
np.mean(x[1:-1])
```

```{code-cell} ipython3
# Remove 2/20 observations from *each* tail before taking the `mean`
np.mean(x[2:-2])
```

Likewise, when the proportion would not result in an integral number of points being removed, how many points are removed? 

For example, in an array with 22 elements, is the number to remove rounded up or down?

```{code-cell} ipython3
x = np.arange(22.)
# Remove `floor(0.1 * 22)` (2) observations from each tail
np.mean(x[2:-2])
```

```{code-cell} ipython3
# Remove `ceil(0.1 * 22)` (2) observations from each tail
np.mean(x[3:-3])
```

What if instead of removing a proportion, we want to remove data below and above certain thresholds. It turns out that we would need to use `scipy.stats.tmean` or (the similarly-named, but curiously relegated to a separate namespace) `scipy.stats.mstats.trimmed_mean`, instead.

```{code-cell} ipython3
stats.tmean(x, limits=(2, 19), inclusive=(True, True))
```

```{code-cell} ipython3
stats.mstats.trimmed_mean(x, limits=(2, 19), relative=False, inclusive=(True, True))
```

What if we want to remove data outside a "confidence interval" estimate. Do any of these functions offer anything to facilitate that?

No, not really.

One could take the mean after using one of the seven functions for trimming data. But which one(s) will do what we need?

- `stats.mstats.trim`
- `stats.mstats.trima`
- `stats.mstats.trimr`
- `stats.mstats.trimtail`
- `stats.mstats.trimboth`
- `stats.trimboth`
- `stats.trim1`

So at the time of writing, SciPy offers at least half a dozen ways of taking a trimmed mean, but none are able to perform *all* variations of the trimmed mean that we might reasonably expect them to.

What about other statistics - should we use `stats.tvar`, `stats.mstats.tvar`, or `stats.mstats.trimmed_var` to compute the trimmed variance? If there are three functions for that and three for the mean, surely there is at least one function for computing the trimmed skewness? (Naw!)

In the spirit of the Zen of Python, we suggest the new "one-- and preferably only one --obvious way" to calculate the trimmed version of many statistics: manually, in two or three lines of code, using more familiar, general-purpose features of NumPy and SciPy. Whether calculating a trimmed or winsorized statistic, there will be three steps:

1. Identify the outlier threshold.
2. Trim or winsorize.
3. Compute the statistic.

+++

## Identifying the threshold

SciPy now offers a `quantile` function with several methods useful for identifying thresholds beyond which data are to be trimmed or winsorized. The most common thresholds are those corresponding with a certain *percentage* of the data in each tail. To set the threshold at $p=10\%$ of the data in each tail, rounding the number of points to trim *down*:

```{code-cell} ipython3
x = np.arange(22.)
p = 0.1
stats.quantile(x, [p, 1-p], method='winsor_less')
```

Rounding the number of points to trim *up*:

```{code-cell} ipython3
stats.quantile(x, [p, 1-p], method='winsor_more')
```

Rounding the number of points to trim to the nearest integer:

```{code-cell} ipython3
stats.quantile(x, [p, 1-p], method='winsor_round')
```

Note that these methods are all subtlely different from estimating the 10th and 90th percentiles of the data, which might be called for in different circumstances. Using the default method:

```{code-cell} ipython3
stats.quantile(x, [p, 1-p])
```

The Harrell-Davis method:

```{code-cell} ipython3
stats.quantile(x, [p, 1-p], method='harrell-davis')
```

Note that `quantile` is fully vectorized. Suppose we have two rows of data, and we wish to set different thresholds for each row.

```{code-cell} ipython3
x2d = np.stack([x, x])
p2d = [[p  , 1 - p  ],  # thresholds for zeroth row
       [2*p, 1 - 2*p]]  # thresholds for next row
stats.quantile(x2d, p2d, method='winsor_less', axis=-1)
```

### Missing Data

+++

When there are missing or ignorable values but omitting them would lead to a ragged array, the two historical options have been to:

- cover the missing/ignorable values with the "mask" of a NumPy masked array and use `stats.mstats` functions for analysis, or
- place NaNs in place of the missing/ignorable values and use `stats` functions with `nan_poliy='omit'` for analysis.

The downsides of these approaches are documented elsewhere; here, we recommend using an [MArray](https://mdhaber.github.io/marray/tutorial.html) with `stats` functions.

```{code-cell} ipython3
from marray import numpy as mxp
mask = ((x2d % 2) == 0)  # ignore odd elements
y = mxp.asarray(x2d, mask=mask)
print(y)
```

`quantile` and many other `scipy.stats` functions natively support `MArray`s, and support grows with each version of SciPy.

```{code-cell} ipython3
stats.quantile(y, mxp.asarray(p2d), method='winsor_less', axis=-1)
```

## Trimming and Winsorizing

+++

### Winsorizing is `clip`ing

+++

Once a numerical threshold has been identified, winsorizing is identical to a routine so fundamental that it is part of the array API standard: `clip`.

Suppose we have $n=22$ observations and we wish to winsorize $p=10\%$ from each side of the data. Rounding down with `method='winsor_less'`, we would expect to winsorize $\left \lfloor{pn}\right \rfloor = 2$ points on each side. In this case, `quantile` returns the *third* data point from the left (at index `2`) as the lower limit. `clip`ing to this limit will winsorize the two lower observations, leaving the third observation unchanged. The story is similar on the right.

```{code-cell} ipython3
x = np.arange(22.)
p = 0.1
low, high = stats.quantile(x, [p, 1-p], method='winsor_less')
x_winsorized = np.clip(x, low, high)
x_winsorized
```

When extracting the `low` and `high` thresholds from the result of quantile, just remember that their positions correspond with the positions of values in the `p` argument.

```{code-cell} ipython3
stats.quantile(x2d, p2d, method='winsor_less', axis=-1)
```

Here, we want the columns, and we want to preserve their shape so they still align with the shape of `x2d`.

```{code-cell} ipython3
low_high = stats.quantile(x2d, p2d, method='winsor_less', axis=-1)
low2d = low_high[:, :1]
high2d = low_high[:, 1:]
low2d, high2d
```

```{code-cell} ipython3
x2d_winsorized = np.clip(x2d, low2d, high2d)
x2d_winsorized
```

### Trimming

+++

For 1-D data, trimming is equivalent to another fundamental routine: indexing.

In the general case that the data are unsorted and threshold values are known, use a boolean mask to pick out the data to be kept.

As in the 1-D winsorization example above, suppose we want to trim two points on either side. `low` is the *third* observation from the left - the first we want to keep. Therefore, we'll form a boolean mask with *non-strict* inequalities to select the data to retain.

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_less')  # get the threshold values
mask = (x >= low) & (x <= high)
x_trimmed = x[mask]
x_trimmed
```

To keep code easier to interpret, use `method='winsor_more'` if you want to round the other way (rather than changing to strict inequalities).

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_more')  # get the threshold values
mask = (x >= low) & (x <= high)
x[mask]
```

In 2-D, this approach of retaining only certain data *may* be possible if the number of elements retained in each row happen to be identical. Note that boolean indexing results in a 1-D array:

```{code-cell} ipython3
mask = (x2d >= low) & (x2d <= high)
x2d[mask]
```

*When the number of elements retained in each row are identical (and only then)*, it is safe to reshape the output to an array with the original number of rows.

```{code-cell} ipython3
x2d[mask].reshape(x2d.shape[0], -1)  # -1 means "infer from the number of elements"
```

The more general option is to keep the trimmed elements in their original locations but to *mask* them.

Note that the mask of elements to *remove* is the logical NOT of the elements to be kept. Since we used *non-strict* inequalities to retain elements, we use *strict* inequalities to trim them.

```{code-cell} ipython3
mask2d = (x2d < low2d) | (x2d > high2d)
x2d_masked = mxp.asarray(x2d, mask=mask2d)
print(x2d_masked)
```

## Computing the statistic

+++

When the data is winsorized or trimmed data has actually been removed, working with it is no different than any other data.

```{code-cell} ipython3
np.mean(x_trimmed)
```

```{code-cell} ipython3
np.mean(x_winsorized)  # result happens to be identical *for this data*
```

```{code-cell} ipython3
np.mean(x2d_winsorized, axis=-1)  # results happen to be identical *for this data*
```

When the trimmed values are represented as masked elements of a masked array, only slightly more care is needed. For simple statistics like `mean` and `var`, use the corresponding function of the MArray namespace, `mxp`, rather than the NumPy namespace.

```{code-cell} ipython3
mxp.mean(x2d_masked, axis=-1)
```

Many SciPy functions also natively support masked arrays (see the documentation of the desired function).

```{code-cell} ipython3
stats.kurtosis(x2d_masked, axis=-1)
```

## Recipes

+++

Here, we reproduce the results of several legacy trimming functions using the preferred approach.

```{code-cell} ipython3
rng = np.random.default_rng(339759635264198322188079188741355346875)
x = rng.random(22)
p = 0.1
```

### `stats.trim_mean`

+++

`trim_mean` accepts a proportion of the total number of elements to cut from *each* side.

```{code-cell} ipython3
stats.trim_mean(x, proportiontocut=p)
```

It always rounds the number of elements to cut *down*, so we choose `winsor_less`.

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_less')
mask = (low <= x) & (x <= high)
np.mean(x[mask])
```

### `stats.tmean`

+++

`tmean` accepts absolute lower and upper `limits`.

```{code-cell} ipython3
limits = (0.1, 0.9)
stats.tmean(x, limits, inclusive=(True, True))
```

The `inclusive` argument has no effect unless one of the limits exactly equals an observation.

```{code-cell} ipython3
stats.tmean(x, limits, inclusive=(False, False))
```

As an example of when it would matter:

```{code-cell} ipython3
limits = (x[9], x[5])
tmean_inclusive = stats.tmean(x, limits, inclusive=(True, True))
tmean_exclusive = stats.tmean(x, limits, inclusive=(False, False))
tmean_inclusive, tmean_exclusive
```

`inclusive` refers to the elements to *keep*, so use *non-strict* inequality to be inclusive:

```{code-cell} ipython3
mask = (limits[0] <= x) & (x <= limits[1])
np.mean(x[mask])
```

and *strict* inequality to be exclusive.

```{code-cell} ipython3
mask = (limits[0] < x) & (x < limits[1])
np.mean(x[mask])
```

### `stats.mstats.trimmed_mean`

+++

`trimmed_mean` is, in some ways, a confusing mash-up of `stats.tmean` and `stats.trim_mean`. It accepts `limits`, like `stats.tmean`, but by default, these are relative, like `stats.trim_mean`.

```{code-cell} ipython3
stats.mstats.trimmed_mean(x, (p, p), inclusive=(True, True), relative=True)
```

Surprisingly, `inclusive` has no effect when `p*len(x)` is integral and `relative` is `True` (default).

```{code-cell} ipython3
stats.mstats.trimmed_mean(x, (p, p), inclusive=(False, False))
```

Rather, it affects the method of rounding `p*len(x)`. So let's consider a case in which it would matter.

```{code-cell} ipython3
rng = np.random.default_rng(324589234958293486245)
x = rng.random(size=19)
stats.mstats.trimmed_mean(x, (p, p), inclusive=(True, True))
```

According to the documentation, `inclusive=(True, True)` means that the number of elements "masked" on each side should be rounded. In our case, `p*len(x) = 1.9`, so we would round to 2 on each side for a total of four.

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_round')
mask = (x < low) | (high < x)
np.count_nonzero(mask)
```

```{code-cell} ipython3
mxp.mean(mxp.asarray(x, mask=mask))
```

What happened? Apparently the documentation has been wrong since the function was added in 2008. If we trim *less* data, we get the answer we expect.

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_less')
mask = (x < low) | (high < x)
np.count_nonzero(mask)
```

```{code-cell} ipython3
mxp.mean(mxp.asarray(x, mask=mask))
```

What about `inclusive=(False, False)`? According to the documentation, the number of elements "masked" should be truncated. In our case, `p*len(x) = 1.9`, so we would trunctate to 1 on each side for a total of two. This should give us the same result as the `winsor_less` calculation just above.

```{code-cell} ipython3
stats.mstats.trimmed_mean(x, (p, p), inclusive=(False, False))
```

But instead, it gives us the same result as the `winsor_round` (or a `winsor_more`) calculation. Apparently, the documentation is backwards.

+++

As it turns out, `stats.mstats.trimmed_mean` is a convenience function equivalent to:

```python3
stats.mstatats.trim(a, limits=limits, inclusive=inclusive, relative=relative).mean(axis=axis)
```

So the behavior is best understood in terms of `stats.mstats.trim`.

+++

### `stats.mstats.trim`

+++

The entire implementation of `stats.mstats.trim` is:

```python3
if relative:
    return trimr(a, limits=limits, inclusive=inclusive, axis=axis)
else:
    return trima(a, limits=limits, inclusive=inclusive)
```

+++

So it is also convenience function that simply wraps the functions `stats.mstats.trimr` and `stats.mstats.trima`.

+++

### `stats.mstats.trima`

`stats.mstats.trima` is itself a convenience function that simply masks values outside of a given interval.

```{code-cell} ipython3
x = np.arange(20)
low, high = (1, 18)
stats.mstats.trima(x, (low, high), inclusive=(True, True))
```

The definition of `inclusive` is ambiguous at best:

> Tuple of (lower flag, upper flag), indicating whether values exactly equal to the lower (upper) limit are allowed.

What does *allowed* mean? Are values equal to the limits included in the data that is *trimmed* or *retained*? Apparently, since only one point has been masked on each side, it refers to elements being *retained*.

Generating and applying the mask manually is much more explicit. It also tends to be more concise.

```{code-cell} ipython3
mxp.asarray(x, mask=(x < low) | (x > high))  # mask values *strictly* less than 1; also mask values *strictly* greater than 18
```

To reproduce the results we'd get with `inclusive=(False, False)`:

```{code-cell} ipython3
stats.mstats.trima(x, (1, 18), inclusive=(False, False))
```

we simply adjust the inequalities.

```{code-cell} ipython3
mxp.asarray(x, mask=(x <= low) | (x >= high))  # mask values 1 or less; also mask values 18 or greater
```

### `stats.mstats.trimr`

`stats.mstats.trimr` is also a convenience function that trims a specified proportion of the *total* number of elements from each tail. `inclusive` has no effect when the product of the proportion and the total number of elements is integral.

```{code-cell} ipython3
p = 0.1
stats.mstats.trimr(x, (p, p), inclusive=(True, True))
```

```{code-cell} ipython3
stats.mstats.trimr(x, (p, p), inclusive=(False, False))
```

The documentation for `inclusive` states that it is a:

> Tuple of flags indicating whether the number of data being masked on the left (right) end should be truncated (True) or rounded (False) to integers.

+++

Suppose we have $n=19$ and the proportion to mask is $p = 0.1$. In that case, $pn = 1.9$. With `inclusive=(False, False)`, we round the number of elements to be masked to `2`.

```{code-cell} ipython3
n = 19.
x = np.arange(n)
stats.mstats.trimr(x, (p, p), inclusive=(False, False))
```

This is equivalent to:

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_round')
mxp.asarray(x, mask=(x < low) | (x > high))
```

With `inclusive=(True, True)`, the number of elements to be masked is truncated to `1`.

```{code-cell} ipython3
stats.mstats.trimr(x, (p, p), inclusive=(True, True))
```

```{code-cell} ipython3
low, high = stats.quantile(x, [p, 1-p], method='winsor_less')
mxp.asarray(x, mask=(x < low) | (x > high))
```

### `stats.mstats.trimboth`, `stats.trimboth`,  `stats.mstats.tail`, and `stats.trim1`.

For your convenience, `stats.mstats.trimboth` passes `proportiontocut` to `trimr` twice!

```{code-cell} ipython3
np.all(stats.mstats.trimboth(x, p) == stats.mstats.trimr(x, (p, p)))
```

`stats.trimboth` does the same job, but it removes (rather than masking) the trimmed data, and shuffles the elements.

```{code-cell} ipython3
a = stats.mstats.trimboth(x, p)
b = stats.trimboth(x, p)
np.all(np.ma.compressed(a) == np.sort(b))
```

`trimtail` performs the arduous task of passing `proportiontocut` and `None` as a tuple.

```{code-cell} ipython3
np.all(stats.mstats.trimtail(x, p, tail='left') == stats.mstats.trimr(x, (p, None)))
```

```{code-cell} ipython3
np.all(stats.mstats.trimtail(x, p, tail='right') == stats.mstats.trimr(x, (None, p)))
```

`stats.trim1` is to `stats.mstats.trimtail` as `stats.trimboth` is to `stats.mstats.trimboth`.

```{code-cell} ipython3
a = stats.mstats.trimtail(x, p, tail='left')
b = stats.trim1(x, p, tail='left')
np.all(np.ma.compressed(a) == np.sort(b))
```

Replicating these results using `quantile` and `MArray` is left as an exercise for the reader.

+++

### stats.mstats.winsorize

+++

In almost all cases, `stats.mstats.winsorize` is equivalent to:

```{code-cell} ipython3
def winsorize(x, limits, inclusive, nan_policy='propagate'):
    method = 'winsor_less' if all(inclusive) else 'winsor_round'
    low_high = stats.quantile(x, [limits[0], 1-limits[1]], method=method, axis=-1, nan_policy=nan_policy)
    return np.clip(x, low_high[..., :1], low_high[..., 1:])
```

```{code-cell} ipython3
stats.mstats.winsorize(x, (p, p), inclusive=(False, False)).data
```

```{code-cell} ipython3
winsorize(x, (p, p), inclusive=(False, False))
```

Minor generalizations would be needed to support `axis` other than `-1`, different left and right values of `inclusive`, and the `inplace` argument.

Results may be different when NaNs are present and `nan_policy='omit'`. This is because `stats.mstats.winsorize` is incorrect. If we have two NaNs to be omitted, we would have 17 elements left, so $np = 1.7$, and `inclusive=(False, False)` would round this number up to 2. This means that two of the values *on each side* would be replaced with the next more central element.

```{code-cell} ipython3
x[9:11] = np.nan
stats.mstats.winsorize(x, (p, p), inclusive=(False, False), nan_policy='omit').data
```

Clearly, no elements in the right tail have been replaced. Our simple implementation fixes the problem.

```{code-cell} ipython3
winsorize(x, (p, p), inclusive=(False, False), nan_policy='omit')
```
