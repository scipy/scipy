---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"tags": ["jupyterlite_sphinx_strip"]}

```{eval-rst}
.. notebooklite:: missing_data.md
   :new_tab: True
```

(missing_data)=

+++

# Working with Missing Data

Suppose a longitudinal study is intended to measure some feature of a subject at a number of points in time. To keep track of which measured value corresponds with each time point, the data is stored in an array: index 0 corresponds with the first time, index 1 corresponds with the second, and so on.

```{code-cell} ipython3
import os
os.environ['SCIPY_ARRAY_API'] = '1'
import numpy as np
# Measurements at each of five times
data = np.asarray([1.2, 2.3, 3.4, 4.5, 5.6])
```

If a measurement is not made at a designated time, *something* must fill the corresponding element of the array. There are a few common approaches.

In the simplest approach, a "sentinel value", which does not appear elsewhere in the valid measurement data, is chosen to represent a missing measurement. Frequently, this is a value that lies outside the possible range of measurement, such as a negative value when the possible range is strictly positive.

```{code-cell} ipython3
# A sentinel value, -1., can be used to represent a
# missing measurement of a positive quantity
data = np.asarray([1.1, 1.2, 1.3, -1., 1.5])
```

Suppose we wish to take the harmonic mean of the valid measurements using {func}`scipy.stats.hmean<scipy.stats.hmean>`. One approach is to eliminate the sentinel values, producing a temporary array of a smaller size, and to pass this temporary array to {func}`scipy.stats.hmean<scipy.stats.hmean>`.

```{code-cell} ipython3
from scipy import stats
temp = data[data > 0]
stats.hmean(temp)
```

But suppose we have more than one subject, with different missing measurements for each subject.

```{code-cell} ipython3
data = np.asarray([[1.1, 1.2, 1.3, -1., 1.5],   # four valid measurements, subject 1
                   [2.9, -1., -1., 2.6, 2.5]])  # three valid measurements, subject 2
```

If all measurements were valid, we could simply pass the whole array into `hmean` and specify that the calculation is to be performed independently for each slice (row) along the last axis; i.e.:

```{code-cell} ipython3
stats.hmean(data, axis=-1)
```

Of course, these negative values are not automatically ignored, so they produce invalid results.
We also cannot simply remove the missing values before calling `hmean`, either:

```{code-cell} ipython3
temp = data[data > 0]
temp
```

This produces a one-dimensional array, so `hmean` would not be able to produce separate harmonic means for each subject.

```{code-cell} ipython3
stats.hmean(temp, axis=-1)
```

One solution is to loop over the rows:

```{code-cell} ipython3
res = []
for row in data:
    temp = row[row > 0]
    res.append(stats.hmean(temp))
res = np.asarray(res)
res
```

This is valid, but cumbersome and potentially slow for datasets with many subjects. To facilitate batched computation with missing data, SciPy provides two approaches.

The simplest requires the user to represent missing data using the floating point value NaN (Not a Number) as the sentinel:

```{code-cell} ipython3
NaN = np.nan
# data[data < 0] = NaN, or more explicitly:
data = np.asarray([[1.1, 1.2, 1.3, NaN, 1.5],   # four valid measurements, subject 1
                   [2.9, NaN, NaN, 2.6, 2.5]])  # three valid measurements, subject 2
```

then use `hmean` with the option `nan_policy='omit'`.

```{code-cell} ipython3
stats.hmean(data, axis=-1, nan_policy='omit')
```

Almost all reducing statistics in {mod}`scipy.stats<scipy.stats>` support `nan_policy='omit'`. Coverage is nearly complete because it is implemented in the generic way: looping over the slices, and eliminating the NaNs before performing the operation for each slice. As discussed, this can be slow when there are many slices, so it is offered merely for batch calculation convenience, not for speed.

The other option is to fill the space of missing values with arbitrary data, and to use a second, boolean array of the same shape - a "mask" - to keep track of which elements are missing.

```{code-cell} ipython3
data = np.asarray([[1.1, 1.2, 1.3, 1.4, 1.5],
                   [2.9, 2.8, 2.7, 2.6, 2.5]])
mask = np.asarray([[False, False, False,  True, False],
                   [False,  True,  True, False, False]])
```

NumPy offers `np.ma.masked_array`, for working with masked data, and functions in the {mod}`scipy.stats<scipy.stats>` were provided to work with these NumPy masked arrays.

In principle, the masked array approach is advantageous because it has the potential to avoid conflating *missing* NaN values with *invalid* NaN values, such as the result of `0 / 0`. It can also be faster in batch calculations with many slices, because batched masked array calculations can be implemented to ignore masked values without introducing Python `for` loops.

However, the {fun}`scipy.stats.mstats.hmean<scipy.stats.mstats.hmean>` function is now deprecated along with the {mod}`scipy.stats.mstats<scipy.stats.mstats>` namespace and all other uses of `np.ma.masked_array` in {mod}`scipy.stats<scipy.stats`.

```{code-cell} ipython3
x = np.ma.masked_array(data, mask=mask)
stats.mstats.hmean(x, axis=-1)
```

There are several reasons.

The first is that NumPy masked arrays do, in fact, conflate invalid and missing values. Consider the following example:

```{code-cell} ipython3
x = np.asarray([0, 1, 2, 3, 4])
np.sum(x / x)
```

Ordinary NumPy arrays warn that `0 / 0` produces NaN, and this invalid value propagates in the sum. It is impossible to get a valid numerical result when a NaN is involved in arithmetic.

Yet NumPy masked arrays seem to provide a number.

```{code-cell} ipython3
y = np.ma.masked_array(x)
np.sum(y / y)
```

This occurs because NaNs arising from invalid numerical calculations involving NumPy masked arrays are  masked without warning (and subsequently ignored).

```{code-cell} ipython3
y / y
```

This can lead to invalid calculations producing bogus but harmless-looking numerical results, which is clearly unsafe in scientific computing. Rather than flagging the problem and allowing the user to *fix* it, it hides the problem and produces seemingly correct but erroneous values. 

+++

The second reason for deprecating `mstats` is that `mstats` function interfaces and implementations were entirely separate from those of `scipy.stats`, and often neglected in terms of maintenance and enhancements. Consider, for instance, the stark difference in documentation thoroughness and feature completeness between [`scipy.stats.mannwhitneyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) and [`scipy.stats.mstats.mannwhitneyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.mannwhitneyu.html#scipy.stats.mstats.mannwhitneyu). Eliminating `mstats` in favor of adding missing data allows SciPy maintainers to provide users with one, complete implementation with a single interface. 

+++

The final reason for deprecating `scipy.stats.mstats` and support for NumPy masked arrays is the rise in support for Array API compatible arrays throughout SciPy. NumPy masked arrays themselves are mostly unmaintained and do not conform to the Array API Standard; so when adding major capabilities in GPU in JIT computing via CuPy, PyTorch, and JAX arrays via the array API, it is difficult to preserve support for the legacy `np.ma.masked_array` type.

Fortunately, as support for the Array API Standard closes this window, it opens the door toward a new way of supporting masked data. Specifically, [MArray](https://mdhaber.github.io/marray/tutorial.html) is a an Array API compatible array type that *wraps* the functionality of other array backends and endows them with support for masks.

```{code-cell} ipython3
from marray import numpy as xp
x = xp.asarray(data, mask=mask)
x
```

```{code-cell} ipython3
stats.hmean(x, axis=-1)
```

Consequently, existing users of `scipy.stats.mstats` and NumPy masked arrays are advised to switch to using corresponding `scipy.stats` functions with NumPy-backed MArrays, where possible, and to using the `scipy.stats` functions with `nan_policy='omit'` otherwise.
