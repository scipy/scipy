## About

This is a Python translation of [Zaikun Zhang](https://www.zhangzk.net)'s [modern-Fortran reference implementation](https://github.com/libprima/prima/tree/main/fortran)
for Powell's derivative-free optimization solvers, which is available at `fortran/` under the root directory.
It is supposed to be a faithful translation, [producing bit-for-bit identical results](https://github.com/scipy/scipy/pull/22350#issue-2795978526).
as the modern-Fortran reference implementation.
If you notice a difference, [raise an issue](https://github.com/libprima/prima/issues/new).

Due to [bug-fixes](https://github.com/libprima/prima#bug-fixes) and [improvements](https://github.com/libprima/prima#improvements),
the modern-Fortran reference implementation by [Zaikun Zhang](https://www.zhangzk.net)
behaves differently from the original Fortran 77 implementation by [M. J. D. Powell](https://www.zhangzk.net/powell.html),
even though the algorithms are essentially the same. Therefore, it is important to point out that you are using
PRIMA rather than the original solvers if you want your results to be reproducible.

Compared to Powell's Fortran 77 implementation, the modern-Fortran implementation and hence any faithful
translation like this one generally [produce better solutions with fewer function evaluations](https://github.com/libprima/prima#improvements),
making them preferable for [applications with expensive function evaluations](https://github.com/orgs/libprima/discussions/145).
However, if function evaluations are not the dominant cost in your application, the Fortran 77
solvers are likely to be faster, as they are more efficient in terms of memory usage and flops
thanks to the careful and ingenious (but unmaintained and unmaintainable) implementation by Powell.

As of April 2025, only the COBYLA solver is available in this Python translation
(many thanks to [Nickolai Belakovski](http://www.nickolai.me/)), and SciPy 1.16.0
integrates it to replace the original Fortran 77 implementation of [COBYLA underlying the
`scipy.optimize.minimize` function](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyla.html).
The other solvers will be translated from the Fortran reference implementation in the future.
If you are interested in doing so, contact [Zaikun Zhang](https://www.zhangzk.net).

## Development notes

To develop, `cd` into the `src` directory and run

```pip install --editable .```

This will install PRIMA locally in an editable fashion. From there you can run the examples/cobyla/cobyla_example.py (from any directory) and go from there.

### Style notes

- Most of the comments are copied from Fortran verbatim, except in cases where they need to modified due to specifics of the Python language. In these cases a note will be made of the difference between Fortran and Python
  - Rationale:
      - The main purpose of this is to keep the Python and Fortran codebases as similar as possible.
- For determining the dimensions of an array, we exclusively use `np.size` instead of `np.shape` or `some_array.shape` or `len`
  - Rationale:
    - Fortran uses `SIZE` so this helps us to be as consistent with the Fortran code as possible.

### A note on Fortran's `maxval`, `maximum`, and `maxval` and their Python equivalents

| Fortran   | Python       | Return value |
|-----------|--------------|--------------|
| `maxval`  | `max`        | scalar       |
| `maximum` | `np.max`     | scalar       |
| `max`     | `np.maximum` | vector       |

The difference between `maxval` and `maximum` is that `maximum` will return NaN if the input contains NaN. Python's `max`
and numpy's `np.max` have a similar distinction.

Fortran's `max` and numpy's `mp.maximum` accept two arguments, either of which can be a scalar or an array,
and returns an elementwise maximum of the two arguments. In the case of a scalar and an array argument it
returns an elementwise maximum of the scalar and each element of the array.

This note applies to `minval`, `minimum`, and `min` as well.


### A note on indices

Consider the following Fortran code

```
do i=0:5
  print *, *
end do
```

It can be easily and automatically translated to Python as

```
for i in range(0, 6):
  print(i)
```

Now consider the following similar loop

```
do i=1:5
  print *, some_array(i)
end do
```

This can be translated to Python as

```
for i in range(1, 6):
  print(some_array[i-1])
```

This leads to awkward Python code, since the more pythonic code would range from 0 to 5, and the indexing would be `some_array[i]`. In order to make the Python code more usable, we will attempt to write more "pythonic" code, even though that makes the translation a little bit more difficult.
