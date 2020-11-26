A Design Specification for `nan_policy`
=======================================

Many functions in `scipy.stats` have a parameter called ``nan_policy``
that determines how the function handles data that contains ``nan``.  In
this section, we provide SciPy developer guidelines for how ``nan_policy``
is intended to be used, to ensure that as this parameter is added to new
functions, we maintain a consistent API.

The basic API
-------------

The parameter ``nan_policy`` accepts three possible strings: ``'omit'``,
``'raise'`` and ``'propagate'``.  The meanings are:

* ``nan_policy='omit'``:
  Ignore occurrences of ``nan`` in the input.  Do not generate a warning
  if the input contains ``nan``. For example, for the simple case of a
  function that accepts a single array (and ignoring the possible use of
  ``axis`` for the moment)::

      func([1.0, 3.0, np.nan, 5.0], nan_policy='omit')

  should behave the same as::

      func([1.0, 3.0, 5.0])

  More generally, ``func(a, nan_policy='omit')`` should behave the same as
  ``func(a[~np.isnan(a)])``.

  Unit tests for this property should be used to test functions that
  handle ``nan_policy``.

  For functions that accept two or more arguments but whose values are
  not related, the same idea applies to each input array.  So::

      func(a, b, nan_policy='omit')

  should behave the same as::

      func(a[~np.isnan(a)], b[~np.isnan(b)])

  For inputs with *related* or *paired* values, the recommended behavior
  is to omit all the values for which any of the related values are ``nan``.
  For a function with two related array inputs, this means::

      y = func(a, b, nan_policy='omit')

  should behave the same as::

      hasnan = np.isnan(a) | np.isnan(b)  # Union of the isnan masks.
      y = func(a[~hasnan], b[~hasnan])

  The docstring for such a function should clearly state this behavior.

* ``nan_policy='raise'``:
  Raise a ``ValueError``.
* ``nan_policy='propagate'``:
  Propagate the ``nan`` value to the output.  Typically, this means just
  execute the function without checking for ``nan``, but see

      https://github.com/scipy/scipy/issues/7818

  for an example where that might lead to unexpected output.


``nan_policy`` combined with an ``axis`` parameter
--------------------------------------------------
There is nothing surprising here--the principle mentioned above still
applies when the function has an ``axis`` parameter.  Suppose, for example,
``func`` reduces a 1-d array to a scalar, and handles n-d arrays as a
collection of 1-d arrays, with the ``axis`` parameter specifying the axis
along which the reduction is to be applied.  If, say::

    func([1, 3, 4])     -> 10.0
    func([2, -3, 8, 2]) ->  4.2
    func([7, 8])        ->  9.5
    func([])            -> -inf

then::

    func([[  1, nan,   3,   4],
          [  2,  -3,   8,   2],
          [nan,   7, nan,   8],
          [nan, nan, nan, nan]], nan_policy='omit', axis=-1)

must give the result::

    np.array([10.0, 4.2, 9.5, -inf])


Edge cases
----------
A function that implements the ``nan_policy`` parameter should gracefully
handle the case where *all* the values in the input array(s) are ``nan``.
The basic principle described above still applies::

    func([nan, nan, nan], nan_policy='omit')

should behave the same as::

    func([])

In practice, when adding ``nan_policy`` to an existing function, it is
not unusual to find that the function doesn't already handle this case
in a well-defined manner, and some thought and design may have to be
applied to ensure that it works.  The correct behavior (whether that be
to return ``nan``, return some other value, raise an exception, or something
else) will be determined on a case-by-case basis.


Why doesn't ``nan_policy`` also apply to ``inf``?
--------------------------------------------------
Although we learn in grade school that "infinity is not a number", the
floating point values ``nan`` and ``inf`` are qualitatively different.
The values ``inf`` and ``-inf`` act much more like regular floating
point values than ``nan``.

* One can compare ``inf`` to other floating point values and it behaves
  as expected, e.g. ``3 < inf`` is True.
* For the most part, arithmetic works "as expected" with ``inf``,
  e.g. ``inf + inf = inf``, ``-2*inf = -inf``, ``1/inf = 0``,
  etc.
* Many existing functions work "as expected" with ``inf``:
  ``np.log(inf) = inf``, ``np.exp(-inf) = 0``,
  ``np.array([1.0, -1.0, np.inf]).min() = -1.0``, etc.

So while ``nan`` almost always means "something went wrong" or "something
is missing", ``inf`` can in many cases be treated as a useful floating
point value.

It is also consistent with the NumPy ``nan`` functions to not ignore
``inf``::

    >>> np.nanmax([1, 2, 3, np.inf, np.nan])
    inf
    >>> np.nansum([1, 2, 3, np.inf, np.nan])
    inf
    >>> np.nanmean([8, -np.inf, 9, 1, np.nan])
    -inf


How *not* to implement ``nan_policy``
-------------------------------------
In the past (and possibly currently), some ``stats`` functions handled
``nan_policy`` by using a masked array to mask the ``nan`` values, and
then computing the result using the functions in the ``mstats`` subpackage.
The problem with this approach is that the masked array code might convert
``inf`` to a masked value, which we don't want to do (see above).  It also
means that, if care is not taken, the return value will be a masked array,
which will likely be a surprise to the user if they passed in regular arrays.
