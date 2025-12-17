.. _scipy_thread_safety:

Thread Safety in SciPy
======================

SciPy supports use in a multithreaded context via the `threading` module in
the standard library. Many SciPy operations release the GIL, as does NumPy (and
a lot of SciPy functionality is implemented as calls to NumPy functions) - so
unlike many situations in Python, it is possible to improve parallel
performance by exploiting multithreaded parallelism in Python.

The easiest performance gains happen when each worker thread owns its own array
or set of array objects, with no data directly shared between threads. Threads
that spend most of their time in low-level code will typically run in parallel.

It is possible to share NumPy arrays between threads, but extreme care must be
taken to avoid creating thread safety issues when mutating arrays that are
shared between multiple threads - please see the NumPy documentation on thread
safety for more details. SciPy functions will not mutate arrays that the user
passes in, unless a function explicitly documents that it will do so (which is
rare). Hence calling SciPy functions in a threaded fashion on the same NumPy array
is safe.

While most of SciPy consists of *functions*, more care has to be taken with
*classes* and *data structures*.

Classes that have state, such as some of the integration and interpolation
objects in `scipy.integrate` and `scipy.interpolate`, are typically robust
against being called in parallel. They either accept parallel calls or raise an
informative error. For example, `scipy.integrate.ode` may raise an
``IntegratorConcurrencyError`` for integration methods that do not support
parallel execution.

SciPy offers a couple of data structures, namely sparse arrays and matrices in
`scipy.sparse`, and k-D trees in `scipy.spatial`. These data structures are
*currently not thread-safe*. Please avoid in particular operations that mutate
a data structure, like using item or slice assignment on sparse arrays, while
the data is shared across multiple threads. That may result in data corruption,
crashes, or other unwanted behavior.

Note that operations that *do not* release the GIL will see no performance
gains from use of the `threading` module, and instead might be better served
with `multiprocessing`.


Free-threaded Python
--------------------

.. versionadded:: 1.15.0

Starting with SciPy 1.15.0 and CPython 3.13, SciPy has experimental support
for Python runtimes with the GIL disabled on all platforms. See
https://py-free-threading.github.io for more information about installing and
using free-threaded Python.

Because free-threaded Python does not have a global interpreter lock (GIL) to
serialize access to Python objects, there are more opportunities for threads to
mutate shared state and create thread safety issues. All SciPy functionality is
tested for usage from parallel threads, however we expect there to be issues
that are as yet undiscovered - if you run into a problem, please check the
`GitHub issues with the free-threading label <https://github.com/scipy/scipy/issues?q=is%3Aopen+is%3Aissue+label%3Afree-threading>`__
and open a new issue if one does not exist yet for the function that is
misbehaving.
