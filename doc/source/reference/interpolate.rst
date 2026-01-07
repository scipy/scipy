.. automodule:: scipy.interpolate
   :no-members:
   :no-inherited-members:
   :no-special-members:


Concurrency Notes
-----------------

Most interpolation objects in ``scipy.interpolate`` are thread-safe for read-only
operations (evaluation at new points). For details on thread-safety, see
:ref:`scipy_thread_safety`.

Interpolator objects created by calling one of the interpolation functions should
not be mutated from multiple threads simultaneously. Read-only operations (e.g.,
calling the interpolator to evaluate at new points) are generally safe.
