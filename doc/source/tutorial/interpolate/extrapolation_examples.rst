.. currentmodule:: scipy.interpolate


.. _tutorial-extrapolation:

=============================
Extrapolation tips and tricks
=============================

Handling of extrapolation—evaluation of the interpolators on query
points outside of the domain of interpolated data— is not fully
consistent among different routines in `scipy.interpolate`. In some
applications you can use the ``extrapolate`` or ``fill_value`` keywords;
this may or may not be adequate to a particular problem.

Special attention needs to be paid to extrapolation of non-linear
interpolants. Very often the extrapolated results quickly stop making
sense with increasing the distance to the data domain. This is of course
to be expected: an interpolant only knows the data within the data
domain.

When the extrapolated results are not adequate, users need to implement
the desired extrapolation mode themselves.

In this tutorial, we consider several worked examples of both
approaches. These examples may or may not be applicable to your
particular problem; they are not necessarily best practices; and they
are deliberately pared down to bare essentials needed to demonstrate the
main ideas, in a hope that they serve as an inspiration for your
handling of your particular problem.


.. _tutorial-extrapolation-left-right:


`interp1d` : replicate `numpy.interp` left and right fill values
====================================================================

TL;DR: Use ``fill_value=(left, right)``

`numpy.interp` uses constant extrapolation, and defaults to extending
the first and last values of the ``y`` array in the interpolation
interval: the output of ``np.interp(xnew, x, y)`` is ``y[0]`` for
``xnew < x[0]`` and ``y[-1]`` for ``xnew > x[-1]``.

By default, `interp1d` refuses to extrapolate, and raises a
``ValueError`` when evaluated on a data point outside of the
interpolation range. This can be switched off by the
``bounds_error=False`` argument: then ``interp1d`` sets the out-of-range
values with the ``fill_value``, which is ``nan`` by default.

To mimic the behavior of `numpy.interp` with `interp1d`, you can use
the fact that it supports a 2-tuple as the ``fill_value``. The tuple
elements are then used to fill for ``xnew < min(x)`` and ``x > max(x)``,
respectively. For multidimensional ``y``, these elements must have the
same shape as ``y`` or be broadcastable to it.

To illustrate:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    
    x = np.linspace(0, 1.5*np.pi, 11)
    y = np.column_stack((np.cos(x), np.sin(x)))   # y.shape is (11, 2)
    
    func = interp1d(x, y,
                    axis=0,  # interpolate along columns
                    bounds_error=False,
                    kind='linear',
                    fill_value=(y[0], y[-1]))
    xnew = np.linspace(-np.pi, 2.5*np.pi, 51)
    ynew = func(xnew)
    
    fix, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    ax1.plot(xnew, ynew[:, 0])
    ax1.plot(x, y[:, 0], 'o')
    
    ax2.plot(xnew, ynew[:, 1])
    ax2.plot(x, y[:, 1], 'o')
    plt.tight_layout()

