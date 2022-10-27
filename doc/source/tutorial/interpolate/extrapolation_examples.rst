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


.. _tutorial-extrapolation-cubicspline-extend:

CubicSpline extend the boundary conditions
==========================================

closes gh-14472

`CubicSpline` needs two extra boundary conditions, which are
controlled by the ``bc_type`` parameter. This parameter can either list
explicit values of derivatives at the edges, or use helpful aliases. For
instance, ``bc_type="clamped"`` sets the first derivatives to zero,
``bc_type="natural"`` sets the second derivatives to zero, and so on.

While the extrapolation is controlled by the boundary condition, the
relation is not very intuitive. For instance, one can expect that for
``bc_type="natural"``, the extrapolation is linear. This expectation is
too strong: each boundary condition sets the derivatives at a single
point, *at the boundary* only. Extrapolation is done from the first and
last polynomial pieces, which — for a natural spline — is a cubic with a
zero second derivative at a given point.

One other way of seeing why this expectation is too strong is to
consider a dataset with only three data points, where the spline has two
polynomial pieces. To extrapolate linearly, both of these pieces need to
be linear. But then, two linear pieces cannot match at a middle point
with a continuous 2nd derivative! (Unless of course, if all three data
points actually lie on a single straight line).

To illustrate the behavior we consider a synthetic dataset and compare
several boundary conditions:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import CubicSpline
    
    xs = [1, 2, 3, 4, 5, 6, 7, 8]
    ys = [4.5, 3.6, 1.6, 0.0, -3.3, -3.1, -1.8, -1.7]
    
    notaknot = CubicSpline(xs, ys, bc_type='not-a-knot')
    natural = CubicSpline(xs, ys, bc_type='natural')
    clamped = CubicSpline(xs, ys, bc_type='clamped')
    xnew = np.linspace(min(xs) - 4, max(xs) + 4, 101)
    
    splines = [notaknot, natural, clamped]
    titles = ['not-a-knot', 'natural', 'clamped']
    
    fig, axs = plt.subplots(3, 3, figsize=(12, 12))
    for i in [0, 1, 2]:
        for j, spline, title in zip(range(3), splines, titles):
            axs[i, j].plot(xs, spline(xs, nu=i),'o')
            axs[i, j].plot(xnew, spline(xnew, nu=i),'-')
            axs[i, j].set_title(f'{title}, deriv={i}')
            
    plt.tight_layout()
    plt.show()

It is clearly seen that e.g. natural spline does have the zero second derivative
at the boundaries, but extrapolation is non-linear. ``bc_type="clamped"``
shows a similar behavior: first derivaties are only equal zero exactly at the
boundary. In all cases, extrapolation is done by extending the first and last
polynomial pieces of the spline, whatever they happen to be.

One possible way to force the extrapolation is to extend the interpolation domain
to add a first and last polynomial pieces which have desired properties.

Here we use ``extend`` method of the `CubicSpline` superclass,
`PPoly`, to add two extra breakpoints and to make sure that the additional
polynomial pieces maintain the values of the derivatives. Then the extrapolation
proceeds using these two additional intervals.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import CubicSpline
    
    def add_boundary_knots(spline):
        """
        Add knots infinitesimally to the left and right with zero 2nd and 3rd derivative.
        
        Maintain the first derivative from whatever boundary condition was selected,
        Modifies the spline in place.
        """
        # determine the slope at the left edge
        leftx = spline.x[0]
        lefty = spline(leftx)
        leftslope = spline(leftx, nu=1)
        
        # add a new breakpoint just to the left and use the
        # known slope to construct the PPoly coefficients.
        leftxnext = np.nextafter(leftx, leftx - 1)
        leftynext = lefty + leftslope*(leftxnext - leftx)
        leftcoeffs = np.array([0, 0, leftslope, leftynext])
        spline.extend(leftcoeffs[..., None], np.r_[leftxnext])
        
        # repeat with additional knots to the right
        rightx = spline.x[-1]
        righty = spline(rightx)
        rightslope = spline(rightx,nu=1)
        rightxnext = np.nextafter(rightx, rightx + 1)
        rightynext = righty + rightslope * (rightxnext - rightx)
        rightcoeffs = np.array([0, 0, rightslope, rightynext])
        spline.extend(rightcoeffs[..., None], np.r_[rightxnext])
    
    xs = [1, 2, 3, 4, 5, 6, 7, 8]
    ys = [4.5, 3.6, 1.6, 0.0, -3.3, -3.1, -1.8, -1.7]
    
    notaknot = CubicSpline(xs,ys, bc_type='not-a-knot')
    # not-a-knot does not require additional boundary conditions for extrapolation
    
    natural = CubicSpline(xs,ys, bc_type='natural')
    # extend the natural natural spline with linear extrapolating knots
    add_boundary_knots(natural)
    
    clamped = CubicSpline(xs,ys, bc_type='clamped')
    # extend the clamped spline with constant extrapolating knots
    add_boundary_knots(clamped)
    
    xnew = np.linspace(min(xs) - 5, max(xs) + 5, 201)
    
    fig, axs = plt.subplots(3, 3,figsize=(12,12))
    
    splines = [notaknot, natural, clamped]
    titles = ['not-a-knot', 'natural', 'clamped']
    
    for i in [0, 1, 2]:
        for j, spline, title in zip(range(3), splines, titles):
            axs[i, j].plot(xs, spline(xs, nu=i),'o')
            axs[i, j].plot(xnew, spline(xnew, nu=i),'-')
            axs[i, j].set_title(f'{title}, deriv={i}')
    
    plt.tight_layout()
    plt.show()


