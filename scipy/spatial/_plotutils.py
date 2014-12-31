from __future__ import division, print_function, absolute_import

import numpy as np
from scipy._lib.decorator import decorator as _decorator

__all__ = ['delaunay_plot_2d', 'convex_hull_plot_2d', 'voronoi_plot_2d']


@_decorator
def _held_figure(func, obj, ax=None, **kw):
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure()
        ax = fig.gca()

    was_held = ax.ishold()
    try:
        ax.hold(True)
        return func(obj, ax=ax, **kw)
    finally:
        ax.hold(was_held)


def _adjust_bounds(ax, points):
    ptp_bound = points.ptp(axis=0)
    ax.set_xlim(points[:,0].min() - 0.1*ptp_bound[0],
                points[:,0].max() + 0.1*ptp_bound[0])
    ax.set_ylim(points[:,1].min() - 0.1*ptp_bound[1],
                points[:,1].max() + 0.1*ptp_bound[1])


@_held_figure
def delaunay_plot_2d(tri, ax=None):
    """
    Plot the given Delaunay triangulation in 2-D

    Parameters
    ----------
    tri : scipy.spatial.Delaunay instance
        Triangulation to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on

    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot

    See Also
    --------
    Delaunay
    matplotlib.pyplot.triplot

    Notes
    -----
    Requires Matplotlib.

    """
    if tri.points.shape[1] != 2:
        raise ValueError("Delaunay triangulation is not 2-D")

    ax.plot(tri.points[:,0], tri.points[:,1], 'o')
    ax.triplot(tri.points[:,0], tri.points[:,1], tri.simplices.copy())

    _adjust_bounds(ax, tri.points)

    return ax.figure


@_held_figure
def convex_hull_plot_2d(hull, ax=None):
    """
    Plot the given convex hull diagram in 2-D

    Parameters
    ----------
    hull : scipy.spatial.ConvexHull instance
        Convex hull to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on

    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot

    See Also
    --------
    ConvexHull

    Notes
    -----
    Requires Matplotlib.

    """
    if hull.points.shape[1] != 2:
        raise ValueError("Convex hull is not 2-D")

    ax.plot(hull.points[:,0], hull.points[:,1], 'o')
    for simplex in hull.simplices:
        ax.plot(hull.points[simplex,0], hull.points[simplex,1], 'k-')

    _adjust_bounds(ax, hull.points)

    return ax.figure


@_held_figure
def voronoi_plot_2d(vor, ax=None):
    """
    Plot the given Voronoi diagram in 2-D

    Parameters
    ----------
    vor : scipy.spatial.Voronoi instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on

    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot

    See Also
    --------
    Voronoi

    Notes
    -----
    Requires Matplotlib.

    """
    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    ax.plot(vor.points[:,0], vor.points[:,1], '.')
    ax.plot(vor.vertices[:,0], vor.vertices[:,1], 'o')

    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            ax.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], 'k-')

    ptp_bound = vor.points.ptp(axis=0)

    center = vor.points.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            ax.plot([vor.vertices[i,0], far_point[0]],
                    [vor.vertices[i,1], far_point[1]], 'k--')

    _adjust_bounds(ax, vor.points)

    return ax.figure
