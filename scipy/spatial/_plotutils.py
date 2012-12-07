import numpy as np

__all__ = ['delaunay_plot_2d', 'convex_hull_plot_2d', 'voronoi_plot_2d']

def delaunay_plot_2d(tri):
    """
    Plot the given Delaunay triangulation in 2-D

    Parameters
    ----------
    tri : scipy.spatial.Delaunay
        Triangulation to plot

    See Also
    --------
    matplotlib.pyplot.triplot

    See Also
    --------
    Delaunay

    Notes
    -----
    Requires Matplotlib.

    """
    import matplotlib.pyplot as plt
    plt.plot(tri.points[:,0], tri.points[:,1], 'o')
    plt.triplot(tri.points[:,0], tri.points[:,1], tri.simplices)

def convex_hull_plot_2d(hull):
    """
    Plot the given convex hull diagram in 2-D

    Parameters
    ----------
    hull : scipy.spatial.ConvexHull
        Convex hull to plot

    See Also
    --------
    ConvexHull

    Notes
    -----
    Requires Matplotlib.

    """
    import matplotlib.pyplot as plt

    plt.plot(hull.points[:,0], hull.points[:,1], 'o')
    for simplex in hull.simplices:
        plt.plot(hull.points[simplex,0], hull.points[simplex,1], 'k-')
    plt.show()

def voronoi_plot_2d(vor):
    """
    Plot the given Voronoi diagram in 2-D

    Parameters
    ----------
    vor : scipy.spatial.Voronoi
        Diagram to plot

    See Also
    --------
    Voronoi

    Notes
    -----
    Requires Matplotlib.

    """
    import matplotlib.pyplot as plt

    plt.plot(vor.points[:,0], vor.points[:,1], '.')
    plt.plot(vor.vertices[:,0], vor.vertices[:,1], 'o')

    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            plt.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], 'k-')

    ptp_bound = vor.points.ptp(axis=0)

    center = vor.points.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0] # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]]) # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * ptp_bound.max()

            plt.plot([vor.vertices[i,0], far_point[0]], 
                     [vor.vertices[i,1], far_point[1]], 'k--')
    plt.xlim(vor.points[:,0].min() - 0.1*ptp_bound[0],
             vor.points[:,0].max() + 0.1*ptp_bound[0])
    plt.ylim(vor.points[:,1].min() - 0.1*ptp_bound[1],
             vor.points[:,1].max() + 0.1*ptp_bound[1])
