"""Lloyd-Max algorithm."""
from __future__ import annotations

import numpy as np
from typing import (
    Optional,
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy._lib._util import (
        DecimalNumber, IntNumber
    )

from scipy.spatial import distance, Voronoi


__all__ = ['lloyd_centroidal_voronoi_tessellation']


def _l1_norm(points: np.ndarray) -> float:
    l1 = distance.cdist(points, points, 'cityblock')
    return np.min(l1[l1.nonzero()])


def _lloyd_centroidal_voronoi_tessellation(
        points: np.ndarray,
        decay: float,
        qhull_options: str,
) -> np.ndarray:
    """Lloyd-Max algorithm iteration.

    Based on the implementation of Stéfan van der Walt:

    https://github.com/stefanv/lloyd

    which is:

        Copyright (c) 2021-04-21 Stéfan van der Walt
        https://github.com/stefanv/lloyd
        MIT License

    Parameters
    ----------
    points : array_like (n, d)
        The points to iterate on.
    decay : float
        Relaxation decay. A positive value would move the points toward
        their centroid, and negative value would move them away.
        1 would move the points to their centroid.
    qhull_options : str
        Additional options to pass to Qhull. See Qhull manual
        for details. (Default: "Qbb Qc Qz Qj Qx" for ndim > 4 and
        "Qbb Qc Qz Qj" otherwise.)

    Returns
    -------
    points : array_like (n, d)
        The points after an iteration of Lloyd's algorithm.

    """
    new_points = np.empty_like(points)

    voronoi = Voronoi(points, qhull_options=qhull_options)

    for ii, idx in enumerate(voronoi.point_region):
        # the region is a series of indices into self.voronoi.vertices
        # remove point at infinity, designated by index -1
        region = [i for i in voronoi.regions[idx] if i != -1]

        # get the vertices for this region
        verts = voronoi.vertices[region]

        # clipping would be wrong, we need to intersect
        # verts = np.clip(verts, 0, 1)

        # move points towards centroids:
        # Centroid in n-D is the mean for uniformly distributed nodes
        # of a geometry.
        centroid = np.mean(verts, axis=0)
        new_points[ii] = points[ii] + (centroid - points[ii]) * decay

    # only update points to centroid within the region
    is_valid = np.all(np.logical_and(new_points >= 0, new_points <= 1), axis=1)
    points[is_valid] = new_points[is_valid]

    return points


def lloyd_centroidal_voronoi_tessellation(
        points: npt.ArrayLike,
        *,
        tol: DecimalNumber = 1e-5,
        maxiter: IntNumber = 10,
        decay: Optional[npt.ArrayLike] = None,
        qhull_options: Optional[str] = None,
) -> np.ndarray:
    """Approximate Centroidal Voronoi Tessellation.

    Perturb points in :math:`[0, 1]^d` using Lloyd-Max algorithm.

    Parameters
    ----------
    points : array_like (n, d)
        The points to iterate on. With ``n`` the number of points and ``d``
        the dimension.
    tol : float, optional
        Tolerance for termination. If the min of the L1-norm over the points
        changes less than `tol`, it stops the algorithm. Default is 1e-5.
    maxiter : int, optional
        Maximum number of iterations. It will stop the algorithm even if
        `tol` is above the threshold.
        Too many iterations tend to cluster the points as a hypersphere.
        Default is 10.
    decay : {float, array_like (maxiter)}, optional
        Relaxation decay. A positive value would move the points toward
        their centroid, and negative value would move them away.
        1 would move the points to their centroid. Default is a varying
        exponential decay starting at 2 and ending at 1 after `maxiter`.
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual
        for details. (Default: "Qbb Qc Qz Qj Qx" for ndim > 4 and
        "Qbb Qc Qz Qj" otherwise.)

    Returns
    -------
    points : array_like (n, d)
        The points after being processed by Lloyd-Max algorithm.

    Notes
    -----
    Lloyd-Max algorithm is an iterative process with the purpose of improving
    the dispersion of points. For given points: (i) compute a Voronoi
    Tessellation; (ii) find the centroid of each Voronoi cell; (iii) move the
    points toward the centroid of their respective cell.

    The process converges to equally spaced points. It implies that measures
    like the discrepancy could suffer from too many iterations. On the other
    hand, L1 and L2 distances should improve. This is especially true with
    QMC methods which tend to favor the discrepancy over other criteria.

    .. warning::

       The Voronoi Tessellation step is expensive and quickly becomes
       intractable with dimensions as low as 10 even for a sample of points
       of size as low as 1000.

    References
    ----------
    .. [1] Lloyd. "Least Squares Quantization in PCM".
       IEEE Transactions on Information Theory, 1982.
    .. [2] Max J. "Quantizing for minimum distortion".
       IEEE Transactions on Information Theory, 1960.

    Examples
    --------
    >>> from scipy.spatial.distance import cdist
    >>> from scipy.spatial import lloyd_centroidal_voronoi_tessellation
    >>> rng = np.random.default_rng()
    >>> points = rng.random((128, 2))

    Compute the quality of the points using the L1 criterion.

    >>> def l1_norm(points):
    ...    l1 = cdist(points, points, 'cityblock')
    ...    return np.min(l1[l1.nonzero()])

    >>> l1_norm(points)
    0.00161...  # random

    Now process the points using Lloyd's algorithm and check the improvement
    on the L1. The value should increase.

    >>> points = lloyd_centroidal_voronoi_tessellation(points)
    >>> l1_norm(points)
    0.0278...  # random

    """
    points = np.asarray(points)

    if not points.ndim == 2:
        raise ValueError('Sample is not a 2D array')

    # Checking that points are within the bounds
    if not ((points.max() <= 1.) and (points.min() >= 0.)):
        raise ValueError('Sample is out of bounds')

    if qhull_options is None:
        qhull_options = 'Qbb Qc Qz QJ'

        if points.shape[1] >= 5:
            qhull_options += ' Qx'

    if decay is None:
        # Fit an exponential to be 2 at 0 and 1 at `maxiter`.
        # The decay is used for relaxation.
        # analytical solution for y=exp(-maxiter/x) - 0.1
        root = -maxiter / np.log(0.1)
        decay = [np.exp(-x / root)+0.9 for x in range(maxiter)]
    else:
        try:
            decay = np.broadcast_to(decay, maxiter)  # type: ignore[arg-type]
        except ValueError as exc:
            msg = ('decay is not a list of float'
                   ' of length maxiter')
            raise ValueError(msg) from exc

    for i in range(maxiter):
        l1_old = _l1_norm(points=points)
        points = _lloyd_centroidal_voronoi_tessellation(
                points=points, decay=decay[i],
                qhull_options=qhull_options,
        )

        l1_new = _l1_norm(points=points)

        if abs(l1_new - l1_old) < tol:
            break

    return points
