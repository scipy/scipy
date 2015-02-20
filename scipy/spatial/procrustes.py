"""
This module provides functions to perform full Procrustes analysis.

This code was originally written by Justin Kucynski and ported over from
scikit-bio by Yoshiki Vazquez-Baeza.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.linalg import orthogonal_procrustes


__all__ = ['procrustes']


def procrustes(data1, data2):
    r"""Procrustes analysis, a similarity test for two data sets

    Each input matrix is a set of points or vectors (the rows of the matrix).
    The dimension of the space is the number of columns of each matrix. Given
    two identically sized matrices, procrustes standardizes both such that:

    - :math:`tr(AA^{T}) = 1`.

    - Both sets of points are centered around the origin.

    Procrustes ([1]_, [2]_) then applies the optimal transform to the second
    matrix (including scaling/dilation, rotations, and reflections) to minimize
    :math:`M^{2}=\sum(data1-data2)^{2}`, or the sum of the squares of the
    pointwise differences between the two input datasets.

    This function was not designed to handle datasets with different numbers of
    datapoints (rows).  If two data sets have different dimensionality
    (different number of columns), simply add columns of zeros the smaller of
    the two.

    Parameters
    ----------
    data1 : array_like
        matrix, n rows represent points in k (columns) space data1 is the
        reference data, after it is standardised, the data from data2 will
        be transformed to fit the pattern in data1 (must have >1 unique
        points).

    data2 : array_like
        n rows of data in k space to be fit to data1.  Must be the  same
        shape (numrows, numcols) as data1 (must have >1 unique points).

    Returns
    -------
    mtx1 : array_like
        a standardized version of data1
    mtx2 : array_like
        the orientation of data2 that best fits data1. Centered, but not
        necessarily :math:`tr(AA^{T}) = 1`.
    disparity : float
        :math:`M^{2}` as defined above.

    See Also
    --------
    orthogonal_procrustes

    Notes
    -----
    - The disparity should not depend on the order of the input matrices, but
      the output matrices will, as only the first output matrix is guaranteed
      to be scaled such that :math:`tr(AA^{T}) = 1`.

    - Duplicate data points are generally ok, duplicating a data point will
      increase its effect on the procrustes fit.

    - The disparity scales as the number of points per input matrix.

    References
    ----------
    .. [1] Krzanowski, W. J. (2000). "Principles of Multivariate analysis".
    .. [2] Gower, J. C. (1975). "Generalized procrustes analysis".

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.spatial import procrustes
    # the b matrix is a rotated, shifted, scaled and mirrored version of a
    >>> a = np.array([[1, 3], [1, 2], [1, 1], [2, 1]], 'd')
    >>> b = np.array([[4, -2], [4, -4], [4, -6], [2, -6]], 'd')
    >>> mtx1, mtx2, disparity = procrustes(a, b)
    >>> round(disparity)
    0.0

    """
    num_rows, num_cols = np.shape(data1)
    if (num_rows, num_cols) != np.shape(data2):
        raise ValueError("Input matrices must be of same shape")
    if num_rows == 0 or num_cols == 0:
        raise ValueError("Input matrices must be >0 rows and >0 cols")

    # standardize each matrix
    mtx1 = _center(data1)
    mtx2 = _center(data2)

    if (not np.any(mtx1)) or (not np.any(mtx2)):
        raise ValueError("input matrices must contain >1 unique points")

    mtx1 = _normalize(mtx1)
    mtx2 = _normalize(mtx2)

    # transform mtx2 to minimize disparity (sum((mtx1[i,j] - mtx2[i,j])^2))
    mtx2, disparity = orthogonal_procrustes(mtx1, mtx2)

    # When mtx1 and mtx2 can match exactly, orthogonal_procrustes returns a
    # matrix with:
    #   |-1 0|
    #   | 0 1|
    # Thus we only really need to return mtx1 and a disparity of zero.
    #
    # Note: we check for shape first, because if the shape is not equal
    # np.allclose will raise an error.
    if mtx2.shape == (2, 2) and np.allclose(np.array([[-1, 0], [0, 1]]), mtx2,
                                            atol=1e-8):
        mtx2 = mtx1
        disparity = 0

    return mtx1, mtx2, disparity


def _center(mtx):
    """Translate all data (rows of the matrix) to center on the origin

    Parameters
    ----------
    mtx : array_like
        Matrix to translate the data for.

    Returns
    -------
    result : array_like ('d') array
        Shifted version of the input data.  The new matrix is such that the
        center of mass of the row vectors is centered at the origin.

    """
    result = np.array(mtx, 'd')
    result -= np.mean(result, 0)
    # subtract each column's mean from each element in that column
    return result


def _normalize(mtx):
    """change scaling of data (in rows) such that trace(mtx*mtx') = 1

    Parameters
    ----------
    mtx : array_like
        Matrix to scale the data for.

    Notes
    -----
    mtx' denotes the transpose of mtx

    """
    mtx = np.asarray(mtx, dtype=float)
    return mtx / np.linalg.norm(mtx)
