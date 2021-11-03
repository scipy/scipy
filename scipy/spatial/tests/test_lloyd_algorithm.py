import pytest
import numpy as np
from numpy.testing import assert_allclose

from scipy.spatial.distance import cdist
from scipy.spatial import lloyd_centroidal_voronoi_tessellation
from scipy.spatial._lloyd_algorithm import _l1_norm


def test_lloyd():
    # mindist
    def l2_norm(points):
        l2 = cdist(points, points)
        return np.min(l2[l2.nonzero()])

    # quite sensible seed as it can go up before going further down
    rng = np.random.RandomState(1809831)
    points = rng.uniform(0, 1, size=(128, 2))
    base_l1 = _l1_norm(points)
    base_l2 = l2_norm(points)

    for _ in range(4):
        points_lloyd = lloyd_centroidal_voronoi_tessellation(
                points, maxiter=1,
        )
        curr_l1 = _l1_norm(points_lloyd)
        curr_l2 = l2_norm(points_lloyd)

        # higher is better for the distance measures
        assert base_l1 < curr_l1
        assert base_l2 < curr_l2

        base_l1 = curr_l1
        base_l2 = curr_l2

        points = points_lloyd


def test_lloyd_non_mutating():
    """
    Verify that the input points are not mutated in place and that they do not
    share memory with the output.
    """
    orig_points = np.array([[0.1, 0.1],
                            [0.1, 0.2],
                            [0.2, 0.1],
                            [0.2, 0.2]])
    points_copy = orig_points.copy()
    new_points = lloyd_centroidal_voronoi_tessellation(points=orig_points)
    assert_allclose(orig_points, points_copy)
    assert not np.may_share_memory(orig_points, new_points)


def test_lloyd_decay():
    rng = np.random.default_rng()
    points = rng.random((20, 2))

    maxiter = 5
    decay_const = lloyd_centroidal_voronoi_tessellation(
            points, decay=[1], maxiter=maxiter
    )
    decay_list = lloyd_centroidal_voronoi_tessellation(
            points, decay=[1, 1, 1, 1, 1], maxiter=maxiter
    )
    assert_allclose(decay_const, decay_list)


def test_lloyd_errors():
    with pytest.raises(ValueError, match=r"Sample is not a 2D array"):
        points = [0, 1, 0.5]
        lloyd_centroidal_voronoi_tessellation(points)

    with pytest.raises(ValueError, match=r"Sample is out of bounds"):
        points = [[-1.1, 0], [0.1, 0.4], [1, 2]]
        lloyd_centroidal_voronoi_tessellation(points)

    with pytest.raises(ValueError, match=r"decay is not a list"):
        points = [[0, 0], [1, 1]]
        decay = [1, 2, 3]
        maxiter = 2
        lloyd_centroidal_voronoi_tessellation(
                points, decay=decay, maxiter=maxiter
        )
