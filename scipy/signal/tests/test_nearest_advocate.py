import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from scipy.signal import nearest_advocate

SEED = 0


def test_nearest_advocate_base():
    N = 1_000
    DEF_DIST = 0.25

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.25, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + np.pi + np.random.normal(
        loc=0, scale=0.1, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=-60, td_max=60, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=True, dist_padding=DEF_DIST)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, np.pi, decimal=1)
    assert_almost_equal(min_mean_dist, 0.07694374, decimal=2)  


def test_nearest_advocate_edge():
    N = 1_000
    DEF_DIST = 0.25

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.5, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + np.pi + np.random.normal(
        loc=0, scale=0.4, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=-60, td_max=60, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=True, dist_padding=DEF_DIST)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, np.pi, decimal=1)
    assert_almost_equal(min_mean_dist, 0.1646458, decimal=2)


def test_nearest_advocate_base_defmax():
    N = 1_000
    DEF_DIST = 0.0

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.25, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + np.pi + np.random.normal(
        loc=0, scale=0.1, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=-60, td_max=60, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=True, dist_padding=DEF_DIST)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, np.pi, decimal=1)
    assert_almost_equal(min_mean_dist, 0.07690712, decimal=2)


def test_nearest_advocate_base_fewoverlap():
    N = 1_000
    DEF_DIST = 0.0
    TIME_SHIFT = 900

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.25, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + TIME_SHIFT + np.random.normal(
        loc=0, scale=0.1, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=850, td_max=950, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=True, dist_padding=DEF_DIST)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, TIME_SHIFT, decimal=1)
    assert_almost_equal(min_mean_dist, 0.07674612, decimal=2)


def test_nearest_advocate_base_nopadding():
    N = 1_000
    DEF_DIST = 0.25

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.25, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + np.pi + np.random.normal(
        loc=0, scale=0.1, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=-60, td_max=60, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=False)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, np.pi, decimal=1)
    assert_almost_equal(min_mean_dist, 0.07690712, decimal=2)


def test_nearest_advocate_base_noverlap():
    N = 100
    DEF_DIST = 0.25
    TIME_SHIFT = 200

    np.random.seed(SEED)
    arr_ref = np.sort(np.cumsum(np.random.normal(
        loc=1, scale=0.25, size=N))).astype(np.float32)
    arr_sig = np.sort(arr_ref + TIME_SHIFT + np.random.normal(
        loc=0, scale=0.1, size=N)).astype(np.float32)

    np_nearest = nearest_advocate(
        arr_ref=arr_ref, arr_sig=arr_sig,
        td_min=-60, td_max=60, sps=20, sparse_factor=1,
        dist_max=DEF_DIST, regulate_paddings=False)
    time_shift, min_mean_dist = np_nearest[np.argmin(np_nearest[:, 1])]

    assert_almost_equal(time_shift, -60, decimal=1)  # each value is the same
    assert_equal(min_mean_dist, DEF_DIST)
