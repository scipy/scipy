import numpy as np
import pytest

from scipy.stats import mannwhitneyu


@pytest.mark.parametrize("offset", [0, 2**53])
@pytest.mark.parametrize("method", ["exact", "auto"])
def test_large_int64_order_preserved(offset, method):
    x = np.array([offset + i for i in range(4)], dtype=np.int64)
    y = np.array([offset + i for i in range(4, 8)], dtype=np.int64)

    assert len(set(map(int, np.concatenate([x, y])))) == 8

    res = mannwhitneyu(x, y, alternative="two-sided", method=method)
    assert res.statistic == 0
    assert res.pvalue == 1 / 35
    assert isinstance(res.statistic, np.floating)
    assert isinstance(res.pvalue, np.floating)


@pytest.mark.parametrize("alternative", ["two-sided", "greater", "less"])
def test_large_int64_alternatives_keep_order(alternative):
    offset = 2**53
    x = np.array([offset + i for i in range(4)], dtype=np.int64)
    y = np.array([offset + i for i in range(4, 8)], dtype=np.int64)

    res = mannwhitneyu(x, y, alternative=alternative, method="exact")

    assert res.statistic == 0
    assert isinstance(res.statistic, np.floating)
    assert isinstance(res.pvalue, np.floating)
    expected = {"two-sided": 1 / 35, "greater": 1.0, "less": 1 / 70}
    assert res.pvalue == expected[alternative]


def test_large_int64_axis_keeps_order():
    offset = 2**53
    x = np.array([[offset + i for i in range(4)] for _ in range(3)], dtype=np.int64)
    y = np.array(
        [[offset + i for i in range(4, 8)] for _ in range(3)], dtype=np.int64
    )

    res = mannwhitneyu(x, y, axis=1, method="exact")

    assert res.statistic.shape == (3,)
    assert np.all(res.statistic == 0)
    assert np.all(res.pvalue == 1 / 35)


def test_large_int64_asymptotic_keeps_order():
    offset = 2**53
    x = np.array([offset + i for i in range(4)], dtype=np.int64)
    y = np.array([offset + i for i in range(4, 8)], dtype=np.int64)

    res = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")

    assert res.statistic == 0
    assert res.pvalue < 0.05


def test_large_int64_axis0_keeps_order():
    offset = 2**53
    x = np.array([[offset + i for i in range(3)] for _ in range(4)], dtype=np.int64)
    y = np.array(
        [[offset + i for i in range(3, 6)] for _ in range(4)], dtype=np.int64
    )

    res = mannwhitneyu(x, y, axis=0, method="exact")

    assert res.statistic.shape == (3,)
    assert np.all(res.statistic == 0)
    assert np.all(res.pvalue == 1 / 35)


def test_large_int64_masked_input_keeps_order():
    offset = 2**53
    x = np.ma.array([offset + i for i in range(4)], mask=[False, False, False, True])
    y = np.ma.array([offset + i for i in range(4, 8)], mask=[True, False, False, False])

    res = mannwhitneyu(x, y, alternative="two-sided", method="exact")

    assert res.statistic == 0
    assert 0 <= res.pvalue <= 1


def test_large_int64_float_input_unchanged():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0, 7.0])

    res = mannwhitneyu(x, y, alternative="two-sided", method="exact")

    assert res.statistic == 0
    assert res.pvalue == 1 / 35


def test_large_uint64_order_preserved():
    offset = 2**53
    x = np.array([offset + i for i in range(4)], dtype=np.uint64)
    y = np.array([offset + i for i in range(4, 8)], dtype=np.uint64)

    res = mannwhitneyu(x, y, alternative="two-sided", method="exact")

    assert res.statistic == 0
    assert res.pvalue == 1 / 35


def test_large_int64_without_continuity_keeps_order():
    offset = 2**53
    x = np.array([offset + i for i in range(4)], dtype=np.int64)
    y = np.array([offset + i for i in range(4, 8)], dtype=np.int64)

    res = mannwhitneyu(
        x, y, alternative="two-sided", method="exact", use_continuity=False
    )

    assert res.statistic == 0
    assert res.pvalue == 1 / 35
