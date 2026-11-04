import sys
from scipy._lib._testutils import _parse_size, _get_mem_available, assert_mem_gain_limit
import numpy as np
import pytest


def test__parse_size():
    expected = {
        '12': 12e6,
        '12 b': 12,
        '12k': 12e3,
        '  12  M  ': 12e6,
        '  12  G  ': 12e9,
        ' 12Tb ': 12e12,
        '12  Mib ': 12 * 1024.0**2,
        '12Tib': 12 * 1024.0**4,
    }

    for inp, outp in sorted(expected.items()):
        if outp is None:
            with pytest.raises(ValueError):
                _parse_size(inp)
        else:
            assert _parse_size(inp) == outp


def test__mem_available():
    # May return None on non-Linux platforms
    available = _get_mem_available()
    if sys.platform.startswith('linux'):
        assert available >= 0
    else:
        assert available is None or available >= 0


@pytest.mark.thread_unsafe(reason="tracemalloc is not threadsafe")
class TestPeakMemoryGain:
    def test_passes(self):
        with assert_mem_gain_limit(400_096):  # 500 x 100 x 8 + 96
            np.zeros((500, 100))

    def test_raises(self):
        with pytest.raises(KeyError):
            with assert_mem_gain_limit(400_096):
                np.zeros((500, 100))
                raise KeyError()

    def test_fails(self):
        with pytest.raises(AssertionError, match=".* is above limit "):
            with assert_mem_gain_limit(400_096):
                np.zeros((500, 101))

    def test_fails_raises(self):
        with pytest.raises(KeyError):
            with pytest.raises(AssertionError, match=".* is above limit "):
                with assert_mem_gain_limit(400_096):
                    np.zeros((500, 101))
                    raise KeyError()


def test_pytest_report_header(monkeypatch):
    from scipy._lib._pytest_header import TEST_ENV_VARS, environment_report_header
    for name in TEST_ENV_VARS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("SCIPY_DEFAULT_DTYPE", "float32")

    header = environment_report_header()
    assert header == (
        "SciPy test environment: SCIPY_ARRAY_API=<unset>, SCIPY_DEVICE=<unset>, "
        "SCIPY_DEFAULT_DTYPE=float32, SCIPY_XSLOW=<unset>, "
        "SCIPY_HYPOTHESIS_PROFILE=<unset>, SCIPY_AVAILABLE_MEM=<unset>"
    )
