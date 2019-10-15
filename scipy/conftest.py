# Pytest customization
from __future__ import division, absolute_import, print_function

import gc
import os
import pytest
import sys
import warnings

from distutils.version import LooseVersion
from scipy._lib._fpumode import get_fpu_mode
from scipy._lib._testutils import FPUModeChangeWarning


def pytest_configure(config):
    config.addinivalue_line("markers",
        "slow: Tests that are very slow.")
    config.addinivalue_line("markers",
        "xslow: mark test as extremely slow (not run unless explicitly requested)")


def pytest_runtest_setup(item):
    if LooseVersion(pytest.__version__) >= LooseVersion("3.6.0"):
        mark = item.get_closest_marker("xslow")
    else:
        mark = item.get_marker("xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('SCIPY_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; set environment variable SCIPY_XSLOW=1 to run it")


@pytest.fixture(scope="function", autouse=True)
def check_fpu_mode(request):
    """
    Check FPU mode was not changed during the test.
    """
    old_mode = get_fpu_mode()
    yield
    new_mode = get_fpu_mode()

    if old_mode != new_mode:
        warnings.warn("FPU mode changed from {0:#x} to {1:#x} during "
                      "the test".format(old_mode, new_mode),
                      category=FPUModeChangeWarning, stacklevel=0)


@pytest.yield_fixture(scope="module")
def pool():
    """Fixture for tests that use multiprocessing.pool."""
    # See https://bugs.python.org/issue38501
    # We only need to decorate tests that do not use `with mapwrapper:`
    # because the problem occurs when close/terminate are not explicitly
    # called, and garbage collection is supposed to take care of it.
    if sys.version_info >= (3, 8):
        pytest.xfail('Python 3.8 hangs when cleaning up a pool')
    yield
    gc.collect()  # clean up
