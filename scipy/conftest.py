# Pytest customization
import json
import os
import warnings

import numpy as np
import numpy.array_api
import numpy.testing as npt
import pytest

from scipy._lib._fpumode import get_fpu_mode
from scipy._lib._testutils import FPUModeChangeWarning
from scipy._lib import _pep440
from scipy._lib._array_api import SCIPY_ARRAY_API, SCIPY_DEVICE


def pytest_configure(config):
    config.addinivalue_line("markers",
        "slow: Tests that are very slow.")
    config.addinivalue_line("markers",
        "xslow: mark test as extremely slow (not run unless explicitly requested)")
    config.addinivalue_line("markers",
        "xfail_on_32bit: mark test as failing on 32-bit platforms")
    try:
        import pytest_timeout  # noqa:F401
    except Exception:
        config.addinivalue_line(
            "markers", 'timeout: mark a test for a non-default timeout')


def _get_mark(item, name):
    if _pep440.parse(pytest.__version__) >= _pep440.Version("3.6.0"):
        mark = item.get_closest_marker(name)
    else:
        mark = item.get_marker(name)
    return mark


def pytest_runtest_setup(item):
    mark = _get_mark(item, "xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('SCIPY_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; set environment variable SCIPY_XSLOW=1 to run it")
    mark = _get_mark(item, 'xfail_on_32bit')
    if mark is not None and np.intp(0).itemsize < 8:
        pytest.xfail(f'Fails on our 32-bit test platform(s): {mark.args[0]}')

    # Older versions of threadpoolctl have an issue that may lead to this
    # warning being emitted, see gh-14441
    with npt.suppress_warnings() as sup:
        sup.filter(pytest.PytestUnraisableExceptionWarning)

        try:
            from threadpoolctl import threadpool_limits

            HAS_THREADPOOLCTL = True
        except Exception:  # observed in gh-14441: (ImportError, AttributeError)
            # Optional dependency only. All exceptions are caught, for robustness
            HAS_THREADPOOLCTL = False

        if HAS_THREADPOOLCTL:
            # Set the number of openmp threads based on the number of workers
            # xdist is using to prevent oversubscription. Simplified version of what
            # sklearn does (it can rely on threadpoolctl and its builtin OpenMP helper
            # functions)
            try:
                xdist_worker_count = int(os.environ['PYTEST_XDIST_WORKER_COUNT'])
            except KeyError:
                # raises when pytest-xdist is not installed
                return

            if not os.getenv('OMP_NUM_THREADS'):
                max_openmp_threads = os.cpu_count() // 2  # use nr of physical cores
                threads_per_worker = max(max_openmp_threads // xdist_worker_count, 1)
                try:
                    threadpool_limits(threads_per_worker, user_api='blas')
                except Exception:
                    # May raise AttributeError for older versions of OpenBLAS.
                    # Catch any error for robustness.
                    return


@pytest.fixture(scope="function", autouse=True)
def check_fpu_mode(request):
    """
    Check FPU mode was not changed during the test.
    """
    old_mode = get_fpu_mode()
    yield
    new_mode = get_fpu_mode()

    if old_mode != new_mode:
        warnings.warn("FPU mode changed from {:#x} to {:#x} during "
                      "the test".format(old_mode, new_mode),
                      category=FPUModeChangeWarning, stacklevel=0)


# Array API backend handling
xp_available_backends = {'numpy': np}

if SCIPY_ARRAY_API and isinstance(SCIPY_ARRAY_API, str):
    # fill the dict of backends with available libraries
    xp_available_backends.update({'numpy.array_api': numpy.array_api})

    try:
        import torch  # type: ignore[import]
        xp_available_backends.update({'pytorch': torch})
        # can use `mps` or `cpu`
        torch.set_default_device(SCIPY_DEVICE)
    except ImportError:
        pass

    try:
        import cupy  # type: ignore[import]
        xp_available_backends.update({'cupy': cupy})
    except ImportError:
        pass

    # by default, use all available backends
    if SCIPY_ARRAY_API.lower() not in ("1", "true"):
        SCIPY_ARRAY_API_ = json.loads(SCIPY_ARRAY_API)

        if 'all' in SCIPY_ARRAY_API_:
            pass  # same as True
        else:
            # only select a subset of backend by filtering out the dict
            try:
                xp_available_backends = {
                    backend: xp_available_backends[backend]
                    for backend in SCIPY_ARRAY_API_
                }
            except KeyError:
                msg = f"'--array-api-backend' must be in {xp_available_backends.keys()}"
                raise ValueError(msg)

if 'cupy' in xp_available_backends:
    SCIPY_DEVICE = 'cuda'

array_api_compatible = pytest.mark.parametrize("xp", xp_available_backends.values())

skip_if_array_api = pytest.mark.skipif(
    SCIPY_ARRAY_API,
    reason="do not run with Array API on",
)

skip_if_array_api_gpu = pytest.mark.skipif(
    SCIPY_ARRAY_API and SCIPY_DEVICE != 'cpu',
    reason="do not run with Array API on and not on CPU",
)


def skip_if_array_api_backend(backend):
    def wrapper(func):
        reason = (
            f"do not run with Array API backend: {backend}"
        )
        # method gets there as a function so we cannot use inspect.ismethod
        if '.' in func.__qualname__:
            def wrapped(self, *args, xp, **kwargs):
                if xp.__name__ == backend:
                    pytest.skip(reason=reason)
                return func(self, *args, xp, **kwargs)
        else:
            def wrapped(*args, xp, **kwargs):  # type: ignore[misc]
                if xp.__name__ == backend:
                    pytest.skip(reason=reason)
                return func(*args, xp, **kwargs)
        return wrapped
    return wrapper
