# Pytest customization
import json
import os
import warnings
import tempfile
from contextlib import contextmanager
from typing import Literal

import numpy as np
import numpy.testing as npt
import pytest
import hypothesis

from scipy._lib._fpumode import get_fpu_mode
from scipy._lib._array_api import (
    SCIPY_ARRAY_API, SCIPY_DEVICE, array_namespace, default_xp, xp_device
)
from scipy._lib._testutils import FPUModeChangeWarning
from scipy._lib.array_api_extra.testing import patch_lazy_xp_functions
from scipy._lib import _pep440

try:
    from scipy_doctest.conftest import dt_config
    HAVE_SCPDT = True
except ModuleNotFoundError:
    HAVE_SCPDT = False

try:
    import pytest_run_parallel  # noqa:F401
    PARALLEL_RUN_AVAILABLE = True
except Exception:
    PARALLEL_RUN_AVAILABLE = False


def pytest_configure(config):
    try:
        import pytest_timeout  # noqa:F401
    except Exception:
        config.addinivalue_line(
            "markers", 'timeout: mark a test for a non-default timeout')
    try:
        # This is a more reliable test of whether pytest_fail_slow is installed
        # When I uninstalled it, `import pytest_fail_slow` didn't fail!
        from pytest_fail_slow import parse_duration  # type: ignore[import-not-found] # noqa:F401,E501
    except Exception:
        config.addinivalue_line(
            "markers", 'fail_slow: mark a test for a non-default timeout failure')

    if not PARALLEL_RUN_AVAILABLE:
        config.addinivalue_line(
            'markers',
            'parallel_threads(n): run the given test function in parallel '
            'using `n` threads.')
        config.addinivalue_line(
            "markers",
            "thread_unsafe: mark the test function as single-threaded",
        )
        config.addinivalue_line(
            "markers",
            "iterations(n): run the given test function `n` times in each thread",
        )


def pytest_runtest_setup(item):
    mark = item.get_closest_marker("xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('SCIPY_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; "
                        "set environment variable SCIPY_XSLOW=1 to run it")
    mark = item.get_closest_marker("xfail_on_32bit")
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
        warnings.warn(f"FPU mode changed from {old_mode:#x} to {new_mode:#x} during "
                      "the test",
                      category=FPUModeChangeWarning, stacklevel=0)


if not PARALLEL_RUN_AVAILABLE:
    @pytest.fixture
    def num_parallel_threads():
        return 1


# Array API backend handling
xp_available_backends = {'numpy': np}

if SCIPY_ARRAY_API and isinstance(SCIPY_ARRAY_API, str):
    # fill the dict of backends with available libraries
    try:
        import array_api_strict
        xp_available_backends.update({'array_api_strict': array_api_strict})
        if _pep440.parse(array_api_strict.__version__) < _pep440.Version('2.0'):
            raise ImportError("array-api-strict must be >= version 2.0")
        array_api_strict.set_array_api_strict_flags(
            api_version='2023.12'
        )
    except ImportError:
        pass

    try:
        import torch  # type: ignore[import-not-found]
        xp_available_backends.update({'torch': torch})
        # can use `mps` or `cpu`
        torch.set_default_device(SCIPY_DEVICE)

        # default to float64 unless explicitly requested
        default = os.getenv('SCIPY_DEFAULT_DTYPE', default='float64')
        if default == 'float64':
            torch.set_default_dtype(torch.float64)
        elif default != "float32":
            raise ValueError(
                "SCIPY_DEFAULT_DTYPE env var, if set, can only be either 'float64' "
               f"or 'float32'. Got '{default}' instead."
            )
    except ImportError:
        pass

    try:
        import cupy  # type: ignore[import-not-found]
        xp_available_backends.update({'cupy': cupy})
    except ImportError:
        pass

    try:
        import jax.numpy  # type: ignore[import-not-found]
        xp_available_backends.update({'jax.numpy': jax.numpy})
        jax.config.update("jax_enable_x64", True)
        jax.config.update("jax_default_device", jax.devices(SCIPY_DEVICE)[0])
    except ImportError:
        pass

    try:
        import dask.array as da
        xp_available_backends.update({'dask.array': da})
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

    # this is annoying in CuPy 13.x
    warnings.filterwarnings(
        'ignore', 'cupyx.jit.rawkernel is experimental', category=FutureWarning
    )
    from cupyx.scipy import signal
    del signal


@pytest.fixture(params=[
    pytest.param(v, id=k, marks=pytest.mark.array_api_backends)
    for k, v in xp_available_backends.items()
])
def xp(request, monkeypatch):
    """Run the test that uses this fixture on each available array API library.

    You can select all and only the tests that use the `xp` fixture by
    passing `-m array_api_backends` to pytest.

    You can select where individual tests run through the `@skip_xp_backends`,
    `@xfail_xp_backends`, and `@skip_xp_invalid_arg` pytest markers.

    Please read: https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#adding-tests
    """
    # Read all @pytest.marks.skip_xp_backends markers that decorate to the test,
    # if any, and raise pytest.skip() if the current xp is in the list.
    skip_or_xfail_xp_backends(request, "skip")
    # Read all @pytest.marks.xfail_xp_backends markers that decorate the test,
    # if any, and raise pytest.xfail() if the current xp is in the list.
    skip_or_xfail_xp_backends(request, "xfail")

    xp = request.param
    # Potentially wrap namespace with array_api_compat
    xp = array_namespace(xp.empty(0))

    if SCIPY_ARRAY_API:
        # If request.param==jax.numpy, wrap tested functions in jax.jit
        patch_lazy_xp_functions(
            xp=request.param, request=request, monkeypatch=monkeypatch
        )

        # Throughout all calls to assert_almost_equal, assert_array_almost_equal, and
        # xp_assert_* functions, test that the array namespace is xp in both the
        # expected and actual arrays. This is to detect the case where both arrays are
        # erroneously just plain numpy while xp is something else.
        with default_xp(xp):
            yield xp
    else:
        yield xp


skip_xp_invalid_arg = pytest.mark.skipif(SCIPY_ARRAY_API,
    reason = ('Test involves masked arrays, object arrays, or other types '
              'that are not valid input when `SCIPY_ARRAY_API` is used.'))


def _backends_kwargs_from_request(request, skip_or_xfail):
    """A helper for {skip,xfail}_xp_backends"""
    # do not allow multiple backends
    args_ = request.keywords[f'{skip_or_xfail}_xp_backends'].args
    if len(args_) > 1:
        # np_only / cpu_only has args=(), otherwise it's ('numpy',)
        # and we do not allow ('numpy', 'cupy')
        raise ValueError(f"multiple backends: {args_}")

    markers = list(request.node.iter_markers(f'{skip_or_xfail}_xp_backends'))
    backends = []
    kwargs = {}
    for marker in markers:
        if marker.kwargs.get('np_only', False):
            kwargs['np_only'] = True
            kwargs['exceptions'] = marker.kwargs.get('exceptions', [])
            kwargs['reason'] = marker.kwargs.get('reason', None)
        elif marker.kwargs.get('cpu_only', False):
            if not kwargs.get('np_only', False):
                # if np_only is given, it is certainly cpu only
                kwargs['cpu_only'] = True
                kwargs['exceptions'] = marker.kwargs.get('exceptions', [])
                kwargs['reason'] = marker.kwargs.get('reason', None)

        # add backends, if any
        if len(marker.args) == 1:
            backend = marker.args[0]
            backends.append(backend)
            kwargs.update(**{backend: marker.kwargs})
            for kwarg in ("cpu_only", "np_only", "exceptions"):
                if marker.kwargs.get(kwarg):
                    raise ValueError(f"{kwarg} is mutually exclusive with {backend}")
        elif len(marker.args) > 1:
            raise ValueError(
                f"Please specify only one backend per marker: {marker.args}"
            )

    return backends, kwargs


def skip_or_xfail_xp_backends(request: pytest.FixtureRequest,
                              skip_or_xfail: Literal['skip', 'xfail']) -> None:
    """
    Helper of the `xp` fixture.
    Skip or xfail based on the ``skip_xp_backends`` or ``xfail_xp_backends`` markers.

    See the "Support for the array API standard" docs page for usage examples.

    Usage
    -----
    ::
        skip_xp_backends = pytest.mark.skip_xp_backends
        xfail_xp_backends = pytest.mark.xfail_xp_backends
        ...

        @skip_xp_backends(backend, *, reason=None)
        @skip_xp_backends(*, cpu_only=True, exceptions=(), reason=None)
        @skip_xp_backends(*, np_only=True, exceptions=(), reason=None)

        @xfail_xp_backends(backend, *, reason=None)
        @xfail_xp_backends(*, cpu_only=True, exceptions=(), reason=None)
        @xfail_xp_backends(*, np_only=True, exceptions=(), reason=None)

    Parameters
    ----------
    backend : str, optional
        Backend to skip/xfail, e.g. ``"torch"``.
        Mutually exclusive with ``cpu_only`` and ``np_only``.
    cpu_only : bool, optional
        When ``True``, the test is skipped/xfailed on non-CPU devices,
        minus exceptions. Mutually exclusive with ``backend``.
    np_only : bool, optional
        When ``True``, the test is skipped/xfailed for all backends other
        than the default NumPy backend and the exceptions.
        Mutually exclusive with ``backend``. Implies ``cpu_only``.
    reason : str, optional
        A reason for the skip/xfail. If omitted, a default reason is used.
    exceptions : list[str], optional
        A list of exceptions for use with ``cpu_only`` or ``np_only``.
        This should be provided when delegation is implemented for some,
        but not all, non-CPU/non-NumPy backends.
    """
    if f"{skip_or_xfail}_xp_backends" not in request.keywords:
        return

    backends, kwargs = _backends_kwargs_from_request(
        request, skip_or_xfail=skip_or_xfail
    )
    xp = request.param
    skip_or_xfail = getattr(pytest, skip_or_xfail)
    np_only = kwargs.get("np_only", False)
    cpu_only = kwargs.get("cpu_only", False)
    exceptions = kwargs.get("exceptions", [])

    if reasons := kwargs.get("reasons"):
        raise ValueError(f"provide a single `reason=` kwarg; got {reasons=} instead")

    # input validation
    if np_only and cpu_only:
        # np_only is a stricter subset of cpu_only
        cpu_only = False

    # Test explicit backends first so that their reason can override
    # those from cpu_only
    for backend in backends:
        if xp.__name__ == backend:
            reason = kwargs[backend].get('reason')
            if not reason:
                reason = f"do not run with array API backend: {backend}"
            skip_or_xfail(reason=reason)

    if np_only:
        reason = kwargs.get("reason")
        if not reason:
            reason = "do not run with non-NumPy backends"

        if xp.__name__ != 'numpy' and xp.__name__ not in exceptions:
            skip_or_xfail(reason=reason)
        return

    if cpu_only:
        reason = kwargs.get("reason")
        if not reason:
            reason = ("no array-agnostic implementation or delegation available "
                      "for this backend and device")

        exceptions = [] if exceptions is None else exceptions
        if SCIPY_ARRAY_API and SCIPY_DEVICE != 'cpu':
            if xp.__name__ == 'cupy' and 'cupy' not in exceptions:
                skip_or_xfail(reason=reason)
            elif xp.__name__ == 'torch' and 'torch' not in exceptions:
                if 'cpu' not in xp.empty(0).device.type:
                    skip_or_xfail(reason=reason)
            elif xp.__name__ == 'jax.numpy' and 'jax.numpy' not in exceptions:
                for d in xp.empty(0).devices():
                    if 'cpu' not in d.device_kind:
                        skip_or_xfail(reason=reason)
            elif xp.__name__ == 'dask.array' and 'dask.array' not in exceptions:
                # dask has no device. 'cpu' is a hack introduced by array-api-compat.
                # Force to revisit this when in the future
                # dask adds proper device support
                assert xp_device(xp.empty(0)) == 'cpu'

# Following the approach of NumPy's conftest.py...
# Use a known and persistent tmpdir for hypothesis' caches, which
# can be automatically cleared by the OS or user.
hypothesis.configuration.set_hypothesis_home_dir(
    os.path.join(tempfile.gettempdir(), ".hypothesis")
)

# We register two custom profiles for SciPy - for details see
# https://hypothesis.readthedocs.io/en/latest/settings.html
# The first is designed for our own CI runs; the latter also
# forces determinism and is designed for use via scipy.test()
hypothesis.settings.register_profile(
    name="nondeterministic", deadline=None, print_blob=True,
)
hypothesis.settings.register_profile(
    name="deterministic",
    deadline=None, print_blob=True, database=None, derandomize=True,
    suppress_health_check=list(hypothesis.HealthCheck),
)

# Profile is currently set by environment variable `SCIPY_HYPOTHESIS_PROFILE`
# In the future, it would be good to work the choice into dev.py.
SCIPY_HYPOTHESIS_PROFILE = os.environ.get("SCIPY_HYPOTHESIS_PROFILE",
                                          "deterministic")
hypothesis.settings.load_profile(SCIPY_HYPOTHESIS_PROFILE)


############################################################################
# doctesting stuff

if HAVE_SCPDT:

    # FIXME: populate the dict once
    @contextmanager
    def warnings_errors_and_rng(test=None):
        """Temporarily turn (almost) all warnings to errors.

        Filter out known warnings which we allow.
        """
        known_warnings = dict()

        # these functions are known to emit "divide by zero" RuntimeWarnings
        divide_by_zero = [
            'scipy.linalg.norm', 'scipy.ndimage.center_of_mass',
        ]
        for name in divide_by_zero:
            known_warnings[name] = dict(category=RuntimeWarning,
                                        message='divide by zero')

        # Deprecated stuff in scipy.signal and elsewhere
        deprecated = [
            'scipy.signal.cwt', 'scipy.signal.morlet', 'scipy.signal.morlet2',
            'scipy.signal.ricker',
            'scipy.integrate.simpson',
            'scipy.interpolate.interp2d',
            'scipy.linalg.kron',
        ]
        for name in deprecated:
            known_warnings[name] = dict(category=DeprecationWarning)

        from scipy import integrate
        # the functions are known to emit IntegrationWarnings
        integration_w = ['scipy.special.ellip_normal',
                         'scipy.special.ellip_harm_2',
        ]
        for name in integration_w:
            known_warnings[name] = dict(category=integrate.IntegrationWarning,
                                        message='The occurrence of roundoff')

        # scipy.stats deliberately emits UserWarnings sometimes
        user_w = ['scipy.stats.anderson_ksamp', 'scipy.stats.kurtosistest',
                  'scipy.stats.normaltest', 'scipy.sparse.linalg.norm']
        for name in user_w:
            known_warnings[name] = dict(category=UserWarning)

        # additional one-off warnings to filter
        dct = {
            'scipy.sparse.linalg.norm':
                dict(category=UserWarning, message="Exited at iteration"),
            # tutorials
            'linalg.rst':
                dict(message='the matrix subclass is not',
                     category=PendingDeprecationWarning),
            'stats.rst':
                dict(message='The maximum number of subdivisions',
                     category=integrate.IntegrationWarning),
        }
        known_warnings.update(dct)

        # these legitimately emit warnings in examples
        legit = set('scipy.signal.normalize')

        # Now, the meat of the matter: filter warnings,
        # also control the random seed for each doctest.

        # XXX: this matches the refguide-check behavior, but is a tad strange:
        # makes sure that the seed the old-fashioned np.random* methods is
        # *NOT* reproducible but the new-style `default_rng()` *IS* repoducible.
        # Should these two be either both repro or both not repro?

        from scipy._lib._util import _fixed_default_rng
        import numpy as np
        with _fixed_default_rng():
            np.random.seed(None)
            with warnings.catch_warnings():
                if test and test.name in known_warnings:
                    warnings.filterwarnings('ignore',
                                            **known_warnings[test.name])
                    yield
                elif test and test.name in legit:
                    yield
                else:
                    warnings.simplefilter('error', Warning)
                    yield

    dt_config.user_context_mgr = warnings_errors_and_rng
    dt_config.skiplist = set([
        'scipy.linalg.LinAlgError',     # comes from numpy
        'scipy.fftpack.fftshift',       # fftpack stuff is also from numpy
        'scipy.fftpack.ifftshift',
        'scipy.fftpack.fftfreq',
        'scipy.special.sinc',           # sinc is from numpy
        'scipy.optimize.show_options',  # does not have much to doctest
        'scipy.signal.normalize',       # manipulates warnings (XXX temp skip)
        'scipy.sparse.linalg.norm',     # XXX temp skip
        # these below test things which inherit from np.ndarray
        # cross-ref https://github.com/numpy/numpy/issues/28019
        'scipy.io.matlab.MatlabObject.strides',
        'scipy.io.matlab.MatlabObject.dtype',
        'scipy.io.matlab.MatlabOpaque.dtype',
        'scipy.io.matlab.MatlabOpaque.strides',
        'scipy.io.matlab.MatlabFunction.strides',
        'scipy.io.matlab.MatlabFunction.dtype'
    ])

    # these are affected by NumPy 2.0 scalar repr: rely on string comparison
    if np.__version__ < "2":
        dt_config.skiplist.update(set([
            'scipy.io.hb_read',
            'scipy.io.hb_write',
            'scipy.sparse.csgraph.connected_components',
            'scipy.sparse.csgraph.depth_first_order',
            'scipy.sparse.csgraph.shortest_path',
            'scipy.sparse.csgraph.floyd_warshall',
            'scipy.sparse.csgraph.dijkstra',
            'scipy.sparse.csgraph.bellman_ford',
            'scipy.sparse.csgraph.johnson',
            'scipy.sparse.csgraph.yen',
            'scipy.sparse.csgraph.breadth_first_order',
            'scipy.sparse.csgraph.reverse_cuthill_mckee',
            'scipy.sparse.csgraph.structural_rank',
            'scipy.sparse.csgraph.construct_dist_matrix',
            'scipy.sparse.csgraph.reconstruct_path',
            'scipy.ndimage.value_indices',
            'scipy.stats.mstats.describe',
    ]))

    # help pytest collection a bit: these names are either private
    # (distributions), or just do not need doctesting.
    dt_config.pytest_extra_ignore = [
        "scipy.stats.distributions",
        "scipy.optimize.cython_optimize",
        "scipy.test",
        "scipy.show_config",
        # equivalent to "pytest --ignore=path/to/file"
        "scipy/special/_precompute",
        "scipy/interpolate/_interpnd_info.py",
        "scipy/_lib/array_api_compat",
        "scipy/_lib/highs",
        "scipy/_lib/unuran",
        "scipy/_lib/_gcutils.py",
        "scipy/_lib/doccer.py",
        "scipy/_lib/_uarray",
    ]

    dt_config.pytest_extra_xfail = {
        # name: reason
        "ND_regular_grid.rst": "ReST parser limitation",
        "extrapolation_examples.rst": "ReST parser limitation",
        "sampling_pinv.rst": "__cinit__ unexpected argument",
        "sampling_srou.rst": "nan in scalar_power",
        "probability_distributions.rst": "integration warning",
    }

    # tutorials
    dt_config.pseudocode = set(['integrate.nquad(func,'])
    dt_config.local_resources = {
        'io.rst': [
            "octave_a.mat",
            "octave_cells.mat",
            "octave_struct.mat"
        ]
    }

    dt_config.strict_check = True

    # ignore Matplotlib's `ax.text`:
    dt_config.stopwords.add('.text(')
############################################################################
