"""
Airspeed Velocity benchmark utilities
"""
import sys
import os
import re
import time
import textwrap
import subprocess
import itertools
import random

from asv_runner.benchmarks.mark import SkipNotImplemented


class Benchmark:
    """
    Base class with sensible options
    """


class XPBenchmark(Benchmark):
    """
    Base class for benchmarks that are run on multiple Array API backends
    and devices. Supports multiple devices, jax.jit, and lazy/asynchronous
    evaluation.

    Basic usage
    -----------
    ::
        def myfunc(x):
            return x + 1

        class MyFunc(XPBenchmark):
            def setup(self, backend):
                super().setup(backend, myfunc)
                x = self.xp.arange(5)
                self.x = self.synchronize(x)
                if self.warmup:
                    self.func(self.x)

            def time_myfunc(self, backend):
                self.func(self.x)

    Adding parameters
    -----------------
    In the below example:
    - We add a `size` asv parameter
    - We add a `plus` function parameter which can't be traced by jax.jit

    ::
        def myfunc(x, plus=True):
            return x + 1 if plus else x - 1

        class MyFunc(XPBenchmark):
            param_names = (*XPBenchmark.param_names, "size")
            params = (*XPBenchmark.params, [5, 10])

            def setup(self, backend, size):
                super().setup(backend, myfunc, static_argnames=("plus",))
                x = self.xp.arange(size)
                self.x = self.synchronize(x)
                if self.warmup:
                    self.func(self.x, plus=True)
                    self.func(self.x, plus=False)

            def time_myfunc_plus(self, backend, size):
                self.func(self.x, plus=True)

            def time_myfunc_minus(self, backend, size):
                self.func(self.x, plus=False)
    """
    backends = ["numpy", "array_api_strict", "cupy", "torch:cpu", "torch:cuda",
                "dask.array", "jax.numpy:cpu", "jax.numpy:cuda"]

    # subclasses can override these
    param_names = ("backend",)
    params = (backends, )

    def setup(self, backend, func, *, static_argnums=None, static_argnames=None):
        """Skip benchmark if backend/device combination is not available.
        Configure namespace.
        Potentially wrap func with jax.jit and ensure timings are correct
        for lazy backends.

        Parameters
        ----------
        backend : str
            backend name from asv parameterization
        func : callable
            function to benchmark
        static_argnums : Sequence[int], optional
            Parameter for jax.jit. Note that, unlike in the unit tests,
            we can't use the automatic parameter and return value wrap/unwrap
            from `array_api_extra.testing.lazy_xp_function`, as it comes with a
            substantial performance overhead.
        static_argnames : Sequence[str], optional
            Parameter for jax.jit

        Sets attributes
        ---------------
        backend : str
            As the parameter (for convenience of helper functions)
        xp : namespace
            array namespace, potentially wrapped by array_api_compat
        func : callable
            function to benchmark, potentially wrapped
        warmup : bool
            Whether setup() should run a warmup iteration
        """
        self.backend = backend
        if ":" in backend:
            backend, device = backend.split(":")
        else:
            device = "cuda" if backend == "cupy" else "cpu"

        with safe_import() as array_api_imports:
            # Requires scipy >=1.16
            from scipy._lib._array_api import array_namespace, xp_capabilities_table
            from scipy.conftest import xp_available_backends, xp_known_backends

            if isinstance(xp_available_backends, dict):  # scipy == 1.16
                backends = xp_available_backends
            else:  # scipy >= 1.17
                backends = {p.id: p.values[0] for p in xp_available_backends}
        if array_api_imports.error:
            # On older scipy versions, disregard SCIPY_ARRAY_API
            import numpy as np
            def array_namespace(*args, **kwargs):
                return np
            xp_capabilities_table = {}
            backends = {"numpy": np}
            xp_known_backends = {"numpy"}

        # If new backends are added to conftest.py, you need to add them here too
        assert not xp_known_backends - set(n.split(":")[0] for n in self.backends)

        try:
            xp = backends[backend]
        except KeyError:
            raise SkipNotImplemented(
                f"{backend} not available or skipped by SCIPY_ARRAY_API")

        if func and func in xp_capabilities_table:
            capabilities = xp_capabilities_table[func]
            skips = {n for n, _ in capabilities["skip_backends"]}
            skips |= {n for n, _ in capabilities["xfail_backends"]}
            if (((capabilities["cpu_only"] and device != "cpu")
                 or (capabilities["np_only"] and backend != "numpy"))
                and backend not in capabilities["exceptions"]):
                skips.add(backend)
            if backend in skips:
                raise SkipNotImplemented(f"{backend} skipped by @xp_capabilities")
        else:
            capabilities = {"jax_jit": False}

        # Potentially wrap namespace with array_api_compat
        xp = array_namespace(xp.empty(0))

        self.xp = xp
        self.func = func
        self.warmup = False

        if backend == "torch":
            import torch

            torch.set_default_dtype(torch.float64)
            try:
                torch.empty(0, device=device)
            except (RuntimeError, AssertionError):
                raise SkipNotImplemented(f"{device=} not available")
            torch.set_default_device(device)

            if device == "cuda":
                def wrapper(*args, **kwargs):
                    res = func(*args, **kwargs)
                    torch.cuda.synchronize()
                    return res

                self.func = wrapper

        elif backend == "jax.numpy":
            import jax

            jax.config.update("jax_enable_x64", True)
            try:
                jax_device = jax.devices(device)[0]
            except RuntimeError:
                raise SkipNotImplemented(f"{device=} not available")
            jax.config.update("jax_default_device", jax_device)

            if capabilities["jax_jit"]:
                func = jax.jit(func, static_argnames=static_argnames,
                               static_argnums=static_argnums)
                self.warmup = True

            def wrapper(*args, **kwargs):
                res = func(*args, **kwargs)
                jax.block_until_ready(res)
                return res

            self.func = wrapper

        elif backend == "dask.array":
            import dask

            def wrapper(*args, **kwargs):
                res = func(*args, **kwargs)
                return dask.compute(res)[0]

            self.func = wrapper

        elif backend == "cupy":
            import cupy
            # The default stream is non-blocking.
            # As of CuPy 13.4.1, explicit non-blocking streams
            # are substantially slower.
            # cupy.cuda.Stream(non_blocking=True).use()

            def wrapper(*args, **kwargs):
                res = func(*args, **kwargs)
                cupy.cuda.get_current_stream().synchronize()
                return res

            self.func = wrapper

        else:
            assert backend in ("numpy", "array_api_strict")

    def synchronize(self, *arrays):
        """Wait until the given arrays have finished generating and return a
        synchronized instance of them.
        You need to call this on all arrays that your setup() function creates.
        """
        if self.backend == "dask.array":
            import dask

            arrays = dask.persist(*arrays)
        elif self.backend in ("jax.numpy:cpu", "jax.numpy:cuda"):
            import jax

            jax.block_until_ready(arrays)
        elif self.backend == "torch:cuda":
            import torch

            torch.cuda.synchronize()
        elif self.backend == "cupy":
            import cupy

            cupy.cuda.get_current_stream().synchronize()
        else:
            assert self.backend in ("numpy", "array_api_strict", "torch:cpu")

        return arrays[0] if len(arrays) == 1 else arrays


def is_xslow():
    try:
        return int(os.environ.get('SCIPY_XSLOW', '0'))
    except ValueError:
        return False


class LimitedParamBenchmark(Benchmark):
    """
    Limits parameter combinations to `max_number` choices, chosen
    pseudo-randomly with fixed seed.
    Raises NotImplementedError (skip) if not in active set.
    """
    num_param_combinations = 0

    def setup(self, *args, **kwargs):
        slow = is_xslow()

        if slow:
            # no need to skip
            return

        param_seed = kwargs.pop('param_seed', None)
        if param_seed is None:
            param_seed = 1

        params = kwargs.pop('params', None)
        if params is None:
            params = self.params

        num_param_combinations = kwargs.pop('num_param_combinations', None)
        if num_param_combinations is None:
            num_param_combinations = self.num_param_combinations

        all_choices = list(itertools.product(*params))

        rng = random.Random(param_seed)
        rng.shuffle(all_choices)
        active_choices = all_choices[:num_param_combinations]

        if args not in active_choices:
            raise NotImplementedError("skipped")


def get_max_rss_bytes(rusage):
    """
    Extract the max RSS value in bytes.
    """
    if not rusage:
        return None

    if sys.platform.startswith('linux'):
        # On Linux getrusage() returns ru_maxrss in kilobytes
        # https://man7.org/linux/man-pages/man2/getrusage.2.html
        return rusage.ru_maxrss * 1024
    elif sys.platform == "darwin":
        # on macOS ru_maxrss is in bytes
        return rusage.ru_maxrss
    else:
        # Unknown, just return whatever is here.
        return rusage.ru_maxrss


def run_monitored_wait4(code):
    """
    Run code in a new Python process, and monitor peak memory usage.

    Returns
    -------
    duration : float
        Duration in seconds (including Python startup time)
    peak_memusage : int
        Peak memory usage in bytes of the child Python process

    Notes
    -----
    Works on Unix platforms (Linux, macOS) that have `os.wait4()`.
    """
    code = textwrap.dedent(code)

    start = time.time()
    process = subprocess.Popen([sys.executable, '-c', code])
    pid, returncode, rusage = os.wait4(process.pid, 0)
    duration = time.time() - start
    max_rss_bytes = get_max_rss_bytes(rusage)

    if returncode != 0:
        raise AssertionError(f"Running failed:\n{code}")

    return duration, max_rss_bytes


def run_monitored_proc(code):
    """
    Run code in a new Python process, and monitor peak memory usage.

    Returns
    -------
    duration : float
        Duration in seconds (including Python startup time)
    peak_memusage : float
        Peak memory usage (rough estimate only) in bytes

    """
    if not sys.platform.startswith('linux'):
        raise RuntimeError("Peak memory monitoring only works on Linux")

    code = textwrap.dedent(code)
    process = subprocess.Popen([sys.executable, '-c', code])

    peak_memusage = -1

    start = time.time()
    while True:
        ret = process.poll()
        if ret is not None:
            break

        with open(f'/proc/{process.pid}/status') as f:
            procdata = f.read()

        m = re.search(r'VmRSS:\s*(\d+)\s*kB', procdata, re.S | re.I)
        if m is not None:
            memusage = float(m.group(1)) * 1e3
            peak_memusage = max(memusage, peak_memusage)

        time.sleep(0.01)

    process.wait()

    duration = time.time() - start

    if process.returncode != 0:
        raise AssertionError(f"Running failed:\n{code}")

    return duration, peak_memusage


def run_monitored(code):
    """
    Run code in a new Python process, and monitor peak memory usage.

    Returns
    -------
    duration : float
        Duration in seconds (including Python startup time)
    peak_memusage : float or int
        Peak memory usage (rough estimate only) in bytes

    """

    if hasattr(os, 'wait4'):
        return run_monitored_wait4(code)
    else:
        return run_monitored_proc(code)


def get_mem_info():
    """Get information about available memory"""
    import psutil
    vm = psutil.virtual_memory()
    return {
        "memtotal": vm.total,
        "memavailable": vm.available,
    }


def set_mem_rlimit(max_mem=None):
    """
    Set address space rlimit
    """
    import resource
    if max_mem is None:
        mem_info = get_mem_info()
        max_mem = int(mem_info['memtotal'] * 0.7)
    cur_limit = resource.getrlimit(resource.RLIMIT_AS)
    if cur_limit[0] > 0:
        max_mem = min(max_mem, cur_limit[0])

    try:
        resource.setrlimit(resource.RLIMIT_AS, (max_mem, cur_limit[1]))
    except ValueError:
        # on macOS may raise: current limit exceeds maximum limit
        pass


def with_attributes(**attrs):
    def decorator(func):
        for key, value in attrs.items():
            setattr(func, key, value)
        return func
    return decorator


class safe_import:

    def __enter__(self):
        self.error = False
        return self

    def __exit__(self, type_, value, traceback):
        if type_ is not None:
            self.error = True
            suppress = not (
                os.getenv('SCIPY_ALLOW_BENCH_IMPORT_ERRORS', '1').lower() in
                ('0', 'false') or not issubclass(type_, ImportError))
            return suppress
