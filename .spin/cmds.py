import os
import sys
import importlib
import importlib.util
import json
import traceback
import warnings
import math
import subprocess
from copy import deepcopy

import click
from spin import util
from spin.cmds import meson

from pathlib import Path

PROJECT_MODULE = "scipy"


# # Check that the meson git submodule is present
# curdir = pathlib.Path(__file__).parent
# meson_import_dir = curdir.parent / 'vendored-meson' / 'meson' / 'mesonbuild'
# if not meson_import_dir.exists():
#     raise RuntimeError(
#         'The `vendored-meson/meson` git submodule does not exist! ' +
#         'Run `git submodule update --init` to fix this problem.'
#     )

def configure_scipy_openblas(self, blas_variant='32'):
    """Create scipy-openblas.pc and scipy/_distributor_init_local.py

    Requires a pre-installed scipy-openblas32 wheel from PyPI.
    """
    basedir = os.getcwd()
    pkg_config_fname = os.path.join(basedir, "scipy-openblas.pc")

    if os.path.exists(pkg_config_fname):
        return None

    module_name = f"scipy_openblas{blas_variant}"
    try:
        openblas = importlib.import_module(module_name)
    except ModuleNotFoundError:
        raise RuntimeError(f"Importing '{module_name}' failed. "
                            "Make sure it is installed and reachable "
                            "by the current Python executable. You can "
                            f"install it via 'pip install {module_name}'.")

    local = os.path.join(basedir, "scipy", "_distributor_init_local.py")
    with open(local, "w", encoding="utf8") as fid:
        fid.write(f"import {module_name}\n")

    with open(pkg_config_fname, "w", encoding="utf8") as fid:
        fid.write(openblas.get_pkg_config())

physical_cores_cache = None

def _cpu_count_cgroup(os_cpu_count):
    # Cgroup CPU bandwidth limit available in Linux since 2.6 kernel
    cpu_max_fname = "/sys/fs/cgroup/cpu.max"
    cfs_quota_fname = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
    cfs_period_fname = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
    if os.path.exists(cpu_max_fname):
        # cgroup v2
        # https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html
        with open(cpu_max_fname) as fh:
            cpu_quota_us, cpu_period_us = fh.read().strip().split()
    elif os.path.exists(cfs_quota_fname) and os.path.exists(cfs_period_fname):
        # cgroup v1
        # https://www.kernel.org/doc/html/latest/scheduler/sched-bwc.html#management
        with open(cfs_quota_fname) as fh:
            cpu_quota_us = fh.read().strip()
        with open(cfs_period_fname) as fh:
            cpu_period_us = fh.read().strip()
    else:
        # No Cgroup CPU bandwidth limit (e.g. non-Linux platform)
        cpu_quota_us = "max"
        cpu_period_us = 100_000  # unused, for consistency with default values

    if cpu_quota_us == "max":
        # No active Cgroup quota on a Cgroup-capable platform
        return os_cpu_count
    else:
        cpu_quota_us = int(cpu_quota_us)
        cpu_period_us = int(cpu_period_us)
        if cpu_quota_us > 0 and cpu_period_us > 0:
            return math.ceil(cpu_quota_us / cpu_period_us)
        else:  # pragma: no cover
            # Setting a negative cpu_quota_us value is a valid way to disable
            # cgroup CPU bandwidth limits
            return os_cpu_count


def _cpu_count_affinity(os_cpu_count):
    # Number of available CPUs given affinity settings
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except NotImplementedError:
            pass

    # On PyPy and possibly other platforms, os.sched_getaffinity does not exist
    # or raises NotImplementedError, let's try with the psutil if installed.
    try:
        import psutil

        p = psutil.Process()
        if hasattr(p, "cpu_affinity"):
            return len(p.cpu_affinity())

    except ImportError:  # pragma: no cover
        if (
            sys.platform == "linux"
            and os.environ.get("LOKY_MAX_CPU_COUNT") is None
        ):
            # PyPy does not implement os.sched_getaffinity on Linux which
            # can cause severe oversubscription problems. Better warn the
            # user in this particularly pathological case which can wreck
            # havoc, typically on CI workers.
            warnings.warn(
                "Failed to inspect CPU affinity constraints on this system. "
                "Please install psutil or explicitly set LOKY_MAX_CPU_COUNT.",
                stacklevel=4
            )

    # This can happen for platforms that do not implement any kind of CPU
    # infinity such as macOS-based platforms.
    return os_cpu_count


def _cpu_count_user(os_cpu_count):
    """Number of user defined available CPUs"""
    cpu_count_affinity = _cpu_count_affinity(os_cpu_count)

    cpu_count_cgroup = _cpu_count_cgroup(os_cpu_count)

    # User defined soft-limit passed as a loky specific environment variable.
    cpu_count_loky = int(os.environ.get("LOKY_MAX_CPU_COUNT", os_cpu_count))

    return min(cpu_count_affinity, cpu_count_cgroup, cpu_count_loky)

def _count_physical_cores():
    """Return a tuple (number of physical cores, exception)

    If the number of physical cores is found, exception is set to None.
    If it has not been found, return ("not found", exception).

    The number of physical cores is cached to avoid repeating subprocess calls.
    """
    exception = None

    # First check if the value is cached
    global physical_cores_cache
    if physical_cores_cache is not None:
        return physical_cores_cache, exception

    # Not cached yet, find it
    try:
        if sys.platform == "linux":
            cpu_info = subprocess.run(
                "lscpu --parse=core".split(), capture_output=True, text=True
            )
            cpu_info = cpu_info.stdout.splitlines()
            cpu_info = {line for line in cpu_info if not line.startswith("#")}
            cpu_count_physical = len(cpu_info)
        elif sys.platform == "win32":
            cpu_info = subprocess.run(
                "wmic CPU Get NumberOfCores /Format:csv".split(),
                capture_output=True,
                text=True,
            )
            cpu_info = cpu_info.stdout.splitlines()
            cpu_info = [
                l.split(",")[1]
                for l in cpu_info
                if (l and l != "Node,NumberOfCores")
            ]
            cpu_count_physical = sum(map(int, cpu_info))
        elif sys.platform == "darwin":
            cpu_info = subprocess.run(
                "sysctl -n hw.physicalcpu".split(),
                capture_output=True,
                text=True,
            )
            cpu_info = cpu_info.stdout
            cpu_count_physical = int(cpu_info)
        else:
            raise NotImplementedError(f"unsupported platform: {sys.platform}")

        # if cpu_count_physical < 1, we did not find a valid value
        if cpu_count_physical < 1:
            raise ValueError(f"found {cpu_count_physical} physical cores < 1")

    except Exception as e:
        exception = e
        cpu_count_physical = "not found"

    # Put the result in cache
    physical_cores_cache = cpu_count_physical

    return cpu_count_physical, exception

def cpu_count(only_physical_cores=False):
    """Return the number of CPUs the current process can use.

    The returned number of CPUs accounts for:
     * the number of CPUs in the system, as given by
       ``multiprocessing.cpu_count``;
     * the CPU affinity settings of the current process
       (available on some Unix systems);
     * Cgroup CPU bandwidth limit (available on Linux only, typically
       set by docker and similar container orchestration systems);
     * the value of the LOKY_MAX_CPU_COUNT environment variable if defined.
    and is given as the minimum of these constraints.

    If ``only_physical_cores`` is True, return the number of physical cores
    instead of the number of logical cores (hyperthreading / SMT). Note that
    this option is not enforced if the number of usable cores is controlled in
    any other way such as: process affinity, Cgroup restricted CPU bandwidth
    or the LOKY_MAX_CPU_COUNT environment variable. If the number of physical
    cores is not found, return the number of logical cores.

    Note that on Windows, the returned number of CPUs cannot exceed 61 (or 60 for
    Python < 3.10), see:
    https://bugs.python.org/issue26903.

    It is also always larger or equal to 1.
    """
    # Note: os.cpu_count() is allowed to return None in its docstring
    os_cpu_count = os.cpu_count() or 1
    if sys.platform == "win32":
        # On Windows, attempting to use more than 61 CPUs would result in a
        # OS-level error. See https://bugs.python.org/issue26903. According to
        # https://learn.microsoft.com/en-us/windows/win32/procthread/processor-groups
        # it might be possible to go beyond with a lot of extra work but this
        # does not look easy.
        os_cpu_count = min(os_cpu_count, _MAX_WINDOWS_WORKERS)

    cpu_count_user = _cpu_count_user(os_cpu_count)
    aggregate_cpu_count = max(min(os_cpu_count, cpu_count_user), 1)

    if not only_physical_cores:
        return aggregate_cpu_count

    if cpu_count_user < os_cpu_count:
        # Respect user setting
        return max(cpu_count_user, 1)

    cpu_count_physical, exception = _count_physical_cores()
    if cpu_count_physical != "not found":
        return cpu_count_physical

    # Fallback to default behavior
    if exception is not None:
        # warns only the first time
        warnings.warn(
            "Could not find the number of physical cores for the "
            f"following reason:\n{exception}\n"
            "Returning the number of logical cores instead. You can "
            "silence this warning by setting LOKY_MAX_CPU_COUNT to "
            "the number of cores you want to use.",
            stacklevel=2
        )
        traceback.print_tb(exception.__traceback__)

    return aggregate_cpu_count

@click.command()
@click.option(
    '--werror', default=False, is_flag=True,
    help="Treat warnings as errors")
@click.option(
    '--gcov', default=False, is_flag=True,
    help="enable C code coverage via gcov (requires GCC)."
            "gcov output goes to build/**/*.gc*")
@click.option(
    '--asan', default=False, is_flag=True,
    help=("Build and run with AddressSanitizer support. "
            "Note: the build system doesn't check whether "
            "the project is already compiled with ASan. "
            "If not, you need to do a clean build (delete "
            "build and build-install directories)."))
@click.option(
    '--debug', '-d', default=False, is_flag=True, help="Debug build")
@click.option(
    '--release', '-r', default=False, is_flag=True, help="Release build")
@click.option(
    '--parallel', '-j', default=None, metavar='N_JOBS',
    help=("Number of parallel jobs for building. "
            "This defaults to the number of available physical CPU cores"))
@click.option(
    '--setup-args', '-C', default=[], multiple=True,
    help=("Pass along one or more arguments to `meson setup` "
            "Repeat the `-C` in case of multiple arguments."))
@click.option(
    '--show-build-log', default=False, is_flag=True,
    help="Show build output rather than using a log file")
@click.option(
    '--with-scipy-openblas', default=False, is_flag=True,
    help=("If set, use the `scipy-openblas32` wheel installed into the "
            "current environment as the BLAS/LAPACK to build against."))
@click.option(
    '--with-accelerate', default=False, is_flag=True,
    help=("If set, use `Accelerate` as the BLAS/LAPACK to build against."
            " Takes precedence over -with-scipy-openblas (macOS only)")
)
@click.option(
    '--tags', default="runtime,python-runtime,tests,devel",
    show_default=True, help="Install tags to be used by meson."
)
@click.argument("meson_args", nargs=-1)
@click.pass_context
def build(ctx, meson_args, with_scipy_openblas, jobs=None, clean=False, verbose=False, quiet=False, *args, **kwargs):
    """🔧 Build package with Meson/ninja and install

    MESON_ARGS are passed through e.g.:

    spin build -- -Dpkg_config_path=/lib64/pkgconfig

    The package is installed to build-install

    By default builds for release, to be able to use a debugger set CFLAGS
    appropriately. For example, for linux use

    CFLAGS="-O0 -g" spin build
    """
    MESON_ARGS = "meson_args"

    meson_args_ = deepcopy(ctx.params[MESON_ARGS])
    ctx.params[MESON_ARGS] = {}
    ctx.params[MESON_ARGS]["compile"] = tuple()
    ctx.params[MESON_ARGS]["install"] = tuple()
    ctx.params[MESON_ARGS]["setup"] = meson_args_

    if sys.platform == "cygwin":
        # Cygwin only has netlib lapack, but can link against
        # OpenBLAS rather than netlib blas at runtime.  There is
        # no libopenblas-devel to enable linking against
        # openblas-specific functions or OpenBLAS Lapack
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + ("-Dlapack=lapack", "-Dblas=blas")

    if ctx.params['werror']:
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + ("--werror", )

    if ctx.params['debug'] or ctx.params['release']:
        if ctx.params['debug'] and ctx.params['release']:
            raise ValueError("Set at most one of `--debug` and `--release`!")
        if ctx.params['debug']:
            buildtype = 'debug'
            cflags_unwanted = ('-O1', '-O2', '-O3')
        elif ctx.params['release']:
            buildtype = 'release'
            cflags_unwanted = ('-O0', '-O1', '-O2')
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + (f"-Dbuildtype={buildtype}", )
        if 'CFLAGS' in os.environ.keys():
            # Check that CFLAGS doesn't contain something that supercedes -O0
            # for a plain debug build (conda envs tend to set -O2)
            cflags = os.environ['CFLAGS'].split()
            for flag in cflags_unwanted:
                if flag in cflags:
                    raise ValueError(f"A {buildtype} build isn't possible, "
                                        f"because CFLAGS contains `{flag}`."
                                        "Please also check CXXFLAGS and FFLAGS.")
    if ctx.params['gcov']:
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + ('-Db_coverage=true', )

    if ctx.params['asan']:
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + ('-Db_sanitize=address,undefined', )

    if ctx.params['setup_args']:
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + tuple([str(arg) for arg in ctx.params['setup_args']])

    if ctx.params['with_accelerate']:
        # on a mac you probably want to use accelerate over scipy_openblas
        ctx.params[MESON_ARGS]["setup"] = ctx.params[MESON_ARGS]["setup"] + ("-Dblas=accelerate", )
    elif ctx.params['with_scipy_openblas']:
        configure_scipy_openblas()
        os.env['PKG_CONFIG_PATH'] = os.pathsep.join([
                os.getcwd(),
                os.env.get('PKG_CONFIG_PATH', '')
                ])

    if ctx.params['parallel'] is None:
        # Use number of physical cores rather than ninja's default of 2N+2,
        # to avoid out of memory issues (see gh-17941 and gh-18443)
        n_cores = cpu_count(only_physical_cores=True)
        ctx.params[MESON_ARGS]["jobs"] = n_cores
    else:
        ctx.params[MESON_ARGS]["jobs"] = ctx.params['parallel']

    ctx.params[MESON_ARGS]["install"] = ctx.params[MESON_ARGS]["install"] + ("--tags=" + ctx.params['tags'], )

    if ctx.params["show_build_log"]:
        ctx.params["verbose"] = ctx.params["show_build_log"]

    for option in ('werror', 'debug', 'release',
                   'gcov', 'asan','setup_args',
                   'with_accelerate', 'with_scipy_openblas',
                   'parallel', 'show_build_log',
                   'tags'):
        ctx.params.pop(option)

    ctx.forward(meson.build)

@click.command()
@click.argument("pytest_args", nargs=-1)
@click.option(
    '--verbose', '-v', default=False, is_flag=True,
    help="more verbosity")
@click.option(
    '--coverage', '-c', default=False, is_flag=True,
    help=("report coverage of project code. "
            "HTML output goes under build/coverage")) # removed doctests as currently not supported by _lib/_testutils.py
                                                    # doctests = Option(['--doctests'], default=False)
@click.option(
    '--durations', '-d', default=None, metavar="NUM_TESTS",
    help="Show timing for the given number of slowest tests"
)
@click.option(
    '--submodule', '-s', default=None, metavar='MODULE_NAME',
    help="Submodule whose tests to run (cluster, constants, ...)")
@click.option(
    '--tests', '-t', default=None, metavar='TESTS',
    help='Specify tests to run')
@click.option(
    '--mode', '-m', default='not slow', metavar='MODE', show_default=True,
    help=("'fast', 'full', or something that could be passed to "
            "`pytest -m` as a marker expression"))
@click.option(
    '--parallel', '-j', default=1, metavar='N_JOBS',
    help="Number of parallel jobs for testing"
)
@click.option(
    '--array-api-backend', '-b', default=None, metavar='ARRAY_BACKEND',
    multiple=True,
    help=(
        "Array API backend "
        "('all', 'numpy', 'torch', 'cupy', 'array_api_strict', 'jax.numpy')."
    )
)
@click.pass_context
def test(ctx, pytest_args, verbose, *args, **kwargs):
    """🔧 Run tests

    PYTEST_ARGS are passed through directly to pytest, e.g.:

      spin test -- --pdb

    To run tests on a directory or file:

     \b
     spin test numpy/linalg
     spin test numpy/linalg/tests/test_linalg.py

    To report the durations of the N slowest tests:

      spin test -- --durations=N

    To run tests that match a given pattern:

     \b
     spin test -- -k "geometric"
     spin test -- -k "geometric and not rgeometric"

    By default, spin will run `-m 'not slow'`. To run the full test suite, use
    `spin test -m full`

    For more, see `pytest --help`.
    """  # noqa: E501
    tests = ctx.params['tests']
    if ctx.params["submodule"]:
        tests = PROJECT_MODULE + "." + ctx.params["submodule"]

    markexpr = ctx.params['mode']
    if (not pytest_args) and (not tests):
        pytest_args = ('scipy',)

    if '-m' not in pytest_args:
        if len(pytest_args) == 1 and not tests:
            tests = pytest_args[0]
            pytest_args = ()
        if markexpr != "full":
            pytest_args = ('-m', markexpr) + pytest_args

    n_jobs = ctx.params['parallel']
    if (n_jobs != 1) and ('-n' not in pytest_args):
        pytest_args = ('-n', str(n_jobs)) + pytest_args

    if tests and '--pyargs' not in pytest_args:
        pytest_args += ('--pyargs', tests)

    if verbose:
        pytest_args = ('-v',) + pytest_args

    ctx.params['pytest_args'] = pytest_args

    if len(ctx.params['array_api_backend']) != 0:
        os.environ['SCIPY_ARRAY_API'] = json.dumps(list(ctx.params['array_api_backend']))

    for extra_param in (
        'verbose', 'coverage', 'durations',
        'submodule', 'array_api_backend',
        'mode', 'parallel', 'tests'):
        del ctx.params[extra_param]

    ctx.forward(meson.test)

@click.command()
@click.option(
    "--clean", is_flag=True,
    default=False,
    help="Clean previously built docs before building"
)
@click.option(
    "--no-build", "-n",
    "first_build",
    default=True,
    is_flag=True,
    help="Build numpy before generating docs",
)
@click.option(
        '--list-targets', '-t', default=False, is_flag=True,
        help='List doc targets',
    )
@click.option(
        '--parallel', '-j', default="auto", metavar='N_JOBS',
        help="Number of parallel jobs"
    )
@click.option(
    '--no-cache', default=False, is_flag=True,
    help="Forces a full rebuild of the docs. Note that this may be " + \
            "needed in order to make docstring changes in C/Cython files " + \
            "show up."
)
@click.pass_context
def docs(ctx, list_targets, clean, first_build, parallel, *args, **kwargs):
    """📖 Build Sphinx documentation

    By default, SPHINXOPTS="-W", raising errors on warnings.
    To build without raising on warnings:

      SPHINXOPTS="" spin docs

    To list all Sphinx targets:

      spin docs targets

    To build another Sphinx target:

      spin docs TARGET

    E.g., to build a zipfile of the html docs for distribution:

      spin docs dist

    """
    meson.docs.ignore_unknown_options = True

    if clean:
        cwd = os.getcwd()
        os.chdir(os.path.join(cwd, "doc"))
        subprocess.call(["make", "clean"], cwd=os.getcwd())
        ctx.params.pop("clean")
        os.chdir(cwd)

    SPHINXOPTS = "-W"
    if "no_cache" in ctx.params:
        if ctx.params["no_cache"]:
            SPHINXOPTS += " -E"
        ctx.params.pop("no_cache")
    if "parallel" in ctx.params:
        ctx.params["jobs"] = ctx.params["parallel"]
        ctx.params.pop("parallel")
    SPHINXOPTS = os.environ.get("SPHINXOPTS", "") + SPHINXOPTS
    os.environ["SPHINXOPTS"] = SPHINXOPTS

    ctx.params.pop("list_targets")
    ctx.params["sphinx_target"] = "html"

    ctx.forward(meson.docs)
