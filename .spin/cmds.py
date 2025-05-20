import contextlib
import os
import sys
import importlib
import importlib.util
import json
import traceback
import warnings
import math
import subprocess
from concurrent.futures.process import _MAX_WINDOWS_WORKERS

import spin
import click
from spin import util
from spin.cmds import meson

from pathlib import Path

PROJECT_MODULE = "scipy"

@click.option(
    '--werror', default=False, is_flag=True,
    help="Treat warnings as errors")
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
    '--use-system-libraries', default=False, is_flag=True,
    help=("If set, use system libraries"
            "if they are available for subprojects."))
@click.option(
    '--tags', default="runtime,python-runtime,tests,devel",
    show_default=True, help="Install tags to be used by meson."
)
@spin.util.extend_command(spin.cmds.meson.build)
def build(*, parent_callback, meson_args, jobs, verbose, werror, asan, debug,
          release, parallel, setup_args, show_build_log,
          with_scipy_openblas, with_accelerate, use_system_libraries,
          tags, **kwargs):
    """🔧 Build package with Meson/ninja and install

    MESON_ARGS are passed through e.g.:

    spin build -- -Dpkg_config_path=/lib64/pkgconfig

    The package is installed to build-install

    By default builds for release, to be able to use a debugger set CFLAGS
    appropriately. For example, for linux use

    CFLAGS="-O0 -g" spin build
    """
    MESON_ARGS = "meson_args"
    MESON_COMPILE_ARGS = "meson_compile_args"
    MESON_INSTALL_ARGS = "meson_install_args"

    meson_compile_args = tuple()
    meson_install_args = tuple()

    if sys.platform == "cygwin":
        # Cygwin only has netlib lapack, but can link against
        # OpenBLAS rather than netlib blas at runtime.  There is
        # no libopenblas-devel to enable linking against
        # openblas-specific functions or OpenBLAS Lapack
        meson_args = meson_args + ("-Dlapack=lapack", "-Dblas=blas")

    if werror:
        meson_args = meson_args + ("--werror", )

    if debug or release:
        if debug and release:
            raise ValueError("Set at most one of `--debug` and `--release`!")
        if debug:
            buildtype = 'debug'
            cflags_unwanted = ('-O1', '-O2', '-O3')
        elif release:
            buildtype = 'release'
            cflags_unwanted = ('-O0', '-O1', '-O2')
        meson_args = meson_args + (f"-Dbuildtype={buildtype}", )
        if 'CFLAGS' in os.environ.keys():
            # Check that CFLAGS doesn't contain something that supercedes -O0
            # for a plain debug build (conda envs tend to set -O2)
            cflags = os.environ['CFLAGS'].split()
            for flag in cflags_unwanted:
                if flag in cflags:
                    raise ValueError(f"A {buildtype} build isn't possible, "
                                        f"because CFLAGS contains `{flag}`."
                                        "Please also check CXXFLAGS and FFLAGS.")

    if asan:
        meson_args = meson_args + ('-Db_sanitize=address,undefined', )

    if setup_args:
        meson_args = meson_args + tuple([str(arg) for arg in setup_args])

    if with_accelerate:
        # on a mac you probably want to use accelerate over scipy_openblas
        meson_args = meson_args + ("-Dblas=accelerate", )
    elif with_scipy_openblas:
        configure_scipy_openblas()
        os.environ['PKG_CONFIG_PATH'] = os.pathsep.join([
                os.getcwd(),
                os.environ.get('PKG_CONFIG_PATH', '')
                ])

    if use_system_libraries:
        meson_args = meson_args + ("-Duse-system-libraries=auto",)

    if parallel is None:
        # Use number of physical cores rather than ninja's default of 2N+2,
        # to avoid out of memory issues (see gh-17941 and gh-18443)
        n_cores = cpu_count(only_physical_cores=True)
        jobs = n_cores
    else:
        jobs = parallel

    meson_install_args = meson_install_args + ("--tags=" + tags, )

    if show_build_log:
        verbose = show_build_log

    parent_callback(**{MESON_ARGS: meson_args,
                       MESON_COMPILE_ARGS: meson_compile_args,
                       MESON_INSTALL_ARGS: meson_install_args,
                       "jobs": jobs,
                       "verbose": verbose,
                       **kwargs})

@click.option(
    '--durations', '-d', default=None, metavar="NUM_TESTS",
    help="Show timing for the given number of slowest tests"
)
@click.option(
    '--submodule', '-s', default=None, metavar='MODULE_NAME',
    help="Submodule whose tests to run (cluster, constants, ...)")
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
        "('all', 'numpy', 'torch', 'cupy', 'array_api_strict', "
        "'jax.numpy', 'dask.array')."
    )
)
@spin.util.extend_command(spin.cmds.meson.test)
def test(*, parent_callback, pytest_args, tests, coverage,
         durations, submodule, mode, parallel,
         array_api_backend, **kwargs):
    """🔧 Run tests

    PYTEST_ARGS are passed through directly to pytest, e.g.:

      spin test -- --pdb

    To run tests on a directory or file:

     \b
     spin test scipy/linalg

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

    build_dir = os.path.abspath(kwargs['build_dir'])
    site_package_dir = get_site_packages(build_dir)

    if coverage:
        if is_editable_install():
            click.secho(
                "Error: cannot generate coverage report for editable installs",
                fg="bright_red",
            )
            raise SystemExit(1)
        elif site_package_dir is None:
            raise FileNotFoundError(
                "SciPy build not found, please execute "
                "``spin build`` before calling ``spin test --coverage``. "
                "We need it to figure out whether ``lcov`` can be called or not.")
        else:
            # Check needed to ensure gcov functions correctly.
            with working_dir(site_package_dir):
                sys.path.insert(0, site_package_dir)
                os.environ['PYTHONPATH'] = os.pathsep.join(
                        (site_package_dir, os.environ.get('PYTHONPATH', '')))
                was_built_with_gcov_flag = len(list(
                    Path(build_dir).rglob("*.gcno"))) > 0
                if was_built_with_gcov_flag:
                    config = importlib.import_module(
                            "scipy.__config__").show(mode='dicts')
                    compilers_config = config['Compilers']
                    cpp = compilers_config['c++']['name']
                    c = compilers_config['c']['name']
                    fortran = compilers_config['fortran']['name']
                    if not (c == 'gcc' and cpp == 'gcc' and fortran == 'gcc'):
                        print("SciPy was built with --gcov flag which requires "
                            "LCOV while running tests.\nFurther, LCOV usage "
                            "requires GCC for C, C++ and Fortran codes in SciPy.\n"
                            "Compilers used currently are:\n"
                            f"  C: {c}\n  C++: {cpp}\n  Fortran: {fortran}\n"
                            "Therefore, exiting without running tests.")
                        exit(1) # Exit because tests will give missing symbol error

    if submodule:
        tests = PROJECT_MODULE + "." + submodule

    markexpr = mode
    if (not pytest_args) and (not tests):
        pytest_args = ('scipy',)

    if '-m' not in pytest_args:
        if len(pytest_args) == 1 and not tests:
            tests = pytest_args[0]
            pytest_args = ()
        if markexpr != "full":
            pytest_args = ('-m', markexpr) + pytest_args

    n_jobs = parallel
    if (n_jobs != 1) and ('-n' not in pytest_args):
        pytest_args = ('-n', str(n_jobs)) + pytest_args

    if durations:
        pytest_args += ('--durations', durations)

    if len(array_api_backend) != 0:
        os.environ['SCIPY_ARRAY_API'] = json.dumps(list(array_api_backend))

    parent_callback(**{"pytest_args": pytest_args, "tests": tests,
                    "coverage": coverage, **kwargs})

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
@spin.util.extend_command(spin.cmds.meson.docs)
def docs(*, parent_callback, sphinx_target, clean, jobs,
         list_targets, parallel, no_cache, **kwargs):
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

    if clean: # SciPy has its own mechanism to clear the previous docs build
        cwd = os.getcwd()
        os.chdir(os.path.join(cwd, "doc"))
        subprocess.call(["make", "clean"], cwd=os.getcwd())
        clean = False
        os.chdir(cwd)

    SPHINXOPTS = "-W"
    if no_cache:
        SPHINXOPTS += " -E"

    jobs = parallel
    SPHINXOPTS = os.environ.get("SPHINXOPTS", "") + SPHINXOPTS
    os.environ["SPHINXOPTS"] = SPHINXOPTS

    sphinx_target = "html"

    parent_callback(**{"sphinx_target": sphinx_target,
                       "clean": clean, "jobs": jobs, **kwargs})

def _set_pythonpath(pythonpath):
    env = os.environ
    env['PYTHONWARNINGS'] = env.get('PYTHONWARNINGS', 'all')

    if pythonpath:
        for p in reversed(pythonpath.split(os.pathsep)):
            sys.path.insert(0, p)

@click.option(
    '--pythonpath', '-p', metavar='PYTHONPATH', default=None,
    help='Paths to prepend to PYTHONPATH')
@spin.util.extend_command(spin.cmds.meson.python)
def python(*, parent_callback, pythonpath, **kwargs):
    """🐍 Launch Python shell with PYTHONPATH set

    OPTIONS are passed through directly to Python, e.g.:

    spin python -c 'import sys; print(sys.path)'
    """
    _set_pythonpath(pythonpath)
    parent_callback(**kwargs)

@click.option(
    '--pythonpath', '-p', metavar='PYTHONPATH', default=None,
    help='Paths to prepend to PYTHONPATH')
@spin.util.extend_command(spin.cmds.meson.ipython)
def ipython(*, parent_callback, pythonpath, **kwargs):
    """💻 Launch IPython shell with PYTHONPATH set

    OPTIONS are passed through directly to IPython, e.g.:

    spin ipython -i myscript.py
    """
    _set_pythonpath(pythonpath)
    parent_callback(**kwargs)

@click.option(
    '--pythonpath', '-p', metavar='PYTHONPATH', default=None,
    help='Paths to prepend to PYTHONPATH')
@spin.util.extend_command(spin.cmds.meson.shell)
def shell(*, parent_callback, pythonpath, **kwargs):
    """💻 Launch shell with PYTHONPATH set

    SHELL_ARGS are passed through directly to the shell, e.g.:

    spin shell -- -c 'echo $PYTHONPATH'

    Ensure that your shell init file (e.g., ~/.zshrc) does not override
    the PYTHONPATH.
    """
    _set_pythonpath(pythonpath)
    parent_callback(**kwargs)

@contextlib.contextmanager
def working_dir(new_dir):
    current_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(current_dir)

@click.command(context_settings={"ignore_unknown_options": True})
@meson.build_dir_option
@click.pass_context
def mypy(ctx, build_dir=None):
    """🦆 Run Mypy tests for SciPy
    """
    if is_editable_install():
        click.secho(
            "Error: Mypy does not work (well) for editable installs",
            fg="bright_red",
        )
        raise SystemExit(1)
    else:
        click.secho(
                "Invoking `build` prior to running mypy tests:",
                bold=True, fg="bright_green"
            )
        ctx.invoke(build)

    try:
        import mypy.api
    except ImportError as e:
        raise RuntimeError(
            "Mypy not found. Please install it by running "
            "pip install -r mypy_requirements.txt from the repo root"
        ) from e

    build_dir = os.path.abspath(build_dir)
    root = Path(build_dir).parent
    config = os.path.join(root, "mypy.ini")
    check_path = PROJECT_MODULE
    install_dir = meson._get_site_packages(build_dir)

    with working_dir(install_dir):
        os.environ['MYPY_FORCE_COLOR'] = '1'
        click.secho(f"mypy.api.run --config-file {config} {check_path}",
                    bold=True, fg="bright_blue")
        report, errors, status = mypy.api.run([
            "--config-file",
            str(config),
            check_path,
        ])
    print(report, end='')
    print(errors, end='', file=sys.stderr)

@spin.util.extend_command(test, doc='')
def smoke_docs(*, parent_callback, pytest_args, **kwargs):
    """🔧 Run doctests of objects in the public API.

    PYTEST_ARGS are passed through directly to pytest, e.g.:

      spin smoke-docs -- --pdb

    To run tests on a directory:

     \b
     spin smoke-docs scipy/linalg

    To report the durations of the N slowest doctests:

      spin smoke-docs -- --durations=N

    To run doctests that match a given pattern:

     \b
     spin smoke-docs -- -k "slogdet"
     spin smoke-docs scipy/linalg -- -k "det and not slogdet"

    \b
    Note:
    -----

    \b
     - This command only runs doctests and skips everything under tests/
     - This command only doctests public objects: those which are accessible
       from the top-level `__init__.py` file.

    """  # noqa: E501
    # prevent obscure error later; cf https://github.com/numpy/numpy/pull/26691/
    if not importlib.util.find_spec("scipy_doctest"):
        raise ModuleNotFoundError("Please install scipy-doctest")

    tests = kwargs["tests"]
    if kwargs["submodule"]:
        tests = PROJECT_MODULE + "." + kwargs["submodule"]

    if not pytest_args and not tests:
        pytest_args = ('scipy', )

    # turn doctesting on:
    doctest_args = (
        '--doctest-modules',
        '--doctest-collect=api'
    )

    if not tests:
        doctest_args += ('--doctest-collect=api', )

    pytest_args = pytest_args + doctest_args

    parent_callback(**{"pytest_args": pytest_args, **kwargs})

@click.command()
@click.option(
    '--verbose', '-v', default=False, is_flag=True,
    help="more verbosity")
@click.option(
    '--submodule', '-s', default=None, metavar='MODULE_NAME',
    help="Submodule whose tests to run (cluster, constants, ...)")
@meson.build_dir_option
@click.pass_context
def refguide_check(ctx, build_dir=None, *args, **kwargs):
    """🔧 Run refguide check."""
    click.secho(
            "Invoking `build` prior to running refguide-check:",
            bold=True, fg="bright_green"
        )
    ctx.invoke(build)

    build_dir = os.path.abspath(build_dir)
    root = Path(build_dir).parent
    install_dir = meson._get_site_packages(build_dir)

    cmd = [f'{sys.executable}',
            os.path.join(root, 'tools', 'refguide_check.py')]

    if ctx.params["verbose"]:
        cmd += ['-vvv']

    if ctx.params["submodule"]:
        cmd += [ctx.params["submodule"]]

    os.environ['PYTHONPATH'] = install_dir
    util.run(cmd)

@click.command()
@click.argument(
    'pytest_args', nargs=-1, metavar='PYTEST-ARGS', required=False
)
@click.option(
    '--tests', '-t', default=None, multiple=True, metavar='TESTS',
    help='Specify *rst files to smoke test')
@click.option(
    '--verbose', '-v', default=False, is_flag=True, help="verbosity")
@meson.build_dir_option
@click.pass_context
def smoke_tutorials(ctx, pytest_args, tests, verbose, build_dir, *args, **kwargs):
    """🔧 Run doctests of user-facing rst tutorials.

    To test all tutorials in the scipy doc/source/tutorial directory, use

      spin smoke-tutorials

    To run tests on a specific RST file:

     \b
     spin smoke-tutorials doc/source/reference/stats.rst
     spin smoke-tutorials -t doc/source/reference/stats.rst

    \b
    Note:
    -----

    \b
     - This command only runs doctests and skips everything under tests/
     - This command only doctests public objects: those which are accessible
       from the top-level `__init__.py` file.

    """  # noqa: E501

    click.secho(
        "Invoking `build` prior to running tests for tutorials:",
        bold=True, fg="bright_green"
    )
    ctx.invoke(build)

    meson._set_pythonpath(build_dir)

    cmd = ['pytest']
    if tests:
        cmd += list(tests)
    else:
        cmd += ['doc/source/tutorial', '--doctest-glob=*rst']
    if verbose:
        cmd += ['-v']

    extra_argv = list(pytest_args[:]) if pytest_args else []
    if extra_argv and extra_argv[0] == '--':
        extra_argv = extra_argv[1:]
    cmd += extra_argv

    cmd_str = ' '.join(cmd)
    click.secho(cmd_str, bold=True, fg="bright_blue")
    util.run(cmd)

@click.command()
@click.option(
    '--fix', default=False, is_flag=True,
    help='Attempt to auto-fix errors')
@click.option("--diff-against", default="main", help="Diff against "
    "this branch and lint modified files. Use either "
    "`--diff-against` or `--files`, but not both.")
@click.option("--files", default="",
    help="Lint these files or directories; "
         "use **/*.py to lint all files")
@click.option("--all", default=False, is_flag=True,
    help="This overrides `--diff-against` and `--files` "
         "to lint all local files (excluding subprojects).")
@click.option("--no-cython", default=True, is_flag=True,
    help="Do not run cython-lint.")
@click.pass_context
def lint(ctx, fix, diff_against, files, all, no_cython):
    """🔦 Run linter on modified files and check for
    disallowed Unicode characters and possibly-invalid test names."""
    cmd_prefix = [sys.executable] if sys.platform == "win32" else []

    cmd_lint = cmd_prefix + [
        os.path.join('tools', 'lint.py'),
        f'--diff-against={diff_against}'
    ]
    if files != "":
        cmd_lint += [f'--files={files}']
    if all:
        cmd_lint += ['--all']
    if no_cython:
        cmd_lint += ['--no-cython']
    if fix:
        cmd_lint += ['--fix']
    util.run(cmd_lint)

    cmd_unicode = cmd_prefix + [
        os.path.join('tools', 'check_unicode.py')
    ]
    util.run(cmd_unicode)

    cmd_check_test_name = cmd_prefix + [
        os.path.join('tools', 'check_test_name.py')
    ]
    util.run(cmd_check_test_name)

# From scipy: benchmarks/benchmarks/common.py
def _set_mem_rlimit(max_mem=None):
    """
    Set address space rlimit
    """
    import resource
    import psutil

    mem = psutil.virtual_memory()

    if max_mem is None:
        max_mem = int(mem.total * 0.7)
    cur_limit = resource.getrlimit(resource.RLIMIT_AS)
    if cur_limit[0] > 0:
        max_mem = min(max_mem, cur_limit[0])

    try:
        resource.setrlimit(resource.RLIMIT_AS, (max_mem, cur_limit[1]))
    except ValueError:
        # on macOS may raise: current limit exceeds maximum limit
        pass

def _run_asv(cmd):
    # Always use ccache, if installed
    PATH = os.environ['PATH']
    EXTRA_PATH = os.pathsep.join([
        '/usr/lib/ccache', '/usr/lib/f90cache',
        '/usr/local/lib/ccache', '/usr/local/lib/f90cache'
    ])
    env = os.environ
    env['PATH'] = f'{EXTRA_PATH}{os.pathsep}{PATH}'

    # Control BLAS/LAPACK threads
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'

    # Limit memory usage
    try:
        _set_mem_rlimit()
    except (ImportError, RuntimeError):
        pass

    util.run(cmd, cwd='benchmarks', env=env)

def _commit_to_sha(commit):
    p = util.run(['git', 'rev-parse', commit], output=False, echo=False)
    if p.returncode != 0:
        raise(
            click.ClickException(
                f'Could not find SHA matching commit `{commit}`'
            )
        )

    return p.stdout.decode('ascii').strip()


def _dirty_git_working_dir():
    # Changes to the working directory
    p0 = util.run(['git', 'diff-files', '--quiet'])

    # Staged changes
    p1 = util.run(['git', 'diff-index', '--quiet', '--cached', 'HEAD'])

    return (p0.returncode != 0 or p1.returncode != 0)

@click.command()
@click.option(
    '--tests', '-t',
    default=None, metavar='TESTS', multiple=True,
    help="Which tests to run"
)
@click.option(
    '--submodule', '-s', default=None, metavar='SUBMODULE',
    help="Submodule whose tests to run (cluster, constants, ...)")
@click.option(
    '--compare', '-c',
    is_flag=True,
    default=False,
    help="Compare benchmarks between the current branch and main "
         "(unless other branches specified). "
         "The benchmarks are each executed in a new isolated "
         "environment."
)
@click.option(
    '--verbose', '-v', is_flag=True, default=False
)
@click.option(
    '--quick', '-q', is_flag=True, default=False,
    help="Run each benchmark only once (timings won't be accurate)"
)
@click.argument(
    'commits', metavar='',
    required=False,
    nargs=-1
)
@meson.build_dir_option
@click.pass_context
def bench(ctx, tests, submodule, compare, verbose, quick,
          commits, build_dir=None, *args, **kwargs):
    """🔧 Run benchmarks.

    \b
    ```python
     Examples:

    $ spin bench -t integrate.SolveBVP
    $ spin bench -t linalg.Norm
    $ spin bench --compare main
    ```
    """
    build_dir = os.path.abspath(build_dir)
    if not commits:
        commits = ('main', 'HEAD')
    elif len(commits) == 1:
        commits = commits + ('HEAD',)
    elif len(commits) > 2:
        raise click.ClickException(
            'Need a maximum of two revisions to compare'
        )

    bench_args = []
    if submodule:
        submodule = (submodule, )
    else:
        submodule = tuple()
    for t in tests + submodule:
        bench_args += ['--bench', t]

    if verbose:
        bench_args = ['-v'] + bench_args

    if quick:
        bench_args = ['--quick'] + bench_args

    if not compare:
        # No comparison requested; we build and benchmark the current version

        click.secho(
            "Invoking `build` prior to running benchmarks:",
            bold=True, fg="bright_green"
        )
        ctx.invoke(build)

        meson._set_pythonpath(build_dir)

        p = util.run(
            ['python', '-c', 'import scipy as sp; print(sp.__version__)'],
            cwd='benchmarks',
            echo=False,
            output=False
        )
        os.chdir('..')

        np_ver = p.stdout.strip().decode('ascii')
        click.secho(
            f'Running benchmarks on SciPy {np_ver}',
            bold=True, fg="bright_green"
        )
        cmd = [
            'asv', 'run', '--dry-run',
            '--show-stderr', '--python=same',
            '--quick'] + bench_args
        _run_asv(cmd)
    else:
        # Ensure that we don't have uncommited changes
        commit_a, commit_b = [_commit_to_sha(c) for c in commits]

        if commit_b == 'HEAD' and _dirty_git_working_dir():
            click.secho(
                "WARNING: you have uncommitted changes --- "
                "these will NOT be benchmarked!",
                fg="red"
            )

        cmd_compare = [
            'asv', 'continuous', '--factor', '1.05', '--quick'
        ] + bench_args + [commit_a, commit_b]
        _run_asv(cmd_compare)


def configure_scipy_openblas(blas_variant='32'):
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

    Note that on Windows, the returned number of CPUs cannot exceed 61, see:
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

def get_site_packages(build_dir):
    """site-packages directory is path to installed in-tree build.

    Returns None if `scipy` wasn't build at all.
    Returns an empty string (from spin.meson call) for an editable install.
    """
    try:
        return meson._get_site_packages(build_dir)
    except FileNotFoundError:
        return None


def is_editable_install():
    return meson._is_editable_install_of_same_source('scipy')
