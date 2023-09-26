#! /usr/bin/env python3

'''
Developer CLI: building (meson), tests, benchmark, etc.

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authoring the development tasks.

REQUIREMENTS:
--------------
- see environment.yml: doit, pydevtool, click, rich-click

# USAGE:

## 1 - click API

Commands can added using default Click API. i.e.

```
@cli.command()
@click.argument('extra_argv', nargs=-1)
@click.pass_obj
def python(ctx_obj, extra_argv):
    """Start a Python shell with PYTHONPATH set"""
```

## 2 - class based Click command definition

`CliGroup` provides an alternative class based API to create Click commands.

Just use the `cls_cmd` decorator. And define a `run()` method

```
@cli.cls_cmd('test')
class Test():
    """Run tests"""

    @classmethod
    def run(cls):
        print('Running tests...')
```

- Command may make use a Click.Group context defining a `ctx` class attribute
- Command options are also define as class attributes

```
@cli.cls_cmd('test')
class Test():
    """Run tests"""
    ctx = CONTEXT

    verbose = Option(
        ['--verbose', '-v'], default=False, is_flag=True, help="verbosity")

    @classmethod
    def run(cls, **kwargs): # kwargs contains options from class and CONTEXT
        print('Running tests...')
```

## 3 - class based interface can be run as a doit task by subclassing from Task

- Extra doit task metadata can be defined as class attribute TASK_META.
- `run()` method will be used as python-action by task

```
@cli.cls_cmd('test')
class Test(Task):   # Task base class, doit will create a task
    """Run tests"""
    ctx = CONTEXT

    TASK_META = {
        'task_dep': ['build'],
    }

    @classmethod
    def run(cls, **kwargs):
        pass
```

## 4 - doit tasks with cmd-action "shell" or dynamic metadata

Define method `task_meta()` instead of `run()`:

```
@cli.cls_cmd('refguide-check')
class RefguideCheck(Task):
    @classmethod
    def task_meta(cls, **kwargs):
        return {
```

'''

import os
import subprocess
import sys
import warnings
import shutil
import json
import datetime
import time
import platform
import importlib.util
import errno
import contextlib
from sysconfig import get_path
import math
import traceback
from concurrent.futures.process import _MAX_WINDOWS_WORKERS

from pathlib import Path
from collections import namedtuple
from types import ModuleType as new_module
from dataclasses import dataclass

import click
from click import Option, Argument
from doit.cmd_base import ModuleTaskLoader
from doit.reporter import ZeroReporter
from doit.exceptions import TaskError
from doit.api import run_tasks
from pydevtool.cli import UnifiedContext, CliGroup, Task
from rich.console import Console
from rich.panel import Panel
from rich.theme import Theme
from rich_click import rich_click

DOIT_CONFIG = {
    'verbosity': 2,
    'minversion': '0.36.0',
}


console_theme = Theme({
    "cmd": "italic gray50",
})

if sys.platform == 'win32':
    class EMOJI:
        cmd = ">"
else:
    class EMOJI:
        cmd = ":computer:"


rich_click.STYLE_ERRORS_SUGGESTION = "yellow italic"
rich_click.SHOW_ARGUMENTS = True
rich_click.GROUP_ARGUMENTS_OPTIONS = False
rich_click.SHOW_METAVARS_COLUMN = True
rich_click.USE_MARKDOWN = True
rich_click.OPTION_GROUPS = {
    "dev.py": [
        {
            "name": "Options",
            "options": [
                "--help", "--build-dir", "--no-build", "--install-prefix"],
        },
    ],

    "dev.py test": [
        {
            "name": "Options",
            "options": ["--help", "--verbose", "--parallel", "--coverage",
                        "--durations"],
        },
        {
            "name": "Options: test selection",
            "options": ["--submodule", "--tests", "--mode"],
        },
    ],
}
rich_click.COMMAND_GROUPS = {
    "dev.py": [
        {
            "name": "build & testing",
            "commands": ["build", "test"],
        },
        {
            "name": "static checkers",
            "commands": ["lint", "mypy"],
        },
        {
            "name": "environments",
            "commands": ["shell", "python", "ipython", "show_PYTHONPATH"],
        },
        {
            "name": "documentation",
            "commands": ["doc", "refguide-check"],
        },
        {
            "name": "release",
            "commands": ["notes", "authors"],
        },
        {
            "name": "benchmarking",
            "commands": ["bench"],
        },
    ]
}


class ErrorOnlyReporter(ZeroReporter):
    desc = """Report errors only"""

    def runtime_error(self, msg):
        console = Console()
        console.print("[red bold] msg")

    def add_failure(self, task, fail_info):
        console = Console()
        if isinstance(fail_info, TaskError):
            console.print(f'[red]Task Error - {task.name}'
                          f' => {fail_info.message}')
        if fail_info.traceback:
            console.print(Panel(
                "".join(fail_info.traceback),
                title=f"{task.name}",
                subtitle=fail_info.message,
                border_style="red",
            ))


CONTEXT = UnifiedContext({
    'build_dir': Option(
        ['--build-dir'], metavar='BUILD_DIR',
        default='build', show_default=True,
        help=':wrench: Relative path to the build directory.'),
    'no_build': Option(
        ["--no-build", "-n"], default=False, is_flag=True,
        help=(":wrench: Do not build the project"
              " (note event python only modification require build).")),
    'install_prefix': Option(
        ['--install-prefix'], default=None, metavar='INSTALL_DIR',
        help=(":wrench: Relative path to the install directory."
              " Default is <build-dir>-install.")),
})


def run_doit_task(tasks):
    """
      :param tasks: (dict) task_name -> {options}
    """
    loader = ModuleTaskLoader(globals())
    doit_config = {
        'verbosity': 2,
        'reporter': ErrorOnlyReporter,
    }
    return run_tasks(loader, tasks, extra_config={'GLOBAL': doit_config})


class CLI(CliGroup):
    context = CONTEXT
    run_doit_task = run_doit_task


@click.group(cls=CLI)
@click.pass_context
def cli(ctx, **kwargs):
    """Developer Tool for SciPy

    \bCommands that require a built/installed instance are marked with :wrench:.


    \b**python dev.py --build-dir my-build test -s stats**

    """  # noqa: E501
    CLI.update_context(ctx, kwargs)


PROJECT_MODULE = "scipy"
PROJECT_ROOT_FILES = ['scipy', 'LICENSE.txt', 'meson.build']


@dataclass
class Dirs:
    """
        root:
            Directory where scr, build config and tools are located
            (and this file)
        build:
            Directory where build output files (i.e. *.o) are saved
        install:
            Directory where .so from build and .py from src are put together.
        site:
            Directory where the built SciPy version was installed.
            This is a custom prefix, followed by a relative path matching
            the one the system would use for the site-packages of the active
            Python interpreter.
    """
    # all paths are absolute
    root: Path
    build: Path
    installed: Path
    site: Path  # <install>/lib/python<version>/site-packages

    def __init__(self, args=None):
        """:params args: object like Context(build_dir, install_prefix)"""
        self.root = Path(__file__).parent.absolute()
        if not args:
            return

        self.build = Path(args.build_dir).resolve()
        if args.install_prefix:
            self.installed = Path(args.install_prefix).resolve()
        else:
            self.installed = self.build.parent / (self.build.stem + "-install")

        if sys.platform == 'win32' and sys.version_info < (3, 10):
            # Work around a pathlib bug; these must be absolute paths
            self.build = Path(os.path.abspath(self.build))
            self.installed = Path(os.path.abspath(self.installed))

        # relative path for site-package with py version
        # i.e. 'lib/python3.10/site-packages'
        self.site = self.get_site_packages()

    def add_sys_path(self):
        """Add site dir to sys.path / PYTHONPATH"""
        site_dir = str(self.site)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))

    def get_site_packages(self):
        """
        Depending on whether we have debian python or not,
        return dist_packages path or site_packages path.
        """
        if sys.version_info >= (3, 12):
            plat_path = Path(get_path('platlib'))
        else:
            # distutils is required to infer meson install path
            # for python < 3.12 in debian patched python
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=DeprecationWarning)
                from distutils import dist
                from distutils.command.install import INSTALL_SCHEMES
            if 'deb_system' in INSTALL_SCHEMES:
                # debian patched python in use
                install_cmd = dist.Distribution().get_command_obj('install')
                install_cmd.select_scheme('deb_system')
                install_cmd.finalize_options()
                plat_path = Path(install_cmd.install_platlib)
            else:
                plat_path = Path(get_path('platlib'))
        return self.installed / plat_path.relative_to(sys.exec_prefix)


@contextlib.contextmanager
def working_dir(new_dir):
    current_dir = os.getcwd()
    try:
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(current_dir)


def import_module_from_path(mod_name, mod_path):
    """Import module with name `mod_name` from file path `mod_path`"""
    spec = importlib.util.spec_from_file_location(mod_name, mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def get_test_runner(project_module):
    """
    get Test Runner from locally installed/built project
    """
    __import__(project_module)
    # scipy._lib._testutils:PytestTester
    test = sys.modules[project_module].test
    version = sys.modules[project_module].__version__
    mod_path = sys.modules[project_module].__file__
    mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))
    return test, version, mod_path


############

@cli.cls_cmd('build')
class Build(Task):
    """:wrench: Build & install package on path.

    \b
    ```shell-session
    Examples:

    $ python dev.py build --asan ;
        ASAN_OPTIONS=detect_leaks=0:symbolize=1:strict_init_order=true
        LD_PRELOAD=$(gcc --print-file-name=libasan.so)
        python dev.py test -v -t
        ./scipy/ndimage/tests/test_morphology.py -- -s
    ```
    """
    ctx = CONTEXT

    werror = Option(
        ['--werror'], default=False, is_flag=True,
        help="Treat warnings as errors")
    gcov = Option(
        ['--gcov'], default=False, is_flag=True,
        help="enable C code coverage via gcov (requires GCC)."
             "gcov output goes to build/**/*.gc*")
    asan = Option(
        ['--asan'], default=False, is_flag=True,
        help=("Build and run with AddressSanitizer support. "
              "Note: the build system doesn't check whether "
              "the project is already compiled with ASan. "
              "If not, you need to do a clean build (delete "
              "build and build-install directories)."))
    debug = Option(
        ['--debug', '-d'], default=False, is_flag=True, help="Debug build")
    parallel = Option(
        ['--parallel', '-j'], default=None, metavar='N_JOBS',
        help=("Number of parallel jobs for building. "
              "This defaults to the number of available physical CPU cores"))
    setup_args = Option(
        ['--setup-args', '-C'], default=[], multiple=True,
        help=("Pass along one or more arguments to `meson setup` "
              "Repeat the `-C` in case of multiple arguments."))
    show_build_log = Option(
        ['--show-build-log'], default=False, is_flag=True,
        help="Show build output rather than using a log file")
    win_cp_openblas = Option(
        ['--win-cp-openblas'], default=False, is_flag=True,
        help=("If set, and on Windows, copy OpenBLAS lib to install directory "
              "after meson install. "
              "Note: this argument may be removed in the future once a "
              "`site.cfg`-like mechanism to select BLAS/LAPACK libraries is "
              "implemented for Meson"))

    @classmethod
    def setup_build(cls, dirs, args):
        """
        Setting up meson-build
        """
        for fn in PROJECT_ROOT_FILES:
            if not (dirs.root / fn).exists():
                print("To build the project, run dev.py in "
                      "git checkout or unpacked source")
                sys.exit(1)

        env = dict(os.environ)
        cmd = ["meson", "setup", dirs.build, "--prefix", dirs.installed]
        build_dir = dirs.build
        run_dir = Path()
        if build_dir.exists() and not (build_dir / 'meson-info').exists():
            if list(build_dir.iterdir()):
                raise RuntimeError("Can't build into non-empty directory "
                                   f"'{build_dir.absolute()}'")

        if sys.platform == "cygwin":
            # Cygwin only has netlib lapack, but can link against
            # OpenBLAS rather than netlib blas at runtime.  There is
            # no libopenblas-devel to enable linking against
            # openblas-specific functions or OpenBLAS Lapack
            cmd.extend(["-Dlapack=lapack", "-Dblas=blas"])

        build_options_file = (
            build_dir / "meson-info" / "intro-buildoptions.json")
        if build_options_file.exists():
            with open(build_options_file) as f:
                build_options = json.load(f)
            installdir = None
            for option in build_options:
                if option["name"] == "prefix":
                    installdir = option["value"]
                    break
            if installdir != str(dirs.installed):
                run_dir = build_dir
                cmd = ["meson", "setup", "--reconfigure",
                       "--prefix", str(dirs.installed)]
            else:
                return
        if args.werror:
            cmd += ["--werror"]
        if args.gcov:
            cmd += ['-Db_coverage=true']
        if args.asan:
            cmd += ['-Db_sanitize=address,undefined']
        if args.setup_args:
            cmd += [str(arg) for arg in args.setup_args]

        # Setting up meson build
        cmd_str = ' '.join([str(p) for p in cmd])
        cls.console.print(f"{EMOJI.cmd} [cmd] {cmd_str}")
        ret = subprocess.call(cmd, env=env, cwd=run_dir)
        if ret == 0:
            print("Meson build setup OK")
        else:
            print("Meson build setup failed!")
            sys.exit(1)
        return env

    @classmethod
    def build_project(cls, dirs, args, env):
        """
        Build a dev version of the project.
        """
        cmd = ["ninja", "-C", str(dirs.build)]
        if args.parallel is None:
            # Use number of physical cores rather than ninja's default of 2N+2,
            # to avoid out of memory issues (see gh-17941 and gh-18443)
            n_cores = cpu_count(only_physical_cores=True)
            cmd += [f"-j{n_cores}"]
        else:
            cmd += ["-j", str(args.parallel)]

        # Building with ninja-backend
        cmd_str = ' '.join([str(p) for p in cmd])
        cls.console.print(f"{EMOJI.cmd} [cmd] {cmd_str}")
        ret = subprocess.call(cmd, env=env, cwd=dirs.root)

        if ret == 0:
            print("Build OK")
        else:
            print("Build failed!")
            sys.exit(1)

    @classmethod
    def install_project(cls, dirs, args):
        """
        Installs the project after building.
        """
        if dirs.installed.exists():
            non_empty = len(os.listdir(dirs.installed))
            if non_empty and not dirs.site.exists():
                raise RuntimeError("Can't install in non-empty directory: "
                                   f"'{dirs.installed}'")
        cmd = ["meson", "install", "-C", args.build_dir, "--only-changed"]
        log_filename = dirs.root / 'meson-install.log'
        start_time = datetime.datetime.now()
        cmd_str = ' '.join([str(p) for p in cmd])
        cls.console.print(f"{EMOJI.cmd} [cmd] {cmd_str}")
        if args.show_build_log:
            ret = subprocess.call(cmd, cwd=dirs.root)
        else:
            print("Installing, see meson-install.log...")
            with open(log_filename, 'w') as log:
                p = subprocess.Popen(cmd, stdout=log, stderr=log,
                                     cwd=dirs.root)

            try:
                # Wait for it to finish, and print something to indicate the
                # process is alive, but only if the log file has grown (to
                # allow continuous integration environments kill a hanging
                # process accurately if it produces no output)
                last_blip = time.time()
                last_log_size = os.stat(log_filename).st_size
                while p.poll() is None:
                    time.sleep(0.5)
                    if time.time() - last_blip > 60:
                        log_size = os.stat(log_filename).st_size
                        if log_size > last_log_size:
                            elapsed = datetime.datetime.now() - start_time
                            print("    ... installation in progress ({} "
                                  "elapsed)".format(elapsed))
                            last_blip = time.time()
                            last_log_size = log_size

                ret = p.wait()
            except:  # noqa: E722
                p.terminate()
                raise
        elapsed = datetime.datetime.now() - start_time

        if ret != 0:
            if not args.show_build_log:
                with open(log_filename) as f:
                    print(f.read())
            print(f"Installation failed! ({elapsed} elapsed)")
            sys.exit(1)

        # ignore everything in the install directory.
        with open(dirs.installed / ".gitignore", "w") as f:
            f.write("*")

        if sys.platform == "cygwin":
            rebase_cmd = ["/usr/bin/rebase", "--database", "--oblivious"]
            rebase_cmd.extend(Path(dirs.installed).glob("**/*.dll"))
            subprocess.check_call(rebase_cmd)

        print("Installation OK")
        return

    @classmethod
    def copy_openblas(cls, dirs):
        """
        Copies OpenBLAS DLL to the SciPy install dir, and also overwrites the
        default `_distributor_init.py` file with the one
        we use for wheels uploaded to PyPI so that DLL gets loaded.

        Assumes pkg-config is installed and aware of OpenBLAS.

        The "dirs" parameter is typically a "Dirs" object with the
        structure as the following, say, if dev.py is run from the
        folder "repo":

        dirs = Dirs(
            root=WindowsPath('C:/.../repo'),
            build=WindowsPath('C:/.../repo/build'),
            installed=WindowsPath('C:/.../repo/build-install'),
            site=WindowsPath('C:/.../repo/build-install/Lib/site-packages'
            )

        """
        # Get OpenBLAS lib path from pkg-config
        cmd = ['pkg-config', '--variable', 'libdir', 'openblas']
        result = subprocess.run(cmd, capture_output=True, text=True)
        # pkg-config does not return any meaningful error message if fails
        if result.returncode != 0:
            print('"pkg-config --variable libdir openblas" '
                  'command did not manage to find OpenBLAS '
                  'succesfully. Try running manually on the '
                  'command prompt for more information.')
            print("OpenBLAS copy failed!")
            sys.exit(result.returncode)

        # Skip the drive letter of the path -> /c to get Windows drive
        # to be appended correctly to avoid "C:\c\..." from stdout.
        openblas_lib_path = Path(result.stdout.strip()[2:]).resolve()
        if not openblas_lib_path.stem == 'lib':
            raise RuntimeError('"pkg-config --variable libdir openblas" '
                               'command did not return a path ending with'
                               ' "lib" folder. Instead it returned '
                               f'"{openblas_lib_path}"')

        # Look in bin subdirectory for OpenBLAS binaries.
        bin_path = openblas_lib_path.parent / 'bin'
        # Locate, make output .libs directory in Scipy install directory.
        scipy_path = dirs.site / 'scipy'
        libs_path = scipy_path / '.libs'
        libs_path.mkdir(exist_ok=True)
        # Copy DLL files from OpenBLAS install to scipy install .libs subdir.
        for dll_fn in bin_path.glob('*.dll'):
            out_fname = libs_path / dll_fn.name
            print(f'Copying {dll_fn} ----> {out_fname}')
            out_fname.write_bytes(dll_fn.read_bytes())

        # Write _distributor_init.py to scipy install dir;
        # this ensures the .libs file is on the DLL search path at run-time,
        # so OpenBLAS gets found
        openblas_support = import_module_from_path(
            'openblas_support',
            dirs.root / 'tools' / 'openblas_support.py'
        )
        openblas_support.make_init(scipy_path)
        print('OpenBLAS copied')

    @classmethod
    def run(cls, add_path=False, **kwargs):
        kwargs.update(cls.ctx.get(kwargs))
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)

        cls.console = Console(theme=console_theme)
        dirs = Dirs(args)
        if args.no_build:
            print("Skipping build")
        else:
            env = cls.setup_build(dirs, args)
            cls.build_project(dirs, args, env)
            cls.install_project(dirs, args)
            if args.win_cp_openblas and platform.system() == 'Windows':
                cls.copy_openblas(dirs)

        # add site to sys.path
        if add_path:
            dirs.add_sys_path()


@cli.cls_cmd('test')
class Test(Task):
    """:wrench: Run tests.

    \b
    ```python
    Examples:

    $ python dev.py test -s {SAMPLE_SUBMODULE}
    $ python dev.py test -t scipy.optimize.tests.test_minimize_constrained
    $ python dev.py test -s cluster -m full --durations 20
    $ python dev.py test -s stats -- --tb=line  # `--` passes next args to pytest
    $ python dev.py test -b numpy -b pytorch -s cluster
    ```
    """  # noqa: E501
    ctx = CONTEXT

    verbose = Option(
        ['--verbose', '-v'], default=False, is_flag=True,
        help="more verbosity")
    # removed doctests as currently not supported by _lib/_testutils.py
    # doctests = Option(['--doctests'], default=False)
    coverage = Option(
        ['--coverage', '-c'], default=False, is_flag=True,
        help=("report coverage of project code. "
              "HTML output goes under build/coverage"))
    durations = Option(
        ['--durations', '-d'], default=None, metavar="NUM_TESTS",
        help="Show timing for the given number of slowest tests"
    )
    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='MODULE_NAME',
        help="Submodule whose tests to run (cluster, constants, ...)")
    tests = Option(
        ['--tests', '-t'], default=None, multiple=True, metavar='TESTS',
        help='Specify tests to run')
    mode = Option(
        ['--mode', '-m'], default='fast', metavar='MODE', show_default=True,
        help=("'fast', 'full', or something that could be passed to "
              "`pytest -m` as a marker expression"))
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='N_JOBS',
        help="Number of parallel jobs for testing"
    )
    array_api_backend = Option(
        ['--array-api-backend', '-b'], default=None, metavar='ARRAY_BACKEND',
        multiple=True,
        help=(
            "Array API backend ('all', 'numpy', 'pytorch', 'cupy', 'numpy.array_api')."
        )
    )
    # Argument can't have `help=`; used to consume all of `-- arg1 arg2 arg3`
    pytest_args = Argument(
        ['pytest_args'], nargs=-1, metavar='PYTEST-ARGS', required=False
    )

    TASK_META = {
        'task_dep': ['build'],
    }

    @classmethod
    def scipy_tests(cls, args, pytest_args):
        dirs = Dirs(args)
        dirs.add_sys_path()
        print(f"SciPy from development installed path at: {dirs.site}")

        # FIXME: support pos-args with doit
        extra_argv = pytest_args[:] if pytest_args else []
        if extra_argv and extra_argv[0] == '--':
            extra_argv = extra_argv[1:]

        if args.coverage:
            dst_dir = dirs.root / args.build_dir / 'coverage'
            fn = dst_dir / 'coverage_html.js'
            if dst_dir.is_dir() and fn.is_file():
                shutil.rmtree(dst_dir)
            extra_argv += ['--cov-report=html:' + str(dst_dir)]
            shutil.copyfile(dirs.root / '.coveragerc',
                            dirs.site / '.coveragerc')

        if args.durations:
            extra_argv += ['--durations', args.durations]

        # convert options to test selection
        if args.submodule:
            tests = [PROJECT_MODULE + "." + args.submodule]
        elif args.tests:
            tests = args.tests
        else:
            tests = None

        if len(args.array_api_backend) != 0:
            os.environ['SCIPY_ARRAY_API'] = json.dumps(list(args.array_api_backend))

        runner, version, mod_path = get_test_runner(PROJECT_MODULE)
        # FIXME: changing CWD is not a good practice
        with working_dir(dirs.site):
            print("Running tests for {} version:{}, installed at:{}".format(
                        PROJECT_MODULE, version, mod_path))
            # runner verbosity - convert bool to int
            verbose = int(args.verbose) + 1
            result = runner(  # scipy._lib._testutils:PytestTester
                args.mode,
                verbose=verbose,
                extra_argv=extra_argv,
                doctests=False,
                coverage=args.coverage,
                tests=tests,
                parallel=args.parallel)
        return result

    @classmethod
    def run(cls, pytest_args, **kwargs):
        """run unit-tests"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        return cls.scipy_tests(args, pytest_args)


@cli.cls_cmd('bench')
class Bench(Task):
    """:wrench: Run benchmarks.

    \b
    ```python
     Examples:

    $ python dev.py bench -t integrate.SolveBVP
    $ python dev.py bench -t linalg.Norm
    $ python dev.py bench --compare main
    ```
    """
    ctx = CONTEXT
    TASK_META = {
        'task_dep': ['build'],
    }
    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")
    tests = Option(
        ['--tests', '-t'], default=None, multiple=True,
        metavar='TESTS', help='Specify tests to run')
    compare = Option(
        ['--compare', '-c'], default=None, metavar='COMPARE', multiple=True,
        help=(
            "Compare benchmark results of current HEAD to BEFORE. "
            "Use an additional --bench COMMIT to override HEAD with COMMIT. "
            "Note that you need to commit your changes first!"))

    @staticmethod
    def run_asv(dirs, cmd):
        EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
                      '/usr/local/lib/ccache', '/usr/local/lib/f90cache']
        bench_dir = dirs.root / 'benchmarks'
        sys.path.insert(0, str(bench_dir))
        # Always use ccache, if installed
        env = dict(os.environ)
        env['PATH'] = os.pathsep.join(EXTRA_PATH +
                                      env.get('PATH', '').split(os.pathsep))
        # Control BLAS/LAPACK threads
        env['OPENBLAS_NUM_THREADS'] = '1'
        env['MKL_NUM_THREADS'] = '1'

        # Limit memory usage
        from benchmarks.common import set_mem_rlimit
        try:
            set_mem_rlimit()
        except (ImportError, RuntimeError):
            pass
        try:
            return subprocess.call(cmd, env=env, cwd=bench_dir)
        except OSError as err:
            if err.errno == errno.ENOENT:
                cmd_str = " ".join(cmd)
                print(f"Error when running '{cmd_str}': {err}\n")
                print("You need to install Airspeed Velocity "
                      "(https://airspeed-velocity.github.io/asv/)")
                print("to run Scipy benchmarks")
                return 1
            raise

    @classmethod
    def scipy_bench(cls, args):
        dirs = Dirs(args)
        dirs.add_sys_path()
        print(f"SciPy from development installed path at: {dirs.site}")
        with working_dir(dirs.site):
            runner, version, mod_path = get_test_runner(PROJECT_MODULE)
            extra_argv = []
            if args.tests:
                extra_argv.append(args.tests)
            if args.submodule:
                extra_argv.append([args.submodule])

            bench_args = []
            for a in extra_argv:
                bench_args.extend(['--bench', ' '.join(str(x) for x in a)])
            if not args.compare:
                print("Running benchmarks for Scipy version %s at %s"
                      % (version, mod_path))
                cmd = ['asv', 'run', '--dry-run', '--show-stderr',
                       '--python=same', '--quick'] + bench_args
                retval = cls.run_asv(dirs, cmd)
                sys.exit(retval)
            else:
                if len(args.compare) == 1:
                    commit_a = args.compare[0]
                    commit_b = 'HEAD'
                elif len(args.compare) == 2:
                    commit_a, commit_b = args.compare
                else:
                    print("Too many commits to compare benchmarks for")
                # Check for uncommitted files
                if commit_b == 'HEAD':
                    r1 = subprocess.call(['git', 'diff-index', '--quiet',
                                          '--cached', 'HEAD'])
                    r2 = subprocess.call(['git', 'diff-files', '--quiet'])
                    if r1 != 0 or r2 != 0:
                        print("*" * 80)
                        print("WARNING: you have uncommitted changes --- "
                              "these will NOT be benchmarked!")
                        print("*" * 80)

                # Fix commit ids (HEAD is local to current repo)
                p = subprocess.Popen(['git', 'rev-parse', commit_b],
                                     stdout=subprocess.PIPE)
                out, err = p.communicate()
                commit_b = out.strip()

                p = subprocess.Popen(['git', 'rev-parse', commit_a],
                                     stdout=subprocess.PIPE)
                out, err = p.communicate()
                commit_a = out.strip()
                cmd_compare = [
                    'asv', 'continuous', '--show-stderr', '--factor', '1.05',
                    '--quick', commit_a, commit_b
                ] + bench_args
                cls.run_asv(dirs, cmd_compare)
                sys.exit(1)

    @classmethod
    def run(cls, **kwargs):
        """run benchmark"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        cls.scipy_bench(args)


###################
# linters

def emit_cmdstr(cmd):
    """Print the command that's being run to stdout

    Note: cannot use this in the below tasks (yet), because as is these command
    strings are always echoed to the console, even if the command isn't run
    (but for example the `build` command is run).
    """
    console = Console(theme=console_theme)
    # The [cmd] square brackets controls the font styling, typically in italics
    # to differentiate it from other stdout content
    console.print(f"{EMOJI.cmd} [cmd] {cmd}")


def task_lint():
    # Lint just the diff since branching off of main using a
    # stricter configuration.
    # emit_cmdstr(os.path.join('tools', 'lint.py') + ' --diff-against main')
    return {
        'basename': 'lint',
        'actions': [str(Dirs().root / 'tools' / 'lint.py') +
                    ' --diff-against=main'],
        'doc': 'Lint only files modified since last commit (stricter rules)',
    }


def task_unicode_check():
    # emit_cmdstr(os.path.join('tools', 'unicode-check.py'))
    return {
        'basename': 'unicode-check',
        'actions': [str(Dirs().root / 'tools' / 'unicode-check.py')],
        'doc': 'Check for disallowed Unicode characters in the SciPy Python '
               'and Cython source code.',
    }


def task_check_test_name():
    # emit_cmdstr(os.path.join('tools', 'check_test_name.py'))
    return {
        "basename": "check-testname",
        "actions": [str(Dirs().root / "tools" / "check_test_name.py")],
        "doc": "Check tests are correctly named so that pytest runs them."
    }


@cli.cls_cmd('lint')
class Lint():
    """:dash: Run linter on modified files and check for
    disallowed Unicode characters and possibly-invalid test names."""
    def run():
        run_doit_task({
            'lint': {},
            'unicode-check': {},
            'check-testname': {},
        })


@cli.cls_cmd('mypy')
class Mypy(Task):
    """:wrench: Run mypy on the codebase."""
    ctx = CONTEXT

    TASK_META = {
        'task_dep': ['build'],
    }

    @classmethod
    def run(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        dirs = Dirs(args)

        try:
            import mypy.api
        except ImportError as e:
            raise RuntimeError(
                "Mypy not found. Please install it by running "
                "pip install -r mypy_requirements.txt from the repo root"
            ) from e

        config = dirs.root / "mypy.ini"
        check_path = PROJECT_MODULE

        with working_dir(dirs.site):
            # By default mypy won't color the output since it isn't being
            # invoked from a tty.
            os.environ['MYPY_FORCE_COLOR'] = '1'
            # Change to the site directory to make sure mypy doesn't pick
            # up any type stubs in the source tree.
            emit_cmdstr(f"mypy.api.run --config-file {config} {check_path}")
            report, errors, status = mypy.api.run([
                "--config-file",
                str(config),
                check_path,
            ])
        print(report, end='')
        print(errors, end='', file=sys.stderr)
        return status == 0


##########################################
# DOC

@cli.cls_cmd('doc')
class Doc(Task):
    """:wrench: Build documentation.

    TARGETS: Sphinx build targets [default: 'html']

    Running `python dev.py doc -j8 html` is equivalent to:
    1. Execute build command (skip by passing the global `-n` option).
    2. Set the PYTHONPATH environment variable
       (query with `python dev.py -n show_PYTHONPATH`).
    3. Run make on `doc/Makefile`, i.e.: `make -C doc -j8 TARGETS`

    To remove all generated documentation do: `python dev.py -n doc clean`
    """
    ctx = CONTEXT

    args = Argument(['args'], nargs=-1, metavar='TARGETS', required=False)
    list_targets = Option(
        ['--list-targets', '-t'], default=False, is_flag=True,
        help='List doc targets',
    )
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='N_JOBS',
        help="Number of parallel jobs"
    )
    no_cache = Option(
        ['--no-cache'], default=False, is_flag=True,
        help="Forces a full rebuild of the docs. Note that this may be " + \
             "needed in order to make docstring changes in C/Cython files " + \
             "show up."
    )

    @classmethod
    def task_meta(cls, list_targets, parallel, no_cache, args, **kwargs):
        if list_targets:  # list MAKE targets, remove default target
            task_dep = []
            targets = ''
        else:
            task_dep = ['build']
            targets = ' '.join(args) if args else 'html'

        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        build_args = Args(**kwargs)
        dirs = Dirs(build_args)

        make_params = [f'PYTHON="{sys.executable}"']
        if parallel or no_cache:
            sphinxopts = ""
            if parallel:
                sphinxopts += f"-j{parallel} "
            if no_cache:
                sphinxopts += "-E"
            make_params.append(f'SPHINXOPTS="{sphinxopts}"')

        return {
            'actions': [
                # move to doc/ so local scipy does not get imported
                (f'cd doc; env PYTHONPATH="{dirs.site}" '
                 f'make {" ".join(make_params)} {targets}'),
            ],
            'task_dep': task_dep,
            'io': {'capture': False},
        }


@cli.cls_cmd('refguide-check')
class RefguideCheck(Task):
    """:wrench: Run refguide check."""
    ctx = CONTEXT

    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")
    verbose = Option(
        ['--verbose', '-v'], default=False, is_flag=True, help="verbosity")

    @classmethod
    def task_meta(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        dirs = Dirs(args)

        cmd = [f'{sys.executable}',
               str(dirs.root / 'tools' / 'refguide_check.py'),
               '--doctests']
        if args.verbose:
            cmd += ['-vvv']
        if args.submodule:
            cmd += [args.submodule]
        cmd_str = ' '.join(cmd)
        return {
            'actions': [f'env PYTHONPATH={dirs.site} {cmd_str}'],
            'task_dep': ['build'],
            'io': {'capture': False},
        }


##########################################
# ENVS

@cli.cls_cmd('python')
class Python:
    """:wrench: Start a Python shell with PYTHONPATH set.

    ARGS: Arguments passed to the Python interpreter.
          If not set, an interactive shell is launched.

    Running `python dev.py shell my_script.py` is equivalent to:
    1. Execute build command (skip by passing the global `-n` option).
    2. Set the PYTHONPATH environment variable
       (query with `python dev.py -n show_PYTHONPATH`).
    3. Run interpreter: `python my_script.py`
    """
    ctx = CONTEXT
    pythonpath = Option(
        ['--pythonpath', '-p'], metavar='PYTHONPATH', default=None,
        help='Paths to prepend to PYTHONPATH')
    extra_argv = Argument(
        ['extra_argv'], nargs=-1, metavar='ARGS', required=False)

    @classmethod
    def _setup(cls, pythonpath, **kwargs):
        vals = Build.opt_defaults()
        vals.update(kwargs)
        Build.run(add_path=True, **vals)
        if pythonpath:
            for p in reversed(pythonpath.split(os.pathsep)):
                sys.path.insert(0, p)

    @classmethod
    def run(cls, pythonpath, extra_argv=None, **kwargs):
        cls._setup(pythonpath, **kwargs)
        if extra_argv:
            # Don't use subprocess, since we don't want to include the
            # current path in PYTHONPATH.
            sys.argv = extra_argv
            with open(extra_argv[0]) as f:
                script = f.read()
            sys.modules['__main__'] = new_module('__main__')
            ns = dict(__name__='__main__', __file__=extra_argv[0])
            exec(script, ns)
        else:
            import code
            code.interact()


@cli.cls_cmd('ipython')
class Ipython(Python):
    """:wrench: Start IPython shell with PYTHONPATH set.

    Running `python dev.py ipython` is equivalent to:
    1. Execute build command (skip by passing the global `-n` option).
    2. Set the PYTHONPATH environment variable
       (query with `python dev.py -n show_PYTHONPATH`).
    3. Run the `ipython` interpreter.
    """
    ctx = CONTEXT
    pythonpath = Python.pythonpath

    @classmethod
    def run(cls, pythonpath, **kwargs):
        cls._setup(pythonpath, **kwargs)
        import IPython
        IPython.embed(user_ns={})


@cli.cls_cmd('shell')
class Shell(Python):
    """:wrench: Start Unix shell with PYTHONPATH set.

    Running `python dev.py shell` is equivalent to:
    1. Execute build command (skip by passing the global `-n` option).
    2. Open a new shell.
    3. Set the PYTHONPATH environment variable in shell
       (query with `python dev.py -n show_PYTHONPATH`).
    """
    ctx = CONTEXT
    pythonpath = Python.pythonpath
    extra_argv = Python.extra_argv

    @classmethod
    def run(cls, pythonpath, extra_argv, **kwargs):
        cls._setup(pythonpath, **kwargs)
        shell = os.environ.get('SHELL', 'sh')
        click.echo(f"Spawning a Unix shell '{shell}' ...")
        os.execv(shell, [shell] + list(extra_argv))
        sys.exit(1)


@cli.cls_cmd('show_PYTHONPATH')
class ShowDirs(Python):
    """:information: Show value of the PYTHONPATH environment variable used in
    this script.

    PYTHONPATH sets the default search path for module files for the
    interpreter. Here, it includes the path to the local SciPy build
    (typically `.../build-install/lib/python3.10/site-packages`).

    Use the global option `-n` to skip the building step, e.g.:
    `python dev.py -n show_PYTHONPATH`
    """
    ctx = CONTEXT
    pythonpath = Python.pythonpath
    extra_argv = Python.extra_argv

    @classmethod
    def run(cls, pythonpath, extra_argv, **kwargs):
        cls._setup(pythonpath, **kwargs)
        py_path = os.environ.get('PYTHONPATH', '')
        click.echo(f"PYTHONPATH={py_path}")


@cli.command()
@click.argument('version_args', nargs=2)
@click.pass_obj
def notes(ctx_obj, version_args):
    """:ledger: Release notes and log generation.

    \b
    ```python
     Example:

    $ python dev.py notes v1.7.0 v1.8.0
    ```
    """
    if version_args:
        sys.argv = version_args
        log_start = sys.argv[0]
        log_end = sys.argv[1]
    cmd = f"python tools/write_release_and_log.py {log_start} {log_end}"
    click.echo(cmd)
    try:
        subprocess.run([cmd], check=True, shell=True)
    except subprocess.CalledProcessError:
        print('Error caught: Incorrect log start or log end version')


@cli.command()
@click.argument('revision_args', nargs=2)
@click.pass_obj
def authors(ctx_obj, revision_args):
    """:ledger: Generate list of authors who contributed within revision
    interval.

    \b
    ```python
    Example:

    $ python dev.py authors v1.7.0 v1.8.0
    ```
    """
    if revision_args:
        sys.argv = revision_args
        start_revision = sys.argv[0]
        end_revision = sys.argv[1]
    cmd = f"python tools/authors.py {start_revision}..{end_revision}"
    click.echo(cmd)
    try:
        subprocess.run([cmd], check=True, shell=True)
    except subprocess.CalledProcessError:
        print('Error caught: Incorrect revision start or revision end')


# The following CPU core count functions were taken from loky/backend/context.py
# See https://github.com/joblib/loky

# Cache for the number of physical cores to avoid repeating subprocess calls.
# It should not change during the lifetime of the program.
physical_cores_cache = None


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
            "the number of cores you want to use."
        )
        traceback.print_tb(exception.__traceback__)

    return aggregate_cpu_count


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
            # cgroup CPU bandwith limits
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
                "Please install psutil or explictly set LOKY_MAX_CPU_COUNT."
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


if __name__ == '__main__':
    cli()
