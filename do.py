#! /usr/bin/env python3

'''
Developer CLI: building (meson), tests, benchmark, etc.

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authring the development tasks.

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
import shutil
import json
import datetime
import time
import platform
import importlib.util
import errno
import contextlib
from sysconfig import get_path
from pathlib import Path
from collections import namedtuple
from types import ModuleType as new_module
from dataclasses import dataclass

import click
from click import Option, Argument
from doit import task_params
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


class EMOJI:
    cmd = ":computer:"


rich_click.STYLE_ERRORS_SUGGESTION = "yellow italic"
rich_click.SHOW_ARGUMENTS = True
rich_click.GROUP_ARGUMENTS_OPTIONS = False
rich_click.SHOW_METAVARS_COLUMN = True
rich_click.USE_MARKDOWN = True
rich_click.OPTION_GROUPS = {
    "do.py": [
        {
            "name": "Options",
            "options": [
                "--help", "--build-dir", "--no-build", "--install-prefix"],
        },
    ],

    "do.py test": [
        {
            "name": "Options",
            "options": ["--help", "--verbose", "--parallel", "--coverage"],
        },
        {
            "name": "Options: test selection",
            "options": ["--submodule", "--tests", "--mode"],
        },
    ],
}
rich_click.COMMAND_GROUPS = {
    "do.py": [
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
            "commands": ["shell", "python", "ipython"],
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
        help=(":wrench: do not build the project"
              " (note event python only modification require build)")),
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

    Commands that require a built/installed instance are marked with :wrench:.



    **python do.py --build-dir my-build test -s stats**
    """
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
        # relative path for site-package with py version
        # i.e. 'lib/python3.10/site-packages'
        py_lib_path = Path(get_path('platlib')).relative_to(sys.exec_prefix)
        self.site = self.installed / py_lib_path

    def add_sys_path(self):
        """Add site dir to sys.path / PYTHONPATH"""
        site_dir = str(self.site)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))


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
    """:wrench: build & install package on path"""
    ctx = CONTEXT

    werror = Option(
        ['--werror'], default=False, is_flag=True,
        help="Treat warnings as errors")
    gcov = Option(
        ['--gcov'], default=False, is_flag=True,
        help="enable C code coverage via gcov (requires GCC)."
             "gcov output goes to build/**/*.gc*")
    debug = Option(
        ['--debug', '-d'], default=False, is_flag=True, help="Debug build")
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='PARALLEL',
        help="Number of parallel jobs for build and testing")
    show_build_log = Option(
        ['--show-build-log'], default=False, is_flag=True,
        help="Show build output rather than using a log file")
    win_cp_openblas = Option(
        ['--win-cp-openblas'], default=False, is_flag=True,
        help=("If set, and on Windows, copy OpenBLAS lib to install directory"
              "after meson install. "
              "Note: this argument may be removed in the future once a "
              "`site.cfg`-like mechanism to select BLAS/LAPACK libraries is"
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
        if build_dir.exists():
            build_options_file = (
                build_dir / "meson-info" / "intro-buildoptions.json")
            with open(build_options_file) as f:
                build_options = json.load(f)
            installdir = None
            for option in build_options:
                if option["name"] == "prefix":
                    installdir = option["value"]
                    break
            if installdir != str(dirs.installed):
                run_dir = build_dir
                cmd = ["meson", "--reconfigure",
                       "--prefix", str(dirs.installed)]
            else:
                return
        if args.werror:
            cmd += ["--werror"]
        if args.gcov:
            cmd += ['-Db_coverage=true']
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
        if args.parallel > 1:
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
        cmd = ["meson", "install", "-C", args.build_dir]
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
                            print("    ... installation in progress ({0} "
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
                with open(log_filename, 'r') as f:
                    print(f.read())
            print(f"Installation failed! ({elapsed} elapsed)")
            sys.exit(1)

        # ignore everything in the install directory.
        with open(dirs.installed / ".gitignore", "w") as f:
            f.write("*")

        print("Installation OK")
        return

    @classmethod
    def copy_openblas(cls, dirs):
        """
        Copies OpenBLAS DLL to the SciPy install dir, and also overwrites the
        default `_distributor_init.py` file with the one
        we use for wheels uploaded to PyPI so that DLL gets loaded.

        Assumes pkg-config is installed and aware of OpenBLAS.
        """
        # Get OpenBLAS lib path from pkg-config
        cmd = ['pkg-config', '--variable', 'libdir', 'openblas']
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stderrr)
            return result.returncode

        openblas_lib_path = Path(result.stdout.strip())
        if not openblas_lib_path.stem == 'lib':
            raise RuntimeError(
                f'Expecting "lib" at end of "{openblas_lib_path}"')

        # Look in bin subdirectory for OpenBLAS binaries.
        bin_path = openblas_lib_path.parent / 'bin'
        # Locate, make output .libs directory in Scipy install directory.
        scipy_path = dirs.site / 'scipy'
        libs_path = scipy_path / '.libs'
        libs_path.mkdir(exist_ok=True)
        # Copy DLL files from OpenBLAS install to scipy install .libs subdir.
        for dll_fn in bin_path.glob('*.dll'):
            out_fname = libs_path / dll_fn.parts[-1]
            print(f'Copying {dll_fn} to {out_fname}')
            out_fname.write_bytes(dll_fn.read_bytes())

        # Write _distributor_init.py to scipy install dir;
        # this ensures the .libs file is on the DLL search path at run-time,
        # so OpenBLAS gets found
        openblas_support = import_module_from_path(
            'openblas_support',
            dirs.root / 'tools' / 'openblas_support.py')
        openblas_support.make_init(scipy_path)
        return 0

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
                if cls.copy_openblas(dirs) == 0:
                    print('OpenBLAS copied')
                else:
                    print("OpenBLAS copy failed!")
                    sys.exit(1)

        # add site to sys.path
        if add_path:
            dirs.add_sys_path()


@cli.cls_cmd('test')
class Test(Task):
    """:wrench: Run tests

    Examples:

    $ python do.py test -s {SAMPLE_SUBMODULE}
    $ python do.py test -t scipy.optimize.tests.test_minimize_constrained
    $ python do.py test -s stats -- --tb=line
    """
    ctx = CONTEXT

    verbose = Option(
        ['--verbose', '-v'], default=False, is_flag=True,
        help="more verbosity")
    # removed doctests as currently not supported by _lib/_testutils.py
    # doctests = Option(['--doctests'], default=False)
    coverage = Option(
        ['--coverage'], default=False, is_flag=True,
        help=("report coverage of project code. "
              "HTML output goes under build/coverage"))
    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")
    tests = Option(
        ['--tests', '-t'], default=None, multiple=True, metavar='TESTS',
        help='Specify tests to run')
    mode = Option(
        ['--mode', '-m'], default='fast', metavar='MODE', show_default=True,
        help=("'fast', 'full', or something that could be passed to "
              "`pytest -m` as a marker expression"))
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='PARALLEL',
        help="Number of parallel jobs for testing"
    )
    pytest_args = Argument(
        ['pytest_args'], nargs=-1, metavar='PYTEST-ARGS', required=False)

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

        # convert options to test selection
        if args.submodule:
            tests = [PROJECT_MODULE + "." + args.submodule]
        elif args.tests:
            tests = args.tests
        else:
            tests = None

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
    """:wrench: Run benchmarks

     Examples:

    $ python do.py bench -t integrate.SolveBVP
    $ python do.py bench -t linalg.Norm
    $ python do.py bench --compare main

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
                    commit_a, commit_b
                ] + bench_args
                cls.run_asv(dirs, cmd_compare)
                sys.exit(1)

    @classmethod
    def run(cls, **kwargs):
        """run benchamark"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        cls.scipy_bench(args)


###################
# linters

@task_params([{'name': 'output_file', 'long': 'output-file', 'default': None,
               'help': 'Redirect report to a file'}])
def task_flake8(output_file):
    """Run flake8 over the code base and benchmarks."""
    opts = ''
    if output_file:
        opts += f'--output-file={output_file}'
    return {
        'actions': [f"flake8 {opts} scipy benchmarks/benchmarks"],
        'doc': 'Lint scipy and benchmarks directory',
    }


def task_pep8diff():
    # Lint just the diff since branching off of main using a
    # stricter configuration.
    return {
        'basename': 'pep8-diff',
        'actions': [str(Dirs().root / 'tools' / 'lint_diff.py')],
        'doc': 'Lint only files modified since last commit (stricker rules)',
    }


@cli.cls_cmd('lint')
class Lint():
    """:dash: run flake8, and check PEP 8 compliance on branch diff."""
    output_file = Option(
        ['--output-file'], default=None, help='Redirect report to a file')

    def run(output_file):
        opts = {'output_file': output_file}
        run_doit_task({'flake8': opts, 'pep8-diff': {}})


@cli.cls_cmd('mypy')
class Mypy(Task):
    """:wrench: Run mypy on the codebase"""
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
    """:wrench: Build documentation

TARGETS: Sphinx build targets [default: 'html-scipyorg']
"""
    ctx = CONTEXT

    args = Argument(['args'], nargs=-1, metavar='TARGETS', required=False)
    list_targets = Option(
        ['--list-targets', '-t'], default=False, is_flag=True,
        help='List doc targets',
    )
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='PARALLEL',
        help="Number of parallel jobs"
    )

    @classmethod
    def task_meta(cls, list_targets, parallel, args, **kwargs):
        if list_targets:  # list MAKE targets, remove default target
            task_dep = []
            targets = ''
        else:
            task_dep = ['build']
            targets = ' '.join(args) if args else 'html-scipyorg'

        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        build_args = Args(**kwargs)
        dirs = Dirs(build_args)

        make_params = [f'PYTHON="{sys.executable}"']
        if parallel:
            make_params.append(f'SPHINXOPTS="-j{parallel}"')

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
    """:wrench: Run refguide check"""
    ctx = CONTEXT

    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")

    @classmethod
    def task_meta(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        dirs = Dirs(args)

        cmd = [str(dirs.root / 'tools' / 'refguide_check.py'), '--doctests']
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
class Python():
    """:wrench: Start a Python shell with PYTHONPATH set"""
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
            with open(extra_argv[0], 'r') as f:
                script = f.read()
            sys.modules['__main__'] = new_module('__main__')
            ns = dict(__name__='__main__', __file__=extra_argv[0])
            exec(script, ns)
        else:
            import code
            code.interact()


@cli.cls_cmd('ipython')
class Ipython(Python):
    """:wrench: Start IPython shell with PYTHONPATH set"""
    ctx = CONTEXT
    pythonpath = Python.pythonpath

    @classmethod
    def run(cls, pythonpath, **kwargs):
        cls._setup(pythonpath, **kwargs)
        import IPython
        IPython.embed(user_ns={})


@cli.cls_cmd('shell')
class Shell(Python):
    """:wrench: Start Unix shell with PYTHONPATH set"""
    ctx = CONTEXT
    pythonpath = Python.pythonpath
    extra_argv = Python.extra_argv

    @classmethod
    def run(cls, pythonpath, extra_argv, **kwargs):
        cls._setup(pythonpath, **kwargs)
        shell = os.environ.get('SHELL', 'sh')
        print("Spawning a Unix shell...")
        os.execv(shell, [shell] + list(extra_argv))
        sys.exit(1)


@cli.command()
@click.argument('version_args', nargs=2)
@click.pass_obj
def notes(ctx_obj, version_args):
    """:ledger: Release notes and log generation

     Example:

    $ python do.py notes v1.7.0 v1.8.0
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
    """:ledger: Generate list of authors who contributed within revision interval

    Example:

    $ python do.py authors v1.7.0 v1.8.0
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


if __name__ == '__main__':
    cli()
