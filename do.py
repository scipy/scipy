#! /usr/bin/env python3

'''
Developer CLI: building (meson), tests, benchmark, etc.

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authring the development tasks.

Note this requires the unreleased doit 0.35.
And also PyPI packages: click, rich, rich-click(1.3.0)


# USAGE:

## 1 - click API

Commands cab added using default Click API. i.e.

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

    verbose = Option(['--verbose', '-v'], default=False, is_flag=True, help="more verbosity")

    @classmethod
    def run(cls, **kwargs): # kwargs contains options from class and CONTEXT
        print('Running tests...')
```

## 3 - This class based interface can be run as a doit task simply by subclassing from Task

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



# TODO:
First milestone replace dev.py

- [ ] support --no-build/-n
- [ ] click does not support '--'. used by test command
      https://github.com/pallets/click/issues/1340

commands:
- [ ] bench
- [ ] --bench-compare

cleanup + doit:
- [ ] move out non-scipy code
- [ ] copy used code from dev.py
- [ ] doit reporter when running under click
- [ ] doit: add io task attribute to control capture/caching of stdout/stderr


BUG:
- [ ] python dev.py -t scipy.optimize.tests.test_minimize_constrained.py::TestTrustRegionConstr::test_args
      unknown marker "slow". seems conftest is not being loaded
- [ ] doc: fail on python3.10 - distutils is deprecated
'''

import os
import subprocess
import sys
import re
import shutil
from sysconfig import get_path
from pathlib import Path
from collections import namedtuple
from types import ModuleType as new_module

import click
from click import Option
from click.globals import get_current_context
from rich_click import RichCommand, RichGroup
from doit import task_params
from doit.task import Task as DoitTask
from doit.cmd_base import ModuleTaskLoader, get_loader
from doit.doit_cmd import DoitMain
from doit.cmdparse import CmdParse, CmdParseError
from doit.exceptions import InvalidDodoFile, InvalidCommand, InvalidTask

# keep compatibility and re-use dev.py code
import dev as dev_module

DOIT_CONFIG = {
    'verbosity': 2,
}


##########################################
### doit / click integration through custom class interface

class Context():
    """Higher level to allow some level of Click.context with doit"""
    def __init__(self, options: dict):
        self.options = options
        self.vals = {}

    def get(self, save=None):
        # try to get from Click
        ctx = get_current_context(silent=True)
        if ctx:
            return ctx.obj
        else:
            if save:
                for name in self.options.keys():
                    if name in save:
                        self.vals[name] = save[name]
            return self.vals



param_type_map = {
    'text': str,
    'boolean': bool,
    'integer': int,
}
def param_click2doit(name: str, val: Option):
    """Convert click param to dict used by doit.cmdparse"""
    pd = {
        'name': name,
        'type': param_type_map[val.type.name], # FIXME: add all types
        'default': val.default,
        'help': val.help or '',
        'metavar': val.metavar,
    }
    for opt in val.opts:
        if opt[:2] == '--':
            pd['long'] = opt[2:]
        elif opt[0] == '-':
            pd['short'] = opt[1:]
    return pd



# convert click.types.ParamType.name to doit param type
CAMEL_PATTERN = re.compile(r'(?!^)([A-Z]+)')
def to_camel(name):
    return CAMEL_PATTERN.sub(r'-\1', name).lower()

def run_as_py_action(cls):
    """used by doit loader to create task instances"""
    if cls is Task:
        return
    task_kwargs = getattr(cls, 'TASK_META', {})
    return DoitTask(
        # convert name to kebab-case
        name=to_camel(cls.__name__),
        doc=cls.__doc__,
        actions=[cls.run],
        params=cls._params,
        **task_kwargs,
    )


class MetaclassDoitTask(type):
    def __new__(meta_cls, name, bases, dct):
        # params/opts from Context and Option attributes
        cls = super().__new__(meta_cls, name, bases, dct)
        params = []
        if ctx := getattr(cls, 'ctx', None):
            for ctx_opt in ctx.options.values():
                params.append(param_click2doit(ctx_opt.name, ctx_opt))
        for attr_name, val in cls.__dict__.items():
            if isinstance(val, Option):
                params.append(param_click2doit(attr_name, val))
        cls._params = params

        task_basename = to_camel(cls.__name__)
        if hasattr(cls, 'task_meta'):
            def creator(**kwargs):
                task_meta = cls.task_meta(**kwargs)
                if 'basename' not in task_meta:
                    task_meta['basename'] = task_basename
                if 'doc' not in task_meta:
                    task_meta['doc'] = cls.__doc__
                return task_meta
            creator._task_creator_params = cls._params
        else:
            def creator():
                return run_as_py_action(cls)
        creator.basename = task_basename
        cls.create_doit_tasks = creator
        return cls


class Task(metaclass=MetaclassDoitTask):
    """Base class to define doit task and/or click command"""

    @classmethod
    def opt_defaults(cls):
        """helper used by another click commands to call this command"""
        return {p['name']:p['default'] for p in cls._params}



class CliGroup(RichGroup):
    COMMAND_CLASS = RichCommand

    def cls_cmd(self, name):
        """class decorator, convert to click.Command"""
        def register_click(cls):
            # get options for class definition
            opts = []
            for attr_name, attr_val in cls.__dict__.items():
                if isinstance(attr_val, Option):
                    opts.append(attr_val)

            if issubclass(cls, Task):
                # run as doit task
                def callback(**kwargs):
                    run_doit_task({name: kwargs})
            else:
                # run as plain function
                def callback(**kwargs):
                    cls.run(**kwargs)

            click_cmd = RichCommand(
                name=name,
                callback=callback,
                help=cls.__doc__,
                params=opts,
            )
            self.add_command(click_cmd)
            return cls
        return register_click




class DoitMainAPI(DoitMain):
    """add new method to run tasks with parsed command line"""

    def run_tasks(self, tasks):
        """
        :params task: str - task name
        """
        # get list of available commands
        sub_cmds = self.get_cmds()
        task_loader = get_loader(self.config, self.task_loader, sub_cmds)

        # execute command
        cmd_name = 'run'
        command = sub_cmds.get_plugin(cmd_name)(
            task_loader=task_loader,
            config=self.config,
            bin_name=self.BIN_NAME,
            cmds=sub_cmds,
            opt_vals={},
            )

        try:
            cmd_opt = CmdParse(command.get_options())
            params, _ = cmd_opt.parse([])
            args = tasks
            command.execute(params, args)
        except (CmdParseError, InvalidDodoFile,
                InvalidCommand, InvalidTask) as err:
            if isinstance(err, InvalidCommand):
                err.cmd_used = cmd_name
                err.bin_name = self.BIN_NAME
            raise err


def run_doit_task(tasks):
    """
      :param tasks: (dict) task_name -> {options}
    """
    loader = ModuleTaskLoader(globals())
    loader.task_opts = tasks # task_opts will be used as @task_param
    doit_main = DoitMainAPI(loader, extra_config={'GLOBAL': {'verbosity': 2}})
    return doit_main.run_tasks(list(tasks.keys()))


###########################################
### SciPY

from rich_click import rich_click

rich_click.STYLE_ERRORS_SUGGESTION = "yellow italic"
rich_click.SHOW_ARGUMENTS = True
rich_click.GROUP_ARGUMENTS_OPTIONS = False
rich_click.SHOW_METAVARS_COLUMN = True
rich_click.USE_MARKDOWN = True
rich_click.OPTION_GROUPS = {
    "do.py": [
        {
            "name": "Options",
            "options": ["--help", "--build-dir", "--install-prefix"],
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
            "commands": ["pep8", "mypy"],
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
    ]
}


CONTEXT = Context({
    'build_dir': Option(
        ['--build-dir'], default='build', metavar='BUILD_DIR', show_default=True,
        help=':wrench: Relative path to the build directory.'),
    'install_prefix': Option(
        ['--install-prefix'], default=None, metavar='INSTALL_DIR',
        help=":wrench: Relative path to the install directory. Default is <build-dir>-install."),
})



@click.group(cls=CliGroup)
@click.pass_context
def cli(ctx, **kwargs):
    """Developer Tool for SciPy

    Commands that require a built/installed instance are marked with :wrench:.



    **python do.py --build-dir my-build test -s stats**
    """
    ctx.ensure_object(dict)
    for opt_name in CONTEXT.options.keys():
        ctx.obj[opt_name] = kwargs.get(opt_name)
cli.params.extend(CONTEXT.options.values())


def get_dirs(args):
    """return (install_dir, site_dir)"""
    build_dir = Path(args.build_dir)
    install_dir = args.install_prefix
    if not install_dir:
        install_dir = build_dir.parent / (build_dir.stem + "-install")

    py_lib_path = Path(get_path('platlib')).relative_to(sys.exec_prefix)
    site_dir = str(Path(install_dir) / py_lib_path)
    return install_dir, site_dir



# FIXME: remove
def set_installed(args):
    """set dev_module.PATH_INSTALLED

    Given install-prefix or <build_dir>-install
    """
    install_dir = get_dirs(args)[0]
    # set dev_module Global
    dev_module.PATH_INSTALLED = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        install_dir
    )
    return str(install_dir) # FIXME use returned value instead of module global

# FIXME: remove
def get_site_dir():
    path_installed = dev_module.PATH_INSTALLED
    # relative path for site-package with py version. i.e. 'lib/python3.10/site-packages'
    py_lib_path = Path(get_path('platlib')).relative_to(sys.exec_prefix)
    return str(Path(path_installed) / py_lib_path)


############

@cli.cls_cmd('build')
class Build(Task):
    """:wrench: build & install package on path"""
    ctx = CONTEXT

    werror = Option(['--werror'], default=False, is_flag=True, help="Treat warnings as errors")
    gcov = Option(
        ['--gcov'], default=False, is_flag=True,
        help="enable C code coverage via gcov (requires GCC)."
             "gcov output goes to build/**/*.gc*")
    debug = Option(['--debug', '-d'], default=False, is_flag=True, help="Debug build")
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='PARALLEL',
        help="Number of parallel jobs for build and testing")
    show_build_log = Option(
        ['--show-build-log'], default=False, is_flag=True,
        help="Show build output rather than using a log file")
    win_cp_openblas = Option(
        ['--win-cp-openblas'], default=False, is_flag=True,
        help="If set, and on Windows, copy OpenBLAS lib to install directory after"
             "meson install. Note: this argument may be removed in the future once a"
             "`site.cfg`-like mechanism to select BLAS/LAPACK libraries is"
             "implemented for Meson")

    @classmethod
    def run(cls, add_path=False, **kwargs):
        """return site_dir"""
        kwargs.update(cls.ctx.get(kwargs))
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        set_installed(args)
        site_dir = dev_module.build_project(args)

        # add site to sys.path
        if add_path:
            sys.path.insert(0, site_dir)
            os.environ['PYTHONPATH'] = \
                os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))
        return site_dir



PROJECT_MODULE = "scipy"
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))

@cli.cls_cmd('test')
class Test(Task):
    """:wrench: Run tests

    Examples:

    $ python do.py -s {SAMPLE_SUBMODULE}

    $ python do.py -s stats
    """
    ctx = CONTEXT

    verbose = Option(['--verbose', '-v'], default=False, is_flag=True, help="more verbosity")
    # removed doctests as currently not supported by _lib/_testutils.py
    # doctests = Option(['--doctests'], default=False, is_flag=True, help="Run doctests in module")
    coverage = Option(
        ['--coverage'], default=False, is_flag=True,
        help="report coverage of project code. HTML output goes under build/coverage")
    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")
    tests = Option(['--tests', '-t'], default=None, multiple=True, metavar='TESTS', help='Specify tests to run')
    mode = Option(
        ['--mode', '-m'], default='fast', metavar='MODE', show_default=True,
        help="'fast', 'full', or something that could be passed to `pytest -m` as a marker expression")
    parallel = Option(
        ['--parallel', '-j'], default=1, metavar='PARALLEL',
        help="Number of parallel jobs for testing"
    )

    TASK_META = {
        'task_dep': ['build'],
    }

    @staticmethod
    def _get_test_runner(path_installed, project_module):
        """
        get Test Runner from locally installed/built project
        """
        __import__(project_module)
        test = sys.modules[project_module].test
        version = sys.modules[project_module].__version__
        mod_path = sys.modules[project_module].__file__
        mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))
        return test, version, mod_path


    @classmethod
    def scipy_tests(cls, args):

        # if NOT BUID:
        #     test_dir = os.path.join(ROOT_DIR, args.build_dir, 'test')
        #     if not os.path.isdir(test_dir):
        #         os.makedirs(test_dir)
        site_dir = get_site_dir()
        # add local installed dir to PYTHONPATH
        print(f"Trying to find scipy from development installed path at: {site_dir}")
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))

        test_dir = site_dir

        # TODO: extra arguments to pytest?
        extra_argv = []
        # extra_argv = args.args[:]
        # if extra_argv and extra_argv[0] == '--':
        #     extra_argv = extra_argv[1:]

        if args.coverage:
            dst_dir = os.path.join(ROOT_DIR, args.build_dir, 'coverage')
            fn = os.path.join(dst_dir, 'coverage_html.js')
            if os.path.isdir(dst_dir) and os.path.isfile(fn):
                shutil.rmtree(dst_dir)
            extra_argv += ['--cov-report=html:' + dst_dir]

            shutil.copyfile(os.path.join(ROOT_DIR, '.coveragerc'),
                            os.path.join(test_dir, '.coveragerc'))

        runner, version, mod_path = cls._get_test_runner(dev_module.PATH_INSTALLED, PROJECT_MODULE)

        # convert options to test selection
        if args.submodule:
            tests = [PROJECT_MODULE + "." + args.submodule]
        elif args.tests:
            tests = args.tests
        else:
            tests = None

        # FIXME: changing CWD is not a good practice and might messed up with other tasks
        cwd = os.getcwd()
        try:
            os.chdir(test_dir)
            print("Running tests for {} version:{}, installed at:{}".format(
                        PROJECT_MODULE, version, mod_path))
            verbose = int(args.verbose) + 1 # runner verbosity - convert bool to int
            result = runner(
                args.mode,
                verbose=verbose,
                extra_argv=extra_argv,
                doctests=False, # args.doctests,
                coverage=args.coverage,
                tests=tests,
                parallel=args.parallel)
        finally:
            os.chdir(cwd)
        return result


    @classmethod
    def run(cls, **kwargs):
        """run unit-tests"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        cls.scipy_tests(args)

"""
**********************Bench taks*************************
Needs more work (WIP)
TODO: Remove redundancy
      Fix bench compare
"""
@cli.cls_cmd('bench')
class Bench(Task):
    """:wrench: Run benchmarks
    """
    ctx = CONTEXT
    TASK_META = {
        'task_dep': ['build'],
    }
    submodule = Option(
        ['--submodule', '-s'],
        default=None,
        metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")
    tests = Option(['--tests', '-t'],
                   default=None,
                   multiple=True,
                   metavar='TESTS',
                   help='Specify tests to run')
    bench_compare = Option(['--bench-compare', '-compare'],
                           default=None,
                           metavar='BENCH-COMPARE',
                           multiple=True,
                           help='Compare benchmark results of current HEAD to BEFORE')

    @staticmethod
    def _get_test_runner(path_installed, project_module):
        """
        get Test Runner from locally installed/built project
        """
        __import__(project_module)
        test = sys.modules[project_module].test
        version = sys.modules[project_module].__version__
        mod_path = sys.modules[project_module].__file__
        mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))
        return test, version, mod_path

    @staticmethod
    def run_asv(cmd):
        EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
                      '/usr/local/lib/ccache', '/usr/local/lib/f90cache']
        cwd = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                           'benchmarks')
        # Always use ccache, if installed
        env = dict(os.environ)
        env['PATH'] = os.pathsep.join(EXTRA_PATH +
                                      env.get('PATH', '').split(os.pathsep))
        # Control BLAS/LAPACK threads
        env['OPENBLAS_NUM_THREADS'] = '1'
        env['MKL_NUM_THREADS'] = '1'

        # Limit memory usage
        sys.path.insert(0, cwd)
        from benchmarks.common import set_mem_rlimit
        try:
            set_mem_rlimit()
        except (ImportError, RuntimeError):
            pass

        # Run
        try:
            return subprocess.call(cmd, env=env, cwd=cwd)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print("Error when running '%s': %s\n" % (" ".join(cmd), str(err),))
                print("You need to install Airspeed Velocity (https://airspeed-velocity.github.io/asv/)")
                print("to run Scipy benchmarks")
                return 1
            raise

    @classmethod
    def scipy_bench(cls, args):
        site_dir = get_site_dir()
        # add local installed dir to PYTHONPATH
        print(f"Trying to find scipy from development installed path at: {site_dir}")
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))
        cwd = os.getcwd()
        os.chdir(site_dir)
        runner, version, mod_path = cls._get_test_runner(dev_module.PATH_INSTALLED, PROJECT_MODULE)
        print("Running tests for {} version:{}, installed at:{}".format(
            PROJECT_MODULE, version, mod_path))

        extra_argv = []
        if args.tests:
            extra_argv.append(args.tests)
        if args.submodule:
            extra_argv.append([args.submodule])

        bench_args = []
        for a in extra_argv:
            bench_args.extend(['--bench', a])
        if not args.bench_compare:
            import scipy

            print("Running benchmarks for Scipy version %s at %s"
                  % (version, mod_path))
            cmd = ['asv', 'run', '--dry-run', '--show-stderr',
                   '--python=same'] + bench_args
            retval = cls.run_asv(cmd)
            sys.exit(retval)
        else:
            if len(args.bench_compare) == 1:
                commit_a = args.bench_compare[0]
                commit_b = 'HEAD'
            elif len(args.bench_compare) == 2:
                commit_a, commit_b = args.bench_compare
            else:
                p.error("Too many commits to compare benchmarks for")

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

            cmd = ['asv', 'continuous', '--show-stderr', '--factor', '1.05',
                   commit_a, commit_b] + bench_args
            run_asv(cmd)
            sys.exit(1)

    @classmethod
    def run(cls, **kwargs):
        """run benchamark"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        cls.scipy_bench(args)

###################
#### linters

@task_params([{'name': 'output_file', 'long': 'output-file', 'default': None,
               'help': 'Redirect report to a file'}])
def task_pep8(output_file):
    """Perform pep8 check with flake8."""
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
        'actions': [os.path.join(ROOT_DIR, 'tools', 'lint_diff.py')],
        'doc': 'Lint only files modified since last commit (with stricker rules)',
        'task_dep': ['pep8'],
    }

@cli.cls_cmd('pep8')
class Pep8():
    """Perform pep8 check with flake8."""
    output_file = Option(['--output-file'], default=None, help= 'Redirect report to a file')
    def run(output_file):
        opts = {'output_file': output_file}
        run_doit_task({'pep8-diff': {}, 'pep8': opts})


@cli.cls_cmd('mypy')
class Mypy(Task):
    """:wrench: Run mypy on the codebase"""
    ctx = CONTEXT

    TASK_META = {
        'task_dep': ['build'],
        # 'stream': 'no-capture',
    }

    @classmethod
    def run(cls):
        try:
            import mypy.api
        except ImportError as e:
            raise RuntimeError(
                "Mypy not found. Please install it by running "
                "pip install -r mypy_requirements.txt from the repo root"
            ) from e

        site_dir = get_site_dir()
        config = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "mypy.ini",
        )

        check_path = PROJECT_MODULE
        # debug
        check_path += '/optimize/minpack2.py'

        runtests = dev_module.import_module_from_path('runtests', Path(ROOT_DIR) / 'runtests.py')
        with runtests.working_dir(site_dir):
            # By default mypy won't color the output since it isn't being
            # invoked from a tty.
            os.environ['MYPY_FORCE_COLOR'] = '1'
            # Change to the site directory to make sure mypy doesn't pick
            # up any type stubs in the source tree.
            report, errors, status = mypy.api.run([
                "--config-file",
                config,
                check_path,
            ])
        print(report, end='')
        print(errors, end='', file=sys.stderr)
        return status == 0



##########################################
### DOC

@cli.cls_cmd('doc')
class Doc(Task):
    """:wrench: Build documentation"""
    ctx = CONTEXT

    # FIXME
    # parser.add_argument("--doc", action="append", nargs="?",
    #                 const="html-scipyorg", help="Build documentation")


    @classmethod
    def task_meta(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        install_dir, site_dir = get_dirs(args)

        # FIXME: support parallel
        # if args.parallel:
        #     cmd.append('SPHINXOPTS="-j{}"'.format(args.parallel))

        return {
            'actions': [
                # move to doc/ so local scipy is not get imported
                f'cd doc; env PYTHONPATH="../{site_dir}" make PYTHON="{sys.executable}" html-scipyorg',
            ],
            'task_dep': ['build'],
        }

@cli.cls_cmd('refguide-check')
class RefguideCheck(Task):
    """:wrench: Run refguide check (do not run regular tests.)"""
    ctx = CONTEXT

    submodule = Option(
        ['--submodule', '-s'], default=None, metavar='SUBMODULE',
        help="Submodule whose tests to run (cluster, constants, ...)")

    @classmethod
    def task_meta(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        install_dir, site_dir = get_dirs(args)

        cmd = [os.path.join(ROOT_DIR, 'tools', 'refguide_check.py'),
               '--doctests']
        if args.submodule:
            cmd += [args.submodule]
        cmd_str = ' '.join(cmd)
        return {
            'actions': [f'env PYTHONPATH={site_dir} {cmd_str}'],
            'task_dep': ['build'],
        }


##########################################
### ENVS
# TODO: it would be nice if click supported sections for commands
# https://stackoverflow.com/questions/57066951/divide-click-commands-into-sections-in-cli-documentation

@cli.command()
@click.argument('extra_argv', nargs=-1)
@click.pass_obj
def python(ctx_obj, extra_argv):
    """:wrench: Start a Python shell with PYTHONPATH set"""
    # not a doit task - manually build
    vals = Build.opt_defaults()
    vals.update(ctx_obj)
    Build.run(add_path=True, **vals)

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


@cli.command()
@click.pass_obj
def ipython(ctx_obj):
    """:wrench: Start IPython shell with PYTHONPATH set"""
    # not a doit task - manually build
    vals = Build.opt_defaults()
    vals.update(ctx_obj)
    Build.run(add_path=True, **vals)

    import IPython
    IPython.embed(user_ns={})


@cli.command()
@click.argument('extra_argv', nargs=-1)
@click.pass_obj
def shell(ctx_obj, extra_argv):
    """:wrench: Start Unix shell with PYTHONPATH set"""
    # not a doit task - manually build
    vals = Build.opt_defaults()
    vals.update(ctx_obj)
    Build.run(add_path=True, **vals)

    shell = os.environ.get('SHELL', 'sh')
    print("Spawning a Unix shell...")
    os.execv(shell, [shell] + extra_argv)
    sys.exit(1)


@cli.command()
@click.argument('version_args', nargs=2)
@click.pass_obj
def notes(ctx_obj, version_args):
    """:ledger: Release notes and log generation"""
    if version_args:
        sys.argv = version_args
        log_start = sys.argv[0]
        log_end = sys.argv[1]
    cmd = f"python tools/write_release_and_log.py v{log_start} v{log_end}"
    click.echo(cmd)
    try:
        subprocess.run([cmd], check = True, shell=True)
    except subprocess.CalledProcessError:
        print('Error caught: Incorrect log start or log end version')


@cli.command()
@click.argument('revision_args', nargs=2)
@click.pass_obj
def authors(ctx_obj, revision_args):
    """:ledger: Task to generate list the authors who contributed within a given revision interval"""
    if revision_args:
        sys.argv = revision_args
        start_revision = sys.argv[0]
        end_revision = sys.argv[1]
    cmd = f"python tools/authors.py v{start_revision}..v{end_revision}"
    click.echo(cmd)
    try:
        subprocess.run([cmd], check = True, shell=True)
    except subprocess.CalledProcessError:
        print('Error caught: Incorrect revision start or revision end')

if __name__ == '__main__':
    cli()
