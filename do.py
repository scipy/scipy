#! /usr/bin/env python3

"""
Developer CLI: building (meson), tests, benchmark, etc.

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authring the development tasks.


Note this requires the unreleased doit 0.35


TODO: First milestone replace dev.py

- [ ] move out non-scipy code
- [ ] document API/reasoning for creating commands/tasks
- [ ] copy used code from dev.py
- [ ] doit reporter when running under click
- [ ] command sections (https://github.com/janluke/cloup or similar)

commands:
- [ ] lcov_html
- [ ] refguide_check
- [ ] bench

BUG:
- [ ] python dev.py -t scipy.optimize.tests.test_minimize_constrained.py::TestTrustRegionConstr::test_args
      unknown marker "slow". seems conftest is not being loaded
- [ ] click does not support '--'. used by test command
- [ ] doc: fail on python3.10 - distutils is deprecated
"""

import os
import sys
import re
import shutil
from sysconfig import get_path
from pathlib import Path
from collections import namedtuple

import click
from click import Option
from click.globals import get_current_context
from doit import task_params
from doit.task import Task as DoitTask
from doit.cmd_base import ModuleTaskLoader, get_loader
from doit.doit_cmd import DoitMain
from doit.cmdparse import DefaultUpdate, CmdParse, CmdParseError
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
    task_kwargs = getattr(cls, 'Meta', {})
    return DoitTask(
        # convert name to kebab-case
        name=to_camel(cls.__name__),
        doc=cls.__doc__,
        actions=[cls.run],
        params=cls._params,
        **task_kwargs,
    )

class MetaTask(type):

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

        if hasattr(cls, 'meta'):
            def creator(**kwargs):
                print('CREATOR got', kwargs)
                task_meta = cls.meta(**kwargs)
                if 'basename' not in task_meta:
                    task_meta['basename'] = to_camel(cls.__name__)
                if 'doc' not in task_meta:
                    task_meta['doc'] = cls.__doc__
                return task_meta
            creator._task_creator_params = cls._params
        else:
            def creator():
                return run_as_py_action(cls)
        cls.create_doit_tasks = creator
        return cls


class Task(metaclass=MetaTask):
    """Base class to define doit task and/or click command"""

    @classmethod
    def opt_defaults(cls):
        """helper used by another click commands to call this command"""
        return {p['name']:p['default'] for p in cls._params}



class CliGroup(click.Group):
    def cmd(self, name):
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

            click_cmd = click.Command(
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


CONTEXT = Context({
    'build_dir': Option(
        ['--build-dir'], default='build', metavar='BUILD_DIR', show_default=True,
        help='Relative path to the build directory.'),
    'install_prefix': Option(
        ['--install-prefix'], default=None, metavar='INSTALL_DIR',
        help="Relative path to the install directory. Default is <build-dir>-install."),
})



@click.group(cls=CliGroup)
@click.pass_context
def cli(ctx, **kwargs):
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

@cli.cmd('build')
class Build(Task):
    """build & install package on path"""
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

@cli.cmd('test')
class Test(Task):
    """Run tests"""
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

    Meta = {
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


    @classmethod
    def run(cls, **kwargs):
        """run unit-tests"""
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)
        cls.scipy_tests(args)



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

@cli.cmd('pep8')
class Pep8():
    """Perform pep8 check with flake8."""
    output_file = Option(['--output-file'], default=None, help= 'Redirect report to a file')
    def run(output_file):
        opts = {'output_file': output_file}
        run_doit_task({'pep8-diff': {}, 'pep8': opts})


@cli.cmd('mypy')
class Mypy(Task):
    """Run mypy on the codebase"""
    ctx = CONTEXT

    Meta = {
        'task_dep': ['build'],
        # 'stream': 'no-capture',
    }

    @classmethod
    def run(cls, **kwargs):
        kwargs.update(cls.ctx.get())
        Args = namedtuple('Args', [k for k in kwargs.keys()])
        args = Args(**kwargs)

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

@cli.cmd('doc')
class doc(Task):
    """Build documentation"""
    ctx = CONTEXT

    # FIXME
    # parser.add_argument("--doc", action="append", nargs="?",
    #                 const="html-scipyorg", help="Build documentation")


    @classmethod
    def meta(cls, **kwargs):
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


##########################################
### ENVS
# TODO: it would be nice if click supported sections for commands
# https://stackoverflow.com/questions/57066951/divide-click-commands-into-sections-in-cli-documentation

@cli.command()
@click.argument('extra_argv', nargs=-1)
@click.pass_obj
def python(ctx_obj, extra_argv):
    """Start a Python shell with PYTHONPATH set"""
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
    """Start IPython shell with PYTHONPATH set"""
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
    """Start Unix shell with PYTHONPATH set"""
    # not a doit task - manually build
    vals = Build.opt_defaults()
    vals.update(ctx_obj)
    Build.run(add_path=True, **vals)

    shell = os.environ.get('SHELL', 'sh')
    print("Spawning a Unix shell...")
    os.execv(shell, [shell] + extra_argv)
    sys.exit(1)



if __name__ == '__main__':
    cli()
