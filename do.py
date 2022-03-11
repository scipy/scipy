#! /usr/bin/env python3

"""
Developer CLI: building (meson), tests, benchmark, etc.

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authring the development tasks.


"""



#####################################
# to be extracted into a doit_click.py.

import sys
import inspect

from doit.cmd_base import ModuleTaskLoader, get_loader
from doit.doit_cmd import DoitMain
from doit.cmdparse import DefaultUpdate, CmdParse, CmdParseError
from doit.exceptions import InvalidDodoFile, InvalidCommand, InvalidTask
import click
from click.globals import push_context

class DoitMainAPI(DoitMain):
    """add new method to run tasks with parsed command line"""

    def run_tasks(self, task):
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
            args = [task]
            command.execute(params, args)
        except (CmdParseError, InvalidDodoFile,
                InvalidCommand, InvalidTask) as err:
            if isinstance(err, InvalidCommand):
                err.cmd_used = cmd_name
                err.bin_name = self.BIN_NAME
            raise err


def run_doit_task(task_name, **kwargs):
    """:param kwargs: contain task_opts"""
    loader = ModuleTaskLoader(globals())
    loader.task_opts = {task_name: kwargs}
    doit_main = DoitMainAPI(loader, extra_config={'GLOBAL': {'verbosity': 2}})
    return doit_main.run_tasks(task_name)

def doit_task_callback(task_name):
    def callback(**kwargs):
        sys.exit(run_doit_task(task_name, **kwargs))
    return callback


def param_doit2click(task_param: dict):
    """converts a doit TaskParam to Click.Parameter"""
    param_decls = [task_param['name']]
    if 'long' in task_param:
        param_decls.append(f"--{task_param['long']}")
    if 'short' in task_param:
        param_decls.append(f"-{task_param['short']}")
    return click.Option(param_decls, default=task_param['default'], help=task_param.get('help'))

class ClickDoit(click.Group):
    """
    subclass with an extra decorator used to create commands from doit task_creators
    (instead of plain function)
    """

    def task_as_cmd(self, name=None, **attrs):
        def decorator(creator):
            task_name = creator.__name__[5:] # 5 is len('task_')
            cmd_name = name if name else task_name
            cmd_help = inspect.getdoc(creator)
            task_params = getattr(creator, '_task_creator_params', None)

            # convert doit task-params to click params
            params = []
            # if param is already define in group, do not add again
            group_options = set(par.name for par in self.params)
            if task_params:
                for tp in task_params:
                    if tp['name'] in group_options:
                        continue
                    params.append(param_doit2click(tp))
            cmd = click.Command(
                name=cmd_name,
                callback=doit_task_callback(task_name),
                help=cmd_help,
                params=params,
            )
            self.add_command(cmd)
            return creator # return original task_creator to be used by doit itself
        return decorator



#######################################
# SciPy tasks / CLI

import os
import sys
from sysconfig import get_path
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from doit import task_params

# keep compatibility and re-use dev.py code
import dev as dev_module

opt_build_dir = {
    'name': 'build_dir',
    'long': 'build-dir',
    'default': 'build',
    'help': "Relative path to build directory. Default is 'build'",
}

opt_install_prefix = {
    'name': 'install_prefix',
    'long': 'install-prefix',
    'default': None,
    'help': "Relative path to the install directory. Default is <build-dir>-install.",
}


# execute tasks through click based CLI
@click.group(cls=ClickDoit)
@click.option('--build-dir', default='build', help=opt_build_dir['help'])
@click.option('--install-prefix', default=None, help=opt_install_prefix['help'])
@click.pass_context
def cli(ctx, build_dir, install_prefix):
    ctx.ensure_object(dict)
    ctx.obj['build_dir'] = build_dir
    ctx.obj['install_prefix'] = install_prefix


@dataclass
class PathArgs:
    """PathArgs are shared by build & all other tasks that used built/installed project"""
    build_dir: str = 'build'
    install_prefix: str = None
    werror: str = None # Treat warnings as errors

def set_installed(args: PathArgs):
    """set dev_module.PATH_INSTALLED

    Given install-prefix or <build_dir>-install
    """
    build_dir = Path(args.build_dir)
    install_dir = args.install_prefix
    if not install_dir:
        install_dir = build_dir.parent / (build_dir.stem + "-install")
    # set dev_module Global
    dev_module.PATH_INSTALLED = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        install_dir
    )
    return str(install_dir) # FIXME use returned value instead of module global


@dataclass
class BuildArgs(PathArgs):
    """build configuration argument options"""
    debug: bool = False
    gcov: bool = False
    parallel: int = 1
    show_build_log: bool = False
    win_cp_openblas: bool = False



@cli.task_as_cmd()
@click.pass_obj
@task_params([opt_build_dir, opt_install_prefix])
def task_build(ctx_obj, build_dir, install_prefix):
    """build & install package on path"""
    args = BuildArgs(build_dir=build_dir, install_prefix=install_prefix)
    args.__dict__.update(**ctx_obj)

    # base = Path('scipy')
    # files = [str(f) for f in base.glob('**/*.py')]
    return {
        'actions': [
            (set_installed, (args,)),
            (dev_module.build_project, (args,)),
        ],
        # 'file_dep': files,
    }


##### Test ####

PROJECT_MODULE = "scipy"
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


@dataclass
class TestArgs(PathArgs):
    parallel: int = 1
    verbose: int = 1 # count
    doctests: bool = False
    coverage: bool = False
    # test selection
    mode: str = 'fast'
    submodule: str = None # -s
    tests: [str] = field(default_factory=list) # -t


def _get_test_runner(path_installed, project_module):
    """
    get Test Runner from locally installed/built project
    """
    # relative path for site-package with py version. i.e. 'lib/python3.10/site-packages'
    py_lib_path = Path(get_path('platlib')).relative_to(sys.exec_prefix)
    site_dir = str(Path(path_installed) / py_lib_path)

    # add local installed dir to PYTHONPATH
    # print(f"Trying to find scipy from development installed path at: {site_dir}")
    sys.path.insert(0, site_dir)
    os.environ['PYTHONPATH'] = \
        os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))

    __import__(project_module)
    test = sys.modules[project_module].test
    version = sys.modules[project_module].__version__
    mod_path = sys.modules[project_module].__file__
    mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))
    return test, version, mod_path


def scipy_tests(args: TestArgs):
    test_dir = os.path.join(ROOT_DIR, args.build_dir, 'test')
    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)

    shutil.copyfile(os.path.join(ROOT_DIR, '.coveragerc'),
                    os.path.join(test_dir, '.coveragerc'))

    runner, version, mod_path = _get_test_runner(dev_module.PATH_INSTALLED, PROJECT_MODULE)

    # TODO: extra arguments to pytest?
    extra_argv = []
    # extra_argv = args.args[:]
    # if extra_argv and extra_argv[0] == '--':
    #     extra_argv = extra_argv[1:]

    if args.submodule:
        tests = [PROJECT_MODULE + "." + args.submodule]
    elif args.tests:
        tests = args.tests
    else:
        tests = None

    print('Testsssss', tests)
    return

    # FIXME: changing CWD is not a good practice and might messed up with other tasks
    cwd = os.getcwd()
    try:
        os.chdir(test_dir)
        print("Running tests for {} version:{}, installed at:{}".format(
                    PROJECT_MODULE, version, mod_path))
        result = runner(
            args.mode,
            verbose=args.verbose,
            extra_argv=extra_argv,
            doctests=args.doctests,
            coverage=args.coverage,
            tests=tests,
            parallel=args.parallel)
    finally:
        os.chdir(cwd)


@cli.task_as_cmd()
@click.pass_obj
@task_params([
    opt_build_dir, opt_install_prefix,
    {'name': 'submodule', 'long': 'submodule', 'short': 's', 'default': None,
     'help': 'Submodule whose tests to run (cluster, constants, ...)',
    },
])
def task_test(ctx_obj, submodule=None, **kwargs):
    """run unit-tests"""
    args = TestArgs(submodule=submodule, **kwargs)
    args.__dict__.update(**ctx_obj)
    return {
        'actions': [
            (set_installed, (args,)),
            (scipy_tests, (args,)),
        ],
        'verbosity': 2,
        'task_dep': ['build'],
    }


#####################
# Example of click command that does not use doit at all
@cli.command()
def pep8():
    """pep8 & lint (not a doit task)"""
    import os
    # Lint the source using the configuration in tox.ini.
    os.system("flake8 scipy benchmarks/benchmarks")
    # Lint just the diff since branching off of main using a
    # stricter configuration.
    lint_diff = os.path.join(ROOT_DIR, 'tools', 'lint_diff.py')
    os.system(lint_diff)


###################
## Example of task / click command being defined separately

@task_params([{'name': 'output_file', 'long': 'output-file', 'default': None,
               'help': 'Redirect report to a file'}])
def task_flake8(output_file):
    """Perform pep8 check with flake8."""
    opts = ''
    if output_file:
        opts += f'--output-file={output_file}'
    return {
        'actions': [f"flake8 {opts} scipy benchmarks/benchmarks"],
    }

@cli.command()
@click.option('output-file', default=None, help= 'Redirect report to a file')
def flake8(output_file):
    """Perform pep8 check with flake8."""
    opts = {'outfile': output_file}
    sys.exit(run_doit_task('flake8', opts))


##################
# DRY exposing simple task as command

@cli.task_as_cmd()
@task_params([{'name': 'output_file', 'long': 'output-file', 'default': None,
               'help': 'Redirect report to a file'}])
def task_flake(output_file):
    """Perform pep8 check with flake8."""
    opts = ''
    if output_file:
        opts += f'--output-file={output_file}'
    return {
        'actions': [f"flake8 {opts} scipy benchmarks/benchmarks"],
    }


# def task_bench():
#     pass

# def task_doc_build():
#     pass



# build using dev.py command line
# def task_builddev():
#     return {
#         'actions': ['python dev.py --build-only'],
#     }


if __name__ == '__main__':
    cli()
else:
    # make sure click.pass_obj decorator still works when running with plain doit
    push_context(click.core.Context(cli, obj={}))
