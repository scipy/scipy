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

from doit.cmd_base import ModuleTaskLoader
from doit.doit_cmd import DoitMain
import click

def run_doit_task(task_name):
    def callback():
        loader = ModuleTaskLoader(globals())
        doit_main = DoitMain(loader, extra_config={'GLOBAL': {'verbosity': 2}})
        sys.exit(doit_main.run([task_name]))
    return callback


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
            cmd = click.Command(name=cmd_name, callback=run_doit_task(task_name), help=cmd_help)
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

# keep compatibility and re-use dev.py code
import dev as dev_module

# execute tasks through click based CLI
cli = ClickDoit()

@dataclass
class PathArgs:
    build_dir: str = 'build'
    install_prefix: str = None

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
def task_build():
    """build & install package on path"""
    args = BuildArgs()
    return {
        'actions': [
            (set_installed, (args,)),
            (dev_module.build_project, (args,)),
        ],
    }


##### Test ####

PROJECT_MODULE = "scipy"
PROJECT_ROOT_FILES = ['scipy', 'LICENSE.txt', 'meson.build']
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
def task_test():
    args = TestArgs()
    return {
        'actions': [
            (set_installed, (args,)),
            (scipy_tests, (args,)),
        ],
        'verbosity': 2,
    }

# def task_bench():
#     pass

# def task_doc_build():
#     pass


if __name__ == '__main__':
    cli()
