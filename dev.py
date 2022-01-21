#!/usr/bin/env python
"""
dev.py [OPTIONS] [-- ARGS]

Run tests, building the project first with Meson

Examples::

    $ python dev.py
    $ python dev.py -s {SAMPLE_SUBMODULE}
    $ python dev.py -t {SAMPLE_TEST}
    $ python dev.py --ipython
    $ python dev.py --python somescript.py
    $ python dev.py --bench
    $ python dev.py --no-build --bench signal.LTI

Run a debugger:

    $ gdb --args python dev.py [...other args...]

Generate C code coverage listing under build/lcov/:
(requires http://ltp.sourceforge.net/coverage/lcov.php)

    $ python dev.py --gcov [...other args...]
    $ python dev.py --lcov-html

"""

#
# This is a generic test runner script for projects using NumPy's test
# framework. Change the following values to adapt to your project:
#

PROJECT_MODULE = "scipy"
PROJECT_ROOT_FILES = ['scipy', 'LICENSE.txt', 'meson.build']
SAMPLE_TEST = "scipy.fftpack.tests.test_real_transforms::TestIDSTIIIInt"
SAMPLE_SUBMODULE = "optimize"

EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
              '/usr/local/lib/ccache', '/usr/local/lib/f90cache']

# ---------------------------------------------------------------------


if __doc__ is None:
    __doc__ = "Run without -OO if you want usage info"
else:
    __doc__ = __doc__.format(**globals())


import sys
import os
# the following multiprocessing import is necessary to prevent tests that use
# multiprocessing from hanging on >= Python3.8 (macOS) using pytest. Just the
# import is enough...
import multiprocessing


# In case we are run from the source directory, we don't want to import the
# project from there:
sys.path.pop(0)
current_sys_path = sys.path.copy()

from argparse import ArgumentParser, REMAINDER
import contextlib
import shutil
import subprocess
import time
import datetime
import importlib.util

try:
    from types import ModuleType as new_module
except ImportError:  # old Python
    from imp import new_module

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))
PATH_INSTALLED = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                   'installdir')

# To import runtests.py
fname = os.path.join(ROOT_DIR, 'runtests.py')
spec = importlib.util.spec_from_file_location('runtests', str(fname))
runtests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtests)

# Reassign sys.path as it is changed by importlib above
sys.path = current_sys_path

def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("--verbose", "-v", action="count", default=1,
                        help="more verbosity")
    parser.add_argument("--no-build", "-n", action="store_true", default=False,
                        help="do not build the project (use system installed version)")
    parser.add_argument("--werror", action="store_true", default=False,
                        help="Treat warnings as errors")
    parser.add_argument("--build-only", "-b", action="store_true", default=False,
                        help="just build, do not run any tests")
    parser.add_argument("--doctests", action="store_true", default=False,
                        help="Run doctests in module")
    parser.add_argument("--refguide-check", action="store_true", default=False,
                        help="Run refguide check (do not run regular tests.)")
    parser.add_argument("--coverage", action="store_true", default=False,
                        help=("report coverage of project code. HTML output"
                              " goes under build/coverage"))
    parser.add_argument("--gcov", action="store_true", default=False,
                        help=("enable C code coverage via gcov (requires GCC)."
                              " gcov output goes to build/**/*.gc*"))
    parser.add_argument("--lcov-html", action="store_true", default=False,
                        help=("produce HTML for C code coverage information "
                              "from a previous run with --gcov. "
                              "HTML output goes to build/lcov/"))
    parser.add_argument("--mode", "-m", default="fast",
                        help="'fast', 'full', or something that could be "
                             "passed to nosetests -A [default: fast]")
    parser.add_argument("--submodule", "-s", default=None,
                        help="Submodule whose tests to run (cluster,"
                             " constants, ...)")
    parser.add_argument("--pythonpath", "-p", default=None,
                        help="Paths to prepend to PYTHONPATH")
    parser.add_argument("--tests", "-t", action='append',
                        help="Specify tests to run")
    parser.add_argument("--python", action="store_true",
                        help="Start a Python shell with PYTHONPATH set")
    parser.add_argument("--ipython", "-i", action="store_true",
                        help="Start IPython shell with PYTHONPATH set")
    parser.add_argument("--shell", action="store_true",
                        help="Start Unix shell with PYTHONPATH set")
    parser.add_argument("--debug", "-g", action="store_true",
                        help="Debug build")
    parser.add_argument("--parallel", "-j", type=int, default=1,
                        help="Number of parallel jobs for build and testing")
    parser.add_argument("--show-build-log", action="store_true",
                        help="Show build output rather than using a log file")
    parser.add_argument("--bench", action="store_true",
                        help="Run benchmark suite instead of test suite")
    parser.add_argument("--bench-compare", action="append", metavar="BEFORE",
                        help=("Compare benchmark results of current HEAD to"
                              " BEFORE. Use an additional "
                              "--bench-compare=COMMIT to override HEAD with"
                              " COMMIT. Note that you need to commit your "
                              "changes first!"
                              ))
    parser.add_argument("args", metavar="ARGS", default=[], nargs=REMAINDER,
                        help="Arguments to pass to Nose, Python or shell")
    parser.add_argument("--pep8", action="store_true", default=False,
                        help="Perform pep8 check with flake8.")
    parser.add_argument("--mypy", action="store_true", default=False,
                        help="Run mypy on the codebase")
    parser.add_argument("--doc", action="append", nargs="?",
                        const="html-scipyorg", help="Build documentation")
    args = parser.parse_args(argv)

    if args.pep8:
        # Lint the source using the configuration in tox.ini.
        os.system("flake8 scipy benchmarks/benchmarks")
        # Lint just the diff since branching off of main using a
        # stricter configuration.
        lint_diff = os.path.join(ROOT_DIR, 'tools', 'lint_diff.py')
        os.system(lint_diff)
        sys.exit(0)

    if args.mypy:
        sys.exit(run_mypy(args))

    if args.bench_compare:
        args.bench = True
        args.no_build = True  # ASV does the building

    if args.lcov_html:
        # generate C code coverage output
        runtests.lcov_generate()
        sys.exit(0)

    if args.pythonpath:
        for p in reversed(args.pythonpath.split(os.pathsep)):
            sys.path.insert(0, p)

    if args.gcov:
        runtests.gcov_reset_counters()

    if args.debug and args.bench:
        print("*** Benchmarks should not be run against debug version; "
              "remove -g flag ***")

    if not args.no_build:
        site_dir = build_project(args)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))

    extra_argv = args.args[:]
    if extra_argv and extra_argv[0] == '--':
        extra_argv = extra_argv[1:]

    if args.python:
        if extra_argv:
            # Don't use subprocess, since we don't want to include the
            # current path in PYTHONPATH.
            sys.argv = extra_argv
            with open(extra_argv[0], 'r') as f:
                script = f.read()
            sys.modules['__main__'] = new_module('__main__')
            ns = dict(__name__='__main__',
                      __file__=extra_argv[0])
            exec(script, ns)
            sys.exit(0)
        else:
            import code
            code.interact()
            sys.exit(0)

    if args.ipython:
        import IPython
        IPython.embed(user_ns={})
        sys.exit(0)

    if args.shell:
        shell = os.environ.get('SHELL', 'sh')
        print("Spawning a Unix shell...")
        os.execv(shell, [shell] + extra_argv)
        sys.exit(1)

    if args.doc:
        cmd = ["make", "-Cdoc", 'PYTHON="{}"'.format(sys.executable)]
        cmd += args.doc
        if args.parallel:
            cmd.append('SPHINXOPTS="-j{}"'.format(args.parallel))
        subprocess.run(cmd, check=True)
        sys.exit(0)

    if args.coverage:
        dst_dir = os.path.join(ROOT_DIR, 'build', 'coverage')
        fn = os.path.join(dst_dir, 'coverage_html.js')
        if os.path.isdir(dst_dir) and os.path.isfile(fn):
            shutil.rmtree(dst_dir)
        extra_argv += ['--cov-report=html:' + dst_dir]

    if args.refguide_check:
        cmd = [os.path.join(ROOT_DIR, 'tools', 'refguide_check.py'),
               '--doctests']
        if args.submodule:
            cmd += [args.submodule]
        os.execv(sys.executable, [sys.executable] + cmd)
        sys.exit(0)

    if args.bench:
        # Run ASV
        items = extra_argv
        if args.tests:
            items += args.tests
        if args.submodule:
            items += [args.submodule]

        bench_args = []
        for a in items:
            bench_args.extend(['--bench', a])

        if not args.bench_compare:
            import scipy
            print("Running benchmarks for Scipy version %s at %s"
                  % (scipy.__version__, scipy.__file__))
            cmd = ['asv', 'run', '--dry-run', '--show-stderr',
                   '--python=same'] + bench_args
            retval = runtests.run_asv(cmd)
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
                    print("*"*80)
                    print("WARNING: you have uncommitted changes --- "
                          "these will NOT be benchmarked!")
                    print("*"*80)

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
            runtests.run_asv(cmd)
            sys.exit(1)

    if args.build_only:
        sys.exit(0)
    else:
        try:
            __import__(PROJECT_MODULE)
            test = sys.modules[PROJECT_MODULE].test
            version = sys.modules[PROJECT_MODULE].__version__
            mod_path = sys.modules[PROJECT_MODULE].__file__
            mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))
        except ImportError:
            current_python_path = os.environ.get('PYTHONPATH', None)
            print("Unable to import {} from: {}".format(PROJECT_MODULE,
                                                       current_python_path))
            from sysconfig import get_path
            py_path = get_path('platlib')
            site_dir = os.path.join(PATH_INSTALLED,
                                    runtests.get_path_suffix(py_path, 3))
            print("Trying to import scipy from development installed path at:",
                  site_dir)
            sys.path.insert(0, site_dir)
            os.environ['PYTHONPATH'] = \
                os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))
            __import__(PROJECT_MODULE)
            test = sys.modules[PROJECT_MODULE].test
            version = sys.modules[PROJECT_MODULE].__version__
            mod_path = sys.modules[PROJECT_MODULE].__file__
            mod_path = os.path.abspath(os.path.join(os.path.dirname(mod_path)))


    if args.submodule:
        tests = [PROJECT_MODULE + "." + args.submodule]
    elif args.tests:
        tests = args.tests
    else:
        tests = None

    # Run the tests

    if not args.no_build:
        test_dir = site_dir
    else:
        test_dir = os.path.join(ROOT_DIR, 'build', 'test')
        if not os.path.isdir(test_dir):
            os.makedirs(test_dir)

    shutil.copyfile(os.path.join(ROOT_DIR, '.coveragerc'),
                    os.path.join(test_dir, '.coveragerc'))

    cwd = os.getcwd()
    try:
        os.chdir(test_dir)
        print("Running tests for {} version:{}, installed at:{}".format(
                    PROJECT_MODULE, version, mod_path))
        result = test(args.mode,
                      verbose=args.verbose,
                      extra_argv=extra_argv,
                      doctests=args.doctests,
                      coverage=args.coverage,
                      tests=tests,
                      parallel=args.parallel)
    finally:
        os.chdir(cwd)

    if isinstance(result, bool):
        sys.exit(0 if result else 1)
    elif result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)


def setup_build(args, env):
    """
    Setting up meson-build
    """
    cmd = ["meson", "setup", "build", "--prefix", PATH_INSTALLED]
    if args.werror:
        cmd += ["--werror"]
    # Setting up meson build
    ret = subprocess.call(cmd, env=env, cwd=ROOT_DIR)
    if ret == 0:
        print("Meson build setup OK")
    else:
        print("Meson build setup failed! ({0} elapsed)")
        sys.exit(1)
    return


def install_project(show_build_log):
    """
    Installs the project after building.
    """
    cmd = ["meson", "install", "-C", "build"]
    log_filename = os.path.join(ROOT_DIR, 'meson-install.log')
    start_time = datetime.datetime.now()
    if show_build_log:
        ret = subprocess.call(cmd, cwd=ROOT_DIR)
    else:
        print("Installing, see meson-install.log...")
        with open(log_filename, 'w') as log:
            p = subprocess.Popen(cmd, stdout=log, stderr=log,
                                 cwd=ROOT_DIR)

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

    if ret == 0:
        print("Installation OK")
    else:
        if not show_build_log:
            with open(log_filename, 'r') as f:
                print(f.read())
        print("Installation failed! ({0} elapsed)".format(elapsed))
        sys.exit(1)
    return


def build_project(args):
    """
    Build a dev version of the project.

    Returns
    -------
    site_dir
        site-packages directory where it was installed

    """

    root_ok = [os.path.exists(os.path.join(ROOT_DIR, fn))
               for fn in PROJECT_ROOT_FILES]
    if not all(root_ok):
        print("To build the project, run dev.py in "
              "git checkout or unpacked source")
        sys.exit(1)

    env = dict(os.environ)

    if args.debug or args.gcov:
        # assume everyone uses gcc/gfortran
        env['OPT'] = '-O0 -ggdb'
        env['FOPT'] = '-O0 -ggdb'
        if args.gcov:
            from sysconfig import get_config_vars
            cvars = get_config_vars()
            env['OPT'] = '-O0 -ggdb'
            env['FOPT'] = '-O0 -ggdb'
            env['CC'] = env.get('CC', cvars['CC']) + ' --coverage'
            env['CXX'] = env.get('CXX', cvars['CXX']) + ' --coverage'
            env['F77'] = 'gfortran --coverage '
            env['F90'] = 'gfortran --coverage '
            env['LDSHARED'] = cvars['LDSHARED'] + ' --coverage'
            env['LDFLAGS'] = " ".join(cvars['LDSHARED'].split()[1:]) +\
                ' --coverage'

    build_dir = os.path.join(ROOT_DIR, 'build')

    # Check if meson is already setup
    if not os.path.exists(os.path.join(build_dir, 'build.ninja')):
        setup_build(args, env)

    from sysconfig import get_path
    py_path = get_path('platlib')
    site_dir = os.path.join(PATH_INSTALLED,
                            runtests.get_path_suffix(py_path, 3))

    cmd = ["ninja", "-C", "build"]
    if args.parallel > 1:
        cmd += ["-j", str(args.parallel)]

    # Building with ninja-backend
    ret = subprocess.call(cmd, env=env, cwd=ROOT_DIR)

    if ret == 0:
        print("Build OK")
    else:
        print("Build failed!")
        sys.exit(1)

    install_project(args.show_build_log)

    return site_dir


def run_mypy(args):
    if args.no_build:
        raise ValueError('Cannot run mypy with --no-build')

    try:
        import mypy.api
    except ImportError as e:
        raise RuntimeError(
            "Mypy not found. Please install it by running "
            "pip install -r mypy_requirements.txt from the repo root"
        ) from e

    site_dir = build_project(args)
    config = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "mypy.ini",
    )
    with runtests.working_dir(site_dir):
        # By default mypy won't color the output since it isn't being
        # invoked from a tty.
        os.environ['MYPY_FORCE_COLOR'] = '1'
        # Change to the site directory to make sure mypy doesn't pick
        # up any type stubs in the source tree.
        report, errors, status = mypy.api.run([
            "--config-file",
            config,
            PROJECT_MODULE,
        ])
    print(report, end='')
    print(errors, end='', file=sys.stderr)
    return status


if __name__ == "__main__":
    main(argv=sys.argv[1:])
