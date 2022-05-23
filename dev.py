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
import warnings  # noqa: E402
from pathlib import Path
import platform
# the following multiprocessing import is necessary to prevent tests that use
# multiprocessing from hanging on >= Python3.8 (macOS) using pytest. Just the
# import is enough...
import multiprocessing

# distutils is required to infer meson install path
# if this needs to be replaced for Python 3.12 support and there's no
# stdlib alternative, use the hack discussed in gh-16058
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    from distutils import dist
    from distutils.command.install import INSTALL_SCHEMES

# In case we are run from the source directory, we don't want to import the
# project from there:
sys.path.pop(0)
current_sys_path = sys.path.copy()

from argparse import ArgumentParser, REMAINDER
import shutil
import subprocess
import time
import datetime
import importlib.util
import json  # noqa: E402
from sysconfig import get_path
from types import ModuleType as new_module  # noqa: E402

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


def import_module_from_path(mod_name, mod_path):
    """Import module with name `mod_name` from file path `mod_path`"""
    spec = importlib.util.spec_from_file_location(mod_name, mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Import runtests.py
runtests = import_module_from_path('runtests', Path(ROOT_DIR) / 'runtests.py')

# Reassign sys.path as it is changed by the `runtests` import above
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
                             "passed to `pytest -m` as a marker expression "
                             "[default: fast]")
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
    parser.add_argument("--win-cp-openblas", action="store_true",
                        help="If set, and on Windows, copy OpenBLAS lib to "
                        "install directory after meson install. "
                        "Note: this argument may be removed in the future "
                        "once a `site.cfg`-like mechanism to select BLAS/LAPACK "
                        "libraries is implemented for Meson")
    parser.add_argument("--build-dir", default="build",
                        help="Relative path to the build directory. "
                             "Default is 'build'")
    parser.add_argument("--install-prefix", default=None,
                        help="Relative path to the install directory. "
                             "Default is <build-dir>-install.")
    args = parser.parse_args(argv)

    global PATH_INSTALLED
    build_dir = Path(args.build_dir)
    install_dir = args.install_prefix
    if not install_dir:
        install_dir = build_dir.parent / (build_dir.stem + "-install")
    PATH_INSTALLED = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        install_dir
    )

    if args.win_cp_openblas and platform.system() != 'Windows':
        raise RuntimeError('--win-cp-openblas only has effect on Windows')

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
        dst_dir = os.path.join(ROOT_DIR, args.build_dir, 'coverage')
        fn = os.path.join(dst_dir, 'coverage_html.js')
        if os.path.isdir(dst_dir) and os.path.isfile(fn):
            shutil.rmtree(dst_dir)
        extra_argv += ['--cov-report=html:' + dst_dir]

    if args.refguide_check:
        cmd = [os.path.join(ROOT_DIR, 'tools', 'refguide_check.py'),
               '--doctests']
        if args.verbose:
            cmd += ['-' + 'v'*args.verbose]
        if args.submodule:
            cmd += [args.submodule]
        os.execv(sys.executable, [sys.executable] + cmd)
        sys.exit(0)

    test, version, mod_path = get_project_info()

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
                  % (version, mod_path))
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
        test_dir = os.path.join(ROOT_DIR, args.build_dir, 'test')
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


def get_project_info():
    """
    Function to import the project module and return its tests, version,
    and path where it is found.
    If the project module is not found, then it tries to find it in the
    development installed path.
    """
    try:
        test, version, mod_path = runtests.import_module()
    except ImportError:
        # this may fail when running with --no-build, so try to detect
        # an installed scipy in a subdir inside a repo
        site_dir = get_site_packages()
        print("Trying to find scipy from development installed "
              "path at:", site_dir)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = \
            os.pathsep.join((site_dir, os.environ.get('PYTHONPATH', '')))
        test, version, mod_path = runtests.import_module()
    return test, version, mod_path


def setup_build(args, env):
    """
    Setting up meson-build
    """
    cmd = ["meson", "setup", args.build_dir, "--prefix", PATH_INSTALLED]
    build_dir = Path(args.build_dir)
    run_dir = os.getcwd()
    if build_dir.exists() and not (build_dir / 'meson-info').exists():
        if list(build_dir.iterdir()):
            raise RuntimeError(
                f"You're using Meson to build in the `{build_dir.absolute()}` directory, "
                "but it looks like that directory is not empty and "
                "was not originally created by Meson. "
                f"Please remove '{build_dir.absolute()}' and try again."
            )
    if os.path.exists(build_dir):
        build_options_file = (build_dir / "meson-info"
                              / "intro-buildoptions.json")
        if build_options_file.exists():
            with open(build_options_file) as f:
                build_options = json.load(f)
            installdir = None
            for option in build_options:
                if option["name"] == "prefix":
                    installdir = option["value"]
                    break
            if installdir != PATH_INSTALLED:
                run_dir = os.path.join(run_dir, build_dir)
                cmd = ["meson", "--reconfigure", "--prefix", PATH_INSTALLED]
            else:
                return
        else:
            run_dir = os.path.join(run_dir, build_dir)
            cmd = ["meson", "--reconfigure", "--prefix", PATH_INSTALLED]

    if args.werror:
        cmd += ["--werror"]
    if args.gcov:
        cmd += ['-Db_coverage=true']
    # Setting up meson build
    ret = subprocess.call(cmd, env=env, cwd=run_dir)
    if ret == 0:
        print("Meson build setup OK")
    else:
        print("Meson build setup failed! ({0} elapsed)")
        sys.exit(1)
    return


def install_project(args):
    """
    Installs the project after building.
    """
    if os.path.exists(PATH_INSTALLED):
        installdir = get_site_packages()
        non_empty = len(os.listdir(PATH_INSTALLED))
        if non_empty and not os.path.exists(installdir):
            raise RuntimeError("Can't install in non-empty directory: "
                               f"'{PATH_INSTALLED}'")
    cmd = ["meson", "install", "-C", args.build_dir]
    log_filename = os.path.join(ROOT_DIR, 'meson-install.log')
    start_time = datetime.datetime.now()
    if args.show_build_log:
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

    if ret != 0:
        if not args.show_build_log:
            with open(log_filename, 'r') as f:
                print(f.read())
        print("Installation failed! ({0} elapsed)".format(elapsed))
        sys.exit(1)

    # ignore everything in the install directory.
    with open(Path(PATH_INSTALLED) / ".gitignore", "w") as f:
        f.write("*")

    print("Installation OK")
    return


def copy_openblas():
    """
    Copies OpenBLAS DLL to the SciPy install dir, and also overwrites the
    default `_distributor_init.py` file with the one we use for wheels uploaded
    to PyPI so that DLL gets loaded.

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
        raise RuntimeError(f'Expecting "lib" at end of "{openblas_lib_path}"')

    # Look in bin subdirectory for OpenBLAS binaries.
    bin_path = openblas_lib_path.parent / 'bin'
    # Locate, make output .libs directory in Scipy install directory.
    scipy_path = Path(get_site_packages()) / 'scipy'
    libs_path = scipy_path / '.libs'
    libs_path.mkdir(exist_ok=True)
    # Copy DLL files from OpenBLAS install to scipy install .libs subdir.
    for dll_fn in bin_path.glob('*.dll'):
        out_fname = libs_path / dll_fn.parts[-1]
        print(f'Copying {dll_fn} to {out_fname}')
        out_fname.write_bytes(dll_fn.read_bytes())

    # Write _distributor_init.py to scipy install dir; this ensures the .libs
    # file is on the DLL search path at run-time, so OpenBLAS gets found
    openblas_support = import_module_from_path(
        'openblas_support',
        Path(ROOT_DIR) / 'tools' / 'openblas_support.py')
    openblas_support.make_init(scipy_path)
    return 0


def get_site_packages():
    """
    Depending on whether we have debian python or not,
    return dist_packages path or site_packages path.
    """
    if 'deb_system' in INSTALL_SCHEMES:
        # Debian patched python in use
        install_cmd = dist.Distribution().get_command_obj('install')
        install_cmd.select_scheme('deb_system')
        install_cmd.finalize_options()
        plat_path = Path(install_cmd.install_platlib)
    else:
        plat_path = Path(get_path('platlib'))
    return str(Path(PATH_INSTALLED) / plat_path.relative_to(sys.exec_prefix))


def build_project(args):
    """
    Build a dev version of the project.

    Returns
    -------
    site_dir
        Directory where the built SciPy version was installed. This is a custom
        prefix, followed by a relative path matching the one the system would
        use for the site-packages of the active Python interpreter.
    """
    root_ok = [os.path.exists(os.path.join(ROOT_DIR, fn))
               for fn in PROJECT_ROOT_FILES]
    if not all(root_ok):
        print("To build the project, run dev.py in "
              "git checkout or unpacked source")
        sys.exit(1)

    env = dict(os.environ)

    setup_build(args, env)

    cmd = ["ninja", "-C", args.build_dir]
    if args.parallel > 1:
        cmd += ["-j", str(args.parallel)]

    # Building with ninja-backend
    ret = subprocess.call(cmd, env=env, cwd=ROOT_DIR)

    if ret == 0:
        print("Build OK")
    else:
        print("Build failed!")
        sys.exit(1)

    install_project(args)

    site_dir = get_site_packages()

    if args.win_cp_openblas and platform.system() == 'Windows':
        if copy_openblas() == 0:
            print('OpenBLAS copied')
        else:
            print("OpenBLAS copy failed!")
            sys.exit(1)

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
