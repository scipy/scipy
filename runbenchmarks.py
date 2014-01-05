#!/usr/bin/env python
"""
runbenchmarks.py [OPTIONS] [-- ARGS]

Run benchmarks, building the project first.
For development build options unrelated to benchmarking, see runtests.py.

Examples::

    $ python runbenchmarks.py
    $ python runbenchmarks.py -s {SAMPLE_SUBMODULE}
    $ python runbenchmarks.py -t {SAMPLE_BENCHMARK}

"""

#
# This is a generic benchmark runner script.
# Change the following values to adapt to your project:
#

PROJECT_MODULE = "scipy"
PROJECT_ROOT_FILES = ['scipy', 'LICENSE.txt', 'setup.py']
SAMPLE_BENCHMARK = "scipy/linalg/benchmarks/bench_decom.py:bench_eigvals"
SAMPLE_SUBMODULE = "linalg"

#TODO this is used only for building
#     so it could possibly be moved out of this script (?)
EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
              '/usr/local/lib/ccache', '/usr/local/lib/f90cache']

# ---------------------------------------------------------------------


if __doc__ is None:
    __doc__ = "Run without -OO if you want usage info"
else:
    __doc__ = __doc__.format(**globals())


import sys
import os

# In case we are run from the source directory, we don't want to import the
# project from there:
sys.path.pop(0)

#TODO these three imports are used only for building,
#     so they could possibly be moved out of this script (?)
import time
import subprocess
import imp

import shutil
from argparse import ArgumentParser, REMAINDER

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("--verbose", "-v", action="count", default=1,
                        help="more verbosity")
    parser.add_argument("--no-build", "-n", action="store_true", default=False,
                        help="do not build the project (use system installed version)")
    parser.add_argument("--mode", "-m", default="fast",
                        help="'fast', 'full', or something that could be "
                             "passed to nosetests -A [default: fast]")
    parser.add_argument("--submodule", "-s", default=None,
                        help="Submodule whose benchmarks to run (cluster, constants, ...)")
    parser.add_argument("--pythonpath", "-p", default=None,
                        help="Paths to prepend to PYTHONPATH")
    parser.add_argument("--benchmarks", "-b", action='append',
                        help="Specify benchmarks to run")
    parser.add_argument("--debug", "-g", action="store_true",
                        help="Debug build")
    parser.add_argument("--show-build-log", action="store_true",
                        help="Show build output rather than using a log file")
    parser.add_argument("args", metavar="ARGS", default=[], nargs=REMAINDER,
                        help="Arguments to pass to Nose, Python or shell")
    args = parser.parse_args(argv)

    if args.pythonpath:
        for p in reversed(args.pythonpath.split(os.pathsep)):
            sys.path.insert(0, p)

    if not args.no_build:
        site_dir = build_project(ROOT_DIR, PROJECT_ROOT_FILES,
                debug=args.debug, show_build_log=args.show_build_log)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = site_dir

    extra_argv = args.args[:]
    if extra_argv and extra_argv[0] == '--':
        extra_argv = extra_argv[1:]

    bench_dir = os.path.join(ROOT_DIR, 'build', 'bench')

    if args.submodule:
        modname = PROJECT_MODULE + '.' + args.submodule
        try:
            __import__(modname)
            bench = sys.modules[modname].bench
        except (ImportError, KeyError, AttributeError) as e:
            print("Cannot run benchmarks for %s (%s)" % (modname, e))
            sys.exit(2)
    elif args.benchmarks:
        def fix_bench_path(x):
            # fix up bench path
            p = x.split(':')
            p[0] = os.path.relpath(os.path.abspath(p[0]), bench_dir)
            return ':'.join(p)

        benchmarks = [fix_bench_path(x) for x in args.benchmarks]

        def bench(*a, **kw):
            extra_argv = kw.pop('extra_argv', ())
            extra_argv = extra_argv + benchmarks[1:]
            kw['extra_argv'] = extra_argv
            from numpy.testing import Tester
            return Tester(benchmarks[0]).bench(*a, **kw)
    else:
        __import__(PROJECT_MODULE)
        bench = sys.modules[PROJECT_MODULE].bench

    # Run the benchmarks under build/bench
    try:
        shutil.rmtree(bench_dir)
    except OSError:
        pass
    try:
        os.makedirs(bench_dir)
    except OSError:
        pass

    success = False
    cwd = os.getcwd()
    try:
        os.chdir(bench_dir)
        success = bench(args.mode,
                      verbose=args.verbose,
                      extra_argv=extra_argv)
    finally:
        os.chdir(cwd)

    if success:
        sys.exit(0)
    else:
        sys.exit(1)


# TODO this function could be moved out of this script (?)
def build_project(root_dir, project_root_files,
        debug=None, show_build_log=None, gcov=None):
    """
    Build a dev version of the project.

    Returns
    -------
    site_dir
        site-packages directory where it was installed

    """

    root_ok = [os.path.exists(os.path.join(root_dir, fn))
               for fn in project_root_files]
    if not all(root_ok):
        print("To build the project, run this script in "
              "git checkout or unpacked source")
        sys.exit(1)

    dst_dir = os.path.join(root_dir, 'build', 'testenv')

    env = dict(os.environ)
    cmd = [sys.executable, 'setup.py']

    # Always use ccache, if installed
    env['PATH'] = os.pathsep.join(
            EXTRA_PATH + env.get('PATH', '').split(os.pathsep))

    if debug or gcov:
        # assume everyone uses gcc/gfortran
        env['OPT'] = '-O0 -ggdb'
        env['FOPT'] = '-O0 -ggdb'
        if gcov:
            import distutils.sysconfig
            cvars = distutils.sysconfig.get_config_vars()
            env['OPT'] = '-O0 -ggdb'
            env['FOPT'] = '-O0 -ggdb'
            env['CC'] = cvars['CC'] + ' --coverage'
            env['CXX'] = cvars['CXX'] + ' --coverage'
            env['F77'] = 'gfortran --coverage '
            env['F90'] = 'gfortran --coverage '
            env['LDSHARED'] = cvars['LDSHARED'] + ' --coverage'
            env['LDFLAGS'] = " ".join(
                    cvars['LDSHARED'].split()[1:]) + ' --coverage'
        cmd += ["build"]

    cmd += ['install', '--prefix=' + dst_dir]

    log_filename = os.path.join(root_dir, 'build.log')

    if show_build_log:
        ret = subprocess.call(cmd, env=env, cwd=root_dir)
    else:
        log_filename = os.path.join(root_dir, 'build.log')
        print("Building, see build.log...")
        with open(log_filename, 'w') as log:
            p = subprocess.Popen(cmd, env=env, stdout=log, stderr=log,
                                 cwd=root_dir)

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
                    print("    ... build in progress")
                    last_blip = time.time()
                    last_log_size = log_size

        ret = p.wait()

    if ret == 0:
        print("Build OK")
    else:
        if not show_build_log:
            with open(log_filename, 'r') as f:
                print(f.read())
            print("Build failed!")
        sys.exit(1)

    from distutils.sysconfig import get_python_lib
    site_dir = get_python_lib(prefix=dst_dir, plat_specific=True)

    return site_dir


#
# Python 3 support
#

if sys.version_info[0] >= 3:
    import builtins
    exec_ = getattr(builtins, "exec")
else:
    def exec_(code, globs=None, locs=None):
        """Execute code in a namespace."""
        if globs is None:
            frame = sys._getframe(1)
            globs = frame.f_globals
            if locs is None:
                locs = frame.f_locals
            del frame
        elif locs is None:
            locs = globs
        exec("""exec code in globs, locs""")

if __name__ == "__main__":
    main(argv=sys.argv[1:])
