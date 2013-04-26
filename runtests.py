#!/usr/bin/env python
"""
runtests.py [OPTIONS] [-- ARGS]

Run Scipy tests, building Scipy first.

Examples::

    $ python runtests.py
    $ python runtests.py -s optimize
    $ python runtests.py -t scipy/special/tests/test_basic.py:test_xlogy
    $ python runtests.py --ipython

"""

import sys
import os

# In case we are run from the source directory, we don't want to import scipy
# from there, we want to import the installed version:
sys.path.pop(0)

import shutil
import subprocess
from argparse import ArgumentParser, REMAINDER

def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("--quieter", "-q", action="count", default=0,
                        help="less verbosity")
    parser.add_argument("--no-build", "-n", action="store_true", default=False,
                        help="do not build Scipy (use system installed version)")
    parser.add_argument("--build-only", "-b", action="store_true", default=False,
                        help="just build, do not run any tests")
    parser.add_argument("--doctests", action="store_true", default=False,
                        help="Run doctests in module")
    parser.add_argument("--coverage", action="store_true", default=False,
                        help=("report coverage of Scipy code. HTML output goes "
                              "under build/coverage"))
    parser.add_argument("--mode", "-m", default="fast",
                        help="'fast', 'full', or something that could be "
                             "passed to nosetests -A [default: fast]")
    parser.add_argument("--submodule", "-s", default=None,
                        help="Submodule whose tests to run (cluster, constants, ...)")
    parser.add_argument("--pythonpath", "-p", default=None,
                        help="Paths to prepend to PYTHONPATH")
    parser.add_argument("--tests", "-t", default=[], nargs='*',
                        help="Specify tests to run")
    parser.add_argument("--python", action="store_true",
                        help="Start a Python shell with PYTHONPATH set")
    parser.add_argument("--ipython", "-i", action="store_true",
                        help="Start IPython shell with PYTHONPATH set")
    parser.add_argument("--shell", action="store_true",
                        help="Start Unix shell with PYTHONPATH set")
    parser.add_argument("--debug", "-g", action="store_true",
                        help="Debug build")
    parser.add_argument("args", metavar="ARGS", default=[], nargs=REMAINDER,
                        help="Arguments to pass to Nose")
    args = parser.parse_args(argv)

    if args.pythonpath:
        for p in reversed(args.pythonpath.split(os.pathsep)):
            sys.path.insert(0, p)

    if not args.no_build:
        site_dir = build_scipy(args)
        sys.path.insert(0, site_dir)
        os.environ['PYTHONPATH'] = site_dir

    if args.python:
        import code
        code.interact()
        sys.exit(0)

    if args.ipython:
        import IPython
        IPython.embed()
        sys.exit(0)

    if args.shell:
        shell = os.environ.get('SHELL', 'sh')
        print("Spawning a Unix shell...")
        os.execv(shell, [shell])
        sys.exit(1)

    extra_argv = args.args

    if args.coverage:
        dst_dir = os.path.join('build', 'coverage')
        fn = os.path.join(dst_dir, 'coverage_html.js')
        if os.path.isdir(dst_dir) and os.path.isfile(fn):
            shutil.rmtree(dst_dir)
        extra_argv += ['--cover-html', 
                       '--cover-html-dir='+dst_dir]

    if args.build_only:
        sys.exit(0)
    elif args.submodule:
        modname = 'scipy.' + args.submodule
        try:
            __import__(modname)
            test = sys.modules[modname].test
        except (ImportError, KeyError, AttributeError):
            print("Cannot run tests for %s" % modname)
            sys.exit(2)
    elif args.tests:
        import nose
        sys.argv = ['nosetests', '--exe']
        if args.mode == "fast":
            sys.argv += ['-A', "not slow"]
        if args.doctests:
            sys.argv += ['--doctest-tests']
        if args.quieter <= 1:
            sys.argv += ['-v']
        sys.argv += args.tests
        nose.main(sys.argv)
        sys.exit(0)
    else:
        import scipy
        test = scipy.test

    result = test(args.mode,
                  verbose=max(0, 2-args.quieter),
                  extra_argv=args.args,
                  doctests=args.doctests,
                  coverage=args.coverage)

    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)

def build_scipy(args):
    """
    Build a dev version of Scipy.

    Returns
    -------
    site_dir
        site-packages directory where it was installed

    """

    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    root_ok = [os.path.exists(os.path.join(root_dir, fn))
               for fn in ['HACKING.rst.txt', 'scipy']]
    if not all(root_ok):
        print("To use the -b option, run the scipy/runtests.py in "
              "Scipy's git checkout or unpacked source")
        sys.exit(1)

    dst_dir = os.path.join(root_dir, 'build', 'testenv')

    env = dict(os.environ)
    cmd = [sys.executable, 'setup.py']

    # Always use ccache, if installed
    env['PATH'] = os.pathsep.join(['/usr/lib/ccache'] 
                                  + env.get('PATH', '').split(os.pathsep))

    if args.debug:
        # assume everyone uses gcc/gfortran
        env['OPT'] = '-O0 -ggdb'
        env['FOPT'] = '-O0 -ggdb'
        cmd += ["build", "--debug"]

    cmd += ['install', '--prefix=' + dst_dir]

    print("Building, see build.log...")
    with open('build.log', 'w') as log:
        ret = subprocess.call(cmd, env=env, stdout=log, stderr=log)

    if ret == 0:
        print("Build OK")
    else:
        with open('build.log', 'r') as f:
            print(f.read())
        print("Build failed!")
        sys.exit(1)

    from distutils.sysconfig import get_python_lib
    site_dir = get_python_lib(prefix=dst_dir)

    return site_dir
        
if __name__ == "__main__":
    main(argv=sys.argv[1:])
