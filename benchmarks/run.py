#!/usr/bin/env python
"""
run.py [options] ASV_COMMAND..

Convenience wrapper around the ``asv`` command; just sets environment
variables and chdirs to the correct place etc.

"""
from __future__ import division, absolute_import, print_function

import os
import sys
import subprocess
import json
import shutil
import argparse
import sysconfig
import errno


EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
              '/usr/local/lib/ccache', '/usr/local/lib/f90cache']

from benchmarks.common import set_mem_rlimit


def main():
    class ASVHelpAction(argparse.Action):
        nargs = 0
        def __call__(self, parser, namespace, values, option_string=None):
            sys.exit(run_asv(['--help']))

    p = argparse.ArgumentParser(usage=__doc__.strip())
    p.add_argument('--help-asv', nargs=0, action=ASVHelpAction,
        help="""show ASV help""")
    p.add_argument('asv_command', nargs=argparse.REMAINDER)
    args = p.parse_args()

    sys.exit(run_asv(args.asv_command))


def run_asv(args):
    cwd = os.path.abspath(os.path.dirname(__file__))

    repo_dir = os.path.join(cwd, 'scipy')

    cmd = ['asv'] + list(args)
    env = dict(os.environ)

    # Inject ccache/f90cache paths
    if sys.platform.startswith('linux'):
        env['PATH'] = os.pathsep.join(EXTRA_PATH + env.get('PATH', '').split(os.pathsep))

    # Control BLAS and CFLAGS
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['CFLAGS'] = drop_bad_flags(sysconfig.get_config_var('CFLAGS'))

    # Limit memory usage
    try:
        set_mem_rlimit()
    except (ImportError, RuntimeError):
        pass

    # Check scipy version if in dev mode; otherwise clone and setup results
    # repository
    if args and (args[0] == 'dev' or '--python=same' in args):
        import scipy
        print("Running benchmarks for Scipy version %s at %s" % (scipy.__version__, scipy.__file__))

    # Override gh-pages
    if 'gh-pages' in args:
        print("gh-pages command is disabled")
        return 1

    # Run
    try:
        return subprocess.call(cmd, env=env, cwd=cwd)
    except OSError as err:
        if err.errno == errno.ENOENT:
            print("Error when running '%s': %s\n" % (" ".join(cmd), str(err),))
            print("You need to install Airspeed Velocity https://spacetelescope.github.io/asv/")
            print("to run Scipy benchmarks")
            return 1
        raise

def drop_bad_flags(flags):
    """
    Drop flags that are problematic for compiling old scipy versions
    """
    if not flags:
        return flags
    return " ".join(x for x in flags.split()
                    if not (x.startswith("-Werror")
                            or x in ("-pedantic-errors",)))


if __name__ == "__main__":
    sys.exit(main())
