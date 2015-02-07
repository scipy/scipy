#!/usr/bin/env python
"""
Convenience wrapper around the ``asv`` command; just sets environment
variables and chdirs to the correct place.
"""
from __future__ import division, absolute_import, print_function

import os
import sys
import subprocess


EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
              '/usr/local/lib/ccache', '/usr/local/lib/f90cache']


from benchmarks.common import set_mem_rlimit


def main():
    sys.exit(run_asv(*sys.argv[1:]))


def run_asv(*args):
    cmd = ['asv'] + list(args)
    cwd = os.path.abspath(os.path.dirname(__file__))
    env = dict(os.environ)

    # Inject ccache/f90cache paths
    if sys.platform.startswith('linux'):
        env['PATH'] = os.pathsep.join(EXTRA_PATH + env.get('PATH', '').split(os.pathsep))

    # Control BLAS config
    env['ATLAS'] = 'None'
    env['OPENBLAS_NUM_THREADS'] = '1'

    # Limit memory usage
    try:
        set_mem_rlimit()
    except (ImportError, RuntimeError):
        pass

    # Check scipy version if in dev mode
    if 'dev' in args:
        import scipy
        print("Running benchmarks for Scipy version %s at %s" % (scipy.__version__, scipy.__file__))

    # Run
    try:
        return subprocess.call(cmd, env=env, cwd=cwd)
    except OSError as err:
        if err.errno == 2:
            print("Error when running '%s': %s\n" % (" ".join(cmd), str(err),))
            print("You need to install Airspeed Velocity https://spacetelescope.github.io/asv/")
            print("to run Scipy benchmarks")
            return 1
        raise


if __name__ == "__main__":
    sys.exit(main())
