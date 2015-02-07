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
import argparse


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
    p.add_argument("--current-repo", action="store_true",
        help="""use current repository as the upstream repository,
        rather than cloning it from the internet; enables running
        benchmarks on e.g. your own branches""")
    p.add_argument('asv_command', nargs=argparse.REMAINDER)
    args = p.parse_args()

    sys.exit(run_asv(args.asv_command, current_repo=args.current_repo))


def run_asv(args, current_repo=False):
    cwd = os.path.abspath(os.path.dirname(__file__))

    if current_repo:
        try:
            from asv.util import load_json, write_json
            conf = load_json(os.path.join(cwd, 'asv.conf.json'))
            conf['repo'] = os.path.normpath(os.path.join(cwd, '..'))
            cfg_fn = os.path.join(cwd, '.asvconf.tmp')
            write_json(cfg_fn, conf)
            args = ['--config', cfg_fn] + args
        except ImportError:
            pass

    repo_dir = os.path.join(cwd, 'scipy')
    if is_git_repo_root(repo_dir):
        if current_repo:
            url = os.path.normpath(os.path.join(cwd, '..'))
        else:
            url = "https://github.com/scipy/scipy.git"
        subprocess.call(['git', 'remote', 'set-url', "origin", url],
                        cwd=repo_dir)

    cmd = ['asv'] + list(args)
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
    if args and args[0] == 'dev':
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


def is_git_repo_root(path):
    try:
        p = subprocess.Popen(['git', '-C', path, 'rev-parse', '--git-dir'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode != 0:
            return False
        return (out.strip() == '.git')
    except OSError:
        return False


if __name__ == "__main__":
    sys.exit(main())
