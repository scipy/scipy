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


EXTRA_PATH = ['/usr/lib/ccache', '/usr/lib/f90cache',
              '/usr/local/lib/ccache', '/usr/local/lib/f90cache']

RESULT_REPO = "git@github.com:pv/scipy-bench.git"

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
    env['CFLAGS'] = drop_bad_flags(sysconfig.get_config_var('CFLAGS'))

    # Limit memory usage
    try:
        set_mem_rlimit()
    except (ImportError, RuntimeError):
        pass

    # Check scipy version if in dev mode; otherwise clone and setup results
    # repository
    if args and args[0] == 'dev':
        import scipy
        print("Running benchmarks for Scipy version %s at %s" % (scipy.__version__, scipy.__file__))
    else:
        setup_sync()

    # Override gh-pages
    if 'gh-pages' in args:
        return do_gh_pages(env)

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


def setup_sync():
    cwd = os.path.abspath(os.path.dirname(__file__))
    sync_dir = os.path.join(cwd, 'scipy-benchmarks')
    results_dir = os.path.join(cwd, 'results')
    results_dir_bak = os.path.join(cwd, 'results.bak')

    if not is_git_repo_root(sync_dir):
        subprocess.check_call(['git', 'clone', RESULT_REPO, sync_dir])

    if os.path.isdir(results_dir) and not os.path.islink(results_dir):
        os.rename(results_dir, results_dir_bak)

    if not os.path.islink(results_dir):
        os.symlink(os.path.join('scipy-benchmarks', 'results'), results_dir)


def do_gh_pages(env):
    cwd = os.path.abspath(os.path.dirname(__file__))
    sync_dir = os.path.join(cwd, 'scipy-benchmarks')
    tmp_clone = os.path.join(cwd, 'scipy-benchmarks.tmp')

    subprocess.check_call(['git', 'checkout', 'master'], cwd=sync_dir)
    subprocess.check_call(['git', 'pull', '--ff-only'], cwd=sync_dir)

    subprocess.check_call(['asv', 'publish'], env=env, cwd=cwd)

    if os.path.exists(tmp_clone):
        shutil.rmtree(tmp_clone)

    subprocess.check_call(['git', 'clone', '--shared', '-b', 'master',
                           sync_dir, tmp_clone])
    assert is_git_repo_root(tmp_clone)

    def git_tmp(*args):
        subprocess.check_call(['git', '-C', tmp_clone] + list(args))

    def git_tmp_failok(*args):
        subprocess.call(['git', '-C', tmp_clone] + list(args))

    git_tmp('remote', 'add', 'upstream', RESULT_REPO)
    git_tmp_failok('branch', '-D', 'gh-pages')
    git_tmp('checkout', '--orphan', 'gh-pages')

    html_root = os.path.join(cwd, 'html')
    for fn in os.listdir(html_root):
        src = os.path.join(html_root, fn)
        dst = os.path.join(tmp_clone, fn)
        if os.path.isdir(src):
            shutil.copytree(src, dst)
        else:
            shutil.copyfile(src, dst)

    with open(os.path.join(tmp_clone, '.nojekyll'), 'wb') as fd:
        fd.write(b'\n')

    git_tmp('add', '-f', '.')
    git_tmp('commit', '-m', 'Generated from sources')
    git_tmp('push', '-f', 'upstream', 'gh-pages')

    shutil.rmtree(tmp_clone)

    return 0


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
