#!/usr/bin/env python3
import os
import subprocess
import argparse


MAJOR = 1
MINOR = 15
MICRO = 0
ISRELEASED = False
IS_RELEASE_BRANCH = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


def get_version_info(source_root):
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists(os.path.join(source_root, '.git')):
        GIT_REVISION, COMMIT_COUNT = git_version(source_root)
    elif os.path.exists('scipy/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load scipy/__init__.py
        import runpy
        ns = runpy.run_path('scipy/version.py')
        GIT_REVISION = ns['git_revision']
        COMMIT_COUNT = ns['git_revision']
    else:
        GIT_REVISION = "Unknown"
        COMMIT_COUNT = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + COMMIT_COUNT + '.' + GIT_REVISION

    return FULLVERSION, GIT_REVISION, COMMIT_COUNT


def write_version_py(source_root, filename='scipy/version.py'):
    cnt = """\
# THIS FILE IS GENERATED DURING THE SCIPY BUILD
# See tools/version_utils.py for details

short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
commit_count = '%(commit_count)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION, COMMIT_COUNT = get_version_info(source_root)

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'commit_count': COMMIT_COUNT,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


# Return the git revision as a string
def git_version(cwd):
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               env=env, cwd=cwd).communicate()[0]
        return out

    try:
        git_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        git_dir = os.path.join(git_dir, ".git")
        out = _minimal_ext_cmd(['git',
                                '--git-dir',
                                git_dir,
                                'rev-parse',
                                'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')[:7]

        # We need a version number that's regularly incrementing for newer commits,
        # so the sort order in a wheelhouse of nightly builds is correct (see
        # https://github.com/MacPython/scipy-wheels/issues/114). It should also be
        # a reproducible version number, so don't rely on date/time but base it on
        # commit history. This gives the commit count since the previous branch
        # point from the current branch (assuming a full `git clone`, it may be
        # less if `--depth` was used - commonly the default in CI):
        prev_version_tag = f'^v{MAJOR}.{MINOR - 2}.0'
        out = _minimal_ext_cmd(['git', '--git-dir', git_dir,
                                'rev-list', 'HEAD', prev_version_tag,
                                '--count'])
        COMMIT_COUNT = out.strip().decode('ascii')
        COMMIT_COUNT = '0' if not COMMIT_COUNT else COMMIT_COUNT
    except OSError:
        GIT_REVISION = "Unknown"
        COMMIT_COUNT = "Unknown"

    return GIT_REVISION, COMMIT_COUNT


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=str, default='.',
                        help="Relative path to the root of the source directory")
    args = parser.parse_args()

    write_version_py(args.source_root)
