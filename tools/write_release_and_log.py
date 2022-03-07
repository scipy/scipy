"""
Standalone script for writing release doc and logs
-------------------------------------------------
Cmd = python write_release_and_log.py

Note:
-----
Work in progress (migrating away from paver implementation)
"""

import os
import sys
import subprocess
import shutil
import warnings
from hashlib import md5
from hashlib import sha256
from pathlib import Path


sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(1, os.path.join(os.path.dirname(__file__), 'tools'))
try:
    version_utils = __import__("version_utils")
    FULLVERSION = version_utils.VERSION
    # This is duplicated from setup.py
    if os.path.exists('.git'):
        GIT_REVISION, _ = version_utils.git_version(
                os.path.join(os.path.dirname(__file__), '..'))
    else:
        GIT_REVISION = "Unknown"

    if not version_utils.ISRELEASED:
        if GIT_REVISION == "Unknown":
            FULLVERSION += '.dev0+Unknown'
        else:
            FULLVERSION += '.dev0+' + GIT_REVISION[:7]
finally:
    sys.path.pop(1)
    sys.path.pop(0)

try:
    # Ensure sensible file permissions
    os.umask(0o022)
except AttributeError:
    # No umask on non-posix
    pass


#-----------------------------------
# Things to be changed for a release
#-----------------------------------

# Source of the release notes
RELEASE = 'doc/release/1.8.0-notes.rst'

# Start/end of the log (from git)
LOG_START = 'v1.8.0'
LOG_END = 'main'


#-------------------------------------------------------
# Hardcoded build/install dirs, virtualenv options, etc.
#-------------------------------------------------------

# Default Python version
PYVER="3.9"


def tarball_name(type_name='gztar'):
    root = 'scipy-%s' % FULLVERSION
    if type_name == 'gztar':
        return root + '.tar.gz'
    elif type_name == 'xztar':
        return root + '.tar.xz'
    elif type_name == 'tar':
        return root + '.tar'
    elif type_name == 'zip':
        return root + '.zip'
    raise ValueError("Unknown type %s" % type_name)


def sdist():
    # First, clean the repo and update submodules (for up-to-date doc html theme
    # and Sphinx extensions)
    os.system('git clean -xdf')
    os.system('git submodule init')
    os.system('git submodule update')

    # Fix file permissions
    os.system('chmod -R a+rX *')

    # To be sure to bypass paver when building sdist... paver + scipy.distutils
    # do not play well together.
    # Cython is run over all Cython files in setup.py, so generated C files
    # will be included.
    os.system('python setup.py sdist --formats=gztar,zip')
    os.system('python setup.py sdist --formats=tar')
    if os.path.exists(os.path.join('dist', tarball_name("xztar"))):
        os.unlink(os.path.join('dist', tarball_name("xztar")))
    # var_path = os.path.join('dist', tarball_name("tar"))
    os.system('xz %s' % os.path.join('dist', tarball_name("tar")))

    # Copy the sdists into installers dir
    if not os.path.exists(os.path.join("release", "installers")):
        os.makedirs(os.path.join("release", "installers"))

    if not os.path.exists(os.path.join('dist', tarball_name("xztar"))):
        warnings.warn("Could not create tar.xz! Do you have xz installed?")
    else:
        t = 'xztar'
        source = os.path.join('dist', tarball_name(t))
        target = os.path.join(os.path.join("release", "installers"), tarball_name(t))
        shutil.copy(source, target)

    for t in ['gztar', 'zip']:
        source = os.path.join('dist', tarball_name(t))
        target = os.path.join(os.path.join("release", "installers"), tarball_name(t))
        shutil.copy(source, target)


#----------------------------
# Release notes and Changelog
#----------------------------

def compute_md5(idirs):
    released = os.listdir(idirs)
    print(released)
    checksums = []
    for fn in sorted(released):
        with open(fn, 'rb') as f:
            m = md5(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums


def compute_sha256(idirs):
    # better checksum so gpg signed README.txt containing the sums can be used
    # to verify the binaries instead of signing all binaries
    released = os.listdir(idirs)
    checksums = []
    for fn in sorted(released):
        with open(fn, 'rb') as f:
            m = sha256(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums


def write_release_task(filename='NOTES.txt'):
    idirs = os.path.join("release", "installers")
    source = Path(RELEASE)
    target = Path(filename)
    if target.exists():
        target.remove()

    # set the file as .rst/.tmp
    tmp_target = Path(filename + '.tmp')
    os.system(f'cp {source} {tmp_target}')

    with open(str(tmp_target), 'a') as ftarget:
        ftarget.writelines("""
Checksums
=========

MD5
~~~

""")
        ftarget.writelines(['%s\n' % c for c in compute_md5(idirs)])
        ftarget.writelines("""
SHA256
~~~~~~

""")
        ftarget.writelines(['%s\n' % c for c in compute_sha256(idirs)])


def write_log_task(filename='Changelog'):
    st = subprocess.Popen(
            ['git', 'log',  '%s..%s' % (LOG_START, LOG_END)],
            stdout=subprocess.PIPE)

    out = st.communicate()[0].decode()
    with open(filename, 'w') as a:
        a.writelines(out)


def release_intermediate():
    write_release_task(os.path.join(os.path.join("release", "installers"), 'README'))
    write_log_task(os.path.join(os.path.join("release", "installers"), 'Changelog'))


def main():
    if not os.path.exists('release') or not os.path.exists(os.path.join('release', 'installers')):
        sdist()
        release_intermediate()
    else:
        release_intermediate()


if __name__ == '__main__':
    main()
