"""
Standalone script for writing release doc and logs
-------------------------------------------------
Cmd -> python write_release_and_log.py <LOG_START> <LOG_END> <VERSION_NOTES>
Example -> python write_release_and_log.py v1.6.0 v1.8.0 1.9.0

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


# ----------------------------
# Release notes and Changelog
# ----------------------------

def compute_md5(idirs):
    released = os.listdir(idirs)
    checksums = []
    for fn in sorted(released):
        fn_updated = os.path.join(os.path.join("release", "installers"), fn)
        with open(fn_updated, 'rb') as f:
            m = md5(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))
    return checksums


def compute_sha256(idirs):
    # better checksum so gpg signed README.txt containing the sums can be used
    # to verify the binaries instead of signing all binaries
    released = os.listdir(idirs)
    checksums = []
    for fn in sorted(released):
        fn_updated = os.path.join(os.path.join("release", "installers"), fn)
        with open(fn_updated, 'rb') as f:
            m = sha256(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums


def write_release_task(filename='NOTES.txt'):
    idirs = os.path.join('release', 'installers')
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
        ['git', 'log', '%s..%s' % (LOG_START, LOG_END)],
        stdout=subprocess.PIPE)

    out = st.communicate()[0].decode()
    with open(filename, 'w') as a:
        a.writelines(out)


def main():
    try:
        if not os.path.exists(os.path.join("release", "installers")):
            os.makedirs(os.path.join("release", "installers"))
        else:
            print('Release/installer directory present, executing release tasks')
        write_release_task(os.path.join(os.path.join("release", "installers"), 'README'))
        write_log_task(os.path.join(os.path.join("release", "installers"), 'Changelog'))
        print("Release Logs and Readme generated successfully")
    except:
        print("Something went wrong")


if __name__ == '__main__':
    if len(sys.argv) == 4:
        LOG_START = str(sys.argv[1])
        LOG_END = str(sys.argv[2])
        VERSION_NOTES = str(sys.argv[3])
        RELEASE = f'doc/release/{VERSION_NOTES}-notes.rst'
        print(len(sys.argv))
    else:
        print("invalid number of arguments, pass LOG_START,LOG_END and VERSION_NOTES")
    main()
