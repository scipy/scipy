"""
Standalone script for writing release doc and logs
-------------------------------------------------
Cmd -> python write_release_and_log.py <LOG_START> <LOG_END>
Example -> python write_release_and_log.py v1.6.0 v1.8.0
"""

import os
import sys
import glob
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


# ----------------------------
# Release notes and Changelog
# ----------------------------

def compute_md5(idirs):
    released = os.listdir(idirs)
    checksums = []
    for fn in sorted(released):
        fn_updated = os.path.join("release", fn)
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
        fn_updated = os.path.join("release", fn)
        with open(fn_updated, 'rb') as f:
            m = sha256(f.read())
        checksums.append('%s  %s' % (m.hexdigest(), os.path.basename(fn)))

    return checksums


def write_release_task(filename='NOTES.txt'):
    idirs = Path('release')
    source = Path(RELEASE)
    target = Path(filename)
    if target.exists():
        target.remove()

    tmp_target = Path(filename + '.txt')
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
        if not os.path.exists("release"):
            os.makedirs("release")
        else:
            print('Release directory present, executing release tasks')
        write_release_task(os.path.join("release", 'README'))
        write_log_task(os.path.join("release", 'Changelog'))
        print("Release Logs and Readme generated successfully")
    except:
        print("Something went wrong")


def get_latest_release_doc(path):
    files_path = os.path.join(path, '*')
    files = sorted(
        glob.iglob(files_path), key=os.path.getmtime, reverse=True)
    return files[0]


if __name__ == '__main__':
    if len(sys.argv) == 3:
        LOG_START = str(sys.argv[1])
        LOG_END = str(sys.argv[2])
    else:
        print("invalid number of arguments, please add LOG_START and LOG_END")
    RELEASE = get_latest_release_doc('doc/release')
    main()
