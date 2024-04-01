"""
Standalone script for writing release doc and logs::

    python tools/write_release_and_log.py <LOG_START> <LOG_END>

Example::

    python tools/write_release_and_log.py v1.7.0 v1.8.0

Needs to be run from the root of the repository.

"""

import os
import sys
import subprocess
from hashlib import md5
from hashlib import sha256
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(1, os.path.join(os.path.dirname(__file__), 'tools'))

try:
    version_utils = __import__("version_utils")
    FULLVERSION = version_utils.VERSION
    # This is duplicated from tools/version_utils.py
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


def get_latest_release_doc(path):
    """
    Method to pick the file from 'doc/release' with the highest
    release number (e.g., `1.9.0-notes.rst`).
    """
    file_paths = os.listdir(path)
    file_paths.sort(key=lambda x: list(map(int, (x.split('-')[0].split('.')))))
    return os.path.join(path, file_paths[-1])


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
        checksums.append(f'{m.hexdigest()}  {os.path.basename(fn)}')
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
        checksums.append(f'{m.hexdigest()}  {os.path.basename(fn)}')

    return checksums


def write_release_task(filename='NOTES.txt'):
    idirs = Path('release')
    source = Path(get_latest_release_doc('doc/source/release'))
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
        print("Release README generated successfully")


def write_log_task(filename='Changelog'):
    st = subprocess.Popen(
        ['git', 'log', f'{LOG_START}..{LOG_END}'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = st.communicate()
    if not st.returncode == 0:
        raise RuntimeError("%s failed" % str(error))

    out = st.communicate()[0].decode()
    with open(filename, 'w') as a:
        a.writelines(out)

    print("Release logs generated successfully")


def main():
    """
    Checks weather release directory is present or not
    and calls the method to generate logs and notes
    """
    if not os.path.exists("release"):
        os.makedirs("release")

    write_release_task(os.path.join("release", 'README'))
    write_log_task(os.path.join("release", 'Changelog'))


if __name__ == '__main__':
    if len(sys.argv) == 3:
        LOG_START = str(sys.argv[1])
        LOG_END = str(sys.argv[2])
    else:
        print("invalid number of arguments, please add LOG_START and LOG_END")

    main()
