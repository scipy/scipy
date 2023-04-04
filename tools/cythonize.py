#!/usr/bin/env python
"""cythonize

Cythonize pyx files into C files as needed.

Usage: cythonize [root_dir]

Default [root_dir] is 'scipy'.

The number of parallel Cython processes is controlled by the
environment variable SCIPY_NUM_CYTHONIZE_JOBS. If not set, determined
from the number of CPUs.

Checks pyx files to see if they have been changed relative to their
corresponding C files.  If they have, then runs cython on these files to
recreate the C files.

The script thinks that the pyx files have changed relative to the C files
by comparing hashes stored in a database file.

Simple script to invoke Cython (and Tempita) on all .pyx (.pyx.in)
files; while waiting for a proper build system. Uses file hashes to
figure out if rebuild is needed.

For now, this script should be run by developers when changing Cython files
only, and the resulting C files checked in, so that end-users (and Python-only
developers) do not get the Cython/Tempita dependencies.

Originally written by Dag Sverre Seljebotn, and copied here from:

https://raw.github.com/dagss/private-scipy-refactor/cythonize/cythonize.py

Note: this script does not check any of the dependent C libraries; it only
operates on the Cython .pyx files.

"""

import os
import re
import sys
import hashlib
import subprocess
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool, Lock
from os.path import dirname, join

HASH_FILE = 'cythonize.dat'
DEFAULT_ROOT = 'scipy'

#
# Rules
#
def process_pyx(fromfile, tofile, cwd):
    try:
        from Cython.Compiler.Version import version as cython_version
        from scipy._lib import _pep440

        # Try to find pyproject.toml
        pyproject_toml = join(dirname(__file__), '..', 'pyproject.toml')
        if not os.path.exists(pyproject_toml):
            raise ImportError()

        # Try to find the minimum version from pyproject.toml
        with open(pyproject_toml) as pt:
            for line in pt:
                if "cython" not in line.lower():
                    continue
                line = ''.join(line.split('=')[1:])  # get rid of "Cython>="
                if ',<' in line:
                    # There's an upper bound as well
                    split_on = ',<'
                    if ',<=' in line:
                        split_on = ',<='
                    min_required_version, max_required_version = line.split(split_on)
                    max_required_version, _ = max_required_version.split('"')
                else:
                    min_required_version, _ = line.split('"')

                break
            else:
                raise ImportError()

        # Note: we only check lower bound, for upper bound we rely on pip
        # respecting pyproject.toml. Reason: we want to be able to build/test
        # with more recent Cython locally or on main, upper bound is for
        # sdist in a release.
        if _pep440.parse(cython_version) < _pep440.Version(min_required_version):
            raise Exception('Building SciPy requires Cython >= {}, found '
                            '{}'.format(min_required_version, cython_version))

    except ImportError:
        pass

    flags = ['--fast-fail', '-3']
    if tofile.endswith('.cxx'):
        flags += ['--cplus']

    try:
        try:
            r = subprocess.call(['cython'] + flags + ["-o", tofile, fromfile], cwd=cwd)
            if r != 0:
                raise Exception('Cython failed')
        except OSError as e:
            # There are ways of installing Cython that don't result in a cython
            # executable on the path, see gh-2397.
            r = subprocess.call([sys.executable, '-c',
                                 'import sys; from Cython.Compiler.Main import '
                                 'setuptools_main as main; sys.exit(main())'] + flags +
                                 ["-o", tofile, fromfile],
                                cwd=cwd)
            if r != 0:
                raise Exception("Cython either isn't installed or it failed.") from e
    except OSError as e:
        raise OSError('Cython needs to be installed') from e

def process_tempita_pyx(fromfile, tofile, cwd):
    try:
        try:
            from Cython import Tempita as tempita
        except ImportError:
            import tempita
    except ImportError as e:
        raise Exception('Building SciPy requires Tempita: '
                        'pip install --user Tempita') from e
    with open(os.path.join(cwd, fromfile), mode='r') as f_in:
        template = f_in.read()
        pyxcontent = tempita.sub(template)
        assert fromfile.endswith('.pyx.in')
        pyxfile = fromfile[:-len('.in')]
        with open(os.path.join(cwd, pyxfile), "w", encoding='utf8') as f_out:
            f_out.write(pyxcontent)
    process_pyx(pyxfile, tofile, cwd)


rules = {
    # fromext : function
    '.pyx': process_pyx,
    '.pyx.in': process_tempita_pyx
    }

#
# Hash db
#
def load_hashes(filename):
    # Return { filename : (sha1 of input, sha1 of output) }
    if os.path.isfile(filename):
        hashes = {}
        with open(filename, 'r') as f:
            for line in f:
                filename, inhash, outhash = line.split()
                if outhash == "None":
                    outhash = None
                hashes[filename] = (inhash, outhash)
    else:
        hashes = {}
    return hashes

def save_hashes(hash_db, filename):
    with open(filename, 'w') as f:
        for key, value in sorted(hash_db.items()):
            f.write("%s %s %s\n" % (key, value[0], value[1]))

def sha1_of_file(filename):
    h = hashlib.sha1()
    with open(filename, "rb") as f:
        h.update(f.read())
    return h.hexdigest()

#
# Main program
#

def normpath(path):
    path = path.replace(os.sep, '/')
    if path.startswith('./'):
        path = path[2:]
    return path

def get_hash(frompath, topath):
    from_hash = sha1_of_file(frompath)
    if topath:
        to_hash = sha1_of_file(topath) if os.path.exists(topath) else None
    else:
        to_hash = None
    return (from_hash, to_hash)

def get_cython_dependencies(fullfrompath):
    fullfromdir = os.path.dirname(fullfrompath)
    deps = set()
    with open(fullfrompath, 'r') as f:
        pxipattern = re.compile(r'include "([a-zA-Z0-9_]+\.pxi)"')
        pxdpattern1 = re.compile(r'from \. cimport ([a-zA-Z0-9_]+)')
        pxdpattern2 = re.compile(r'from \.([a-zA-Z0-9_]+) cimport')

        for line in f:
            m = pxipattern.match(line)
            if m:
                deps.add(os.path.join(fullfromdir, m.group(1)))
            m = pxdpattern1.match(line)
            if m:
                deps.add(os.path.join(fullfromdir, m.group(1) + '.pxd'))
            m = pxdpattern2.match(line)
            if m:
                deps.add(os.path.join(fullfromdir, m.group(1) + '.pxd'))
    return list(deps)

def process(path, fromfile, tofile, processor_function, hash_db,
            dep_hashes, lock):
    with lock:
        fullfrompath = os.path.join(path, fromfile)
        fulltopath = os.path.join(path, tofile)
        current_hash = get_hash(fullfrompath, fulltopath)
        if current_hash == hash_db.get(normpath(fullfrompath), None):
            file_changed = False
        else:
            file_changed = True

        deps_changed = False
        deps = get_cython_dependencies(fullfrompath)
        for dep in deps:
            dep_hash = get_hash(dep, None)
            if dep_hash == hash_db.get(normpath(dep), None):
                continue
            else:
                dep_hashes[normpath(dep)] = dep_hash
                deps_changed = True

        if not file_changed and not deps_changed:
            print('%s has not changed' % fullfrompath)
            sys.stdout.flush()
            return

        print('Processing %s' % fullfrompath)
        sys.stdout.flush()

    processor_function(fromfile, tofile, cwd=path)

    with lock:
        # changed target file, recompute hash
        current_hash = get_hash(fullfrompath, fulltopath)
        # store hash in db
        hash_db[normpath(fullfrompath)] = current_hash

def process_generate_pyx(path, lock):
    with lock:
        print('Running {}'.format(path))
    ret = subprocess.call([sys.executable, path])
    with lock:
        if ret != 0:
            raise RuntimeError("Running {} failed".format(path))

def find_process_files(root_dir):
    lock = Lock()

    try:
        num_proc = int(os.environ.get('SCIPY_NUM_CYTHONIZE_JOBS', cpu_count()))
        pool = Pool(processes=num_proc)
    except ImportError as e:
        # Allow building (single-threaded) on GNU/Hurd, which does not
        # support semaphores so Pool cannot initialize.
        pool = type('', (), {'imap_unordered': lambda self, func,
                iterable: map(func, iterable)})()
    except ValueError:
        pool = Pool()

    hash_db = load_hashes(HASH_FILE)
    # Keep changed pxi/pxd hashes in a separate dict until the end
    # because if we update hash_db and multiple files include the same
    # .pxi file the changes won't be detected.
    dep_hashes = {}

    # Run any _generate_pyx.py scripts
    jobs = []
    for cur_dir, dirs, files in os.walk(root_dir):
        generate_pyx = os.path.join(cur_dir, '_generate_pyx.py')
        if os.path.exists(generate_pyx):
            jobs.append(generate_pyx)

    for result in pool.imap_unordered(lambda fn: process_generate_pyx(fn, lock), jobs):
        pass

    # Process pyx files
    jobs = []
    for cur_dir, dirs, files in os.walk(root_dir):
        for filename in files:
            in_file = os.path.join(cur_dir, filename + ".in")
            if filename.endswith('.pyx') and os.path.isfile(in_file):
                continue
            for fromext, function in rules.items():
                if filename.endswith(fromext):
                    toext = ".c"
                    with open(os.path.join(cur_dir, filename), 'rb') as f:
                        data = f.read()
                        m = re.search(br"^\s*#\s*distutils:\s*language\s*=\s*c\+\+\s*$", data, re.I | re.M)
                        if m:
                            toext = ".cxx"
                    fromfile = filename
                    tofile = filename[:-len(fromext)] + toext
                    jobs.append((cur_dir, fromfile, tofile, function,
                                 hash_db, dep_hashes, lock))

    for result in pool.imap_unordered(lambda args: process(*args), jobs):
        pass

    hash_db.update(dep_hashes)
    save_hashes(hash_db, HASH_FILE)

def main():
    try:
        root_dir = sys.argv[1]
    except IndexError:
        root_dir = DEFAULT_ROOT
    find_process_files(root_dir)


if __name__ == '__main__':
    main()
