#!/usr/bin/env python3
# -*- python -*-
"""
%prog SUBMODULE...

Hack to pipe submodules of Numpy through 2to3 and build them in-place
one-by-one.

Example usage:

    python3 tools/py3tool.py testing distutils core

This will copy files to _py3k/numpy, add a dummy __init__.py and
version.py on the top level, and copy and 2to3 the files of the three
submodules.

When running py3tool again, only changed files are re-processed, which
makes the test-bugfix cycle faster.

"""
from optparse import OptionParser
import shutil
import os
import sys
import re
import subprocess
import fnmatch

if os.environ.get('USE_2TO3CACHE'):
    import lib2to3cache

BASE = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
TEMP = os.path.normpath(os.path.join(BASE, '_py3k'))

SCRIPT_2TO3 = os.path.join(BASE, 'tools', '2to3.py')

EXTRA_2TO3_FLAGS = {
    '*/setup.py': '-x import',
    #'numpy/core/code_generators/generate_umath.py': '-x import',
    #'numpy/core/code_generators/generate_numpy_api.py': '-x import',
    #'numpy/core/code_generators/generate_ufunc_api.py': '-x import',
    #'numpy/core/defchararray.py': '-x unicode',
    #'numpy/compat/py3k.py': '-x unicode',
    #'numpy/ma/timer_comparison.py': 'skip',
    #'numpy/distutils/system_info.py': '-x reduce',
    #'numpy/f2py/auxfuncs.py': '-x reduce',
    #'numpy/lib/arrayterator.py': '-x reduce',
    #'numpy/lib/tests/test_arrayterator.py': '-x reduce',
    #'numpy/ma/core.py': '-x reduce',
    #'numpy/ma/tests/test_core.py': '-x reduce',
    #'numpy/ma/tests/test_old_ma.py': '-x reduce',
    #'numpy/ma/timer_comparison.py': '-x reduce',
    #'numpy/oldnumeric/ma.py': '-x reduce',
}

def main():
    p = OptionParser(usage=__doc__.strip())
    p.add_option("--clean", "-c", action="store_true",
                 help="clean source directory")
    options, args = p.parse_args()

    if not args:
        p.error('no submodules given')
    else:
        dirs = ['scipy/%s' % x for x in map(os.path.basename, args)]

    # Prepare
    if not os.path.isdir(TEMP):
        os.makedirs(TEMP)

    # Set up dummy files (for building only submodules)
    dummy_files = {
        '__init__.py': 'from scipy.version import version as __version__',
        'version.py': 'version = "0.8.0.dev"'
    }

    for fn, content in dummy_files.items():
        fn = os.path.join(TEMP, 'scipy', fn)
        if not os.path.isfile(fn):
            try:
                os.makedirs(os.path.dirname(fn))
            except OSError:
                pass
            f = open(fn, 'wb+')
            f.write(content.encode('ascii'))
            f.close()

    # Environment
    pp = [os.path.abspath(TEMP)]
    def getenv():
        env = dict(os.environ)
        env.update({'PYTHONPATH': ':'.join(pp)})
        return env

    # Copy
    for d in dirs:
        src = os.path.join(BASE, d)
        dst = os.path.join(TEMP, d)

        # Run 2to3
        sync_2to3(dst=dst,
                  src=src,
                  patchfile=os.path.join(TEMP, os.path.basename(d) + '.patch'),
                  clean=options.clean)

        # Run setup.py, falling back to Pdb post-mortem on exceptions
        setup_py = os.path.join(dst, 'setup.py')
        if os.path.isfile(setup_py):
            code = """\
import pdb, sys, traceback
p = pdb.Pdb()
try:
    import __main__
    __main__.__dict__.update({
        "__name__": "__main__", "__file__": "setup.py",
        "__builtins__": __builtins__})
    fp = open("setup.py", "rb")
    try:
        exec(compile(fp.read(), "setup.py", 'exec'))
    finally:
        fp.close()
except SystemExit:
    raise
except:
    traceback.print_exc()
    t = sys.exc_info()[2]
    p.interaction(None, t)
"""
            ret = subprocess.call([sys.executable, '-c', code,
                                   'build_ext', '-i'],
                                  cwd=dst,
                                  env=getenv())
            if ret != 0:
                raise RuntimeError("Build failed.")

        # Run nosetests
        subprocess.call(['nosetests3', '-v', d], cwd=TEMP)

def custom_mangling(filename):
    import_mangling = [
        os.path.join('cluster', '__init__.py'),
        os.path.join('cluster', 'hierarchy.py'),
        os.path.join('cluster', 'vq.py'),
        os.path.join('fftpack', 'basic.py'),
        os.path.join('fftpack', 'pseudo_diffs.py'),
        os.path.join('integrate', 'odepack.py'),
        os.path.join('integrate', 'quadpack.py'),
        os.path.join('integrate', '_ode.py'),
        os.path.join('interpolate', 'fitpack.py'),
        os.path.join('interpolate', 'fitpack2.py'),
        os.path.join('interpolate', 'interpolate.py'),
        os.path.join('interpolate', 'interpolate_wrapper.py'),
        os.path.join('interpolate', 'ndgriddata.py'),
        os.path.join('io', 'array_import.py'),
        os.path.join('io', '__init__.py'),
        os.path.join('io', 'matlab', 'miobase.py'),
        os.path.join('io', 'matlab', 'mio4.py'),
        os.path.join('io', 'matlab', 'mio5.py'),
        os.path.join('io', 'matlab', 'mio5_params.py'),
        os.path.join('linalg', 'basic.py'),
        os.path.join('linalg', 'decomp.py'),
        os.path.join('linalg', 'lapack.py'),
        os.path.join('linalg', 'flinalg.py'),
        os.path.join('linalg', 'iterative.py'),
        os.path.join('linalg', 'misc.py'),
        os.path.join('lib', 'blas', '__init__.py'),
        os.path.join('lib', 'lapack', '__init__.py'),
        os.path.join('ndimage', 'filters.py'),
        os.path.join('ndimage', 'fourier.py'),
        os.path.join('ndimage', 'interpolation.py'),
        os.path.join('ndimage', 'measurements.py'),
        os.path.join('ndimage', 'morphology.py'),
        os.path.join('optimize', 'minpack.py'),
        os.path.join('optimize', 'zeros.py'),
        os.path.join('optimize', 'lbfgsb.py'),
        os.path.join('optimize', 'cobyla.py'),
        os.path.join('optimize', 'slsqp.py'),
        os.path.join('optimize', 'nnls.py'),
        os.path.join('signal', '__init__.py'),
        os.path.join('signal', 'bsplines.py'),
        os.path.join('signal', 'signaltools.py'),
        os.path.join('signal', 'fir_filter_design.py'),
        os.path.join('special', '__init__.py'),
        os.path.join('special', 'add_newdocs.py'),
        os.path.join('special', 'basic.py'),
        os.path.join('special', 'orthogonal.py'),
        os.path.join('spatial', '__init__.py'),
        os.path.join('spatial', 'distance.py'),
        os.path.join('sparse', 'linalg', 'isolve', 'iterative.py'),
        os.path.join('sparse', 'linalg', 'dsolve', 'linsolve.py'),
        os.path.join('sparse', 'linalg', 'dsolve', 'umfpack', 'umfpack.py'),
        os.path.join('sparse', 'linalg', 'eigen', 'arpack', 'arpack.py'),
        os.path.join('sparse', 'linalg', 'eigen', 'arpack', 'speigs.py'),
        os.path.join('sparse', 'linalg', 'iterative', 'isolve', 'iterative.py'),
        os.path.join('sparse', 'csgraph', '__init__.py'),
        os.path.join('sparse', 'csgraph', '_validation.py'),
        os.path.join('stats', 'stats.py'),
        os.path.join('stats', 'distributions.py'),
        os.path.join('stats', 'morestats.py'),
        os.path.join('stats', 'kde.py'),
        os.path.join('stats', 'mstats_basic.py'),
    ]

    if any(filename.endswith(x) for x in import_mangling):
        print(filename)
        f = open(filename, 'r', encoding='utf-8')
        text = f.read()
        f.close()
        for mod in ['_vq', '_hierarchy_wrap', '_fftpack', 'convolve',
                    '_flinalg', 'fblas', 'flapack', 'cblas', 'clapack',
                    'calc_lwork', '_cephes', 'specfun', 'orthogonal_eval',
                    'lambertw', 'ckdtree', '_distance_wrap', '_logit',
                    '_minpack', '_zeros', '_lbfgsb', '_cobyla', '_slsqp',
                    '_nnls',
                    'sigtools', 'spline', 'spectral',
                    '_fitpack', 'dfitpack', '_interpolate',
                    '_odepack', '_quadpack', 'vode', '_dop',
                    'vonmises_cython', '_rank',
                    'futil', 'mvn',
                    '_nd_image',
                    'numpyio',
                    '_superlu', '_arpack', '_iterative', '_umfpack',
                    'interpnd',
                    'mio_utils', 'mio5_utils', 'streams',
                    '_min_spanning_tree', '_shortest_path', '_tools', '_traversal'
                    ]:
            text = re.sub(r'^(\s*)import %s' % mod,
                          r'\1from . import %s' % mod,
                          text, flags=re.M)
            text = re.sub(r'^(\s*)from %s import' % mod,
                          r'\1from .%s import' % mod,
                          text, flags=re.M)
        #text = text.replace('from matrixlib', 'from .matrixlib')
        f = open(filename, 'w', encoding='utf-8')
        f.write(text)
        f.close()

def walk_sync(dir1, dir2, _seen=None):
    if _seen is None:
        seen = {}
    else:
        seen = _seen

    if not dir1.endswith(os.path.sep):
        dir1 = dir1 + os.path.sep

    # Walk through stuff (which we haven't yet gone through) in dir1
    for root, dirs, files in os.walk(dir1):
        sub = root[len(dir1):]
        if sub in seen:
            dirs = [x for x in dirs if x not in seen[sub][0]]
            files = [x for x in files if x not in seen[sub][1]]
            seen[sub][0].extend(dirs)
            seen[sub][1].extend(files)
        else:
            seen[sub] = (dirs, files)
        if not dirs and not files:
            continue
        yield os.path.join(dir1, sub), os.path.join(dir2, sub), dirs, files

    if _seen is None:
        # Walk through stuff (which we haven't yet gone through) in dir2
        for root2, root1, dirs, files in walk_sync(dir2, dir1, _seen=seen):
            yield root1, root2, dirs, files

def sync_2to3(src, dst, patchfile=None, clean=False):
    import lib2to3.main
    from io import StringIO

    to_convert = []

    for src_dir, dst_dir, dirs, files in walk_sync(src, dst):
        for fn in dirs + files:
            src_fn = os.path.join(src_dir, fn)
            dst_fn = os.path.join(dst_dir, fn)

            # skip temporary etc. files
            if fn.startswith('.#') or fn.endswith('~'):
                continue

            # remove non-existing
            if os.path.exists(dst_fn) and not os.path.exists(src_fn):
                if clean:
                    if os.path.isdir(dst_fn):
                        shutil.rmtree(dst_fn)
                    else:
                        os.unlink(dst_fn)
                continue

            # make directories
            if os.path.isdir(src_fn):
                if not os.path.isdir(dst_fn):
                    os.makedirs(dst_fn)
                continue

            dst_dir = os.path.dirname(dst_fn)
            if os.path.isfile(dst_fn) and not os.path.isdir(dst_dir):
                os.makedirs(dst_dir)

            # don't replace up-to-date files
            try:
                if os.path.isfile(dst_fn) and \
                       os.stat(dst_fn).st_mtime >= os.stat(src_fn).st_mtime:
                    continue
            except OSError:
                pass

            # copy file
            shutil.copyfile(src_fn, dst_fn)

            # add .py files to 2to3 list
            if dst_fn.endswith('.py'):
                to_convert.append((src_fn, dst_fn))

    # run 2to3
    flag_sets = {}
    for fn, dst_fn in to_convert:
        flag = ''
        for pat, opt in EXTRA_2TO3_FLAGS.items():
            if fnmatch.fnmatch(fn, pat):
                flag = opt
                break
        flag_sets.setdefault(flag, []).append(dst_fn)

    if patchfile:
        p = open(patchfile, 'wb+')
    else:
        p = open(os.devnull, 'wb')

    for flags, filenames in flag_sets.items():
        if flags == 'skip':
            continue

        _old_stdout = sys.stdout
        try:
            sys.stdout = StringIO()
            lib2to3.main.main("lib2to3.fixes", ['-w', '-n'] + flags.split()+filenames)
        finally:
            sys.stdout = _old_stdout

    for fn, dst_fn in to_convert:
        # perform custom mangling
        custom_mangling(dst_fn)

    p.close()

if __name__ == "__main__":
    main()
