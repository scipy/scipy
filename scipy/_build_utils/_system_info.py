#!/bin/env python
"""
Disable LAPACK library from OSX Accelerate, as it causes problems with
ARPACK.

The problems apparently come from the LAPACK portion of the library
--- the BLAS parts are apparently OK. However, because -Wl,-framework
enables both, we need to disable that. Accelerate's BLAS part is still
detected via /usr/lib/libblas.dylib, but we refuse to use the LAPACK
portion.

"""
from __future__ import division, absolute_import, print_function

import sys
import os
import re
import warnings

import numpy.distutils.system_info

from numpy.distutils.system_info import (
    get_info, dict_append,
    AtlasNotFoundError, LapackNotFoundError, LapackSrcNotFoundError,
    BlasNotFoundError, BlasSrcNotFoundError,
    lapack_opt_info as _lapack_opt_info,
    blas_opt_info as _blas_opt_info,
    lapack_info as _lapack_info,
    blas_info as _blas_info)

def monkeypatch_numpy_distutils():
    numpy.distutils.system_info.blas_opt_info = blas_opt_info
    numpy.distutils.system_info.lapack_opt_info = lapack_opt_info

    numpy.distutils.system_info.blas_info = blas_info
    numpy.distutils.system_info.lapack_info = lapack_info


ACCELERATE_WARNING = """
******************************************************************************
NOTE: Accelerate library detected --- it will not be used.

Accelerate on OSX 10.7+ is known to cause invalid results when used
with ARPACK: http://forge.scilab.org/index.php/p/arpack-ng/issues/1259/

Try setting LAPACK and BLAS environment variables to point to a
directory that contains *.a or *.dylib files of a non-Accelerate
LAPACK/BLAS library.
******************************************************************************
"""

class lapack_opt_info(_lapack_opt_info):
    def calc_info(self):
        if sys.platform == 'darwin':
            # Don't look for Accelerate
            pass

        lapack_mkl_info = get_info('lapack_mkl')
        if lapack_mkl_info:
            self.set_info(**lapack_mkl_info)
            return

        atlas_info = get_info('atlas_threads')
        if not atlas_info:
            atlas_info = get_info('atlas')
        #atlas_info = {} ## uncomment for testing
        need_lapack = 0
        need_blas = 0
        info = {}
        if atlas_info:
            l = atlas_info.get('define_macros', [])
            if ('ATLAS_WITH_LAPACK_ATLAS', None) in l \
                   or ('ATLAS_WITHOUT_LAPACK', None) in l:
                need_lapack = 1
            info = atlas_info

        else:
            warnings.warn(AtlasNotFoundError.__doc__)
            need_blas = 1
            need_lapack = 1
            dict_append(info, define_macros=[('NO_ATLAS_INFO', 1)])

        if need_lapack:
            lapack_info = get_info('lapack')
            #lapack_info = {} ## uncomment for testing
            if lapack_info:
                dict_append(info, **lapack_info)
            else:
                warnings.warn(LapackNotFoundError.__doc__)
                lapack_src_info = get_info('lapack_src')
                if not lapack_src_info:
                    warnings.warn(LapackSrcNotFoundError.__doc__)
                    return
                dict_append(info, libraries=[('flapack_src', lapack_src_info)])

        if need_blas:
            blas_info = get_info('blas')
            #blas_info = {} ## uncomment for testing
            if blas_info:
                dict_append(info, **blas_info)
            else:
                warnings.warn(BlasNotFoundError.__doc__)
                blas_src_info = get_info('blas_src')
                if not blas_src_info:
                    warnings.warn(BlasSrcNotFoundError.__doc__)
                    return
                dict_append(info, libraries=[('fblas_src', blas_src_info)])

        self.set_info(**info)
        return


class blas_opt_info(_blas_opt_info):
    def calc_info(self):
        if sys.platform == 'darwin':
            # Don't look for Accelerate
            pass

        blas_mkl_info = get_info('blas_mkl')
        if blas_mkl_info:
            self.set_info(**blas_mkl_info)
            return

        atlas_info = get_info('atlas_blas_threads')
        if not atlas_info:
            atlas_info = get_info('atlas_blas')
        need_blas = 0
        info = {}
        if atlas_info:
            info = atlas_info
        else:
            warnings.warn(AtlasNotFoundError.__doc__)
            need_blas = 1
            dict_append(info, define_macros=[('NO_ATLAS_INFO', 1)])

        if need_blas:
            blas_info = get_info('blas')
            if blas_info:
                dict_append(info, **blas_info)
            else:
                warnings.warn(BlasNotFoundError.__doc__)
                blas_src_info = get_info('blas_src')
                if not blas_src_info:
                    warnings.warn(BlasSrcNotFoundError.__doc__)
                    return
                dict_append(info, libraries=[('fblas_src', blas_src_info)])

        self.set_info(**info)
        return


class blas_info(_blas_info):
    def calc_info(self):
        lib_dirs = self.get_lib_dirs()

        blas_libs = self.get_libs('blas_libs', self._lib_names)
        info = self.check_libs(lib_dirs, blas_libs, [])

        if is_link_to_accelerate(info):
            # Using Accelerate BLAS component...
            dict_append(info, define_macros=[('HAVE_ACCELERATE', 1)])

        info['language'] = 'f77'  # XXX: is it generally true?
        
        self.set_info(**info)


class lapack_info(_lapack_info):
    def calc_info(self):
        lib_dirs = self.get_lib_dirs()

        lapack_libs = self.get_libs('lapack_libs', self._lib_names)

        if is_link_to_accelerate(dict(libraries=lapack_libs,
                                      library_dirs=lib_dirs)):
            # Don't use Accelerate by default
            warnings.warn(ACCELERATE_WARNING)
            lib_dirs = [x for x in lib_dirs if x != '/usr/lib']

        info = self.check_libs(lib_dirs, lapack_libs, [])
        if info is None:
            return
        info['language'] = 'f77'

        self.set_info(**info)


def is_link_to_accelerate(info):
    if sys.platform != 'darwin':
        return False

    if ('-Wl,Accelerate' in info.get('extra_link_args', [])
            or '-Wl,vecLib' in info.get('extra_link_args', [])):
        return True

    path = None
    if info.get('libraries') == ['blas'] or info.get('libraries') == ['lapack']:
        for libdir in info.get('library_dirs', []):
            path = os.path.join(libdir, 'lib' + info.get('libraries')[0] + '.dylib')
            if os.path.exists(path):
                break
            path = os.path.join(libdir, 'lib' + info.get('libraries')[0] + '.a')
            if os.path.exists(path):
                break

    if not path:
        return False

    if os.path.islink(path):
        path = os.readlink(path)

    if 'Accelerate.framework' in path or 'vecLib.framework' in path:
        return True
    else:
        return False
