#!/bin/env python
"""
Disable LAPACK library from OSX Accelerate, as it causes problems with
ARPACK.

The problems apparently come from the LAPACK portion of the library
--- the BLAS parts are apparently OK. However, because -Wl,-framework
enables both, we need to disable that. Accelerate's BLAS part is still
detected via /usr/lib/libblas.dylib, but we refuse to use the LAPACK
portion.

We also add Fortran ABI checks to ensure that the LAPACK and BLAS
libraries are ABI compatible with the compiler.

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
            if check_lapack_fortran_abi('Intel MKL', lapack_mkl_info):
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

        lapack_name = 'ATLAS/LAPACK' if not need_lapack else 'LAPACK'
        if not check_lapack_fortran_abi(lapack_name, info):
            return

        self.set_info(**info)
        return


class blas_opt_info(_blas_opt_info):
    def calc_info(self):
        if sys.platform == 'darwin':
            # Don't look for Accelerate
            pass

        blas_mkl_info = get_info('blas_mkl')
        if blas_mkl_info:
            if check_blas_fortran_abi("Intel MKL", blas_mkl_info):
                self.set_info(**blas_mkl_info)
                return

        atlas_info = get_info('atlas_blas_threads')
        if not atlas_info:
            atlas_info = get_info('atlas_blas')
        if not check_blas_fortran_abi("ATLAS", atlas_info):
            atlas_info = None
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

        blas_name = 'ATLAS' if not need_blas else 'Generic BLAS'
        if not check_blas_fortran_abi(blas_name, info):
            return

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

            if not check_blas_fortran_abi("Accelerate/BLAS", info):
                # But the Fortran flags are wrong, so don't bother;
                # the LAPACK supplied has probably an incompatible
                # ABI, and in any case we don't have ABI wrappers for
                # any of the LAPACK routines that need them.
                lib_dirs = [x for x in lib_dirs if x != '/usr/lib']
                info = self.check_libs(lib_dirs, blas_libs, [])
                if not info:
                    return

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


def check_blas_fortran_abi(blas_name, info):
    """
    Compile and run a test program to check whether the given BLAS
    conforms to the ABI of the Fortran compiler.

    The ABI is checked for the main suspect functions: SDOT, DDOT,
    CDOTU, ZDOTU.

    """

    body = """\
      program main
      external sdot, ddot, cdotu, zdotu
      real sx(1), sy(1), sa, sdot
      double precision dx(1), dy(1), da, ddot
      complex cx(1), cy(1), ca, cdotu
      double complex zx(1), zy(1), za, zdotu

      sx(1) = 1e0
      sy(1) = 2e0
      sa = sdot(1, sx, 1, sy, 1)
      if (sa.ne.sx(1)*sy(1)) stop 1

      dx(1) = 1d0
      dy(1) = 2d0
      da = ddot(1, dx, 1, dy, 1)
      if (da.ne.dx(1)*dy(1)) stop 2

      cx(1) = (1e0, 2e0)
      cy(1) = (3e0, 4e0)
      ca = cdotu(1, cx, 1, cy, 1)
      if (ca.ne.cx(1)*cy(1)) stop 3

      zx(1) = (1d0, 2d0)
      zy(1) = (3d0, 4d0)
      za = zdotu(1, zx, 1, zy, 1)
      if (za.ne.zx(1)*zy(1)) stop 4

      write(*,*) 'OK', 'SUCCESS'
      end
    """

    if not info:
        return False

    is_abi_compatible = check_fortran_run(body, info)

    if not is_abi_compatible:
        import textwrap
        msg = textwrap.dedent("""

        ***********************************************************************
        WARNING:

        BLAS library (%s) detected, but its
        Fortran ABI is incompatible with the selected Fortran compiler.
        It is therefore not used now.

        If you are using GNU Fortran Compiler on OSX, setting
        the environment variable FFLAGS=\"-arch i386 -arch x86_64 -fPIC -ff2c\"
        may fix this issue. If you compiled the LAPACK library yourself,
        make sure you use compatible compiler flags here.
        ***********************************************************************

        """
        % (blas_name,))
        warnings.warn(msg)

    return is_abi_compatible


def check_lapack_fortran_abi(lapack_name, info):
    """
    Compile and run a test program to check whether the given LAPACK
    conforms to the ABI of the Fortran compiler.

    The ABI is checked for a main indicator function: CLADIV.

    """

    body = """\
      program main
      external cladiv
      complex cx, cy, ca, cladiv

      cx = (1e0, 1e0)
      cy = (2e0, 2e0)
      ca = cladiv(cx, cy)
      if (abs(ca-cx/cy).gt.1e-7) stop 1

      write(*,*) 'OK', 'SUCCESS'
      end
    """

    if not info:
        return False

    is_abi_compatible = check_fortran_run(body, info)

    if not is_abi_compatible:
        import textwrap
        msg = textwrap.dedent("""

        ***********************************************************************
        WARNING:

        LAPACK library (%s) detected, but its
        Fortran ABI is incompatible with the selected Fortran compiler.
        It is therefore not used now.

        If you are using GNU Fortran Compiler on OSX, setting
        the environment variable FFLAGS=\"-arch i386 -arch x86_64 -fPIC -ff2c\"
        may fix this issue. If you compiled the LAPACK library yourself,
        make sure you use compatible compiler flags here.
        ***********************************************************************

        """
        % (lapack_name,))
        warnings.warn(msg)

    return is_abi_compatible


def check_fortran_run(body, info):
    """
    Compile and run a test program to check whether it runs a line
    write(*,*) 'OK', 'SUCCESS'

    """

    from numpy.distutils.core import get_distribution
    from numpy.distutils.command.config import config as Config

    if not info:
        return False

    dist = get_distribution(True)
    config = Config(dist)
    options = dist.command_options.get('config')
    if options:
        dist._set_command_options('config', config, options)

    libraries = info.get('libraries', [])
    library_dirs = info.get('library_dirs', [])
    extra_compile_args = info.get('extra_compile_args', [])
    extra_link_args = info.get('extra_link_args', [])

    # The distutils config API does not offer a way to pass
    # extra_*_args to the compiler. Therefore, we monkeypatch the
    # active compiler to inject the arguments. (The Fortran compiler
    # classes originate from numpy.distutils so that we are not
    # monkeypatching another library.)

    def new_compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
        if extra_postargs:
            extra_postargs += extra_compile_args
        else:
            extra_postargs = extra_compile_args
        return old_compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts)

    def new_link(self, target_desc, objects,
                 output_filename, output_dir=None, libraries=None,
                 library_dirs=None, runtime_library_dirs=None,
                 export_symbols=None, debug=0, extra_preargs=None,
                 extra_postargs=None, build_temp=None, target_lang=None):
        if extra_postargs:
            extra_postargs += extra_link_args
        else:
            extra_postargs = extra_link_args
        return old_link(self, target_desc, objects,
                        output_filename, output_dir, libraries,
                        library_dirs, runtime_library_dirs,
                        export_symbols, debug, extra_preargs,
                        extra_postargs, build_temp, target_lang)

    config._check_compiler()

    if config.fcompiler is None:
        # No Fortran compiler, so checking the ABI is not needed.
        return True

    old_compile = config.fcompiler.__class__._compile
    old_link = config.fcompiler.__class__.link
    try:
        config.fcompiler.__class__._compile = new_compile
        config.fcompiler.__class__.link = new_link

        # Run the test program
        exitcode, output = config.get_output(body,
                                             libraries=libraries,
                                             library_dirs=library_dirs,
                                             lang="f77")
    finally:
        config.fcompiler.__class__._compile = old_compile
        config.fcompiler.__class__.link = old_link

    # Note: get_output includes also `body` in the output, so be careful
    # in checking the success status. Also, Fortran program exit codes 
    # are undefined.
    success = output and re.search(r'OK\s*SUCCESS', output, re.S)

    return success
