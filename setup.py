#!/usr/bin/env python
"""SciPy: Scientific Library for Python

SciPy (pronounced "Sigh Pie") is open-source software for mathematics,
science, and engineering. The SciPy library
depends on NumPy, which provides convenient and fast N-dimensional
array manipulation. The SciPy library is built to work with NumPy
arrays, and provides many user-friendly and efficient numerical
routines such as routines for numerical integration and optimization.
Together, they run on all popular operating systems, are quick to
install, and are free of charge.  NumPy and SciPy are easy to use,
but powerful enough to be depended upon by some of the world's
leading scientists and engineers. If you need to manipulate
numbers on a computer and display or publish the results,
give SciPy a try!

"""

DOCLINES = (__doc__ or '').split("\n")

import os
import sys
import subprocess
import textwrap
import warnings
import sysconfig
from tools.version_utils import write_version_py, get_version_info
from tools.version_utils import IS_RELEASE_BRANCH
import importlib


if sys.version_info[:2] < (3, 8):
    raise RuntimeError("Python version >= 3.8 required.")

import builtins


CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Topic :: Software Development :: Libraries
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX :: Linux
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""


# BEFORE importing setuptools, remove MANIFEST. Otherwise it may not be
# properly updated when the contents of directories change (true for distutils,
# not sure about setuptools).
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# scipy __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__SCIPY_SETUP__ = True


def check_submodules():
    """ verify that the submodules are checked out and clean
        use `git submodule update --init`; on failure
    """
    if not os.path.exists('.git'):
        return
    with open('.gitmodules') as f:
        for l in f:
            if 'path' in l:
                p = l.split('=')[-1].strip()
                if not os.path.exists(p):
                    raise ValueError('Submodule %s missing' % p)


    proc = subprocess.Popen(['git', 'submodule', 'status'],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    for line in status.splitlines():
        if line.startswith('-') or line.startswith('+'):
            raise ValueError('Submodule not clean: %s' % line)


class concat_license_files():
    """Merge LICENSE.txt and LICENSES_bundled.txt for sdist creation

    Done this way to keep LICENSE.txt in repo as exact BSD 3-clause (see
    NumPy gh-13447).  This makes GitHub state correctly how SciPy is licensed.
    """
    def __init__(self):
        self.f1 = 'LICENSE.txt'
        self.f2 = 'LICENSES_bundled.txt'

    def __enter__(self):
        """Concatenate files and remove LICENSES_bundled.txt"""
        with open(self.f1, 'r') as f1:
            self.bsd_text = f1.read()

        with open(self.f1, 'a') as f1:
            with open(self.f2, 'r') as f2:
                self.bundled_text = f2.read()
                f1.write('\n\n')
                f1.write(self.bundled_text)

    def __exit__(self, exception_type, exception_value, traceback):
        """Restore content of both files"""
        with open(self.f1, 'w') as f:
            f.write(self.bsd_text)


from distutils.command.sdist import sdist
class sdist_checked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        check_submodules()
        with concat_license_files():
            sdist.run(self)


def get_build_ext_override():
    """
    Custom build_ext command to tweak extension building.
    """
    from numpy.distutils.command.build_ext import build_ext as npy_build_ext
    if int(os.environ.get('SCIPY_USE_PYTHRAN', 1)):
        try:
            import pythran
            from pythran.dist import PythranBuildExt
        except ImportError:
            BaseBuildExt = npy_build_ext
        else:
            BaseBuildExt = PythranBuildExt[npy_build_ext]
            _pep440 = importlib.import_module('scipy._lib._pep440')
            if _pep440.parse(pythran.__version__) < _pep440.Version('0.11.0'):
                raise RuntimeError("The installed `pythran` is too old, >= "
                                   "0.11.0 is needed, {} detected. Please "
                                   "upgrade Pythran, or use `export "
                                   "SCIPY_USE_PYTHRAN=0`.".format(
                                   pythran.__version__))
    else:
        BaseBuildExt = npy_build_ext

    class build_ext(BaseBuildExt):
        def finalize_options(self):
            super().finalize_options()

            # Disable distutils parallel build, due to race conditions
            # in numpy.distutils (Numpy issue gh-15957)
            if self.parallel:
                print("NOTE: -j build option not supported. Set NPY_NUM_BUILD_JOBS=4 "
                      "for parallel build.")
            self.parallel = None

        def build_extension(self, ext):
            # When compiling with GNU compilers, use a version script to
            # hide symbols during linking.
            if self.__is_using_gnu_linker(ext):
                export_symbols = self.get_export_symbols(ext)
                text = '{global: %s; local: *; };' % (';'.join(export_symbols),)

                script_fn = os.path.join(self.build_temp, 'link-version-{}.map'.format(ext.name))
                with open(script_fn, 'w') as f:
                    f.write(text)
                    # line below fixes gh-8680
                    ext.extra_link_args = [arg for arg in ext.extra_link_args if not "version-script" in arg]
                    ext.extra_link_args.append('-Wl,--version-script=' + script_fn)

            # Allow late configuration
            hooks = getattr(ext, '_pre_build_hook', ())
            _run_pre_build_hooks(hooks, (self, ext))

            super().build_extension(ext)

        def __is_using_gnu_linker(self, ext):
            if not sys.platform.startswith('linux'):
                return False

            # Fortran compilation with gfortran uses it also for
            # linking. For the C compiler, we detect gcc in a similar
            # way as distutils does it in
            # UnixCCompiler.runtime_library_dir_option
            if ext.language == 'f90':
                is_gcc = (self._f90_compiler.compiler_type in ('gnu', 'gnu95'))
            elif ext.language == 'f77':
                is_gcc = (self._f77_compiler.compiler_type in ('gnu', 'gnu95'))
            else:
                is_gcc = False
                if self.compiler.compiler_type == 'unix':
                    cc = sysconfig.get_config_var("CC")
                    if not cc:
                        cc = ""
                    compiler_name = os.path.basename(cc.split(" ")[0])
                    is_gcc = "gcc" in compiler_name or "g++" in compiler_name
            return is_gcc and sysconfig.get_config_var('GNULD') == 'yes'

    return build_ext


def get_build_clib_override():
    """
    Custom build_clib command to tweak library building.
    """
    from numpy.distutils.command.build_clib import build_clib as old_build_clib

    class build_clib(old_build_clib):
        def finalize_options(self):
            super().finalize_options()

            # Disable parallelization (see build_ext above)
            self.parallel = None

        def build_a_library(self, build_info, lib_name, libraries):
            # Allow late configuration
            hooks = build_info.get('_pre_build_hook', ())
            _run_pre_build_hooks(hooks, (self, build_info))
            old_build_clib.build_a_library(self, build_info, lib_name, libraries)

    return build_clib


def _run_pre_build_hooks(hooks, args):
    """Call a sequence of pre-build hooks, if any"""
    if hooks is None:
        hooks = ()
    elif not hasattr(hooks, '__iter__'):
        hooks = (hooks,)
    for hook in hooks:
        hook(*args)


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                         os.path.join(cwd, 'tools', 'cythonize.py'),
                         'scipy'],
                        cwd=cwd)
    if p != 0:
        # Could be due to a too old pip version and build isolation, check that
        try:
            # Note, pip may not be installed or not have been used
            import pip
        except (ImportError, ModuleNotFoundError):
            raise RuntimeError("Running cythonize failed!")
        else:
            _pep440 = importlib.import_module('scipy._lib._pep440')
            if _pep440.parse(pip.__version__) < _pep440.Version('18.0.0'):
                raise RuntimeError("Cython not found or too old. Possibly due "
                                   "to `pip` being too old, found version {}, "
                                   "needed is >= 18.0.0.".format(
                                   pip.__version__))
            else:
                raise RuntimeError("Running cythonize failed!")


def parse_setuppy_commands():
    """Check the commands and respond appropriately.  Disable broken commands.

    Return a boolean value for whether or not to run the build or not (avoid
    parsing Cython and template files if False).
    """
    args = sys.argv[1:]

    if not args:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands = ['--help-commands', '--name', '--version', '-V',
                     '--fullname', '--author', '--author-email',
                     '--maintainer', '--maintainer-email', '--contact',
                     '--contact-email', '--url', '--license', '--description',
                     '--long-description', '--platforms', '--classifiers',
                     '--keywords', '--provides', '--requires', '--obsoletes']

    for command in info_commands:
        if command in args:
            return False

    # Note that 'alias', 'saveopts' and 'setopt' commands also seem to work
    # fine as they are, but are usually used together with one of the commands
    # below and not standalone.  Hence they're not added to good_commands.
    good_commands = ('develop', 'sdist', 'build', 'build_ext', 'build_py',
                     'build_clib', 'build_scripts', 'bdist_wheel', 'bdist_rpm',
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg')

    for command in good_commands:
        if command in args:
            return True

    # The following commands are supported, but we need to show more
    # useful messages to the user
    if 'install' in args:
        print(textwrap.dedent("""
            Note: for reliable uninstall behaviour and dependency installation
            and uninstallation, please use pip instead of using
            `setup.py install`:

              - `pip install .`       (from a git repo or downloaded source
                                       release)
              - `pip install scipy`   (last SciPy release on PyPI)

            """))
        return True

    if '--help' in args or '-h' in sys.argv[1]:
        print(textwrap.dedent("""
            SciPy-specific help
            -------------------

            To install SciPy from here with reliable uninstall, we recommend
            that you use `pip install .`. To install the latest SciPy release
            from PyPI, use `pip install scipy`.

            For help with build/installation issues, please ask on the
            scipy-user mailing list.  If you are sure that you have run
            into a bug, please report it at https://github.com/scipy/scipy/issues.

            Setuptools commands help
            ------------------------
            """))
        return False


    # The following commands aren't supported.  They can only be executed when
    # the user explicitly adds a --force command-line argument.
    bad_commands = dict(
        test="""
            `setup.py test` is not supported.  Use one of the following
            instead:

              - `python runtests.py`              (to build and test)
              - `python runtests.py --no-build`   (to test installed scipy)
              - `>>> scipy.test()`           (run tests for installed scipy
                                              from within an interpreter)
            """,
        upload="""
            `setup.py upload` is not supported, because it's insecure.
            Instead, build what you want to upload and upload those files
            with `twine upload -s <filenames>` instead.
            """,
        upload_docs="`setup.py upload_docs` is not supported",
        easy_install="`setup.py easy_install` is not supported",
        clean="""
            `setup.py clean` is not supported, use one of the following instead:

              - `git clean -xdf` (cleans all files)
              - `git clean -Xdf` (cleans all versioned files, doesn't touch
                                  files that aren't checked into the git repo)
            """,
        check="`setup.py check` is not supported",
        register="`setup.py register` is not supported",
        bdist_dumb="`setup.py bdist_dumb` is not supported",
        bdist="`setup.py bdist` is not supported",
        flake8="`setup.py flake8` is not supported, use flake8 standalone",
        build_sphinx="`setup.py build_sphinx` is not supported, see doc/README.md",
        )
    bad_commands['nosetests'] = bad_commands['test']
    for command in ('upload_docs', 'easy_install', 'bdist', 'bdist_dumb',
                     'register', 'check', 'install_data', 'install_headers',
                     'install_lib', 'install_scripts', ):
        bad_commands[command] = "`setup.py %s` is not supported" % command

    for command in bad_commands.keys():
        if command in args:
            print(textwrap.dedent(bad_commands[command]) +
                  "\nAdd `--force` to your command to use it anyway if you "
                  "must (unsupported).\n")
            sys.exit(1)

    # Commands that do more than print info, but also don't need Cython and
    # template parsing.
    other_commands = ['egg_info', 'install_egg_info', 'rotate']
    for command in other_commands:
        if command in args:
            return False

    # If we got here, we didn't detect what setup.py command was given
    warnings.warn("Unrecognized setuptools command ('{}'), proceeding with "
                  "generating Cython sources and expanding templates".format(
                  ' '.join(sys.argv[1:])))
    return True

def check_setuppy_command():
    run_build = parse_setuppy_commands()
    if run_build:
        try:
            pkgname = 'numpy'
            import numpy
            pkgname = 'pybind11'
            import pybind11
        except ImportError as exc:  # We do not have our build deps installed
            print(textwrap.dedent(
                    """Error: '%s' must be installed before running the build.
                    """
                    % (pkgname,)))
            sys.exit(1)

    return run_build

def configuration(parent_package='', top_path=None):
    from scipy._build_utils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        if sys.platform == "darwin":
            msg = ('No BLAS/LAPACK libraries found. '
                   'Note: Accelerate is no longer supported.')
        else:
            msg = 'No BLAS/LAPACK libraries found.'
        msg += ("\n"
                "To build Scipy from sources, BLAS & LAPACK libraries "
                "need to be installed.\n"
                "See site.cfg.example in the Scipy source directory and\n"
                "https://docs.scipy.org/doc/scipy/reference/building/index.html "
                "for details.")
        raise NotFoundError(msg)

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('scipy')
    config.add_data_files(('scipy', '*.txt'))

    config.get_version('scipy/version.py')

    return config


def setup_package():
    # In maintenance branch, change np_maxversion to N+3 if numpy is at N
    # Update here, in pyproject.toml, and in scipy/__init__.py
    # Rationale: SciPy builds without deprecation warnings with N; deprecations
    #            in N+1 will turn into errors in N+3
    # For Python versions, if releases is (e.g.) <=3.9.x, set bound to 3.10
    np_minversion = '1.18.5'
    np_maxversion = '1.26.0'
    python_minversion = '3.8'
    python_maxversion = '3.12'
    if IS_RELEASE_BRANCH:
        req_np = 'numpy>={},<{}'.format(np_minversion, np_maxversion)
        req_py = '>={},<{}'.format(python_minversion, python_maxversion)
    else:
        req_np = 'numpy>={}'.format(np_minversion)
        req_py = '>={}'.format(python_minversion)

    # Rewrite the version file every time
    write_version_py('.')

    cmdclass = {'sdist': sdist_checked}

    metadata = dict(
        name='scipy',
        maintainer="SciPy Developers",
        maintainer_email="scipy-dev@python.org",
        description=DOCLINES[0],
        long_description="\n".join(DOCLINES[2:]),
        url="https://www.scipy.org",
        download_url="https://github.com/scipy/scipy/releases",
        project_urls={
            "Bug Tracker": "https://github.com/scipy/scipy/issues",
            "Documentation": "https://docs.scipy.org/doc/scipy/reference/",
            "Source Code": "https://github.com/scipy/scipy",
        },
        license='BSD',
        cmdclass=cmdclass,
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        install_requires=[req_np],
        python_requires=req_py,
        zip_safe=False,
    )

    if "--force" in sys.argv:
        run_build = True
        sys.argv.remove('--force')
    else:
        # Raise errors for unsupported commands, improve help output, etc.
        run_build = check_setuppy_command()

    # Disable OSX Accelerate, it has too old LAPACK
    os.environ['ACCELERATE'] = 'None'

    # This import is here because it needs to be done before importing setup()
    # from numpy.distutils, but after the MANIFEST removing and sdist import
    # higher up in this file.
    from setuptools import setup

    if run_build:
        from numpy.distutils.core import setup

        # Customize extension building
        cmdclass['build_ext'] = get_build_ext_override()
        cmdclass['build_clib'] = get_build_clib_override()

        if not 'sdist' in sys.argv:
            # Generate Cython sources, unless we're creating an sdist
            # Cython is a build dependency, and shipping generated .c files
            # can cause problems (see gh-14199)
            generate_cython()

        metadata['configuration'] = configuration
    else:
        # Don't import numpy here - non-build actions are required to succeed
        # without NumPy for example when pip is used to install Scipy when
        # NumPy is not yet present in the system.

        # Version number is added to metadata inside configuration() if build
        # is run.
        metadata['version'] = get_version_info('.')[0]

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
