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

DOCLINES = __doc__.split("\n")

import os
import sys
import subprocess
import textwrap
import warnings
import sysconfig


if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[:2] < (3, 4):
    raise RuntimeError("Python version 2.7 or >= 3.4 required.")

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3
Programming Language :: Python :: 3.4
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""

MAJOR = 1
MINOR = 2
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
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
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


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


def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('scipy/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load scipy/__init__.py
        import imp
        version = imp.load_source('scipy.version', 'scipy/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='scipy/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SCIPY SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


try:
    from sphinx.setup_command import BuildDoc
    HAVE_SPHINX = True
except Exception:
    HAVE_SPHINX = False

if HAVE_SPHINX:
    class ScipyBuildDoc(BuildDoc):
        """Run in-place build before Sphinx doc build"""
        def run(self):
            ret = subprocess.call([sys.executable, sys.argv[0], 'build_ext', '-i'])
            if ret != 0:
                raise RuntimeError("Building Scipy failed!")
            BuildDoc.run(self)


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


from distutils.command.sdist import sdist
class sdist_checked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        check_submodules()
        sdist.run(self)


def get_build_ext_override():
    """
    Custom build_ext command to tweak extension building.
    """
    from numpy.distutils.command.build_ext import build_ext as old_build_ext

    class build_ext(old_build_ext):
        def build_extension(self, ext):
            # When compiling with GNU compilers, use a version script to
            # hide symbols during linking.
            if self.__is_using_gnu_linker(ext):
                export_symbols = self.get_export_symbols(ext)
                text = '{global: %s; local: *; };' % (';'.join(export_symbols),)

                script_fn = os.path.join(self.build_temp, 'link-version-{}.map'.format(ext.name))
                with open(script_fn, 'w') as f:
                    f.write(text)
                    ext.extra_link_args.append('-Wl,--version-script=' + script_fn)

            old_build_ext.build_extension(self, ext)

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
                    compiler_name = os.path.basename(cc)
                    is_gcc = "gcc" in compiler_name or "g++" in compiler_name
            return is_gcc and sysconfig.get_config_var('GNULD') == 'yes'

    return build_ext


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                         os.path.join(cwd, 'tools', 'cythonize.py'),
                         'scipy'],
                        cwd=cwd)
    if p != 0:
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
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg',
                     'build_sphinx')

    for command in good_commands:
        if command in args:
            return True

    # The following commands are supported, but we need to show more
    # useful messages to the user
    if 'install' in args:
        print(textwrap.dedent("""
            Note: if you need reliable uninstall behavior, then install
            with pip instead of using `setup.py install`:

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
    warnings.warn("Unrecognized setuptools command, proceeding with "
                  "generating Cython sources and expanding templates")
    return True


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info, NotFoundError, numpy_info
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import (get_g77_abi_wrappers, split_fortran_files)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        msg = 'No lapack/blas resources found.'
        if sys.platform == "darwin":
            msg = ('No lapack/blas resources found. '
                   'Note: Accelerate is no longer supported.')
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
    # Rewrite the version file every time
    write_version_py()

    cmdclass = {'sdist': sdist_checked}
    if HAVE_SPHINX:
        cmdclass['build_sphinx'] = ScipyBuildDoc

    # Figure out whether to add ``*_requires = ['numpy']``.
    # We don't want to do that unconditionally, because we risk updating
    # an installed numpy which fails too often.  Just if it's not installed, we
    # may give it a try.  See gh-3379.
    try:
        import numpy
    except ImportError:  # We do not have numpy installed
        build_requires = ['numpy>=1.8.2']
    else:
        # If we're building a wheel, assume there already exist numpy wheels
        # for this platform, so it is safe to add numpy to build requirements.
        # See gh-5184.
        build_requires = (['numpy>=1.8.2'] if 'bdist_wheel' in sys.argv[1:]
                          else [])

    metadata = dict(
        name='scipy',
        maintainer="SciPy Developers",
        maintainer_email="scipy-dev@python.org",
        description=DOCLINES[0],
        long_description="\n".join(DOCLINES[2:]),
        url="https://www.scipy.org",
        download_url="https://github.com/scipy/scipy/releases",
        license='BSD',
        cmdclass=cmdclass,
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        test_suite='nose.collector',
        setup_requires=build_requires,
        install_requires=build_requires,
        python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*',
    )

    if "--force" in sys.argv:
        run_build = True
        sys.argv.remove('--force')
    else:
        # Raise errors for unsupported commands, improve help output, etc.
        run_build = parse_setuppy_commands()

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

        cwd = os.path.abspath(os.path.dirname(__file__))
        if not os.path.exists(os.path.join(cwd, 'PKG-INFO')):
            # Generate Cython sources, unless building from source release
            generate_cython()

        metadata['configuration'] = configuration
    else:
        # Don't import numpy here - non-build actions are required to succeed
        # without Numpy for example when pip is used to install Scipy when
        # Numpy is not yet present in the system.

        # Version number is added to metadata inside configuration() if build
        # is run.
        metadata['version'] = get_version_info()[0]

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
