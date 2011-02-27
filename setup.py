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
import warnings
import subprocess
import shutil
import re

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""

MAJOR               = 0
MINOR               = 9
MICRO               = 0
ISRELEASED          = True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

os.environ['NO_SCIPY_IMPORT']='SciPy/setup.py'

# Return the svn version as a string, raise a ValueError otherwise
def svn_version():
    from numpy.compat import asstr

    env = os.environ.copy()
    env['LC_ALL'] = 'C'
    try:
        out = subprocess.Popen(['svn', 'info'], stdout=subprocess.PIPE,
                env=env).communicate()[0]
    except OSError:
        warnings.warn(" --- Could not run svn info --- ")
        return ""

    r = re.compile('Revision: ([0-9]+)')
    svnver = None
    for line in asstr(out).split('\n'):
        m = r.match(line)
        if m:
            svnver = m.group(1)

    if not svnver:
        raise ValueError("Error while parsing svn version ?")
    return svnver

FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.dev'
    # If in git or something, bypass the svn rev
    if os.path.exists('.svn'):
        FULLVERSION += svn_version()

def write_version_py(filename='scipy/version.py'):
    cnt = """\
# THIS FILE IS GENERATED FROM SCIPY SETUP.PY
short_version='%(version)s'
version='%(version)s'
release=%(isrelease)s

if not release:
    version += '.dev'
    import os
    svn_version_file = os.path.join(os.path.dirname(__file__),
                                    '__svn_version__.py')
    if os.path.isfile(svn_version_file):
        import imp
        svn = imp.load_module('scipy.__svn_version__',
                              open(svn_version_file),
                              svn_version_file,
                              ('.py','U',1))
        version += svn.version
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION, 'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('scipy')
    config.add_data_files(('scipy','*.txt'))

    config.get_version('scipy/version.py')

    return config

def setup_package():

    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    src_path = local_path
    if sys.version_info[0] == 3:
        src_path = os.path.join(local_path, 'build', 'py3k')
        sys.path.insert(0, os.path.join(local_path, 'tools'))
        import py3tool
        print("Converting to Python3 via 2to3...")
        py3tool.sync_2to3('scipy', os.path.join(src_path, 'scipy'))

        site_cfg = os.path.join(local_path, 'site.cfg')
        if os.path.isfile(site_cfg):
            shutil.copy(site_cfg, src_path)

    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'scipy')) # to retrive version

    # Run build
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    if os.path.exists('scipy/version.py'): os.remove('scipy/version.py')
    write_version_py()

    try:
        setup(
            name = 'scipy',
            maintainer = "SciPy Developers",
            maintainer_email = "scipy-dev@scipy.org",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "http://www.scipy.org",
            download_url = "http://sourceforge.net/project/showfiles.php?group_id=27747&package_id=19531",
            license = 'BSD',
            classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
