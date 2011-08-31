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


CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""


# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

os.environ['NO_SCIPY_IMPORT'] = 'SciPy/setup.py'

sys.path.insert(0, os.path.dirname(__file__))
try:
    setup_py = __import__("setup")
    write_version_py = setup_py.write_version_py
finally:
    sys.path.pop(0)


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path, setup_name = "setupscons.py")
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('scipy')
    config.add_data_files(('scipy','*.txt'))

    config.get_version('scipy/version.py') # sets config.version

    return config


def setup_package():

    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'scipy')) # to retrive version

    # Rewrite the version file everytime
    if os.path.exists('scipy/version.py'):
        os.remove('scipy/version.py')

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
            classifiers=filter(None, CLASSIFIERS.split('\n')),
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return


if __name__ == '__main__':
    setup_package()
