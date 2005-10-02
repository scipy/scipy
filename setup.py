#!/usr/bin/env python
"""
Installing scipy:
  python setup.py install
  python setup.py fc_config --fcompiler=Gnu install

Building extensions in-place:
  python setup.py build_src build_ext --inplace

Creating scipy distribution:
  python setup.py sdist -f            #-> Scipy_complete
  python setup.py sdist_packagers -f  #-> Scipy and Scipy_core
"""

import os
import sys
from glob import glob

# Use 'setup.py sdist_packagers -f' to create separate tar-balls
# for Scipy and Scipy_core.
if 'sdist_packagers' in sys.argv:
    sys.argv[sys.argv.index('sdist_packagers')] = 'sdist'
    command_sdist = 1
elif 'bdist_rpm' in sys.argv:
    command_sdist = 1
else:
    command_sdist = 0

# Try to detect if we are building in the scipy source directory (the most
# common case).  We do this by checking for a scipy_core subdir.  If this is
# the case, we add the current dir+/scipy_core to sys.path and to the
# environment's PYTHONPATH, so that bdist_rpm works without requiring
# scipy.distutils to be previously installed. However, such a situation
# is abnormal because building scipy requires f2py, and f2py in turn
# requires scipy.distutils (though, f2py can be installed without
# scipy.distutils).

scipy_core_dir = os.path.join(os.getcwd(),'scipy_core')
if os.path.isdir(scipy_core_dir):
    ppath = os.environ.setdefault('PYTHONPATH',scipy_core_dir)
    if ppath != scipy_core_dir:
        ppath = '%s%s%s' % (scipy_core_dir,os.pathsep,ppath)
    os.environ['PYTHONPATH'] = ppath
sys.path.insert(0,scipy_core_dir)

try:
    import scipy.distutils
    # Declare all scipy.distutils related imports here.
    from scipy.distutils.misc_util import default_config_dict
    from scipy.distutils.misc_util import get_path, merge_config_dicts
    from scipy.distutils.misc_util import get_subpackages, generate_config_py
    from scipy.distutils.core import setup, Extension
finally:
    del sys.path[0]

#-------------------------------

def setup_package(ignore_packages=[]):

    if command_sdist and os.path.isdir('scipy_core'):
        # Applying the same commands to scipy_core.
        # Results can be found in scipy_core directory.
        c = '%s %s %s' % (sys.executable,
                          os.path.abspath(os.path.join('scipy_core','setup.py')),
                          ' '.join(sys.argv[1:]))
        print c
        s = os.system(c)
        assert not s,'failed on scipy_core'

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))

    os.chdir(local_path)
    sys.path.insert(0,os.path.join(local_path,'Lib'))
    try:
        from scipy_version import scipy_version

        # If a minor version number is odd then this indicates
        # development version from CVS. Otherwise, its release version.
        # Uncomment when making releases:
        #if not command_sdist: scipy_version = '0.3.3'

        packages_path = ['Lib']
        if command_sdist:
            name = 'scipy'
        else:
            name = 'scipy_complete'
            packages_path.append('scipy_core')
        config_list = [{'name': name,
                        'packages':['scipy','scipy.tests'],
                        'package_dir':
                        {'scipy':'Lib',
                         'scipy.tests':os.path.join('Lib','tests')}}]

        for d in packages_path:
            config_list += get_subpackages(os.path.join(local_path,d),
                                           parent = 'scipy',
                                           parent_path=local_path,
                                           ignore_packages = ignore_packages)
        config_dict = merge_config_dicts(config_list)

        ext = Extension(name='scipy.config',
                        sources=[generate_config_py])
        config_dict['ext_modules'].append(ext)

        print 'SciPy Version %s' % scipy_version
        setup( version = scipy_version,
               maintainer = "SciPy Developers",
               maintainer_email = "scipy-dev@scipy.org",
               description = "Scientific Algorithms Library for Python",
               license = "SciPy License (BSD Style)",
               url = "http://www.scipy.org",
               **config_dict
              )

    finally:
        del sys.path[0]
        os.chdir(old_path)

if __name__ == "__main__":
    ignore_packages = [
        # list of package names that will not be built.
        ]
    setup_package(ignore_packages)
