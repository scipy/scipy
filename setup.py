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

# Note that when running bdist_rpm, scipy_core directory does
# not exist. So, scipy_distutils must be installed before
# running bdist_rpm.
sys.path.insert(0,'scipy_core')
try:
    import scipy_distutils
    # Declare all scipy_distutils related imports here.
    from scipy_distutils.misc_util import default_config_dict
    from scipy_distutils.misc_util import get_path, merge_config_dicts
    from scipy_distutils.misc_util import get_subpackages, generate_config_py
    from scipy_distutils.core import setup, Extension
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
    # setup files of subpackages require scipy_core:
    sys.path.insert(0,os.path.join(local_path,'scipy_core'))
    try:
        
        #from scipy_version import scipy_version
        scipy_version = '0.3.2'

        packages_path = ['Lib']
        if command_sdist:
            name = 'SciPy'
        else:
            name = 'SciPy_complete'
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
        del sys.path[0]
        os.chdir(old_path)

if __name__ == "__main__":
    ignore_packages = [
        # list of package names that will not be built.
        ]
    setup_package(ignore_packages)
