#!/usr/bin/env python
"""
Installing scipy:
  python setup.py install
  python setup.py build build_flib install --fcompiler=Gnu

Creating scipy distribution:
  python setup.py sdist -f

"""

import os
import sys

sys.path.insert(0,'scipy_core')
try:
    # Declare all scipy_distutils related imports here.
    from scipy_distutils.misc_util import default_config_dict
    from scipy_distutils.misc_util import get_path, merge_config_dicts
    from scipy_distutils.core import setup
    if sys.platform == 'win32':
        # This forces g77 for windows platforms:
        from scipy_distutils.mingw32_support import *
finally:
    del sys.path[0]

#-------------------------------

standard_packages = ['io','linalg',
                     'special','signal','stats',
                     'interpolate','integrate','optimize',
                     'cluster','cow','ga','fftpack']
standard_packages = [os.path.join('Lib',p) for p in standard_packages]

graphics_packages = ['plt','gplt','xplt']
graphics_packages = [os.path.join('Lib',p) for p in graphics_packages]

chaco_packages = ['chaco','kiva','traits','freetype']
chaco_packages = [os.path.join('Lib_chaco',p) for p in chaco_packages]

core_packages = ['scipy_distutils','scipy_test','scipy_base']
core_packages = [os.path.join('scipy_core',p) for p in core_packages]

#---------------

parent_package = 'scipy'
scipy_packages = standard_packages + graphics_packages
#scipy_packages += chaco_packages

#---------------

# these packages aren't nested under scipy
separate_packages = ['gui_thread','weave']
separate_packages = [os.path.join('Lib',p) for p in separate_packages]
separate_packages += core_packages

#-------------------------------

def get_package_config(name, parent=parent_package):
    sys.path.insert(0,name)
    try:
        mod = __import__('setup_'+os.path.basename(name))
        config = mod.configuration(parent)
    finally:
        del sys.path[0]
    return config

def get_separate_package_config(name):
    return get_package_config(name,'')

def setup_package():
    old_path = os.getcwd()
    path = get_path(__name__)
    os.chdir(path)
    sys.path.insert(0,path)
    try:
        config_list = [{'packages':['scipy','scipy.tests'],
                        'package_dir':
                        {'scipy':'Lib',
                         'scipy.tests':os.path.join('Lib','tests')}}]
        config_list += map(get_package_config,scipy_packages)
        config_list += map(get_separate_package_config,separate_packages)

        config_dict = merge_config_dicts(config_list)

        sys.path.insert(0,'Lib')
        from scipy_version import scipy_version
        del sys.path[0]

        print 'SciPy Version %s' % scipy_version
        setup (name = "SciPy",
               version = scipy_version,
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
    setup_package()
