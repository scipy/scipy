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
from glob import glob

# Note that when running bdist_rpm, scipy_core directory does
# not exist. So, scipy_distutils must be installed before
# running bdist_rpm.
sys.path.insert(0,'scipy_core')
try:
    import scipy_distutils
    # Declare all scipy_distutils related imports here.
    from scipy_distutils.misc_util import default_config_dict
    from scipy_distutils.misc_util import get_path, merge_config_dicts
    from scipy_distutils.core import setup
finally:
    del sys.path[0]

#------ drop-to-Lib packages --------

def get_packages(path,ignore_packages=[],
                 parent='scipy',parent_path=None):

    config_list = []

    for info_file in glob(os.path.join(path,'*','info_*.py')):

        package_name = os.path.basename(os.path.dirname(info_file))
        if package_name != os.path.splitext(os.path.basename(info_file))[0][5:]:
            print '  !! Mismatch of package name %r and %s' \
                  % (package_name, info_file)
            continue

        if package_name in ignore_packages:
            continue

        sys.path.insert(0,os.path.dirname(info_file))
        try:
            exec 'import %s as info_module' \
                 % (os.path.splitext(os.path.basename(info_file))[0])
            if not getattr(info_module,'ignore',0):
                exec 'import setup_%s as setup_module' % (package_name)
                if getattr(info_module,'standalone',0):
                    args = ('',)
                else:
                    args = (parent,)
                if setup_module.configuration.func_code.co_argcount>1:
                    args = args + (parent_path,)
                config = setup_module.configuration(*args)
                config_list.append(config)
        finally:
            del sys.path[0]

    return config_list
#-------------------------------

def setup_package(ignore_packages=[]):

    if os.path.isdir('scipy_core'):
        # Applying the same commands to scipy_core.
        # Results can be found in scipy_core directory.
        c = '%s %s %s' % (sys.executable,
                          os.path.join('scipy_core','setup.py'),
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
        from scipy_version import scipy_version

        config_list = [{'name': 'SciPy',
                        'packages':['scipy','scipy.tests'],
                        'package_dir':
                        {'scipy':'Lib',
                         'scipy.tests':os.path.join('Lib','tests')}}]

        for d in ['Lib']:
            config_list += get_packages(os.path.join(local_path,d),
                                        ignore_packages,
                                        parent_path=local_path)

        config_dict = merge_config_dicts(config_list)

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
        #'sparse',
        #'kiva','freetype','chaco','traits',
        ]
    if sys.platform=='win32':
        # XXX: Fix xplt for windows! It is failing due
        # to missing m.lib file while running pygist_config.
        ignore_packages.append('xplt')
        ignore_packages.append('sparse')
    setup_package(ignore_packages)
