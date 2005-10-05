import os
import sys

def configuration(parent_package='',top_path=None):
    from scipy.distutils.misc_util import Configuration
    config = Configuration()
    config.add_subpackage('Lib')
    return config.todict()

def setup_package():
    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))

    os.chdir(local_path)
    sys.path.insert(0,os.path.join(local_path,'Lib'))
    try:
        config_dict = configuration(top_path='')
        from scipy_version import scipy_version
        from scipy.distutils.core import setup
        setup( version = scipy_version,
               name = 'scipy',
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
    return

if __name__ == '__main__':
    setup_package()
