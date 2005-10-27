import os
import sys

def setup_package():
    os.environ['NO_SCIPY_IMPORT']='1'

    from scipy.distutils.core import setup
    from scipy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'Lib'))

    try:
        config = Configuration(\
                               maintainer = "SciPy Developers",
                               maintainer_email = "scipy-dev@scipy.org",
                               description = "Scientific Algorithms Library for Python",
                               url = "http://www.scipy.org",
                               license = 'BSD',
                               )
        # Force scipy to be a package (its __init__.py file comes from scipy_core)
        config.packages.append('scipy')
        config.package_dir['scipy'] = os.path.join(config.local_path,'Lib')
        
        config.add_subpackage('Lib')
        config.name = 'SciPy'
        from scipy_version import scipy_version as version
        config.dict_append(version=version)

        print config.name,'version',config.version

        setup( **config.todict() )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return

if __name__ == '__main__':
    setup_package()
