import os
import sys

os.environ['NO_SCIPY_IMPORT']='SciPy/setup.py'

def configuration(parent_package='',top_path=None):
      from numpy.distutils.misc_util import Configuration
      config = Configuration(None, parent_package, top_path)
      config.set_options(ignore_setup_xxx_py=True,
                         assume_default_configuration=True,
                         delegate_options_to_subpackages=True,
                         quiet=True)

      config.add_subpackage('Lib')
      config.add_data_files(('scipy','*.txt'))

      config.get_version('Lib/version.py') # sets config.version

      return config

def setup_package():

    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'Lib')) # to retrive version

    try:
        from version import version as version
        setup(
            name = 'scipy',
            version = version, # will be overwritten by configuration version
            maintainer = "SciPy Developers",
            maintainer_email = "scipy-dev@scipy.org",
            description = "Scientific Algorithms Library for Python",
            url = "http://www.scipy.org",
            license = 'BSD',
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
