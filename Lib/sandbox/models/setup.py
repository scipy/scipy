import os, glob, string, shutil
import sys
sys.path.insert(0,"lib")
from distutils.core import setup
from neuroimaging import packages, __version__, __doc__

packages = (
    'models',
    'models.family',
    'models.robust',
    'models.tests',
    )

def main():
        
    setup( name = 'models',
           version = __version__,
           description = """
           This is a package to try to add statistical models to scipy,
           originally in neuroimaging.scipy.org repository.
           """,
           author = 'Jonathan Taylor + hopefully others',
           author_email = 'jonathan.taylor@stanford.edu',
           ext_package = 'models',
           packages=['']+list(packages),
           package_dir = {'': 'lib'},
           url = 'http://neuroimaging.scipy.org',
           long_description = __doc__)

if __name__ == "__main__": main()
