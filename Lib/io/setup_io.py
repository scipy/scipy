import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path

def configuration(parent_package=''):
    local_path = get_path(__name__)
    if parent_package:
        parent_package += '.'
    packages = []
    ext_modules = []
    
    from scipy_distutils.core import Extension
    packages.append(parent_package+'io')
    
    sources = ['numpyiomodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'io.numpyio',sources)
    ext_modules.append(ext)
    #packages.append(parent_package+'io.tests') 
    results = {'packages': packages,
               'ext_modules': ext_modules,
              }
    return results          
