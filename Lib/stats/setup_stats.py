import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)

    config = default_config_dict()
    config['packages'].append(parent_package+'stats')
    #config['packages'].append(parent_package+'stats.tests') 

    # Extension
    sources = ['randmodule.c','ranlib_all.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,'stats.rand'),sources)
    config['ext_modules'].append(ext)

    return config
