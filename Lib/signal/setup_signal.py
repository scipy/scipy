import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):    
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)
    
    config = default_config_dict()

    config['packages'].append(parent_package+'signal')
    
    sources = ['sigtoolsmodule.c','firfilter.c','medianfilter.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'signal.sigtools', sources)
    config['ext_modules'].append(ext)
    
    sources = ['splinemodule.c','S_bspline_util.c','D_bspline_util.c',
               'C_bspline_util.c','Z_bspline_util.c','bspline_util.c']
    sources = [os.path.join(local_path,x) for x in sources]               
    ext = Extension(parent_package+'signal.spline',sources)
    config['ext_modules'].append(ext)

    return config
