import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)
    
    packages = []
    packages.append(parent_package+'cow')
    #packages.append(parent_package+'cow.tests') 
    return {'packages': packages}