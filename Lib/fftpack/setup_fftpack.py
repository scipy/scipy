#!/usr/bin/env python
# Created by Pearu Peterson, August 2002

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path,dot_join,\
     default_config_dict,dict_append
from scipy_distutils.system_info import get_info,FFTWNotFoundError,\
     DJBFFTNotFoundError

def configuration(parent_package=''):
    package_name = 'fftpack'
    fftw_info = get_info('fftw') or get_info('dfftw')
    if not fftw_info:
        print FFTWNotFoundError.__doc__
    djbfft_info = get_info('djbfft')
    if not djbfft_info:
        print DJBFFTNotFoundError.__doc__
    #djbfft_info = None
    #fftw_info = None
    config = default_config_dict(package_name,parent_package)
    local_path = get_path(__name__)
    test_path = os.path.join(local_path,'tests')
    config['packages'].append(dot_join(parent_package,package_name,'tests'))
    config['package_dir'][package_name+'.tests'] = test_path

    dfftpack = glob(os.path.join(local_path,'dfftpack','*.f'))
    config['fortran_libraries'].append(('dfftpack',{'sources':dfftpack}))
    
    sources = ['fftpack.pyf','src/zfft.c','src/drfft.c','src/zrfft.c',
               'src/zfftnd.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext_args = {
        'name': dot_join(parent_package,package_name,'fftpack'),
        'sources': sources,
        'libraries': ['dfftpack'],
        }
    if fftw_info:
        dict_append(ext_args,**fftw_info)
    if djbfft_info:
        dict_append(ext_args,**djbfft_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['convolve.pyf','src/convolve.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext_args = {
        'name': dot_join(parent_package,package_name,'convolve'),
        'sources': sources,
        'libraries': ['dfftpack'],
        }
    if fftw_info:
        dict_append(ext_args,**fftw_info)
    if djbfft_info:
        dict_append(ext_args,**djbfft_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    import f2py2e
    if f2py2e.__version__.version < '2.23.190-1367':
        import warnings
        warnings.warn(\
        '\n'+'WARNING!'*9+'\n\n'
        '  fftpack2 requires F2PY version 2.23.190-1367 or higher but got %s\n\n'%
        f2py2e.__version__.version
        +'WARNING!'*9+'\n\n')
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(version='0.3',
          description='fftpack - Discrete Fourier Transform package',
          author='Pearu Peterson',
          author_email = 'pearu@cens.ioc.ee',
          license = 'SciPy License (BSD Style)',
          **configuration())
