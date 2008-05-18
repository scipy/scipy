#!/usr/bin/env python
from os.path import join

from numpy.distutils.system_info import get_info

SRC = {}
SRC['fftw'] = [join('fftw/src', i) for i in 
                   ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
SRC['fftw3'] = [join('fftw3/src', i) for i in 
                ['zfft.cxx', 'drfft.cxx', 'zfftnd.cxx']]
SRC['mkl'] = [join('mkl/src', i) for i in 
                ['zfft.cxx', 'zfftnd.cxx']]

def build_backend(config, name, config_name = None):
    if config_name is None:
        config_name = name

    info = get_info(config_name)
    if info:
        config.add_library("%s_backend" % name,
                           sources = SRC[name],
                           include_dirs = ["common",
                                           info['include_dirs']])
        config.add_extension("%s._%s" % (name, name), 
                sources = ["%s/%s.pyf.src" % (name, name)], 
                extra_info = info, libraries = ["%s_backend" % name])

    config.add_subpackage(name)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('backends', parent_package, top_path)
        
    # It is a bit of an hack: we build the extensions here, but we add
    # subpackages. The reason why we build the extensions here is because we
    # share files between backends, and scons cannot use sources outside its
    # top directory.
    build_backend(config, 'fftw', 'fftw2')
    build_backend(config, 'fftw3')
    build_backend(config, 'mkl')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
