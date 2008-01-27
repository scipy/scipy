#!/usr/bin/env python

from os.path import join

def configuration(parent_package = '', top_path = None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cluster', parent_package, top_path)

    config.add_data_dir('tests')

    #config.add_extension('_vq',
    #    sources=[join('src', 'vq_module.c'), join('src', 'vq.c')],
    #    include_dirs = [get_numpy_include_dirs()])
    config.add_sconscript('SConstruct')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer = "SciPy Developers",
          author = "Eric Jones",
          maintainer_email = "scipy-dev@scipy.org",
          description = "Clustering Algorithms (Information Theory)",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration(top_path='').todict()
          )
