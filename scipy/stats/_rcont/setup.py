from os.path import join
import numpy as np


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('_rcont', parent_package, top_path)

    config.add_extension(
        'rcont',
        sources=['rcont.c', '_rcont.c', 'logfactorial.c'],
        include_dirs=[np.get_include()],
        library_dirs=[join(np.get_include(),
                      '..', '..', 'random', 'lib'),
                      join(np.get_include(),
                      '..', 'lib')],
        libraries=['npyrandom', 'npymath']
    )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
