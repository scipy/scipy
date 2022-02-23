from os.path import join


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('_levy_stable', parent_package, top_path)

    config.add_library(
        '_levyst',
        sources=[join('c_src', 'levyst.c')],
        headers=[join('c_src', 'levyst.h')]
    )
    config.add_extension(
        'levyst',
        libraries=['_levyst'],
        sources=['levyst.c']
    )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
