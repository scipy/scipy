def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('_rcont', parent_package, top_path)

    config.add_library(
        '_rcont',
        sources=['rcont.c']
    )
    config.add_extension(
        'rcont',
        libraries=['_rcont'],
        sources=['rcont.c']
    )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
