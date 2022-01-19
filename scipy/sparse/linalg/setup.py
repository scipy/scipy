
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('linalg', parent_package, top_path)

    config.add_subpackage('_isolve')
    config.add_subpackage('_dsolve')
    config.add_subpackage('_eigen')

    config.add_data_dir('tests')

    # PROPACK
    config.add_subpackage('_propack')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
