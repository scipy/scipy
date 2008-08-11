
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('models',parent_package,top_path)

    config.add_subpackage('family')
    config.add_subpackage('robust')

    config.add_data_dir('tests')

    config.add_extension('_hbspline',
                         sources=['src/bspline_ext.c',
                                  'src/bspline_impl.c'],
    )
    return config

if __name__ == '__main__':

    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
