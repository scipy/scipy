
def configuration(parent_package='',top_path=None, setup_name = 'setupscons.py'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('scipy',parent_package,top_path, setup_name = 'setupscons.py')
    config.add_subpackage('cluster')
    config.add_subpackage('fftpack')
    config.add_subpackage('integrate')
    config.add_subpackage('interpolate')
    config.add_subpackage('io')
    config.add_subpackage('lib')
    config.add_subpackage('linalg')
    config.add_subpackage('linsolve')
    config.add_subpackage('maxentropy')
    config.add_subpackage('misc')
    config.add_subpackage('odr')
    config.add_subpackage('optimize')
    config.add_subpackage('sandbox')
    config.add_subpackage('signal')
    config.add_subpackage('sparse')
    config.add_subpackage('special')
    config.add_subpackage('splinalg')
    config.add_subpackage('stats')
    config.add_subpackage('ndimage')
    config.add_subpackage('stsci')
    config.add_subpackage('weave')
    config.add_subpackage('testing')
    config.make_svn_version_py()  # installs __svn_version__.py
    config.scons_make_config_py() # installs __config__.py
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
