
def configuration(parent_package='',top_path=None, package_name='models'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path)

    config.add_subpackage('*')

    config.add_data_dir('tests')

    try:
        from scipy.sandbox.models.bspline_module import mod
        n, s, d = weave_ext(mod)
        config.add_extension(n, s, **d)
    except ImportError: pass

    return config

def weave_ext(mod):
    d = mod.setup_extension().__dict__
    n = d['name']; del(d['name'])
    s = d['sources']; del(d['sources'])
    return n, s, d

if __name__ == '__main__':

    from numpy.distutils.core import setup
    setup(**configuration(top_path='', package_name='scipy.sandbox.models').todict())
