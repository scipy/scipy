from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('numexpr', parent_package, top_path)
    config.add_extension('interpreter',
                         sources = ['interpreter.c'],
                         depends = ['interp_body.c',
                                    'complex_functions.inc'],
                         extra_compile_args=['-O2', '-funroll-all-loops'],
                         )
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(author='David M. Cooke',
          author_email='cookedm@physics.mcmaster.ca',
          version='0.1',
          zip_safe=False,
          **configuration(top_path='').todict())
