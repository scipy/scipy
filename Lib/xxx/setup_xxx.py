#!/usr/bin/env python

import os

def configuration(parent_package='', parent_path=None):

    # The following three lines constitute minimal contents
    # of configuration(..) that is suitable for pure Python
    # packages.
    package = 'xxx'
    from scipy_distutils.misc_util import default_config_dict, dot_join, get_path
    config = default_config_dict(package,parent_package)

    # local_path should be joined with local source files
    # when specifying sources below
    local_path = get_path(__name__)

    # xxx contains pure Python sub-package yyy
    config['packages'].append(dot_join(parent_package,package,'yyy'))
    # And don't forget adding sub-packages of yyy:
    config['packages'].append(dot_join(parent_package,package,'yyy.tests'))

    # xxx has extension module ..
    from scipy_distutils.core import Extension

    
    # xxx has f2py generated extension module ..
    # xxx has swig generated extension module ..

    # xxx generates source code
    from scipy_distutils.misc_util import SourceGenerator
    def generate_spam_pyf(target, sources):
        fin = open(sources[0])
        body = fin.read()
        fin.close()
        fout = open(target,'w')
        fout.write('python module spam\n%s\nend python module spam' % body)
        fout.close()
    pyf_file = SourceGenerator(generate_spam_pyf,
                               target='spam.pyf',
                               sources=[os.path.join(local_path,'spam_src.pyf')])
    ext = Extension(name=dot_join(parent_package,package,'spam'),
                    sources=[pyf_file],
                    depends = pyf_file.sources)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
