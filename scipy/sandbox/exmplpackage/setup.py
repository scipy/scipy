#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):
    # The following two lines with `return config` constitutes a
    # minimal contents of configuration(..) that is suitable for pure
    # Python packages.
    from numpy.distutils.misc_util import Configuration
    config = Configuration('exmplpackage',parent_package,top_path)

    # include test scripts from tests
    config.add_data_dir('tests')

    # exmplpackage contains Python sub-package yyy
    config.add_subpackage('yyy')

    # exmplpackage generates source code, that will be processed with f2py
    def generate_spam_pyf(ext, build_dir):
        from distutils.dep_util import newer
        target = join(build_dir,'spam.pyf')
        source = ext.depends[0]
        if newer(source,target):
            fin = open(source)
            body = fin.read()
            fin.close()
            fout = open(target,'w')
            fout.write('python module spam\n%s\nend python module spam' % body)
            fout.close()
        return target
    config.add_extension('spam',
                         sources = [generate_spam_pyf],
                         depends = ['spam_src.pyf'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
