#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('image', parent_package, top_path)
    config.add_data_files("""ciexyz31_1.txt ciexyz64_1.txt ciexyzjv.txt
        linss2_10e_1.txt sbrgb2.txt""".split())

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration())
