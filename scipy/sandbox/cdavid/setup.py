#! /usr/bin/env python
# Last Change: Tue Nov 28 03:00 PM 2006 J

""" toolbox of Cournapeau David. Implements various things such as
autocorrelation with lag argument, lpc coefficients computation.

2006, David Cournapeau

LICENSE: the license of pyaudio is the same than scipy"""

from os.path import join
import os

from info import version as cdavid_version

DISTNAME        = 'cdavid'
VERSION         = cdavid_version
DESCRIPTION     ='A scipy package for various speech processing tools lpc'
MAINTAINER      ='David Cournapeau',
MAINTAINER_EMAIL='david@ar.media.kyoto-u.ac.jp',
URL             ='http://ar.media.kyoto-u.ac.jp/members/david',
LICENSE         = 'BSD'

def configuration(parent_package='',top_path=None, package_name=DISTNAME):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path,
             version     = VERSION)
    config.add_data_dir('tests')
    #config.add_data_dir('profile_data')
    config.add_extension('gabsig',
                         #define_macros=[('LIBSVM_EXPORTS', None),
                         #               ('LIBSVM_DLL', None)],
                         sources=[join('src', 'levinson.c'),
                         join('src', 'autocorr_nofft.c'),
                         join('src', 'lpc.c')])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
