#!/usr/bin/env python
# !!! This is not used at all by the actual setup.py package
# !!! It's just nice rebuilding cluster in place.

import os, sys, string, re
from glob import glob

# Check for an advanced enough Distutils.
import distutils
from distutils.core import setup, Extension

setup (name = "cluster",
       maintainer = "SciPy Developers",
       author = "eric jones",
       maintainer_email = "scipy-devel@scipy.org",
       description = "Clustering Algorithms (Information Theory)",
       url = "http://www.scipy.org",

       packages = ['cluster'],
       package_dir = {'cluster':'.'},
       include_dirs = ['src'],
       ext_modules =  [ Extension('cluster._vq',['src/vq_wrap.cpp']) ]
       )
