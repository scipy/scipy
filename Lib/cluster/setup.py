#!/usr/bin/env python
import os, sys, string, re
from glob import glob

# Check for an advanced enough Distutils.
import distutils
from distutils.core import setup, Extension

setup (name = "_vq",
       include_dirs = ['src'],
       ext_modules =  [ Extension('_vq',['src/vq_wrap.cpp']) ]
       )