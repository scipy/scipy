#!/usr/bin/env python

from setup import setup_package

ignore_packages = ['chaco','kiva','freetype','traits']
if sys.platform=='win32':
    ignore_packages.append('xplt')

if __name__ == "__main__":
    setup_package(ignore_packages)
