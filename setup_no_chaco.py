#!/usr/bin/env python

import sys
from setup import setup_package

ignore_packages = ['chaco','kiva','freetype','traits']
if sys.platform=='win32':
    ignore_packages.append('xplt')
    ignore_packages.append('sparse')

if __name__ == "__main__":
    setup_package(ignore_packages)
