#!/usr/bin/env python

from setup import setup_package

ignore_packages = ['chaco','kiva','freetype','traits']

if __name__ == "__main__":
    setup_package(ignore_packages)
