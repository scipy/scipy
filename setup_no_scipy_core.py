#!/usr/bin/env python

import setup

if __name__ == '__main__':
    ignore_packages = ['scipy_test','scipy_base','scipy_distutils']
    setup.setup_package(ignore_packages)
