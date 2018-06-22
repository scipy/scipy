"""
Find scipy functions whose docstrings do not contain "Examples".

This looks for *functions* only, not all callables.
"""

from __future__ import print_function

import importlib
import types
import numpy as np
import scipy

modules = ['cluster.hierarchy', 'cluster.vq', 'constants', 'fftpack',
           'integrate', 'interpolate', 'io', 'linalg', 'misc', 'ndimage',
           'odr', 'optimize', 'signal', 'sparse', 'sparse.linalg',
           'sparse.csgraph', 'spatial', 'special', 'stats',
           'stats.contingency', 'stats.mstats']

skip = ['misc.lena']

print("scipy version %s" % scipy.__version__)
print()

total = 0
for module_name in modules:
    mod = importlib.import_module('.' + module_name, package='scipy')
    objects = [(name, getattr(mod, name))
                   for name in getattr(mod, '__all__', dir(mod))
                       if not name.startswith('_')]
    funcs = [item for item in objects
                      if isinstance(item[1], (types.FunctionType,
                                              types.BuiltinFunctionType,
                                              np.ufunc))]
    noex = [item for item in funcs
                     if ((module_name + '.' + item[0]) not in skip and
                         (item[1].__doc__ is None or
                          ("is deprecated" not in item[1].__doc__) and
                          ("Examples" not in item[1].__doc__)))]
    if len(noex) > 0:
        total += len(noex)
        print(module_name, "(%d)" % len(noex))
        for name, func in noex:
            print("   ", name, end='')
            if func.__doc__ is None:
                print(" \t[no docstring]")
            else:
                print()

print()
