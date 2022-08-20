"""
Find scipy functions and methods whose docstrings contain "Examples" and
use 'np.' but do not have 'import numpy as np'.
"""

import importlib
import types
import numpy as np
import scipy
from scipy._lib.uarray import _Function


def is_missing_import_np(docstring):
    if docstring is None:
        return False
    examples_start = docstring.find('Examples\n')
    if examples_start == -1:
        return False
    examples_section = docstring[examples_start:]
    return ('np.' in examples_section and
            'import numpy as np' not in examples_section)


modules = ['cluster.hierarchy', 'cluster.vq', 'constants', 'fft', 'fftpack',
           'integrate', 'interpolate', 'io', 'io.arff', 'io.wavfile',
           'linalg', 'misc', 'ndimage', 'odr', 'optimize',
           'signal', 'sparse', 'sparse.linalg', 'sparse.csgraph', 'spatial',
           'special', 'stats', 'stats.contingency', 'stats.mstats']

skip = []

print(f"scipy version {scipy.__version__}")
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
                                     np.ufunc,
                                     _Function))]
    no_np = [item for item in funcs
             if (((module_name + '.' + item[0]) not in skip)
                 and is_missing_import_np(item[1].__doc__))]

    method_no_np = []
    classes = [item for item in objects
               if isinstance(item[1], type)]
    for cls_item in classes:
        name, cls = cls_item
        for cls_attr in dir(cls):
            cls_obj = getattr(cls, cls_attr)
            if (callable(cls_obj)
                    and not cls_attr.startswith('_')
                    and not isinstance(cls_obj, types.MemberDescriptorType)
                    and is_missing_import_np(cls_obj.__doc__)):
                method_no_np.append((name, cls_attr))

    num_found = len(no_np) + len(method_no_np)
    if num_found > 0:
        total += num_found
        print(module_name, f"({num_found})")
        no_np.sort()
        for name, func in no_np:
            print("   ", name,)

        prev_name = None
        for name, cls_attr in method_no_np:
            if name != prev_name:
                print(f'    {name} (class)')
                prev_name = name
            print(f'        .{cls_attr}')

print()
print(f"Found {total} objects missing 'import numpy as np'")