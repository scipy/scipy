import os
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

# Define common include directories
include_dirs = [get_numpy_include_dirs()]

# List of extensions to be added
extensions = [
    ('_vq', ['_vq.c']),
    ('_hierarchy', ['_hierarchy.c']),
    ('_optimal_leaf_ordering', ['_optimal_leaf_ordering.c'])
]

# Define macros
DEFINE_MACROS = [("SCIPY_PY3K", None)]

def configuration(parent_package='', top_path=None):
    config = Configuration('cluster', parent_package, top_path)

    # Add data directory
    config.add_data_dir('tests')

    # Add extensions using a loop
    for ext_name, sources in extensions:
        config.add_extension(ext_name, sources=sources, include_dirs=include_dirs, define_macros=DEFINE_MACROS)

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
