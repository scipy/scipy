def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.compiler_helper import set_cxx_flags_hook

    config = Configuration('csgraph', parent_package, top_path)

    config.add_data_dir('tests')

    ext = config.add_extension('_shortest_path',
         sources=['_shortest_path.cxx'],
         include_dirs=[numpy.get_include()])
    ext._pre_build_hook = set_cxx_flags_hook

    config.add_extension('_traversal',
         sources=['_traversal.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_min_spanning_tree',
         sources=['_min_spanning_tree.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_matching',
         sources=['_matching.c'],
         include_dirs=[numpy.get_include()])
    
    config.add_extension('_flow',
         sources=['_flow.c'],
         include_dirs=[numpy.get_include()])
    
    config.add_extension('_reordering',
         sources=['_reordering.c'],
         include_dirs=[numpy.get_include()])

    config.add_extension('_tools',
         sources=['_tools.c'],
         include_dirs=[numpy.get_include()])

    return config
