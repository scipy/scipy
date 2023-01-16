
def pre_build_hook(build_ext, ext):
    from scipy._build_utils.compiler_helper import get_cxx_std_flag
    std_flag = get_cxx_std_flag(build_ext._cxx_compiler)
    if std_flag is not None:
        ext.extra_compile_args.append(std_flag)


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('csgraph', parent_package, top_path)

    config.add_data_dir('tests')

    ext = config.add_extension('_shortest_path',
                         sources=['_shortest_path.cxx'],
                         include_dirs=[numpy.get_include()],
                         language='c++')
    ext._pre_build_hook = pre_build_hook

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
