pyx_files_for_c = [
  ['_flow', '_flow.pyx'],
  ['_matching', '_matching.pyx'],
  ['_min_spanning_tree', '_min_spanning_tree.pyx'],
  ['_reordering', '_reordering.pyx'],
  ['_tools', '_tools.pyx'],
  ['_traversal', '_traversal.pyx']
]
pyx_files_for_cpp = [
  ['_shortest_path', '_shortest_path.pyx'],
]

cython_gen_csgraph_for_c = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, fs.copyfile('parameters.pxi')],
)

foreach pyx_file: pyx_files_for_c
  py3.extension_module(pyx_file[0],
    cython_gen_csgraph_for_c.process(pyx_file[1]),
    c_args: cython_c_args,
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/sparse/csgraph',
  )
endforeach

cython_gen_csgraph_for_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree],
)

foreach pyx_file: pyx_files_for_cpp
  py3.extension_module(pyx_file[0],
    cython_gen_csgraph_for_cpp.process(pyx_file[1]),
    cpp_args: cython_cpp_args,
    include_directories: inc_np,
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/sparse/csgraph',
  )
endforeach

python_sources = [
  '__init__.py',
  '_laplacian.py',
  '_validation.py',
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/sparse/csgraph'
)

subdir('tests')
