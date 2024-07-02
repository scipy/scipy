pyx_files = [
  ['_flow', '_flow.pyx'],
  ['_matching', '_matching.pyx'],
  ['_min_spanning_tree', '_min_spanning_tree.pyx'],
  ['_reordering', '_reordering.pyx'],
  ['_tools', '_tools.pyx'],
  ['_traversal', '_traversal.pyx']
]

cython_gen_csgraph = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, fs.copyfile('parameters.pxi')],
)

foreach pyx_file: pyx_files
  py3.extension_module(pyx_file[0],
    cython_gen_csgraph.process(pyx_file[1]),
    c_args: cython_c_args,
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/sparse/csgraph'
  )
endforeach

py3.extension_module('_shortest_path',
  cython_gen_cpp.process('_shortest_path.pyx'),
  cpp_args: cython_cpp_args,
  include_directories: inc_np,
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/sparse/csgraph'
)

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
