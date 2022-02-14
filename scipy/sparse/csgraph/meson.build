pyx_files = [
  ['_flow', '_flow.pyx'],
  ['_matching', '_matching.pyx'],
  ['_min_spanning_tree', '_min_spanning_tree.pyx'],
  ['_reordering', '_reordering.pyx'],
  ['_shortest_path', '_shortest_path.pyx'],
  ['_tools', '_tools.pyx'],
  ['_traversal', '_traversal.pyx']
]

foreach pyx_file: pyx_files
  py3.extension_module(pyx_file[0],
    cython_gen.process(pyx_file[1]),
    c_args: cython_c_args,
    include_directories: inc_np,
    dependencies: py3_dep,
    install: true,
    subdir: 'scipy/sparse/csgraph'
  )
endforeach

python_sources = [
  '__init__.py',
  '_laplacian.py',
  '_validation.py',
  'setup.py'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/sparse/csgraph'
)

subdir('tests')
