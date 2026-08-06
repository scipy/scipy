project('qhull_r',
  'c',
  version : '8.0.2',
  meson_version: '>= 1.5.0',
)

cc = meson.get_compiler('c')

libqhull_r_sources = [
  'libqhull_r/geom2_r.c',
  'libqhull_r/geom_r.c',
  'libqhull_r/global_r.c',
  'libqhull_r/io_r.c',
  'libqhull_r/libqhull_r.c',
  'libqhull_r/mem_r.c',
  'libqhull_r/merge_r.c',
  'libqhull_r/poly2_r.c',
  'libqhull_r/poly_r.c',
  'libqhull_r/qset_r.c',
  'libqhull_r/random_r.c',
  'libqhull_r/rboxlib_r.c',
  'libqhull_r/stat_r.c',
  'libqhull_r/user_r.c',
  'libqhull_r/usermem_r.c',
  'libqhull_r/userprintf_r.c',
  'libqhull_r/userprintf_rbox_r.c'
]

libqhull_r = static_library(
    'libqhull_r',
    [libqhull_r_sources],
    c_args: cc.get_supported_arguments('-Wno-unused-but-set-variable'),
    include_directories: 'libqhull_r',
    # Ensure that if we link a static library into a shared library,
    # private symbols don't get re-exported.
    gnu_symbol_visibility: 'inlineshidden',
    install: false,
)

qhull_r_dep = declare_dependency(
    link_with: libqhull_r,
    include_directories: '.',
)
