project(
  'uarray',
  'cpp',
)

# NOTE: downstream projects may rely on this variable name when vendoring,
# do not change it without consulting with downstream projects
python_sources = {
  'uarray': files(
    'uarray/uarray/__init__.py',
    'uarray/uarray/_backend.py',
    'uarray/uarray/_typing.pyi',
    'uarray/uarray/_typing.pyi',
    'uarray/uarray/conftest.py',
  )
}

cpp_sources = files(
  'uarray/_uarray_dispatch.cxx',
  'uarray/vectorcall.cxx',
)

include_dirs = include_directories('uarray')
