python_sources = uarray.get_variable('python_sources')
cpp_sources = uarray.get_variable('cpp_sources')
include_dirs = uarray.get_variable('include_dirs')

foreach prefix, files : python_sources
  py3.install_sources(files, subdir: external_dir / prefix)
endforeach

py3.extension_module('_uarray',
  cpp_sources,
  cpp_args: [_cpp_Wno_terminate, _cpp_Wno_unused_function],
  link_args: version_link_args,
  include_directories: include_dirs,
  install: true,
  subdir: external_dir / 'uarray',
)
