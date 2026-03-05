# Generate _version.py

uarray_dir = py.get_install_dir() / 'uarray'
meson.add_dist_script(
   ['../../tools/gitversion.py']
)
if not fs.exists('_version.py')
  custom_target(
    'generate-version',
    install: true,
    build_always_stale: true,
    build_by_default: true,
    output: '_version.py',
    input: '../../tools/gitversion.py',
    command: ['../../tools/gitversion.py'],
    install_dir: uarray_dir,
    install_tag: 'python-runtime',
  )
endif

# Installation #

if not meson.is_subproject()
  foreach prefix, files : python_sources
    py.install_sources(files, subdir: prefix)
  endforeach
endif

py.install_sources([
    'conftest.py',
  ],
  subdir: 'uarray',
  install_tag: 'tests',
)

subdir('tests')
