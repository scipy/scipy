project('ilp64_test_package', 'c', 'cython')

fs = import('fs')
py = import('python').find_installation(pure: false)
cy = meson.get_compiler('cython')

if not cy.version().version_compare('>=3.2.0')
  error('Requires Cython >= 3.2.0')
endif

cython_args = [
  '-Xfreethreading_compatible=True',
  '-I' + meson.project_build_root() / 'src' / 'ilp64_test_package',
]

subdir('src/ilp64_test_package')
