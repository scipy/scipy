cmake = import('cmake')
opt_var = cmake.subproject_options()

opt_var.add_cmake_defines(
  {'EIGEN_BUILD_DOC': false,
   'EIGEN_BUILD_BTL': false,
   'EIGEN_BUILD_BLAS': false,
   'EIGEN_BUILD_LAPACK': false,
   'EIGEN_BUILD_TESTING': false,
   'EIGEN_BUILD_PKGCONFIG': false,
   'EIGEN_LEAVE_TEST_IN_ALL_TARGET': false})

# Configure Eigen as a CMake subproject from the git submodule.
eigen_sp = cmake.subproject('eigen', options: opt_var)

message('CMake targets:\n - ' + '\n - '.join(eigen_sp.target_list()))
# Eigen’s exported target is usually Eigen3::Eigen
eigen3_dep = eigen_sp.dependency('eigen', include_type: 'system')

# Riccaticpp is header-only. Its CMakeLists.txt would FetchContent-pull
# a second eigen (3.4.0) and try to build its Fortran BLAS/LAPACK,
# which triggers a meson+cmake gfortran '-J' issue and is wasted work
# anyway. Bypass CMake entirely by referencing the headers directly.
riccaticpp_dep = declare_dependency(
  include_directories: include_directories('../../../subprojects/riccaticpp/include'))
riccaticpp_pymain_sources = files('_pymain.cpp')

eigen_deps = [np_dep, pybind11_dep, eigen3_dep, riccaticpp_dep]

py3.extension_module(
  'pyriccaticpp',
  riccaticpp_pymain_sources,
  dependencies: eigen_deps,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate/_ivp',
  cpp_args: ['-std=c++17'],
)

py3.install_sources([
    '__init__.py',
    '_osc.py',
    'base.py',
    'bdf.py',
    'common.py',
    'dop853_coefficients.py',
    'ivp.py',
    'lsoda.py',
    'radau.py',
    'rk.py'
  ],
  subdir: 'scipy/integrate/_ivp'
)

subdir('tests')
