fs = import('fs')
if not fs.exists('../../../subprojects/highs/README.md')
    error(
        'Missing the `highs` submodule! Run `git submodule update --init` to fix this.',
    )
endif

highs_opts = {
    'default_library': 'static',
    'with_pybind11': false,  # We don't want to build and install highspy
    'use_zlib': 'disabled',
    'do_install': false,
}

highs_proj = subproject('highs', default_options: highs_opts)

highs_args = highs_proj.get_variable('_args')
highs_deps = highs_proj.get_variable('_deps')
highslib_srcs = highs_proj.get_variable('highslib_srcs')
highs_incdirs = highs_proj.get_variable('_incdirs')

_highslib = static_library(
    'highs',
    highslib_srcs,
    dependencies: highs_deps,
    cpp_args: highs_args,
    c_args: highs_args,
    include_directories: highs_incdirs,
    # Ensure that if we link a static library into a shared library,
    # private symbols don't get re-exported.
    gnu_symbol_visibility: 'inlineshidden',
    pic: true,
    install: true,
    install_dir: py3.get_install_dir() / 'scipy/optimize/_highspy',
)

_highs_dep = declare_dependency(
    link_with: _highslib,
    dependencies: highs_deps,
    include_directories: highs_incdirs,
)

highspy_cpp = highs_proj.get_variable('highspy_cpp')
highs_conf_file = highs_proj.get_variable('highs_conf_file')

scipy_highspy_dep = [py3_dep, pybind11_dep, _highs_dep, thread_dep, atomic_dep]

py3.install_sources(
    ['__init__.py', '_highs_wrapper.py'],
    subdir: 'scipy/optimize/_highspy',
)

py3.extension_module(
    '_core',
    sources: highspy_cpp + highs_conf_file,
    dependencies: scipy_highspy_dep,
    c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
    cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
    link_args: version_link_args,
    gnu_symbol_visibility: 'hidden',  # default
    install: true,
    subdir: 'scipy/optimize/_highspy',
)

py3.extension_module(
    '_highs_options',
    sources: [highs_conf_file, 'highs_options.cpp'],
    dependencies: scipy_highspy_dep,
    c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
    cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
    link_args: version_link_args,
    gnu_symbol_visibility: 'hidden',  # default
    install: true,
    subdir: 'scipy/optimize/_highspy',
)
