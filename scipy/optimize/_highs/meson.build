fs = import('fs')
if not fs.exists('../../../subprojects/highs/README.md')
    error(
        'Missing the `highs` submodule! Run `git submodule update --init` to fix this.',
    )
endif

highs_opts = {
    'default_library': 'static',
    'with_pybind11': false, # We don't want to build and install highspy
    'use_zlib': 'disabled',
}

highs_proj = subproject('highs', default_options: highs_opts)

highs_dep = highs_proj.get_variable('highs_dep')
highspy_cpp = highs_proj.get_variable('highspy_cpp')
highs_conf_file = highs_proj.get_variable('highs_conf_file')
highsoptions_cpp = highs_proj.get_variable('highsoptions_cpp')

scipy_highspy_dep = [py3_dep, pybind11_dep, highs_dep, thread_dep, atomic_dep]

py3.install_sources(
    ['__init__.py', '_highs_wrapper.py'],
    subdir: 'scipy/optimize/_highspy',
)

py3.extension_module(
    '_core',
    sources: highspy_cpp + highs_conf_file,
    dependencies: scipy_highspy_dep,
    install: true,
    subdir: 'scipy/optimize/_highspy',
)

py3.extension_module(
    '_highs_options',
    sources: highsoptions_cpp + highs_conf_file,
    dependencies: scipy_highspy_dep,
    install: true,
    subdir: 'scipy/optimize/_highspy',
)
