project('repro_gh_11577', 'c')

openblas_dep = dependency('openblas')

executable('repro_c',
    'ggev_repro_gh_11577.c',
    dependencies: openblas_dep
)
