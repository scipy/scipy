project('ducc', 'cpp', version: '0.38.0')

duccfft_sources = files(
    'ducc/src/ducc0/infra/mav.cc',
    'ducc/src/ducc0/infra/string_utils.cc',
    'ducc/src/ducc0/infra/threading.cc',
)

duccfft_dep = declare_dependency(
    sources: duccfft_sources,
    include_directories: 'ducc/src',
)
