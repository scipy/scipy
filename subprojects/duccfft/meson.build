project('duccfft', version: '0.40.0')

duccfft_sources = files(
    'ducc0/infra/mav.cc',
    'ducc0/infra/string_utils.cc',
    'ducc0/infra/threading.cc',
)

duccfft_dep = declare_dependency(
    sources: duccfft_sources,
    include_directories: '.',
)
