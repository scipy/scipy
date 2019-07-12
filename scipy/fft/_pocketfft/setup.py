
def try_compile(compiler, code=None, flags=[], ext='.cpp'):
    """Returns True if the compiler is able to compile the given code"""
    import tempfile
    from distutils.errors import CompileError
    import os

    code = code or 'int main (int argc, char **argv) { return 0; }'

    with tempfile.TemporaryDirectory() as temp_dir:
        fname = os.path.join(temp_dir, 'main'+ext)
        with open(fname, 'w') as f:
            f.write(code)

        try:
            compiler.compile([fname], output_dir=temp_dir, extra_postargs=flags)
        except CompileError:
            return False
    return True


def has_flag(compiler, flag):
    return try_compile(compiler, flags=[flag])


def get_std_flag(compiler):
    # Test the compiler for the highest available c++ standard flag
    gnu_flags = ['--std=c++14', '--std=c++11']
    flags_by_cc = {
        'msvc': ['/std:c++14', None],
        'intelw': ['/Qstd=c++14', '/Qstd=c++11']
    }
    flags = flags_by_cc.get(compiler.compiler_type, gnu_flags)

    for flag in flags:
        if flag is None:
            return None

        if has_flag(compiler, flag):
            return flag

    from numpy.distutils import log
    log.warn('Could not detect c++ standard flag')
    return None


def try_add_flag(args, compiler, flag):
    """Appends flag to the list of arguments if supported by the compiler"""
    if try_compile(compiler, flags=args+[flag]):
        args.append(flag)


def pre_build_hook(build_ext, ext):
    cc = build_ext._cxx_compiler
    args = ext.extra_compile_args

    std_flag = get_std_flag(build_ext._cxx_compiler)
    if std_flag is not None:
        args.append(std_flag)

    if cc.compiler_type == 'msvc':
        args.append('/EHsc')
    else:
        try_add_flag(args, cc, '-fvisibility=hidden')

        import sys
        if sys.platform == 'darwin':
            args.append('-mmacosx-version-min=10.7')
            try_add_flag(args, cc, '-stdlib=libc++')


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import pybind11
    include_dirs = [pybind11.get_include(True), pybind11.get_include(False)]

    config = Configuration('_pocketfft', parent_package, top_path)
    ext = config.add_extension('pypocketfft',
                               sources=['pypocketfft.cxx'],
                               depends=['pocketfft_hdronly.h'],
                               include_dirs=include_dirs,
                               language='c++')
    ext._pre_build_hook = pre_build_hook

    config.add_data_files('LICENSE.md')
    config.add_data_dir('tests')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
