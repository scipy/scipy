import re
import os

import numpy as np

__all__ = ['blas_ilp64_pre_build_hook', 'generic_pre_build_hook',
           'write_file_content']


# See _fflag_lp64/ilp64 in scipy/meson.build for these
# def get_fcompiler_ilp64_flags():
# def get_fcompiler_macro_include_flags(path): 

# Done in scipy/meson.build
# def get_f2py_int64_options():


def _blas_ilp64_pre_build_hook(outdir, prefix='', suffix=''):
    if suffix:
        if not suffix.endswith('_'):
            # Symbol suffix has to end with '_' to be Fortran-compatible
            raise RuntimeError("BLAS/LAPACK has incompatible symbol suffix: "
                               "{!r}".format(suffix))

        suffix = suffix[:-1]

    # When symbol prefix/suffix is present, we have to patch sources

    # Create name-mapping include files
    include_name_f = 'blas64-prefix-defines.inc'
    include_name_c = 'blas64-prefix-defines.h'
    include_fn_f = os.path.join(outdir, include_name_f)
    include_fn_c = os.path.join(outdir, include_name_c)

    text = ""
    for symbol in get_blas_lapack_symbols():
        text += f'#define {symbol} {prefix}{symbol}_{suffix}\n'
        text += f'#define {symbol.upper()} {prefix}{symbol}_{suffix}\n'

        # Code generation may give source codes with mixed-case names
        for j in (1, 2):
            s = symbol[:j].lower() + symbol[j:].upper()
            text += f'#define {s} {prefix}{symbol}_{suffix}\n'
            s = symbol[:j].upper() + symbol[j:].lower()
            text += f'#define {s} {prefix}{symbol}_{suffix}\n'

    write_file_content(include_fn_f, text)

    ctext = re.sub(r'^#define (.*) (.*)$', r'#define \1_ \2_', text, flags=re.M)
    write_file_content(include_fn_c, text + "\n" + ctext)

    # Patch sources to include it
    def patch_source(filename, old_text):
        text = f'#include "{include_name_f}"\n'
        text += old_text
        return text

    # TODO: update source file names
    #return generic_pre_build_hook(patch_source_func=patch_source,
    #                             source_fnpart="_blas64")


def generic_pre_build_hook(cmd, ext, fcompiler_flags, patch_source_func=None,
                           source_fnpart=None):
    """
    Pre-build hook for adding compiler flags and patching sources.

    Parameters
    ----------
    cmd : distutils.core.Command
        Hook input. Current distutils command (build_clib or build_ext).
    ext : dict or numpy.distutils.extension.Extension
        Hook input. Configuration information for library (dict, build_clib)
        or extension (numpy.distutils.extension.Extension, build_ext).
    fcompiler_flags : dict
        Dictionary of ``{'compiler_name': ['-flag1', ...]}`` containing
        compiler flags to set.
    patch_source_func : callable, optional
        Function patching sources, see `_generic_patch_sources` below.
    source_fnpart : str, optional
        String to append to the modified file basename before extension.

    """
    # Mangle sources
    if patch_source_func is not None:
        build_info.setdefault('depends', []).extend(build_info['sources'])
        new_sources = _generic_patch_sources(build_info['sources'], patch_source_func,
                                             source_fnpart)
        build_info['sources'][:] = new_sources



def _generic_patch_sources(filenames, patch_source_func, source_fnpart, root_dir=None):
    """
    Patch Fortran sources, creating new source files.

    Parameters
    ----------
    filenames : list
        List of Fortran source files to patch.
        Files not ending in ``.f`` or ``.f90`` are left unaltered.
    patch_source_func : callable(filename, old_contents) -> new_contents
        Function to apply to file contents, returning new file contents
        as a string.
    source_fnpart : str
        String to append to the modified file basename before extension.
    root_dir : str, optional
        Source root directory. Default: cwd

    Returns
    -------
    new_filenames : list
        List of names of the newly created patched sources.

    """
    new_filenames = []

    if root_dir is None:
        root_dir = os.getcwd()

    root_dir = os.path.abspath(root_dir)
    src_dir = os.path.join(root_dir, _get_build_src_dir())

    for src in filenames:
        base, ext = os.path.splitext(os.path.basename(src))

        if ext not in ('.f', '.f90'):
            new_filenames.append(src)
            continue

        with open(src) as fsrc:
            text = patch_source_func(src, fsrc.read())

        # Generate useful target directory name under src_dir
        src_path = os.path.abspath(os.path.dirname(src))

        for basedir in [src_dir, root_dir]:
            if os.path.commonpath([src_path, basedir]) == basedir:
                rel_path = os.path.relpath(src_path, basedir)
                break
        else:
            raise ValueError(f"{src!r} not under {root_dir!r}")

        dst = os.path.join(src_dir, rel_path, base + source_fnpart + ext)
        write_file_content(dst, text)

        new_filenames.append(dst)

    return new_filenames


def write_file_content(filename, content):
    """
    Write content to file, but only if it differs from the current one.
    """
    if os.path.isfile(filename):
        with open(filename) as f:
            old_content = f.read()

        if old_content == content:
            return

    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    with open(filename, 'w') as f:
        f.write(content)


def get_blas_lapack_symbols():
    cached = getattr(get_blas_lapack_symbols, 'cached', None)
    if cached is not None:
        return cached

    # Obtain symbol list from Cython Blas/Lapack interface
    srcdir = os.path.join(os.path.dirname(__file__), os.pardir, 'linalg')

    symbols = []

    # Get symbols from the generated files
    for fn in ['cython_blas_signatures.txt', 'cython_lapack_signatures.txt']:
        with open(os.path.join(srcdir, fn)) as f:
            for line in f:
                m = re.match(r"^\s*[a-z]+\s+([a-z0-9]+)\(", line)
                if m:
                    symbols.append(m.group(1))

    # Get the rest from the generator script
    # (we cannot import it directly here, so use exec)
    sig_fn = os.path.join(srcdir, '_cython_signature_generator.py')
    with open(sig_fn) as f:
        code = f.read()
    ns = {'__name__': '<module>'}
    exec(code, ns)
    symbols.extend(ns['blas_exclusions'])
    symbols.extend(ns['lapack_exclusions'])

    get_blas_lapack_symbols.cached = tuple(sorted(set(symbols)))
    return get_blas_lapack_symbols.cached
