"""
Show which shared libraries are loaded by numpy and scipy imports.

Reads /proc/self/maps before and after each import to determine which .so
files were pulled in. Libraries are categorized into:

- Shared libraries: "real" shared libs (BLAS, compilers, system libs)
- Stdlib extensions: CPython extension modules from lib-dynload/
- numpy/scipy extensions: extension modules from site-packages or build dirs

By default, only shared libraries are shown. Use -s and -e to include the
other categories.

Usage::

    # Via spin check:
    spin check --loaded-sharedlibs

    # Via pixi task (builds automatically if needed):
    pixi run check-mkl-lp64-sharedlibs
    pixi run check-mkl-ilp64-sharedlibs

    # Direct invocation with extra options:
    spin python -- tools/verify_loaded_sharedlibs.py [-s] [-e] [IMPORT ...]

Positional arguments (direct invocation only)::

    IMPORT  Dotted import names to check (default: numpy). Each import is
            performed in order; new shared libraries are reported per import.

Options (direct invocation only)::

    -s, --show-stdlib      Show Python stdlib extension modules (lib-dynload/)
    -e, --show-extensions  Show numpy/scipy extension modules

Examples::

    # Check with MKL LP64 build:
    pixi run -e mkl-lp64 spin check --no-build --build-dir=build-mkl-lp64 \
            --loaded-sharedlibs

    # Or use the pixi task shorthand:
    pixi run check-mkl-lp64-sharedlibs

    # Direct invocation for numpy then scipy.linalg:
    spin python -- tools/verify_loaded_sharedlibs.py numpy scipy.linalg

    # Single scipy submodule:
    spin python -- tools/verify_loaded_sharedlibs.py scipy.sparse

Sample output from ``pixi run check-mkl-lp64-sharedlibs numpy scipy.linalg``::

    numpy pulled in:
      .pixi/envs/mkl-lp64/lib/libffi.so.8.2.0
      .pixi/envs/mkl-lp64/lib/libgcc_s.so.1
      .pixi/envs/mkl-lp64/lib/libmkl_avx512.so.2
      .pixi/envs/mkl-lp64/lib/libmkl_core.so.2
      .pixi/envs/mkl-lp64/lib/libmkl_gf_lp64.so.2
      .pixi/envs/mkl-lp64/lib/libmkl_intel_thread.so.2
      .pixi/envs/mkl-lp64/lib/libmkl_rt.so.2
      .pixi/envs/mkl-lp64/lib/libmkl_vml_avx512.so.2
      .pixi/envs/mkl-lp64/lib/libomp.so
      .pixi/envs/mkl-lp64/lib/libstdc++.so.6.0.34
      /usr/lib/librt.so.1
      (4 stdlib + 2 numpy/scipy extension modules hidden)

    scipy.linalg pulled in:
      .pixi/envs/build-mkl-lp64/lib/libmkl_intel_lp64.so.2
      .pixi/envs/build-mkl-lp64/lib/libmkl_sequential.so.2
      .pixi/envs/mkl-lp64/lib/libcrypto.so.3
      (14 stdlib + 23 numpy/scipy extension modules hidden)
"""
import argparse
import importlib
import re


def loaded_libs():
    with open("/proc/self/maps") as f:
        return {
            line.split()[-1]
            for line in f
            if ".so" in line and line.split()[-1].startswith("/")
        }


def strip_prefix(path):
    """Strip the repo root prefix, keeping paths starting from .pixi/ or build-*."""
    m = re.search(r'(\.pixi/|build-)', path)
    if m:
        return path[m.start():]
    return path


def classify_libs(libs):
    """Classify libraries into shared libs, stdlib extensions, and package extensions"""
    shared = []
    stdlib = []
    extensions = []
    for lib in libs:
        if '/lib-dynload/' in lib:
            stdlib.append(lib)
        elif '/site-packages/' in lib or re.search(r'/build-[^/]*-install/', lib):
            extensions.append(lib)
        else:
            shared.append(lib)
    return {
        'shared': sorted(shared, key=strip_prefix),
        'stdlib': sorted(stdlib, key=strip_prefix),
        'extensions': sorted(extensions, key=strip_prefix),
    }


def print_section(name, libs, args):
    """Print one import section with categorized libraries."""
    categorized = classify_libs(libs)

    print(f"\n{name} pulled in:")
    for lib in categorized['shared']:
        print(f"  {strip_prefix(lib)}")

    if args.show_stdlib and categorized['stdlib']:
        print(f"\n  stdlib extensions ({len(categorized['stdlib'])}):")
        for lib in categorized['stdlib']:
            print(f"    {strip_prefix(lib)}")

    if args.show_extensions and categorized['extensions']:
        print(f"\n  numpy/scipy extensions ({len(categorized['extensions'])}):")
        for lib in categorized['extensions']:
            print(f"    {strip_prefix(lib)}")

    hidden = []
    if not args.show_stdlib and categorized['stdlib']:
        hidden.append(f"{len(categorized['stdlib'])} stdlib")
    if not args.show_extensions and categorized['extensions']:
        hidden.append(f"{len(categorized['extensions'])} numpy/scipy")
    if hidden:
        print(f"  ({' + '.join(hidden)} extension modules hidden)")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Show shared libraries loaded by numpy and scipy imports."
    )
    parser.add_argument(
        '-s', '--show-stdlib', action='store_true',
        help="Show Python stdlib extension modules (lib-dynload/)"
    )
    parser.add_argument(
        '-e', '--show-extensions', action='store_true',
        help="Show numpy/scipy extension modules"
    )
    parser.add_argument(
        'imports', nargs='*', default=['numpy'], metavar='IMPORT',
        help="Dotted import names to check (default: numpy)"
    )
    args = parser.parse_args()

    prev = loaded_libs()
    for name in args.imports:
        importlib.import_module(name)
        curr = loaded_libs()
        print_section(name, curr - prev, args)
        prev = curr


if __name__ == "__main__":
    main()
