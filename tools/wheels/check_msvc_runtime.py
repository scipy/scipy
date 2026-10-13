#!/usr/bin/env python
"""
check_msvc_runtime.py [MODULE]

Check that a Windows wheel vendors the MSVC C++ runtime.

The extension modules are built against the MSVC C++ standard library. Unlike
vcruntime140.dll, msvcp140.dll is not part of Windows and is not shipped with
CPython, so delvewheel vendors it into ``scipy.libs`` under a mangled name.
When that silently stops happening, importing SciPy can fail with "DLL load
failed" on machines without the VC++ redistributable installed - see
scipy/scipy#14998 and scipy/scipy#17191.

This is a no-op on non-Windows platforms.

"""
import pathlib
import sys
import argparse


def main():
    p = argparse.ArgumentParser(usage=__doc__.rstrip())
    p.add_argument("module", nargs="?", default="scipy")
    args = p.parse_args()

    if sys.platform != "win32":
        sys.exit(0)

    # Drop '' from sys.path
    sys.path.pop(0)

    # Find module path
    __import__(args.module)
    mod = sys.modules[args.module]

    # The vendored DLLs are installed alongside the package, in <module>.libs
    libs_path = pathlib.Path(mod.__file__).parent.parent / f"{args.module}.libs"

    vendored = sorted(p.name for p in libs_path.glob("msvcp140*.dll"))
    if not vendored:
        print(
            f"ERROR: no msvcp140*.dll found in {libs_path}. The MSVC C++ "
            "runtime is not being vendored, so importing this wheel may fail "
            "on machines without the VC++ redistributable installed.\n"
        )
        sys.exit(1)

    print(f"Found vendored MSVC C++ runtime: {vendored}")
    sys.exit(0)


if __name__ == "__main__":
    main()
