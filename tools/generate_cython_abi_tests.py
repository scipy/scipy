#!/usr/bin/env python3
"""
generate_cython_abi_tests.py

Run against a known-good SciPy release to generate a regression test that
asserts the Cython public API binary interfaces are unchanged in future builds.

Usage:
    python generate_cython_abi_tests.py linalg
    python generate_cython_abi_tests.py special
    python generate_cython_abi_tests.py optimize

This writes two files per submodule:
    scipy/<submodule>/tests/test_cython_abi.py        — the test code
    scipy/<submodule>/tests/cython_abi_signatures.json — expected signatures
"""

import argparse
import ctypes
import json
import sys
import textwrap
from datetime import date
from pathlib import Path

# ---------------------------------------------------------------------------
# Submodule → Cython module mapping
# ---------------------------------------------------------------------------

SUBMODULE_MAP = {
    "linalg": {
        "modules": [
            ("scipy.linalg.cython_blas", "cython_blas"),
            ("scipy.linalg.cython_lapack", "cython_lapack"),
        ],
        # The cython_blas/cython_lapack ABI depends on whether ILP64 BLAS
        # is used, so skip when it doesn't match what the signatures were
        # generated against.
        "skip_condition": "CYTHON_BLAS_ILP64",
        "skip_reason": "cython_blas/cython_lapack ABI depends on ILP64 BLAS",
    },
    "special": {
        "modules": [
            ("scipy.special.cython_special", "cython_special"),
        ],
    },
    "optimize": {
        "modules": [
            ("scipy.optimize.cython_optimize._zeros", "cython_optimize"),
        ],
    },
}

# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------

_get_name = ctypes.pythonapi.PyCapsule_GetName
_get_name.restype = ctypes.c_char_p
_get_name.argtypes = [ctypes.py_object]


def extract_capi(module):
    """
    Return a sorted dict {function_name: capsule_signature_string} for every
    entry in module.__pyx_capi__.

    Entries for which PyCapsule_GetName returns NULL (should never happen for
    Cython-generated capsules, but be defensive) are omitted with a warning.
    """
    result = {}
    for name, capsule in module.__pyx_capi__.items():
        raw = _get_name(capsule)
        if raw is None:
            print(f"WARNING: PyCapsule_GetName returned NULL for {name!r}",
                  file=sys.stderr)
            continue
        result[name] = raw.decode("utf-8")
    return dict(sorted(result.items()))


# ---------------------------------------------------------------------------
# Code generation
# ---------------------------------------------------------------------------

def generate(submodule, config):
    """
    Generate a test file and JSON signatures file for the given submodule.

    submodule: name of the scipy submodule (e.g. "linalg", "special")
    config: dict with keys "modules" (list of (import_path, attr) pairs),
            and optional "skip_condition" / "skip_reason"
    """
    import importlib
    import scipy

    modules = config["modules"]
    skip_condition = config.get("skip_condition")
    skip_reason = config.get("skip_reason", "")

    # Extract signatures from all modules
    signatures = {}
    for import_path, attr in modules:
        mod = importlib.import_module(import_path)
        signatures[import_path] = extract_capi(mod)

    module_names = " / ".join(attr for _, attr in modules)

    # --- Write JSON signatures file ---
    repo_root = Path(__file__).resolve().parent.parent
    tests_dir = repo_root / "scipy" / submodule / "tests"
    json_path = tests_dir / "cython_abi_signatures.json"
    test_path = tests_dir / "test_cython_abi.py"

    json_path.write_text(json.dumps(signatures, indent=2) + "\n")
    print(f"Wrote {json_path}")

    # --- Write test file ---
    header = f"""\
# ---------------------------------------------------------------------------
# AUTO-GENERATED - do not edit by hand.
# Re-generate with:  python tools/generate_cython_abi_tests.py {submodule}
#
# Generated against SciPy {scipy.__version__} on {date.today().isoformat()}
# Python {sys.version.split()[0]}
# ---------------------------------------------------------------------------

\"\"\"
Regression test for the {module_names} binary interface.

The __pyx_capi__ dict of each module contains PyCapsule objects whose
*name* string encodes the C-level function signature that Cython checks
at import time when another module does ``cimport``.  Any change to
those strings is a binary incompatibility that breaks every downstream
package (e.g., scikit-learn, statsmodels) compiled against the previous release,
without any compile-time warning.

The expected signatures are stored in cython_abi_signatures.json.
The test fails if:
  * an existing name disappears, or
  * an existing name's signature string changes.
Adding new names is explicitly allowed (it is not a break).
\"\"\"

import ctypes
import importlib
import importlib.resources
import json

import pytest
import scipy

CYTHON_BLAS_ILP64 = (
    scipy.show_config(mode='dicts')['Build Dependencies']['blas']['cython blas ilp64']
)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

_get_capsule_name = ctypes.pythonapi.PyCapsule_GetName
_get_capsule_name.restype = ctypes.c_char_p
_get_capsule_name.argtypes = [ctypes.py_object]


def _extract_capi(module):
    \"\"\"Return {{name: signature}} for every entry in module.__pyx_capi__.\"\"\"
    result = {{}}
    for name, capsule in module.__pyx_capi__.items():
        raw = _get_capsule_name(capsule)
        if raw is not None:
            result[name] = raw.decode('utf-8')
    return result


_sig_file = importlib.resources.files('scipy.{submodule}.tests').joinpath(
    'cython_abi_signatures.json'
)
with importlib.resources.as_file(_sig_file) as _f:
    _EXPECTED = json.loads(_f.read_text())

"""

    # Emit one test function per module
    test_funcs = []
    for import_path, attr in modules:
        # Build optional skip decorator
        skip_decorator = ""
        if skip_condition:
            skip_decorator = (
                f"@pytest.mark.skipif(\n"
                f"    {skip_condition},\n"
                f"    reason={skip_reason!r},\n"
                f")\n"
            )

        # Use a regular string with .format() to avoid nested f-string escaping
        test_funcs.append(textwrap.dedent("""
            {skip_decorator}def test_{attr}_abi_stability():
                \"\"\"No existing {attr} signature may change or disappear.\"\"\"
                expected = _EXPECTED[{import_path!r}]
                mod = importlib.import_module({import_path!r})
                actual = _extract_capi(mod)
                errors = []
                for name, expected_sig in expected.items():
                    if name not in actual:
                        errors.append(
                            f"REMOVED  {{name!r}} (was {{expected_sig!r}})"
                        )
                    elif actual[name] != expected_sig:
                        errors.append(
                            f"CHANGED  {{name!r}}\\n"
                            f"  expected: {{expected_sig!r}}\\n"
                            f"  actual:   {{actual[name]!r}}"
                        )
                if errors:
                    joined = '\\n'.join(errors)
                    pytest.fail(
                        f"{import_path}.__pyx_capi__ has {{len(errors)}} "
                        f"ABI breakage(s):\\n{{joined}}"
                    )
        """).format(
            attr=attr,
            import_path=import_path,
            skip_decorator=skip_decorator,
        ))

    test_path.write_text(header + "\n".join(test_funcs))
    print(f"Wrote {test_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Cython ABI stability tests for a SciPy submodule."
    )
    parser.add_argument(
        "submodule",
        choices=sorted(SUBMODULE_MAP),
        help="SciPy submodule to generate tests for",
    )
    args = parser.parse_args()

    generate(args.submodule, config=SUBMODULE_MAP[args.submodule])
