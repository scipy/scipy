"""Generate flat tables showing Array API capabilities for use in docs.

These tables are intended for presenting Array API capabilities across
a wide number of functions at once. Rows correspond to functions and
columns correspond to library/device/option combinations.
"""

import inspect
import types

from collections import defaultdict
from importlib import import_module

from scipy._lib._array_api import xp_capabilities_table
from scipy._lib._array_api import _make_sphinx_capabilities


def _process_capabilities_table_entry(entry: dict | None) -> dict[str, str]:
    """Returns flat dict showing what is and isn't supported in entry."""
    if entry is None:
        # If there is no entry, assume no alternative backends are supported.
        # If the list of supported backends will grows, this hard-coded dict
        # will need to be updated.
        return {
            "cupy": False, "torch (cpu)": False,
            "torch (gpu)": False, "jax.numpy (cpu)": False,
            "jax.numpy (gpu)": False, "jax.numpy (jit)": False,
            "dask.array": False, "dask.array (lazy)": False
        }
    row = {}
    # For now, use _make_sphinx_capabilities because that's where
    # the relevant logic for determining what is and isn't
    # supported based on xp_capabilities_table entries lives.
    # Perhaps this logic should be decoupled from sphinx.
    backends = _make_sphinx_capabilities(**entry).keys()
    for backend, capabilities in _make_sphinx_capabilities(**entry).items():
        if backend in {"array_api_strict", "numpy"}:
            continue
        cpu, gpu = capabilities.cpu, capabilities.gpu
        if cpu is None:
            row[backend] = capabilities.gpu
        elif gpu is None:
            row[backend] = capabilities.cpu
        else:
            row[f"{backend} (cpu)"] = capabilities.cpu
            row[f"{backend} (gpu)"] = capabilities.gpu
        if backend == "jax.numpy":
            row["jax.numpy (jit)"] = entry["jax_jit"] and row["jax.numpy (cpu)"]
        if backend == "dask.array":
            support_lazy = not entry["allow_dask_compute"] and row["dask.array"]
            row["dask.array (lazy)"] = support_lazy
    return row


def make_flat_capabilities_table(
        modules: str | list[str],
        /,
        *,
        capabilities_table: list[str] | None = None,
) -> list[dict[str, str]]:
    """Generate full table of array api capabilities across public functions.

    Parameters
    ----------
    modules : str | list[str]
        A string containing single SciPy module, (e.g `scipy.stats`, `scipy.fft`)
        or a list of such strings.

    Returns
    -------
    output : list[dict[str, str]]
        `output` is a table in dict format
        (keys corresponding to column names). If `modules` is a list, then the
        first column will be "module" If `modules` is a string, then the
        "module" column will be omitted. The other colums are currently
    
        "cupy", "torch (cpu)", "torch (gpu)", "jax.numpy (cpu)", "jax.numpy (gpu)",
        "jax.numpy (jit)", "dask.array", "dask.array (lazy)".

        "numpy" is omitted because it will always be supported.
    """
    multiple_modules = True
    if isinstance(modules, str):
        modules = [modules]
        multiple_modules = False
    if capabilities_table is None:
        capabilities_table = xp_capabilities_table

    output = []

    for module_name in modules:
        module = import_module(module_name)
        public_things = module.__all__
        for name in public_things:
            thing = getattr(module, name)
            if isinstance(thing, type) and issubclass(thing, Exception):
                # Skip exception types
                continue
            if isinstance(thing, types.ModuleType):
                # Skip modules
                continue
            if isinstance(thing, type):
                # Skip classes for now, but we may want to handle these in some
                # way later, so giving them their own branch.
                continue
            elif callable(thing) and hasattr(thing, "__name__"):
                entry = xp_capabilities_table.get(thing, None)
                capabilities = _process_capabilities_table_entry(entry)
                # If a list of multiple modules is passed in, add the module
                # as an entry of the table.
                row = {"module": module_name} if multiple_modules else {}
                row.update({"function": thing.__name__})
                row.update(capabilities)
                output.append(row)
            else:
                # Skip anything else which isn't callable and also unnamed
                # callables. The only unnamed callables are distributions
                # from stats.distributions.
                continue
    return output
