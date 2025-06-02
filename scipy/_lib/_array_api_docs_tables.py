"""Generate flat tables showing Array API capabilities for use in docs.

These tables are intended for presenting Array API capabilities across
a wide number of functions at once. Rows correspond to functions and
columns correspond to library/device/option combinations.
"""

import types

from collections import defaultdict
from importlib import import_module

from scipy._lib._array_api import xp_capabilities_table
from scipy._lib._array_api import _make_sphinx_capabilities

# For undocumented aliases of public functions which are kept around for
# backwards compatibility reasons. These should be excluded from the
# tables since they would be redundant. There are also no docs pages to
# link entries to.
ALIASES = {
    "scipy.linalg": {
        # Alias of scipy.linalg.solve_continuous_lyapunov
        "solve_lyapunov",
    },
    "scipy.ndimage": {
        # Alias of scipy.ndimage.sum_labels
        "sum",
    },
    "scipy.special": {
        # Alias of scipy.special.jv
        "jn",
        # Alias of scipy.special.roots_legendre
        "p_roots",
        # Alias of scipy.special.roots_chebyt
        "t_roots",
        # Alias of scipy.special.roots_chebyu
        "u_roots",
        # Alias of scipy.special.roots_chebyc
        "c_roots",
        # Alias of scipy.special.roots_chebys
        "s_roots",
        # Alias of scipy.special.roots_jacobi
        "j_roots",
        # Alias of scipy.special.roots_laguerre
        "l_roots",
        # Alias of scipy.special.roots_genlaguerre
        "la_roots",
        # Alias of scipy.special.roots_hermite
        "h_roots",
        # Alias of scipy.special.roots_hermitenorm
        "he_roots",
        # Alias of scipy.special.roots_gegenbauer
        "cg_roots",
        # Alias of scipy.special.roots_sh_legendre
        "ps_roots",
        # Alias of scipy.special.roots_sh_chebyt
        "ts_roots",
        # Alias of scipy.special.roots_chebyu
        "us_roots",
        # Alias of scipy.special.roots_sh_jacobi
        "js_roots",
    }
}

# Shortened names for use in table.
BACKEND_NAMES_MAP = {
    "jax.numpy": "jax",
    "dask.array": "dask",
}


def _process_capabilities_table_entry(entry: dict | None) -> dict[str, dict[str, bool]]:
    """Returns dict showing alternative backend support in easy to consume form.

    Parameters
    ----------
    entry : Optional[dict]
       A dict with the structure of the values of the dict
       scipy._lib._array_api.xp_capabilities_table. If None, it is
       assumped that no alternative backends are supported.
       Default: None.

    Returns
    -------
    dict[str, dict[str, bool]]
        The output dict currently has keys "cpu", "gpu", "jit" and "lazy".
        The value associated to each key is itself a dict. The keys of
        the inner dicts correspond to backends, with bool values stating
        whether or not the backend is supported with a given device or
        mode. Inapplicable backends do not appear in the inner dicts
        (e.g. since cupy is gpu-only, it does not appear in the inner
        dict keyed on "cpu"). Only alternative backends to NumPy are
        included since NumPY support should be guaranteed.

    """
    output = {
        "cpu": {"torch": False, "jax": False, "dask": False},
        "gpu": {"cupy": False, "torch": False, "jax": False},
        "jit": {"jax": False},
        "lazy": {"dask": False},
    }
    if entry is None:
        # If there is no entry, assume no alternative backends are supported.
        # If the list of supported backends will grows, this hard-coded dict
        # will need to be updated.
        return output

    # For now, use _make_sphinx_capabilities because that's where
    # the relevant logic for determining what is and isn't
    # supported based on xp_capabilities_table entries lives.
    # Perhaps this logic should be decoupled from sphinx.
    for backend, capabilities in _make_sphinx_capabilities(**entry).items():
        if backend in {"array_api_strict", "numpy"}:
            continue
        backend = BACKEND_NAMES_MAP.get(backend, backend)
        cpu, gpu = capabilities.cpu, capabilities.gpu
        if cpu is not None:
            if backend not in output["cpu"]:
                raise ValueError(
                    "Input capabilities table entry contains unhandled"
                    f" backend {backend} on cpu."
                )
            output["cpu"][backend] = capabilities.cpu
        if gpu is not None:
            if backend not in output["gpu"]:
                raise ValueError(
                    "Input capabilities table entry contains unhandled"
                    f" backend {backend} on gpu."
                )
            output["gpu"][backend] = capabilities.gpu
        if backend == "jax":
            output["jit"]["jax"] = entry["jax_jit"] and output["cpu"]["jax"]
        if backend == "dask.array":
            support_lazy = not entry["allow_dask_compute"] and output["dask"]
            output["lazy"]["dask"] = support_lazy
    return output


def make_flat_capabilities_table(
        modules: str | list[str],
        backend_type: str,
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

    backend_type : {'cpu', 'gpu', 'jit', 'lazy'}

    capabilities_table : Optional[list[str]]
        Table in the form of `scipy._lib._array_api.xp_capabilities_table`.
        If None, uses `scipy._lib._array_api.xp_capabilities_table`.
        Default: None.

    Returns
    -------
    output : list[dict[str, str]]
        `output` is a table in dict format
        (keys corresponding to column names). If `modules` is a list, then the
        first column will be "module" If `modules` is a string, then the
        "module" column will be omitted. The other columns correspond to
        supported backends for the given `backend_type`,
        e.g. jax, torch, and dask on cpu. numpy is excluded because it should
         always be supported. See the helper function
        `_process_capabilities_table_entry` above).

    """
    if backend_type not in {"cpu", "gpu", "jit", "lazy"}:
        raise ValueError(f"Received unhandled backend type {backend_type}")

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
            if name in ALIASES.get(module_name, {}):
                # Skip undocumented aliases that are kept
                # for backwards compatibility reasons.
                continue
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
            if callable(thing) and hasattr(thing, "__name__"):
                entry = xp_capabilities_table.get(thing, None)
                capabilities = _process_capabilities_table_entry(entry)[backend_type]
                # If a list of multiple modules is passed in, add the module
                # as an entry of the table.
                row = {"module": module_name} if multiple_modules else {}
                row.update({"function": name})
                row.update(capabilities)
                output.append(row)
            else:
                # Skip anything else which isn't a callable. Also skip unnamed
                # callables. Currently the only unnamed callables are distributions
                # from the old distribution infrastructure.
                continue
    return output


def calculate_table_statistics(
    flat_table: list[dict[str, str]]
) -> dict[str, str] | dict[str, dict[str, str]]:
    """Get counts of what is supported per module.

    Parameters
    ----------
    flat_table : list[dict[str, str]]
        A table as returned by `make_flat_capabilities_table`

    Returns
    -------
    dict[str, str] | dict[str, dict[str, str]]
        If `flat_table` is only for one module (and has no module column)
        then the output is a dictionary with a key "total" and a key for
        each backend column of the flat capabilities table. The value
        corresponding to total is the total count of functions in the given
        module, and the value associated to the other keys is the count of
        functions that support that particular backend/details.

        If `flat_table` has a module column, then the output will be a
        dictionary with keys for each included module and values of type
        dict[str, str] of the kind returned by this function when `flat_table`
        has no module column.
    """
    if not flat_table:
        return []
    table_contains_modules = "module" in flat_table[0]
    if table_contains_modules:
        counter = defaultdict(lambda: defaultdict(int))
    else:
        counter = defaultdict(int)

    for entry in flat_table:
        entry = entry.copy()
        entry.pop("function")
        if table_contains_modules:
            module = entry.pop("module")
            current_counter = counter[module]
        else:
            current_counter = counter
        current_counter["total"] += 1
        for key, value in entry.items():
            current_counter[key] += value
    # Strip away dangerous defaultdictness
    if table_contains_modules:
        counter = {key: dict(value) for key, value in counter.items()}
    else:
        counter = dict(counter)
    return counter
