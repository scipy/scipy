from functools import partial
import inspect
import sys
from importlib import import_module
from importlib.util import find_spec
import os
from pathlib import Path
from types import FunctionType, ModuleType
from typing import Dict, List

from scipy.sparse import linalg

from .array_api_version import __array_api_version__

def not_implemented(*args, **kwargs):
    raise NotImplementedError(
        f"{inspect.stack()[1][3]} is not implemented likely "
        f"because the operation does not make sense for sparse data."
    )

def fill_with_not_implemented(object_to_be_filled, expected_methods):
    is_dict = isinstance(object_to_be_filled, dict)
    for name, implementation in expected_methods:
        if is_dict:
            has_method = name in object_to_be_filled
        else:
            has_method = hasattr(object_to_be_filled, name)
        if not has_method:
            not_implemented_with_signature = partial(not_implemented)
            not_implemented_with_signature.__signature__ = inspect.signature(implementation)
            not_implemented_with_signature.__name__ = name
            if is_dict:
                object_to_be_filled[name] = not_implemented_with_signature
            else:
                setattr(object_to_be_filled, name, not_implemented_with_signature)

__all__ = [
    "name_to_func",
    "array_methods",
    "array_attributes",
    "category_to_funcs",
    "EXTENSIONS",
    "extension_to_funcs",
]

spec_module = "_" + __array_api_version__.replace('.', '_')
try:
    array_api_repo = Path(os.environ["ARRAY_API_REPO_PATH"])
    spec_dir = array_api_repo / "spec" / __array_api_version__ / "API_specification"
    assert spec_dir.exists(), f"{spec_dir} not found - the array api needs to be checked out to the top level"
    sigs_dir = array_api_repo / "src" / "array_api_stubs" / spec_module
    assert sigs_dir.exists()

    sigs_abs_path: str = str(sigs_dir.parent.parent.resolve())
    sys.path.append(sigs_abs_path)
    assert find_spec(f"array_api_stubs.{spec_module}") is not None

    name_to_mod: Dict[str, ModuleType] = {}
    for path in sigs_dir.glob("*.py"):
        name = path.name.replace(".py", "")
        name_to_mod[name] = import_module(f"array_api_stubs.{spec_module}.{name}")

    array = name_to_mod["array_object"].array
    array_methods = [
        f for n, f in inspect.getmembers(array, predicate=inspect.isfunction)
        if n != "__init__"  # probably exists for Sphinx
    ]

    array_methods_with_names_without___init__ = [
        (n, f) for n, f in inspect.getmembers(array, predicate=inspect.isfunction)
        if n != "__init__"  # probably exists for Sphinx
    ]

    def fill_array_with_not_implemented(array):
        return fill_with_not_implemented(array, array_methods_with_names_without___init__)


    top_level_functions = []
    for name, mod in name_to_mod.items():
        if name.endswith("_functions"):
            top_level_functions.extend([
                (n, f) for n, f in inspect.getmembers(mod, predicate=inspect.isfunction)
                if n in mod.__all__
            ])

    EXTENSIONS: List[str] = ["linalg"]  # TODO: add "fft" once stubs available
    extension_to_funcs: Dict[str, List[FunctionType]] = {}
    for ext in EXTENSIONS:
        mod = name_to_mod[ext]
        extension_to_funcs[ext] = [
            (n, f) for n, f in inspect.getmembers(mod, predicate=inspect.isfunction)
            if n in mod.__all__
        ]
    def fill_linalg_with_not_implemented():
        return fill_with_not_implemented(linalg, extension_to_funcs["linalg"])

    def fill_array_module_with_not_implemented(module):
        return fill_with_not_implemented(module, top_level_functions)

    # for funcs in extension_to_funcs.values():
    #     for func in funcs:
    #         if func.__name__ not in name_to_func.keys():
    #             name_to_func[func.__name__] = func

    # # sanity check public attributes are not empty
    # for attr in __all__:
    #     assert len(locals()[attr]) != 0, f"{attr} is empty"
except KeyError:
 
    def fill_array_with_not_implemented(array):
        pass
    def fill_linalg_with_not_implemented():
        pass
    def fill_array_module_with_not_implemented(module):
        pass