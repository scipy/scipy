"""Generate a single module with ufuncs for Boost functions."""

from typing import List
import json
from textwrap import dedent
import itertools
import pathlib
import argparse

from _common import _UFuncDef
from _distribution_codegen import _handle_distribution, _add_distribution_defaults
from _function_codegen import _handle_function

THIS_DIR = pathlib.Path(__file__).parent


def _make_boost_module(ufuncdefs: List[_UFuncDef]):
    all_includes = set(itertools.chain(*[d.includes for d in ufuncdefs]))
    NL = "\n        "
    src = dedent(f"""\
        #define PY_SSIZE_T_CLEAN
        #include <Python.h>
        #include "numpy/ndarraytypes.h"
        #include "numpy/ufuncobject.h"
        #include "numpy/halffloat.h"
        #include <math.h>
        
        {NL.join(f'#include "{inc}"' for inc in all_includes)}
        
        static PyMethodDef BoostMethods[] = {{
            {{NULL, NULL, 0, NULL}}
        }};
    
    """)

    for d in ufuncdefs:
        for loop_fun in d.loop_funs:
            src += loop_fun + "\n\n"

        src += f"PyUFuncGenericFunction {d.name}_funcs[{len(d.loop_funs)}] = {{"
        for name in d.loop_fun_names:
            src += f"\n    &{name},"
        src += "\n};\n\n"

        src += f"static char {d.name}_types[] = {{"
        for t in d.types:
            NPY_TYPE = {
                "npy_half": "NPY_HALF",
                "float": "NPY_FLOAT",
                "double": "NPY_DOUBLE",
                "long double": "NPY_LONGDOUBLE",
            }[t]
            src += f"\n    " + f"{NPY_TYPE}, "*(1 + d.num_inputs)
        src += "\n};\n\n"

    src += dedent("""\
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_boostgen",
        NULL,
        -1,
        BoostMethods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    
    static void *data[1] = {NULL};  // see https://github.com/numpy/numpy/issues/14613
    
    PyMODINIT_FUNC PyInit__boostgen(void)
    {
        PyObject *m, *d;
        m = PyModule_Create(&moduledef);
        if (!m) {
            return NULL;
        }
    
        import_array();
        import_umath();
        d = PyModule_GetDict(m);

    """)

    pyobj_names = []
    for d in ufuncdefs:
        pyobj_name = f"{d.name}_ufunc"
        pyobj_names.append(pyobj_name)
        src += f'    PyObject *{pyobj_name};\n'

    src += "\n"
    for pyobj_name, d in zip(pyobj_names, ufuncdefs):
        src += f'    {pyobj_name} = PyUFunc_FromFuncAndData({d.name}_funcs, data, {d.name}_types, {len(d.types)}, {d.num_inputs}, 1, PyUFunc_None, "{d.name}", "docstring placeholder", 0);\n'
        src += f'    PyDict_SetItemString(d, "{d.name}", {pyobj_name});\n'

    src += "\n"
    for pyobj_name in pyobj_names:
        src += f"    Py_DECREF({pyobj_name});\n"

    src += dedent("""\
        return m;
    }
    """)

    return src


def main(outdir: str):
    # Boost "things" we want to consider:
    # - distributions
    # - special functions
    ufuncdefs = []

    # start with statistical distributions:
    with open(THIS_DIR / "shortcut_distributions.json", "r") as fp:
        shortcut_distribution_list = json.load(fp)
    with open(THIS_DIR / "distributions.json", "r") as fp:
        distribution_list = json.load(fp)
    # fill-in defaults
    _add_distribution_defaults(shortcut_distribution_list, distribution_list)
    for f in distribution_list:
        ufuncdefs += _handle_distribution(f)

    # collect and generate special functions now
    with open(THIS_DIR / "functions.json", "r") as fp:
        function_list = json.load(fp)
    for f in function_list:
        ufuncdefs += _handle_function(f)

    # now create ufuncs in a single boost module
    src = _make_boost_module(ufuncdefs=ufuncdefs)
    outfile = (pathlib.Path(outdir) / "boost_ufuncs.cpp").resolve().absolute()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with open(str(outfile), "w") as fp:
        fp.write(src)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--outdir", type=str, help="Location to write generated files.")
    args = parser.parse_args()
    main(outdir=args.outdir)
