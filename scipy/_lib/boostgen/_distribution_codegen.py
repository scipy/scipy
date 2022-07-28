"""Code generation for Boost Statistical Distributions."""

from textwrap import dedent

from _common import _UFuncDef


def _handle_distribution(f):
    """Generate specialized ufunc loops for Boost Statistical Distributions."""
    ufuncdefs = []
    for method in f["methods"]:
        loop_fun_names = []
        loop_types = []
        loop_funs = []
        ninputs = f["numCtorArgs"] + method["numMethodArgs"]
        for t in f["specializedTypes"]:
            loop_fun_name = f"{t.replace(' ', '_')}_{method['scipyFunctionName']}"
            boost_call = method["boostMethodName"]
            if "complement" in method and method["complement"]:
                boost_call += "(boost::math::complement"
            else:
                method["complement"] = False
            loop_fun = dedent(f"""\
                static void {loop_fun_name}(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
                {{
                    npy_intp ii;
                    npy_intp n = dimensions[0];
                    char {', '.join(f'*in{ii + 1} = args[{ii}]' for ii in range(ninputs))};
                    char *out = args[{ninputs}];
                    npy_intp {', '.join(f'in{ii + 1}_step = steps[{ii}]' for ii in range(ninputs))}, out_step = steps[{ninputs}];

                    for (ii = 0; ii < n; ++ii) {{
                        {f["boostClassName"]}<{t}> obj({', '.join(f'*({t} *)in{ii + 1}' for ii in range(f["numCtorArgs"]))});
                        *(({t} *)out) = {boost_call}(obj{', ' * method['numMethodArgs']}{', '.join(f'*({t} *)in{ii + 1}' for ii in range(f["numCtorArgs"], ninputs))}){')' * method['complement']};

                        {' '.join(f'in{ii + 1} += in{ii + 1}_step;' for ii in range(ninputs))}
                        out += out_step;
                    }}
                }}""")
            loop_fun_names.append(loop_fun_name)
            loop_funs.append(loop_fun)
            loop_types.append(t)
        ufuncdefs.append(_UFuncDef(
            name=method["scipyFunctionName"],
            loop_fun_names=loop_fun_names,
            loop_funs=loop_funs,
            types=loop_types,
            num_inputs=ninputs,
            includes=[f["boostInclude"]]))
    return ufuncdefs


def _add_distribution_defaults(shortcut_distribution_list, distribution_list):
    """Generate complete distribution list entries using defaults given basic distribution information."""
    for d in shortcut_distribution_list:
        distribution_list.append({
            "boostInclude": d["boostInclude"],
            "boostClassName": d["boostClassName"],
            "specializedTypes": ["float", "double", "long double"],
            "numCtorArgs": d["numCtorArgs"],
            "methods": [
                {
                    "boostMethodName": "boost::math::pdf",
                    "scipyFunctionName": f"{d['scipyPrefix']}_pdf",
                    "numMethodArgs": 1
                },
                {
                    "boostMethodName": "boost::math::cdf",
                    "scipyFunctionName": f"{d['scipyPrefix']}_cdf",
                    "numMethodArgs": 1
                },
                {
                    "boostMethodName": "boost::math::cdf",
                    "complement": True,
                    "scipyFunctionName": f"{d['scipyPrefix']}_sf",
                    "numMethodArgs": 1
                },
                {
                    "boostMethodName": "boost::math::quantile",
                    "scipyFunctionName": f"{d['scipyPrefix']}_ppf",
                    "numMethodArgs": 1
                },
                {
                    "boostMethodName": "boost::math::quantile",
                    "complement": True,
                    "scipyFunctionName": f"{d['scipyPrefix']}_isf",
                    "numMethodArgs": 1
                },
                {
                    "boostMethodName": "boost::math::mean",
                    "scipyFunctionName": f"{d['scipyPrefix']}_mean",
                    "numMethodArgs": 0
                },
                {
                    "boostMethodName": "boost::math::variance",
                    "scipyFunctionName": f"{d['scipyPrefix']}_variance",
                    "numMethodArgs": 0
                },
                {
                    "boostMethodName": "boost::math::skewness",
                    "scipyFunctionName": f"{d['scipyPrefix']}_skewness",
                    "numMethodArgs": 0
                },
                {
                    "boostMethodName": "boost::math::kurtosis_excess",
                    "scipyFunctionName": f"{d['scipyPrefix']}_kurtosis_excess",
                    "numMethodArgs": 0
                }
            ]
        })
