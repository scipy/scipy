"""Code generation for Boost Special Functions."""

from textwrap import dedent

from _common import _UFuncDef


def _handle_function(f):
    """Generate specialized ufunc loops for Boost Special Functions."""
    ufuncdefs = []
    loop_fun_names = []
    loop_types = []
    loop_funs = []
    ninputs = f["numArgs"]
    for t in f["specializedTypes"]:
        loop_fun_name = f"{t.replace(' ', '_')}_{f['scipyFunctionName']}"
        tmpl_params = ""
        if f["numTemplateParameters"] > 0:
            tmpl_params = f"<{', '.join([t]*f['numTemplateParameters'])}>"
        extra_policy_arg = ""
        if "extraPolicyArg" in f and f["extraPolicyArg"]:
            extra_policy_arg = ", " + f["extraPolicyArg"]
        loop_fun = dedent(f"""\
            static void {loop_fun_name}(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data)
            {{
                npy_intp ii;
                npy_intp n = dimensions[0];
                char {', '.join(f'*in{ii + 1} = args[{ii}]' for ii in range(ninputs))};
                char *out = args[{ninputs}];
                npy_intp {', '.join(f'in{ii + 1}_step = steps[{ii}]' for ii in range(ninputs))}, out_step = steps[{ninputs}];
    
                for (ii = 0; ii < n; ++ii) {{
                    *(({t} *)out) = {f["boostFunctionName"]}{tmpl_params}({', '.join(f'*({t} *)in{ii + 1}' for ii in range(ninputs))}{extra_policy_arg});
    
                    {' '.join(f'in{ii + 1} += in{ii + 1}_step;' for ii in range(ninputs))}
                    out += out_step;
                }}
            }}""")
        loop_fun_names.append(loop_fun_name)
        loop_funs.append(loop_fun)
        loop_types.append(t)
    ufuncdefs.append(_UFuncDef(
        name=f["scipyFunctionName"],
        loop_fun_names=loop_fun_names,
        loop_funs=loop_funs,
        types=loop_types,
        num_inputs=ninputs,
        includes=[f["boostInclude"]]))
    return ufuncdefs
