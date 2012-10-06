#!/usr/bin/python
"""
generate_ufuncs.py

Generate Ufunc definition source files for scipy.special.

"""

#---------------------------------------------------------------------------------
# Ufunc listing
#---------------------------------------------------------------------------------

# Ufuncs without C++
UFUNCS = """
_lambertw -- lambertw_scalar: Dld->D                       -- lambertw.pxd
logit -- logitf: f->f, logit: d->d, logitl: g->g           -- _logit.h
expit -- expitf: f->f, expit: d->d, expitl: g->g           -- _logit.h
"""

# Ufuncs with C++
UFUNCS_CXX = """
"""

#---------------------------------------------------------------------------------
# Code generation
#---------------------------------------------------------------------------------

import subprocess
import re
import textwrap
import add_newdocs

C_TYPES = {
    'f': 'float',
    'd': 'double',
    'g': 'long double',
    'F': 'float complex',
    'D': 'double complex',
    'G': 'long double complex',
    'i': 'int',
    'l': 'long',
    'v': 'void',
}

TYPE_NAMES = {
    'f': 'np.NPY_FLOAT',
    'd': 'np.NPY_DOUBLE',
    'g': 'np.NPY_LONGDOUBLE',
    'F': 'np.NPY_CFLOAT',
    'D': 'np.NPY_CDOUBLE',
    'G': 'np.NPY_CLONGDOUBLE',
    'i': 'np.NPY_INT',
    'l': 'np.NPY_LONG',
}

def generate_loop(func_inputs, func_outputs, func_retval,
                  ufunc_inputs=None, ufunc_outputs=None):
    """
    Generate a UFunc loop function that calls a function given as its
    data parameter with the specified input and output arguments and
    return value.

    This function can be passed to PyUFunc_FromFuncAndData.

    Parameters
    ----------
    func_inputs, func_outputs, func_retval : str
        Signature of the function to call, given as type codes of the
        input, output and return value arguments. These 1-character
        codes are given according to the C_TYPES and TYPE_NAMES lists
        above.

        The corresponding C function signature to be called is:

            retval func(intype1 iv1, intype2 iv2, ..., outtype1 *ov1, ...);

        If the function does not have any output arguments, its return
        value is considered to be the output. Otherwise, the return
        value is ignored.

    ufunc_inputs, ufunc_outputs : str
        Ufunc input and output signature.

        This does not have to exactly match the function signature,
        as long as the type casts work out on the C level.

    Returns
    -------
    loop_name
        Name of the generated loop function.
    loop_body
        Generated C code for the loop.

    """
    if ufunc_inputs is None:
        ufunc_inputs = func_inputs
        if not func_outputs:
            ufunc_outputs = func_retval
        else:
            ufunc_outputs = func_outputs

    if len(func_inputs) != len(ufunc_inputs):
        raise ValueError("Function and ufunc have different number of inputs")
    if not func_outputs:
        if len(func_retval) != 1 or len(ufunc_outputs) != 1:
            raise ValueError("Function retval and ufunc outputs don't match")
    else:
        if len(func_outputs) != len(ufunc_outputs):
            raise ValueError("Function and ufunc have different number of outputs")

    name = "loop_%s_%s_%s_As_%s_%s" % (
        func_retval, func_inputs, func_outputs, ufunc_inputs, ufunc_outputs
        )
    body = "cdef void %s(char **args, np.npy_intp *dims, np.npy_intp *steps, void *func) nogil:\n" % name
    body += "    cdef np.npy_intp i, n = dims[0]\n"

    pointers = []
    for j in range(len(ufunc_inputs)):
        pointers.append("*ip%d = args[%d]" % (j, j))
    for j in range(len(ufunc_outputs)):
        pointers.append("*op%d = args[%d]" % (j, j + len(ufunc_inputs)))
    body += "    cdef char %s\n" % ", ".join(pointers)

    ftypes = []
    fvars = []
    outtypecodes = []
    for j in range(len(func_inputs)):
        ftypes.append(C_TYPES[func_inputs[j]])
        fvars.append("(<%s*>ip%d)[0]" % (C_TYPES[ufunc_inputs[j]], j))
    for j, outtype in enumerate(func_outputs):
        body += "    cdef %s ov%d\n" % (C_TYPES[outtype], j)
        ftypes.append("%s *" % C_TYPES[outtype])
        fvars.append("&ov%d" % j)
        outtypecodes.append(outtype)

    if not func_outputs:
        outtypecodes.append(func_retval)
        body += "    cdef %s ov0\n" % (C_TYPES[func_retval],)

    body += "    for i in range(n):\n"
    if not func_outputs:
        body += "        ov0 = (<%s(*)(%s) nogil>func)(%s)\n" % (
            C_TYPES[func_retval], ", ".join(ftypes), ", ".join(fvars))
    else:
        body += "        (<%s(*)(%s) nogil>func)(%s)\n" % (
            C_TYPES[func_retval], ", ".join(ftypes), ", ".join(fvars))

    for j, (outtype, fouttype) in enumerate(zip(ufunc_outputs, outtypecodes)):
        if outtype in 'FDG' and outtype != fouttype:
            if fouttype in 'FDG':
                body += "        (<%s *>op%d)[0].real = <%s>ov%d.real\n" % (
                    C_TYPES[outtype], j, C_TYPES[outtype.lower()], j)
                body += "        (<%s *>op%d)[0].imag = <%s>ov%d.imag\n" % (
                    C_TYPES[outtype], j, C_TYPES[outtype.lower()], j)
            else:
                body += "        (<%s *>op%d)[0].real = <%s>ov%d\n" % (
                    C_TYPES[outtype], j, j)
                body += "        (<%s *>op%d)[0].imag = 0\n" % (
                    C_TYPES[outtype], j, j)
        else:
            body += "        (<%s *>op%d)[0] = <%s>ov%d\n" % (
                C_TYPES[outtype], j, C_TYPES[outtype], j)
    for j in range(len(ufunc_inputs)):
        body += "        ip%d += steps[%d]\n" % (j, j)
    for j in range(len(ufunc_outputs)):
        body += "        op%d += steps[%d]\n" % (j, j + len(ufunc_inputs))

    return name, body

def iter_variants(inputs, outputs):
    """
    Generate variants of UFunc signatures, by changing variable types,
    within the limitation that the corresponding C types casts still
    work out.

    This does not generate all possibilities, just the ones required
    for the ufunc to work properly with the most common data types.

    Parameters
    ----------
    inputs, outputs : str
        UFunc input and output signature strings

    Yields
    ------
    new_input, new_output : str
        Modified input and output strings.
        Also the original input/output pair is yielded.

    """
    yield inputs, outputs
    yield inputs.replace('d', 'f'), outputs.replace('d', 'f')
    yield inputs.replace('D', 'F'), outputs.replace('D', 'F')
    yield inputs.replace('i', 'l'), outputs.replace('i', 'l')

class Ufunc(object):
    """
    Ufunc signature, restricted format suitable for special functions.

    Parameters
    ----------
    name
        Name of the ufunc to create
    signature
        String of form 'func: fff->f, func2: ddd->d' describing
        the C-level functions and types of their input arguments
        and return values.
    
    """
    def __init__(self, name, signatures):
        self.name = name
        self.signatures = self._parse_signatures(signatures)
        self.doc = add_newdocs.get("scipy.special." + name)
        if self.doc is None:
            raise ValueError("No docstring for ufunc %r" % name)
        self.doc = textwrap.dedent(self.doc).strip()

    def _parse_signatures(self, sigs):
        return [self._parse_signature(x) for x in sigs.split(",")
                if x.strip()]

    def _parse_signature(self, sig):
        m = re.match("\s*(.*):\s*([fdgFDGil]*)\\*([fdgFDGil]*)->([fdgFDGil]*)\s*$", sig)
        if m:
            func, inarg, outarg, ret = map(lambda x: x.strip(), m.groups())
            return (func, inarg, outarg, ret)
            raise ValueError("Invalid signature: %r" % sig)
        m = re.match("\s*(.*):\s*([fdgFDGil]*)->([fdgFDGil]?)\s*$", sig)
        if m:
            func, inarg, ret = map(lambda x: x.strip(), m.groups())
            return (func, inarg, "", ret)
        raise ValueError("Invalid signature: %r" % sig)

    def _get_signatures_and_loops(self, all_loops):
        inarg_num = None
        outarg_num = None
        
        variants = {}
        for func_name, inarg, outarg, ret in self.signatures:
            base_input = inarg
            base_output = outarg or ret

            if inarg_num is None:
                inarg_num = len(base_input)
                outarg_num = len(base_output)

            for inp, outp in iter_variants(base_input, base_output):
                sig = (func_name, inp, outp)
                if (inp, outp) in variants:
                    continue
                if "v" in outp:
                    raise ValueError("%s: void signature %r" % (self.name, sig))
                if len(inp) != inarg_num or len(outp) != outarg_num:
                    raise ValueError("%s: signature %r does not have %d/%d input/output args" % (
                        self.name, sig,
                        inarg_num, outarg_num))

                loop_name, loop = generate_loop(inarg, outarg, ret,
                                                inp, outp)
                all_loops[loop_name] = loop
                variants[(inp, outp)] = (func_name, loop_name, inp, outp)
        return variants, inarg_num, outarg_num

    def generate(self, all_loops):
        toplevel = ""

        variants, inarg_num, outarg_num = self._get_signatures_and_loops(all_loops)
        variants = variants.items()
        variants.sort()

        loops = []
        datas = []
        types = []

        for (inp, outp), (func_name, loop_name, inputs, outputs) in variants:
            for x in inputs:
                types.append(TYPE_NAMES[x])
            for x in outputs:
                types.append(TYPE_NAMES[x])
            loops.append(loop_name)
            datas.append(func_name)

        toplevel += "cdef np.PyUFuncGenericFunction ufunc_%s_loops[%d]\n" % (self.name, len(loops))
        toplevel += "cdef void *ufunc_%s_data[%d]\n" % (self.name, len(datas))
        toplevel += "cdef char ufunc_%s_types[%d]\n" % (self.name, len(types))
        toplevel += 'cdef char *ufunc_%s_doc = (\n    "%s")\n' % (
            self.name,
            self.doc.replace('"', '\\"').replace('\n', '\\n\"\n    "')
            )

        for j, function in enumerate(loops):
            toplevel += "ufunc_%s_loops[%d] = <np.PyUFuncGenericFunction>%s\n" % (self.name, j, function)
        for j, type in enumerate(types):
            toplevel += "ufunc_%s_types[%d] = <char>%s\n" % (self.name, j, type)
        for j, data in enumerate(datas):
            toplevel += "ufunc_%s_data[%d] = <void*>_func_%s\n" % (self.name, j, data)

        toplevel += ('@ = np.PyUFunc_FromFuncAndData(ufunc_@_loops, '
                     'ufunc_@_data, ufunc_@_types, %d, %d, %d, 0, '
                     '"@", ufunc_@_doc, 0)\n' % (len(types)/(inarg_num+outarg_num),
                                                 inarg_num, outarg_num)
                     ).replace('@', self.name)

        return toplevel, list(set(datas))

def generate(filename, ufunc_str):
    ufuncs = []
    headers = {}

    for line in ufunc_str.splitlines():
        line = line.strip()
        if not line:
            continue
        m = re.match("^([a-z0-9_]+)\s*--\s*(.*?)\s*(--.*)?$", line)
        if not m:
            raise ValueError("Unparseable line %r" % line)
        ufuncs.append(Ufunc(m.group(1), m.group(2)))
        if m.group(3):
            headers[ufuncs[-1].name] = m.group(3)[2:].strip()

    toplevel = ""
    defs = ""
    all_loops = {}

    ufuncs.sort(key=lambda u: u.name)
    for ufunc in ufuncs:
        t, cfuncs = ufunc.generate(all_loops)
        toplevel += t + "\n"

        header = headers.get(ufunc.name, 'cephes.h')
        if cfuncs:
            if header.endswith('.pxd'):
                for cfunc in cfuncs:
                    defs += "from %s cimport %s as _func_%s\n" % (header[:-4], cfunc, cfunc)
            else:
                defs += "cdef extern from \"%s\":\n" % header
                for cfunc in cfuncs:
                    defs += "    void _func_%s \"%s\"()\n" % (cfunc, cfunc)

    toplevel = "\n".join(all_loops.values() + [defs, toplevel])

    f = open(filename, 'wb')
    f.write("""\
# This file is automatically generated by generate_ufuncs.py.
# Do not edit manually!

cdef extern from "complex.h":
    pass

cimport numpy as np
cimport libc

np.import_array()
np.import_ufunc()

""")

    f.write(toplevel)

    f.write("\n\n")
    for ufunc in ufuncs:
        f.write("# %s = _ufuncs.%s\n" % (ufunc.name, ufunc.name))

    f.close()

def main():
    generate("_ufuncs.pyx", UFUNCS)
    generate("_ufuncs_cxx.pyx", UFUNCS_CXX)

    subprocess.call(['cython', '_ufuncs.pyx'])

if __name__ == "__main__":
    main()
