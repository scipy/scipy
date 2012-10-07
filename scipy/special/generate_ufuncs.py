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
bdtrc -- bdtrc: iid->d                                     -- cephes.h
bdtr -- bdtr: iid->d                                       -- cephes.h
bdtri -- bdtri: iid->d                                     -- cephes.h
btdtr -- btdtr: ddd->d                                     -- cephes.h
btdtri -- incbi: ddd->d                                    -- cephes.h
fdtrc -- fdtrc: ddd->d                                     -- cephes.h
fdtr -- fdtr: ddd->d                                       -- cephes.h
fdtri -- fdtri: ddd->d                                     -- cephes.h
gdtrc -- gdtrc: ddd->d                                     -- cephes.h
gdtr -- gdtr: ddd->d                                       -- cephes.h
hyp2f1 -- hyp2f1: dddd->d, chyp2f1_wrap: dddD->D           -- cephes.h, specfun_wrappers.h
hyp1f1 -- chyp1f1_wrap: ddD->D, hyp1f1_wrap: ddd->d        -- specfun_wrappers.h
hyperu -- hypU_wrap: ddd->d                                -- specfun_wrappers.h
hyp2f0 -- hyp2f0: dddi*dd->i                               -- cephes.h
hyp1f2 -- onef2: dddd*dd->i                                -- cephes.h
hyp3f0 -- threef0: dddd*dd->i                              -- cephes.h
betainc -- incbet: ddd->d                                  -- cephes.h
betaincinv -- incbi: ddd->d                                -- cephes.h
nbdtrc -- nbdtrc: iid->d                                   -- cephes.h
nbdtr -- nbdtr: iid->d                                     -- cephes.h
nbdtri -- nbdtri: iid->d                                   -- cephes.h
beta -- beta: dd->d                                        -- cephes.h
betaln -- lbeta: dd->d                                     -- cephes.h
cbrt -- cbrt: d->d                                         -- cephes.h
chdtrc -- chdtrc: dd->d                                    -- cephes.h
chdtr -- chdtr: dd->d                                      -- cephes.h
chdtri -- chdtri: dd->d                                    -- cephes.h
dawsn -- dawsn: d->d                                       -- cephes.h
ellipeinc -- ellie: dd->d                                  -- cephes.h
ellipkinc -- ellik: dd->d                                  -- cephes.h
ellipe -- ellpe: d->d                                      -- cephes.h
ellipkm1 -- ellpk: d->d                                    -- cephes.h
exp10 -- exp10: d->d                                       -- cephes.h
exp2 -- exp2: d->d                                         -- cephes.h
gamma -- Gamma: d->d, cgamma_wrap: D->D                    -- cephes.h, specfun_wrappers.h
gammaln -- lgam: d->d, clngamma_wrap: D->D                 -- cephes.h, specfun_wrappers.h
i0 -- i0: d->d                                             -- cephes.h
i0e -- i0e: d->d                                           -- cephes.h
i1 -- i1: d->d                                             -- cephes.h
i1e -- i1e: d->d                                           -- cephes.h
gammaincc -- igamc: dd->d                                  -- cephes.h
gammainc -- igam: dd->d                                    -- cephes.h
gammainccinv -- igami: dd->d                               -- cephes.h
iv -- iv: dd->d, cbesi_wrap: dD->D                         -- cephes.h, amos_wrappers.h
ive -- cbesi_wrap_e_real: dd->d, cbesi_wrap_e: dD->D       -- amos_wrappers.h
ellipj -- ellpj: dd*dddd->i                                -- cephes.h
expn -- expn: id->d                                        -- cephes.h
exp1 -- cexp1_wrap: D->D, exp1_wrap: d->d                  -- specfun_wrappers.h
expi -- cexpi_wrap: D->D, expi_wrap: d->d                  -- specfun_wrappers.h
kn -- kn: id->d                                            -- cephes.h
pdtrc -- pdtrc: id->d                                      -- cephes.h
pdtr -- pdtr: id->d                                        -- cephes.h
pdtri -- pdtri: id->d                                      -- cephes.h
yn -- yn: id->d                                            -- cephes.h
smirnov -- smirnov: id->d                                  -- cephes.h
smirnovi -- smirnovi: id->d                                -- cephes.h
airy -- airy: d*dddd->i, cairy_wrap: D*DDDD->i             -- cephes.h, amos_wrappers.h
itairy -- itairy_wrap: d*dddd->i                           -- specfun_wrappers.h
airye -- cairy_wrap_e: D*DDDD->i, cairy_wrap_e_real: d*dddd->i -- amos_wrappers.h
fresnel -- fresnl: d*dd->i, cfresnl_wrap: D*DD->i          -- cephes.h, specfun_wrappers.h
shichi -- shichi: d*dd->i                                  -- cephes.h
sici -- sici: d*dd->i                                      -- cephes.h
itj0y0 -- it1j0y0_wrap: d*dd->i                            -- specfun_wrappers.h
it2j0y0 -- it2j0y0_wrap: d*dd->i                           -- specfun_wrappers.h
iti0k0 -- it1i0k0_wrap: d*dd->i                            -- specfun_wrappers.h
it2i0k0 -- it2i0k0_wrap: d*dd->i                           -- specfun_wrappers.h
j0 -- j0: d->d                                             -- cephes.h
y0 -- y0: d->d                                             -- cephes.h
j1 -- j1: d->d                                             -- cephes.h
y1 -- y1: d->d                                             -- cephes.h
jv -- jv: dd->d, cbesj_wrap: dD->D                         -- cephes.h, amos_wrappers.h
jn -- jv: dd->d, cbesj_wrap: dD->D                         -- cephes.h, amos_wrappers.h
jve -- cbesj_wrap_e: dD->D, cbesj_wrap_e_real: dd->d       -- amos_wrappers.h
yv -- yv: dd->d, cbesy_wrap: dD->D                         -- cephes.h, amos_wrappers.h
yve -- cbesy_wrap_e_real: dd->d, cbesy_wrap_e: dD->D       -- amos_wrappers.h
k0 -- k0: d->d                                             -- cephes.h
k0e -- k0e: d->d                                           -- cephes.h
k1 -- k1: d->d                                             -- cephes.h
k1e -- k1e: d->d                                           -- cephes.h
kv -- cbesk_wrap: dD->D, cbesk_wrap_real: dd->d            -- amos_wrappers.h
kve -- cbesk_wrap_e_real: dd->d, cbesk_wrap_e: dD->D       -- amos_wrappers.h
hankel1 -- cbesh_wrap1: dD->D                              -- amos_wrappers.h
hankel1e -- cbesh_wrap1_e: dD->D                           -- amos_wrappers.h
hankel2 -- cbesh_wrap2: dD->D                              -- amos_wrappers.h
hankel2e -- cbesh_wrap2_e: dD->D                           -- amos_wrappers.h
ndtr -- ndtr: d->d                                         -- cephes.h
log_ndtr -- log_ndtr: d->d                                 -- cephes.h
erfc -- erfc: d->d                                         -- cephes.h
erf -- erf: d->d, cerf_wrap: D->D                          -- cephes.h, specfun_wrappers.h
ndtri -- ndtri: d->d                                       -- cephes.h
psi -- psi: d->d, cpsi_wrap: D->D                          -- cephes.h, specfun_wrappers.h
rgamma -- rgamma: d->d, crgamma_wrap: D->D                 -- cephes.h, specfun_wrappers.h
round -- round: d->d                                       -- cephes.h
sindg -- sindg: d->d                                       -- cephes.h
cosdg -- cosdg: d->d                                       -- cephes.h
radian -- radian: ddd->d                                   -- cephes.h
tandg -- tandg: d->d                                       -- cephes.h
cotdg -- cotdg: d->d                                       -- cephes.h
log1p -- log1p: d->d                                       -- cephes.h
expm1 -- expm1: d->d                                       -- cephes.h
cosm1 -- cosm1: d->d                                       -- cephes.h
spence -- spence: d->d                                     -- cephes.h
zetac -- zetac: d->d                                       -- cephes.h
struve -- struve_wrap: dd->d                               -- specfun_wrappers.h
modstruve -- modstruve_wrap: dd->d                         -- specfun_wrappers.h
itstruve0 -- itstruve0_wrap: d->d                          -- specfun_wrappers.h
it2struve0 -- it2struve0_wrap: d->d                        -- specfun_wrappers.h
itmodstruve0 -- itmodstruve0_wrap: d->d                    -- specfun_wrappers.h
kelvin -- kelvin_wrap: d*DDDD->i                           -- specfun_wrappers.h
ber -- ber_wrap: d->d                                      -- specfun_wrappers.h
bei -- bei_wrap: d->d                                      -- specfun_wrappers.h
ker -- ker_wrap: d->d                                      -- specfun_wrappers.h
kei -- kei_wrap: d->d                                      -- specfun_wrappers.h
berp -- berp_wrap: d->d                                    -- specfun_wrappers.h
beip -- beip_wrap: d->d                                    -- specfun_wrappers.h
kerp -- kerp_wrap: d->d                                    -- specfun_wrappers.h
keip -- keip_wrap: d->d                                    -- specfun_wrappers.h
zeta -- zeta: dd->d                                        -- cephes.h
kolmogorov -- kolmogorov: d->d                             -- cephes.h
kolmogi -- kolmogi: d->d                                   -- cephes.h
wofz -- cwofz_wrap: D->D                                   -- toms_wrappers.h
besselpoly -- besselpoly: ddd->d                           -- c_misc/misc.h
btdtria -- cdfbet3_wrap: ddd->d                            -- cdf_wrappers.h
btdtrib -- cdfbet4_wrap: ddd->d                            -- cdf_wrappers.h
bdtrik -- cdfbin2_wrap: ddd->d                             -- cdf_wrappers.h
bdtrin -- cdfbin3_wrap: ddd->d                             -- cdf_wrappers.h
chdtriv -- cdfchi3_wrap: dd->d                             -- cdf_wrappers.h
chndtr -- cdfchn1_wrap: ddd->d                             -- cdf_wrappers.h
chndtrix -- cdfchn2_wrap: ddd->d                           -- cdf_wrappers.h
chndtridf -- cdfchn3_wrap: ddd->d                          -- cdf_wrappers.h
chndtrinc -- cdfchn4_wrap: ddd->d                          -- cdf_wrappers.h
fdtridfd -- cdff4_wrap: ddd->d                             -- cdf_wrappers.h
ncfdtr -- cdffnc1_wrap: dddd->d                            -- cdf_wrappers.h
ncfdtri -- cdffnc2_wrap: dddd->d                           -- cdf_wrappers.h
ncfdtridfn -- cdffnc3_wrap: dddd->d                        -- cdf_wrappers.h
ncfdtridfd -- cdffnc4_wrap: dddd->d                        -- cdf_wrappers.h
ncfdtrinc -- cdffnc5_wrap: dddd->d                         -- cdf_wrappers.h
gdtrix -- cdfgam2_wrap: ddd->d                             -- cdf_wrappers.h
gdtrib -- cdfgam3_wrap: ddd->d                             -- cdf_wrappers.h
gdtria -- cdfgam4_wrap: ddd->d                             -- cdf_wrappers.h
nbdtrik -- cdfnbn2_wrap: ddd->d                            -- cdf_wrappers.h
nbdtrin -- cdfnbn3_wrap: ddd->d                            -- cdf_wrappers.h
nrdtrimn -- cdfnor3_wrap: ddd->d                           -- cdf_wrappers.h
nrdtrisd -- cdfnor4_wrap: ddd->d                           -- cdf_wrappers.h
pdtrik -- cdfpoi2_wrap: dd->d                              -- cdf_wrappers.h
stdtr -- cdft1_wrap: dd->d                                 -- cdf_wrappers.h
stdtrit -- cdft2_wrap: dd->d                               -- cdf_wrappers.h
stdtridf -- cdft3_wrap: dd->d                              -- cdf_wrappers.h
nctdtr -- cdftnc1_wrap: ddd->d                             -- cdf_wrappers.h
nctdtrit -- cdftnc2_wrap: ddd->d                           -- cdf_wrappers.h
nctdtridf -- cdftnc3_wrap: ddd->d                          -- cdf_wrappers.h
nctdtrinc -- cdftnc4_wrap: ddd->d                          -- cdf_wrappers.h
tklmbda -- tukeylambdacdf: dd->d                           -- cdf_wrappers.h
mathieu_a -- cem_cva_wrap: dd->d                           -- specfun_wrappers.h
mathieu_b -- sem_cva_wrap: dd->d                           -- specfun_wrappers.h
mathieu_cem -- cem_wrap: ddd*dd->i                         -- specfun_wrappers.h
mathieu_sem -- sem_wrap: ddd*dd->i                         -- specfun_wrappers.h
mathieu_modcem1 -- mcm1_wrap: ddd*dd->i                    -- specfun_wrappers.h
mathieu_modcem2 -- mcm2_wrap: ddd*dd->i                    -- specfun_wrappers.h
mathieu_modsem1 -- msm1_wrap: ddd*dd->i                    -- specfun_wrappers.h
mathieu_modsem2 -- msm2_wrap: ddd*dd->i                    -- specfun_wrappers.h
lpmv -- pmv_wrap: ddd->d                                   -- specfun_wrappers.h
pbwa -- pbwa_wrap: dd*dd->i                                -- specfun_wrappers.h
pbdv -- pbdv_wrap: dd*dd->i                                -- specfun_wrappers.h
pbvv -- pbvv_wrap: dd*dd->i                                -- specfun_wrappers.h
pro_cv -- prolate_segv_wrap: ddd->d                        -- specfun_wrappers.h
obl_cv -- oblate_segv_wrap: ddd->d                         -- specfun_wrappers.h
pro_ang1_cv -- prolate_aswfa_wrap: ddddd*dd->i             -- specfun_wrappers.h
pro_rad1_cv -- prolate_radial1_wrap: ddddd*dd->i           -- specfun_wrappers.h
pro_rad2_cv -- prolate_radial2_wrap: ddddd*dd->i           -- specfun_wrappers.h
obl_ang1_cv -- oblate_aswfa_wrap: ddddd*dd->i              -- specfun_wrappers.h
obl_rad1_cv -- oblate_radial1_wrap: ddddd*dd->i            -- specfun_wrappers.h
obl_rad2_cv -- oblate_radial2_wrap: ddddd*dd->i            -- specfun_wrappers.h
pro_ang1 -- prolate_aswfa_nocv_wrap: dddd*dd->i            -- specfun_wrappers.h
pro_rad1 -- prolate_radial1_nocv_wrap: dddd*dd->i          -- specfun_wrappers.h
pro_rad2 -- prolate_radial2_nocv_wrap: dddd*dd->i          -- specfun_wrappers.h
obl_ang1 -- oblate_aswfa_nocv_wrap: dddd*dd->i             -- specfun_wrappers.h
obl_rad1 -- oblate_radial1_nocv_wrap: dddd*dd->i           -- specfun_wrappers.h
obl_rad2 -- oblate_radial2_nocv_wrap: dddd*dd->i           -- specfun_wrappers.h
modfresnelp -- modified_fresnel_plus_wrap: d*DD->i         -- specfun_wrappers.h
modfresnelm -- modified_fresnel_minus_wrap: d*DD->i        -- specfun_wrappers.h
"""

# Ufuncs with C++
UFUNCS_CXX = """
"""

#---------------------------------------------------------------------------------
# Extra code (error handling system)
#---------------------------------------------------------------------------------

EXTRA_CODE = """

cdef extern int scipy_special_print_error_messages

cdef extern from "stdarg.h":
     ctypedef struct va_list:
         pass
     ctypedef struct fake_type:
         pass
     void va_start(va_list, void* arg) nogil
     void* va_arg(va_list, fake_type) nogil
     void va_end(va_list) nogil
     fake_type int_type "int"

cdef extern from "Python.h":
    object PyErr_Warn(object, char *)
    int PyOS_vsnprintf(char *, size_t, char *, va_list va) nogil

cdef public void scipy_special_raise_warning(char *fmt, ...) nogil:
    cdef char msg[1024]
    cdef va_list ap

    va_start(ap, fmt)
    PyOS_vsnprintf(msg, 1024, fmt, ap)
    va_end(ap)

    with gil:
        from scipy.special import SpecialFunctionWarning
        PyErr_Warn(SpecialFunctionWarning, msg)

def errprint(inflag=None):
    \"\"\"
    errprint(flag)

    Sets the error printing flag for special functions (from the
    cephesmodule). The output is the previous state.  With errprint(0)
    no error messages are shown; the default is errprint(1).  If no
    argument is given the current state of the flag is returned and no
    change occurs.

    \"\"\"
    global scipy_special_print_error_messages
    cdef int oldflag

    oldflag = scipy_special_print_error_messages
    if inflag is not None:
        scipy_special_print_error_messages = int(bool(inflag))

    return oldflag

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
    yield inputs.replace('d', 'f').replace('D', 'F'), outputs.replace('d', 'f').replace('D', 'F')
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

    lines = ufunc_str.splitlines()
    lines.sort()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        m = re.match("^([a-z0-9_]+)\s*--\s*(.*?)\s*(--.*)?$", line)
        if not m:
            raise ValueError("Unparseable line %r" % line)
        ufuncs.append(Ufunc(m.group(1), m.group(2)))
        if m.group(3):
            headers[ufuncs[-1].name] = [x.strip() for x in m.group(3)[2:].split(",")]

    toplevel = ""
    defs = ""
    all_loops = {}

    ufuncs.sort(key=lambda u: u.name)
    for ufunc in ufuncs:
        t, cfuncs = ufunc.generate(all_loops)
        toplevel += t + "\n"


        hdrs = headers.get(ufunc.name, ['cephes.h'])
        if len(hdrs) == 1:
            hdrs = [hdrs[0]] * len(cfuncs)
        elif len(hdrs) != len(cfuncs):
            raise ValueError("%s: wrong number of headers" % ufunc.name)

        for cfunc, header in zip(cfuncs, hdrs):
            if header.endswith('.pxd'):
                defs += "from %s cimport %s as _func_%s\n" % (header[:-4], cfunc, cfunc)
            else:
                defs += "cdef extern from \"%s\":\n" % header
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

    f.write(EXTRA_CODE)

    f.close()

def main():
    generate("_ufuncs.pyx", UFUNCS)
    generate("_ufuncs_cxx.pyx", UFUNCS_CXX)

    subprocess.call(['cython', '_ufuncs.pyx'])

if __name__ == "__main__":
    main()
