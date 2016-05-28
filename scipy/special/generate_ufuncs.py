#!/usr/bin/python
"""
generate_ufuncs.py

Generate Ufunc definition source files for scipy.special.  Produces
files '_ufuncs.c' and '_ufuncs_cxx.c' by first producing Cython.

This will generate both calls to PyUFunc_FromFuncAndData and the
required ufunc inner loops.

The syntax in the ufunc signature list is

    <line>:           <ufunc_name> '--' <kernels> '--' <headers>
    <kernels>:        <function> [',' <function>]*
    <function>:       <name> ':' <input> '*' <output>
                        '->' <retval> '*' <ignored_retval>
    <input>:          <typecode>*
    <output>:         <typecode>*
    <retval>:         <typecode>?
    <ignored_retval>: <typecode>?
    <headers>:        <header_name> [',' <header_name>]*

The input parameter types are denoted by single character type
codes, according to

   'f': 'float'
   'd': 'double'
   'g': 'long double'
   'F': 'float complex'
   'D': 'double complex'
   'G': 'long double complex'
   'i': 'int'
   'l': 'long'
   'v': 'void'

If multiple kernel functions are given for a single ufunc, the one
which is used is determined by the standard ufunc mechanism. Kernel
functions that are listed first are also matched first against the
ufunc input types, so functions listed earlier take precedence.

In addition, versions with casted variables, such as d->f,D->F and
i->d are automatically generated.

There should be either a single header that contains all of the kernel
functions listed, or there should be one header for each kernel
function. Cython pxd files are allowed in addition to .h files.

Cython functions may use fused types, but the names in the list
should be the specialized ones, such as 'somefunc[float]'.

Function coming from C++ should have ``++`` appended to the name of
the header.

Floating-point exceptions inside these Ufuncs are converted to
special function errors --- which are separately controlled by the
user, and off by default, as they are usually not especially useful
for the user.


The C++ module
--------------
In addition to ``_ufuncs`` module, a second module ``_ufuncs_cxx`` is
generated. This module only exports function pointers that are to be
used when constructing some of the ufuncs in ``_ufuncs``. The function
pointers are exported via Cython's standard mechanism.

This mainly avoids build issues --- Python distutils has no way to
figure out what to do if you want to link both C++ and Fortran code in
the same shared library.

"""

from __future__ import division, print_function, absolute_import

#---------------------------------------------------------------------------------
# Special function listing
#---------------------------------------------------------------------------------

#
#

# Ufuncs without C++
FUNCS = """
sph_harm -- sph_harmonic: iidd->D, sph_harmonic_unsafe: dddd->D -- sph_harm.pxd, _legacy.pxd
_lambertw -- lambertw_scalar: Dld->D                       -- lambertw.pxd
_ellip_harm -- ellip_harmonic: ddiiddd->d, ellip_harmonic_unsafe: ddddddd->d --_ellip_harm.pxd, _legacy.pxd
logit -- logitf: f->f, logit: d->d, logitl: g->g           -- _logit.h
expit -- expitf: f->f, expit: d->d, expitl: g->g           -- _logit.h
bdtrc -- bdtrc: iid->d, bdtrc_unsafe: ddd->d               -- cephes.h, _legacy.pxd
bdtr -- bdtr: iid->d, bdtr_unsafe: ddd->d                  -- cephes.h, _legacy.pxd
bdtri -- bdtri: iid->d, bdtri_unsafe: ddd->d               -- cephes.h, _legacy.pxd
binom -- binom: dd->d                                      -- orthogonal_eval.pxd
btdtr -- btdtr: ddd->d                                     -- cephes.h
btdtri -- incbi: ddd->d                                    -- cephes.h
fdtrc -- fdtrc: ddd->d                                     -- cephes.h
fdtr -- fdtr: ddd->d                                       -- cephes.h
fdtri -- fdtri: ddd->d                                     -- cephes.h
gdtrc -- gdtrc: ddd->d                                     -- cephes.h
gdtr -- gdtr: ddd->d                                       -- cephes.h
hyp0f1 -- _hyp0f1_real: dd->d, _hyp0f1_cmplx: dD->D        -- _hyp0f1.pxd
hyp2f1 -- hyp2f1: dddd->d, chyp2f1_wrap: dddD->D           -- cephes.h, specfun_wrappers.h
hyp1f1 -- hyp1f1_wrap: ddd->d, chyp1f1_wrap: ddD->D        -- specfun_wrappers.h
hyperu -- hypU_wrap: ddd->d                                -- specfun_wrappers.h
hyp2f0 -- hyp2f0: dddi*d->d, hyp2f0_unsafe: dddd*d->d      -- cephes.h, _legacy.pxd
hyp1f2 -- onef2: dddd*d->d                                 -- cephes.h
hyp3f0 -- threef0: dddd*d->d                               -- cephes.h
betainc -- incbet: ddd->d                                  -- cephes.h
betaincinv -- incbi: ddd->d                                -- cephes.h
nbdtrc -- nbdtrc: iid->d, nbdtrc_unsafe: ddd->d            -- cephes.h, _legacy.pxd
nbdtr -- nbdtr: iid->d, nbdtr_unsafe: ddd->d               -- cephes.h, _legacy.pxd
nbdtri -- nbdtri: iid->d, nbdtri_unsafe: ddd->d            -- cephes.h, _legacy.pxd
beta -- beta: dd->d                                        -- cephes.h
betaln -- lbeta: dd->d                                     -- cephes.h
cbrt -- cbrt: d->d                                         -- cephes.h
chdtrc -- chdtrc: dd->d                                    -- cephes.h
chdtr -- chdtr: dd->d                                      -- cephes.h
chdtri -- chdtri: dd->d                                    -- cephes.h
ellipeinc -- ellie: dd->d                                  -- cephes.h
ellipkinc -- ellik: dd->d                                  -- cephes.h
ellipe -- ellpe: d->d                                      -- cephes.h
ellipkm1 -- ellpk: d->d                                    -- cephes.h
eval_jacobi --      eval_jacobi[double]: dddd->d,     eval_jacobi[double complex]: dddD->D,     eval_jacobi_l: lddd->d -- orthogonal_eval.pxd
eval_sh_jacobi --   eval_sh_jacobi[double]: dddd->d,  eval_sh_jacobi[double complex]: dddD->D,  eval_sh_jacobi_l: lddd->d -- orthogonal_eval.pxd
eval_gegenbauer --  eval_gegenbauer[double]: ddd->d,  eval_gegenbauer[double complex]: ddD->D,  eval_gegenbauer_l: ldd->d -- orthogonal_eval.pxd
eval_chebyt --      eval_chebyt[double]: dd->d,       eval_chebyt[double complex]: dD->D,       eval_chebyt_l: ld->d -- orthogonal_eval.pxd
eval_chebyu --      eval_chebyu[double]: dd->d,       eval_chebyu[double complex]: dD->D,       eval_chebyu_l: ld->d -- orthogonal_eval.pxd
eval_chebyc --      eval_chebyc[double]: dd->d,       eval_chebyc[double complex]: dD->D,       eval_chebyc_l: ld->d -- orthogonal_eval.pxd
eval_chebys --      eval_chebys[double]: dd->d,       eval_chebys[double complex]: dD->D,       eval_chebys_l: ld->d -- orthogonal_eval.pxd
eval_sh_chebyt --   eval_sh_chebyt[double]: dd->d,    eval_sh_chebyt[double complex]: dD->D,    eval_sh_chebyt_l:ld->d -- orthogonal_eval.pxd
eval_sh_chebyu --   eval_sh_chebyu[double]: dd->d,    eval_sh_chebyu[double complex]: dD->D,    eval_sh_chebyu_l:ld->d -- orthogonal_eval.pxd
eval_legendre --    eval_legendre[double]: dd->d,     eval_legendre[double complex]: dD->D,     eval_legendre_l: ld->d -- orthogonal_eval.pxd
eval_sh_legendre -- eval_sh_legendre[double]: dd->d,  eval_sh_legendre[double complex]: dD->D,  eval_sh_legendre_l:ld->d -- orthogonal_eval.pxd
eval_genlaguerre -- eval_genlaguerre[double]: ddd->d, eval_genlaguerre[double complex]: ddD->D, eval_genlaguerre_l:ldd->d -- orthogonal_eval.pxd
eval_laguerre --    eval_laguerre[double]: dd->d,     eval_laguerre[double complex]: dD->D,     eval_laguerre_l:ld->d -- orthogonal_eval.pxd
eval_hermite  -- eval_hermite: ld->d                       -- orthogonal_eval.pxd
eval_hermitenorm -- eval_hermitenorm: ld->d                -- orthogonal_eval.pxd
exp10 -- exp10: d->d                                       -- cephes.h
exp2 -- exp2: d->d                                         -- cephes.h
gamma -- Gamma: d->d, cgamma: D->D                         -- cephes.h, _loggamma.pxd
_gammaln -- lgam: d->d, clngamma_wrap: D->D                -- cephes.h, specfun_wrappers.h
gammasgn -- gammasgn: d->d                                 -- c_misc/misc.h
i0 -- i0: d->d                                             -- cephes.h
i0e -- i0e: d->d                                           -- cephes.h
i1 -- i1: d->d                                             -- cephes.h
i1e -- i1e: d->d                                           -- cephes.h
gammaincc -- igamc: dd->d                                  -- cephes.h
gammainc -- igam: dd->d                                    -- cephes.h
gammaincinv -- gammaincinv: dd->d                          -- cephes.h
gammainccinv -- igami: dd->d                               -- cephes.h
iv -- iv: dd->d, cbesi_wrap: dD->D                         -- cephes.h, amos_wrappers.h
ive -- cbesi_wrap_e_real: dd->d, cbesi_wrap_e: dD->D       -- amos_wrappers.h
ellipj -- ellpj: dd*dddd->*i                               -- cephes.h
expn -- expn: id->d, expn_unsafe: dd->d                    -- cephes.h, _legacy.pxd
exp1 -- exp1_wrap: d->d, cexp1_wrap: D->D                  -- specfun_wrappers.h
expi -- expi_wrap: d->d, cexpi_wrap: D->D                  -- specfun_wrappers.h
kn -- cbesk_wrap_real_int: id->d, kn_unsafe: dd->d         -- cephes.h, _legacy.pxd
pdtrc -- pdtrc: id->d, pdtrc_unsafe: dd->d                 -- cephes.h, _legacy.pxd
pdtr -- pdtr: id->d, pdtr_unsafe: dd->d                    -- cephes.h, _legacy.pxd
pdtri -- pdtri: id->d, pdtri_unsafe: dd->d                 -- cephes.h, _legacy.pxd
yn -- yn: id->d, yn_unsafe: dd->d                          -- cephes.h, _legacy.pxd
smirnov -- smirnov: id->d, smirnov_unsafe: dd->d           -- cephes.h, _legacy.pxd
smirnovi -- smirnovi: id->d, smirnovi_unsafe: dd->d        -- cephes.h, _legacy.pxd
airy -- airy_wrap: d*dddd->*i, cairy_wrap: D*DDDD->*i      -- amos_wrappers.h
itairy -- itairy_wrap: d*dddd->*i                          -- specfun_wrappers.h
airye -- cairy_wrap_e_real: d*dddd->*i, cairy_wrap_e: D*DDDD->*i -- amos_wrappers.h
fresnel -- fresnl: d*dd->*i, cfresnl_wrap: D*DD->*i        -- cephes.h, specfun_wrappers.h
shichi -- shichi: d*dd->*i                                 -- cephes.h
sici -- sici: d*dd->*i                                     -- cephes.h
itj0y0 -- it1j0y0_wrap: d*dd->*i                           -- specfun_wrappers.h
it2j0y0 -- it2j0y0_wrap: d*dd->*i                          -- specfun_wrappers.h
iti0k0 -- it1i0k0_wrap: d*dd->*i                           -- specfun_wrappers.h
it2i0k0 -- it2i0k0_wrap: d*dd->*i                          -- specfun_wrappers.h
j0 -- j0: d->d                                             -- cephes.h
y0 -- y0: d->d                                             -- cephes.h
j1 -- j1: d->d                                             -- cephes.h
y1 -- y1: d->d                                             -- cephes.h
jv -- cbesj_wrap_real: dd->d, cbesj_wrap: dD->D            -- amos_wrappers.h
jve -- cbesj_wrap_e_real: dd->d, cbesj_wrap_e: dD->D       -- amos_wrappers.h
yv -- cbesy_wrap_real: dd->d, cbesy_wrap: dD->D            -- amos_wrappers.h
yve -- cbesy_wrap_e_real: dd->d, cbesy_wrap_e: dD->D       -- amos_wrappers.h
k0 -- k0: d->d                                             -- cephes.h
k0e -- k0e: d->d                                           -- cephes.h
k1 -- k1: d->d                                             -- cephes.h
k1e -- k1e: d->d                                           -- cephes.h
kv -- cbesk_wrap_real: dd->d, cbesk_wrap: dD->D            -- amos_wrappers.h
kve -- cbesk_wrap_e_real: dd->d, cbesk_wrap_e: dD->D       -- amos_wrappers.h
hankel1 -- cbesh_wrap1: dD->D                              -- amos_wrappers.h
hankel1e -- cbesh_wrap1_e: dD->D                           -- amos_wrappers.h
hankel2 -- cbesh_wrap2: dD->D                              -- amos_wrappers.h
hankel2e -- cbesh_wrap2_e: dD->D                           -- amos_wrappers.h
ndtr -- ndtr: d->d, faddeeva_ndtr: D->D                    -- cephes.h, _faddeeva.h++
log_ndtr -- log_ndtr: d->d, faddeeva_log_ndtr: D->D        -- cephes.h, _faddeeva.h++
ndtri -- ndtri: d->d                                       -- cephes.h
psi -- digamma: d->d, cdigamma: D->D                       -- _digamma.pxd, _digamma.pxd
rgamma -- rgamma: d->d, crgamma: D->D                      -- cephes.h, _loggamma.pxd
round -- round: d->d                                       -- cephes.h
sindg -- sindg: d->d                                       -- cephes.h
cosdg -- cosdg: d->d                                       -- cephes.h
radian -- radian: ddd->d                                   -- cephes.h
tandg -- tandg: d->d                                       -- cephes.h
cotdg -- cotdg: d->d                                       -- cephes.h
log1p -- log1p: d->d, clog1p: D->D                         -- cephes.h, _cunity.pxd
expm1 -- expm1: d->d, cexpm1: D->D                         -- cephes.h, _cunity.pxd
cosm1 -- cosm1: d->d                                       -- cephes.h
spence -- spence: d->d, cspence: D-> D                     -- cephes.h, _spence.pxd
zetac -- zetac: d->d                                       -- cephes.h
struve -- struve_h: dd->d                                  -- misc.h
modstruve -- struve_l: dd->d                               -- misc.h
_struve_power_series -- struve_power_series:  ddi*d->d     -- misc.h
_struve_asymp_large_z -- struve_asymp_large_z: ddi*d->d    -- misc.h
_struve_bessel_series -- struve_bessel_series: ddi*d->d    -- misc.h
itstruve0 -- itstruve0_wrap: d->d                          -- specfun_wrappers.h
it2struve0 -- it2struve0_wrap: d->d                        -- specfun_wrappers.h
itmodstruve0 -- itmodstruve0_wrap: d->d                    -- specfun_wrappers.h
kelvin -- kelvin_wrap: d*DDDD->*i                          -- specfun_wrappers.h
ber -- ber_wrap: d->d                                      -- specfun_wrappers.h
bei -- bei_wrap: d->d                                      -- specfun_wrappers.h
ker -- ker_wrap: d->d                                      -- specfun_wrappers.h
kei -- kei_wrap: d->d                                      -- specfun_wrappers.h
berp -- berp_wrap: d->d                                    -- specfun_wrappers.h
beip -- beip_wrap: d->d                                    -- specfun_wrappers.h
kerp -- kerp_wrap: d->d                                    -- specfun_wrappers.h
keip -- keip_wrap: d->d                                    -- specfun_wrappers.h
_zeta -- zeta: dd->d                                       -- cephes.h
kolmogorov -- kolmogorov: d->d                             -- cephes.h
kolmogi -- kolmogi: d->d                                   -- cephes.h
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
mathieu_cem -- cem_wrap: ddd*dd->*i                        -- specfun_wrappers.h
mathieu_sem -- sem_wrap: ddd*dd->*i                        -- specfun_wrappers.h
mathieu_modcem1 -- mcm1_wrap: ddd*dd->*i                   -- specfun_wrappers.h
mathieu_modcem2 -- mcm2_wrap: ddd*dd->*i                   -- specfun_wrappers.h
mathieu_modsem1 -- msm1_wrap: ddd*dd->*i                   -- specfun_wrappers.h
mathieu_modsem2 -- msm2_wrap: ddd*dd->*i                   -- specfun_wrappers.h
lpmv -- pmv_wrap: ddd->d                                   -- specfun_wrappers.h
pbwa -- pbwa_wrap: dd*dd->*i                               -- specfun_wrappers.h
pbdv -- pbdv_wrap: dd*dd->*i                               -- specfun_wrappers.h
pbvv -- pbvv_wrap: dd*dd->*i                               -- specfun_wrappers.h
pro_cv -- prolate_segv_wrap: ddd->d                        -- specfun_wrappers.h
obl_cv -- oblate_segv_wrap: ddd->d                         -- specfun_wrappers.h
pro_ang1_cv -- prolate_aswfa_wrap: ddddd*dd->*i            -- specfun_wrappers.h
pro_rad1_cv -- prolate_radial1_wrap: ddddd*dd->*i          -- specfun_wrappers.h
pro_rad2_cv -- prolate_radial2_wrap: ddddd*dd->*i          -- specfun_wrappers.h
obl_ang1_cv -- oblate_aswfa_wrap: ddddd*dd->*i             -- specfun_wrappers.h
obl_rad1_cv -- oblate_radial1_wrap: ddddd*dd->*i           -- specfun_wrappers.h
obl_rad2_cv -- oblate_radial2_wrap: ddddd*dd->*i           -- specfun_wrappers.h
pro_ang1 -- prolate_aswfa_nocv_wrap: dddd*d->d             -- specfun_wrappers.h
pro_rad1 -- prolate_radial1_nocv_wrap: dddd*d->d           -- specfun_wrappers.h
pro_rad2 -- prolate_radial2_nocv_wrap: dddd*d->d           -- specfun_wrappers.h
obl_ang1 -- oblate_aswfa_nocv_wrap: dddd*d->d              -- specfun_wrappers.h
obl_rad1 -- oblate_radial1_nocv_wrap: dddd*d->d            -- specfun_wrappers.h
obl_rad2 -- oblate_radial2_nocv_wrap: dddd*d->d            -- specfun_wrappers.h
modfresnelp -- modified_fresnel_plus_wrap: d*DD->*i        -- specfun_wrappers.h
modfresnelm -- modified_fresnel_minus_wrap: d*DD->*i       -- specfun_wrappers.h
wofz -- faddeeva_w: D->D                                   -- _faddeeva.h++
erfc -- erfc: d->d, faddeeva_erfc: D->D                    -- cephes.h, _faddeeva.h++
erf -- erf: d->d, faddeeva_erf: D->D                       -- cephes.h, _faddeeva.h++
dawsn -- faddeeva_dawsn: d->d, faddeeva_dawsn_complex: D->D -- _faddeeva.h++
erfcx -- faddeeva_erfcx: d->d, faddeeva_erfcx_complex: D->D -- _faddeeva.h++
erfi -- faddeeva_erfi: d->d, faddeeva_erfi_complex: D->D   -- _faddeeva.h++
xlogy -- xlogy[double]: dd->d, xlogy[double_complex]: DD->D -- _xlogy.pxd
xlog1py -- xlog1py[double]: dd->d, xlog1py[double_complex]: DD->D   -- _xlogy.pxd
poch -- poch: dd->d                                        -- c_misc/misc.h
boxcox -- boxcox: dd->d                                    -- _boxcox.pxd
boxcox1p -- boxcox1p: dd->d                                -- _boxcox.pxd
inv_boxcox -- inv_boxcox: dd->d                            -- _boxcox.pxd
inv_boxcox1p -- inv_boxcox1p: dd->d                        -- _boxcox.pxd
entr -- entr: d->d                                         -- _convex_analysis.pxd
kl_div -- kl_div: dd->d                                    -- _convex_analysis.pxd
rel_entr -- rel_entr: dd->d                                -- _convex_analysis.pxd
huber -- huber: dd->d                                      -- _convex_analysis.pxd
pseudo_huber -- pseudo_huber: dd->d                        -- _convex_analysis.pxd
exprel -- exprel: d->d                                     -- _exprel.pxd
_spherical_yn -- spherical_yn_real: ld->d, spherical_yn_complex: lD->D -- _spherical_bessel.pxd
_spherical_jn -- spherical_jn_real: ld->d, spherical_jn_complex: lD->D -- _spherical_bessel.pxd
_spherical_in -- spherical_in_real: ld->d, spherical_in_complex: lD->D -- _spherical_bessel.pxd
_spherical_kn -- spherical_kn_real: ld->d, spherical_kn_complex: lD->D -- _spherical_bessel.pxd
_spherical_yn_d -- spherical_yn_d_real: ld->d, spherical_yn_d_complex: lD->D -- _spherical_bessel.pxd
_spherical_jn_d -- spherical_jn_d_real: ld->d, spherical_jn_d_complex: lD->D -- _spherical_bessel.pxd
_spherical_in_d -- spherical_in_d_real: ld->d, spherical_in_d_complex: lD->D -- _spherical_bessel.pxd
_spherical_kn_d -- spherical_kn_d_real: ld->d, spherical_kn_d_complex: lD->D -- _spherical_bessel.pxd
loggamma -- loggamma: D->D                                 -- _loggamma.pxd
_sinpi -- sinpi[double]: d->d, sinpi[double_complex]: D->D -- _trig.pxd
_cospi -- cospi[double]: d->d, cospi[double_complex]: D->D -- _trig.pxd
"""

#---------------------------------------------------------------------------------
# Extra code
#---------------------------------------------------------------------------------

UFUNCS_EXTRA_CODE_COMMON = """\
# This file is automatically generated by generate_ufuncs.py.
# Do not edit manually!

cdef extern from "_complexstuff.h":
    # numpy/npy_math.h doesn't have correct extern "C" declarations,
    # so we must include a wrapped version first
    pass

cdef extern from "numpy/npy_math.h":
    double NPY_NAN

cimport numpy as np
from numpy cimport (
    npy_float, npy_double, npy_longdouble,
    npy_cfloat, npy_cdouble, npy_clongdouble,
    npy_int, npy_long,
    NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
    NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
    NPY_INT, NPY_LONG)

ctypedef double complex double_complex

cdef extern from "numpy/ufuncobject.h":
    int PyUFunc_getfperr() nogil

cdef public int wrap_PyUFunc_getfperr() nogil:
    \"\"\"
    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
    this avoids messing with the UNIQUE_SYMBOL #defines
    \"\"\"
    return PyUFunc_getfperr()

cimport libc

cimport sf_error

np.import_array()
np.import_ufunc()

cdef int _set_errprint(int flag) nogil:
    return sf_error.set_print(flag)
"""

UFUNCS_EXTRA_CODE = """
cimport scipy.special._ufuncs_cxx

def errprint(inflag=None):
    \"\"\"
    errprint(inflag=None)

    Sets or returns the error printing flag for special functions.

    Parameters
    ----------
    inflag : bool, optional
        Whether warnings concerning evaluation of special functions in
        scipy.special are shown. If omitted, no change is made to the
        current setting.

    Returns
    -------
    old_flag
        Previous value of the error flag

    \"\"\"
    if inflag is not None:
        scipy.special._ufuncs_cxx._set_errprint(int(bool(inflag)))
        return sf_error.set_print(int(bool(inflag)))
    else:
        return sf_error.get_print()
"""

UFUNCS_EXTRA_CODE_BOTTOM = """\
#
# Aliases
#
jn = jv
"""

CYTHON_SPECIAL_PXD = """\
# This file is automatically generated by generate_ufuncs.py.
# Do not edit manually!
"""

CYTHON_SPECIAL_PYX = '''\
# This file is automatically generated by generate_ufuncs.py.
# Do not edit manually!
"""
================================
Cython API for Special Functions
================================

Scalar, typed versions of many of the functions in ``scipy.special``
can be accessed directly from Cython; the complete list is given
below. The list is in the form::

    <function>: <inargs>-><outargs> [, <inargs>-><outargs> ...]

Here ``<inargs>`` and ``<outargs>`` are specified as lists of
typecodes which correspond to the following cython types:

- f: float
- d: double
- g: long double
- F: float complex
- D: double complex
- G: long double complex
- i: int
- l: long
- v: void.

Multiple pairs of ``<inargs>-><outargs>`` mean the function has
several possible signatures. For example, ``jv`` is listed as::

    jv: dd->d, dD->D,

which means that in Cython the following functions are available:

- ``double jv(double, double)``
- ``double complex jv(double, double complex)``.

The module is usable from Cython via::

    cimport scipy.special.cython_special

Available Functions
===================

FUNCLIST
"""
cimport numpy as np
from numpy cimport (
    npy_float, npy_double, npy_longdouble,
    npy_cfloat, npy_cdouble, npy_clongdouble,
    npy_int, npy_long,
    NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
    NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
    NPY_INT, NPY_LONG)

cdef extern from "numpy/ufuncobject.h":
    int PyUFunc_getfperr() nogil

cdef public int wrap_PyUFunc_getfperr() nogil:
    \"\"\"
    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
    this avoids messing with the UNIQUE_SYMBOL #defines
    \"\"\"
    return PyUFunc_getfperr()

cdef extern from "numpy/npy_math.h":
    double NPY_NAN

cimport sf_error
cimport _complexstuff
cimport scipy.special._ufuncs_cxx
from . import _ufuncs

ctypedef long double long_double
ctypedef float complex float_complex
ctypedef double complex double_complex
ctypedef long double complex long_double_complex
'''

CYTHON_SPECIAL_TEST = """\
# This file is automatically generated by generate_ufuncs.py.
# Do not edit manually!

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from scipy import special
from scipy.special import MODNAME

real_points = [-10, -1, 1, 10]
complex_points = [complex(*tup) for tup in product(real_points, repeat=2)]

TEST_POINTS = {
    'f': list(np.array(real_points, dtype=np.float32)),
    'd': list(np.array(real_points, dtype=np.float64)),
    'g': list(np.array(real_points, dtype=np.longdouble)),
    'F': list(np.array(complex_points, dtype=np.complex64)),
    'D': list(np.array(complex_points, dtype=np.complex128)),
    'G': list(np.array(complex_points, dtype=np.longcomplex)),
    'i': list(np.array(real_points, dtype=np.int)),
    'l': list(np.array(real_points, dtype=np.long)),
}


def _generate_test_points(typecodes):
    axes = tuple(map(lambda x: TEST_POINTS[x], typecodes))
    pts = list(product(*axes))
    return pts


def test_cython_api():
    params = [
PARAMS
    ]

    def check(param):
        pyfunc, cyfunc, specializations = param
        for typecodes in specializations:
            if "f" in typecodes or "F" in typecodes or "g" in typecodes or "G" in typecodes:
                # Don't test np.float32, np.complex64, np.longdouble,
                # or np.longcomplex
                continue
            pts = _generate_test_points(typecodes)
            pyres, cyres = [], []
            for pt in pts:
                pyres.append(pyfunc(*pt))
                cyres.append(cyfunc(*pt))
            pyres = np.asarray(pyres)
            cyres = np.asarray(cyres)
            assert_allclose(pyres, cyres)

    for param in params:
        yield check, param
"""

#---------------------------------------------------------------------------------
# Code generation
#---------------------------------------------------------------------------------

import os
import optparse
import re
import textwrap

add_newdocs = __import__('add_newdocs')

CY_TYPES = {
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

CY_UNDERSCORE_TYPES = {
    'f': 'float',
    'd': 'double',
    'g': 'long_double',
    'F': 'float_complex',
    'D': 'double_complex',
    'G': 'long_double_complex',
    'i': 'int',
    'l': 'long',
    'v': 'void',
}    

C_TYPES = {
    'f': 'npy_float',
    'd': 'npy_double',
    'g': 'npy_longdouble',
    'F': 'npy_cfloat',
    'D': 'npy_cdouble',
    'G': 'npy_clongdouble',
    'i': 'npy_int',
    'l': 'npy_long',
    'v': 'void',
}

TYPE_NAMES = {
    'f': 'NPY_FLOAT',
    'd': 'NPY_DOUBLE',
    'g': 'NPY_LONGDOUBLE',
    'F': 'NPY_CFLOAT',
    'D': 'NPY_CDOUBLE',
    'G': 'NPY_CLONGDOUBLE',
    'i': 'NPY_INT',
    'l': 'NPY_LONG',
}

CYTHON_SPECIAL_BENCH = {
    'jv': [('dd', 1, 1), ('dD', 1, 1)]
}


def cast_order(c):
    return ['ilfdgFDG'.index(x) for x in c]

# These downcasts will cause the function to return NaNs, unless the
# values happen to coincide exactly.
DANGEROUS_DOWNCAST = set([
    ('F', 'i'), ('F', 'l'), ('F', 'f'), ('F', 'd'), ('F', 'g'),
    ('D', 'i'), ('D', 'l'), ('D', 'f'), ('D', 'd'), ('D', 'g'),
    ('G', 'i'), ('G', 'l'), ('G', 'f'), ('G', 'd'), ('G', 'g'),
    ('f', 'i'), ('f', 'l'),
    ('d', 'i'), ('d', 'l'),
    ('g', 'i'), ('g', 'l'),
    ('l', 'i'),
])

NAN_VALUE = {
    'f': 'NPY_NAN',
    'd': 'NPY_NAN',
    'g': 'NPY_NAN',
    'F': 'NPY_NAN',
    'D': 'NPY_NAN',
    'G': 'NPY_NAN',
    'i': '0xbad0bad0',
    'l': '0xbad0bad0',
}


def generate_loop(func_inputs, func_outputs, func_retval,
                  ufunc_inputs, ufunc_outputs):
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
        codes are given according to the CY_TYPES and TYPE_NAMES
        lists above.

        The corresponding C function signature to be called is:

            retval func(intype1 iv1, intype2 iv2, ..., outtype1 *ov1, ...);

        If len(ufunc_outputs) == len(func_outputs)+1, the return value
        is treated as the first output argument. Otherwise, the return
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
    if len(func_inputs) != len(ufunc_inputs):
        raise ValueError("Function and ufunc have different number of inputs")

    if len(func_outputs) != len(ufunc_outputs) and not (
            func_retval != "v" and len(func_outputs)+1 == len(ufunc_outputs)):
        raise ValueError("Function retval and ufunc outputs don't match")

    name = "loop_%s_%s_%s_As_%s_%s" % (
        func_retval, func_inputs, func_outputs, ufunc_inputs, ufunc_outputs
        )
    body = "cdef void %s(char **args, np.npy_intp *dims, np.npy_intp *steps, void *data) nogil:\n" % name
    body += "    cdef np.npy_intp i, n = dims[0]\n"
    body += "    cdef void *func = (<void**>data)[0]\n"
    body += "    cdef char *func_name = <char*>(<void**>data)[1]\n"

    for j in range(len(ufunc_inputs)):
        body += "    cdef char *ip%d = args[%d]\n" % (j, j)
    for j in range(len(ufunc_outputs)):
        body += "    cdef char *op%d = args[%d]\n" % (j, j + len(ufunc_inputs))

    ftypes = []
    fvars = []
    outtypecodes = []
    for j in range(len(func_inputs)):
        ftypes.append(CY_TYPES[func_inputs[j]])
        fvars.append("<%s>(<%s*>ip%d)[0]" % (
            CY_TYPES[func_inputs[j]],
            CY_TYPES[ufunc_inputs[j]], j))

    if len(func_outputs)+1 == len(ufunc_outputs):
        func_joff = 1
        outtypecodes.append(func_retval)
        body += "    cdef %s ov0\n" % (CY_TYPES[func_retval],)
    else:
        func_joff = 0

    for j, outtype in enumerate(func_outputs):
        body += "    cdef %s ov%d\n" % (CY_TYPES[outtype], j+func_joff)
        ftypes.append("%s *" % CY_TYPES[outtype])
        fvars.append("&ov%d" % (j+func_joff))
        outtypecodes.append(outtype)

    body += "    for i in range(n):\n"
    if len(func_outputs)+1 == len(ufunc_outputs):
        rv = "ov0 = "
    else:
        rv = ""

    funcall = "        %s(<%s(*)(%s) nogil>func)(%s)\n" % (
        rv, CY_TYPES[func_retval], ", ".join(ftypes), ", ".join(fvars))

    # Cast-check inputs and call function
    input_checks = []
    for j in range(len(func_inputs)):
        if (ufunc_inputs[j], func_inputs[j]) in DANGEROUS_DOWNCAST:
            chk = "<%s>(<%s*>ip%d)[0] == (<%s*>ip%d)[0]" % (
                CY_TYPES[func_inputs[j]], CY_TYPES[ufunc_inputs[j]], j,
                CY_TYPES[ufunc_inputs[j]], j)
            input_checks.append(chk)

    if input_checks:
        body += "        if %s:\n" % (" and ".join(input_checks))
        body += "    " + funcall
        body += "        else:\n"
        body += "            sf_error.error(func_name, sf_error.DOMAIN, \"invalid input argument\")\n"
        for j, outtype in enumerate(outtypecodes):
            body += "            ov%d = <%s>%s\n" % (
                j, CY_TYPES[outtype], NAN_VALUE[outtype])
    else:
        body += funcall

    # Assign and cast-check output values
    for j, (outtype, fouttype) in enumerate(zip(ufunc_outputs, outtypecodes)):
        if (fouttype, outtype) in DANGEROUS_DOWNCAST:
            body += "        if ov%d == <%s>ov%d:\n" % (j, CY_TYPES[outtype], j)
            body += "            (<%s *>op%d)[0] = <%s>ov%d\n" % (
                CY_TYPES[outtype], j, CY_TYPES[outtype], j)
            body += "        else:\n"
            body += "            sf_error.error(func_name, sf_error.DOMAIN, \"invalid output\")\n"
            body += "            (<%s *>op%d)[0] = <%s>%s\n" % (
                CY_TYPES[outtype], j, CY_TYPES[outtype], NAN_VALUE[outtype])
        else:
            body += "        (<%s *>op%d)[0] = <%s>ov%d\n" % (
                CY_TYPES[outtype], j, CY_TYPES[outtype], j)
    for j in range(len(ufunc_inputs)):
        body += "        ip%d += steps[%d]\n" % (j, j)
    for j in range(len(ufunc_outputs)):
        body += "        op%d += steps[%d]\n" % (j, j + len(ufunc_inputs))

    body += "    sf_error.check_fpe(func_name)\n"

    return name, body


def generate_fused_type(typecodes):
    """
    Generate name of and cython code for a fused type.

    Parameters
    ----------
    typecodes : str
        Valid inputs to CY_TYPES (i.e. f, d, g, ...).

    """
    cytypes = map(lambda x: CY_TYPES[x], typecodes)
    name = typecodes + "_number_t"
    declaration = ["ctypedef fused " + name + ":"]
    for cytype in cytypes:
        declaration.append("    " + cytype)
    declaration = "\n".join(declaration)
    return name, cytypes, declaration


def generate_bench_func(name, pt):
    tab = " "*4
    typecodes, args = pt[0], pt[1:]
    template = ["def _bench_{}_{}_{}(int N):"]
    template.append(tab + "cdef int n")
    argnames = []
    for n, typecode in enumerate(typecodes):
        argname = "x{}".format(n)
        argnames.append(argname)
        line = tab
        line += "cdef {} {} = {}".format(CY_TYPES[typecode], argname, args[n])
        template.append(line)
    template.append(tab + "for n in range(N):")
    template.append(2*tab + "{}(" + ", ".join(argnames) + ")")
    template = "\n".join(template)
    py_bench = template.format(name, typecodes, "py", "_ufuncs." + name)
    cy_bench = template.format(name, typecodes, "cy", name)
    return py_bench, cy_bench


def iter_variants(inputs, outputs, always_long=True):
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
    if always_long:
        maps = [
            # always use long instead of int (more common type on 64-bit)
            ('i', 'l'),
        ]
    else:
        maps = [
            ('i', 'i'),
        ]

    # float32-preserving signatures
    if not ('i' in inputs or 'l' in inputs):
        # Don't add float32 versions of ufuncs with integer arguments, as this
        # can lead to incorrect dtype selection if the integer arguments are
        # arrays, but float arguments are scalars.
        # For instance sph_harm(0,[0],0,0).dtype == complex64
        # This may be a Numpy bug, but we need to work around it.
        # cf. gh-4895, https://github.com/numpy/numpy/issues/5895
        maps = maps + [(a + 'dD', b + 'fF') for a, b in maps]

    # do the replacements
    for src, dst in maps:
        new_inputs = inputs
        new_outputs = outputs
        for a, b in zip(src, dst):
            new_inputs = new_inputs.replace(a, b)
            new_outputs = new_outputs.replace(a, b)
        yield new_inputs, new_outputs


class Func(object):
    """
    Base class for Ufunc and FusedFunc.

    """
    def __init__(self, name, signatures, headers):
        self.name = name
        self.signatures = self._parse_signatures(signatures, headers)
        self.function_name_overrides = {}

    def _parse_signatures(self, sigs_str, headers_str):
        sigs = [x.strip() for x in sigs_str.split(",") if x.strip()]
        headers = [x.strip() for x in headers_str.split(",") if x.strip()]
        if len(headers) == 1:
            headers = headers * len(sigs)
        if len(headers) != len(sigs):
            raise ValueError("%s: Number of headers and signatures doesn't match: %r -- %r" % (
                self.name, sigs_str, headers_str))
        return [self._parse_signature(x) + (h,) for x, h in zip(sigs, headers)]

    def _parse_signature(self, sig):
        m = re.match("\s*(.*):\s*([fdgFDGil]*)\s*\\*\s*([fdgFDGil]*)\s*->\s*([*fdgFDGil]*)\s*$", sig)
        if m:
            func, inarg, outarg, ret = [x.strip() for x in m.groups()]
            if ret.count('*') > 1:
                raise ValueError("%s: Invalid signature: %r" % (self.name, sig))
            return (func, inarg, outarg, ret)
        m = re.match("\s*(.*):\s*([fdgFDGil]*)\s*->\s*([fdgFDGil]?)\s*$", sig)
        if m:
            func, inarg, ret = [x.strip() for x in m.groups()]
            return (func, inarg, "", ret)
        raise ValueError("%s: Invalid signature: %r" % (self.name, sig))

    def get_prototypes(self, nptypes_for_h=False):
        prototypes = []
        for func_name, inarg, outarg, ret, header in self.signatures:
            ret = ret.replace('*', '')
            c_args = ([C_TYPES[x] for x in inarg]
                      + [C_TYPES[x] + ' *' for x in outarg])
            cy_args = ([CY_TYPES[x] for x in inarg]
                       + [CY_TYPES[x] + ' *' for x in outarg])
            c_proto = "%s (*)(%s)" % (C_TYPES[ret], ", ".join(c_args))
            if header.endswith("h") and nptypes_for_h:
                cy_proto = c_proto + "nogil"
            else:
                cy_proto = "%s (*)(%s) nogil" % (CY_TYPES[ret], ", ".join(cy_args))
            prototypes.append((func_name, c_proto, cy_proto, header))
        return prototypes

    def cython_func_name(self, c_name, specialized=False, prefix="_func_",
                         override=True):
        # act on function name overrides
        if override and c_name in self.function_name_overrides:
            c_name = self.function_name_overrides[c_name]
            prefix = ""

        # support fused types
        m = re.match(r'^(.*?)(\[.*\])$', c_name)
        if m:
            c_base_name, fused_part = m.groups()
        else:
            c_base_name, fused_part = c_name, ""
        if specialized:
            return "%s%s%s" % (prefix, c_base_name, fused_part.replace(' ', '_'))
        else:
            return "%s%s" % (prefix, c_base_name,)

    @classmethod
    def parse_all(cls, ufunc_str):
        ufuncs = []

        lines = ufunc_str.splitlines()
        lines.sort()

        for line in lines:
            line = line.strip()
            if not line:
                continue
            m = re.match("^([a-z0-9_]+)\s*--\s*(.*?)\s*--(.*)$", line)
            if not m:
                raise ValueError("Unparseable line %r" % line)
            ufuncs.append(cls(m.group(1), m.group(2), m.group(3)))
        return ufuncs


class Ufunc(Func):
    """
    Ufunc signature, restricted format suitable for special functions.

    Parameters
    ----------
    name
        Name of the ufunc to create
    signature
        String of form 'func: fff*ff->f, func2: ddd->*i' describing
        the C-level functions and types of their input arguments
        and return values.

        The syntax is 'function_name: inputparams*outputparams->output_retval*ignored_retval'

    Attributes
    ----------
    name : str
        Python name for the Ufunc
    signatures : list of (func_name, inarg_spec, outarg_spec, ret_spec, header_name)
        List of parsed signatures
    doc : str
        Docstring, obtained from add_newdocs
    function_name_overrides : dict of str->str
        Overrides for the function names in signatures

    """
    def __init__(self, name, signatures, headers):
        super(Ufunc, self).__init__(name, signatures, headers)
        self.doc = add_newdocs.get("scipy.special." + name)
        if self.doc is None:
            raise ValueError("No docstring for ufunc %r" % name)
        self.doc = textwrap.dedent(self.doc).strip()

    def _get_signatures_and_loops(self, all_loops):
        inarg_num = None
        outarg_num = None

        seen = set()
        variants = []

        def add_variant(func_name, inarg, outarg, ret, inp, outp):
            if inp in seen:
                return
            seen.add(inp)

            sig = (func_name, inp, outp)
            if "v" in outp:
                raise ValueError("%s: void signature %r" % (self.name, sig))
            if len(inp) != inarg_num or len(outp) != outarg_num:
                raise ValueError("%s: signature %r does not have %d/%d input/output args" % (
                    self.name, sig,
                    inarg_num, outarg_num))

            loop_name, loop = generate_loop(inarg, outarg, ret, inp, outp)
            all_loops[loop_name] = loop
            variants.append((func_name, loop_name, inp, outp))

        # First add base variants
        for func_name, inarg, outarg, ret, header in self.signatures:
            outp = re.sub(r'\*.*', '', ret) + outarg
            ret = ret.replace('*', '')
            if inarg_num is None:
                inarg_num = len(inarg)
                outarg_num = len(outp)

            inp, outp = list(iter_variants(inarg, outp))[0]
            add_variant(func_name, inarg, outarg, ret, inp, outp)

        # Then the supplementary ones
        for func_name, inarg, outarg, ret, header in self.signatures:
            outp = re.sub(r'\*.*', '', ret) + outarg
            ret = ret.replace('*', '')
            for inp, outp in iter_variants(inarg, outp):
                add_variant(func_name, inarg, outarg, ret, inp, outp)

        # Then sort variants to input argument cast order
        # -- the sort is stable, so functions earlier in the signature list
        #    are still preferred
        variants.sort(key=lambda v: cast_order(v[2]))

        return variants, inarg_num, outarg_num

    def generate(self, all_loops):
        toplevel = ""

        variants, inarg_num, outarg_num = self._get_signatures_and_loops(all_loops)

        loops = []
        funcs = []
        types = []

        for func_name, loop_name, inputs, outputs in variants:
            for x in inputs:
                types.append(TYPE_NAMES[x])
            for x in outputs:
                types.append(TYPE_NAMES[x])
            loops.append(loop_name)
            funcs.append(func_name)

        toplevel += "cdef np.PyUFuncGenericFunction ufunc_%s_loops[%d]\n" % (self.name, len(loops))
        toplevel += "cdef void *ufunc_%s_ptr[%d]\n" % (self.name, 2*len(funcs))
        toplevel += "cdef void *ufunc_%s_data[%d]\n" % (self.name, len(funcs))
        toplevel += "cdef char ufunc_%s_types[%d]\n" % (self.name, len(types))
        toplevel += 'cdef char *ufunc_%s_doc = (\n    "%s")\n' % (
            self.name,
            self.doc.replace("\\", "\\\\").replace('"', '\\"').replace('\n', '\\n\"\n    "')
            )

        for j, function in enumerate(loops):
            toplevel += "ufunc_%s_loops[%d] = <np.PyUFuncGenericFunction>%s\n" % (self.name, j, function)
        for j, type in enumerate(types):
            toplevel += "ufunc_%s_types[%d] = <char>%s\n" % (self.name, j, type)
        for j, func in enumerate(funcs):
            toplevel += "ufunc_%s_ptr[2*%d] = <void*>%s\n" % (self.name, j,
                                                              self.cython_func_name(func, specialized=True))
            toplevel += "ufunc_%s_ptr[2*%d+1] = <void*>(<char*>\"%s\")\n" % (self.name, j,
                                                                             self.name)
        for j, func in enumerate(funcs):
            toplevel += "ufunc_%s_data[%d] = &ufunc_%s_ptr[2*%d]\n" % (
                self.name, j, self.name, j)

        toplevel += ('@ = np.PyUFunc_FromFuncAndData(ufunc_@_loops, '
                     'ufunc_@_data, ufunc_@_types, %d, %d, %d, 0, '
                     '"@", ufunc_@_doc, 0)\n' % (len(types)/(inarg_num+outarg_num),
                                                 inarg_num, outarg_num)
                     ).replace('@', self.name)

        return toplevel


class FusedFunc(Func):
    """
    Generate code for a fused-type special function that can be
    cimported in cython.

    """
    def __init__(self, name, signatures, headers):
        super(FusedFunc, self).__init__(name, signatures, headers)
        self.doc = "See the documentation for scipy.special." + self.name
        self.inargs, self.outargs = self._get_args()

    def _get_args(self):
        inarg_num, outarg_num = None, None
        all_inp, all_outp = [], []
        for func_name, inarg, outarg, ret, header in self.signatures:
            outp = re.sub(r'\*.*', '', ret) + outarg
            ret = ret.replace('*', '')
            if inarg_num is None:
                inarg_num = len(inarg)
                outarg_num = len(outp)
            inp, outp = list(iter_variants(inarg, outp, always_long=False))[0]
            all_inp.append(inp)
            all_outp.append(outp)

        inargs = []
        for n in range(inarg_num):
            types = unique(map(lambda x: x[n], all_inp))
            types.sort()
            inargs.append(''.join(types))
        outargs = []
        for n in range(outarg_num):
            types = unique(map(lambda x: x[n], all_outp))
            types.sort()
            outargs.append(''.join(types))

        return tuple(inargs), tuple(outargs)

    def generate(self):
        tab = " "*4

        # Get the types of the arguments and return values of the
        # function and create new fused types as necessary.
        fused_types = set()
        if len(self.outargs) > 1:
            raise ValueError("Only functions with one return value allowed.")
        else:
            outarg, = self.outargs
        if len(outarg) == 1:
            # It's not a fused type
            outtype = CY_TYPES[outarg]
        else:
            # It's a fused type
            outtype, _, type_declaration = generate_fused_type(outarg)
            fused_types.add(type_declaration)

        intypes = []
        for inarg in self.inargs:
            if len(inarg) == 1:
                cytype = CY_TYPES[inarg]
                intypes.append((cytype, [cytype]))
            else:
                intype, cytypes, type_declaration = generate_fused_type(inarg)
                fused_types.add(type_declaration)
                intypes.append((intype, cytypes))

        # Generate the body of the function. Do this first because we
        # won't know if we need a tmp variable until all the
        # signatures are parsed.
        needtmp = False
        multiple_specializations = False
        specializations = []
        body = []
        for n, (func_name, inargs, outargs, ret, header) in enumerate(self.signatures):
            # Generate the if-elif-else statements that select
            # specializations of the fused types.
            clauses = []
            seen = set()
            for (intype, cytypes), inarg in zip(intypes, inargs):
                if len(cytypes) == 1:
                    continue
                elif intype not in seen:
                    clauses.append("{} is {}".format(intype, CY_UNDERSCORE_TYPES[inarg]))
                    seen.add(intype)

            if clauses:
                multiple_specializations = True
                if n == 0:
                    line = tab + "if " + " and ".join(clauses) + ":"
                else:
                    line = tab + "elif " + " and ".join(clauses) + ":"
                body.append(line)
                sp = 2*tab
            else:
                sp = tab

            # If the function is from _ufuncs_cxx then it is exported
            # as a void * pointer; cast it to the correct type.
            ret = re.sub(r'\*.*', '', ret)
            if header.endswith("++"):
                call = "scipy.special._ufuncs_cxx._export_{}".format(func_name)
                argtypes = map(lambda x: CY_TYPES[x], inargs + outargs)
                rettype = CY_TYPES[ret]
                call = "(<{}(*)(".format(rettype) + ", ".join(argtypes) + ") nogil>" + call + ")"
            else:
                call = self.cython_func_name(func_name, specialized=True)
            call += "({})"

            # Get the arguments of the function, casting cython double
            # complex to npy_cdouble as necessary
            args = []
            for n, arg in enumerate(inargs):
                if header.endswith("h") and arg == "D":
                    args.append("_complexstuff.npy_cdouble_from_double_complex(x{})".format(n))
                else:
                    args.append("x{}".format(n))

            # Generate the call to a particular specialization
            if outargs and ret:
                raise ValueError("Only functions with one return value allowed.")
            elif ret:
                specializations.append(inargs + "->" + ret)
                if header.endswith("h") and ret == "D":
                    # Cast npy_cdouble to cython double complex
                    call = "_complexstuff.double_complex_from_npy_cdouble(" + call + ")"
                line = sp + "res = " + call.format(", ".join(args))
            elif outargs:
                specializations.append(inargs + "->" + outargs)
                if header.endswith("h") and outargs == "D":
                    needtmp = True
                    args += "&tmp"
                else:
                    args += "&res"
                line = sp + func_name + "(" + ", ".join(args) + ")"
            else:
                raise ValueError("Must return at least one value.")
            body.append(line)
            if needtmp:
                body.append(sp + "res = _complexstuff.npy_cdouble_from_double_complex(tmp)")
        if multiple_specializations:
            # If there was a specialization that we don't have a
            # signature for, fall back by returning NAN.
            body.append(tab + "else:")
            outarg, = self.outargs
            if len(outarg) == 1:
                body.append(2*tab + "res = {}".format(NAN_VALUE[outarg]))
            else:
                for n, t in enumerate(outarg):
                    if n == 0:
                        line = 2*tab + "if {} is {}:".format(outtype, CY_UNDERSCORE_TYPES[t])
                    elif n == len(outarg) - 1:
                        line = 2*tab + "else:"
                    else:
                        line = 2*tab + "elif {} is {}:".format(outtype, CY_UNDERSCORE_TYPES[t])
                    body.append(line)
                    body.append(3*tab + "res = {}".format(NAN_VALUE[t]))
        body.append(tab + "return res")

        # Generate the head of the function
        args, head = [], []
        for n, (intype, _) in enumerate(intypes):
            args.append("{} x{}".format(intype, n))
        declaration = "cpdef {} {}(".format(outtype, self.name) + ", ".join(args) + ") nogil"
        head.append(declaration + ":")
        head.append(tab + '"""' + self.doc + '"""')
        head.append(tab + "cdef {} res".format(outtype))
        if needtmp:
            head.append(tab + "cdef npy_cdouble tmp")

        source = "\n".join(head + body)
        return declaration, source, fused_types, specializations


def get_declaration(ufunc, c_name, c_proto, cy_proto, header, proto_h_filename):
    """
    Construct a Cython declaration of a function coming either from a
    pxd or a header file. Do sufficient tricks to enable compile-time
    type checking against the signature expected by the ufunc.
    """

    defs = []
    defs_h = []

    var_name = c_name.replace('[', '_').replace(']', '_').replace(' ', '_')

    if header.endswith('.pxd'):
        defs.append("from %s cimport %s as %s" % (
            header[:-4], ufunc.cython_func_name(c_name, prefix=""),
            ufunc.cython_func_name(c_name)))

        # check function signature at compile time
        proto_name = '_proto_%s_t' % var_name
        defs.append("ctypedef %s" % (cy_proto.replace('(*)', proto_name)))
        defs.append("cdef %s *%s_var = &%s" % (
            proto_name, proto_name, ufunc.cython_func_name(c_name, specialized=True)))
    else:
        # redeclare the function, so that the assumed
        # signature is checked at compile time
        new_name = "%s \"%s\"" % (ufunc.cython_func_name(c_name), c_name)
        defs.append("cdef extern from \"%s\":" % proto_h_filename)
        defs.append("    cdef %s" % (cy_proto.replace('(*)', new_name)))
        defs_h.append("#include \"%s\"" % header)
        defs_h.append("%s;" % (c_proto.replace('(*)', c_name)))

    return defs, defs_h, var_name


def generate_ufuncs(fn_prefix, cxx_fn_prefix, ufuncs):
    filename = fn_prefix + ".pyx"
    proto_h_filename = fn_prefix + '_defs.h'

    cxx_proto_h_filename = cxx_fn_prefix + '_defs.h'
    cxx_pyx_filename = cxx_fn_prefix + ".pyx"
    cxx_pxd_filename = cxx_fn_prefix + ".pxd"

    toplevel = ""

    # for _ufuncs*
    defs = []
    defs_h = []
    all_loops = {}

    # for _ufuncs_cxx*
    cxx_defs = []
    cxx_pxd_defs = ["cdef int _set_errprint(int flag) nogil"]
    cxx_defs_h = []

    ufuncs.sort(key=lambda u: u.name)

    for ufunc in ufuncs:
        # generate function declaration and type checking snippets
        cfuncs = ufunc.get_prototypes()
        for c_name, c_proto, cy_proto, header in cfuncs:
            if header.endswith('++'):
                header = header[:-2]

                # for the CXX module
                item_defs, item_defs_h, var_name = get_declaration(ufunc, c_name, c_proto, cy_proto,
                                                                   header, cxx_proto_h_filename)
                cxx_defs.extend(item_defs)
                cxx_defs_h.extend(item_defs_h)

                cxx_defs.append("cdef void *_export_%s = <void*>%s" % (
                    var_name, ufunc.cython_func_name(c_name, specialized=True, override=False)))
                cxx_pxd_defs.append("cdef void *_export_%s" % (var_name,))

                # let cython grab the function pointer from the c++ shared library
                ufunc.function_name_overrides[c_name] = "scipy.special._ufuncs_cxx._export_" + var_name
            else:
                # usual case
                item_defs, item_defs_h, _ = get_declaration(ufunc, c_name, c_proto, cy_proto, header,
                                                            proto_h_filename)
                defs.extend(item_defs)
                defs_h.extend(item_defs_h)

        # ufunc creation code snippet
        t = ufunc.generate(all_loops)
        toplevel += t + "\n"

    # Produce output
    toplevel = "\n".join(sorted(all_loops.values()) + defs + [toplevel])

    with open(filename, 'w') as f:
        f.write(UFUNCS_EXTRA_CODE_COMMON)
        f.write(UFUNCS_EXTRA_CODE)
        f.write(toplevel)
        f.write(UFUNCS_EXTRA_CODE_BOTTOM)

    defs_h = unique(defs_h)
    with open(proto_h_filename, 'w') as f:
        f.write("#ifndef UFUNCS_PROTO_H\n#define UFUNCS_PROTO_H 1\n")
        f.write("\n".join(defs_h))
        f.write("\n#endif\n")

    cxx_defs_h = unique(cxx_defs_h)
    with open(cxx_proto_h_filename, 'w') as f:
        f.write("#ifndef UFUNCS_PROTO_H\n#define UFUNCS_PROTO_H 1\n")
        f.write("\n".join(cxx_defs_h))
        f.write("\n#endif\n")

    with open(cxx_pyx_filename, 'w') as f:
        f.write(UFUNCS_EXTRA_CODE_COMMON)
        f.write("\n".join(cxx_defs))
        f.write("\n# distutils: language = c++\n")

    with open(cxx_pxd_filename, 'w') as f:
        f.write("\n".join(cxx_pxd_defs))


def generate_fused_funcs(modname, ufunc_fn_prefix, fused_funcs):
    pxdfile = modname + ".pxd"
    pyxfile = modname + ".pyx"
    proto_h_filename = ufunc_fn_prefix + '_defs.h'
    testfile = os.path.join(os.path.dirname(__file__), "tests",
                            "test_" + modname + ".py")

    sources = []
    declarations = []
    # Code for benchmarks
    bench = []
    fused_types = set()
    # Parameters for the tests
    params = []
    doc = []
    defs = []

    for func in fused_funcs:
        if func.name.startswith("_"):
            # Don't try to deal with functions that have extra layers
            # of wrappers.
            continue
        elif len(func.outargs) > 1:
            # Don't try to deal multiple return values
            continue
        # Get the function declaration for the .pxd and the source
        # code for the .pyx
        dec, source, func_fused_types, specs = func.generate()
        declarations.append(dec)
        sources.append(source)
        fused_types.update(func_fused_types)
        # Declare the specializations
        cfuncs = func.get_prototypes(nptypes_for_h=True)
        for c_name, c_proto, cy_proto, header in cfuncs:
            if header.endswith('++'):
                # We grab the c++ functions from the c++ module
                continue
            item_defs, _, _ = get_declaration(func, c_name, c_proto,
                                              cy_proto, header,
                                              proto_h_filename)
            defs.extend(item_defs)
        # Get a list of all the functions for the test generator
        pyfunc = "special." + func.name
        cyfunc = modname + "." + func.name
        inargs = str(tuple(map(lambda x: x.split("->")[0], specs)))
        param = [pyfunc, cyfunc, inargs]
        params.append("(" + ", ".join(param) + ")")
        # Add a line to the documentation
        doc.append('- :py:func:`~scipy.special.' + func.name + '`: ' + ", ".join(specs))
        # Add auxillary functions for benchmarking
        if func.name in CYTHON_SPECIAL_BENCH:
            for pt in CYTHON_SPECIAL_BENCH[func.name]:
                py_bench, cy_bench = generate_bench_func(func.name, pt)
                bench.extend([py_bench, cy_bench])
            
    fused_types = list(fused_types)
    fused_types.sort()

    with open(pxdfile, 'w') as f:
        f.write(CYTHON_SPECIAL_PXD)
        f.write("\n")
        f.write("\n\n".join(fused_types))
        f.write("\n\n")
        f.write("\n".join(declarations))
    with open(pyxfile, 'w') as f:
        header = CYTHON_SPECIAL_PYX
        header = header.replace("FUNCLIST", "\n".join(doc))
        f.write(header)
        f.write("\n\n")
        f.write("\n".join(defs))
        f.write("\n\n")
        f.write("\n\n".join(sources))
        f.write("\n\n")
        f.write("\n\n".join(bench))
    with open(testfile, 'w') as f:
        params = map(lambda x: " "*8 + x, params)
        params = ",\n".join(params)
        test = CYTHON_SPECIAL_TEST
        test = test.replace("MODNAME", modname).replace("PARAMS", params)
        f.write(test)


def unique(lst):
    """
    Return a list without repeated entries (first occurrence is kept),
    preserving order.
    """
    seen = set()
    new_lst = []
    for item in lst:
        if item in seen:
            continue
        seen.add(item)
        new_lst.append(item)
    return new_lst


def main():
    p = optparse.OptionParser(usage=__doc__.strip())
    options, args = p.parse_args()
    if len(args) != 0:
        p.error('invalid number of arguments')

    ufuncs = Ufunc.parse_all(FUNCS)
    generate_ufuncs("_ufuncs", "_ufuncs_cxx", ufuncs)
    fused_funcs = FusedFunc.parse_all(FUNCS)
    generate_fused_funcs("cython_special", "_ufuncs", fused_funcs)


if __name__ == "__main__":
    main()
