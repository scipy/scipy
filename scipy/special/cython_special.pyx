"""
.. highlight:: cython

Cython API for special functions
================================

Scalar, typed versions of many of the functions in ``scipy.special``
can be accessed directly from Cython; the complete list is given
below. Functions are overloaded using Cython fused types so their
names match their Python counterpart. The module follows the following
conventions:

- If a function's Python counterpart returns multiple values, then the
  function returns its outputs via pointers in the final arguments.
- If a function's Python counterpart returns a single value, then the
  function's output is returned directly.

The module is usable from Cython via::

    cimport scipy.special.cython_special

Error handling
--------------

Functions can indicate an error by returning ``nan``; however they
cannot emit warnings like their counterparts in ``scipy.special``.

Available functions
-------------------

- :py:func:`~scipy.special.voigt_profile`::

        double voigt_profile(double, double, double)

- :py:func:`~scipy.special.agm`::

        double agm(double, double)

- :py:func:`~scipy.special.airy`::

        void airy(double, double *, double *, double *, double *)
        void airy(double complex, double complex *, double complex *, double complex *, double complex *)

- :py:func:`~scipy.special.airye`::

        void airye(double complex, double complex *, double complex *, double complex *, double complex *)
        void airye(double, double *, double *, double *, double *)

- :py:func:`~scipy.special.bdtr`::

        double bdtr(double, double, double)
        double bdtr(double, long, double)

- :py:func:`~scipy.special.bdtrc`::

        double bdtrc(double, double, double)
        double bdtrc(double, long, double)

- :py:func:`~scipy.special.bdtri`::

        double bdtri(double, double, double)
        double bdtri(double, long, double)

- :py:func:`~scipy.special.bdtrik`::

        double bdtrik(double, double, double)

- :py:func:`~scipy.special.bdtrin`::

        double bdtrin(double, double, double)

- :py:func:`~scipy.special.bei`::

        double bei(double)

- :py:func:`~scipy.special.beip`::

        double beip(double)

- :py:func:`~scipy.special.ber`::

        double ber(double)

- :py:func:`~scipy.special.berp`::

        double berp(double)

- :py:func:`~scipy.special.besselpoly`::

        double besselpoly(double, double, double)

- :py:func:`~scipy.special.beta`::

        double beta(double, double)

- :py:func:`~scipy.special.betainc`::

        float betainc(float, float, float)
        double betainc(double, double, double)

- :py:func:`~scipy.special.betaincc`::

        float betaincc(float, float, float)
        double betaincc(double, double, double)

- :py:func:`~scipy.special.betaincinv`::

        float betaincinv(float, float, float)
        double betaincinv(double, double, double)

- :py:func:`~scipy.special.betainccinv`::

        float betainccinv(float, float, float)
        double betainccinv(double, double, double)

- :py:func:`~scipy.special.betaln`::

        double betaln(double, double)

- :py:func:`~scipy.special.binom`::

        double binom(double, double)

- :py:func:`~scipy.special.boxcox`::

        double boxcox(double, double)

- :py:func:`~scipy.special.boxcox1p`::

        double boxcox1p(double, double)

- :py:func:`~scipy.special.btdtr`::

        double btdtr(double, double, double)

- :py:func:`~scipy.special.btdtri`::

        double btdtri(double, double, double)

- :py:func:`~scipy.special.btdtria`::

        double btdtria(double, double, double)

- :py:func:`~scipy.special.btdtrib`::

        double btdtrib(double, double, double)

- :py:func:`~scipy.special.cbrt`::

        double cbrt(double)

- :py:func:`~scipy.special.chdtr`::

        double chdtr(double, double)

- :py:func:`~scipy.special.chdtrc`::

        double chdtrc(double, double)

- :py:func:`~scipy.special.chdtri`::

        double chdtri(double, double)

- :py:func:`~scipy.special.chdtriv`::

        double chdtriv(double, double)

- :py:func:`~scipy.special.chndtr`::

        double chndtr(double, double, double)

- :py:func:`~scipy.special.chndtridf`::

        double chndtridf(double, double, double)

- :py:func:`~scipy.special.chndtrinc`::

        double chndtrinc(double, double, double)

- :py:func:`~scipy.special.chndtrix`::

        double chndtrix(double, double, double)

- :py:func:`~scipy.special.cosdg`::

        double cosdg(double)

- :py:func:`~scipy.special.cosm1`::

        double cosm1(double)

- :py:func:`~scipy.special.cotdg`::

        double cotdg(double)

- :py:func:`~scipy.special.dawsn`::

        double dawsn(double)
        double complex dawsn(double complex)

- :py:func:`~scipy.special.ellipe`::

        double ellipe(double)

- :py:func:`~scipy.special.ellipeinc`::

        double ellipeinc(double, double)

- :py:func:`~scipy.special.ellipj`::

        void ellipj(double, double, double *, double *, double *, double *)

- :py:func:`~scipy.special.ellipkinc`::

        double ellipkinc(double, double)

- :py:func:`~scipy.special.ellipkm1`::

        double ellipkm1(double)

- :py:func:`~scipy.special.ellipk`::

        double ellipk(double)

- :py:func:`~scipy.special.elliprc`::

        double elliprc(double, double)
        double complex elliprc(double complex, double complex)

- :py:func:`~scipy.special.elliprd`::

        double elliprd(double, double, double)
        double complex elliprd(double complex, double complex, double complex)

- :py:func:`~scipy.special.elliprf`::

        double elliprf(double, double, double)
        double complex elliprf(double complex, double complex, double complex)

- :py:func:`~scipy.special.elliprg`::

        double elliprg(double, double, double)
        double complex elliprg(double complex, double complex, double complex)

- :py:func:`~scipy.special.elliprj`::

        double elliprj(double, double, double, double)
        double complex elliprj(double complex, double complex, double complex, double complex)

- :py:func:`~scipy.special.entr`::

        double entr(double)

- :py:func:`~scipy.special.erf`::

        double complex erf(double complex)
        double erf(double)

- :py:func:`~scipy.special.erfc`::

        double complex erfc(double complex)
        double erfc(double)

- :py:func:`~scipy.special.erfcx`::

        double erfcx(double)
        double complex erfcx(double complex)

- :py:func:`~scipy.special.erfi`::

        double erfi(double)
        double complex erfi(double complex)

- :py:func:`~scipy.special.erfinv`::

        float erfinv(float)
        double erfinv(double)

- :py:func:`~scipy.special.erfcinv`::

        double erfcinv(double)

- :py:func:`~scipy.special.eval_chebyc`::

        double complex eval_chebyc(double, double complex)
        double eval_chebyc(double, double)
        double eval_chebyc(long, double)

- :py:func:`~scipy.special.eval_chebys`::

        double complex eval_chebys(double, double complex)
        double eval_chebys(double, double)
        double eval_chebys(long, double)

- :py:func:`~scipy.special.eval_chebyt`::

        double complex eval_chebyt(double, double complex)
        double eval_chebyt(double, double)
        double eval_chebyt(long, double)

- :py:func:`~scipy.special.eval_chebyu`::

        double complex eval_chebyu(double, double complex)
        double eval_chebyu(double, double)
        double eval_chebyu(long, double)

- :py:func:`~scipy.special.eval_gegenbauer`::

        double complex eval_gegenbauer(double, double, double complex)
        double eval_gegenbauer(double, double, double)
        double eval_gegenbauer(long, double, double)

- :py:func:`~scipy.special.eval_genlaguerre`::

        double complex eval_genlaguerre(double, double, double complex)
        double eval_genlaguerre(double, double, double)
        double eval_genlaguerre(long, double, double)

- :py:func:`~scipy.special.eval_hermite`::

        double eval_hermite(long, double)

- :py:func:`~scipy.special.eval_hermitenorm`::

        double eval_hermitenorm(long, double)

- :py:func:`~scipy.special.eval_jacobi`::

        double complex eval_jacobi(double, double, double, double complex)
        double eval_jacobi(double, double, double, double)
        double eval_jacobi(long, double, double, double)

- :py:func:`~scipy.special.eval_laguerre`::

        double complex eval_laguerre(double, double complex)
        double eval_laguerre(double, double)
        double eval_laguerre(long, double)

- :py:func:`~scipy.special.eval_legendre`::

        double complex eval_legendre(double, double complex)
        double eval_legendre(double, double)
        double eval_legendre(long, double)

- :py:func:`~scipy.special.eval_sh_chebyt`::

        double complex eval_sh_chebyt(double, double complex)
        double eval_sh_chebyt(double, double)
        double eval_sh_chebyt(long, double)

- :py:func:`~scipy.special.eval_sh_chebyu`::

        double complex eval_sh_chebyu(double, double complex)
        double eval_sh_chebyu(double, double)
        double eval_sh_chebyu(long, double)

- :py:func:`~scipy.special.eval_sh_jacobi`::

        double complex eval_sh_jacobi(double, double, double, double complex)
        double eval_sh_jacobi(double, double, double, double)
        double eval_sh_jacobi(long, double, double, double)

- :py:func:`~scipy.special.eval_sh_legendre`::

        double complex eval_sh_legendre(double, double complex)
        double eval_sh_legendre(double, double)
        double eval_sh_legendre(long, double)

- :py:func:`~scipy.special.exp1`::

        double complex exp1(double complex)
        double exp1(double)

- :py:func:`~scipy.special.exp10`::

        double exp10(double)

- :py:func:`~scipy.special.exp2`::

        double exp2(double)

- :py:func:`~scipy.special.expi`::

        double complex expi(double complex)
        double expi(double)

- :py:func:`~scipy.special.expit`::

        double expit(double)
        float expit(float)
        long double expit(long double)

- :py:func:`~scipy.special.expm1`::

        double complex expm1(double complex)
        double expm1(double)

- :py:func:`~scipy.special.expn`::

        double expn(double, double)
        double expn(long, double)

- :py:func:`~scipy.special.exprel`::

        double exprel(double)

- :py:func:`~scipy.special.fdtr`::

        double fdtr(double, double, double)

- :py:func:`~scipy.special.fdtrc`::

        double fdtrc(double, double, double)

- :py:func:`~scipy.special.fdtri`::

        double fdtri(double, double, double)

- :py:func:`~scipy.special.fdtridfd`::

        double fdtridfd(double, double, double)

- :py:func:`~scipy.special.fresnel`::

        void fresnel(double, double *, double *)
        void fresnel(double complex, double complex *, double complex *)

- :py:func:`~scipy.special.gamma`::

        double complex gamma(double complex)
        double gamma(double)

- :py:func:`~scipy.special.gammainc`::

        double gammainc(double, double)

- :py:func:`~scipy.special.gammaincc`::

        double gammaincc(double, double)

- :py:func:`~scipy.special.gammainccinv`::

        double gammainccinv(double, double)

- :py:func:`~scipy.special.gammaincinv`::

        double gammaincinv(double, double)

- :py:func:`~scipy.special.gammaln`::

        double gammaln(double)

- :py:func:`~scipy.special.gammasgn`::

        double gammasgn(double)

- :py:func:`~scipy.special.gdtr`::

        double gdtr(double, double, double)

- :py:func:`~scipy.special.gdtrc`::

        double gdtrc(double, double, double)

- :py:func:`~scipy.special.gdtria`::

        double gdtria(double, double, double)

- :py:func:`~scipy.special.gdtrib`::

        double gdtrib(double, double, double)

- :py:func:`~scipy.special.gdtrix`::

        double gdtrix(double, double, double)

- :py:func:`~scipy.special.hankel1`::

        double complex hankel1(double, double complex)

- :py:func:`~scipy.special.hankel1e`::

        double complex hankel1e(double, double complex)

- :py:func:`~scipy.special.hankel2`::

        double complex hankel2(double, double complex)

- :py:func:`~scipy.special.hankel2e`::

        double complex hankel2e(double, double complex)

- :py:func:`~scipy.special.huber`::

        double huber(double, double)

- :py:func:`~scipy.special.hyp0f1`::

        double complex hyp0f1(double, double complex)
        double hyp0f1(double, double)

- :py:func:`~scipy.special.hyp1f1`::

        double hyp1f1(double, double, double)
        double complex hyp1f1(double, double, double complex)

- :py:func:`~scipy.special.hyp2f1`::

        double hyp2f1(double, double, double, double)
        double complex hyp2f1(double, double, double, double complex)

- :py:func:`~scipy.special.hyperu`::

        double hyperu(double, double, double)

- :py:func:`~scipy.special.i0`::

        double i0(double)

- :py:func:`~scipy.special.i0e`::

        double i0e(double)

- :py:func:`~scipy.special.i1`::

        double i1(double)

- :py:func:`~scipy.special.i1e`::

        double i1e(double)

- :py:func:`~scipy.special.inv_boxcox`::

        double inv_boxcox(double, double)

- :py:func:`~scipy.special.inv_boxcox1p`::

        double inv_boxcox1p(double, double)

- :py:func:`~scipy.special.it2i0k0`::

        void it2i0k0(double, double *, double *)

- :py:func:`~scipy.special.it2j0y0`::

        void it2j0y0(double, double *, double *)

- :py:func:`~scipy.special.it2struve0`::

        double it2struve0(double)

- :py:func:`~scipy.special.itairy`::

        void itairy(double, double *, double *, double *, double *)

- :py:func:`~scipy.special.iti0k0`::

        void iti0k0(double, double *, double *)

- :py:func:`~scipy.special.itj0y0`::

        void itj0y0(double, double *, double *)

- :py:func:`~scipy.special.itmodstruve0`::

        double itmodstruve0(double)

- :py:func:`~scipy.special.itstruve0`::

        double itstruve0(double)

- :py:func:`~scipy.special.iv`::

        double complex iv(double, double complex)
        double iv(double, double)

- :py:func:`~scipy.special.ive`::

        double complex ive(double, double complex)
        double ive(double, double)

- :py:func:`~scipy.special.j0`::

        double j0(double)

- :py:func:`~scipy.special.j1`::

        double j1(double)

- :py:func:`~scipy.special.jv`::

        double complex jv(double, double complex)
        double jv(double, double)

- :py:func:`~scipy.special.jve`::

        double complex jve(double, double complex)
        double jve(double, double)

- :py:func:`~scipy.special.k0`::

        double k0(double)

- :py:func:`~scipy.special.k0e`::

        double k0e(double)

- :py:func:`~scipy.special.k1`::

        double k1(double)

- :py:func:`~scipy.special.k1e`::

        double k1e(double)

- :py:func:`~scipy.special.kei`::

        double kei(double)

- :py:func:`~scipy.special.keip`::

        double keip(double)

- :py:func:`~scipy.special.kelvin`::

        void kelvin(double, double complex *, double complex *, double complex *, double complex *)

- :py:func:`~scipy.special.ker`::

        double ker(double)

- :py:func:`~scipy.special.kerp`::

        double kerp(double)

- :py:func:`~scipy.special.kl_div`::

        double kl_div(double, double)

- :py:func:`~scipy.special.kn`::

        double kn(double, double)
        double kn(long, double)

- :py:func:`~scipy.special.kolmogi`::

        double kolmogi(double)

- :py:func:`~scipy.special.kolmogorov`::

        double kolmogorov(double)

- :py:func:`~scipy.special.kv`::

        double complex kv(double, double complex)
        double kv(double, double)

- :py:func:`~scipy.special.kve`::

        double complex kve(double, double complex)
        double kve(double, double)

- :py:func:`~scipy.special.log1p`::

        double complex log1p(double complex)
        double log1p(double)

- :py:func:`~scipy.special.log_expit`::

        double log_expit(double)
        float log_expit(float)
        long double log_expit(long double)

- :py:func:`~scipy.special.log_ndtr`::

        double log_ndtr(double)
        double complex log_ndtr(double complex)

- :py:func:`~scipy.special.loggamma`::

        double loggamma(double)
        double complex loggamma(double complex)

- :py:func:`~scipy.special.logit`::

        double logit(double)
        float logit(float)
        long double logit(long double)

- :py:func:`~scipy.special.lpmv`::

        double lpmv(double, double, double)

- :py:func:`~scipy.special.mathieu_a`::

        double mathieu_a(double, double)

- :py:func:`~scipy.special.mathieu_b`::

        double mathieu_b(double, double)

- :py:func:`~scipy.special.mathieu_cem`::

        void mathieu_cem(double, double, double, double *, double *)

- :py:func:`~scipy.special.mathieu_modcem1`::

        void mathieu_modcem1(double, double, double, double *, double *)

- :py:func:`~scipy.special.mathieu_modcem2`::

        void mathieu_modcem2(double, double, double, double *, double *)

- :py:func:`~scipy.special.mathieu_modsem1`::

        void mathieu_modsem1(double, double, double, double *, double *)

- :py:func:`~scipy.special.mathieu_modsem2`::

        void mathieu_modsem2(double, double, double, double *, double *)

- :py:func:`~scipy.special.mathieu_sem`::

        void mathieu_sem(double, double, double, double *, double *)

- :py:func:`~scipy.special.modfresnelm`::

        void modfresnelm(double, double complex *, double complex *)

- :py:func:`~scipy.special.modfresnelp`::

        void modfresnelp(double, double complex *, double complex *)

- :py:func:`~scipy.special.modstruve`::

        double modstruve(double, double)

- :py:func:`~scipy.special.nbdtr`::

        double nbdtr(double, double, double)
        double nbdtr(long, long, double)

- :py:func:`~scipy.special.nbdtrc`::

        double nbdtrc(double, double, double)
        double nbdtrc(long, long, double)

- :py:func:`~scipy.special.nbdtri`::

        double nbdtri(double, double, double)
        double nbdtri(long, long, double)

- :py:func:`~scipy.special.nbdtrik`::

        double nbdtrik(double, double, double)

- :py:func:`~scipy.special.nbdtrin`::

        double nbdtrin(double, double, double)

- :py:func:`~scipy.special.ncfdtr`::

        double ncfdtr(double, double, double, double)

- :py:func:`~scipy.special.ncfdtri`::

        double ncfdtri(double, double, double, double)

- :py:func:`~scipy.special.ncfdtridfd`::

        double ncfdtridfd(double, double, double, double)

- :py:func:`~scipy.special.ncfdtridfn`::

        double ncfdtridfn(double, double, double, double)

- :py:func:`~scipy.special.ncfdtrinc`::

        double ncfdtrinc(double, double, double, double)

- :py:func:`~scipy.special.nctdtr`::

        double nctdtr(double, double, double)

- :py:func:`~scipy.special.nctdtridf`::

        double nctdtridf(double, double, double)

- :py:func:`~scipy.special.nctdtrinc`::

        double nctdtrinc(double, double, double)

- :py:func:`~scipy.special.nctdtrit`::

        double nctdtrit(double, double, double)

- :py:func:`~scipy.special.ndtr`::

        double complex ndtr(double complex)
        double ndtr(double)

- :py:func:`~scipy.special.ndtri`::

        double ndtri(double)

- :py:func:`~scipy.special.nrdtrimn`::

        double nrdtrimn(double, double, double)

- :py:func:`~scipy.special.nrdtrisd`::

        double nrdtrisd(double, double, double)

- :py:func:`~scipy.special.obl_ang1`::

        void obl_ang1(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.obl_ang1_cv`::

        void obl_ang1_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.obl_cv`::

        double obl_cv(double, double, double)

- :py:func:`~scipy.special.obl_rad1`::

        void obl_rad1(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.obl_rad1_cv`::

        void obl_rad1_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.obl_rad2`::

        void obl_rad2(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.obl_rad2_cv`::

        void obl_rad2_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.owens_t`::

        double owens_t(double, double)

- :py:func:`~scipy.special.pbdv`::

        void pbdv(double, double, double *, double *)

- :py:func:`~scipy.special.pbvv`::

        void pbvv(double, double, double *, double *)

- :py:func:`~scipy.special.pbwa`::

        void pbwa(double, double, double *, double *)

- :py:func:`~scipy.special.pdtr`::

        double pdtr(double, double)

- :py:func:`~scipy.special.pdtrc`::

        double pdtrc(double, double)

- :py:func:`~scipy.special.pdtri`::

        double pdtri(double, double)
        double pdtri(long, double)

- :py:func:`~scipy.special.pdtrik`::

        double pdtrik(double, double)

- :py:func:`~scipy.special.poch`::

        double poch(double, double)

- :py:func:`~scipy.special.powm1`::

        float powm1(float, float)
        double powm1(double, double)

- :py:func:`~scipy.special.pro_ang1`::

        void pro_ang1(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pro_ang1_cv`::

        void pro_ang1_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pro_cv`::

        double pro_cv(double, double, double)

- :py:func:`~scipy.special.pro_rad1`::

        void pro_rad1(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pro_rad1_cv`::

        void pro_rad1_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pro_rad2`::

        void pro_rad2(double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pro_rad2_cv`::

        void pro_rad2_cv(double, double, double, double, double, double *, double *)

- :py:func:`~scipy.special.pseudo_huber`::

        double pseudo_huber(double, double)

- :py:func:`~scipy.special.psi`::

        double complex psi(double complex)
        double psi(double)

- :py:func:`~scipy.special.radian`::

        double radian(double, double, double)

- :py:func:`~scipy.special.rel_entr`::

        double rel_entr(double, double)

- :py:func:`~scipy.special.rgamma`::

        double complex rgamma(double complex)
        double rgamma(double)

- :py:func:`~scipy.special.round`::

        double round(double)

- :py:func:`~scipy.special.shichi`::

        void shichi(double complex, double complex *, double complex *)
        void shichi(double, double *, double *)

- :py:func:`~scipy.special.sici`::

        void sici(double complex, double complex *, double complex *)
        void sici(double, double *, double *)

- :py:func:`~scipy.special.sindg`::

        double sindg(double)

- :py:func:`~scipy.special.smirnov`::

        double smirnov(double, double)
        double smirnov(long, double)

- :py:func:`~scipy.special.smirnovi`::

        double smirnovi(double, double)
        double smirnovi(long, double)

- :py:func:`~scipy.special.spence`::

        double complex spence(double complex)
        double spence(double)

- :py:func:`~scipy.special.sph_harm`::

        double complex sph_harm(double, double, double, double)
        double complex sph_harm(long, long, double, double)

- :py:func:`~scipy.special.stdtr`::

        double stdtr(double, double)

- :py:func:`~scipy.special.stdtridf`::

        double stdtridf(double, double)

- :py:func:`~scipy.special.stdtrit`::

        double stdtrit(double, double)

- :py:func:`~scipy.special.struve`::

        double struve(double, double)

- :py:func:`~scipy.special.tandg`::

        double tandg(double)

- :py:func:`~scipy.special.tklmbda`::

        double tklmbda(double, double)

- :py:func:`~scipy.special.wofz`::

        double complex wofz(double complex)

- :py:func:`~scipy.special.wrightomega`::

        double complex wrightomega(double complex)
        double wrightomega(double)

- :py:func:`~scipy.special.xlog1py`::

        double xlog1py(double, double)
        double complex xlog1py(double complex, double complex)

- :py:func:`~scipy.special.xlogy`::

        double xlogy(double, double)
        double complex xlogy(double complex, double complex)

- :py:func:`~scipy.special.y0`::

        double y0(double)

- :py:func:`~scipy.special.y1`::

        double y1(double)

- :py:func:`~scipy.special.yn`::

        double yn(double, double)
        double yn(long, double)

- :py:func:`~scipy.special.yv`::

        double complex yv(double, double complex)
        double yv(double, double)

- :py:func:`~scipy.special.yve`::

        double complex yve(double, double complex)
        double yve(double, double)

- :py:func:`~scipy.special.zetac`::

        double zetac(double)

- :py:func:`~scipy.special.wright_bessel`::

        double wright_bessel(double, double, double)

- :py:func:`~scipy.special.log_wright_bessel`::

        double log_wright_bessel(double, double, double)

- :py:func:`~scipy.special.ndtri_exp`::

        double ndtri_exp(double)


Custom functions
----------------

Some functions in ``scipy.special`` which are not ufuncs have custom
Cython wrappers.

Spherical Bessel functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional ``derivative`` boolean argument is replaced with an
optional Cython ``bint``, leading to the following signatures.

- :py:func:`~scipy.special.spherical_jn`::

        double complex spherical_jn(long, double complex)
        double complex spherical_jn(long, double complex, bint)
        double spherical_jn(long, double)
        double spherical_jn(long, double, bint)

- :py:func:`~scipy.special.spherical_yn`::

        double complex spherical_yn(long, double complex)
        double complex spherical_yn(long, double complex, bint)
        double spherical_yn(long, double)
        double spherical_yn(long, double, bint)

- :py:func:`~scipy.special.spherical_in`::

        double complex spherical_in(long, double complex)
        double complex spherical_in(long, double complex, bint)
        double spherical_in(long, double)
        double spherical_in(long, double, bint)

- :py:func:`~scipy.special.spherical_kn`::

        double complex spherical_kn(long, double complex)
        double complex spherical_kn(long, double complex, bint)
        double spherical_kn(long, double)
        double spherical_kn(long, double, bint)

"""

from libc.math cimport NAN

from numpy cimport npy_float, npy_double, npy_longdouble, npy_cdouble, npy_int, npy_long

cdef extern from "numpy/ufuncobject.h":
    int PyUFunc_getfperr() nogil

cdef public int wrap_PyUFunc_getfperr() noexcept nogil:
    """
    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
    this avoids messing with the UNIQUE_SYMBOL #defines
    """
    return PyUFunc_getfperr()

from . cimport _complexstuff
cimport scipy.special._ufuncs_cxx
from scipy.special import _ufuncs

ctypedef long double long_double
ctypedef float complex float_complex
ctypedef double complex double_complex
ctypedef long double complex long_double_complex

cdef extern from r"xsf_wrappers.h":
    double _func_gammaln_wrap "gammaln_wrap"(double) nogil

    double special_bei(double) nogil
    double special_beip(double) nogil
    double special_ber(double) nogil
    double special_berp(double) nogil
    npy_double special_kei(npy_double) nogil
    npy_double special_keip(npy_double) nogil
    void special_ckelvin(npy_double, npy_cdouble *, npy_cdouble *, npy_cdouble *, npy_cdouble *) nogil
    npy_double special_ker(npy_double) nogil
    double special_kerp(double) nogil
    npy_double _func_cem_cva_wrap "cem_cva_wrap"(npy_double, npy_double) nogil
    npy_double _func_sem_cva_wrap "sem_cva_wrap"(npy_double, npy_double) nogil
    void _func_cem_wrap "cem_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_mcm1_wrap "mcm1_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_mcm2_wrap "mcm2_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_msm1_wrap "msm1_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_msm2_wrap "msm2_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_sem_wrap "sem_wrap"(npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_modified_fresnel_minus_wrap "modified_fresnel_minus_wrap"(npy_double, npy_cdouble *, npy_cdouble *) nogil
    void _func_modified_fresnel_plus_wrap "modified_fresnel_plus_wrap"(npy_double, npy_cdouble *, npy_cdouble *) nogil
    npy_double _func_oblate_aswfa_nocv_wrap "oblate_aswfa_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_oblate_aswfa_wrap "oblate_aswfa_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_double _func_oblate_segv_wrap "oblate_segv_wrap"(npy_double, npy_double, npy_double) nogil
    npy_double _func_oblate_radial1_nocv_wrap "oblate_radial1_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_oblate_radial1_wrap "oblate_radial1_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_double _func_oblate_radial2_nocv_wrap "oblate_radial2_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_oblate_radial2_wrap "oblate_radial2_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_double _func_prolate_aswfa_nocv_wrap "prolate_aswfa_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_prolate_aswfa_wrap "prolate_aswfa_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_double _func_prolate_segv_wrap "prolate_segv_wrap"(npy_double, npy_double, npy_double) nogil
    npy_double _func_prolate_radial1_nocv_wrap "prolate_radial1_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_prolate_radial1_wrap "prolate_radial1_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_double _func_prolate_radial2_nocv_wrap "prolate_radial2_nocv_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double *) nogil
    void _func_prolate_radial2_wrap "prolate_radial2_wrap"(npy_double, npy_double, npy_double, npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_cdouble special_cexp1(npy_cdouble) nogil
    npy_double special_exp1(npy_double) nogil
    npy_cdouble special_cexpi(npy_cdouble) nogil
    npy_double special_expi(npy_double) nogil
    void _func_it2i0k0_wrap "it2i0k0_wrap"(npy_double, npy_double *, npy_double *) nogil
    void _func_it2j0y0_wrap "it2j0y0_wrap"(npy_double, npy_double *, npy_double *) nogil
    npy_double special_it2struve0(npy_double) nogil
    void special_itairy(npy_double, npy_double *, npy_double *, npy_double *, npy_double *) nogil
    void _func_it1i0k0_wrap "it1i0k0_wrap"(npy_double, npy_double *, npy_double *) nogil
    void _func_it1j0y0_wrap "it1j0y0_wrap"(npy_double, npy_double *, npy_double *) nogil
    npy_double special_itmodstruve0(npy_double) nogil
    npy_double special_itstruve0(npy_double) nogil
    void _func_pbdv_wrap "pbdv_wrap"(npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_pbvv_wrap "pbvv_wrap"(npy_double, npy_double, npy_double *, npy_double *) nogil
    void _func_pbwa_wrap "pbwa_wrap"(npy_double, npy_double, npy_double *, npy_double *) nogil
    npy_int _func_cfresnl_wrap "cfresnl_wrap"(npy_cdouble, npy_cdouble *, npy_cdouble *) nogil

    void special_airy(npy_double, npy_double *, npy_double *, npy_double *, npy_double *) nogil
    void special_cairy(npy_cdouble, npy_cdouble *, npy_cdouble *, npy_cdouble *, npy_cdouble *) nogil
    void special_airye(npy_double, npy_double *, npy_double *, npy_double *, npy_double *) nogil
    void special_cairye(npy_cdouble, npy_cdouble *, npy_cdouble *, npy_cdouble *, npy_cdouble *) nogil

    npy_cdouble special_ccyl_hankel_1(npy_double, npy_cdouble) nogil
    npy_cdouble special_ccyl_hankel_1e(npy_double, npy_cdouble) nogil
    npy_cdouble special_ccyl_hankel_2(npy_double, npy_cdouble) nogil
    npy_cdouble special_ccyl_hankel_2e(npy_double, npy_cdouble) nogil

    npy_double special_binom(npy_double, npy_double) nogil

    npy_double special_digamma(npy_double) nogil
    npy_cdouble special_cdigamma(npy_cdouble) nogil

    npy_double special_cyl_bessel_j(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_j(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_je(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_je(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_y(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_y(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_ye(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_ye(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_k_int(npy_int, npy_double) nogil
    npy_double special_cyl_bessel_k(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_k(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_ke(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_ke(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_i(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_i(npy_double, npy_cdouble) nogil

    npy_double special_cyl_bessel_ie(npy_double, npy_double) nogil
    npy_cdouble special_ccyl_bessel_ie(npy_double, npy_cdouble) nogil

    npy_double special_exprel(npy_double) nogil

    npy_double special_gamma(npy_double) nogil
    npy_cdouble special_cgamma(npy_cdouble) nogil

    npy_float special_expitf(npy_float) nogil
    npy_double special_expit(npy_double) nogil
    npy_longdouble special_expitl(npy_longdouble) nogil

    npy_float special_log_expitf(npy_float) nogil
    npy_double special_log_expit(npy_double) nogil
    npy_longdouble special_log_expitl(npy_longdouble) nogil

    npy_float special_logitf(npy_float) nogil
    npy_double special_logit(npy_double) nogil
    npy_longdouble special_logitl(npy_longdouble) nogil

    npy_double special_loggamma(npy_double) nogil
    npy_cdouble special_cloggamma(npy_cdouble) nogil

    npy_double special_hyp2f1(npy_double, npy_double, npy_double, npy_double) nogil
    npy_cdouble special_chyp2f1(npy_double, npy_double, npy_double, npy_cdouble) nogil

    npy_double special_rgamma(npy_double) nogil
    npy_cdouble special_crgamma(npy_cdouble) nogil

    npy_long special_sph_bessel_j(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_j(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_j_jac(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_j_jac(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_y(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_y(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_y_jac(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_y_jac(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_i(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_i(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_i_jac(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_i_jac(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_k(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_k(npy_long, npy_cdouble) nogil

    npy_long special_sph_bessel_k_jac(npy_long, npy_double) nogil
    npy_cdouble special_csph_bessel_k_jac(npy_long, npy_cdouble) nogil

    npy_cdouble special_sph_harm(npy_long, npy_long, npy_double, npy_double) nogil
    npy_cdouble special_sph_harm_unsafe(npy_double, npy_double, npy_double, npy_double) nogil
    double _func_cephes_iv_wrap "cephes_iv_wrap"(double, double) nogil

    npy_double special_wright_bessel(npy_double, npy_double, npy_double) nogil
    npy_double special_log_wright_bessel(npy_double, npy_double, npy_double) nogil
    double special_ellipk(double m) nogil

    double cephes_besselpoly(double a, double lmbda, double nu) nogil
    double cephes_beta(double a, double b) nogil
    double cephes_chdtr(double df, double x) nogil
    double cephes_chdtrc(double df, double x) nogil
    double cephes_chdtri(double df, double y) nogil
    double cephes_lbeta(double a, double b) nogil
    double cephes_sinpi(double x) nogil
    double cephes_cospi(double x) nogil
    double cephes_cbrt(double x) nogil
    double cephes_Gamma(double x) nogil
    double cephes_gammasgn(double x) nogil
    double cephes_hyp2f1(double a, double b, double c, double x) nogil
    double cephes_i0(double x) nogil
    double cephes_i0e(double x) nogil
    double cephes_i1(double x) nogil
    double cephes_i1e(double x) nogil
    double cephes_iv(double v, double x) nogil
    double cephes_j0(double x) nogil
    double cephes_j1(double x) nogil
    double cephes_k0(double x) nogil
    double cephes_k0e(double x) nogil
    double cephes_k1(double x) nogil
    double cephes_k1e(double x) nogil
    double cephes_y0(double x) nogil
    double cephes_y1(double x) nogil
    double cephes_yn(int n, double x) nogil
    double cephes_igam(double a, double x) nogil
    double cephes_igamc(double a, double x) nogil
    double cephes_igami(double a, double p) nogil
    double cephes_igamci(double a, double p) nogil
    double cephes_igam_fac(double a, double x) nogil
    double cephes_lanczos_sum_expg_scaled(double x) nogil
    double cephes_kolmogorov(double x) nogil
    double cephes_kolmogc(double x) nogil
    double cephes_kolmogi(double x) nogil
    double cephes_kolmogci(double x) nogil
    double cephes_kolmogp(double x) nogil
    double cephes_smirnov(int n, double x) nogil
    double cephes_smirnovc(int n, double x) nogil
    double cephes_smirnovi(int n, double x) nogil
    double cephes_smirnovci(int n, double x) nogil
    double cephes_smirnovp(int n, double x) nogil
    double cephes_ndtr(double x) nogil
    double cephes_erf(double x) nogil
    double cephes_erfc(double x) nogil
    double cephes_poch(double x, double m) nogil
    double cephes_rgamma(double x) nogil
    double cephes_zeta(double x, double q) nogil
    double cephes_zetac(double x) nogil
    double cephes_riemann_zeta(double x) nogil
    double cephes_log1p(double x) nogil
    double cephes_log1pmx(double x) nogil
    double cephes_lgam1p(double x) nogil
    double cephes_expm1(double x) nogil
    double cephes_cosm1(double x) nogil
    double cephes_expn(int n, double x) nogil
    double cephes_ellpe(double x) nogil
    double cephes_ellpk(double x) nogil
    double cephes_ellie(double phi, double m) nogil
    double cephes_ellik(double phi, double m) nogil
    double cephes_sindg(double x) nogil
    double cephes_cosdg(double x) nogil
    double cephes_tandg(double x) nogil
    double cephes_cotdg(double x) nogil
    double cephes_radian(double d, double m, double s) nogil
    double cephes_ndtri(double x) nogil
    double cephes_bdtr(double k, int n, double p) nogil
    double cephes_bdtri(double k, int n, double y) nogil
    double cephes_bdtrc(double k, int n, double p) nogil
    double cephes_btdtri(double aa, double bb, double yy0) nogil
    double cephes_btdtr(double a, double b, double x) nogil
    double cephes_erfcinv(double y) nogil
    double cephes_exp10(double x) nogil
    double cephes_exp2(double x) nogil
    double cephes_fdtr(double a, double b, double x) nogil
    double cephes_fdtrc(double a, double b, double x) nogil
    double cephes_fdtri(double a, double b, double y) nogil
    double cephes_gdtr(double a, double b, double x) nogil
    double cephes_gdtrc(double a, double b, double x) nogil
    double cephes_owens_t(double h, double a) nogil
    double cephes_nbdtr(int k, int n, double p) nogil
    double cephes_nbdtrc(int k, int n, double p) nogil
    double cephes_nbdtri(int k, int n, double p) nogil
    double cephes_pdtr(double k, double m) nogil
    double cephes_pdtrc(double k, double m) nogil
    double cephes_pdtri(int k, double y) nogil
    double cephes_round(double x) nogil
    double cephes_spence(double x) nogil

    double cephes_tukeylambdacdf(double x, double lmbda) nogil
    double cephes_struve_h(double v, double z) nogil
    double cephes_struve_l(double v, double z) nogil

from ._agm cimport agm as _func_agm
ctypedef double _proto_agm_t(double, double) noexcept nogil
cdef _proto_agm_t *_proto_agm_t_var = &_func_agm
from ._legacy cimport bdtr_unsafe as _func_bdtr_unsafe
ctypedef double _proto_bdtr_unsafe_t(double, double, double) noexcept nogil
cdef _proto_bdtr_unsafe_t *_proto_bdtr_unsafe_t_var = &_func_bdtr_unsafe
from ._legacy cimport bdtrc_unsafe as _func_bdtrc_unsafe
ctypedef double _proto_bdtrc_unsafe_t(double, double, double) noexcept nogil
cdef _proto_bdtrc_unsafe_t *_proto_bdtrc_unsafe_t_var = &_func_bdtrc_unsafe
from ._legacy cimport bdtri_unsafe as _func_bdtri_unsafe
ctypedef double _proto_bdtri_unsafe_t(double, double, double) noexcept nogil
cdef _proto_bdtri_unsafe_t *_proto_bdtri_unsafe_t_var = &_func_bdtri_unsafe
from ._cdflib_wrappers cimport bdtrik as _func_bdtrik
ctypedef double _proto_bdtrik_t(double, double, double) noexcept nogil
cdef _proto_bdtrik_t *_proto_bdtrik_t_var = &_func_bdtrik
from ._cdflib_wrappers cimport bdtrin as _func_bdtrin
ctypedef double _proto_bdtrin_t(double, double, double) noexcept nogil
cdef _proto_bdtrin_t *_proto_bdtrin_t_var = &_func_bdtrin
from ._boxcox cimport boxcox as _func_boxcox
ctypedef double _proto_boxcox_t(double, double) noexcept nogil
cdef _proto_boxcox_t *_proto_boxcox_t_var = &_func_boxcox
from ._boxcox cimport boxcox1p as _func_boxcox1p
ctypedef double _proto_boxcox1p_t(double, double) noexcept nogil
cdef _proto_boxcox1p_t *_proto_boxcox1p_t_var = &_func_boxcox1p
from ._cdflib_wrappers cimport btdtria as _func_btdtria
ctypedef double _proto_btdtria_t(double, double, double) noexcept nogil
cdef _proto_btdtria_t *_proto_btdtria_t_var = &_func_btdtria
from ._cdflib_wrappers cimport btdtrib as _func_btdtrib
ctypedef double _proto_btdtrib_t(double, double, double) noexcept nogil
cdef _proto_btdtrib_t *_proto_btdtrib_t_var = &_func_btdtrib
from ._cdflib_wrappers cimport chdtriv as _func_chdtriv
ctypedef double _proto_chdtriv_t(double, double) noexcept nogil
cdef _proto_chdtriv_t *_proto_chdtriv_t_var = &_func_chdtriv
from ._cdflib_wrappers cimport chndtr as _func_chndtr
ctypedef double _proto_chndtr_t(double, double, double) noexcept nogil
cdef _proto_chndtr_t *_proto_chndtr_t_var = &_func_chndtr
from ._cdflib_wrappers cimport chndtridf as _func_chndtridf
ctypedef double _proto_chndtridf_t(double, double, double) noexcept nogil
cdef _proto_chndtridf_t *_proto_chndtridf_t_var = &_func_chndtridf
from ._cdflib_wrappers cimport chndtrinc as _func_chndtrinc
ctypedef double _proto_chndtrinc_t(double, double, double) noexcept nogil
cdef _proto_chndtrinc_t *_proto_chndtrinc_t_var = &_func_chndtrinc
from ._cdflib_wrappers cimport chndtrix as _func_chndtrix
ctypedef double _proto_chndtrix_t(double, double, double) noexcept nogil
cdef _proto_chndtrix_t *_proto_chndtrix_t_var = &_func_chndtrix
cdef extern from r"_ufuncs_defs.h":
    cdef npy_int _func_cephes_ellpj_wrap "cephes_ellpj_wrap"(npy_double, npy_double, npy_double *, npy_double *, npy_double *, npy_double *)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_ellik "ellik"(npy_double, npy_double)nogil
from ._ellipk cimport ellipk as _func_ellipk
ctypedef double _proto_ellipk_t(double) noexcept nogil
cdef _proto_ellipk_t *_proto_ellipk_t_var = &_func_ellipk
from ._convex_analysis cimport entr as _func_entr
ctypedef double _proto_entr_t(double) noexcept nogil
cdef _proto_entr_t *_proto_entr_t_var = &_func_entr
from .orthogonal_eval cimport eval_chebyc as _func_eval_chebyc
ctypedef double complex _proto_eval_chebyc_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_chebyc_double_complex__t *_proto_eval_chebyc_double_complex__t_var = &_func_eval_chebyc[double_complex]
from .orthogonal_eval cimport eval_chebyc as _func_eval_chebyc
ctypedef double _proto_eval_chebyc_double__t(double, double) noexcept nogil
cdef _proto_eval_chebyc_double__t *_proto_eval_chebyc_double__t_var = &_func_eval_chebyc[double]
from .orthogonal_eval cimport eval_chebyc_l as _func_eval_chebyc_l
ctypedef double _proto_eval_chebyc_l_t(long, double) noexcept nogil
cdef _proto_eval_chebyc_l_t *_proto_eval_chebyc_l_t_var = &_func_eval_chebyc_l
from .orthogonal_eval cimport eval_chebys as _func_eval_chebys
ctypedef double complex _proto_eval_chebys_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_chebys_double_complex__t *_proto_eval_chebys_double_complex__t_var = &_func_eval_chebys[double_complex]
from .orthogonal_eval cimport eval_chebys as _func_eval_chebys
ctypedef double _proto_eval_chebys_double__t(double, double) noexcept nogil
cdef _proto_eval_chebys_double__t *_proto_eval_chebys_double__t_var = &_func_eval_chebys[double]
from .orthogonal_eval cimport eval_chebys_l as _func_eval_chebys_l
ctypedef double _proto_eval_chebys_l_t(long, double) noexcept nogil
cdef _proto_eval_chebys_l_t *_proto_eval_chebys_l_t_var = &_func_eval_chebys_l
from .orthogonal_eval cimport eval_chebyt as _func_eval_chebyt
ctypedef double complex _proto_eval_chebyt_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_chebyt_double_complex__t *_proto_eval_chebyt_double_complex__t_var = &_func_eval_chebyt[double_complex]
from .orthogonal_eval cimport eval_chebyt as _func_eval_chebyt
ctypedef double _proto_eval_chebyt_double__t(double, double) noexcept nogil
cdef _proto_eval_chebyt_double__t *_proto_eval_chebyt_double__t_var = &_func_eval_chebyt[double]
from .orthogonal_eval cimport eval_chebyt_l as _func_eval_chebyt_l
ctypedef double _proto_eval_chebyt_l_t(long, double) noexcept nogil
cdef _proto_eval_chebyt_l_t *_proto_eval_chebyt_l_t_var = &_func_eval_chebyt_l
from .orthogonal_eval cimport eval_chebyu as _func_eval_chebyu
ctypedef double complex _proto_eval_chebyu_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_chebyu_double_complex__t *_proto_eval_chebyu_double_complex__t_var = &_func_eval_chebyu[double_complex]
from .orthogonal_eval cimport eval_chebyu as _func_eval_chebyu
ctypedef double _proto_eval_chebyu_double__t(double, double) noexcept nogil
cdef _proto_eval_chebyu_double__t *_proto_eval_chebyu_double__t_var = &_func_eval_chebyu[double]
from .orthogonal_eval cimport eval_chebyu_l as _func_eval_chebyu_l
ctypedef double _proto_eval_chebyu_l_t(long, double) noexcept nogil
cdef _proto_eval_chebyu_l_t *_proto_eval_chebyu_l_t_var = &_func_eval_chebyu_l
from .orthogonal_eval cimport eval_gegenbauer as _func_eval_gegenbauer
ctypedef double complex _proto_eval_gegenbauer_double_complex__t(double, double, double complex) noexcept nogil
cdef _proto_eval_gegenbauer_double_complex__t *_proto_eval_gegenbauer_double_complex__t_var = &_func_eval_gegenbauer[double_complex]
from .orthogonal_eval cimport eval_gegenbauer as _func_eval_gegenbauer
ctypedef double _proto_eval_gegenbauer_double__t(double, double, double) noexcept nogil
cdef _proto_eval_gegenbauer_double__t *_proto_eval_gegenbauer_double__t_var = &_func_eval_gegenbauer[double]
from .orthogonal_eval cimport eval_gegenbauer_l as _func_eval_gegenbauer_l
ctypedef double _proto_eval_gegenbauer_l_t(long, double, double) noexcept nogil
cdef _proto_eval_gegenbauer_l_t *_proto_eval_gegenbauer_l_t_var = &_func_eval_gegenbauer_l
from .orthogonal_eval cimport eval_genlaguerre as _func_eval_genlaguerre
ctypedef double complex _proto_eval_genlaguerre_double_complex__t(double, double, double complex) noexcept nogil
cdef _proto_eval_genlaguerre_double_complex__t *_proto_eval_genlaguerre_double_complex__t_var = &_func_eval_genlaguerre[double_complex]
from .orthogonal_eval cimport eval_genlaguerre as _func_eval_genlaguerre
ctypedef double _proto_eval_genlaguerre_double__t(double, double, double) noexcept nogil
cdef _proto_eval_genlaguerre_double__t *_proto_eval_genlaguerre_double__t_var = &_func_eval_genlaguerre[double]
from .orthogonal_eval cimport eval_genlaguerre_l as _func_eval_genlaguerre_l
ctypedef double _proto_eval_genlaguerre_l_t(long, double, double) noexcept nogil
cdef _proto_eval_genlaguerre_l_t *_proto_eval_genlaguerre_l_t_var = &_func_eval_genlaguerre_l
from .orthogonal_eval cimport eval_hermite as _func_eval_hermite
ctypedef double _proto_eval_hermite_t(long, double) noexcept nogil
cdef _proto_eval_hermite_t *_proto_eval_hermite_t_var = &_func_eval_hermite
from .orthogonal_eval cimport eval_hermitenorm as _func_eval_hermitenorm
ctypedef double _proto_eval_hermitenorm_t(long, double) noexcept nogil
cdef _proto_eval_hermitenorm_t *_proto_eval_hermitenorm_t_var = &_func_eval_hermitenorm
from .orthogonal_eval cimport eval_jacobi as _func_eval_jacobi
ctypedef double complex _proto_eval_jacobi_double_complex__t(double, double, double, double complex) noexcept nogil
cdef _proto_eval_jacobi_double_complex__t *_proto_eval_jacobi_double_complex__t_var = &_func_eval_jacobi[double_complex]
from .orthogonal_eval cimport eval_jacobi as _func_eval_jacobi
ctypedef double _proto_eval_jacobi_double__t(double, double, double, double) noexcept nogil
cdef _proto_eval_jacobi_double__t *_proto_eval_jacobi_double__t_var = &_func_eval_jacobi[double]
from .orthogonal_eval cimport eval_jacobi_l as _func_eval_jacobi_l
ctypedef double _proto_eval_jacobi_l_t(long, double, double, double) noexcept nogil
cdef _proto_eval_jacobi_l_t *_proto_eval_jacobi_l_t_var = &_func_eval_jacobi_l
from .orthogonal_eval cimport eval_laguerre as _func_eval_laguerre
ctypedef double complex _proto_eval_laguerre_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_laguerre_double_complex__t *_proto_eval_laguerre_double_complex__t_var = &_func_eval_laguerre[double_complex]
from .orthogonal_eval cimport eval_laguerre as _func_eval_laguerre
ctypedef double _proto_eval_laguerre_double__t(double, double) noexcept nogil
cdef _proto_eval_laguerre_double__t *_proto_eval_laguerre_double__t_var = &_func_eval_laguerre[double]
from .orthogonal_eval cimport eval_laguerre_l as _func_eval_laguerre_l
ctypedef double _proto_eval_laguerre_l_t(long, double) noexcept nogil
cdef _proto_eval_laguerre_l_t *_proto_eval_laguerre_l_t_var = &_func_eval_laguerre_l
from .orthogonal_eval cimport eval_legendre as _func_eval_legendre
ctypedef double complex _proto_eval_legendre_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_legendre_double_complex__t *_proto_eval_legendre_double_complex__t_var = &_func_eval_legendre[double_complex]
from .orthogonal_eval cimport eval_legendre as _func_eval_legendre
ctypedef double _proto_eval_legendre_double__t(double, double) noexcept nogil
cdef _proto_eval_legendre_double__t *_proto_eval_legendre_double__t_var = &_func_eval_legendre[double]
from .orthogonal_eval cimport eval_legendre_l as _func_eval_legendre_l
ctypedef double _proto_eval_legendre_l_t(long, double) noexcept nogil
cdef _proto_eval_legendre_l_t *_proto_eval_legendre_l_t_var = &_func_eval_legendre_l
from .orthogonal_eval cimport eval_sh_chebyt as _func_eval_sh_chebyt
ctypedef double complex _proto_eval_sh_chebyt_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_sh_chebyt_double_complex__t *_proto_eval_sh_chebyt_double_complex__t_var = &_func_eval_sh_chebyt[double_complex]
from .orthogonal_eval cimport eval_sh_chebyt as _func_eval_sh_chebyt
ctypedef double _proto_eval_sh_chebyt_double__t(double, double) noexcept nogil
cdef _proto_eval_sh_chebyt_double__t *_proto_eval_sh_chebyt_double__t_var = &_func_eval_sh_chebyt[double]
from .orthogonal_eval cimport eval_sh_chebyt_l as _func_eval_sh_chebyt_l
ctypedef double _proto_eval_sh_chebyt_l_t(long, double) noexcept nogil
cdef _proto_eval_sh_chebyt_l_t *_proto_eval_sh_chebyt_l_t_var = &_func_eval_sh_chebyt_l
from .orthogonal_eval cimport eval_sh_chebyu as _func_eval_sh_chebyu
ctypedef double complex _proto_eval_sh_chebyu_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_sh_chebyu_double_complex__t *_proto_eval_sh_chebyu_double_complex__t_var = &_func_eval_sh_chebyu[double_complex]
from .orthogonal_eval cimport eval_sh_chebyu as _func_eval_sh_chebyu
ctypedef double _proto_eval_sh_chebyu_double__t(double, double) noexcept nogil
cdef _proto_eval_sh_chebyu_double__t *_proto_eval_sh_chebyu_double__t_var = &_func_eval_sh_chebyu[double]
from .orthogonal_eval cimport eval_sh_chebyu_l as _func_eval_sh_chebyu_l
ctypedef double _proto_eval_sh_chebyu_l_t(long, double) noexcept nogil
cdef _proto_eval_sh_chebyu_l_t *_proto_eval_sh_chebyu_l_t_var = &_func_eval_sh_chebyu_l
from .orthogonal_eval cimport eval_sh_jacobi as _func_eval_sh_jacobi
ctypedef double complex _proto_eval_sh_jacobi_double_complex__t(double, double, double, double complex) noexcept nogil
cdef _proto_eval_sh_jacobi_double_complex__t *_proto_eval_sh_jacobi_double_complex__t_var = &_func_eval_sh_jacobi[double_complex]
from .orthogonal_eval cimport eval_sh_jacobi as _func_eval_sh_jacobi
ctypedef double _proto_eval_sh_jacobi_double__t(double, double, double, double) noexcept nogil
cdef _proto_eval_sh_jacobi_double__t *_proto_eval_sh_jacobi_double__t_var = &_func_eval_sh_jacobi[double]
from .orthogonal_eval cimport eval_sh_jacobi_l as _func_eval_sh_jacobi_l
ctypedef double _proto_eval_sh_jacobi_l_t(long, double, double, double) noexcept nogil
cdef _proto_eval_sh_jacobi_l_t *_proto_eval_sh_jacobi_l_t_var = &_func_eval_sh_jacobi_l
from .orthogonal_eval cimport eval_sh_legendre as _func_eval_sh_legendre
ctypedef double complex _proto_eval_sh_legendre_double_complex__t(double, double complex) noexcept nogil
cdef _proto_eval_sh_legendre_double_complex__t *_proto_eval_sh_legendre_double_complex__t_var = &_func_eval_sh_legendre[double_complex]
from .orthogonal_eval cimport eval_sh_legendre as _func_eval_sh_legendre
ctypedef double _proto_eval_sh_legendre_double__t(double, double) noexcept nogil
cdef _proto_eval_sh_legendre_double__t *_proto_eval_sh_legendre_double__t_var = &_func_eval_sh_legendre[double]
from .orthogonal_eval cimport eval_sh_legendre_l as _func_eval_sh_legendre_l
ctypedef double _proto_eval_sh_legendre_l_t(long, double) noexcept nogil
cdef _proto_eval_sh_legendre_l_t *_proto_eval_sh_legendre_l_t_var = &_func_eval_sh_legendre_l
from ._cunity cimport cexpm1 as _func_cexpm1
ctypedef double complex _proto_cexpm1_t(double complex) noexcept nogil
cdef _proto_cexpm1_t *_proto_cexpm1_t_var = &_func_cexpm1
from ._legacy cimport expn_unsafe as _func_expn_unsafe
ctypedef double _proto_expn_unsafe_t(double, double) noexcept nogil
cdef _proto_expn_unsafe_t *_proto_expn_unsafe_t_var = &_func_expn_unsafe
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_expn "expn"(npy_int, npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_fdtr "fdtr"(npy_double, npy_double, npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_fdtrc "fdtrc"(npy_double, npy_double, npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_fdtri "fdtri"(npy_double, npy_double, npy_double)nogil
from ._cdflib_wrappers cimport fdtridfd as _func_fdtridfd
ctypedef double _proto_fdtridfd_t(double, double, double) noexcept nogil
cdef _proto_fdtridfd_t *_proto_fdtridfd_t_var = &_func_fdtridfd
cdef extern from r"_ufuncs_defs.h":
    cdef npy_int _func_cephes_fresnl_wrap "cephes_fresnl_wrap"(npy_double, npy_double *, npy_double *)nogil
from ._cdflib_wrappers cimport gdtria as _func_gdtria
ctypedef double _proto_gdtria_t(double, double, double) noexcept nogil
cdef _proto_gdtria_t *_proto_gdtria_t_var = &_func_gdtria
from ._cdflib_wrappers cimport gdtrib as _func_gdtrib
ctypedef double _proto_gdtrib_t(double, double, double) noexcept nogil
cdef _proto_gdtrib_t *_proto_gdtrib_t_var = &_func_gdtrib
from ._cdflib_wrappers cimport gdtrix as _func_gdtrix
ctypedef double _proto_gdtrix_t(double, double, double) noexcept nogil
cdef _proto_gdtrix_t *_proto_gdtrix_t_var = &_func_gdtrix
from ._convex_analysis cimport huber as _func_huber
ctypedef double _proto_huber_t(double, double) noexcept nogil
cdef _proto_huber_t *_proto_huber_t_var = &_func_huber
from ._hyp0f1 cimport _hyp0f1_cmplx as _func__hyp0f1_cmplx
ctypedef double complex _proto__hyp0f1_cmplx_t(double, double complex) noexcept nogil
cdef _proto__hyp0f1_cmplx_t *_proto__hyp0f1_cmplx_t_var = &_func__hyp0f1_cmplx
from ._hyp0f1 cimport _hyp0f1_real as _func__hyp0f1_real
ctypedef double _proto__hyp0f1_real_t(double, double) noexcept nogil
cdef _proto__hyp0f1_real_t *_proto__hyp0f1_real_t_var = &_func__hyp0f1_real
cdef extern from r"_ufuncs_defs.h":
    cdef npy_cdouble _func_chyp1f1_wrap "chyp1f1_wrap"(npy_double, npy_double, npy_cdouble)nogil
from ._hypergeometric cimport hyperu as _func_hyperu
ctypedef double _proto_hyperu_t(double, double, double) noexcept nogil
cdef _proto_hyperu_t *_proto_hyperu_t_var = &_func_hyperu
from ._boxcox cimport inv_boxcox as _func_inv_boxcox
ctypedef double _proto_inv_boxcox_t(double, double) noexcept nogil
cdef _proto_inv_boxcox_t *_proto_inv_boxcox_t_var = &_func_inv_boxcox
from ._boxcox cimport inv_boxcox1p as _func_inv_boxcox1p
ctypedef double _proto_inv_boxcox1p_t(double, double) noexcept nogil
cdef _proto_inv_boxcox1p_t *_proto_inv_boxcox1p_t_var = &_func_inv_boxcox1p
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_j0 "j0"(npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_j1 "j1"(npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_k0 "k0"(npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_k0e "k0e"(npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_k1 "k1"(npy_double)nogil
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_k1e "k1e"(npy_double)nogil
from ._convex_analysis cimport kl_div as _func_kl_div
ctypedef double _proto_kl_div_t(double, double) noexcept nogil
cdef _proto_kl_div_t *_proto_kl_div_t_var = &_func_kl_div
from ._legacy cimport kn_unsafe as _func_kn_unsafe
ctypedef double _proto_kn_unsafe_t(double, double) noexcept nogil
cdef _proto_kn_unsafe_t *_proto_kn_unsafe_t_var = &_func_kn_unsafe
from ._cunity cimport clog1p as _func_clog1p
ctypedef double complex _proto_clog1p_t(double complex) noexcept nogil
cdef _proto_clog1p_t *_proto_clog1p_t_var = &_func_clog1p
cdef extern from r"_ufuncs_defs.h":
    cdef npy_double _func_pmv_wrap "pmv_wrap"(npy_double, npy_double, npy_double)nogil
from ._legacy cimport nbdtr_unsafe as _func_nbdtr_unsafe
ctypedef double _proto_nbdtr_unsafe_t(double, double, double) noexcept nogil
cdef _proto_nbdtr_unsafe_t *_proto_nbdtr_unsafe_t_var = &_func_nbdtr_unsafe
from ._legacy cimport nbdtrc_unsafe as _func_nbdtrc_unsafe
ctypedef double _proto_nbdtrc_unsafe_t(double, double, double) noexcept nogil
cdef _proto_nbdtrc_unsafe_t *_proto_nbdtrc_unsafe_t_var = &_func_nbdtrc_unsafe
from ._legacy cimport nbdtri_unsafe as _func_nbdtri_unsafe
ctypedef double _proto_nbdtri_unsafe_t(double, double, double) noexcept nogil
cdef _proto_nbdtri_unsafe_t *_proto_nbdtri_unsafe_t_var = &_func_nbdtri_unsafe
from ._cdflib_wrappers cimport nbdtrik as _func_nbdtrik
ctypedef double _proto_nbdtrik_t(double, double, double) noexcept nogil
cdef _proto_nbdtrik_t *_proto_nbdtrik_t_var = &_func_nbdtrik
from ._cdflib_wrappers cimport nbdtrin as _func_nbdtrin
ctypedef double _proto_nbdtrin_t(double, double, double) noexcept nogil
cdef _proto_nbdtrin_t *_proto_nbdtrin_t_var = &_func_nbdtrin
from ._cdflib_wrappers cimport ncfdtr as _func_ncfdtr
ctypedef double _proto_ncfdtr_t(double, double, double, double) noexcept nogil
cdef _proto_ncfdtr_t *_proto_ncfdtr_t_var = &_func_ncfdtr
from ._cdflib_wrappers cimport ncfdtri as _func_ncfdtri
ctypedef double _proto_ncfdtri_t(double, double, double, double) noexcept nogil
cdef _proto_ncfdtri_t *_proto_ncfdtri_t_var = &_func_ncfdtri
from ._cdflib_wrappers cimport ncfdtridfd as _func_ncfdtridfd
ctypedef double _proto_ncfdtridfd_t(double, double, double, double) noexcept nogil
cdef _proto_ncfdtridfd_t *_proto_ncfdtridfd_t_var = &_func_ncfdtridfd
from ._cdflib_wrappers cimport ncfdtridfn as _func_ncfdtridfn
ctypedef double _proto_ncfdtridfn_t(double, double, double, double) noexcept nogil
cdef _proto_ncfdtridfn_t *_proto_ncfdtridfn_t_var = &_func_ncfdtridfn
from ._cdflib_wrappers cimport ncfdtrinc as _func_ncfdtrinc
ctypedef double _proto_ncfdtrinc_t(double, double, double, double) noexcept nogil
cdef _proto_ncfdtrinc_t *_proto_ncfdtrinc_t_var = &_func_ncfdtrinc
from ._cdflib_wrappers cimport nctdtr as _func_nctdtr
ctypedef double _proto_nctdtr_t(double, double, double) noexcept nogil
cdef _proto_nctdtr_t *_proto_nctdtr_t_var = &_func_nctdtr
from ._cdflib_wrappers cimport nctdtridf as _func_nctdtridf
ctypedef double _proto_nctdtridf_t(double, double, double) noexcept nogil
cdef _proto_nctdtridf_t *_proto_nctdtridf_t_var = &_func_nctdtridf
from ._cdflib_wrappers cimport nctdtrinc as _func_nctdtrinc
ctypedef double _proto_nctdtrinc_t(double, double, double) noexcept nogil
cdef _proto_nctdtrinc_t *_proto_nctdtrinc_t_var = &_func_nctdtrinc
from ._cdflib_wrappers cimport nctdtrit as _func_nctdtrit
ctypedef double _proto_nctdtrit_t(double, double, double) noexcept nogil
cdef _proto_nctdtrit_t *_proto_nctdtrit_t_var = &_func_nctdtrit
from ._cdflib_wrappers cimport nrdtrimn as _func_nrdtrimn
ctypedef double _proto_nrdtrimn_t(double, double, double) noexcept nogil
cdef _proto_nrdtrimn_t *_proto_nrdtrimn_t_var = &_func_nrdtrimn
from ._cdflib_wrappers cimport nrdtrisd as _func_nrdtrisd
ctypedef double _proto_nrdtrisd_t(double, double, double) noexcept nogil
cdef _proto_nrdtrisd_t *_proto_nrdtrisd_t_var = &_func_nrdtrisd
from ._legacy cimport pdtri_unsafe as _func_pdtri_unsafe
ctypedef double _proto_pdtri_unsafe_t(double, double) noexcept nogil
cdef _proto_pdtri_unsafe_t *_proto_pdtri_unsafe_t_var = &_func_pdtri_unsafe
from ._cdflib_wrappers cimport pdtrik as _func_pdtrik
ctypedef double _proto_pdtrik_t(double, double) noexcept nogil
cdef _proto_pdtrik_t *_proto_pdtrik_t_var = &_func_pdtrik
from ._convex_analysis cimport pseudo_huber as _func_pseudo_huber
ctypedef double _proto_pseudo_huber_t(double, double) noexcept nogil
cdef _proto_pseudo_huber_t *_proto_pseudo_huber_t_var = &_func_pseudo_huber
from ._convex_analysis cimport rel_entr as _func_rel_entr
ctypedef double _proto_rel_entr_t(double, double) noexcept nogil
cdef _proto_rel_entr_t *_proto_rel_entr_t_var = &_func_rel_entr
from ._sici cimport cshichi as _func_cshichi
ctypedef int _proto_cshichi_t(double complex, double complex *, double complex *) noexcept nogil
cdef _proto_cshichi_t *_proto_cshichi_t_var = &_func_cshichi
cdef extern from r"_ufuncs_defs.h":
    cdef npy_int _func_cephes_shichi_wrap "cephes_shichi_wrap"(npy_double, npy_double *, npy_double *)nogil
from ._sici cimport csici as _func_csici
ctypedef int _proto_csici_t(double complex, double complex *, double complex *) noexcept nogil
cdef _proto_csici_t *_proto_csici_t_var = &_func_csici
cdef extern from r"_ufuncs_defs.h":
    cdef npy_int _func_cephes_sici_wrap "cephes_sici_wrap"(npy_double, npy_double *, npy_double *)nogil
from ._legacy cimport smirnov_unsafe as _func_smirnov_unsafe
ctypedef double _proto_smirnov_unsafe_t(double, double) noexcept nogil
cdef _proto_smirnov_unsafe_t *_proto_smirnov_unsafe_t_var = &_func_smirnov_unsafe
from ._legacy cimport smirnovi_unsafe as _func_smirnovi_unsafe
ctypedef double _proto_smirnovi_unsafe_t(double, double) noexcept nogil
cdef _proto_smirnovi_unsafe_t *_proto_smirnovi_unsafe_t_var = &_func_smirnovi_unsafe
from ._spence cimport cspence as _func_cspence
ctypedef double complex _proto_cspence_t(double complex) noexcept nogil
cdef _proto_cspence_t *_proto_cspence_t_var = &_func_cspence
from ._cdflib_wrappers cimport stdtr as _func_stdtr
ctypedef double _proto_stdtr_t(double, double) noexcept nogil
cdef _proto_stdtr_t *_proto_stdtr_t_var = &_func_stdtr
from ._cdflib_wrappers cimport stdtridf as _func_stdtridf
ctypedef double _proto_stdtridf_t(double, double) noexcept nogil
cdef _proto_stdtridf_t *_proto_stdtridf_t_var = &_func_stdtridf
from ._cdflib_wrappers cimport stdtrit as _func_stdtrit
ctypedef double _proto_stdtrit_t(double, double) noexcept nogil
cdef _proto_stdtrit_t *_proto_stdtrit_t_var = &_func_stdtrit
from ._xlogy cimport xlog1py as _func_xlog1py
ctypedef double _proto_xlog1py_double__t(double, double) noexcept nogil
cdef _proto_xlog1py_double__t *_proto_xlog1py_double__t_var = &_func_xlog1py[double]
from ._xlogy cimport xlog1py as _func_xlog1py
ctypedef double complex _proto_xlog1py_double_complex__t(double complex, double complex) noexcept nogil
cdef _proto_xlog1py_double_complex__t *_proto_xlog1py_double_complex__t_var = &_func_xlog1py[double_complex]
from ._xlogy cimport xlogy as _func_xlogy
ctypedef double _proto_xlogy_double__t(double, double) noexcept nogil
cdef _proto_xlogy_double__t *_proto_xlogy_double__t_var = &_func_xlogy[double]
from ._xlogy cimport xlogy as _func_xlogy
ctypedef double complex _proto_xlogy_double_complex__t(double complex, double complex) noexcept nogil
cdef _proto_xlogy_double_complex__t *_proto_xlogy_double_complex__t_var = &_func_xlogy[double_complex]
from ._legacy cimport yn_unsafe as _func_yn_unsafe
ctypedef double _proto_yn_unsafe_t(double, double) noexcept nogil
cdef _proto_yn_unsafe_t *_proto_yn_unsafe_t_var = &_func_yn_unsafe
from ._ndtri_exp cimport ndtri_exp as _func_ndtri_exp
ctypedef double _proto_ndtri_exp_t(double) noexcept nogil
cdef _proto_ndtri_exp_t *_proto_ndtri_exp_t_var = &_func_ndtri_exp

cpdef double voigt_profile(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.voigt_profile"""
    return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_voigt_profile)(x0, x1, x2)

cpdef double agm(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.agm"""
    return _func_agm(x0, x1)

cdef void airy(Dd_number_t x0, Dd_number_t *y0, Dd_number_t *y1, Dd_number_t *y2, Dd_number_t *y3) noexcept nogil:
    """See the documentation for scipy.special.airy"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    cdef npy_cdouble tmp2
    cdef npy_cdouble tmp3
    if Dd_number_t is double:
        special_airy(x0, y0, y1, y2, y3)
    elif Dd_number_t is double_complex:
        special_cairy(_complexstuff.npy_cdouble_from_double_complex(x0), &tmp0, &tmp1, &tmp2, &tmp3)
        y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
        y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)
        y2[0] = _complexstuff.double_complex_from_npy_cdouble(tmp2)
        y3[0] = _complexstuff.double_complex_from_npy_cdouble(tmp3)
    else:
        if Dd_number_t is double_complex:
            y0[0] = NAN
            y1[0] = NAN
            y2[0] = NAN
            y3[0] = NAN
        else:
            y0[0] = NAN
            y1[0] = NAN
            y2[0] = NAN
            y3[0] = NAN

def _airy_pywrap(Dd_number_t x0):
    cdef Dd_number_t y0
    cdef Dd_number_t y1
    cdef Dd_number_t y2
    cdef Dd_number_t y3
    airy(x0, &y0, &y1, &y2, &y3)
    return y0, y1, y2, y3

cdef void airye(Dd_number_t x0, Dd_number_t *y0, Dd_number_t *y1, Dd_number_t *y2, Dd_number_t *y3) noexcept nogil:
    """See the documentation for scipy.special.airye"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    cdef npy_cdouble tmp2
    cdef npy_cdouble tmp3
    if Dd_number_t is double_complex:
        special_cairye(_complexstuff.npy_cdouble_from_double_complex(x0), &tmp0, &tmp1, &tmp2, &tmp3)
        y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
        y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)
        y2[0] = _complexstuff.double_complex_from_npy_cdouble(tmp2)
        y3[0] = _complexstuff.double_complex_from_npy_cdouble(tmp3)
    elif Dd_number_t is double:
        special_airye(x0, y0, y1, y2, y3)
    else:
        if Dd_number_t is double_complex:
            y0[0] = NAN
            y1[0] = NAN
            y2[0] = NAN
            y3[0] = NAN
        else:
            y0[0] = NAN
            y1[0] = NAN
            y2[0] = NAN
            y3[0] = NAN

def _airye_pywrap(Dd_number_t x0):
    cdef Dd_number_t y0
    cdef Dd_number_t y1
    cdef Dd_number_t y2
    cdef Dd_number_t y3
    airye(x0, &y0, &y1, &y2, &y3)
    return y0, y1, y2, y3

cpdef double bdtr(double x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.bdtr"""
    if dl_number_t is double:
        return _func_bdtr_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_bdtr(x0, x1, x2)
    else:
        return NAN

cpdef double bdtrc(double x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.bdtrc"""
    if dl_number_t is double:
        return _func_bdtrc_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_bdtrc(x0, x1, x2)
    else:
        return NAN

cpdef double bdtri(double x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.bdtri"""
    if dl_number_t is double:
        return _func_bdtri_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_bdtri(x0, x1, x2)
    else:
        return NAN

cpdef double bdtrik(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.bdtrik"""
    return _func_bdtrik(x0, x1, x2)

cpdef double bdtrin(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.bdtrin"""
    return _func_bdtrin(x0, x1, x2)

cpdef double bei(double x0) noexcept nogil:
    """See the documentation for scipy.special.bei"""
    return special_bei(x0)

cpdef double beip(double x0) noexcept nogil:
    """See the documentation for scipy.special.beip"""
    return special_beip(x0)

cpdef double ber(double x0) noexcept nogil:
    """See the documentation for scipy.special.ber"""
    return special_ber(x0)

cpdef double berp(double x0) noexcept nogil:
    """See the documentation for scipy.special.berp"""
    return special_berp(x0)

cpdef double besselpoly(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.besselpoly"""
    return cephes_besselpoly(x0, x1, x2)

cpdef double beta(double x0, double x1) noexcept nogil:
    return cephes_beta(x0, x1)

cpdef df_number_t betainc(df_number_t x0, df_number_t x1, df_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.betainc"""
    if df_number_t is float:
        return (<float(*)(float, float, float) noexcept nogil>scipy.special._ufuncs_cxx._export_ibeta_float)(x0, x1, x2)
    elif df_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_ibeta_double)(x0, x1, x2)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cpdef df_number_t betaincc(df_number_t x0, df_number_t x1, df_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.betaincc"""
    if df_number_t is float:
        return (<float(*)(float, float, float) noexcept nogil>scipy.special._ufuncs_cxx._export_ibetac_float)(x0, x1, x2)
    elif df_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_ibetac_double)(x0, x1, x2)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cpdef df_number_t betaincinv(df_number_t x0, df_number_t x1, df_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.betaincinv"""
    if df_number_t is float:
        return (<float(*)(float, float, float) noexcept nogil>scipy.special._ufuncs_cxx._export_ibeta_inv_float)(x0, x1, x2)
    elif df_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_ibeta_inv_double)(x0, x1, x2)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cpdef df_number_t betainccinv(df_number_t x0, df_number_t x1, df_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.betainccinv"""
    if df_number_t is float:
        return (<float(*)(float, float, float) noexcept nogil>scipy.special._ufuncs_cxx._export_ibetac_inv_float)(x0, x1, x2)
    elif df_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_ibetac_inv_double)(x0, x1, x2)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cpdef double betaln(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.betaln"""
    return cephes_lbeta(x0, x1)

cpdef double binom(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.binom"""
    return special_binom(x0, x1)

cpdef double boxcox(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.boxcox"""
    return _func_boxcox(x0, x1)

cpdef double boxcox1p(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.boxcox1p"""
    return _func_boxcox1p(x0, x1)

cpdef double btdtr(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.btdtr"""
    return cephes_btdtr(x0, x1, x2)

cpdef double btdtri(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.btdtri"""
    return cephes_btdtri(x0, x1, x2)

cpdef double btdtria(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.btdtria"""
    return _func_btdtria(x0, x1, x2)

cpdef double btdtrib(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.btdtrib"""
    return _func_btdtrib(x0, x1, x2)

cpdef double cbrt(double x0) noexcept nogil:
    """See the documentation for scipy.special.cbrt"""
    return cephes_cbrt(x0)

cpdef double chdtr(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.chdtr"""
    return cephes_chdtr(x0, x1)

cpdef double chdtrc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.chdtrc"""
    return cephes_chdtrc(x0, x1)

cpdef double chdtri(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.chdtri"""
    return cephes_chdtri(x0, x1)

cpdef double chdtriv(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.chdtriv"""
    return _func_chdtriv(x0, x1)

cpdef double chndtr(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.chndtr"""
    return _func_chndtr(x0, x1, x2)

cpdef double chndtridf(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.chndtridf"""
    return _func_chndtridf(x0, x1, x2)

cpdef double chndtrinc(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.chndtrinc"""
    return _func_chndtrinc(x0, x1, x2)

cpdef double chndtrix(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.chndtrix"""
    return _func_chndtrix(x0, x1, x2)

cpdef double cosdg(double x0) noexcept nogil:
    """See the documentation for scipy.special.cosdg"""
    return cephes_cosdg(x0)

cpdef double cosm1(double x0) noexcept nogil:
    """See the documentation for scipy.special.cosm1"""
    return cephes_cosm1(x0)

cpdef double cotdg(double x0) noexcept nogil:
    """See the documentation for scipy.special.cotdg"""
    return cephes_cotdg(x0)

cpdef Dd_number_t dawsn(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.dawsn"""
    if Dd_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_dawsn)(x0)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_dawsn_complex)(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double ellipe(double x0) noexcept nogil:
    """See the documentation for scipy.special.ellipe"""
    return cephes_ellpe(x0)

cpdef double ellipeinc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.ellipeinc"""
    return cephes_ellie(x0, x1)

cdef void ellipj(double x0, double x1, double *y0, double *y1, double *y2, double *y3) noexcept nogil:
    """See the documentation for scipy.special.ellipj"""
    _func_cephes_ellpj_wrap(x0, x1, y0, y1, y2, y3)

def _ellipj_pywrap(double x0, double x1):
    cdef double y0
    cdef double y1
    cdef double y2
    cdef double y3
    ellipj(x0, x1, &y0, &y1, &y2, &y3)
    return y0, y1, y2, y3

cpdef double ellipkinc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.ellipkinc"""
    return cephes_ellik(x0, x1)

cpdef double ellipkm1(double x0) noexcept nogil:
    """See the documentation for scipy.special.ellipkm1"""
    return cephes_ellpk(x0)

cpdef double ellipk(double x0) noexcept nogil:
    """See the documentation for scipy.special.ellipk"""
    return special_ellipk(x0)

cpdef Dd_number_t elliprc(Dd_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.elliprc"""
    if Dd_number_t is double:
        return (<double(*)(double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_fellint_RC)(x0, x1)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex, double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_cellint_RC)(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t elliprd(Dd_number_t x0, Dd_number_t x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.elliprd"""
    if Dd_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_fellint_RD)(x0, x1, x2)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex, double complex, double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_cellint_RD)(x0, x1, x2)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t elliprf(Dd_number_t x0, Dd_number_t x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.elliprf"""
    if Dd_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_fellint_RF)(x0, x1, x2)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex, double complex, double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_cellint_RF)(x0, x1, x2)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t elliprg(Dd_number_t x0, Dd_number_t x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.elliprg"""
    if Dd_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_fellint_RG)(x0, x1, x2)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex, double complex, double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_cellint_RG)(x0, x1, x2)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t elliprj(Dd_number_t x0, Dd_number_t x1, Dd_number_t x2, Dd_number_t x3) noexcept nogil:
    """See the documentation for scipy.special.elliprj"""
    if Dd_number_t is double:
        return (<double(*)(double, double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_fellint_RJ)(x0, x1, x2, x3)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex, double complex, double complex, double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_cellint_RJ)(x0, x1, x2, x3)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double entr(double x0) noexcept nogil:
    """See the documentation for scipy.special.entr"""
    return _func_entr(x0)

cpdef Dd_number_t erf(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.erf"""
    if Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erf)(x0)
    elif Dd_number_t is double:
        return cephes_erf(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t erfc(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.erfc"""
    if Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erfc_complex)(x0)
    elif Dd_number_t is double:
        return cephes_erfc(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t erfcx(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.erfcx"""
    if Dd_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erfcx)(x0)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erfcx_complex)(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t erfi(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.erfi"""
    if Dd_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erfi)(x0)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_erfi_complex)(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef df_number_t erfinv(df_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.erfinv"""
    if df_number_t is float:
        return (<float(*)(float) noexcept nogil>scipy.special._ufuncs_cxx._export_erfinv_float)(x0)
    elif df_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_erfinv_double)(x0)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cpdef double erfcinv(double x0) noexcept nogil:
    """See the documentation for scipy.special.erfcinv"""
    return cephes_erfcinv(x0)

cpdef Dd_number_t eval_chebyc(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_chebyc"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_chebyc[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_chebyc[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_chebyc_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_chebys(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_chebys"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_chebys[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_chebys[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_chebys_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_chebyt(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_chebyt"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_chebyt[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_chebyt[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_chebyt_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_chebyu(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_chebyu"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_chebyu[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_chebyu[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_chebyu_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_gegenbauer(dl_number_t x0, double x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.eval_gegenbauer"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_gegenbauer[double_complex](x0, x1, x2)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_gegenbauer[double](x0, x1, x2)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_gegenbauer_l(x0, x1, x2)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_genlaguerre(dl_number_t x0, double x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.eval_genlaguerre"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_genlaguerre[double_complex](x0, x1, x2)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_genlaguerre[double](x0, x1, x2)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_genlaguerre_l(x0, x1, x2)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double eval_hermite(long x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.eval_hermite"""
    return _func_eval_hermite(x0, x1)

cpdef double eval_hermitenorm(long x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.eval_hermitenorm"""
    return _func_eval_hermitenorm(x0, x1)

cpdef Dd_number_t eval_jacobi(dl_number_t x0, double x1, double x2, Dd_number_t x3) noexcept nogil:
    """See the documentation for scipy.special.eval_jacobi"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_jacobi[double_complex](x0, x1, x2, x3)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_jacobi[double](x0, x1, x2, x3)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_jacobi_l(x0, x1, x2, x3)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_laguerre(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_laguerre"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_laguerre[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_laguerre[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_laguerre_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_legendre(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_legendre"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_legendre[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_legendre[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_legendre_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_sh_chebyt(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_sh_chebyt"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_sh_chebyt[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_sh_chebyt[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_sh_chebyt_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_sh_chebyu(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_sh_chebyu"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_sh_chebyu[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_sh_chebyu[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_sh_chebyu_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_sh_jacobi(dl_number_t x0, double x1, double x2, Dd_number_t x3) noexcept nogil:
    """See the documentation for scipy.special.eval_sh_jacobi"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_sh_jacobi[double_complex](x0, x1, x2, x3)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_sh_jacobi[double](x0, x1, x2, x3)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_sh_jacobi_l(x0, x1, x2, x3)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t eval_sh_legendre(dl_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.eval_sh_legendre"""
    if dl_number_t is double and Dd_number_t is double_complex:
        return _func_eval_sh_legendre[double_complex](x0, x1)
    elif dl_number_t is double and Dd_number_t is double:
        return _func_eval_sh_legendre[double](x0, x1)
    elif dl_number_t is long and Dd_number_t is double:
        return _func_eval_sh_legendre_l(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t exp1(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.exp1"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_cexp1(_complexstuff.npy_cdouble_from_double_complex(x0)))
    elif Dd_number_t is double:
        return special_exp1(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double exp10(double x0) noexcept nogil:
    """See the documentation for scipy.special.exp10"""
    return cephes_exp10(x0)

cpdef double exp2(double x0) noexcept nogil:
    """See the documentation for scipy.special.exp2"""
    return cephes_exp2(x0)

cpdef Dd_number_t expi(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.expi"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_cexpi(_complexstuff.npy_cdouble_from_double_complex(x0)))
    elif Dd_number_t is double:
        return special_expi(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef dfg_number_t expit(dfg_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.expit"""
    if dfg_number_t is double:
        return special_expit(x0)
    elif dfg_number_t is float:
        return special_expitf(x0)
    elif dfg_number_t is long_double:
        return special_expitl(x0)
    else:
        if dfg_number_t is double:
            return NAN
        elif dfg_number_t is float:
            return NAN
        else:
            return NAN

cpdef Dd_number_t expm1(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.expm1"""
    if Dd_number_t is double_complex:
        return _func_cexpm1(x0)
    elif Dd_number_t is double:
        return cephes_expm1(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double expn(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.expn"""
    if dl_number_t is double:
        return _func_expn_unsafe(x0, x1)
    elif dl_number_t is long:
        return cephes_expn(x0, x1)
    else:
        return NAN

cpdef double exprel(double x0) noexcept nogil:
    """See the documentation for scipy.special.exprel"""
    return special_exprel(x0)

cpdef double fdtr(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.fdtr"""
    return cephes_fdtr(x0, x1, x2)

cpdef double fdtrc(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.fdtrc"""
    return cephes_fdtrc(x0, x1, x2)

cpdef double fdtri(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.fdtri"""
    return cephes_fdtri(x0, x1, x2)

cpdef double fdtridfd(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.fdtridfd"""
    return _func_fdtridfd(x0, x1, x2)

cdef void fresnel(Dd_number_t x0, Dd_number_t *y0, Dd_number_t *y1) noexcept nogil:
    """See the documentation for scipy.special.fresnel"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    if Dd_number_t is double:
        _func_cephes_fresnl_wrap(x0, y0, y1)
    elif Dd_number_t is double_complex:
        _func_cfresnl_wrap(_complexstuff.npy_cdouble_from_double_complex(x0), &tmp0, &tmp1)
        y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
        y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)
    else:
        if Dd_number_t is double_complex:
            y0[0] = NAN
            y1[0] = NAN
        else:
            y0[0] = NAN
            y1[0] = NAN

def _fresnel_pywrap(Dd_number_t x0):
    cdef Dd_number_t y0
    cdef Dd_number_t y1
    fresnel(x0, &y0, &y1)
    return y0, y1

cpdef Dd_number_t gamma(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.gamma"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_cgamma(_complexstuff.npy_cdouble_from_double_complex(x0)))
    elif Dd_number_t is double:
        return special_gamma(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double gammainc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.gammainc"""
    return cephes_igam(x0, x1)

cpdef double gammaincc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.gammaincc"""
    return cephes_igamc(x0, x1)

cpdef double gammainccinv(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.gammainccinv"""
    return cephes_igamci(x0, x1)

cpdef double gammaincinv(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.gammaincinv"""
    return cephes_igami(x0, x1)

cpdef double gammaln(double x0) noexcept nogil:
    """See the documentation for scipy.special.gammaln"""
    return _func_gammaln_wrap(x0)

cpdef double gammasgn(double x0) noexcept nogil:
    """See the documentation for scipy.special.gammasgn"""
    return cephes_gammasgn(x0)

cpdef double gdtr(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.gdtr"""
    return cephes_gdtr(x0, x1, x2)

cpdef double gdtrc(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.gdtrc"""
    return cephes_gdtrc(x0, x1, x2)

cpdef double gdtria(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.gdtria"""
    return _func_gdtria(x0, x1, x2)

cpdef double gdtrib(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.gdtrib"""
    return _func_gdtrib(x0, x1, x2)

cpdef double gdtrix(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.gdtrix"""
    return _func_gdtrix(x0, x1, x2)

cpdef double complex hankel1(double x0, double complex x1) noexcept nogil:
    """See the documentation for scipy.special.hankel1"""
    return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_hankel_1(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))

cpdef double complex hankel1e(double x0, double complex x1) noexcept nogil:
    """See the documentation for scipy.special.hankel1e"""
    return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_hankel_1e(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))

cpdef double complex hankel2(double x0, double complex x1) noexcept nogil:
    """See the documentation for scipy.special.hankel2"""
    return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_hankel_2(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))

cpdef double complex hankel2e(double x0, double complex x1) noexcept nogil:
    """See the documentation for scipy.special.hankel2e"""
    return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_hankel_2e(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))

cpdef double huber(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.huber"""
    return _func_huber(x0, x1)

cpdef Dd_number_t hyp0f1(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.hyp0f1"""
    if Dd_number_t is double_complex:
        return _func__hyp0f1_cmplx(x0, x1)
    elif Dd_number_t is double:
        return _func__hyp0f1_real(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t hyp1f1(double x0, double x1, Dd_number_t x2) noexcept nogil:
    """See the documentation for scipy.special.hyp1f1"""
    if Dd_number_t is double:
        return (<double(*)(double, double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_hyp1f1_double)(x0, x1, x2)
    elif Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(_func_chyp1f1_wrap(x0, x1, _complexstuff.npy_cdouble_from_double_complex(x2)))
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t hyp2f1(double x0, double x1, double x2, Dd_number_t x3) noexcept nogil:
    """See the documentation for scipy.special.hyp2f1"""
    if Dd_number_t is double:
        return special_hyp2f1(x0, x1, x2, x3)
    elif Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_chyp2f1(x0, x1, x2, _complexstuff.npy_cdouble_from_double_complex(x3)))
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double hyperu(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.hyperu"""
    return _func_hyperu(x0, x1, x2)

cpdef double i0(double x0) noexcept nogil:
    """See the documentation for scipy.special.i0"""
    return cephes_i0(x0)

cpdef double i0e(double x0) noexcept nogil:
    """See the documentation for scipy.special.i0e"""
    return cephes_i0e(x0)

cpdef double i1(double x0) noexcept nogil:
    """See the documentation for scipy.special.i1"""
    return cephes_i1(x0)

cpdef double i1e(double x0) noexcept nogil:
    """See the documentation for scipy.special.i1e"""
    return cephes_i1e(x0)

cpdef double inv_boxcox(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.inv_boxcox"""
    return _func_inv_boxcox(x0, x1)

cpdef double inv_boxcox1p(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.inv_boxcox1p"""
    return _func_inv_boxcox1p(x0, x1)

cdef void it2i0k0(double x0, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.it2i0k0"""
    _func_it2i0k0_wrap(x0, y0, y1)

def _it2i0k0_pywrap(double x0):
    cdef double y0
    cdef double y1
    it2i0k0(x0, &y0, &y1)
    return y0, y1

cdef void it2j0y0(double x0, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.it2j0y0"""
    _func_it2j0y0_wrap(x0, y0, y1)

def _it2j0y0_pywrap(double x0):
    cdef double y0
    cdef double y1
    it2j0y0(x0, &y0, &y1)
    return y0, y1

cpdef double it2struve0(double x0) noexcept nogil:
    """See the documentation for scipy.special.it2struve0"""
    return special_it2struve0(x0)

cdef void itairy(double x0, double *y0, double *y1, double *y2, double *y3) noexcept nogil:
    """See the documentation for scipy.special.itairy"""
    special_itairy(x0, y0, y1, y2, y3)

def _itairy_pywrap(double x0):
    cdef double y0
    cdef double y1
    cdef double y2
    cdef double y3
    itairy(x0, &y0, &y1, &y2, &y3)
    return y0, y1, y2, y3

cdef void iti0k0(double x0, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.iti0k0"""
    _func_it1i0k0_wrap(x0, y0, y1)

def _iti0k0_pywrap(double x0):
    cdef double y0
    cdef double y1
    iti0k0(x0, &y0, &y1)
    return y0, y1

cdef void itj0y0(double x0, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.itj0y0"""
    _func_it1j0y0_wrap(x0, y0, y1)

def _itj0y0_pywrap(double x0):
    cdef double y0
    cdef double y1
    itj0y0(x0, &y0, &y1)
    return y0, y1

cpdef double itmodstruve0(double x0) noexcept nogil:
    """See the documentation for scipy.special.itmodstruve0"""
    return special_itmodstruve0(x0)

cpdef double itstruve0(double x0) noexcept nogil:
    """See the documentation for scipy.special.itstruve0"""
    return special_itstruve0(x0)

cpdef Dd_number_t iv(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.iv"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_i(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_i(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t ive(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.ive"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_ie(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_ie(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double j0(double x0) noexcept nogil:
    """See the documentation for scipy.special.j0"""
    return cephes_j0(x0)

cpdef double j1(double x0) noexcept nogil:
    """See the documentation for scipy.special.j1"""
    return cephes_j1(x0)

cpdef Dd_number_t jv(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.jv"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_j(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_j(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t jve(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.jve"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_je(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_je(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double k0(double x0) noexcept nogil:
    """See the documentation for scipy.special.k0"""
    return cephes_k0(x0)

cpdef double k0e(double x0) noexcept nogil:
    """See the documentation for scipy.special.k0e"""
    return cephes_k0e(x0)

cpdef double k1(double x0) noexcept nogil:
    """See the documentation for scipy.special.k1"""
    return cephes_k1(x0)

cpdef double k1e(double x0) noexcept nogil:
    """See the documentation for scipy.special.k1e"""
    return cephes_k1e(x0)

cpdef double kei(double x0) noexcept nogil:
    """See the documentation for scipy.special.kei"""
    return special_kei(x0)

cpdef double keip(double x0) noexcept nogil:
    """See the documentation for scipy.special.keip"""
    return special_keip(x0)

cdef void kelvin(double x0, double complex *y0, double complex *y1, double complex *y2, double complex *y3) noexcept nogil:
    """See the documentation for scipy.special.kelvin"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    cdef npy_cdouble tmp2
    cdef npy_cdouble tmp3
    special_ckelvin(x0, &tmp0, &tmp1, &tmp2, &tmp3)
    y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
    y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)
    y2[0] = _complexstuff.double_complex_from_npy_cdouble(tmp2)
    y3[0] = _complexstuff.double_complex_from_npy_cdouble(tmp3)

def _kelvin_pywrap(double x0):
    cdef double complex y0
    cdef double complex y1
    cdef double complex y2
    cdef double complex y3
    kelvin(x0, &y0, &y1, &y2, &y3)
    return y0, y1, y2, y3

cpdef double ker(double x0) noexcept nogil:
    """See the documentation for scipy.special.ker"""
    return special_ker(x0)

cpdef double kerp(double x0) noexcept nogil:
    """See the documentation for scipy.special.kerp"""
    return special_kerp(x0)

cpdef double kl_div(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.kl_div"""
    return _func_kl_div(x0, x1)

cpdef double kn(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.kn"""
    if dl_number_t is double:
        return _func_kn_unsafe(x0, x1)
    elif dl_number_t is long:
        return special_cyl_bessel_k_int(x0, x1)
    else:
        return NAN

cpdef double kolmogi(double x0) noexcept nogil:
    """See the documentation for scipy.special.kolmogi"""
    return cephes_kolmogi(x0)

cpdef double kolmogorov(double x0) noexcept nogil:
    """See the documentation for scipy.special.kolmogorov"""
    return cephes_kolmogorov(x0)

cpdef Dd_number_t kv(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.kv"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_k(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_k(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t kve(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.kve"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_ke(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_ke(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t log1p(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.log1p"""
    if Dd_number_t is double_complex:
        return _func_clog1p(x0)
    elif Dd_number_t is double:
        return cephes_log1p(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef dfg_number_t log_expit(dfg_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.log_expit"""
    if dfg_number_t is double:
        return special_log_expit(x0)
    elif dfg_number_t is float:
        return special_log_expitf(x0)
    elif dfg_number_t is long_double:
        return special_log_expitl(x0)
    else:
        if dfg_number_t is double:
            return NAN
        elif dfg_number_t is float:
            return NAN
        else:
            return NAN

cpdef Dd_number_t log_ndtr(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.log_ndtr"""
    if Dd_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_log_ndtr)(x0)
    elif Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_log_ndtr_complex)(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t loggamma(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.loggamma"""
    if Dd_number_t is double:
        return special_loggamma(x0)
    elif Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_cloggamma(_complexstuff.npy_cdouble_from_double_complex(x0)))
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef dfg_number_t logit(dfg_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.logit"""
    if dfg_number_t is double:
        return special_logit(x0)
    elif dfg_number_t is float:
        return special_logitf(x0)
    elif dfg_number_t is long_double:
        return special_logitl(x0)
    else:
        if dfg_number_t is double:
            return NAN
        elif dfg_number_t is float:
            return NAN
        else:
            return NAN

cpdef double lpmv(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.lpmv"""
    return _func_pmv_wrap(x0, x1, x2)

cpdef double mathieu_a(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_a"""
    return _func_cem_cva_wrap(x0, x1)

cpdef double mathieu_b(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_b"""
    return _func_sem_cva_wrap(x0, x1)

cdef void mathieu_cem(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_cem"""
    _func_cem_wrap(x0, x1, x2, y0, y1)

def _mathieu_cem_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_cem(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void mathieu_modcem1(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_modcem1"""
    _func_mcm1_wrap(x0, x1, x2, y0, y1)

def _mathieu_modcem1_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_modcem1(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void mathieu_modcem2(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_modcem2"""
    _func_mcm2_wrap(x0, x1, x2, y0, y1)

def _mathieu_modcem2_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_modcem2(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void mathieu_modsem1(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_modsem1"""
    _func_msm1_wrap(x0, x1, x2, y0, y1)

def _mathieu_modsem1_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_modsem1(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void mathieu_modsem2(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_modsem2"""
    _func_msm2_wrap(x0, x1, x2, y0, y1)

def _mathieu_modsem2_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_modsem2(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void mathieu_sem(double x0, double x1, double x2, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.mathieu_sem"""
    _func_sem_wrap(x0, x1, x2, y0, y1)

def _mathieu_sem_pywrap(double x0, double x1, double x2):
    cdef double y0
    cdef double y1
    mathieu_sem(x0, x1, x2, &y0, &y1)
    return y0, y1

cdef void modfresnelm(double x0, double complex *y0, double complex *y1) noexcept nogil:
    """See the documentation for scipy.special.modfresnelm"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    _func_modified_fresnel_minus_wrap(x0, &tmp0, &tmp1)
    y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
    y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)

def _modfresnelm_pywrap(double x0):
    cdef double complex y0
    cdef double complex y1
    modfresnelm(x0, &y0, &y1)
    return y0, y1

cdef void modfresnelp(double x0, double complex *y0, double complex *y1) noexcept nogil:
    """See the documentation for scipy.special.modfresnelp"""
    cdef npy_cdouble tmp0
    cdef npy_cdouble tmp1
    _func_modified_fresnel_plus_wrap(x0, &tmp0, &tmp1)
    y0[0] = _complexstuff.double_complex_from_npy_cdouble(tmp0)
    y1[0] = _complexstuff.double_complex_from_npy_cdouble(tmp1)

def _modfresnelp_pywrap(double x0):
    cdef double complex y0
    cdef double complex y1
    modfresnelp(x0, &y0, &y1)
    return y0, y1

cpdef double modstruve(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.modstruve"""
    return cephes_struve_l(x0, x1)

cpdef double nbdtr(dl_number_t x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nbdtr"""
    if dl_number_t is double:
        return _func_nbdtr_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_nbdtr(x0, x1, x2)
    else:
        return NAN

cpdef double nbdtrc(dl_number_t x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nbdtrc"""
    if dl_number_t is double:
        return _func_nbdtrc_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_nbdtrc(x0, x1, x2)
    else:
        return NAN

cpdef double nbdtri(dl_number_t x0, dl_number_t x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nbdtri"""
    if dl_number_t is double:
        return _func_nbdtri_unsafe(x0, x1, x2)
    elif dl_number_t is long:
        return cephes_nbdtri(x0, x1, x2)
    else:
        return NAN

cpdef double nbdtrik(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nbdtrik"""
    return _func_nbdtrik(x0, x1, x2)

cpdef double nbdtrin(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nbdtrin"""
    return _func_nbdtrin(x0, x1, x2)

cpdef double ncfdtr(double x0, double x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.ncfdtr"""
    return _func_ncfdtr(x0, x1, x2, x3)

cpdef double ncfdtri(double x0, double x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.ncfdtri"""
    return _func_ncfdtri(x0, x1, x2, x3)

cpdef double ncfdtridfd(double x0, double x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.ncfdtridfd"""
    return _func_ncfdtridfd(x0, x1, x2, x3)

cpdef double ncfdtridfn(double x0, double x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.ncfdtridfn"""
    return _func_ncfdtridfn(x0, x1, x2, x3)

cpdef double ncfdtrinc(double x0, double x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.ncfdtrinc"""
    return _func_ncfdtrinc(x0, x1, x2, x3)

cpdef double nctdtr(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nctdtr"""
    return _func_nctdtr(x0, x1, x2)

cpdef double nctdtridf(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nctdtridf"""
    return _func_nctdtridf(x0, x1, x2)

cpdef double nctdtrinc(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nctdtrinc"""
    return _func_nctdtrinc(x0, x1, x2)

cpdef double nctdtrit(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nctdtrit"""
    return _func_nctdtrit(x0, x1, x2)

cpdef Dd_number_t ndtr(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.ndtr"""
    if Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_ndtr)(x0)
    elif Dd_number_t is double:
        return cephes_ndtr(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double ndtri(double x0) noexcept nogil:
    """See the documentation for scipy.special.ndtri"""
    return cephes_ndtri(x0)

cpdef double nrdtrimn(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nrdtrimn"""
    return _func_nrdtrimn(x0, x1, x2)

cpdef double nrdtrisd(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.nrdtrisd"""
    return _func_nrdtrisd(x0, x1, x2)

cdef void obl_ang1(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_ang1"""
    y0[0] = _func_oblate_aswfa_nocv_wrap(x0, x1, x2, x3, y1)

def _obl_ang1_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    obl_ang1(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void obl_ang1_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_ang1_cv"""
    _func_oblate_aswfa_wrap(x0, x1, x2, x3, x4, y0, y1)

def _obl_ang1_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    obl_ang1_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cpdef double obl_cv(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.obl_cv"""
    return _func_oblate_segv_wrap(x0, x1, x2)

cdef void obl_rad1(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_rad1"""
    y0[0] = _func_oblate_radial1_nocv_wrap(x0, x1, x2, x3, y1)

def _obl_rad1_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    obl_rad1(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void obl_rad1_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_rad1_cv"""
    _func_oblate_radial1_wrap(x0, x1, x2, x3, x4, y0, y1)

def _obl_rad1_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    obl_rad1_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cdef void obl_rad2(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_rad2"""
    y0[0] = _func_oblate_radial2_nocv_wrap(x0, x1, x2, x3, y1)

def _obl_rad2_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    obl_rad2(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void obl_rad2_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.obl_rad2_cv"""
    _func_oblate_radial2_wrap(x0, x1, x2, x3, x4, y0, y1)

def _obl_rad2_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    obl_rad2_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cpdef double owens_t(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.owens_t"""
    return cephes_owens_t(x0, x1)

cdef void pbdv(double x0, double x1, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pbdv"""
    _func_pbdv_wrap(x0, x1, y0, y1)

def _pbdv_pywrap(double x0, double x1):
    cdef double y0
    cdef double y1
    pbdv(x0, x1, &y0, &y1)
    return y0, y1

cdef void pbvv(double x0, double x1, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pbvv"""
    _func_pbvv_wrap(x0, x1, y0, y1)

def _pbvv_pywrap(double x0, double x1):
    cdef double y0
    cdef double y1
    pbvv(x0, x1, &y0, &y1)
    return y0, y1

cdef void pbwa(double x0, double x1, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pbwa"""
    _func_pbwa_wrap(x0, x1, y0, y1)

def _pbwa_pywrap(double x0, double x1):
    cdef double y0
    cdef double y1
    pbwa(x0, x1, &y0, &y1)
    return y0, y1

cpdef double pdtr(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.pdtr"""
    return cephes_pdtr(x0, x1)

cpdef double pdtrc(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.pdtrc"""
    return cephes_pdtrc(x0, x1)

cpdef double pdtri(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.pdtri"""
    if dl_number_t is double:
        return _func_pdtri_unsafe(x0, x1)
    elif dl_number_t is long:
        return cephes_pdtri(x0, x1)
    else:
        return NAN

cpdef double pdtrik(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.pdtrik"""
    return _func_pdtrik(x0, x1)

cpdef double poch(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.poch"""
    return cephes_poch(x0, x1)

cpdef df_number_t powm1(df_number_t x0, df_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.powm1"""
    if df_number_t is float:
        return (<float(*)(float, float) noexcept nogil>scipy.special._ufuncs_cxx._export_powm1_float)(x0, x1)
    elif df_number_t is double:
        return (<double(*)(double, double) noexcept nogil>scipy.special._ufuncs_cxx._export_powm1_double)(x0, x1)
    else:
        if df_number_t is double:
            return NAN
        else:
            return NAN

cdef void pro_ang1(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_ang1"""
    y0[0] = _func_prolate_aswfa_nocv_wrap(x0, x1, x2, x3, y1)

def _pro_ang1_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    pro_ang1(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void pro_ang1_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_ang1_cv"""
    _func_prolate_aswfa_wrap(x0, x1, x2, x3, x4, y0, y1)

def _pro_ang1_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    pro_ang1_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cpdef double pro_cv(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.pro_cv"""
    return _func_prolate_segv_wrap(x0, x1, x2)

cdef void pro_rad1(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_rad1"""
    y0[0] = _func_prolate_radial1_nocv_wrap(x0, x1, x2, x3, y1)

def _pro_rad1_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    pro_rad1(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void pro_rad1_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_rad1_cv"""
    _func_prolate_radial1_wrap(x0, x1, x2, x3, x4, y0, y1)

def _pro_rad1_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    pro_rad1_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cdef void pro_rad2(double x0, double x1, double x2, double x3, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_rad2"""
    y0[0] = _func_prolate_radial2_nocv_wrap(x0, x1, x2, x3, y1)

def _pro_rad2_pywrap(double x0, double x1, double x2, double x3):
    cdef double y0
    cdef double y1
    pro_rad2(x0, x1, x2, x3, &y0, &y1)
    return y0, y1

cdef void pro_rad2_cv(double x0, double x1, double x2, double x3, double x4, double *y0, double *y1) noexcept nogil:
    """See the documentation for scipy.special.pro_rad2_cv"""
    _func_prolate_radial2_wrap(x0, x1, x2, x3, x4, y0, y1)

def _pro_rad2_cv_pywrap(double x0, double x1, double x2, double x3, double x4):
    cdef double y0
    cdef double y1
    pro_rad2_cv(x0, x1, x2, x3, x4, &y0, &y1)
    return y0, y1

cpdef double pseudo_huber(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.pseudo_huber"""
    return _func_pseudo_huber(x0, x1)

cpdef Dd_number_t psi(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.psi"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_cdigamma(_complexstuff.npy_cdouble_from_double_complex(x0)))
    elif Dd_number_t is double:
        return special_digamma(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double radian(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.radian"""
    return cephes_radian(x0, x1, x2)

cpdef double rel_entr(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.rel_entr"""
    return _func_rel_entr(x0, x1)

cpdef Dd_number_t rgamma(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.rgamma"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_crgamma(_complexstuff.npy_cdouble_from_double_complex(x0))) 
    elif Dd_number_t is double:
        return special_rgamma(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double round(double x0) noexcept nogil:
    """See the documentation for scipy.special.round"""
    return cephes_round(x0)

cdef void shichi(Dd_number_t x0, Dd_number_t *y0, Dd_number_t *y1) noexcept nogil:
    """See the documentation for scipy.special.shichi"""
    if Dd_number_t is double_complex:
        _func_cshichi(x0, y0, y1)
    elif Dd_number_t is double:
        _func_cephes_shichi_wrap(x0, y0, y1)
    else:
        if Dd_number_t is double_complex:
            y0[0] = NAN
            y1[0] = NAN
        else:
            y0[0] = NAN
            y1[0] = NAN

def _shichi_pywrap(Dd_number_t x0):
    cdef Dd_number_t y0
    cdef Dd_number_t y1
    shichi(x0, &y0, &y1)
    return y0, y1

cdef void sici(Dd_number_t x0, Dd_number_t *y0, Dd_number_t *y1) noexcept nogil:
    """See the documentation for scipy.special.sici"""
    if Dd_number_t is double_complex:
        _func_csici(x0, y0, y1)
    elif Dd_number_t is double:
        _func_cephes_sici_wrap(x0, y0, y1)
    else:
        if Dd_number_t is double_complex:
            y0[0] = NAN
            y1[0] = NAN
        else:
            y0[0] = NAN
            y1[0] = NAN

def _sici_pywrap(Dd_number_t x0):
    cdef Dd_number_t y0
    cdef Dd_number_t y1
    sici(x0, &y0, &y1)
    return y0, y1

cpdef double sindg(double x0) noexcept nogil:
    """See the documentation for scipy.special.sindg"""
    return cephes_sindg(x0)

cpdef double smirnov(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.smirnov"""
    if dl_number_t is double:
        return _func_smirnov_unsafe(x0, x1)
    elif dl_number_t is long:
        return cephes_smirnov(x0, x1)
    else:
        return NAN

cpdef double smirnovi(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.smirnovi"""
    if dl_number_t is double:
        return _func_smirnovi_unsafe(x0, x1)
    elif dl_number_t is long:
        return cephes_smirnovi(x0, x1)
    else:
        return NAN

cpdef Dd_number_t spence(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.spence"""
    if Dd_number_t is double_complex:
        return _func_cspence(x0)
    elif Dd_number_t is double:
        return cephes_spence(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double complex sph_harm(dl_number_t x0, dl_number_t x1, double x2, double x3) noexcept nogil:
    """See the documentation for scipy.special.sph_harm"""
    if dl_number_t is double:
        return _complexstuff.double_complex_from_npy_cdouble(special_sph_harm_unsafe(x0, x1, x2, x3))
    elif dl_number_t is long:
        return _complexstuff.double_complex_from_npy_cdouble(special_sph_harm(x0, x1, x2, x3))
    else:
        return NAN

cpdef double stdtr(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.stdtr"""
    return _func_stdtr(x0, x1)

cpdef double stdtridf(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.stdtridf"""
    return _func_stdtridf(x0, x1)

cpdef double stdtrit(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.stdtrit"""
    return _func_stdtrit(x0, x1)

cpdef double struve(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.struve"""
    return cephes_struve_h(x0, x1)

cpdef double tandg(double x0) noexcept nogil:
    """See the documentation for scipy.special.tandg"""
    return cephes_tandg(x0)

cpdef double tklmbda(double x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.tklmbda"""
    return cephes_tukeylambdacdf(x0, x1)

cpdef double complex wofz(double complex x0) noexcept nogil:
    """See the documentation for scipy.special.wofz"""
    return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_faddeeva_w)(x0)

cpdef Dd_number_t wrightomega(Dd_number_t x0) noexcept nogil:
    """See the documentation for scipy.special.wrightomega"""
    if Dd_number_t is double_complex:
        return (<double complex(*)(double complex) noexcept nogil>scipy.special._ufuncs_cxx._export_wrightomega)(x0)
    elif Dd_number_t is double:
        return (<double(*)(double) noexcept nogil>scipy.special._ufuncs_cxx._export_wrightomega_real)(x0)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t xlog1py(Dd_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.xlog1py"""
    if Dd_number_t is double:
        return _func_xlog1py[double](x0, x1)
    elif Dd_number_t is double_complex:
        return _func_xlog1py[double_complex](x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t xlogy(Dd_number_t x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.xlogy"""
    if Dd_number_t is double:
        return _func_xlogy[double](x0, x1)
    elif Dd_number_t is double_complex:
        return _func_xlogy[double_complex](x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double y0(double x0) noexcept nogil:
    """See the documentation for scipy.special.y0"""
    return cephes_y0(x0)

cpdef double y1(double x0) noexcept nogil:
    """See the documentation for scipy.special.y1"""
    return cephes_y1(x0)

cpdef double yn(dl_number_t x0, double x1) noexcept nogil:
    """See the documentation for scipy.special.yn"""
    if dl_number_t is double:
        return _func_yn_unsafe(x0, x1)
    elif dl_number_t is long:
        return cephes_yn(x0, x1)
    else:
        return NAN

cpdef Dd_number_t yv(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.yv"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_y(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_y(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef Dd_number_t yve(double x0, Dd_number_t x1) noexcept nogil:
    """See the documentation for scipy.special.yve"""
    if Dd_number_t is double_complex:
        return _complexstuff.double_complex_from_npy_cdouble(special_ccyl_bessel_ye(x0, _complexstuff.npy_cdouble_from_double_complex(x1)))
    elif Dd_number_t is double:
        return special_cyl_bessel_ye(x0, x1)
    else:
        if Dd_number_t is double_complex:
            return NAN
        else:
            return NAN

cpdef double zetac(double x0) noexcept nogil:
    """See the documentation for scipy.special.zetac"""
    return cephes_zetac(x0)

cpdef double wright_bessel(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.wright_bessel"""
    return special_wright_bessel(x0, x1, x2)

cpdef double log_wright_bessel(double x0, double x1, double x2) noexcept nogil:
    """See the documentation for scipy.special.log_wright_bessel"""
    return special_log_wright_bessel(x0, x1, x2)

cpdef double ndtri_exp(double x0) noexcept nogil:
    """See the documentation for scipy.special.ndtri_exp"""
    return _func_ndtri_exp(x0)

cpdef number_t spherical_jn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_jn"""
    if derivative:
        if number_t is double:
            return special_sph_bessel_j_jac(n, z)
        else:
            return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_j_jac(n, _complexstuff.npy_cdouble_from_double_complex(z)))

    if number_t is double:
        return special_sph_bessel_j(n, z)
    else:
        return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_j(n, _complexstuff.npy_cdouble_from_double_complex(z)))

cpdef number_t spherical_yn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_yn"""
    if derivative:
        if number_t is double:
            return special_sph_bessel_y_jac(n, z)
        else:
            return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_y_jac(n, _complexstuff.npy_cdouble_from_double_complex(z)))

    if number_t is double:
        return special_sph_bessel_y(n, z)
    else:
        return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_y(n, _complexstuff.npy_cdouble_from_double_complex(z)))

cpdef number_t spherical_in(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_in"""
    if derivative:
        if number_t is double:
            return special_sph_bessel_i_jac(n, z)
        else:
            return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_i_jac(n, _complexstuff.npy_cdouble_from_double_complex(z)))

    if number_t is double:
        return special_sph_bessel_i(n, z)
    else:
        return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_i(n, _complexstuff.npy_cdouble_from_double_complex(z)))

cpdef number_t spherical_kn(long n, number_t z, bint derivative=0) noexcept nogil:
    """See the documentation for scipy.special.spherical_kn"""
    if derivative:
        if number_t is double:
            return special_sph_bessel_k_jac(n, z)
        else:
            return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_k_jac(n, _complexstuff.npy_cdouble_from_double_complex(z)))

    if number_t is double:
        return special_sph_bessel_k(n, z)
    else:
        return _complexstuff.double_complex_from_npy_cdouble(special_csph_bessel_k(n, _complexstuff.npy_cdouble_from_double_complex(z)))

def _bench_airy_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.airy(x0)

def _bench_airy_d_cy(int N, double x0):
    cdef int n
    cdef double y0
    cdef double y1
    cdef double y2
    cdef double y3
    for n in range(N):
        airy(x0, &y0, &y1, &y2, &y3)

def _bench_airy_D_py(int N, double complex x0):
    cdef int n
    for n in range(N):
        _ufuncs.airy(x0)

def _bench_airy_D_cy(int N, double complex x0):
    cdef int n
    cdef double complex y0
    cdef double complex y1
    cdef double complex y2
    cdef double complex y3
    for n in range(N):
        airy(x0, &y0, &y1, &y2, &y3)

def _bench_beta_dd_py(int N, double x0, double x1):
    cdef int n
    for n in range(N):
        _ufuncs.beta(x0, x1)

def _bench_beta_dd_cy(int N, double x0, double x1):
    cdef int n
    for n in range(N):
        beta(x0, x1)

def _bench_erf_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.erf(x0)

def _bench_erf_d_cy(int N, double x0):
    cdef int n
    for n in range(N):
        erf(x0)

def _bench_erf_D_py(int N, double complex x0):
    cdef int n
    for n in range(N):
        _ufuncs.erf(x0)

def _bench_erf_D_cy(int N, double complex x0):
    cdef int n
    for n in range(N):
        erf(x0)

def _bench_exprel_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.exprel(x0)

def _bench_exprel_d_cy(int N, double x0):
    cdef int n
    for n in range(N):
        exprel(x0)

def _bench_gamma_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.gamma(x0)

def _bench_gamma_d_cy(int N, double x0):
    cdef int n
    for n in range(N):
        gamma(x0)

def _bench_gamma_D_py(int N, double complex x0):
    cdef int n
    for n in range(N):
        _ufuncs.gamma(x0)

def _bench_gamma_D_cy(int N, double complex x0):
    cdef int n
    for n in range(N):
        gamma(x0)

def _bench_jv_dd_py(int N, double x0, double x1):
    cdef int n
    for n in range(N):
        _ufuncs.jv(x0, x1)

def _bench_jv_dd_cy(int N, double x0, double x1):
    cdef int n
    for n in range(N):
        jv(x0, x1)

def _bench_jv_dD_py(int N, double x0, double complex x1):
    cdef int n
    for n in range(N):
        _ufuncs.jv(x0, x1)

def _bench_jv_dD_cy(int N, double x0, double complex x1):
    cdef int n
    for n in range(N):
        jv(x0, x1)

def _bench_loggamma_D_py(int N, double complex x0):
    cdef int n
    for n in range(N):
        _ufuncs.loggamma(x0)

def _bench_loggamma_D_cy(int N, double complex x0):
    cdef int n
    for n in range(N):
        loggamma(x0)

def _bench_logit_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.logit(x0)

def _bench_logit_d_cy(int N, double x0):
    cdef int n
    for n in range(N):
        logit(x0)

def _bench_psi_d_py(int N, double x0):
    cdef int n
    for n in range(N):
        _ufuncs.psi(x0)

def _bench_psi_d_cy(int N, double x0):
    cdef int n
    for n in range(N):
        psi(x0)

def _bench_psi_D_py(int N, double complex x0):
    cdef int n
    for n in range(N):
        _ufuncs.psi(x0)

def _bench_psi_D_cy(int N, double complex x0):
    cdef int n
    for n in range(N):
        psi(x0)
