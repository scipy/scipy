#pragma once

/* Use this header to include functions from cephes in compiled special
 * function kernels written in C++.

 * cephes_names.h defines aliases airy -> cephes_airy etc.
 * which causes trouble when one tries to use functions from cephes in the
 * same translation unit where boost is used, due to name clashes. We undef
 * all of these aliases and disambiguate the cephes functions by putting them
 * in a special::cephes namespace.
 */


// Namespace the include to avoid polluting global namespace.
namespace special {
    namespace cephes {
	namespace detail {

#include "cephes.h"

#undef airy
#undef bdtrc
#undef bdtr
#undef bdtri
#undef besselpoly
#undef beta
#undef lbeta
#undef btdtr
#undef cbrt
#undef chdtrc
#undef chbevl
#undef chdtr
#undef chdtri
#undef dawsn
#undef ellie
#undef ellik
#undef ellpe
#undef ellpj
#undef ellpk
#undef exp10
#undef exp2
#undef expn
#undef fdtrc
#undef fdtr
#undef fdtri
#undef fresnl
#undef Gamma
#undef lgam
#undef lgam_sgn
#undef gammasgn
#undef gdtr
#undef gdtrc
#undef gdtri
#undef hyp2f1
#undef hyperg
#undef i0
#undef i0e
#undef i1
#undef i1e
#undef igamc
#undef igam
#undef igami
#undef incbet
#undef incbi
#undef iv
#undef j0
#undef y0
#undef j1
#undef y1
#undef jn
#undef jv
#undef k0
#undef k0e
#undef k1
#undef k1e
#undef kn
#undef nbdtrc
#undef nbdtr
#undef nbdtri
#undef ndtr
#undef erfc
#undef erf
#undef erfinv
#undef erfcinv
#undef ndtri
#undef pdtrc
#undef pdtr
#undef pdtri
#undef poch
#undef psi
#undef rgamma
#undef riemann_zeta
#undef round
#undef shichi
#undef sici
#undef radian
#undef sindg
#undef sinpi
#undef cosdg
#undef cospi
#undef sincos
#undef spence
#undef stdtr
#undef stdtri
#undef struve_h
#undef struve_l
#undef struve_power_series
#undef struve_asymp_large_z
#undef struve_bessel_series
#undef yv
#undef tandg
#undef cotdg
#undef log1p
#undef expm1
#undef cosm1
#undef yn
#undef zeta
#undef zetac
#undef smirnov
#undef smirnovc
#undef smirnovi
#undef smirnovci
#undef smirnovp
#undef kolmogorov
#undef kolmogi
#undef kolmogp
#undef kolmogc
#undef kolmogci
#undef owens_t
	}
    }
}

namespace special {
    namespace cephes {
	// Functions are being added as needed.

	inline double beta(double a, double b) {
	    return detail::cephes_beta(a, b);
	}

	inline double lbeta(double a, double b) {
	    return detail::cephes_lbeta(a, b);
	}

	inline double Gamma(double x) {
	    return detail::cephes_Gamma(x);
	}

    }
}
