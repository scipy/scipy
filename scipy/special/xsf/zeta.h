/* Complex riemann-zeta function implementation based on Python implementation
 * written by Matt Haberland (@mdhaber) in:
 * https://colab.research.google.com/drive/1zMDSAJlXCLRqMMtJ0e9nDGjQ8iZCnmn5?usp=sharing
 */

#pragma once

#include "config.h"
#include "error.h"
#include "gamma.h"
#include "trig.h"

#include "cephes/const.h"
#include "cephes/zeta.h"
#include "cephes/zetac.h"

namespace xsf {

namespace detail {

    /* Log of absolute value of expansion coefficients for Euler-Maclaurin
     * summation formula. log(|B2k / (2k)!|)
     *
     * See https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms
     *
     * Generated with the script
     *
     * import numpy as np
     * from mpmath import mp
     *
     * mp.dps = 10000
     *
     * results = []
     * for k in range(51):
     *     results.append(
     *         float(
     *             mp.log(abs(mp.bernoulli(2*k)/(mp.factorial(2*k))))
     *         )
     *     )
     */
    constexpr double zeta_em_log_abs_coeff_lookup[] = {
	0.0,
	-2.4849066497880004,
	-6.579251212010101,
	-10.31692083029347,
	-14.005800284407405,
	-17.68462940266784,
	-21.361131560073222,
	-25.037070502911423,
	-28.712870599846948,
	-32.388636197522295,
	-36.06439319366539,
	-39.740148041995184,
	-43.41590235365616,
	-47.091656531181485,
	-50.76741067517639,
	-54.44316481078909,
	-58.11891894430628,
	-61.79467307729959,
	-65.47042721016194,
	-69.14618134299154,
	-72.82193547581296,
	-76.49768960863234,
	-80.1734437414512,
	-83.84919787426993,
	-87.52495200708863,
	-91.20070613990733,
	-94.87646027272602,
	-98.55221440554472,
	-102.2279685383634,
	-105.9037226711821,
	-109.57947680400078,
	-113.25523093681947,
	-116.93098506963817,
	-120.60673920245685,
	-124.28249333527555,
	-127.95824746809424,
	-131.63400160091294,
	-135.30975573373163,
	-138.9855098665503,
	-142.661263999369,
	-146.3370181321877,
	-150.0127722650064,
	-153.6885263978251,
	-157.36428053064375,
	-161.04003466346245,
	-164.71578879628115,
	-168.39154292909984,
	-172.06729706191854,
	-175.74305119473723,
	-179.4188053275559,
	-183.0945594603746
    };

    // Complex log of expansion coefficients for Euler-Maclaurin summation formula.
    XSF_HOST_DEVICE inline std::complex<double> zeta_em_log_coeff(std::size_t n) {
	std::complex<double> J(0.0, 1.0);
	std::complex<double> result;
	if (n < 50) {
	    result = zeta_em_log_abs_coeff_lookup[n];
	} else {
	    /* Asymptotic formula
	     * Uses https://dlmf.nist.gov/24.11#E1 to approximate B_{2n} and
	     * Stirling's approximation for (2n)!.
	     */
	    result = std::log(2.0) - 2.0*n*std::log(2*M_PI);
	}
	if (n % 2 == 0) {
	    /* B_{2n}/(2n)! is negative for even n. This contributes a term
	     * pi*i when taking the log. */
	    result += M_PI * J;
	}
	return result;
    }

    /* Compute riemann_zeta for complex input z using the Euler-Maclaurin formula.
     * Computation of individual terms in expansion are logarithmized to avoid
     * overflow. TODO: only logarithmize when necessary. */
    XSF_HOST_DEVICE inline std::complex<double> zeta_euler_maclaurin(std::complex<double> z) {
	if (z == 1.0) {
	    /* Return NaN at pole since value depends on how z approaches 1.0. */
	    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
	}
	std::size_t n = static_cast<std::size_t>(std::max(std::abs(z.imag()) / 4.0, 50.0));
	std::size_t m = n;
	std::complex<double> result = 0.0;
	for (std::size_t i = 1; i < n; i++) {
	    std::complex<double> term = std::pow(static_cast<double>(i), -z);
	    result += term;
	    // When z.real() > 1, series converges and we can consider early termination
	    if (z.real() > 1 && std::abs(term) / std::abs(result) <= std::numeric_limits<double>::epsilon()) {
	     	return result;
	    }
	}
	double N = static_cast<double>(n);
	std::complex<double> b = std::pow(n, -z);
	result += b * (0.5 + N / (z - 1.0));
	/* The terms of the Euler-Maclaurin
	 * expansion below are T(k, n) = B2k/(2k)! * n^(1 - z - 2k) * z(z+1)...(z+2k-2).
	 * We work with logarithms to avoid overflow in all cases at the expense of
	 * some accuracy. At the start of iteration k:
	 *     log_poch will equal log(z(z+1)...(z+2k-2))
	 *     log_factor will equal log(n^(1 - z - 2k))
	 * These are updated one extra time after the loop completes for use in the
	 * Euler-Maclaurin error estimate.
	 */
	std::complex<double> log_poch = std::log(z);
	std::complex<double> log_factor = -(z + 1.0) * std::log(N);
	for (std::size_t k = 1; k <= m; k++) {
	    std::complex<double> term = std::exp(zeta_em_log_coeff(k) + log_factor + log_poch);
	    result += term;
	    if (std::abs(term)/std::abs(result) <= std::numeric_limits<double>::epsilon()) {
		return result;
	    }
	    log_poch += std::log(z + static_cast<double>(2*k - 1)) + std::log(z + static_cast<double>(2*k));
	    log_factor -= 2*std::log(N);
	}
	/* Euler-maclaurin absolute error estimate.
	 * The error is bounded above by |(z + 2m + 1)/(z.real + 2m + 1) * T(m+1, n)|
	 * See https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms
	 */
	double error;
	error = std::abs(std::exp(zeta_em_log_coeff(m + 1) + log_factor + log_poch));
	error *= std::abs((z + 2.0*m + 1.0)/(z.real() + 2.0*m + 1.0));
	// convert to relative error estimate
	error /= std::abs(result);
	if (error > 1e-8) {
	    if (error > 1e-1) {
		/* If error estimate predicts we don't even get 1 digit of precision, return NaN
		 * and signal no result */
		set_error("zeta", SF_ERROR_NO_RESULT, NULL);
		return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
	    }
	    // Signal reduced precision.
	    set_error("zeta", SF_ERROR_LOSS, NULL);
	}
	return result;
    }

    /* Lookup table of coefficients for Algorithm 2 from Borwein 1995
     * Borwein, Peter B.. “An efficient algorithm for the Riemann zeta function.” (1995).
     *
     * Stores coefficients as dk / dn, where dn is the final coefficient.
     *
     * Generated with the Python script:
     *
     * import numpy as np
     * import math
     * from mpmath import mp
     *
     * mp.dps = 1000
     * n = 50
     *
     * coeffs = []
     * S = mp.zero
     * for i in range(n + 1):
     *     num = math.factorial(n + i - 1) * 4**i
     *     den = math.factorial(n - i) * math.factorial(2*i)
     *     S += mp.mpf(num) / mp.mpf(den)
     *     coeffs.append(S*n)
     *
     * dn = coeffs[-1]
     * coeffs = [float(dk/dn) for dk in coeffs[:-1]]
     * coeffs = np.asarray(coeffs)
     */
    constexpr double zeta_borwein_coeff[] = {
	1.0555078361382878e-38,
	5.278594688527578e-35,
	4.4014687322044963e-32,
	1.467453546497519e-29,
	2.617862196688831e-27,
	2.900097799958025e-25,
	2.184440361492933e-23,
	1.1890977312913296e-21,
	4.8871396166872276e-20,
	1.5672253698802734e-18,
	4.022931234264572e-17,
	8.435973533351745e-16,
	1.469296379914116e-14,
	2.1548747054571902e-13,
	2.6919530537535124e-12,
	2.8925409162768484e-11,
	2.6957505693699856e-10,
	2.194772239130839e-09,
	1.57078229370057e-08,
	9.936187220749133e-08,
	5.581721578217702e-07,
	2.796271112037765e-06,
	1.253886254275813e-05,
	5.049261002939051e-05,
	0.00018312884459703666,
	0.0005997690328552426,
	0.0017780501082460968,
	0.004781802283665968,
	0.011690432287131671,
	0.026034302929535964,
	0.052922982472754856,
	0.09842471411648858,
	0.16789610796540344,
	0.26350429194768626,
	0.3819442810600665,
	0.5137731385068898,
	0.645292538541866,
	0.7625449322050362,
	0.8556063057019102,
	0.921056062886525,
	0.9616100580028147,
	0.9835905431606953,
	0.9939187229336753,
	0.9980782525166648,
	0.9994930141459889,
	0.999891478844585,
	0.9999819091990203,
	0.9999977981288319,
	0.9999998260580315,
	0.9999999933099243
    };

    /* Compute riemann_zeta for complex input z using Algorithm 2 from Borwein 1995. */
    XSF_HOST_DEVICE inline std::complex<double> zeta_borwein(std::complex<double> z) {
	std::complex<double> result = 0.0;
	// Sum in reverse order because smaller terms come later.
	for (int k = 49; k >= 0; k--) {
	    double sign = std::pow(-1.0, k);
	    std::complex<double> den = std::pow(k + 1, z);
	    std::complex<double> term = sign * (zeta_borwein_coeff[k] - 1.0) / den;
	    result += term;
	}
	return result * -1.0/(1.0 - std::pow(2.0, 1.0 - z));
    }

    /* Compute riemann zeta for complex z and real part >= 0 */
    XSF_HOST_DEVICE inline std::complex<double> zeta_right_halfplane(std::complex<double> z) {
	if (z == 1.0) {
	    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
	}
	/* Cutoff for using Euler-MacLaurin chosen based on cursory empirical search.
	 * TODO: Choose cutoffs in a more principled way. */
	if (z.real() < 50.0 && std::abs(z.imag()) > 50.0) {
	    if (z.real() >= 0.0 && z.real() < 2.5 && std::abs(z.imag()) > 1e9) {
		/* Euler-MacLaurin summation starts to take an unreasonable amount of time in this
		 * region, so just give up and return NaN instead. */
		set_error("zeta", SF_ERROR_NO_RESULT, NULL);
		return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
	    }
	    return zeta_euler_maclaurin(z);
	}
	return zeta_borwein(z);
    }

    XSF_HOST_DEVICE inline std::complex<double> exppi(std::complex<double> z) {
	    // exp(pi*z) for complex z.
	    double x = z.real();
	    double y = z.imag();
	    std::complex<double> factor1(xsf::cospi(y), xsf::sinpi(y));
	    double factor2 = std::exp(M_PI*x);
	    return factor1 * factor2;
	}

    XSF_HOST_DEVICE inline std::complex<double> logsinpi(std::complex<double> z) {
	/* log(sinpi(z)) using sin(z) = (exp(i*pi*z) - exp(-i*pi*z)) / 2i
	 *
	 * No attempt is made to choose any particular branch of the logarithm.
	 * This is an internal function and the intent is that this that the
	 * result of log(sinpi(z)) will be added to other terms, and the sum
	 * will then be exponentiated, making the choice of a specific branch
	 * unnecessary.
	 */
	std::complex<double> result = std::log(xsf::sinpi(z));
	// If it doesn't overflow, just do the regular calculation.
	if (std::isfinite(result.real()) && !std::isfinite(result.imag())) {
	    return result;
	}
	/* Otherwise factor before taking log. This is where we may end up
	 * taking a branch other than the principal branch. */
	std::complex<double> J(0.0, 1.0);
	/* Calculating log((exp(i*pi*z) - exp(-i*pi*z)) / 2i). Factor out term
	 * with larger magnitude before taking log. */
	if (z.imag() > 0 ) {
	    /* if z.imag() > 0 then, exp(-i*pi*z) has greatest magnitude. Factor it
	     * out to get:
	     * log(exp(-i*pi*z)*((exp(2*i*pi*z) - 1.0)/(2i)) =
	     * log(exp(-i*pi*z)) + log((exp(2*i*pi*z) - 1.0)/(2i)) =
	     * -i*pi*z + log((exp(2*i*pi*z) - 1.0)/(2i)) */
	    return -J * M_PI * z + std::log((exppi(2.0 * z * J) - 1.0) / (2.0*J));
	}
	/* if z.imag() < 0 then, exp(i*pi*z) has greatest magnitude. Factor similarly
	 * to above */
	return J * M_PI * z + std::log((1.0 - exppi(-2.0 * z * J)) / (2.0*J));
    }

    /* Leading factor in reflection formula for zeta function.
     * zeta(z) = 2^z * pi^(z-1) * sin(pi*z/2) * gamma(1 - z) * zeta(1 - z)
     * This computes 2^z * pi^(z - 1) * sin(pi*z/2) * gamma(1 - z)
     *
     * Computation is logarithimized to prevent overflow.
     * TODO: Complexify the cephes zeta_reflection implementation, which uses
     * the lanczos approximation for the gamma function. */
    XSF_HOST_DEVICE inline std::complex<double> zeta_reflection_factor_with_logs(std::complex<double> z) {
 	std::complex<double> t1 = z * M_LN2;
	std::complex<double> t2 = (z - 1.0) * xsf::cephes::detail::LOGPI;
	std::complex<double> t3 = logsinpi(z / 2.0);
	std::complex<double> t4 = xsf::loggamma(1.0 - z);
	std::complex<double> factor = std::exp(t1 + t2 + t3 + t4);
	return factor;
    }

    XSF_HOST_DEVICE inline std::complex<double> zeta_reflection(std::complex<double> z) {
	std::complex<double> factor = 2.0 * std::pow(2*M_PI, z - 1.0) * xsf::sinpi(z/2.0) * xsf::gamma(1.0 - z);
	if (!std::isfinite(factor.real()) || !std::isfinite(factor.imag())) {
	    // Try again with logs if standard calculation had overflow.
	    factor = zeta_reflection_factor_with_logs(z);
	}
	std::complex<double> result = zeta_right_halfplane(1.0 - z);
	/* zeta tends to 1.0 as real part tends to +inf. In cases where
	 * the real part of zeta tends to -inf, then zeta(1 - z) in the
	 * reflection formula will tend to 1.0. Factor overflows then,
	 * factor * result below will become NaN. In this case, we just
	 * return factor to preserve complex infinity. Only zeta(1 - z) == 1.0
	 * is handled because this is the only practical case where we should
	 *  expect zeta(1 - z) == x for a real number x when z is not on the
	 * real line. */
	return (result == 1.0) ? factor : factor * result;
    }
}

XSF_HOST_DEVICE inline std::complex<double> riemann_zeta(std::complex<double> z) {
    if (z.imag() == 0.0) {
	return cephes::riemann_zeta(z.real());
    }
    if (z.real() >= 0.5) {
	return detail::zeta_right_halfplane(z);
    }
    return detail::zeta_reflection(z);
}

XSF_HOST_DEVICE inline std::complex<float> riemann_zeta(std::complex<float> z) {
    return static_cast<std::complex<float>>(riemann_zeta(static_cast<std::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double riemann_zeta(double x) { return cephes::riemann_zeta(x); }

XSF_HOST_DEVICE inline float riemann_zeta(float x) { return riemann_zeta(static_cast<double>(x)); }

XSF_HOST_DEVICE inline double zeta(double x, double q) { return cephes::zeta(x, q); }

XSF_HOST_DEVICE inline float zeta(float x, float q) { return zeta(static_cast<double>(x), static_cast<double>(q)); }

XSF_HOST_DEVICE inline std::complex<double> zeta(std::complex<double> z, double q) {
    if (z.imag() == 0.0) {
	return zeta(z.real(), q);
    }
    // Complex input for Hurwitz Zeta is not currently supported.
    set_error("zeta", SF_ERROR_DOMAIN, NULL);
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
}

XSF_HOST_DEVICE inline std::complex<float> zeta(std::complex<float> z, float q) {
    return static_cast<std::complex<float>>(zeta(static_cast<std::complex<double>>(z), static_cast<float>(q)));
}

XSF_HOST_DEVICE inline double zetac(double x) { return cephes::zetac(x); }

XSF_HOST_DEVICE inline float zetac(float x) { return zetac(static_cast<double>(x)); }

} // namespace xsf
