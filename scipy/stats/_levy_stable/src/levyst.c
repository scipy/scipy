/*
 * Implements special functions for stable distribution calculations.
 *
 * A function g appears in the integrand in Nolan's method for calculating
 * stable densities and distribution functions. It takes a different form for
 * alpha = 1 vs alpha ≠ 1. See [NO] for more info.
 *
 * References
 * [NO] John P. Nolan (1997) Numerical calculation of stable densities and
 *      distribution functions.
 */
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "levyst.h"

/* M_PI et al. are not defined in math.h in C99, even with _USE_MATH_DEFINES */
#ifndef M_PI_2
# define M_PI_2  1.57079632679489661923  /* pi/2 */
# define M_1_PI  0.31830988618379067154  /* 1/pi */
# define M_2_PI  0.63661977236758134308  /* 2/pi */
#endif

double
g_alpha_ne_one(struct nolan_precanned *sp, double theta)
{
    if (theta == -sp->xi) {
        if (sp->alpha < 1) {
            return 0;
        }
        else {
            return INFINITY;
        }
    }
    if (theta == M_PI_2) {
        if (sp->alpha < 1) {
            return INFINITY;
        }
        else {
            return 0;
        }
    }

    double cos_theta = cos(theta);
    return (
        sp->zeta_prefactor
        * pow(
            cos_theta
            / sin(sp->alpha_xi + sp->alpha * theta)
            * sp->zeta_offset, sp->alpha_exp)
        * cos(sp->alpha_xi + (sp->alpha - 1) * theta)
        / cos_theta
    );
}

double
g_alpha_eq_one(struct nolan_precanned *sp, double theta)
{
    if (theta == -sp->xi) {
        return 0;
    }

    if (theta == M_PI_2) {
        return INFINITY;
    }

    return (
        (1 + theta * sp->two_beta_div_pi)
        * exp((sp->pi_div_two_beta + theta) * tan(theta) - sp->x0_div_term)
        / cos(theta)
    );
}

struct nolan_precanned *
nolan_precan(double alpha, double beta, double x0)
{
    /* Stores results of intermediate computations so they need not be
     * recomputed when g is called many times during numerical integration
     * through QUADPACK.
     */
    struct nolan_precanned *sp = malloc(sizeof(struct nolan_precanned));
    if (!sp) {
        abort();
    }
    sp->alpha = alpha;
    sp->zeta = -beta * tan(M_PI_2 * alpha);

    if (alpha != 1.) {
        sp->xi = atan(-sp->zeta) / alpha;
        sp->zeta_prefactor = pow(
            pow(sp->zeta, 2.) + 1., -1. / (2. * (alpha - 1.)));
        sp->alpha_exp = alpha / (alpha - 1.);
        sp->alpha_xi = atan(-sp->zeta);
        sp->zeta_offset = x0 - sp->zeta;
        if (alpha < 1.) {
            sp->c1 = 0.5 - sp->xi * M_1_PI;
            sp->c3 = M_1_PI;
        }
        else {
            sp->c1 = 1.;
            sp->c3 = -M_1_PI;
        }
        sp->c2 = alpha * M_1_PI / fabs(alpha - 1.) / (x0 - sp->zeta);
        sp->g = &g_alpha_ne_one;
    }
    else {
        sp->xi = M_PI_2;
        sp->two_beta_div_pi = beta * M_2_PI;
        sp->pi_div_two_beta = M_PI_2 / beta;
        sp->x0_div_term = x0 / sp->two_beta_div_pi;
        sp->c1 = 0.;
        sp->c2 = .5 / fabs(beta);
        sp->c3 = M_1_PI;
        sp->g = &g_alpha_eq_one;
    }
    return sp;
}
