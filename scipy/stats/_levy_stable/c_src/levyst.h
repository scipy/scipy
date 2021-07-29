#ifndef LEVYST_H
#define LEVYST_H

struct nolan_precanned
{
    double (*g)(struct nolan_precanned *, double);
    double alpha;
    double zeta;
    double xi;
    double zeta_prefactor;
    double alpha_exp;
    double alpha_xi;
    double zeta_offset;
    double two_beta_div_pi;
    double pi_div_two_beta;
    double x0_div_term;
    double c1;
    double c2;
    double c3;
};

typedef double (*g_callback)(struct nolan_precanned *, double);

extern struct nolan_precanned *
nolan_precan(double, double, double);
 
#endif
