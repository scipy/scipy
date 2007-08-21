#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void compute_root_from_lambda(double, double *, double *);


void
compute_root_from_lambda(lambda, r, omega)
     double lambda;
     double *r;
     double *omega;
{
    double xi;
    double tmp, tmp2;

    tmp = sqrt(3 + 144*lambda);
    xi = 1 - 96*lambda + 24*lambda * tmp;
    *omega = atan(sqrt((144*lambda - 1.0)/xi));
    tmp2 = sqrt(xi);
    *r = (24*lambda - 1 - tmp2)/(24*lambda) \
	* sqrt((48*lambda + 24*lambda*tmp))/tmp2;
    return;
}
