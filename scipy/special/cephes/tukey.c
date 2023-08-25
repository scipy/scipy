
/* Compute the CDF of the Tukey-Lambda distribution 
 * using a bracketing search with special checks
 * 
 * The PPF of the Tukey-lambda distribution is 
 * G(p) = (p**lam + (1-p)**lam) / lam
 * 
 * Author:  Travis Oliphant 
 */

#include <Python.h>
#include <math.h>

#define SMALLVAL 1e-4
#define EPS 1.0e-14
#define MAXCOUNT 60

double tukeylambdacdf(double x, double lmbda)
{
    double pmin, pmid, pmax, plow, phigh, xeval;
    int count;

    if (isnan(x) || isnan(lmbda)) {
        return NAN;
    }

    xeval = 1.0 / lmbda;
    if (lmbda > 0.0) {
        if (x <= (-xeval)) {
            return 0.0;
        }
        if (x >= xeval) {
            return 1.0;
        }
    }

    if ((-SMALLVAL < lmbda) && (lmbda < SMALLVAL)) {
        if (x >= 0) {
            return 1.0 / (1.0 + exp(-x));
        }
        else {
            return exp(x) / (1.0 + exp(x));
        }
    }

    pmin = 0.0;
    pmid = 0.5;
    pmax = 1.0;
    plow = pmin;
    phigh = pmax;
    count = 0;

    while ((count < MAXCOUNT) && (fabs(pmid - plow) > EPS)) {
        xeval = (pow(pmid, lmbda) - pow(1.0 - pmid, lmbda)) / lmbda;
        if (xeval == x) {
            return pmid;
        }
        if (xeval > x) {
            phigh = pmid;
            pmid = (pmid + plow) / 2.0;
        }
        else {
            plow = pmid;
            pmid = (pmid + phigh) / 2.0;
        }
        count++;
    }
    return pmid;
}
