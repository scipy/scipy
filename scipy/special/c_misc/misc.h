#ifndef C_MISC_MISC_H
#define C_MISC_MISC_H

typedef enum {
  /* An exact solution was found, in which case the first point
     on the interval is the value */
  FSOLVE_EXACT,
  /* Interval width is less than the tolerance */
  FSOLVE_CONVERGED,
  /* Not a bracket */
  FSOLVE_NOT_BRACKET,
  /* Root-finding didn't converge in a set number of iterations. */
  FSOLVE_MAX_ITERATIONS
} fsolve_result_t;

typedef double (*objective_function)(double, void *);

fsolve_result_t false_position(double *a, double *fa, double *b, double *fb,
                       objective_function f, void *f_extra,
                       double abserr, double relerr, double bisect_til,
                       double *best_x, double *best_f, double *errest);

double besselpoly(double a, double lambda, double nu);
double gammaincinv(double a, double x);

#define gammaincinv_doc """gammaincinv(a, y) returns x such that gammainc(a, x) = y."""

#endif /* C_MISC_MISC_H */
