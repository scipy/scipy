#include <stdio.h>

main()
{
    /* Local variables */
    double base, emin, prec, emax, rmin, rmax, t, sfmin;
    extern double dlamch_(char *);
    double rnd, eps;

    eps = dlamch_("Epsilon");
    sfmin = dlamch_("Safe minimum");
    base = dlamch_("Base");
    prec = dlamch_("Precision");
    t = dlamch_("Number of digits in mantissa");
    rnd = dlamch_("Rounding mode");
    emin = dlamch_("Minnimum exponent");
    rmin = dlamch_("Underflow threshold");
    emax = dlamch_("Largest exponent");
    rmax = dlamch_("Overflow threshold");

    printf(" Epsilon                      = %e\n", eps);
    printf(" Safe minimum                 = %e\n", sfmin);
    printf(" Base                         = %.0f\n", base);
    printf(" Precision                    = %e\n", prec);
    printf(" Number of digits in mantissa = %.0f\n", t);
    printf(" Rounding mode                = %.0f\n", rnd);
    printf(" Minimum exponent             = %.0f\n", emin);
    printf(" Underflow threshold          = %e\n", rmin);
    printf(" Largest exponent             = %.0f\n", emax);
    printf(" Overflow threshold           = %e\n", rmax);
    printf(" Reciprocal of safe minimum   = %e\n", 1./sfmin);

    return 0;
}
