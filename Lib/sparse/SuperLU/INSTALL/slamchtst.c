#include <stdio.h>

main()
{
    /* Local variables */
    float base, emin, prec, emax, rmin, rmax, t, sfmin;
    extern double slamch_(char *);
    float rnd, eps;

    eps = slamch_("Epsilon");
    sfmin = slamch_("Safe minimum");
    base = slamch_("Base");
    prec = slamch_("Precision");
    t = slamch_("Number of digits in mantissa");
    rnd = slamch_("Rounding mode");
    emin = slamch_("Minnimum exponent");
    rmin = slamch_("Underflow threshold");
    emax = slamch_("Largest exponent");
    rmax = slamch_("Overflow threshold");

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
