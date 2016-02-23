#include <float.h>
#include <math.h>
#include <stdio.h>

double dmach(char *cmach)
{
/*  -- SuperLU auxiliary routine (version 5.0) --   
    This uses C99 standard constants, and is thread safe.

    Must be compiled with -std=c99 flag.


    Purpose   
    =======   

    DMACH returns double precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DMACH:   
            = 'E' or 'e',   DMACH := eps   
            = 'S' or 's ,   DMACH := sfmin   
            = 'B' or 'b',   DMACH := base   
            = 'P' or 'p',   DMACH := eps*base   
            = 'N' or 'n',   DMACH := t   
            = 'R' or 'r',   DMACH := rnd   
            = 'M' or 'm',   DMACH := emin   
            = 'U' or 'u',   DMACH := rmin   
            = 'L' or 'l',   DMACH := emax   
            = 'O' or 'o',   DMACH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/

    double sfmin, small, rmach;

    if (strncmp(cmach, "E", 1)==0) {
	rmach = DBL_EPSILON * 0.5;
    } else if (strncmp(cmach, "S", 1)==0) {
	sfmin = DBL_MIN;
	small = 1. / DBL_MAX;
	if (small >= sfmin) {
	    /* Use SMALL plus a bit, to avoid the possibility of rounding   
	       causing overflow when computing  1/sfmin. */
	    sfmin = small * (DBL_EPSILON*0.5 + 1.);
	}
	rmach = sfmin;
    } else if (strncmp(cmach, "B", 1)==0) {
	rmach = FLT_RADIX;
    } else if (strncmp(cmach, "P", 1)==0) {
	rmach = DBL_EPSILON * 0.5 * FLT_RADIX;
    } else if (strncmp(cmach, "N", 1)==0) {
	rmach = DBL_MANT_DIG;
    } else if (strncmp(cmach, "R", 1)==0) {
	rmach = FLT_ROUNDS;
    } else if (strncmp(cmach, "M", 1)==0) {
	rmach = DBL_MIN_EXP;
    } else if (strncmp(cmach, "U", 1)==0) {
	rmach = DBL_MIN;
    } else if (strncmp(cmach, "L", 1)==0) {
	rmach = DBL_MAX_EXP;
    } else if (strncmp(cmach, "O", 1)==0) {
	rmach = DBL_MAX;
    }

    return rmach;

} /* end dmach */
