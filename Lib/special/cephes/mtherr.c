/*							mtherr.c
 *
 *	Library common error handling routine
 *
 *
 *
 * SYNOPSIS:
 *
 * char *fctnam;
 * int code;
 * int mtherr();
 *
 * mtherr( fctnam, code );
 *
 *
 *
 * DESCRIPTION:
 *
 * This routine may be called to report one of the following
 * error conditions (in the include file mconf.h).
 *  
 *   Mnemonic        Value          Significance
 *
 *    DOMAIN            1       argument domain error
 *    SING              2       function singularity
 *    OVERFLOW          3       overflow range error
 *    UNDERFLOW         4       underflow range error
 *    TLOSS             5       total loss of precision
 *    PLOSS             6       partial loss of precision
 *    EDOM             33       Unix domain error code
 *    ERANGE           34       Unix range error code
 *
 * The default version of the file prints the function name,
 * passed to it by the pointer fctnam, followed by the
 * error condition.  The display is directed to the standard
 * output device.  The routine then returns to the calling
 * program.  Users may wish to modify the program to abort by
 * calling exit() under severe error conditions such as domain
 * errors.
 *
 * Since all error conditions pass control to this function,
 * the display may be easily changed, eliminated, or directed
 * to an error logging device.
 *
 * SEE ALSO:
 *
 * mconf.h
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <stdio.h>
#include "mconf.h"

int merror = 0;

/* Notice: the order of appearance of the following
 * messages is bound to the error codes defined
 * in mconf.h.
 */
static char *ermsg[7] = {
"unknown",      /* error code 0 */
"domain",       /* error code 1 */
"singularity",  /* et seq.      */
"overflow",
"underflow",
"total loss of precision",
"partial loss of precision"
};


int mtherr( name, code )
char *name;
int code;
{

/* Display string passed by calling program,
 * which is supposed to be the name of the
 * function in which the error occurred:
 */
printf( "\n%s ", name );

/* Set global error message word */
merror = code;

/* Display error message defined
 * by the code argument.
 */
if( (code <= 0) || (code >= 7) )
	code = 0;
printf( "%s error\n", ermsg[code] );

/* Return to calling
 * program
 */
return( 0 );
}
