/*							fabs.c
 *
 *		Absolute value
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y;
 *
 * y = fabs( x );
 *
 *
 *
 * DESCRIPTION:
 * 
 * Returns the absolute value of the argument.
 *
 */


#include "mconf.h"
/* Avoid using UNK if possible.  */
#ifdef UNK
#if BIGENDIAN
#define MIEEE 1
#else
#define IBMPC 1
#endif
#endif

double fabs(double x)
{
union
  {
    double d;
    short i[4];
  } u;

u.d = x;
#ifdef IBMPC
    u.i[3] &= 0x7fff;
#endif
#ifdef MIEEE
    u.i[0] &= 0x7fff;
#endif
#ifdef DEC
    u.i[3] &= 0x7fff;
#endif
#ifdef UNK
if( u.d < 0 )
   u.d = -u.d;
#endif
return( u.d );
}
