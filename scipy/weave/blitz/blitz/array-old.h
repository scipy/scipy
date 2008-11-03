/***************************************************************************
 * blitz/array-old.h  Maximal include version of Array<P_numtype, N_rank>
 *                    Note: see <blitz/array-impl.h> for the class def.
 *
 * $Id: array-old.h 1414 2005-11-01 22:04:59Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

#ifndef BZ_ARRAY_OLD_H
#define BZ_ARRAY_OLD_H

/*
 * <blitz/array.h> used to include most of the Blitz++ library
 * functionality, totally ~ 120000 lines of source code.  This
 * made for extremely slow compile times; processing #include <blitz/array.h>
 * took gcc about 25 seconds on a 500 MHz pentium box.
 *
 * Much of this compile time was due to the old vector expression templates
 * implementation.  Since this is not really needed for the Array<T,N>
 * class, the headers were redesigned so that:
 *
 * #include <blitz/array-old.h>   is the old-style include, pulls in most
 *                                of Blitz++ including vector e.t. 
 * #include <blitz/array.h>       pulls in much less of the library, and
 *                                in particular excludes the vector e.t. code
 *
 * With <blitz/array-old.h>, one gets TinyVector expressions automatically.
 * With <blitz/array.h>, one must now also include <blitz/tinyvec-et.h> 
 * to get TinyVector expressions.
 *
 * The implementation of Array<T,N> has been moved to <blitz/array-impl.h>.
 */

#include <blitz/tinyvec-et.h>
#include <blitz/array-impl.h>

#endif // BZ_ARRAY_OLD_H

