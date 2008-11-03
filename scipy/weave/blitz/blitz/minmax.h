/***************************************************************************
 * blitz/minmax.h  Declaration of min and max functions
 *
 * $Id: minmax.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_MINMAX_H
#define BZ_MINMAX_H

#include <blitz/promote.h>

BZ_NAMESPACE(blitz)

/*
 * These functions are in their own namespace (blitz::minmax) to avoid
 * conflicts with the array reduction operations min and max.
 */

BZ_NAMESPACE(minmax)

template<typename T1, typename T2>
BZ_PROMOTE(T1,T2) min(const T1& a, const T2& b)
{
    typedef BZ_PROMOTE(T1,T2) T_promote;

    if (a <= b)
        return T_promote(a);
    else
        return T_promote(b);
}

template<typename T1, typename T2>
BZ_PROMOTE(T1,T2) max(const T1& a, const T2& b)
{
    typedef BZ_PROMOTE(T1,T2) T_promote;

    if (a >= b)
        return T_promote(a);
    else
        return T_promote(b);
}

BZ_NAMESPACE_END

BZ_NAMESPACE_END

#endif
