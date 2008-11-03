/***************************************************************************
 * blitz/tvcross.h      Cross product of TinyVector<N,3>'s
 *
 * $Id: tvcross.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_TVCROSS_H
#define BZ_TVCROSS_H

#ifndef BZ_TINYVEC_H
 #error <blitz/tvcross.h> must be included via <blitz/tinyvec.h>
#endif

BZ_NAMESPACE(blitz)

/*
 * cross product.
 *
 * NEEDS_WORK: - cross product of two different vector types
 *             - cross product involving expressions
 */

template<typename T_numtype>
TinyVector<T_numtype,3> cross(const TinyVector<T_numtype,3>& x, 
    const TinyVector<T_numtype,3>& y)
{
    return TinyVector<T_numtype,3>(x[1]*y[2] - y[1]*x[2],
        y[0]*x[2] - x[0]*y[2], x[0]*y[1] - y[0]*x[1]);
}


BZ_NAMESPACE_END

#endif // BZ_TVCROSS_H
