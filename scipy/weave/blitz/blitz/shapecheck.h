/***************************************************************************
 * blitz/shapecheck.h    Functions for checking conformability of arrays
 *
 * $Id: shapecheck.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_SHAPECHECK_H
#define BZ_SHAPECHECK_H

BZ_NAMESPACE(blitz)

/*
 * The function areShapesConformable(A,B) checks that the shapes 
 * A and B are conformable (i.e. the same size/geometry).  Typically 
 * the A and B parameters are of type TinyVector<int,N_rank> and represent 
 * the extent of the arrays.  It's possible that in the future jagged-edged
 * arrays will be supported, in which case shapes may be lists
 * of subdomains.
 */

template<typename T_shape1, typename T_shape2>
inline bool areShapesConformable(const T_shape1&, const T_shape2&)
{
    // If the shape objects are different types, this means
    // that the arrays are different ranks, or one is jagged
    // edged, etc.  In this case the two arrays are not
    // conformable.
    return false;
}

template<typename T_shape>
inline bool areShapesConformable(const T_shape& a, const T_shape& b)
{
    // The shape objects are the same type, so compare them.

    // NEEDS_WORK-- once the "all" reduction is implemented, should
    // use it.
    // return all(a == b);

    for (unsigned i=0; i < a.length(); ++i)
    {
        if (a[i] != b[i])
        {
            BZ_DEBUG_MESSAGE("Incompatible shapes detected: " << endl 
                 << a << endl << b << endl);
            return false;
        }
    }

    return true;
}

BZ_NAMESPACE_END

#endif
