/***************************************************************************
 * blitz/vector-et.h      Vector<P_numtype> class + expression templates
 *
 * $Id: vector-et.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_VECTOR_ET_H
#define BZ_VECTOR_ET_H

#include <blitz/vector.h>

// These are compile-time expensive things not included
// by <blitz/vector.h>, but needed if we want vector expressions.

#include <blitz/vecbops.cc>         // Operators with two operands
#include <blitz/vecuops.cc>         // Functions with one argument
#include <blitz/vecbfn.cc>          // Functions with two arguments

#endif  // BZ_VECTOR_ET_H

