/***************************************************************************
 * blitz/array/et.h  Include expression templates implementation for arrays
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
 ****************************************************************************/
#ifndef BZ_ARRAY_ET_H
#define BZ_ARRAY_ET_H

#ifdef BZ_NEW_EXPRESSION_TEMPLATES
 #include <blitz/array/newet.h>     // Expression templates
#else
 #include <blitz/array/bops.cc>     // Expression templates, two operands
 #include <blitz/array/uops.cc>     // Expression templates, math functions
 #include <blitz/array/misc.cc>     // Expression templates, miscellaneous
#endif

#endif
