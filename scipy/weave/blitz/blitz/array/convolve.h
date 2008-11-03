/***************************************************************************
 * blitz/array/convolve.h   One-dimensional convolution
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
#ifndef BZ_ARRAY_CONVOLVE_H
#define BZ_ARRAY_CONVOLVE_H

#ifndef BZ_ARRAY_H
 #error <blitz/array/convolve.h> must be included after <blitz/array.h>
#endif

BZ_NAMESPACE(blitz)

template<typename T>
Array<T,1> convolve(const Array<T,1>& B, const Array<T,1>& C);

BZ_NAMESPACE_END

#include <blitz/array/convolve.cc>

#endif // BZ_ARRAY_CONVOLVE_H
