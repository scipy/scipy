// -*- C++ -*-
/***************************************************************************
 * blitz/meta/metaprog.h   Useful metaprogram declarations
 *
 * $Id: metaprog.h 1413 2005-11-01 22:04:15Z cookedm $
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

#ifndef BZ_META_METAPROG_H
#define BZ_META_METAPROG_H

BZ_NAMESPACE(blitz)

// Null Operand

class _bz_meta_nullOperand {
public:
    _bz_meta_nullOperand() { }
};

template<typename T> inline T operator+(const T& a, _bz_meta_nullOperand)
{ return a; }
template<typename T> inline T operator*(const T& a, _bz_meta_nullOperand)
{ return a; }

// MetaMax

template<int N1, int N2>
class _bz_meta_max {
public:
    static const int max = (N1 > N2) ? N1 : N2;
};

// MetaMin

template<int N1, int N2>
class _bz_meta_min {
public:
    static const int min = (N1 < N2) ? N1 : N2;
};

BZ_NAMESPACE_END 

#endif // BZ_META_METAPROG_H
