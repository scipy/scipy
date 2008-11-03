/***************************************************************************
 * blitz/randref.h      Random number generators, expression templates
 *                      wrapper
 *
 * $Id: randref.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_RANDREF_H
#define BZ_RANDREF_H

#ifndef BZ_RANDOM_H
 #error <blitz/randref.h> must be included via <blitz/random.h>
#endif // BZ_RANDOM_H

BZ_NAMESPACE(blitz)

template<typename P_distribution>
class _bz_VecExprRandom {

public:
    typedef _bz_typename Random<P_distribution>::T_numtype T_numtype;

    _bz_VecExprRandom(Random<P_distribution>& random)
        : random_(random)
    { }

#ifdef BZ_MANUAL_VECEXPR_COPY_CONSTRUCTOR
    _bz_VecExprRandom(_bz_VecExprRandom<P_distribution>& x)
        : random_(x.random_)
    { }
#endif

    T_numtype operator[](unsigned) const
    { return random_.random(); }

    T_numtype operator()(unsigned) const
    { return random_.random(); }

    unsigned length(unsigned recommendedLength) const
    { return recommendedLength; }

    unsigned _bz_suggestLength() const
    { return 0; }

    bool _bz_hasFastAccess() const
    { return 1; }

    T_numtype _bz_fastAccess(unsigned) const
    { return random_.random(); }

private:
    _bz_VecExprRandom() : random_( Random<P_distribution>() ) { }

    Random<P_distribution>& random_;
};

BZ_NAMESPACE_END

#endif // BZ_RANDREF_H

