/***************************************************************************
 * blitz/random.h       Random number generator wrapper class
 *
 * $Id: random.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_RANDOM_H
#define BZ_RANDOM_H

#ifndef BZ_BLITZ_H
 #include <blitz/blitz.h>
#endif

BZ_NAMESPACE(blitz)

template<typename P_distribution>
class Random {

public:
    typedef P_distribution T_distribution;
    typedef _bz_typename T_distribution::T_numtype T_numtype;

    Random(double parm1=0.0, double parm2=1.0, double parm3=0.0)
        : generator_(parm1, parm2, parm3)
    { }

    void randomize()
    { generator_.randomize(); }
   
    T_numtype random()
    { return generator_.random(); }

    operator T_numtype()
    { return generator_.random(); }

protected: 
    T_distribution generator_;
};

BZ_NAMESPACE_END

#include <blitz/randref.h>

#endif // BZ_RANDOM_H

