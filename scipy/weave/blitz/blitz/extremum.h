/***************************************************************************
 * blitz/extremum.h      Declaration of the Extremum<T_numtype, T_index> class
 *
 * $Id: extremum.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_EXTREMUM_H
#define BZ_EXTREMUM_H

#ifndef BZ_BLITZ_H
 #include <blitz/blitz.h>
#endif

BZ_NAMESPACE(blitz)

// The Extremum class is used for returning extreme values and their
// locations in a numeric container.  It's a simple 2-tuple, with the
// first element being the extreme value, and the send its location.
// An object of type Extremum can be automatically converted to
// the numeric type via operator T_numtype().
template<typename P_numtype, typename P_index>
class Extremum {
public:
    typedef P_numtype T_numtype;
    typedef P_index   T_index;

    Extremum(T_numtype value, T_index index)
        : value_(value), index_(index)
    { }

    T_numtype value() const
    { return value_; }

    T_index index() const
    { return index_; }

    void setValue(T_numtype value)
    { value_ = value; }

    void setIndex(T_index index)
    { index_ = index; }

    operator T_numtype() const
    { return value_; }

protected:
    T_numtype value_;
    T_index index_;
};

BZ_NAMESPACE_END

#endif // BZ_EXTREMUM_H

