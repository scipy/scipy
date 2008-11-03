// -*- C++ -*-
/***************************************************************************
 * blitz/etbase.h    Declaration of the ETBase<T> class
 *
 * $Id: etbase.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_ETBASE_H
#define BZ_ETBASE_H

BZ_NAMESPACE(blitz)

template<typename T>
class ETBase { 
public:
    ETBase() 
    { }

    ETBase(const ETBase<T>&)
    { }
    
    T& unwrap() { return static_cast<T&>(*this); }
    
    const T& unwrap() const { return static_cast<const T&>(*this); }
};

BZ_NAMESPACE_END

#endif // BZ_ETBASE_H

