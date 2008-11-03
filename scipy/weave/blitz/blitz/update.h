// -*- C++ -*-
/***************************************************************************
 * blitz/update.h      Declaration of the _bz_XXXX updater classes
 *
 * $Id: update.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_UPDATE_H
#define BZ_UPDATE_H

#include <blitz/blitz.h>

BZ_NAMESPACE(blitz)

class _bz_updater_base { };

#define BZ_DECL_UPDATER(name,op,symbol)                     \
  template<typename X, typename Y>                          \
  class name : public _bz_updater_base {                    \
  public:                                                   \
    static inline void update(X& restrict x, Y y)           \
    { x op y; }                                             \
    static void prettyPrint(BZ_STD_SCOPE(string) &str)      \
    { str += symbol; }                                      \
  }

template<typename X, typename Y>
class _bz_update : public _bz_updater_base {
  public:
    static inline void update(X& restrict x, Y y)
    { x = (X)y; }

    static void prettyPrint(BZ_STD_SCOPE(string) &str)
    { str += "="; }
};

BZ_DECL_UPDATER(_bz_plus_update, +=, "+=");
BZ_DECL_UPDATER(_bz_minus_update, -=, "-=");
BZ_DECL_UPDATER(_bz_multiply_update, *=, "*=");
BZ_DECL_UPDATER(_bz_divide_update, /=, "/=");
BZ_DECL_UPDATER(_bz_mod_update, %=, "%=");
BZ_DECL_UPDATER(_bz_xor_update, ^=, "^=");
BZ_DECL_UPDATER(_bz_bitand_update, &=, "&=");
BZ_DECL_UPDATER(_bz_bitor_update, |=, "|=");
BZ_DECL_UPDATER(_bz_shiftl_update, <<=, "<<=");
BZ_DECL_UPDATER(_bz_shiftr_update, >>=, ">>=");

BZ_NAMESPACE_END

#endif // BZ_UPDATE_H

