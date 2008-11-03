// -*- C++ -*-
/***************************************************************************
 * blitz/tinymatexpr.h   Tiny Matrix Expressions
 *
 * $Id: tinymatexpr.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_TINYMATEXPR_H
#define BZ_TINYMATEXPR_H

#ifndef BZ_TINYMAT_H
 #error <blitz/tinymatexpr.h> must be included via <blitz/tinymat.h>
#endif

BZ_NAMESPACE(blitz)

template<typename T_expr>
class _bz_tinyMatExpr {
public:
    typedef _bz_typename T_expr::T_numtype T_numtype;

    static const int
        rows = T_expr::rows,
        columns = T_expr::columns;

    _bz_tinyMatExpr(T_expr expr)
        : expr_(expr)
    { }

    _bz_tinyMatExpr(const _bz_tinyMatExpr<T_expr>& x)
        : expr_(x.expr_)
    { }

    T_numtype operator()(int i, int j) const
    { return expr_(i,j); }

protected:
    T_expr expr_;
};

BZ_NAMESPACE_END

#endif // BZ_TINYMATEXPR_H

