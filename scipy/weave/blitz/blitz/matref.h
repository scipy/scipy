/***************************************************************************
 * blitz/matref.h      Declaration of the _bz_MatrixRef<P_numtype, P_structure>
 *                     class.
 *
 * $Id: matref.h 1414 2005-11-01 22:04:59Z cookedm $
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

#ifndef BZ_MATREF_H
#define BZ_MATREF_H

#ifndef BZ_MATEXPR_H
 #error <blitz/matref.h> must be included via <blitz/matexpr.h>
#endif // BZ_MATEXPR_H

BZ_NAMESPACE(blitz)

template<typename P_numtype, typename P_structure>
class _bz_MatrixRef {

public:
    typedef P_numtype T_numtype;

    _bz_MatrixRef(const Matrix<P_numtype, P_structure>& m)
        : matrix_(&m)
    { }

    T_numtype operator()(unsigned i, unsigned j) const
    { return (*matrix_)(i,j); }

    unsigned rows(unsigned) const
    { return matrix_->rows(); }

    unsigned cols(unsigned) const
    { return matrix_->cols(); }

private:
    _bz_MatrixRef() { } 

    const Matrix<P_numtype, P_structure>* matrix_;
};

BZ_NAMESPACE_END

#endif // BZ_MATREF_H
