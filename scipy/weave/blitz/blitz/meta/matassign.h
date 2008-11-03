// -*- C++ -*-
/***************************************************************************
 * blitz/meta/matassign.h   TinyMatrix assignment metaprogram
 *
 * $Id: matassign.h 1413 2005-11-01 22:04:15Z cookedm $
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


#ifndef BZ_META_MATASSIGN_H
#define BZ_META_MATASSIGN_H

BZ_NAMESPACE(blitz)

template<int N_rows, int N_columns, int I, int J>
class _bz_meta_matAssign2 {
public:
    static const int go = (J < N_columns - 1) ? 1 : 0;

    template<typename T_matrix, typename T_expr, typename T_updater>
    static inline void f(T_matrix& mat, T_expr expr, T_updater u)
    {
        u.update(mat(I,J), expr(I,J));
        _bz_meta_matAssign2<N_rows * go, N_columns * go, I * go, (J+1) * go>
            ::f(mat, expr, u);
    }
};

template<>
class _bz_meta_matAssign2<0,0,0,0> {
public:
    template<typename T_matrix, typename T_expr, typename T_updater>
    static inline void f(T_matrix&, T_expr, T_updater)
    { }
};

template<int N_rows, int N_columns, int I> 
class _bz_meta_matAssign {
public:
    static const int go = (I < N_rows-1) ? 1 : 0;

    template<typename T_matrix, typename T_expr, typename T_updater>
    static inline void f(T_matrix& mat, T_expr expr, T_updater u)
    {
        _bz_meta_matAssign2<N_rows, N_columns, I, 0>::f(mat, expr, u);
        _bz_meta_matAssign<N_rows * go, N_columns * go, (I+1) * go>
            ::f(mat, expr, u);
    }
};

template<>
class _bz_meta_matAssign<0,0,0> {
public:
    template<typename T_matrix, typename T_expr, typename T_updater>
    static inline void f(T_matrix&, T_expr, T_updater)
    { }
};


BZ_NAMESPACE_END

#endif // BZ_META_ASSIGN_H
