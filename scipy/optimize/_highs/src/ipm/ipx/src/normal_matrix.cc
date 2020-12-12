// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "normal_matrix.h"
#include <cassert>
#include "timer.h"
#include "utils.h"

namespace ipx {

// Method for computing matrix-vector products:
//  1 one pass over A
//  2 two passes over A, first by cols, second by rows
//  3 two passes over A, first by cols, second by cols
// Method 2 and 3 are implemented only for W != NULL.
// For W == NULL method 1 is used always.
//
// The three variants of matrix-vector products were implemented and profiled
// during the development phase of IPX. It turned out that the one-pass variant
// is the fastest on average (about 20% better than the best two-pass variant),
// and also the fastest on most LP models. Therefore, it is used for
// matrix-vector products of the form AA' here and in SplittedNormalMatrix.
#define MATVECMETHOD 1

NormalMatrix::NormalMatrix(const Model& model) : model_(model) {
    #if MATVECMETHOD > 1
    // The two-pass variants require n+m workspace to store the intermediate
    // result W*AI'*rhs.
    work_.resize(model.rows() + model.cols());
    #endif
}

void NormalMatrix::Prepare(const double* W) {
    W_ = W;
    prepared_ = true;
}

double NormalMatrix::time() const {
    return time_;
}

void NormalMatrix::reset_time() {
    time_ = 0.0;
}

void NormalMatrix::_Apply(const Vector& rhs, Vector& lhs,
                           double* rhs_dot_lhs) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Int* Ap = model_.AI().colptr();
    const Int* Ai = model_.AI().rowidx();
    const double* Ax = model_.AI().values();
    #if MATVECMETHOD == 2
    const Int* Atp = model_.AIt().colptr();
    const Int* Ati = model_.AIt().rowidx();
    const double* Atx = model_.AIt().values();
    #endif
    Timer timer;

    assert(prepared_);
    assert((int)lhs.size() == m);
    assert((int)rhs.size() == m);

    if (W_) {
        #if MATVECMETHOD == 1
        for (Int i = 0; i < m; i++)
            lhs[i] = rhs[i] * W_[n+i];
        for (Int j = 0; j < n; j++) {
            Int begin = Ap[j], end = Ap[j+1];
            double d = 0.0;
            for (Int p = begin; p < end; p++)
                d += rhs[Ai[p]] * Ax[p];
            d *= W_[j];
            for (Int p = begin; p < end; p++)
                lhs[Ai[p]] += d * Ax[p];
        }
        #elif MATVECMETHOD == 2
        for (Int j = 0; j < n; j++) {
            Int begin = Ap[j], end = Ap[j+1];
            double d = 0.0;
            for (Int p = begin; p < end; p++)
                d += rhs[Ai[p]] * Ax[p];
            work_[j] = d * W_[j];
        }
        for (Int i = 0; i < m; i++)
            lhs[i] = rhs[i] * W_[n+i];
        for (Int i = 0; i < m; i++) {
            Int begin = Atp[i], end = Atp[i+1]-1; // skip identity entry
            double d = 0.0;
            for (Int p = begin; p < end; p++)
                d += work_[Ati[p]] * Atx[p];
            lhs[i] += d;
        }
        #elif MATVECMETHOD == 3
        for (Int j = 0; j < n; j++) {
            Int begin = Ap[j], end = Ap[j+1];
            double d = 0.0;
            for (Int p = begin; p < end; p++)
                d += rhs[Ai[p]] * Ax[p];
            work_[j] = d * W_[j];
        }
        for (Int i = 0; i < m; i++)
            lhs[i] = rhs[i] * W_[n+i];
        for (Int j = 0; j < n; j++) {
            Int begin = Ap[j], end = Ap[j+1];
            double d = work_[j];
            for (Int p = begin; p < end; p++)
                lhs[Ai[p]] += d * Ax[p];
        }
        #else
        #error "invalid MATHVECMETHOD"
        #endif
    } else {
        lhs = 0.0;
        for (Int j = 0; j < n; j++) {
            Int begin = Ap[j], end = Ap[j+1];
            double d = 0.0;
            for (Int p = begin; p < end; p++)
                d += rhs[Ai[p]] * Ax[p];
            for (Int p = begin; p < end; p++)
                lhs[Ai[p]] += d * Ax[p];
        }
    }
    if (rhs_dot_lhs)
        *rhs_dot_lhs = Dot(rhs,lhs);
    time_ += timer.Elapsed();
}

}  // namespace ipx
