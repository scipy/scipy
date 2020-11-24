// Copyright (c) 2018-2019 ERGO-Code. See license.txt for license.

#include "diagonal_precond.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "timer.h"

namespace ipx {

DiagonalPrecond::DiagonalPrecond(const Model& model) : model_(model) {
    const Int m = model_.rows();
    diagonal_.resize(m);
}

void DiagonalPrecond::Factorize(const double* W, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const SparseMatrix& AI = model_.AI();

    factorized_ = false;

    // Build diagonal of normal matrix.
    if (W) {
        for (Int i = 0; i < m; i++)
            diagonal_[i] = W[n+i];
        for (Int j = 0; j < n; j++) {
            double w = W[j];
            for (Int p = AI.begin(j); p < AI.end(j); p++)
                diagonal_[AI.index(p)] += AI.value(p) * w * AI.value(p);
        }
    } else {
        diagonal_ = 0.0;        // rightmost m columns have weight zero
        for (Int j = 0; j < n; j++) {
            for (Int p = AI.begin(j); p < AI.end(j); p++)
                diagonal_[AI.index(p)] += AI.value(p) * AI.value(p);
        }
    }

    factorized_ = true;
}

double DiagonalPrecond::time() const {
    return time_;
}

void DiagonalPrecond::reset_time() {
    time_ = 0.0;
}

void DiagonalPrecond::_Apply(const Vector& rhs, Vector& lhs,
                              double* rhs_dot_lhs) {
    const Int m = model_.rows();
    double rldot = 0.0;
    Timer timer;

    assert(factorized_);
    assert((int)lhs.size() == m);
    assert((int)rhs.size() == m);

    for (Int i = 0; i < m; i++) {
        lhs[i] = rhs[i] / diagonal_[i];
        rldot += lhs[i] * rhs[i];
    }
    if (rhs_dot_lhs)
        *rhs_dot_lhs = rldot;
    time_ += timer.Elapsed();
}

}  // namespace ipx
