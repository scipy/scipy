#include "ipm/ipx/indexed_vector.h"

namespace ipx {

IndexedVector::IndexedVector(Int dim) : elements_(dim), pattern_(dim) { }

bool IndexedVector::sparse() const {
    return nnz() >= 0 && nnz() <= kHypersparseThreshold * dim();
}

void IndexedVector::set_to_zero() {
    if (sparse()) {
        for (Int p = 0; p < nnz_; p++)
            elements_[pattern_[p]] = 0.0;
    } else {
        elements_ = 0.0;
    }
    nnz_ = 0;
}

double Dot(const IndexedVector& x, const Vector& y) {
    double d = 0.0;
    auto add = [&](Int p, double f) {
        d += y[p] * f;
    };
    for_each_nonzero(x, add);
    return d;
}

}  // namespace ipx
