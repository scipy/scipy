#ifndef IPX_INDEXED_VECTOR_H_
#define IPX_INDEXED_VECTOR_H_

#include <vector>
#include "ipm/ipx/ipx_internal.h"

namespace ipx {

// IndexedVector stores a vector and optionally the indices of its nonzero
// entries (its pattern). To iterate over the entries in an IndexedVector v:
//
//   double vmax = 0.0;
//   Int imax = 0;
//   auto update_max = [&vmax,&imax](Int i, double x) {
//       if (x > vmax) {
//           vmax = x;
//           imax = i;
//       }
//   };
//   for_each_nonzero(v, update_max);
//
// The loop iterates with indirect addressing over the nonzero entries if the
// pattern is known and sparse, and with direct addressing over the whole vector
// otherwise.
//
// When modifying the vector changes its pattern (e.g. by writing to v[i] for an
// arbitrary index i), you have to invalidate the pattern or provide the new one.

class IndexedVector {
public:
    // Constructs a vector of dimension @dim. Entries are initialized to zero
    // and pattern becomes empty.
    explicit IndexedVector(Int dim = 0);

    Int dim() const { return elements_.size(); }

    // Accesses entry 0 <= @i < dim() by value or reference.
    double operator[](Int i) const { return elements_[i]; }
    double& operator[](Int i) { return elements_[i]; }

    // Returns true if pattern is known and sparse.
    bool sparse() const;

    // Returns # nonzeros, or a negative value if the pattern is unknown.
    Int nnz() const { return nnz_; }

    // Invalidates the pattern.
    void InvalidatePattern() { nnz_ = -1; }

    // Changes the # indices in the pattern to new_nnz <= dim(). A negative
    // value invalidates the pattern.
    void set_nnz(Int new_nnz) { nnz_ = new_nnz; }

    // Returns the underlying array that stores the vector.
    const double* elements() const { return &elements_[0]; }
    double* elements() { return &elements_[0]; }

    // Returns the underlying array that stores the pattern.
    const Int* pattern() const { return pattern_.data(); }
    Int* pattern() { return pattern_.data(); }

    // Sets all entries to zero and the pattern to empty.
    void set_to_zero();

private:
    Vector elements_;
    std::vector<Int> pattern_;
    // If nnz_ >= 0, then pattern_[0..nnz_-1] are the indices of (possible)
    // nonzeros in elements_. If nnz_ < 0, then the pattern is unknown.
    Int nnz_{0};
};

template <typename C>
void for_each_nonzero(IndexedVector& v, C& c) {
    if (v.sparse()) {
        const Int* pattern = v.pattern();
        Int nnz = v.nnz();
        for (Int p = 0; p < nnz; p++) {
            const Int i = pattern[p];
            c(i, v[i]);
        }
    } else {
        Int dim = v.dim();
        for (Int i = 0; i < dim; i++) {
            const Int ii = i;   // make sure that caller does not change i
            c(ii, v[i]);
        }
    }
}

template <typename C>
void for_each_nonzero(const IndexedVector& v, C& c) {
    if (v.sparse()) {
        const Int* pattern = v.pattern();
        Int nnz = v.nnz();
        for (Int p = 0; p < nnz; p++) {
            const Int i = pattern[p];
            c(i, v[i]);
        }
    } else {
        Int dim = v.dim();
        for (Int i = 0; i < dim; i++) {
            const Int ii = i;   // make sure that caller does not change i
            c(ii, v[i]);
        }
    }
}

double Dot(const IndexedVector& x, const Vector& y);

}  // namespace ipx

#endif  // IPX_INDEXED_VECTOR_H_
