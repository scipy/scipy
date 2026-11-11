#ifndef IPX_SPARSE_MATRIX_H_
#define IPX_SPARSE_MATRIX_H_

#include <vector>
#include "ipm/ipx/ipx_internal.h"

namespace ipx {

// Sparse matrix in CSC format.

class SparseMatrix {
public:
    SparseMatrix();
    SparseMatrix(Int nrow, Int ncol);
    SparseMatrix(Int nrow, Int ncol, Int min_capacity);

    Int rows() const { return nrow_; }
    Int cols() const { return colptr_.size()-1; }
    Int entries() const { return colptr_.back(); }

    // # entries in column j
    Int entries(Int j) const { return end(j)-begin(j); }

    // Maximum # entries that can be stored in the matrix.
    Int capacity() const { return rowidx_.size(); }

    // Increases capacity if necessary such that capacity() >= min_capacity.
    // Matrix remains unchanged, pointers are invalidated.
    void reserve(Int min_capacity);

    // Changes matrix dimensions. Matrix becomes empty, pointers are
    // invalidated.
    void resize(Int nrow, Int ncol, Int min_capacity = 0);

    // Identical to resize(0,0).
    void clear();

    // Builds matrix from data in compressed column format. The matrix data must
    // be valid (no duplicates, no out of range indices, no inf/nan values);
    // this is not checked. The row indices in the input matrix need not be
    // sorted, but those in the output matrix will be.
    void LoadFromArrays(Int nrow, Int ncol, const Int* Abegin, const Int* Aend,
                        const Int* Ai, const double* Ax);

    Int begin(Int j) const { return colptr_[j]; }
    Int end(Int j) const { return colptr_[j+1]; }

    // Accesses entry at position @pos by value.
    Int index(Int pos) const { return rowidx_[pos]; }
    double value(Int pos) const { return values_[pos]; }

    // Accesses entry at position @pos by reference.
    Int& index(Int pos) { return rowidx_[pos]; }
    double& value(Int pos) { return values_[pos]; }

    // Accesses underlying arrays.
    const Int *colptr() const { return colptr_.data(); }
    const Int *rowidx() const { return rowidx_.data(); }
    const double *values() const { return values_.data(); }
    Int *colptr() { return colptr_.data(); }
    Int *rowidx() { return rowidx_.data(); }
    double *values() { return values_.data(); }

    // Stores the entries in each column in increasing order of index.
    void SortIndices();

    // The following methods provide a queue for adding new columns to the
    // matrix. Entries in the queue are not part of the matrix (so do not
    // contribute to entries()).

    // Appends an entry to the end of the queue.
    void push_back(Int i, double x) {
        rowidx_queue_.push_back(i);
        values_queue_.push_back(x);
    }

    // Returns # entries in the queue.
    Int queue_size() const { return rowidx_queue_.size(); }

    // Accesses entry at position @pos in the queue by value.
    Int qindex(Int pos) const { return rowidx_queue_[pos]; }
    double qvalue(Int pos) const { return values_queue_[pos]; }

    // Accesses entry at position @pos in the queue by reference.
    Int& qindex(Int pos) { return rowidx_queue_[pos]; }
    double& qvalue(Int pos) { return values_queue_[pos]; }

    // Makes new column from queue. The queue becomes empty and cols()
    // increases by 1.
    void add_column();

    // Discards queue.
    void clear_queue();

private:
    // Returns true if row indices are sorted.
    bool IsSorted() const;

    Int nrow_;
    std::vector<Int> colptr_;
    std::vector<Int> rowidx_;
    std::vector<double> values_;
    std::vector<Int> rowidx_queue_;
    std::vector<double> values_queue_;
};

// Builds transpose of matrix.
SparseMatrix Transpose(const SparseMatrix& A);

// Resizes @AT as necessary and fills with the tranpose of A.
void Transpose(const SparseMatrix& A, SparseMatrix& AT);

// Returns a copy of A[:,cols].
SparseMatrix CopyColumns(const SparseMatrix& A, const std::vector<Int>& cols);

// Permutes rows in place so that row i becomes row perm[i].
void PermuteRows(SparseMatrix& A, const std::vector<Int>& perm);

// Multiplies column j by s.
inline void ScaleColumn(SparseMatrix& A, Int j, double s) {
    Int p1 = A.begin(j);
    Int p2 = A.end(j);
    for (Int p = p1; p < p2; p++)
        A.value(p) *= s;
}

// Removes diagonal entries from A. If @diag is not NULL, then it must be an
// array of dimension A.cols() that holds the diagonal of A on return. Diagonal
// entries that are not present in A are set to zero in @diag. Returns the #
// entries removed from A.
Int RemoveDiagonal(SparseMatrix& A, double* diag);

// Returns dot(A[:,j], rhs).
inline double DotColumn(const SparseMatrix& A, Int j, const Vector& rhs) {
    Int p1 = A.begin(j);
    Int p2 = A.end(j);
    double d = 0.0;
    for (Int p = p1; p < p2; p++)
        d += rhs[A.index(p)] * A.value(p);
    return d;
}

// Updates lhs := lhs + alpha * A[:,j].
inline void ScatterColumn(const SparseMatrix& A, Int j, double alpha,
                          Vector& lhs) {
    Int p1 = A.begin(j);
    Int p2 = A.end(j);
    for (Int p = p1; p < p2; p++)
        lhs[A.index(p)] += alpha * A.value(p);
}

// Updates lhs := lhs + alpha*A*rhs or lhs := lhs + alpha*A'*rhs.
// @trans: 't' or 'T' for transposed product.
void MultiplyAdd(const SparseMatrix& A, const Vector& rhs, double alpha,
                 Vector& lhs, char trans);

// Updates lhs := lhs + A*A'*rhs or lhs := lhs + A*D*D*A'*rhs,
// where D is diagonal matrix if @D != NULL.
void AddNormalProduct(const SparseMatrix& A, const double* D, const Vector& rhs,
                      Vector& lhs);

// Triangular solve with sparse matrix.
// @x: right-hand side on entry, left-hand side on return.
// @trans: 't' or 'T' for transposed system.
// @uplo: must have *uplo == 'u' or *uplo == 'U' if A is upper triangular;
//        otherwise A is lower triangular.
// @unitdiag: nonzero if A has a unit diagonal that is not stored.
// Returns the # nonzeros in the solution.
Int TriangularSolve(const SparseMatrix& A, Vector& x, char trans,
                    const char* uplo, int unitdiag);

// Solves (L*U) x = x.
// L unit lower triangular, stored without diagonal.
// U upper triangular, diagonal element at end of column.
void ForwardSolve(const SparseMatrix& L, const SparseMatrix& U, Vector& x);

// Solves (L*U)' x = x.
// L unit lower triangular, stored without diagonal.
// U upper triangular, diagonal element at end of column.
void BackwardSolve(const SparseMatrix& L, const SparseMatrix& U, Vector& x);

// Returns the 1-norm and infinity-norm of A.
double Onenorm(const SparseMatrix& A);
double Infnorm(const SparseMatrix& A);

// Estimates the 1-norm of inverse(A).
// @A must be square and lower or upper triangular.
// @uplo: must have *uplo == 'u' or *uplo == 'U' if A is upper triangular;
//        otherwise A is lower triangular.
// @unitdiag: nonzero if A has a unit diagonal that is not stored.
double NormestInverse(const SparseMatrix& A, const char* uplo, int unitdiag);

}  // namespace ipx

#endif  // IPX_SPARSE_MATRIX_H_
