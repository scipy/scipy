#include "ipm/ipx/sparse_matrix.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include "ipm/ipx/utils.h"
#include "../extern/pdqsort/pdqsort.h"

namespace ipx {

SparseMatrix::SparseMatrix() {
    resize(0,0);
}

SparseMatrix::SparseMatrix(Int nrow, Int ncol) {
    resize(nrow, ncol);
}

SparseMatrix::SparseMatrix(Int nrow, Int ncol, Int min_capacity) {
    resize(nrow, ncol, min_capacity);
}

void SparseMatrix::reserve(Int min_capacity) {
    if (min_capacity > capacity()) {
        rowidx_.resize(min_capacity);
        values_.resize(min_capacity);
    }
}

void SparseMatrix::resize(Int nrow, Int ncol, Int min_capacity) {
    assert(nrow >= 0);
    assert(ncol >= 0);
    assert(min_capacity >= 0);
    nrow_ = nrow;
    colptr_.resize(ncol+1);
    colptr_.shrink_to_fit();
    std::fill(colptr_.begin(), colptr_.end(), 0);
    rowidx_.resize(min_capacity);
    rowidx_.shrink_to_fit();
    values_.resize(min_capacity);
    values_.shrink_to_fit();
}

void SparseMatrix::clear() {
    resize(0,0);
}

void SparseMatrix::LoadFromArrays(Int nrow, Int ncol, const Int* Abegin,
                                  const Int* Aend, const Int* Ai,
                                  const double* Ax) {
    Int nz = 0;
    for (Int j = 0; j < ncol; j++)
        nz += Aend[j]-Abegin[j];
    resize(nrow, ncol, nz);
    Int put = 0;
    for (Int j = 0; j < ncol; j++) {
        colptr_[j] = put;
        for (Int p = Abegin[j]; p < Aend[j]; p++) {
            if (Ax[p] != 0.0) {
                rowidx_[put] = Ai[p];
                values_[put] = Ax[p];
                put++;
            }
        }
    }
    colptr_[ncol] = put;
    SortIndices();
}

void SparseMatrix::SortIndices() {
    if (IsSorted())
        return;
    std::vector<std::pair<Int,double>> work(rows());
    for (Int j = 0; j < cols(); j++) {
        Int nz = 0;             // # entries in column j
        for (Int p = begin(j); p < end(j); p++) {
            work[nz].first = index(p);
            work[nz].second = value(p);
            nz++;
        }
        pdqsort(work.begin(), work.begin() + nz);
        for (Int k = 0, p = begin(j); p < end(j); k++, p++) {
            index(p) = work[k].first;
            value(p) = work[k].second;
        }
    }
}

void SparseMatrix::add_column() {
    Int nz = entries();
    Int nznew = nz + queue_size();
    reserve(nznew);
    std::copy(rowidx_queue_.begin(), rowidx_queue_.end(), rowidx_.begin() + nz);
    std::copy(values_queue_.begin(), values_queue_.end(), values_.begin() + nz);
    colptr_.push_back(nznew);
    clear_queue();
}

void SparseMatrix::clear_queue() {
    rowidx_queue_.clear();
    values_queue_.clear();
}

bool SparseMatrix::IsSorted() const {
    for (Int j = 0; j < cols(); j++) {
        for (Int p = begin(j); p < end(j)-1; p++)
            if (index(p) > index(p+1))
                return false;
    }
    return true;
}

SparseMatrix Transpose(const SparseMatrix& A) {
    SparseMatrix AT;
    Transpose(A, AT);
    return AT;
}

void Transpose(const SparseMatrix& A, SparseMatrix& AT) {
    const Int m = A.rows();
    const Int n = A.cols();
    const Int nz = A.entries();
    AT.resize(n, m, nz);

    // Compute row counts of A in workspace.
    std::vector<Int> work(m);
    for (Int p = 0; p < nz; p++)
        work[A.index(p)]++;

    // Set column pointers for AT.
    Int* ATp = AT.colptr();
    Int sum = 0;
    for (Int i = 0; i < m; i++) {
        ATp[i] = sum;
        sum += work[i];
        work[i] = ATp[i];
    }
    assert(sum == nz);
    ATp[m] = sum;

    // Fill AT with one column of A at a time.
    // work[i] is the next free slot in column i of AT.
    for (Int j = 0; j < n; j++) {
        for (Int p = A.begin(j); p < A.end(j); p++) {
            Int put = work[A.index(p)]++;
            AT.index(put) = j;
            AT.value(put) = A.value(p);
        }
    }
}

SparseMatrix CopyColumns(const SparseMatrix& A, const std::vector<Int>& cols) {
    SparseMatrix A2(A.rows(), 0);
    for (Int j : cols) {
        for (Int p = A.begin(j); p < A.end(j); p++)
            A2.push_back(A.index(p), A.value(p));
        A2.add_column();
    }
    return A2;
}

void PermuteRows(SparseMatrix& A, const std::vector<Int>& perm) {
    for (Int p = 0; p < A.entries(); p++)
        A.index(p) = perm[A.index(p)];
}

Int RemoveDiagonal(SparseMatrix& A, double* diag) {
    Int ncol = A.cols();
    Int* Ap = A.colptr();
    Int* Ai = A.rowidx();
    double* Ax = A.values();
    Int get = 0;
    Int put = 0;
    for (Int j = 0; j < ncol; j++) {
        if (diag)
            diag[j] = 0.0;      // if no diagonal entry in column j
        Ap[j] = put;
        for (; get < Ap[j+1]; get++) {
            if (Ai[get] == j) {
                if (diag)
                    diag[j] = Ax[get];
            } else {
                Ai[put] = Ai[get];
                Ax[put] = Ax[get];
                put++;
            }
        }
    }
    Ap[ncol] = put;
    return get-put;
}

void MultiplyAdd(const SparseMatrix& A, const Vector& rhs, double alpha,
                 Vector& lhs, char trans) {
    const Int m = A.rows();
    const Int n = A.cols();
    if (trans == 't' || trans == 'T') {
        assert((Int)rhs.size() == m);
        assert((Int)lhs.size() == n);
        for (Int j = 0; j < n; j++)
            lhs[j] += alpha * DotColumn(A, j, rhs);
    } else {
        assert((Int)rhs.size() == n);
        assert((Int)lhs.size() == m);
        for (Int j = 0; j < n; j++)
            ScatterColumn(A, j, alpha*rhs[j], lhs);
    }
}

void AddNormalProduct(const SparseMatrix& A, const double* D, const Vector& rhs,
                      Vector& lhs) {
    const Int m = A.rows();
    const Int n = A.cols();
    assert((Int)rhs.size() == m);
    assert((Int)lhs.size() == m);
    for (Int j = 0; j < n; j++) {
        double temp = DotColumn(A, j, rhs);
        if (D) temp *= D[j]*D[j];
        ScatterColumn(A, j, temp, lhs);
    }
}

Int TriangularSolve(const SparseMatrix& A, Vector& x, char trans,
                    const char* uplo, int unitdiag) {
    const Int ncol = A.cols();
    const Int* Ap = A.colptr();
    const Int* Ai = A.rowidx();
    const double* Ax = A.values();
    Int nz = 0;
    if (trans == 't' || trans == 'T') {
        if (*uplo == 'u' || *uplo == 'U') {
            // transposed solve with upper triangular matrix
            for (Int i = 0; i < ncol; i++) {
                Int begin = Ap[i];
                Int end = Ap[i+1] - (unitdiag ? 0 : 1);
                double d = 0.0;
                for (Int p = begin; p < end; p++)
                    d += x[Ai[p]] * Ax[p];
                x[i] -= d;
                if (!unitdiag) {
                    assert(Ai[end] == i);
                    x[i] /= Ax[end];
                }
                if (x[i] != 0.0)
                    nz++;
            }
        } else {
            // transposed solve with lower triangular matrix
            for (Int i = ncol-1; i >= 0; i--) {
                Int begin = Ap[i] + (unitdiag ? 0 : 1);
                Int end = Ap[i+1];
                double d = 0.0;
                for (Int p = begin; p < end; p++)
                    d += x[Ai[p]] * Ax[p];
                x[i] -= d;
                if (!unitdiag) {
                    assert(Ai[begin-1] == i);
                    x[i] /= Ax[begin-1];
                }
                if (x[i] != 0.0)
                    nz++;
            }
        }
    } else {
        if (*uplo == 'u' || *uplo == 'U') {
            // forward solve with upper triangular matrix
            for (Int j = ncol-1; j >= 0; j--) {
                Int begin = Ap[j];
                Int end = Ap[j+1] - (unitdiag ? 0 : 1);
                if (!unitdiag) {
                    assert(Ai[end] == j);
                    x[j] /= Ax[end];
                }
                double temp = x[j];
                if (temp != 0.0) {
                    for (Int p = begin; p < end; p++)
                        x[Ai[p]] -= Ax[p] * temp;
                    nz++;
                }
            }
        } else {
            // forward solve with lower triangular matrix
            for (Int j = 0; j < ncol; j++) {
                Int begin = Ap[j] + (unitdiag ? 0 : 1);
                Int end = Ap[j+1];
                if (!unitdiag) {
                    assert(Ai[begin-1] == j);
                    x[j] /= Ax[begin-1];
                }
                double temp = x[j];
                if (temp != 0.0) {
                    for (Int p = begin; p < end; p++)
                        x[Ai[p]] -= Ax[p] * temp;
                    nz++;
                }
            }
        }
    }
    return nz;
}

void ForwardSolve(const SparseMatrix& L, const SparseMatrix& U, Vector& x) {
    TriangularSolve(L, x, 'n', "lower", 1);
    TriangularSolve(U, x, 'n', "upper", 0);
}

void BackwardSolve(const SparseMatrix& L, const SparseMatrix& U, Vector& x) {
    TriangularSolve(U, x, 't', "upper", 0);
    TriangularSolve(L, x, 't', "lower", 1);
}

double Onenorm(const SparseMatrix& A) {
    double norm = 0.0;
    for (Int j = 0; j < A.cols(); j++) {
        double colsum = 0.0;
        for (Int p = A.begin(j); p < A.end(j); p++)
            colsum += std::abs(A.value(p));
        norm = std::max(norm, colsum);
    }
    return norm;
}

double Infnorm(const SparseMatrix& A) {
    Vector rowsum(A.rows());
    for (Int j = 0; j < A.cols(); j++) {
        for (Int p = A.begin(j); p < A.end(j); p++)
            rowsum[A.index(p)] += std::abs(A.value(p));
    }
    return Infnorm(rowsum);
}

double NormestInverse(const SparseMatrix& A, const char* uplo, int unitdiag) {
    const Int m = A.rows();
    Vector x(m);
    assert(A.rows() == A.cols());

    // Solve A'x=b, where the entries of b are +/-1 chosen dynamically to make x
    // large.
    if (*uplo == 'u' || *uplo == 'U') {
        for (Int j = 0; j < m; j++) {
            Int begin = A.begin(j);
            Int end = A.end(j);
            if (!unitdiag)
                end--;
            double temp = 0.0;
            for (Int p = begin; p < end; p++)
                temp -= x[A.index(p)] * A.value(p);
            temp += temp >= 0.0 ? 1.0 : -1.0; // choose b[j] = 1 or b[j] = -1
            if (!unitdiag) {
                assert(A.index(end) == j);
                temp /= A.value(end);
            }
            x[j] = temp;
        }
    } else {
        for (Int j = m-1; j >= 0; j--) {
            Int begin = A.begin(j);
            Int end = A.end(j);
            if (!unitdiag)
                begin++;
            double temp = 0.0;
            for (Int p = begin; p < end; p++)
                temp -= x[A.index(p)] * A.value(p);
            temp += temp >= 0.0 ? 1.0 : -1.0; // choose b[j] = 1 or b[j] = -1
            if (!unitdiag) {
                assert(A.index(begin-1) == j);
                temp /= A.value(begin-1);
            }
            x[j] = temp;
        }
    }
    double x1norm = Onenorm(x);
    double xinfnorm = Infnorm(x);

    // Solve Ay=x, solution overwrites x.
    TriangularSolve(A, x, 'n', uplo, unitdiag);
    double y1norm = Onenorm(x);

    return std::max(y1norm/x1norm, xinfnorm);
}

} // namespace ipx
