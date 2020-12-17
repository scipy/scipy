// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_LU_UPDATE_H_
#define IPX_LU_UPDATE_H_

#include <vector>
#include "indexed_vector.h"
#include "sparse_matrix.h"

namespace ipx {

// Interface class for LU factorization + update implementations.

class LuUpdate {
public:
    virtual ~LuUpdate() {}
    LuUpdate& operator=(const LuUpdate&) = delete;
    LuUpdate& operator=(LuUpdate&&) = delete;

    // Factorizes matrix given in compressed column format (4-array notation).
    // The matrix dimension (denoted dim in the documentation of the following
    // methods) is set by the constructor of the derived class. If during the
    // course of the factorization a column of the active submatrix becomes
    // (numerically) zero, it is replaced by a unit column in L*U.
    //
    // @strict_abs_pivottol: If true, then the implementation is asked to use
    //                       kLuDependencyTol as absolute pivot tolerance and to
    //                       remove columns from the active submatrix
    //                       immediately when all entries became smaller than
    //                       the abolute pivot tolerance. Need not be supported
    //                       by the implementation.
    //
    // Factorize() cannot fail other than for out of memory, in which case
    // std::bad_alloc is thrown. Returns
    // 0 OK
    // 1 OK, but the factorization is numerically unstable; suggests tightening
    //   the pivot tolerance (see below) and to refactorize.
    // 2 OK, but singularities occured and were replaced by unit columns.
    // 3 = 1 and 2
    Int Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                  const double* Bx, bool strict_abs_pivottol);

    // Exports LU factors. The method can only be called after a fresh
    // factorization. (Otherwise an assertion will fail.) It returns L, U and
    // permutations such that
    //
    //   B[rowperm,colperm] = (L+I)*U,
    //
    // where B is the matrix given to Factorize() after replacing dependent
    // columns by unit columns. The indices in L and U are sorted. Any of the
    // arguments can be NULL, in which case the quantity is not returned.
    //
    // @L: returns matrix L, will be resized as necessary.
    // @U: returns matrix U, will be resized as necessary.
    // @rowperm: size dim array, returns row permutation.
    // @colperm: size dim array, returns column permutation.
    // @dependent_cols: returns indices k for which U[:,k] was replaced by a
    //                  unit column in Factorize().
    void GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                    Int* colperm, std::vector<Int>* dependent_cols);

    // Solves linear system with dense right-hand side and solution.
    // @rhs, @lhs: size dim vectors, may refer to the same object.
    // @trans: 't' or 'T' for transposed system
    void SolveDense(const Vector& rhs, Vector& lhs, char trans);

    // Solves B*x=b in preparation for replacing a column of B by b.
    // @nz, @bi, @bx: b as compressed sparse vector.
    // @lhs: size dim vector returning x.
    // The 3-argument version can save the triangular solve with U.
    void FtranForUpdate(Int nz, const Int* bi, const double* bx);
    void FtranForUpdate(Int nz, const Int* bi, const double* bx,
                        IndexedVector& lhs);

    // Solves B'*y=ej in preparation for replacing column j of B, where ej
    // denotes the j-th unit vector.
    // @lhs: size dim vector returning y.
    // The 1-argument version can save the triangular solve with L.
    void BtranForUpdate(Int j);
    void BtranForUpdate(Int j, IndexedVector& lhs);

    // Updates factorization. The column to be replaced was defined in the last
    // call to BtranForUpdate(). The new column was given in the last call to
    // FtranForUpdate(). Both methods must have been called prior to Update().
    // (Otherwise an assertion will fail.)
    // @pivot: p-th entry of B^{-1}*b, if column p of B is to be replaced by b.
    // Returns: < 0 if the updated matrix is (numerically) singular. The update
    //              may or may not have been performed. The behaviour is
    //              undefined if any of the solve routines is called before
    //              computing a fresh factorization.
    //          > 0 if the updated factorization looks numerically unstable.
    //            0 otherwise
    Int Update(double pivot);

    // Returns true if refactorization is required or recommended for memory
    // or speed. No stability test is done.
    bool NeedFreshFactorization();

    // Returns (nnz(L)+nnz(U))/nnz(B) from the last factorization.
    double fill_factor() const;

    // Gets/sets the tolerance for partial threshold pivoting in Factorize().
    double pivottol() const;
    void pivottol(double new_pivottol);

    // Returns the number of updates since the last factorization.
    Int updates() const;

private:
    virtual Int _Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                           const double* Bx, bool strict_abs_pivottol) = 0;
    virtual void _GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                             Int* colperm, std::vector<Int>* dep_cols) = 0;
    virtual void _SolveDense(const Vector& rhs, Vector& lhs, char trans) = 0;
    virtual void _FtranForUpdate(Int nz, const Int* bi, const double* bx) = 0;
    virtual void _FtranForUpdate(Int nz, const Int* bi, const double* bx,
                                IndexedVector& lhs) = 0;
    virtual void _BtranForUpdate(Int p) = 0;
    virtual void _BtranForUpdate(Int p, IndexedVector& lhs) = 0;
    virtual Int _Update(double pivot) = 0;
    virtual bool _NeedFreshFactorization() = 0;
    virtual double _fill_factor() const = 0;
    virtual double _pivottol() const = 0;
    virtual void _pivottol(double new_pivottol) = 0;

    Int updates_{0};            // counts updates since factorization
};

}  // namespace ipx

#endif  // IPX_LU_UPDATE_H_
