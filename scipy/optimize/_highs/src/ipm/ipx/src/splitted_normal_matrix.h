// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_SPLITTED_NORMAL_MATRIX_H_
#define IPX_SPLITTED_NORMAL_MATRIX_H_

#include <vector>
#include "basis.h"
#include "linear_operator.h"
#include "model.h"
#include "sparse_matrix.h"

namespace ipx {

// SplittedNormalMatrix provides matrix-vector products with
//
//   C = inv(B)*AI*AI'*inv(B') = I + inv(B)*N*N'*inv(B'),
//
// where AI is the m-by-(n+m) matrix defined by the model and [B N] is the
// partitioning of AI into basic and nonbasic columns defined by the basis.
// The columns of B and N are scaled by the interior point scaling factors
// provided in the call to Prepare().
//
// When a variable has status BASIC_FREE, the row and column of C become a unit
// vector. When a variable has status NONBASIC_FIXED, it is dropped from N.

class SplittedNormalMatrix : public LinearOperator {
public:
    // Constructor stores a reference to the model. No data is copied. The model
    // must be valid as long as the object is used.
    explicit SplittedNormalMatrix(const Model& model);

    // Prepares object for subsequent calls to Apply(). @colscale must hold n+m
    // scaling factors for the columns of AI. The scaling factors are copied.
    void Prepare(const Basis& basis, const double* colscale);

    // Returns the column permutation from the LU factorization of the basis
    // matrix. The permutation was stored in the object by Prepare().
    const Int* colperm() const;

    // Returns computation times for operations since last call to reset_time().
    double time_B() const;
    double time_Bt() const;
    double time_NNt() const;
    void reset_time();

private:
    void _Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs) override;

    const Model& model_;
    SparseMatrix L_;           // lower triangular factor without unit diagonal
    SparseMatrix U_;           // upper triangular factor with scaled columns
    SparseMatrix N_;           // N with scaled columns and permuted row indices
    std::vector<Int> free_positions_; // positions corresponding to free vars
    std::vector<Int> colperm_;        // column permutation from LU factor
    std::vector<Int> rowperm_inv_;    // inverse row permutation from LU factor
    Vector work_;                     // size m workspace
    bool prepared_{false};            // operator prepared?
    double time_B_{0.0};              // time solves with B
    double time_Bt_{0.0};             // time solves with B'
    double time_NNt_{0.0};            // time matrix-vector products with NN'
};

}  // namespace ipx

#endif  // IPX_SPLITTED_NORMAL_MATRIX_H_
