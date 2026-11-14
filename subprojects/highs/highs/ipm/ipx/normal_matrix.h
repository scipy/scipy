#ifndef IPX_NORMAL_MATRIX_H_
#define IPX_NORMAL_MATRIX_H_

#include "ipm/ipx/linear_operator.h"
#include "ipm/ipx/model.h"

namespace ipx {

// NormalMatrix provides matrix-vector operations with the matrix
//
//   AI*W*AI',
//
// where AI is the m-by-(n+m) matrix defined by the model, and W is a diagonal
// (weight) matrix defined by the user.

class NormalMatrix : public LinearOperator {
public:
    // Constructor stores a reference to the model. No data is copied. The model
    // must be valid as long as the object is used.
    explicit NormalMatrix(const Model& model);

    // Prepares normal matrix for subsequent calls to Apply(). If W is not NULL,
    // then W must hold n+m entries. No data is copied. The array must be valid
    // in each subsequent call to Apply(). If W is NULL, then the first n
    // entries are assumed 1.0 and the last m entries are assumed 0.0.
    void Prepare(const double* W);

    // Returns computation time for calls to Apply() since last reset_time().
    double time() const;
    void reset_time();

private:
    void _Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs) override;

    const Model& model_;
    const double* W_{nullptr};
    bool prepared_{false};
    Vector work_;            // size n+m workspace (2-pass matvec products only)
    double time_{0.0};
};

}  // namespace ipx

#endif  // IPX_NORMAL_MATRIX_H_
