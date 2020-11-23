// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_LINEAR_OPERATOR_H_
#define IPX_LINEAR_OPERATOR_H_

#include "ipx_internal.h"

namespace ipx {

class LinearOperator {
public:
    LinearOperator& operator=(const LinearOperator&) = delete;
    LinearOperator& operator=(LinearOperator&&) = delete;
    virtual ~LinearOperator() {}

    // Computes lhs = F(rhs), where F is a linear function. If rhs_dot_lhs is
    // not NULL, then the argument returns dot(rhs,lhs). The implementation can
    // assume that rhs and lhs do not refer to the same object.
    void Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs);

private:
    virtual void _Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs)
        = 0;
};

}  // namespace ipx

#endif  // IPX_LINEAR_OPERATOR_H_
