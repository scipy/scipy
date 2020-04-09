// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "linear_operator.h"

namespace ipx {

void LinearOperator::Apply(const Vector& rhs, Vector& lhs,
                           double* rhs_dot_lhs) {
    _Apply(rhs, lhs, rhs_dot_lhs);
}

}  // namespace ipx
