// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_STARTING_BASIS_H_
#define IPX_STARTING_BASIS_H_

#include "basis.h"
#include "iterate.h"

namespace ipx {

// Constructs a basis with the following properties:
//
// If lb[j]=-inf and ub[j]=inf, then the variable becomes either
// - basic with status BASIC_FREE, or
// - nonbasic with status NONBASIC_FIXED.                         (1)
//
// If lb[j]==ub[j] and j is a slack variable, then it becomes either
// - basic with status BASIC_FREE, or                             (2)
// - nonbasic with status NONBASIC_FIXED.
//
// If lb[j]==ub[j] and j is not a slack variable, then it becomes
// - nonbasic with status NONBASIC_FIXED.
//
// All other variables get status BASIC or NONBASIC.
//
// In case (1) the columns corresponding to free variables are linearly
// dependent. In case (2) the rows to equality constraints are linearly
// dependent. In each case the dependent variables are moved to zero without
// altering the primal or dual residual.
// TODO: we need to check for primal/dual infeasibility here.
//
// The method calls ConstructBasisFromWeights() using the interior point
// scaling factors as column weights. If a variable gets status BASIC_FREE or
// NONBASIC_FIXED, then its state in @iterate is changed accordingly to free or
// fixed. On return info->errflag is nonzero if an error occured.
//
void StartingBasis(Iterate* iterate, Basis* basis, Info* info);

}  // namespace ipx

#endif  // IPX_STARTING_BASIS_H_
