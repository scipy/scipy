// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_INTERNAL_H_
#define IPX_INTERNAL_H_

#include "ipx_config.h"
#include "ipx_info.h"
#include "ipx_parameters.h"
#include "ipx_status.h"
#include <valarray>

namespace ipx {

using Int = ipxint;
using Info = ipx_info;
using Parameters = ipx_parameters;
using Vector = std::valarray<double>;

// A vector is treated sparse if it has no more than kHypersparseThreshold * dim
// nonzeros.
static constexpr double kHypersparseThreshold = 0.1;

// When LU factorization is used for rank detection, columns of the active
// submatrix whose maximum entry is <= kLuDependencyTol are removed immediately
// without choosing a pivot.
static constexpr double kLuDependencyTol = 1e-3;

// A fresh LU factorization is considered unstable if
//   ||b-Bx|| / (||b||+||B||*||x||) > kLuStabilityThreshold,
// where x=B\b is computed from the LU factors, b has components +/- 1 that are
// chosen to make x large, and ||.|| is the 1-norm. An unstable factorization
// triggers tightening of the pivot tolerance and refactorization.
static constexpr double kLuStabilityThreshold = 1e-12;

// A Forrest-Tomlin LU update is declared numerically unstable if the relative
// error in the new diagonal entry of U is larger than kFtDiagErrorTol.
static constexpr double kFtDiagErrorTol = 1e-8;

}  // namespace ipx

#endif // IPX_INTERNAL_H_
