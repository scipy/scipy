// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_UTILS_H_
#define IPX_UTILS_H_

#include <vector>
#include "ipx_internal.h"

namespace ipx {

bool AllFinite(const Vector& x);

double Onenorm(const Vector& x);
double Twonorm(const Vector& x);
double Infnorm(const Vector& x);
double Dot(const Vector& x, const Vector& y);

// Returns the index of an entry of maximum absolute value.
Int FindMaxAbs(const Vector& x);

// lhs[permuted_index] = rhs
void Permute(const std::vector<Int>& permuted_index, const Vector& rhs,
             Vector& lhs);

// lhs = rhs[permuted_index]
void PermuteBack(const std::vector<Int>& permuted_index, const Vector& rhs,
                 Vector& lhs);

// Returns the inverse permutation to @perm.
std::vector<Int> InversePerm(const std::vector<Int>& perm);

// Returns the permutation that puts values[0..m-1] in increasing (if reverse is
// false) or in decreasing (if reverse is true) order. If values==NULL, returns
// the identity permutation.
std::vector<Int> Sortperm(Int m, const double* values, bool reverse);

}  // namespace ipx

#endif  // IPX_UTILS_H_
