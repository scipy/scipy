#include "ipm/ipx/utils.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include "../extern/pdqsort/pdqsort.h"

namespace ipx {

bool AllFinite(const Vector& x) {
    for (double xi : x)
        if (!std::isfinite(xi))
            return false;
    return true;
}

double Onenorm(const Vector& x) {
    double norm = 0.0;
    for (double xi : x)
        norm += std::abs(xi);
    return norm;
}

double Twonorm(const Vector& x) {
    double norm = 0.0;
    for (double xi : x)
        norm += xi * xi;
    return std::sqrt(norm);
}

double Infnorm(const Vector& x) {
    double norm = 0.0;
    for (double xi : x)
        norm = std::max(norm, std::abs(xi));
    return norm;
}

double Dot(const Vector& x, const Vector& y) {
    assert(x.size() == y.size());
    double d = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        d += x[i]*y[i];
    return d;
}

Int FindMaxAbs(const Vector& x) {
    double xmax = 0.0;
    size_t imax = 0;
    for (size_t i = 0; i < x.size(); i++) {
        if (std::abs(x[i]) > xmax) {
            xmax = std::abs(x[i]);
            imax = i;
        }
    }
    return static_cast<Int>(imax);
}

void Permute(const std::vector<Int>& permuted_index, const Vector& rhs,
             Vector& lhs) {
    for (size_t i = 0; i < permuted_index.size(); i++)
        lhs[permuted_index[i]] = rhs[i];
}

void PermuteBack(const std::vector<Int>& permuted_index, const Vector& rhs,
                 Vector& lhs) {
    for (size_t i = 0; i < permuted_index.size(); i++)
        lhs[i] = rhs[permuted_index[i]];
}

std::vector<Int> InversePerm(const std::vector<Int>& perm) {
    std::vector<Int> invperm(perm.size());
    // at() throws an exception if the index is out of range
    for (size_t i = 0; i < perm.size(); i++)
        invperm.at(perm[i]) = static_cast<Int>(i);
    return invperm;
}

std::vector<Int> Sortperm(Int m, const double* values, bool reverse) {
    std::vector<Int> perm(m);
    for (Int i = 0; i < m; i++) perm[i] = i;
    if (!values) return perm;

    if (reverse)
        pdqsort(perm.begin(), perm.end(), [&](Int i, Int j) {
            return std::make_pair(values[i], i) > std::make_pair(values[j], j);
        });
    else
        pdqsort(perm.begin(), perm.end(), [&](Int i, Int j) {
            return std::make_pair(values[i], i) < std::make_pair(values[j], j);
        });

    return perm;
}

}  // namespace ipx
