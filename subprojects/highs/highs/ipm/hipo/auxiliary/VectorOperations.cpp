#include "VectorOperations.h"

#include <cassert>
#include <cmath>

namespace hipo {

void vectorAdd(std::vector<double>& v1, double alpha,
               const std::vector<double>& v2, double beta) {
  for (Int i = 0; i < static_cast<Int>(v1.size()); ++i) {
    v1[i] = alpha * v1[i] + beta * v2[i];
  }
}

void vectorAdd(std::vector<double>& v1, const double alpha) {
  for (Int i = 0; i < static_cast<Int>(v1.size()); ++i) {
    v1[i] += alpha;
  }
}

void vectorDivide(std::vector<double>& v1, const std::vector<double>& v2) {
  for (Int i = 0; i < static_cast<Int>(v1.size()); ++i) {
    v1[i] /= v2[i];
  }
}

void vectorScale(std::vector<double>& v1, double alpha) {
  for (Int i = 0; i < static_cast<Int>(v1.size()); ++i) {
    v1[i] *= alpha;
  }
}

double dotProd(const std::vector<double>& v1, const std::vector<double>& v2) {
  double result{};
  for (Int i = 0; i < static_cast<Int>(v1.size()); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

double norm2(const std::vector<double>& x) {
  double norm{};
  for (Int i = 0; i < static_cast<Int>(x.size()); ++i) {
    norm += (x[i] * x[i]);
  }
  return std::sqrt(norm);
}

double infNorm(const std::vector<double>& x) {
  double norm{};
  for (Int i = 0; i < static_cast<Int>(x.size()); ++i) {
    norm = std::max(norm, std::fabs(x[i]));
  }
  return norm;
}

bool isNanVector(const std::vector<double>& x) {
  for (Int i = 0; i < static_cast<Int>(x.size()); ++i) {
    if (std::isnan(x[i])) return true;
  }
  return false;
}

bool isInfVector(const std::vector<double>& x) {
  for (Int i = 0; i < static_cast<Int>(x.size()); ++i) {
    if (std::isinf(x[i])) return true;
  }
  return false;
}

}  // namespace hipo
