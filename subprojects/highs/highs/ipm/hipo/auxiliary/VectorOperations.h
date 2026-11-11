#ifndef HIPO_VECTOR_OPERATIONS_H
#define HIPO_VECTOR_OPERATIONS_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// =======================================================================
// COMPONENT-WISE VECTOR OPERATIONS
// =======================================================================

// v1[i] = alpha * v1[i] + beta * v2[i]
void vectorAdd(std::vector<double>& v1, double alpha,
               const std::vector<double>& v2, double beta);

// v1[i] + alpha
void vectorAdd(std::vector<double>& v1, const double alpha);

// v1[i] / v2[i]
void vectorDivide(std::vector<double>& v1, const std::vector<double>& v2);

// v1[i] * alpha
void vectorScale(std::vector<double>& v1, double alpha);

// =======================================================================

double dotProd(const std::vector<double>& v1, const std::vector<double>& v2);

double norm2(const std::vector<double>& x);

double infNorm(const std::vector<double>& x);

bool isNanVector(const std::vector<double>& x);

bool isInfVector(const std::vector<double>& x);

}  // namespace hipo

#endif