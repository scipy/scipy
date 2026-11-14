/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/linalg.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_LINALG_HPP
#define PDLP_HIPDLP_LINALG_HPP

#include <vector>

#include "Highs.h"

namespace linalg {
double projectBox(double x, double l, double u);
double projectNonNegative(double y);
void projectBounds(const HighsLp& lp, std::vector<double>& x);

// Function to compute A*x for a given HighsLp and vector x
void ax(const HighsLp& lp, const std::vector<double>& x,
        std::vector<double>& result);

// Function to compute A^T*y for a given HighsLp and vector y
void aTy(const HighsLp& lp, const std::vector<double>& y,
         std::vector<double>& result);

double norm2(const std::vector<double>& vec);
void scale(std::vector<double>& vec, double factor);

void normalize(std::vector<double>& vec);

double dot(const std::vector<double>& a, const std::vector<double>& b);

double diffTwoNorm(const std::vector<double>& v1,
                   const std::vector<double>& v2);

// General norm functions
double vectorNorm(const std::vector<double>& vec, double p = 2.0);
double vectorNorm(const double* values, size_t size, double p = 2.0);
double vectorNormSquared(const std::vector<double>& vec);

// LP-specific norm calculations
double computeCostNorm(const HighsLp& lp, double p = 2.0);
double computeRhsNorm(const HighsLp& lp, double p = 2.0);

// Matrix column/row norm calculations
std::vector<double> computeColumnNorms(
    const HighsLp& lp, double p = std::numeric_limits<double>::infinity());
std::vector<double> computeRowNorms(
    const HighsLp& lp, double p = std::numeric_limits<double>::infinity());

std::vector<double> vectorSubtract(const std::vector<double>& a,
                                   const std::vector<double>& b);

}  // namespace linalg

#endif  // PDLP_HIPDLP_LINALG_HPP
