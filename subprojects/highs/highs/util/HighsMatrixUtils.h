/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef UTIL_HIGHSMATRIXUTILS_H_
#define UTIL_HIGHSMATRIXUTILS_H_

#include <cassert>
#include <vector>

// #include "lp_data/HighsStatus.h"
#include "lp_data/HighsOptions.h"

using std::vector;

// Entry point for HighsHessian
HighsStatus assessMatrix(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt vec_dim,
    const HighsInt num_vec, vector<HighsInt>& matrix_start,
    vector<HighsInt>& matrix_index, vector<double>& matrix_value,
    const double small_matrix_value, const double large_matrix_value,
    const bool sum_duplicates, const std::string* col_names = nullptr,
    const std::string* row_names = nullptr);

// Entry point for HighsSparseMatrix
HighsStatus assessMatrix(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt vec_dim,
    const HighsInt num_vec, const bool partitioned,
    vector<HighsInt>& matrix_start, vector<HighsInt>& matrix_p_end,
    vector<HighsInt>& matrix_index, vector<double>& matrix_value,
    const double small_matrix_value, const double large_matrix_value,
    const bool sum_duplicates, const std::string* col_names = nullptr,
    const std::string* row_names = nullptr);

bool legalMatrixDimensions(const HighsLogOptions& log_options,
                           const std::string& class_name,
                           const std::string& matrix_name,
                           const HighsInt num_vec, const bool partitioned,
                           const vector<HighsInt>& start,
                           const vector<HighsInt>& p_end,
                           const vector<HighsInt>& index,
                           const vector<double>& value);

HighsStatus trimMatrixDimensions(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt num_vec,
    const bool partitioned, vector<HighsInt>& start, vector<HighsInt>& p_end,
    vector<HighsInt>& index, vector<double>& value);

#endif  // UTIL_HIGHSMATRIXUTILS_H_
