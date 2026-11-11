/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file model/HighsHessianUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef MODEL_HIGHSHESSIANUTILS_H_
#define MODEL_HIGHSHESSIANUTILS_H_

#include <vector>

#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"
#include "model/HighsHessian.h"

// class HighsHessian;
// class HighsOptions;

using std::vector;

HighsStatus assessHessian(HighsHessian& hessian, const HighsOptions& options,
                          const std::string* col_names = nullptr);
void completeHessianDiagonal(const HighsOptions& options,
                             HighsHessian& hessian);
bool okHessianDiagonal(const HighsOptions& options, HighsHessian& hessian,
                       const ObjSense sense = ObjSense::kMinimize);
HighsStatus normaliseHessian(const HighsOptions& options,
                             HighsHessian& hessian);
HighsStatus extractTriangularHessian(const HighsOptions& options,
                                     HighsHessian& hessian);
void triangularToSquareHessian(const HighsHessian& hessian,
                               vector<HighsInt>& start, vector<HighsInt>& index,
                               vector<double>& value);
void completeHessian(const HighsInt full_dim, HighsHessian& hessian);

void reportHessian(const HighsLogOptions& log_options, const HighsInt dim,
                   const HighsInt num_nz, const HighsInt* start,
                   const HighsInt* index, const double* value);

void userScaleHessian(HighsHessian& hessian, HighsUserScaleData& data,
                      const bool apply = true);
#endif  // MODEL_HIGHSHESSIANUTILS_H_
