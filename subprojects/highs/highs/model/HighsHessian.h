/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file model/HighsHessian.h
 * @brief
 */
#ifndef MODEL_HIGHS_HESSIAN_H_
#define MODEL_HIGHS_HESSIAN_H_

#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsCDouble.h"
#include "util/HighsUtils.h"

// class HighsHessian;

class HighsHessian {
 public:
  HighsHessian() { clear(); }
  std::string name_;
  HighsInt dim_;
  HessianFormat format_;
  std::vector<HighsInt> start_;
  std::vector<HighsInt> index_;
  std::vector<double> value_;
  HessianOracle oracle_;
  bool operator==(const HighsHessian& hessian) const;
  bool isOracle() const { return oracle_.call_ != nullptr; }
  void formFromOracle();
  void product(const std::vector<double>& solution,
               std::vector<double>& product) const;
  void alphaProductPlusY(const double alpha, const std::vector<double>& x,
                         std::vector<double>& y) const;
  double objectiveValue(const std::vector<double>& solution) const;
  HighsCDouble objectiveCDoubleValue(const std::vector<double>& solution) const;
  void exactResize();
  void deleteCols(const HighsIndexCollection& index_collection);
  void clear();
  bool formatOk() const {
    return (this->format_ == HessianFormat::kTriangular ||
            this->format_ == HessianFormat::kSquare);
  };
  bool scaleOk(const HighsInt cost_scale, const double small_matrix_value,
               const double large_matrix_value) const;
  HighsInt dim() const;
  HighsInt numNz() const;
  bool empty() const { return isEmpty(); }
  bool isEmpty() const { return this->dim() <= 0; }
  bool isDiagonal() const;
  double diag(HighsInt i) const;
  double diagonal(HighsInt i) const;
  HighsHessian toSquare() const;
  HighsStatus checkOracle(const HighsLogOptions& log_options,
                          const bool exit_on_first_error) const;
  void print(const std::string& message = "") const;
  bool legalDimensions(const HighsLogOptions& log_options) const;
  HighsStatus trimDimensions(const HighsLogOptions& log_options);
  HighsStatus assess(const HighsOptions& options,
                     const std::string* col_names = nullptr);
};

#endif
