/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
 */
#include "Highs.h"

#include <algorithm>
#include <cassert>
#include <csignal>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>

#include "HighsExternalApi.h"
#include "io/Filereader.h"
#include "io/LoadOptions.h"
#include "ipm/IpxWrapper.h"
#include "lp_data/HighsCallbackStruct.h"
#include "lp_data/HighsInfoDebug.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsSolve.h"
#include "mip/HighsMipSolver.h"
#include "model/HighsHessianUtils.h"
#include "parallel/HighsParallel.h"
#include "presolve/ICrashX.h"
#include "simplex/HSimplex.h"
#include "simplex/HSimplexDebug.h"
#include "util/HighsMatrixPic.h"
#include "util/HighsSort.h"

#define STRINGFY(s) STRINGFY0(s)
#define STRINGFY0(s) #s
const char* highsVersion() {
  return STRINGFY(HIGHS_VERSION_MAJOR) "." STRINGFY(
      HIGHS_VERSION_MINOR) "." STRINGFY(HIGHS_VERSION_PATCH);
}
HighsInt highsVersionMajor() { return HIGHS_VERSION_MAJOR; }
HighsInt highsVersionMinor() { return HIGHS_VERSION_MINOR; }
HighsInt highsVersionPatch() { return HIGHS_VERSION_PATCH; }
const char* highsGithash() { return HIGHS_GITHASH; }
const char* highsCompilationDate() { return "deprecated"; }

Highs::Highs() : callback_(this) {}

std::string Highs::getThirdPartyNotice() const {
  return HighsExternalApi::getThirdPartyNotice<HighsExtras::all>();
}

HighsStatus Highs::clear() {
  resetOptions();
  return clearModel();
}

HighsStatus Highs::clearModel() {
  model_.clear();
  multi_linear_objective_.clear();
  saved_objective_and_solution_.clear();
  return clearSolver();
}

HighsStatus Highs::clearSolver() {
  clearDerivedModelProperties();
  invalidateSolverData();
  ekk_instance_.clear();
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::clearSolverDualData() {
  clearDerivedModelProperties();
  invalidateSolverDualData();
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::releaseMemory() {
  // 1. Clear all model and solver state (same as clear()).
  clearModel();

  // 2. Clear and shrink vectors that clearModel/clearSolver only
  //    invalidate or leave with residual capacity.  After this block the
  //    instance holds zero allocated heap memory for solver state —
  //    equivalent to a freshly constructed Highs.
  saved_objective_and_solution_.shrink_to_fit();
  solution_.clear();
  solution_.col_value.shrink_to_fit();
  solution_.col_dual.shrink_to_fit();
  solution_.row_value.shrink_to_fit();
  solution_.row_dual.shrink_to_fit();
  basis_.clear();
  basis_.col_status.shrink_to_fit();
  basis_.row_status.shrink_to_fit();
  ranging_.clear();
  ranging_.col_cost_up.value_.shrink_to_fit();
  ranging_.col_cost_up.objective_.shrink_to_fit();
  ranging_.col_cost_up.in_var_.shrink_to_fit();
  ranging_.col_cost_up.ou_var_.shrink_to_fit();
  ranging_.col_cost_dn.value_.shrink_to_fit();
  ranging_.col_cost_dn.objective_.shrink_to_fit();
  ranging_.col_cost_dn.in_var_.shrink_to_fit();
  ranging_.col_cost_dn.ou_var_.shrink_to_fit();
  ranging_.col_bound_up.value_.shrink_to_fit();
  ranging_.col_bound_up.objective_.shrink_to_fit();
  ranging_.col_bound_up.in_var_.shrink_to_fit();
  ranging_.col_bound_up.ou_var_.shrink_to_fit();
  ranging_.col_bound_dn.value_.shrink_to_fit();
  ranging_.col_bound_dn.objective_.shrink_to_fit();
  ranging_.col_bound_dn.in_var_.shrink_to_fit();
  ranging_.col_bound_dn.ou_var_.shrink_to_fit();
  ranging_.row_bound_up.value_.shrink_to_fit();
  ranging_.row_bound_up.objective_.shrink_to_fit();
  ranging_.row_bound_up.in_var_.shrink_to_fit();
  ranging_.row_bound_up.ou_var_.shrink_to_fit();
  ranging_.row_bound_dn.value_.shrink_to_fit();
  ranging_.row_bound_dn.objective_.shrink_to_fit();
  ranging_.row_bound_dn.in_var_.shrink_to_fit();
  ranging_.row_bound_dn.ou_var_.shrink_to_fit();
  iis_.clear();

  // 3. Shrink standard form LP vectors.
  standard_form_cost_.shrink_to_fit();
  standard_form_rhs_.shrink_to_fit();

  // 4. Reset timer clocks (accumulated timing data).
  timer_.zeroAllClocks();

  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::setOptionValue(const std::string& option, const bool value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const HighsInt value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const double value) {
  if (setLocalOptionValue(options_.log_options, option, options_.records,
                          value) == OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const std::string& value) {
  HighsLogOptions report_log_options = options_.log_options;
  if (setLocalOptionValue(report_log_options, option, options_.log_options,
                          options_.records, value) == OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::setOptionValue(const std::string& option,
                                  const char* value) {
  HighsLogOptions report_log_options = options_.log_options;
  if (setLocalOptionValue(report_log_options, option, options_.log_options,
                          options_.records, value) == OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::readOptions(const std::string& filename) {
  if (filename.size() <= 0) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Empty file name so not reading options\n");
    return HighsStatus::kWarning;
  }
  HighsLogOptions report_log_options = options_.log_options;
  switch (loadOptionsFromFile(report_log_options, options_, filename)) {
    case HighsLoadOptionsStatus::kError:
    case HighsLoadOptionsStatus::kEmpty:
      return HighsStatus::kError;
    default:
      break;
  }
  return optionChangeAction();
}

HighsStatus Highs::passOptions(const HighsOptions& options) {
  if (passLocalOptions(options_.log_options, options, options_) ==
      OptionStatus::kOk)
    return optionChangeAction();
  return HighsStatus::kError;
}

HighsStatus Highs::resetOptions() {
  this->logHeader();
  resetLocalOptions(options_.records);
  return optionChangeAction();
}

HighsStatus Highs::writeOptions(const std::string& filename,
                                const bool report_only_deviations) const {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  FILE* file;
  HighsFileType file_type;
  return_status = interpretCallStatus(
      options_.log_options,
      openWriteFile(filename, "writeOptions", file, file_type), return_status,
      "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  if (filename == "") file_type = HighsFileType::kMinimal;
  // Report to user that options are being written to a file
  if (filename != "")
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Writing the option values to %s\n", filename.c_str());
  return_status = interpretCallStatus(
      options_.log_options,
      writeOptionsToFile(file, options_.log_options, options_.records,
                         report_only_deviations, file_type),
      return_status, "writeOptionsToFile");
  if (file != stdout) fclose(file);
  return return_status;
}

// HighsStatus Highs::getOptionType(const char* option, HighsOptionType* type)
// const { return getOptionType(option, type);}

HighsStatus Highs::getOptionName(const HighsInt index,
                                 std::string* name) const {
  HighsInt num_option = this->options_.records.size();
  if (index < 0 || index >= num_option) {
    highsLogUser(
        options_.log_options, HighsLogType::kInfo,
        "Index %d passed to getOptionName is not in the valid range [0, %d)\n",
        int(index), int(num_option));
    return HighsStatus::kError;
  }
  *name = this->options_.records[index]->name;
  return HighsStatus::kOk;
}

HighsStatus Highs::getOptionType(const std::string& option,
                                 HighsOptionType* type) const {
  if (getLocalOptionType(options_.log_options, option, options_.records,
                         type) == OptionStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getBoolOptionValues(const std::string& option,
                                       bool* current_value,
                                       bool* default_value) const {
  if (getLocalOptionValues(options_.log_options, option, options_.records,
                           current_value, default_value) != OptionStatus::kOk)
    return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::getIntOptionValues(const std::string& option,
                                      HighsInt* current_value,
                                      HighsInt* min_value, HighsInt* max_value,
                                      HighsInt* default_value) const {
  if (getLocalOptionValues(options_.log_options, option, options_.records,
                           current_value, min_value, max_value,
                           default_value) != OptionStatus::kOk)
    return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::getDoubleOptionValues(const std::string& option,
                                         double* current_value,
                                         double* min_value, double* max_value,
                                         double* default_value) const {
  if (getLocalOptionValues(options_.log_options, option, options_.records,
                           current_value, min_value, max_value,
                           default_value) != OptionStatus::kOk)
    return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::getStringOptionValues(const std::string& option,
                                         std::string* current_value,
                                         std::string* default_value) const {
  if (getLocalOptionValues(options_.log_options, option, options_.records,
                           current_value, default_value) != OptionStatus::kOk)
    return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::getRunDataValue(const std::string& run_data,
                                   HighsInt& value) const {
  this->logHeader();
  RunDataStatus status =
      getLocalRunDataValue(options_.log_options, run_data, run_data_.valid,
                           run_data_.records, value);
  if (status == RunDataStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == RunDataStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

#ifndef HIGHSINT64
HighsStatus Highs::getRunDataValue(const std::string& run_data,
                                   int64_t& value) const {
  this->logHeader();
  RunDataStatus status =
      getLocalRunDataValue(options_.log_options, run_data, run_data_.valid,
                           run_data_.records, value);
  if (status == RunDataStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == RunDataStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}
#endif

HighsStatus Highs::getRunDataType(const std::string& run_data,
                                  HighsRunDataType& type) const {
  this->logHeader();
  if (getLocalRunDataType(options_.log_options, run_data, run_data_.records,
                          type) == RunDataStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getRunDataValue(const std::string& run_data,
                                   double& value) const {
  this->logHeader();
  RunDataStatus status =
      getLocalRunDataValue(options_.log_options, run_data, run_data_.valid,
                           run_data_.records, value);
  if (status == RunDataStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == RunDataStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

HighsStatus Highs::writeRunData(const std::string& filename) const {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  FILE* file;
  HighsFileType file_type;
  return_status = interpretCallStatus(
      options_.log_options,
      openWriteFile(filename, "writeRunData", file, file_type), return_status,
      "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  // Report to user that options are being written to a file
  if (filename != "")
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Writing the run_data values to %s\n", filename.c_str());
  return_status = interpretCallStatus(
      options_.log_options,
      writeRunDataToFile(file, run_data_.valid, run_data_.records, file_type),
      return_status, "writeRunDataToFile");
  if (file != stdout) fclose(file);
  return return_status;
}

HighsStatus Highs::getInfoValue(const std::string& info,
                                HighsInt& value) const {
  this->logHeader();
  InfoStatus status = getLocalInfoValue(options_.log_options, info, info_.valid,
                                        info_.records, value);
  if (status == InfoStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == InfoStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

#ifndef HIGHSINT64
HighsStatus Highs::getInfoValue(const std::string& info, int64_t& value) const {
  this->logHeader();
  InfoStatus status = getLocalInfoValue(options_.log_options, info, info_.valid,
                                        info_.records, value);
  if (status == InfoStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == InfoStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}
#endif

HighsStatus Highs::getInfoType(const std::string& info,
                               HighsInfoType& type) const {
  this->logHeader();
  if (getLocalInfoType(options_.log_options, info, info_.records, type) ==
      InfoStatus::kOk)
    return HighsStatus::kOk;
  return HighsStatus::kError;
}

HighsStatus Highs::getInfoValue(const std::string& info, double& value) const {
  this->logHeader();
  InfoStatus status = getLocalInfoValue(options_.log_options, info, info_.valid,
                                        info_.records, value);
  if (status == InfoStatus::kOk) {
    return HighsStatus::kOk;
  } else if (status == InfoStatus::kUnavailable) {
    return HighsStatus::kWarning;
  } else {
    return HighsStatus::kError;
  }
}

HighsStatus Highs::writeInfo(const std::string& filename) const {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  FILE* file;
  HighsFileType file_type;
  return_status =
      interpretCallStatus(options_.log_options,
                          openWriteFile(filename, "writeInfo", file, file_type),
                          return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  // Report to user that options are being written to a file
  if (filename != "")
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Writing the info values to %s\n", filename.c_str());
  return_status = interpretCallStatus(
      options_.log_options,
      writeInfoToFile(file, info_.valid, info_.records, file_type),
      return_status, "writeInfoToFile");
  if (file != stdout) fclose(file);
  return return_status;
}

/**
 * @brief Get the size of HighsInt
 */
// HighsInt getSizeofHighsInt() {

// Methods below change the incumbent model or solver information
// associated with it. Hence returnFromHighs is called at the end of
// each
HighsStatus Highs::passModel(HighsModel model) {
  this->logHeader();
  // This is the "master" Highs::passModel, in that all the others
  // (and readModel) eventually call it
  //
  HighsLogOptions& log_options = options_.log_options;
  // Possibly analyse the LP data
  if (kHighsAnalysisLevelModelData & options_.highs_analysis_level)
    analyseLp(log_options, model.lp_);
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsHessian& hessian = model_.hessian_;
  // Clear the incumbent model and any associated data
  clearModel();
  // Move the model's LP and Hessian to the internal LP and Hessian
  lp = std::move(model.lp_);
  hessian = std::move(model.hessian_);

  auto clearModelAndReturnError = [&]() {
    // Clear the incumbent model and any associated data
    this->clearModel();
    return HighsStatus::kError;
  };

  if (lp.num_col_ == 0 || lp.num_row_ == 0) {
    // Model constraint matrix has either no columns or no
    // rows. Clearly the matrix is empty, so may have no orientation
    // or starts assigned. HiGHS assumes that such a model will have
    // null starts, so make it column-wise
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Model has either no columns or no rows, so ignoring user "
                 "constraint matrix data and initialising empty matrix\n");
    lp.a_matrix_.setEmptyColwise(lp.num_col_);
  } else {
    // Matrix has rows and columns, so a_matrix format must be valid
    if (!lp.a_matrix_.formatOk()) return clearModelAndReturnError();
  }
  lp.origin_name_ = "User";
  lp.a_matrix_.name_ = "User matrix";
  assert(lp.a_matrix_.formatOk());
  // Dimensions in a_matrix_ may not be set, so take them from lp.
  lp.setMatrixDimensions();
  // Residual scale factors may be present. ToDo Allow user-defined
  // scale factors to be passed
  assert(!lp.is_scaled_);
  assert(!lp.is_moved_);
  lp.resetScale();
  // Check that the LP array dimensions are legal (ie at least the
  // size required by num_col_ and num_row_)
  if (!lp.legalDimensions("passModel", log_options))
    return clearModelAndReturnError();
  // Check that the Hessian format is valid
  if (!hessian.formatOk()) return clearModelAndReturnError();
  // Check validity of the LP, normalising its values
  return_status = interpretCallStatus(log_options, assessLp(lp, options_),
                                      return_status, "assessLp");
  if (return_status == HighsStatus::kError) return clearModelAndReturnError();
  // Now legality of matrix is established, ensure that it is
  // column-wise
  lp.ensureColwise();
  // First check raw consistency of LP and Hessian (column dimensions)
  if (this->model_.lpHessianDimensionsCompatible(log_options) ==
      HighsStatus::kError)
    return clearModelAndReturnError();

  // Check validity of any Hessian, normalising its entries
  return_status =
      interpretCallStatus(log_options, assessHessian(hessian, options_),
                          return_status, "assessHessian");
  if (return_status == HighsStatus::kError) return clearModelAndReturnError();
  if (hessian.dim_) {
    // Clear any zero Hessian
    if (hessian.numNz() == 0) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Hessian has dimension %" HIGHSINT_FORMAT
                   " but no nonzeros, so is ignored\n",
                   hessian.dim_);
      hessian.clear();
    }
  }
  // Ensure that any non-zero Hessian of dimension less than the
  // number of columns in the model is completed
  if (hessian.dim_) completeHessian(this->model_.lp_.num_col_, hessian);
  // Possibly create image files of the constraint matrix and Hessian
  if (options_.write_matrix_image)
    writeLpMatrixPicToFile(options_, "LpMatrix", model_.lp_);
  if (options_.write_hessian_image)
    writeHessianPicToFile(options_, "Hessian", model_.hessian_);

  // Clear solver status, solution, basis and info associated with any
  // previous model
  return_status = interpretCallStatus(log_options, clearSolver(), return_status,
                                      "clearSolver");
  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(HighsLp lp) {
  this->logHeader();
  HighsModel model;
  model.lp_ = std::move(lp);
  return passModel(std::move(model));
}

HighsStatus Highs::passModel(
    const HighsInt num_col, const HighsInt num_row, const HighsInt a_num_nz,
    const HighsInt q_num_nz, const HighsInt a_format, const HighsInt q_format,
    const HighsInt sense, const double offset, const double* costs,
    const double* col_lower, const double* col_upper, const double* row_lower,
    const double* row_upper, const HighsInt* a_start, const HighsInt* a_index,
    const double* a_value, const HighsInt* q_start, const HighsInt* q_index,
    const double* q_value, const HighsInt* integrality) {
  this->logHeader();
  HighsModel model;
  HighsLp& lp = model.lp_;
  // Check that the formats of the constraint matrix and Hessian are valid
  if (!aFormatOk(a_num_nz, a_format)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model has illegal constraint matrix format\n");
    return HighsStatus::kError;
  }
  if (!qFormatOk(q_num_nz, q_format)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model has illegal Hessian matrix format\n");
    return HighsStatus::kError;
  }
  const bool a_rowwise =
      a_num_nz > 0 ? a_format == (HighsInt)MatrixFormat::kRowwise : false;
  //  if (num_nz) a_rowwise = a_format == (HighsInt)MatrixFormat::kRowwise;

  lp.num_col_ = num_col;
  lp.num_row_ = num_row;
  if (num_col > 0) {
    assert(costs != NULL);
    assert(col_lower != NULL);
    assert(col_upper != NULL);
    lp.col_cost_.assign(costs, costs + num_col);
    lp.col_lower_.assign(col_lower, col_lower + num_col);
    lp.col_upper_.assign(col_upper, col_upper + num_col);
  }
  if (num_row > 0) {
    assert(row_lower != NULL);
    assert(row_upper != NULL);
    lp.row_lower_.assign(row_lower, row_lower + num_row);
    lp.row_upper_.assign(row_upper, row_upper + num_row);
  }
  if (a_num_nz > 0) {
    assert(num_col > 0);
    assert(num_row > 0);
    assert(a_start != NULL);
    assert(a_index != NULL);
    assert(a_value != NULL);
    if (a_rowwise) {
      lp.a_matrix_.start_.assign(a_start, a_start + num_row);
    } else {
      lp.a_matrix_.start_.assign(a_start, a_start + num_col);
    }
    lp.a_matrix_.index_.assign(a_index, a_index + a_num_nz);
    lp.a_matrix_.value_.assign(a_value, a_value + a_num_nz);
  }
  if (a_rowwise) {
    lp.a_matrix_.start_.resize(num_row + 1);
    lp.a_matrix_.start_[num_row] = a_num_nz;
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  } else {
    lp.a_matrix_.start_.resize(num_col + 1);
    lp.a_matrix_.start_[num_col] = a_num_nz;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
  }
  if (sense == (HighsInt)ObjSense::kMaximize) {
    lp.sense_ = ObjSense::kMaximize;
  } else {
    lp.sense_ = ObjSense::kMinimize;
  }
  lp.offset_ = offset;
  if (num_col > 0 && integrality != NULL) {
    lp.integrality_.resize(num_col);
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      HighsInt integrality_status = integrality[iCol];
      const bool legal_integrality_status =
          integrality_status == HighsInt(HighsVarType::kContinuous) ||
          integrality_status == HighsInt(HighsVarType::kInteger) ||
          integrality_status == HighsInt(HighsVarType::kSemiContinuous) ||
          integrality_status == HighsInt(HighsVarType::kSemiInteger) ||
          integrality_status == HighsInt(HighsVarType::kImplicitInteger);
      if (!legal_integrality_status) {
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Model has illegal integer value of %d (type %s) for "
                     "integrality[%d]\n",
                     int(integrality_status),
                     highsVarTypeToString(integrality_status).c_str(), iCol);
        return HighsStatus::kError;
      }
      lp.integrality_[iCol] = (HighsVarType)integrality_status;
    }
  }
  if (q_num_nz > 0) {
    assert(num_col > 0);
    assert(q_start != NULL);
    assert(q_index != NULL);
    assert(q_value != NULL);
    HighsHessian& hessian = model.hessian_;
    hessian.dim_ = num_col;
    hessian.format_ = static_cast<HessianFormat>(q_format);
    hessian.start_.assign(q_start, q_start + num_col);
    hessian.start_.resize(num_col + 1);
    hessian.start_[num_col] = q_num_nz;
    hessian.index_.assign(q_index, q_index + q_num_nz);
    hessian.value_.assign(q_value, q_value + q_num_nz);
  }
  return passModel(std::move(model));
}

HighsStatus Highs::passModel(const HighsInt num_col, const HighsInt num_row,
                             const HighsInt num_nz, const HighsInt a_format,
                             const HighsInt sense, const double offset,
                             const double* costs, const double* col_lower,
                             const double* col_upper, const double* row_lower,
                             const double* row_upper, const HighsInt* a_start,
                             const HighsInt* a_index, const double* a_value,
                             const HighsInt* integrality) {
  this->logHeader();
  return passModel(num_col, num_row, num_nz, 0, a_format, 0, sense, offset,
                   costs, col_lower, col_upper, row_lower, row_upper, a_start,
                   a_index, a_value, NULL, NULL, NULL, integrality);
}

HighsStatus Highs::passHessian(HighsHessian hessian_) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  HighsHessian& hessian = this->model_.hessian_;
  hessian = std::move(hessian_);

  auto clearHessianAndReturnError = [&]() {
    // Clear the incumbent hessian and any associated data
    hessian.clear();
    return HighsStatus::kError;
  };

  // First check raw consistency of LP and Hessian (column dimensions)
  if (this->model_.lpHessianDimensionsCompatible(options_.log_options) ==
      HighsStatus::kError)
    return clearHessianAndReturnError();

  // Check validity of any Hessian, normalising its entries
  return_status = interpretCallStatus(options_.log_options,
                                      assessHessian(hessian, options_),
                                      return_status, "assessHessian");
  if (return_status == HighsStatus::kError) return clearHessianAndReturnError();
  if (hessian.dim_) {
    // Clear any zero Hessian
    if (hessian.numNz() == 0) {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Hessian has dimension %" HIGHSINT_FORMAT
                   " but no nonzeros, so is ignored\n",
                   hessian.dim_);
      hessian.clear();
    }
  }
  // Ensure that any non-zero Hessian of dimension less than the
  // number of columns in the model is completed
  if (hessian.dim_) completeHessian(this->model_.lp_.num_col_, hessian);

  return_status = interpretCallStatus(options_.log_options, clearSolver(),
                                      return_status, "clearSolver");
  return returnFromHighs(return_status);
}

HighsStatus Highs::passHessian(const HighsInt dim, const HighsInt num_nz,
                               const HighsInt format, const HighsInt* start,
                               const HighsInt* index, const double* value) {
  this->logHeader();
  HighsHessian hessian;
  if (!qFormatOk(num_nz, format)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model has illegal Hessian matrix format\n");
    return HighsStatus::kError;
  }
  HighsInt num_col = model_.lp_.num_col_;
  if (dim != num_col) return HighsStatus::kError;
  hessian.dim_ = num_col;
  hessian.format_ = static_cast<HessianFormat>(format);
  if (dim > 0) {
    assert(start != NULL);
    hessian.start_.assign(start, start + num_col);
    hessian.start_.resize(num_col + 1);
    hessian.start_[num_col] = num_nz;
  }
  if (num_nz > 0) {
    assert(index != NULL);
    assert(value != NULL);
    hessian.index_.assign(index, index + num_nz);
    hessian.value_.assign(value, value + num_nz);
  }
  return passHessian(hessian);
}

HighsStatus Highs::passHessian(const HighsInt dim,
                               HighsHessianFunctionType oracleCall,
                               void* oracle_data,
                               HighsCHessianFunctionType c_oracleCall) {
  this->logHeader();
  if (dim <= 0) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Ignoring Hessian oracle data since dimension is %d\n",
                 int(dim));
    return HighsStatus::kWarning;
  }
  HighsHessian test_hessian;
  HessianOracle& oracle = test_hessian.oracle_;
  oracle.dim_ = dim;
  if (c_oracleCall) {
    oracle.call_ = [c_oracleCall](
                       const HighsInt type, const HighsInt* x_num_entries,
                       const HighsInt* x_index, const double* x_value,
                       HighsInt* hessian_x_num_entries,
                       HighsInt* hessian_x_index, double* hessian_x_value,
                       void* data) {
      return c_oracleCall(type, x_num_entries, x_index, x_value,
                          hessian_x_num_entries, hessian_x_index,
                          hessian_x_value, data);
    };
  } else {
    oracle.call_ = std::move(oracleCall);
  }
  oracle.data_ = oracle_data;
  // Check whether the new oracle is valid
  if (!oracle.isValid()) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot solve QP when Hessian oracle has no product call\n");
    return HighsStatus::kError;
  }
  // Update the incumbent Hessian
  this->model_.hessian_ = test_hessian;
  return HighsStatus::kOk;
}

HighsStatus Highs::checkHessianOracle(const bool exit_on_first_error) const {
  this->logHeader();
  return this->model_.hessian_.checkOracle(options_.log_options,
                                           exit_on_first_error);
}

HighsStatus Highs::passLinearObjectives(
    const HighsInt num_linear_objective,
    const HighsLinearObjective* linear_objective) {
  this->logHeader();
  if (num_linear_objective < 0) return HighsStatus::kOk;
  this->multi_linear_objective_.clear();
  for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++)
    if (this->addLinearObjective(linear_objective[iObj], iObj) !=
        HighsStatus::kOk)
      return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus Highs::addLinearObjective(
    const HighsLinearObjective& linear_objective, const HighsInt iObj) {
  this->logHeader();
  if (model_.isQp()) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot define additional linear objective for QP\n");
    return HighsStatus::kError;
  }
  if (!this->validLinearObjective(linear_objective, iObj))
    return HighsStatus::kError;
  this->multi_linear_objective_.push_back(linear_objective);
  return HighsStatus::kOk;
}

HighsStatus Highs::clearLinearObjectives() {
  this->logHeader();
  this->multi_linear_objective_.clear();
  return HighsStatus::kOk;
}

HighsStatus Highs::passColName(const HighsInt col, const std::string& name) {
  this->logHeader();
  HighsLp& lp = this->model_.lp_;
  const HighsInt num_col = lp.num_col_;
  if (col < 0 || col >= num_col) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Index %d for column name \"%s\" is outside the range [0, "
                 "num_col = %d)\n",
                 int(col), name.c_str(), int(num_col));
    return HighsStatus::kError;
  }
  if (int(name.length()) <= 0) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot define empty column names\n");
    return HighsStatus::kError;
  }
  lp.col_names_.resize(num_col);
  if (!lp.col_hash_.size()) lp.col_hash_.form(lp.col_names_);
  const bool is_duplicate =
      lp.col_hash_.updateFindsDuplicate(col, lp.col_names_[col], name);
  if (is_duplicate)
    highsLogUser(
        options_.log_options, HighsLogType::kWarning,
        "Name \"%s\" passed for column %d is a duplicate, but accepted\n",
        name.c_str(), int(col));
  lp.col_names_[col] = name;
  assert(lp.col_hash_.ok(lp.col_names_));
  return is_duplicate ? HighsStatus::kWarning : HighsStatus::kOk;
}

HighsStatus Highs::passRowName(const HighsInt row, const std::string& name) {
  this->logHeader();
  HighsLp& lp = this->model_.lp_;
  const HighsInt num_row = lp.num_row_;
  if (row < 0 || row >= num_row) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Index %d for row name \"%s\" is outside the range [0, num_row = %d)\n",
        int(row), name.c_str(), int(num_row));
    return HighsStatus::kError;
  }
  if (int(name.length()) <= 0) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot define empty row names\n");
    return HighsStatus::kError;
  }
  lp.row_names_.resize(num_row);
  if (!lp.row_hash_.size()) lp.row_hash_.form(lp.row_names_);
  const bool is_duplicate =
      lp.row_hash_.updateFindsDuplicate(row, lp.row_names_[row], name);
  if (is_duplicate)
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Name \"%s\" passed for row %d is a duplicate, but accepted\n",
                 name.c_str(), int(row));
  lp.row_names_[row] = name;
  assert(lp.row_hash_.ok(lp.row_names_));
  return is_duplicate ? HighsStatus::kWarning : HighsStatus::kOk;
}

HighsStatus Highs::passModelName(const std::string& name) {
  this->logHeader();
  this->model_.lp_.model_name_ = name;
  return HighsStatus::kOk;
}

HighsStatus Highs::readModel(const std::string& filename) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  Filereader* reader =
      Filereader::getFilereader(options_.log_options, filename);
  if (reader == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model file %s not supported\n", filename.c_str());
    return HighsStatus::kError;
  }

  HighsModel model;
  FilereaderRetcode call_code =
      reader->readModelFromFile(options_, filename, model);
  delete reader;
  if (call_code != FilereaderRetcode::kOk) {
    interpretFilereaderRetcode(options_.log_options, filename, call_code);
    const HighsStatus call_status = call_code == FilereaderRetcode::kWarning
                                        ? HighsStatus::kWarning
                                        : HighsStatus::kError;
    return_status = interpretCallStatus(options_.log_options, call_status,
                                        return_status, "readModelFromFile");
    if (return_status == HighsStatus::kError) return return_status;
  }
  model.lp_.model_name_ = extractModelName(filename);
  const bool remove_rows_of_count_1 = false;
  if (remove_rows_of_count_1) {
    // .lp files from PWSC (notably st-test23.lp) have bounds for
    // semi-continuous variables in the constraints section. By default,
    // these are interpreted as constraints, so the semi-continuous
    // variables are not set up correctly. Fix is to remove all rows of
    // count 1, interpreting their bounds as bounds on the corresponding
    // variable.
    removeRowsOfCountOne(options_.log_options, model.lp_);
  }
  return_status =
      interpretCallStatus(options_.log_options, passModel(std::move(model)),
                          return_status, "passModel");
  return returnFromHighs(return_status);
}

HighsStatus Highs::matrixImage(
    const std::string& matrix_image_filename,
    const std::string& hessian_image_filename) const {
  this->logHeader();
  HighsStatus status = HighsStatus::kOk;
  if (matrix_image_filename != "") {
    status =
        writeLpMatrixPicToFile(options_, matrix_image_filename, model_.lp_);
    if (status != HighsStatus::kOk) return status;
  }
  if (hessian_image_filename != "")
    status = writeHessianPicToFile(options_, hessian_image_filename,
                                   model_.hessian_);
  return status;
}

HighsStatus Highs::readBasis(const std::string& filename) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  // Try to read basis file into read_basis
  HighsBasis read_basis = basis_;
  return_status = interpretCallStatus(
      options_.log_options,
      readBasisFile(options_.log_options, model_.lp_, read_basis, filename),
      return_status, "readBasis");
  if (return_status != HighsStatus::kOk) return return_status;
  // Basis read OK: check whether it's consistent with the LP
  if (!isBasisConsistent(model_.lp_, read_basis)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "readBasis: invalid basis\n");
    return HighsStatus::kError;
  }
  // Update the HiGHS basis and invalidate any simplex basis for the model
  basis_ = read_basis;
  basis_.valid = true;
  basis_.useful = true;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
}

HighsStatus Highs::writeModel(const std::string& filename) {
  this->logHeader();
  return writeLocalModel(model_, filename);
}

HighsStatus Highs::writePresolvedModel(const std::string& filename) {
  this->logHeader();
  return writeLocalModel(presolved_model_, filename);
}

HighsStatus Highs::writeIisModel(const std::string& filename) {
  this->logHeader();
  HighsStatus return_status = this->getIisInterface();
  if (return_status == HighsStatus::kError) return return_status;
  return writeLocalModel(iis_.model_, filename);
}

HighsStatus Highs::writeLocalModel(HighsModel& model,
                                   const std::string& filename) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsLp& lp = model.lp_;

  // Dimensions in a_matrix_ may not be set, so take them from lp
  lp.setMatrixDimensions();

  // Normalise names according to file type
  HighsFileType file_type = getFileType(filename);
  call_status = normaliseNames(this->options_.log_options, lp, file_type);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "normaliseNames");
  assert(call_status != HighsStatus::kError);

  // Ensure that the LP is column-wise
  lp.ensureColwise();

  // Ensure that the dimensions are OK
  if (!lp.legalDimensions("writeLocalModel", options_.log_options))
    return HighsStatus::kError;

  if (model.hessian_.dim_ > 0 &&
      !model.hessian_.legalDimensions(options_.log_options))
    return HighsStatus::kError;

  // Check that the matrix starts are OK
  call_status = lp.a_matrix_.assessStart(options_.log_options);
  if (call_status == HighsStatus::kError) return call_status;

  // Check that the matrix indices are within bounds
  call_status = lp.a_matrix_.assessIndexBounds(options_.log_options);
  if (call_status == HighsStatus::kError) return call_status;

  // Check for repeated column or row names that would corrupt the file
  if (model.lp_.col_hash_.hasDuplicate(model.lp_.col_names_)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model has repeated column names\n");
    return returnFromHighs(HighsStatus::kError);
  }
  if (model.lp_.row_hash_.hasDuplicate(model.lp_.row_names_)) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Model has repeated row names\n");
    return returnFromHighs(HighsStatus::kError);
  }
  if (filename == "") {
    // Empty file name: report model on logging stream
    reportModel(model);
  } else {
    if (model.hessian_.isOracle()) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Cannot write QP to an %s file with Hessian represented via "
                   "an oracle\n",
                   getFilenameExt(filename).c_str());
      return returnFromHighs(HighsStatus::kError);
    }
    assert(lp.a_matrix_.isColwise());
    if (model.hessian_.dim_)
      assert(model.hessian_.format_ == HessianFormat::kTriangular);
    Filereader* writer =
        Filereader::getFilereader(options_.log_options, filename);
    if (writer == NULL) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Model file %s not supported\n", filename.c_str());
      return HighsStatus::kError;
    }
    // Report to user that model is being written
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Writing the model to %s\n", filename.c_str());
    return_status =
        interpretCallStatus(options_.log_options,
                            writer->writeModelToFile(options_, filename, model),
                            return_status, "writeModelToFile");
    delete writer;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::writeBasis(const std::string& filename) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  FILE* file;
  HighsFileType file_type;
  call_status = openWriteFile(filename, "writeBasis", file, file_type);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  call_status = normaliseNames(this->options_.log_options, this->model_.lp_);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "normaliseNames");
  assert(call_status != HighsStatus::kError);

  // Report to user that basis is being written
  if (filename != "") {
    if (!basis_.valid) {
      highsLogUser(options_.log_options, HighsLogType::kWarning,
                   "No basis to write: generated null basis file %s\n",
                   filename.c_str());
      return_status = HighsStatus::kWarning;
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Writing the basis to %s\n", filename.c_str());
    }
  }
  writeBasisFile(file, options_, model_.lp_, basis_);
  if (file != stdout) fclose(file);
  return return_status;
}

HighsStatus Highs::presolve() {
  this->logHeader();
  const HighsLogOptions& log_options = options_.log_options;
  if (model_.needsMods(options_.infinite_cost)) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Model contains infinite costs or semi-variables, so cannot "
                 "be presolved independently\n");
    return HighsStatus::kError;
  }
  if (model_.isQp()) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Model is a QP, for which no presolve techniques are implemented\n");
    return HighsStatus::kError;
  }
  HighsStatus return_status = HighsStatus::kOk;
  this->reportModelStats();
  clearPresolve();
  if (model_.isEmpty()) {
    model_presolve_status_ = HighsPresolveStatus::kNotReduced;
  } else {
    const bool force_presolve = true;
    // make sure scheduler is initialized before calling presolve, since
    // MIP presolve may use parallelism
    //
    return_status = this->initializeMultiThreading();
    if (return_status != HighsStatus::kOk) return return_status;
    HighsProfiling profiling;
    const bool already_profiling = this->profiling_;
    if (!already_profiling) this->initializeProfiling(&profiling);
    // If problem is a MIP and solve_relaxation is true, it's natural
    // to force LP presolve. It means that someone wanting the
    // presolved relaxation of a MIP just has to set solve_relaxation
    // to true, rather than clearing the integrality of the model (cf
    // #2234)
    const bool force_lp_presolve = options_.solve_relaxation;
    model_presolve_status_ = runPresolve(force_lp_presolve, force_presolve);
    if (!already_profiling) this->clearProfiling();
  }

  bool using_reduced_lp = false;
  reportPresolveReductions(log_options, model_presolve_status_, model_.lp_,
                           presolve_.getReducedProblem());
  switch (model_presolve_status_) {
    case HighsPresolveStatus::kNotPresolved: {
      // Shouldn't happen
      assert(model_presolve_status_ != HighsPresolveStatus::kNotPresolved);
      return_status = HighsStatus::kError;
      break;
    }
    case HighsPresolveStatus::kNotReduced:
    case HighsPresolveStatus::kInfeasible:
    case HighsPresolveStatus::kReduced:
    case HighsPresolveStatus::kReducedToEmpty:
    case HighsPresolveStatus::kUnboundedOrInfeasible: {
      // All OK
      if (model_presolve_status_ == HighsPresolveStatus::kInfeasible) {
        // Infeasible model, so indicate that the incumbent model is
        // known as such
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kInfeasible);
      } else if (model_presolve_status_ == HighsPresolveStatus::kNotReduced) {
        // No reduction, so fill Highs presolved model with the
        // incumbent model
        presolved_model_ = model_;
      } else if (model_presolve_status_ == HighsPresolveStatus::kReduced ||
                 model_presolve_status_ ==
                     HighsPresolveStatus::kReducedToEmpty) {
        // Nontrivial reduction, so fill Highs presolved model with the
        // presolved model
        using_reduced_lp = true;
      }
      return_status = HighsStatus::kOk;
      break;
    }
    case HighsPresolveStatus::kTimeout: {
      // Timeout, so assume that it's OK to fill the Highs presolved model with
      // the presolved model, but return warning.
      using_reduced_lp = true;
      return_status = HighsStatus::kWarning;
      break;
    }
    default: {
      assert(model_presolve_status_ == HighsPresolveStatus::kOutOfMemory ||
             model_presolve_status_ == HighsPresolveStatus::kException ||
             model_presolve_status_ == HighsPresolveStatus::kMatrixValueError);
      HighsModelStatus model_status;
      if (model_presolve_status_ == HighsPresolveStatus::kOutOfMemory) {
        model_status = HighsModelStatus::kMemoryLimit;
      } else {
        model_status = HighsModelStatus::kPresolveError;
      }
      setHighsModelStatusAndClearSolutionAndBasis(model_status);
      return_status = HighsStatus::kError;
    }
  }
  if (using_reduced_lp) {
    presolved_model_.lp_ = presolve_.getReducedProblem();
    presolved_model_.lp_.setMatrixDimensions();
  }

  highsLogUser(log_options, HighsLogType::kInfo, "Presolve status: %s\n",
               presolveStatusToString(model_presolve_status_).c_str());
  return returnFromHighs(return_status);
}

HighsStatus Highs::run() {
  this->logHeader();
  // Level 0 of Highs::run()
  // Action the file operations associated with running HiGHS, and
  // call Highs::runFromUserScaling()
  HighsStatus status = HighsStatus::kOk;
  if (this->options_.read_solution_file != "")
    status = this->readSolution(this->options_.read_solution_file);
  if (this->options_.read_basis_file != "")
    status = this->readBasis(this->options_.read_basis_file);
  if (this->options_.write_model_file != "")
    status = this->writeModel(this->options_.write_model_file);

  if (status != HighsStatus::kOk) return status;

  this->reportModelStats();

  // Apply any user objective and bound scaling
  HighsUserScaleData user_scale_data;
  if (this->userScale(user_scale_data) == HighsStatus::kError)
    return HighsStatus::kError;

  // Determine coefficient ranges and possibly warn the user about
  // excessive values, obtaining suggested values for user_objective_scale
  // and user_bound_scale
  this->model_.lp_.a_matrix_.ensureColwise();
  assessExcessiveObjectiveBoundScaling(this->options_, this->model_,
                                       user_scale_data);

  status = model_.assess(options_);
  if (status != HighsStatus::kOk) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Highs::run: model_.assess returns status %s\n",
                 this->highsStatusToString(status).c_str());
    return HighsStatus::kError;
  }

  // Optimize the model in the Highs instance
  status = optimizeHighs();
  if (status == HighsStatus::kError) return status;

  // Assume that run data are valid
  this->run_data_.valid = true;

  // Undo any user objective and bound scaling
  if (this->userUnscale(user_scale_data) == HighsStatus::kError)
    return HighsStatus::kError;

  if (this->options_.write_iis_model_file != "")
    status = this->writeIisModel(this->options_.write_iis_model_file);
  if (this->options_.solution_file != "")
    status = this->writeSolution(this->options_.solution_file,
                                 this->options_.write_solution_style);
  if (this->options_.write_basis_file != "")
    status = this->writeBasis(this->options_.write_basis_file);

  return status;
}

HighsStatus Highs::optimizeHighs() {
  this->logHeader();
  // Level 1 of Highs::run()
  //
  // Move the "mods" to here
  return this->multi_linear_objective_.size() ? this->multiobjectiveSolve()
                                              : this->optimizeModel();
}

HighsStatus Highs::optimizeLp() {
  this->logHeader();
  // Solve what's in the HighsLp instance Highs::model_.lp_
  assert(!this->model_.isQp());
  assert(!this->model_.lp_.hasSemiVariables());
  assert(!this->multi_linear_objective_.size());
  return this->optimizeModelTryCatch();
}

HighsStatus Highs::optimizeModel() {
  this->logHeader();
  // Level 2a of Highs::run()
  //
  HighsStatus status = this->initializeMultiThreading();
  if (status != HighsStatus::kOk) return status;

  const bool already_profiling = this->profiling_;
  HighsProfiling profiling;
  if (!already_profiling) this->initializeProfiling(&profiling);
  status = this->optimizeModelTryCatch();
  if (!already_profiling) {
    this->reportProfiling();
    this->clearProfiling();
  }
  return status;
}

HighsStatus Highs::optimizeModelTryCatch() {
  this->logHeader();
  // Wrap the call to Highs::calledOptimizeModel() in a try-catch
  // block
  HighsStatus status;
  auto handleCatch = [&]() {
    // Clear all solver data, since there may be nothing useful
    this->clearSolver();
    status = HighsStatus::kError;
  };
  try {
    status = calledOptimizeModel();
  } catch (const std::exception& exception) {
    model_status_ = handleExceptionIsOom(options_.log_options,
                                         "calledOptimizeModel", exception)
                        ? HighsModelStatus::kMemoryLimit
                        : HighsModelStatus::kSolveError;
    handleCatch();
  } catch (const HighsTask::Interrupt&) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "HighsTask interrupt in calledOptimizeModel\n");
    model_status_ = HighsModelStatus::kSolveError;
    handleCatch();
  }
  return status;
}

// Checks the options calls presolve and postsolve if needed.
//
// LP solvers are called with callSolveLp(..)
HighsStatus Highs::calledOptimizeModel() {
  this->logHeader();
  // Level 2b of Highs::run()
  //
  HighsInt min_highs_debug_level = kHighsDebugLevelMin;
  // kHighsDebugLevelCostly;
  // kHighsDebugLevelMax;
  if (options_.highs_debug_level < min_highs_debug_level)
    options_.highs_debug_level = min_highs_debug_level;

  const bool possibly_use_log_dev_level_2 = false;
  const HighsInt log_dev_level = options_.log_dev_level;
  const bool output_flag = options_.output_flag;
  HighsInt use_log_dev_level = log_dev_level;
  bool use_output_flag = output_flag;
  const HighsInt check_debug_run_call_num = -103757;
  const HighsInt check_num_col = -317;
  const HighsInt check_num_row = -714;
  if (possibly_use_log_dev_level_2) {
    if (this->debug_optimize_call_num_ == check_debug_run_call_num &&
        model_.lp_.num_col_ == check_num_col &&
        model_.lp_.num_row_ == check_num_row) {
      std::string message =
          "Entering Highs::optimizeModel(): run/col/row matching check ";
      highsLogDev(options_.log_options, HighsLogType::kInfo,
                  "%s: run %d: LP(%6d, %6d)\n", message.c_str(),
                  int(this->debug_optimize_call_num_), int(model_.lp_.num_col_),
                  int(model_.lp_.num_row_));
      // highsPause(true, message);
      use_log_dev_level = 2;
      use_output_flag = true;
    }
  }

  // As exposed by #1819, retaining the name hashes generated when
  // adding names using addColName/addRowName can lead to serious
  // performance regression, so clear them now
  this->model_.lp_.clearAllNameHash();

  if (!options_.use_warm_start) this->clearSolver();
  if (ekk_instance_.status_.has_nla)
    assert(ekk_instance_.lpFactorRowCompatible(model_.lp_.num_row_));

  // returnFromOptimizeModel() is a common exit method to ensure
  // consistency of values set by optimizeModel() and many other
  // things. It's important to be able to check that it's been called,
  // and this is done with this->called_return_from_optimize_model
  //
  // Make sure here that returnFromOptimizeModel() has been called
  // after any previous call to optimizeModel()

  assert(called_return_from_optimize_model);
  if (!called_return_from_optimize_model) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "Highs::optimizeModel() called with "
                "called_return_from_optimize_model false\n");
    return HighsStatus::kError;
  }

  // HiGHS solvers require models with no infinite costs, and no semi-variables
  //
  // Since completeSolutionFromDiscreteAssignment() may require a call
  // to optimizeModel() - with initial check that
  // called_return_from_optimize_model is true
  // - called_return_from_optimize_model cannot yet be set false.
  //
  // This possible call to optimizeModel() means that any need to modify the
  // problem to remove infinite costs must be done first.
  //
  // Set undo_mods = false so that returnFromOptimizeModel() doesn't undo any
  // mods that must be preserved - such as when solving a MIP node
  bool undo_mods = false;
  if (model_.lp_.has_infinite_cost_) {
    // If the model has infinite costs, then try to remove them. The
    // return_status indicates the success of this operation and, if
    // it's unsuccessful, the model will not have been modified and
    // optimizeModel() can simply return an error with model status
    // HighsModelStatus::kUnknown
    assert(model_.lp_.hasInfiniteCost(options_.infinite_cost));
    HighsStatus return_status = handleInfCost();
    if (return_status != HighsStatus::kOk) {
      assert(return_status == HighsStatus::kError);
      setHighsModelStatusAndClearSolutionAndBasis(HighsModelStatus::kUnknown);
      return return_status;
    }
    // Modifications have been performed, so must be undone before
    // this call to optimizeModel() returns
    assert(!model_.lp_.has_infinite_cost_);
    undo_mods = true;
  } else {
    assert(!model_.lp_.hasInfiniteCost(options_.infinite_cost));
  }

  // Ensure that all vectors in the model have exactly the right size
  exactResizeModel();

  if (model_.isMip() && solution_.value_valid) {
    // Determine whether the current (partial) solution of a MIP is
    // feasible and, if not, try to complete the assignment with
    // integer values (if necessary) and continuous values (if
    // necessary) to achieve a feasible solution. Valuable in the case
    // where users make a heuristic (partial) assignment of discrete variables
    HighsStatus call_status = completeSolutionFromDiscreteAssignment();
    if (call_status != HighsStatus::kOk) return HighsStatus::kError;
  }

  // Set this so that calls to returnFromOptimizeModel() can be checked: from
  // here all return statements execute returnFromOptimizeModel()
  called_return_from_optimize_model = false;
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Initialise the HiGHS model status
  model_status_ = HighsModelStatus::kNotset;
  // Clear the run info and data
  invalidateInfo();
  invalidateRunData();
  // Zero the iteration counts
  zeroIterationCounts();
  // Start the HiGHS run clock
  timer_.start();
  // Return immediately if the model has no columns
  if (!model_.lp_.num_col_) {
    setHighsModelStatusAndClearSolutionAndBasis(HighsModelStatus::kModelEmpty);
    return returnFromOptimizeModel(HighsStatus::kOk, undo_mods);
  }
  // Process the bounds, performing inward rounding of bounds on any
  // integer variables and modifying any bounds with inconsistencies
  // less than the primal feasibility tolerance, but returning
  // infeasibility if the model is has any bounds with larger
  // inconsistencies
  if (processedBoundsFindsInfeasibility()) {
    setHighsModelStatusAndClearSolutionAndBasis(HighsModelStatus::kInfeasible);
    return returnFromOptimizeModel(return_status, undo_mods);
  }
  // Ensure that the HighsLp has the matrix column-wise
  model_.lp_.ensureColwise();
  // Shouldn't have to check validity of the model since this is done
  // when it is loaded or modified
  assert((model_.assess(options_) == HighsStatus::kOk) &&
         "model_.assess not HighsStatus::kOk");

  if (options_.highs_debug_level > min_highs_debug_level) {
    // Shouldn't have to check validity of the LP since this is done when it is
    // loaded or modified
    call_status = model_.lp_.assess(options_);
    // If any errors have been found or normalisation carried out,
    // call_status will be kError or kWarning, so only valid return is OK.
    assert(call_status == HighsStatus::kOk);
    return_status = interpretCallStatus(options_.log_options, call_status,
                                        return_status, "assessLp");
    if (return_status == HighsStatus::kError)
      return returnFromOptimizeModel(return_status, undo_mods);
    // Shouldn't have to check that the options settings are legal,
    // since they are checked when modified
    if (checkOptions(options_.log_options, options_.records) !=
        OptionStatus::kOk) {
      return_status = HighsStatus::kError;
      return returnFromOptimizeModel(return_status, undo_mods);
    }
  }

  if (model_.lp_.model_name_.compare(""))
    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Solving model: %s\n", model_.lp_.model_name_.c_str());

  if (!options_.solve_relaxation) {
    // Not solving the relaxation, so check validity of any
    // integrality, keeping a record of any bound and type
    // modifications for semi-variables
    bool made_semi_variable_mods = false;
    call_status =
        assessSemiVariables(model_.lp_, options_, made_semi_variable_mods);
    undo_mods = undo_mods || made_semi_variable_mods;
    if (call_status == HighsStatus::kError) {
      setHighsModelStatusAndClearSolutionAndBasis(
          HighsModelStatus::kSolveError);
      return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
    }
  }
  // Choose method according to model class
  if (model_.isQp()) {
    if (model_.isMip()) {
      if (options_.solve_relaxation) {
        // Relax any semi-variables
        bool made_semi_variable_mods = false;
        relaxSemiVariables(model_.lp_, made_semi_variable_mods);
        undo_mods = undo_mods || made_semi_variable_mods;
      } else {
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Cannot solve MIQP problems with HiGHS\n");
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
    }
    // Ensure that its diagonal entries are OK in the context of the
    // objective sense. It's OK to be semi-definite
    if (!okHessianDiagonal(options_, model_.hessian_, model_.lp_.sense_)) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Cannot solve non-convex QP problems with HiGHS\n");
      return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
    }
    if (!solverValidForQp(options_.solver)) warnSolverInvalid(options_, "QP");
    try {
      call_status = callSolveQp(this->model_, "Solve incumbent QP");
    } catch (const std::exception& exception) {
      model_status_ =
          handleExceptionIsOom(options_.log_options, "callSolveQp", exception)
              ? HighsModelStatus::kMemoryLimit
              : HighsModelStatus::kSolveError;
      call_status = HighsStatus::kError;
    }
    return_status = interpretCallStatus(options_.log_options, call_status,
                                        return_status, "callSolveQp");
    return returnFromOptimizeModel(return_status, undo_mods);
  } else if (model_.isMip()) {
    // Model is MIP
    if (options_.solve_relaxation) {
      // Solving just the relaxation, so relax any semi-variables
      bool made_semi_variable_mods = false;
      relaxSemiVariables(model_.lp_, made_semi_variable_mods);
      undo_mods = undo_mods || made_semi_variable_mods;
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Solving LP relaxation since solve_relaxation is true\n");
    } else {
      // Solve model as a MIP
      if (!solverValidForMip(options_.solver))
        warnSolverInvalid(options_, "MIP");
      call_status = callSolveMip(model_.lp_, "Solve incumbent MIP");
      return_status = interpretCallStatus(options_.log_options, call_status,
                                          return_status, "callSolveMip");
      return returnFromOptimizeModel(return_status, undo_mods);
    }
  }
  // Solve the model as an LP
  HighsLp& incumbent_lp = model_.lp_;
  HighsLogOptions& log_options = options_.log_options;
  bool no_incumbent_lp_solution_or_basis = false;

  if (!solverValidForLp(options_.solver)) warnSolverInvalid(options_, "LP");

  //
  // Record the initial time and set the component times and postsolve
  // iteration count to -1 to identify whether they are not required
  double initial_time = timer_.read();
  double this_presolve_time = -1;
  double this_solve_presolved_lp_time = -1;
  double this_postsolve_time = -1;
  double this_solve_original_lp_time = -1;
  HighsInt postsolve_iteration_count = -1;
  const bool may_require_basis_postsolve = mayRequireBasisPostsolve(&options_);
  if (options_.icrash) {
    ICrashStrategy strategy = ICrashStrategy::kICA;
    bool strategy_ok = parseICrashStrategy(options_.icrash_strategy, strategy);
    if (!strategy_ok) {
      // std::cout << "ICrash error: unknown strategy." << std::endl;
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "ICrash error: unknown strategy.\n");
      return HighsStatus::kError;
    }
    ICrashOptions icrash_options{
        options_.icrash_dualize,         strategy,
        options_.icrash_starting_weight, options_.icrash_iterations,
        options_.icrash_approx_iter,     options_.icrash_exact,
        options_.icrash_breakpoints,     options_.log_options};

    HighsStatus icrash_status =
        callICrash(model_.lp_, icrash_options, icrash_info_);

    if (icrash_status != HighsStatus::kOk)
      return returnFromOptimizeModel(icrash_status, undo_mods);

    // for now set the solution_.col_value
    solution_.col_value = icrash_info_.x_values;
    // Better not to use Highs::crossover
    const bool use_highs_crossover = false;
    if (use_highs_crossover) {
      crossover(solution_);
      // loops:
      called_return_from_optimize_model = true;

      options_.icrash = false;  // to avoid loop
    } else {
      HighsStatus crossover_status =
          callCrossover(options_, model_.lp_, basis_, solution_, model_status_,
                        info_, callback_);
      // callCrossover can return HighsStatus::kWarning due to
      // imprecise dual values. Ignore this since primal simplex will
      // be called to clean up duals
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Crossover following iCrash has return status of %s, and "
                   "problem status is %s\n",
                   highsStatusToString(crossover_status).c_str(),
                   modelStatusToString(model_status_).c_str());
      if (crossover_status == HighsStatus::kError)
        return returnFromOptimizeModel(crossover_status, undo_mods);
      assert(options_.simplex_strategy == kSimplexStrategyPrimal);
    }
    // timer_.stop();
    // run();

    // todo: add "dual" values
    // return HighsStatus::kOk;
  }

  // Even if options_.solver == kHighsChooseString in isolation will,
  // ultimately lead to a choice between simplex and IPM, if a basis
  // is available, simplex should surely be chosen.
  const bool solver_will_use_basis = options_.solver == kSimplexString ||
                                     options_.solver == kHighsChooseString;

  if (solver_will_use_basis) {
    if (!basis_.valid && solution_.value_valid) {
      // There is no valid basis, but there is a valid solution, so use
      // it to construct a basis
      return_status =
          interpretCallStatus(options_.log_options, basisForSolution(),
                              return_status, "basisForSolution");
      if (return_status == HighsStatus::kError)
        return returnFromOptimizeModel(return_status, undo_mods);
      assert(basis_.valid);
    }
  } else {
    // The basis won't be used, so clear it to ensure that, after any
    // presolve the solver choice won't be over-ruled by choosing
    // simplex due to the existence of a basis - which mustn't be used
    // after a strict reduction due to presolve!
    basis_.clear();
  }

  // lambda for Lp solving
  auto timedSolveLp = [&](HighsLp& lp, const std::string& lpSolveDescription,
                          double& time) {
    time = -timer_.read(timer_.solve_clock);
    if (possibly_use_log_dev_level_2) {
      options_.log_dev_level = use_log_dev_level;
      options_.output_flag = use_output_flag;
    }
    timer_.start(timer_.solve_clock);
    const bool as_mip = false;
    call_status = callSolveLp(lp, lpSolveDescription);
    timer_.stop(timer_.solve_clock);
    if (possibly_use_log_dev_level_2) {
      options_.log_dev_level = log_dev_level;
      options_.output_flag = output_flag;
    }
    time += timer_.read(timer_.solve_clock);
  };

  const bool unconstrained_lp = incumbent_lp.numNz() == 0;
  assert(incumbent_lp.num_row_ || unconstrained_lp);
  const bool has_basis = basis_.useful;
  if (has_basis) {
    assert(basis_.col_status.size() ==
           static_cast<size_t>(incumbent_lp.num_col_));
    assert(basis_.row_status.size() ==
           static_cast<size_t>(incumbent_lp.num_row_));
  }
  if (basis_.valid) assert(basis_.useful);

  const bool without_presolve = options_.presolve == kHighsOffString;
  if ((has_basis && solver_will_use_basis) || without_presolve ||
      unconstrained_lp) {
    // There is a valid basis for the problem and the solver will use
    // it, or presolve is off, or LP has no constraint matrix
    //
    // Determine a coherent message about how the LP is being solved
    std::stringstream lp_solve_ss;
    if (unconstrained_lp) {
      lp_solve_ss << "Solving unconstrained LP";
    } else if (has_basis) {
      if (without_presolve) {
        lp_solve_ss << "Solving LP with useful basis";
      } else {
        lp_solve_ss << "Solving LP with useful basis so presolve not used";
      }
    } else {
      // One of unconstrained_lp, has_basis and without_presolve must
      // be true, and the first two aren't
      assert(without_presolve);
      lp_solve_ss << "Solving LP without presolve or useful basis";
    }
    std::string lp_solve = lp_solve_ss.str();
    ekk_instance_.lp_name_ = lp_solve;
    // If there is a valid HiGHS basis, refine any status values that
    // are simply HighsBasisStatus::kNonbasic
    if (basis_.useful) refineBasis(incumbent_lp, solution_, basis_);
    timedSolveLp(incumbent_lp, lp_solve, this_solve_original_lp_time);
    return_status = interpretCallStatus(options_.log_options, call_status,
                                        return_status, "callSolveLp");
    this->run_data_.solve_time = this_solve_original_lp_time;
    if (return_status == HighsStatus::kError)
      return returnFromOptimizeModel(return_status, undo_mods);
  } else {
    // Possibly presolve - according to option_.presolve
    //
    // If solving the relaxation of a MIP, make sure that LP presolve
    // is used - so that rules assuming MIP properties are not applied.
    const double from_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time = -from_presolve_time;
    timer_.start(timer_.presolve_clock);
    const bool force_lp_presolve = true;
    model_presolve_status_ = runPresolve(force_lp_presolve);
    timer_.stop(timer_.presolve_clock);
    const double to_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time += to_presolve_time;
    presolve_.info_.presolve_time = this_presolve_time;
    this->run_data_.presolve_time = this_presolve_time;

    // Set an illegal local pivot threshold value that's updated after
    // solving the presolved LP - if simplex is used
    double factor_pivot_threshold = -1;

    // Run solver.
    bool have_optimal_solution = false;
    // ToDo Put solution of presolved problem in a separate method

    // Record the number of iterations when PDLP is used to solve the
    // presolved problem, since the iteration count is reset to zero
    // if PDLP is used to clean up after postsolve
    HighsInt presolved_lp_pdlp_iteration_count = 0;
    // Log the presolve reductions
    reportPresolveReductions(log_options, model_presolve_status_, incumbent_lp,
                             presolve_.getReducedProblem());
    this->run_data_.presolved_model_num_col =
        presolve_.getReducedProblem().num_col_;
    this->run_data_.presolved_model_num_row =
        presolve_.getReducedProblem().num_row_;
    this->run_data_.presolved_model_num_nz =
        presolve_.getReducedProblem().a_matrix_.numNz();
    switch (model_presolve_status_) {
      case HighsPresolveStatus::kNotPresolved: {
        ekk_instance_.lp_name_ = "Original LP";
        timedSolveLp(incumbent_lp, "Not presolved: solving the LP",
                     this_solve_original_lp_time);
        this->run_data_.solve_time = this_solve_original_lp_time;
        return_status = interpretCallStatus(options_.log_options, call_status,
                                            return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromOptimizeModel(return_status, undo_mods);
        break;
      }
      case HighsPresolveStatus::kNotReduced: {
        ekk_instance_.lp_name_ = "Unreduced LP";
        // Log the presolve reductions
        timedSolveLp(incumbent_lp,
                     "Problem not reduced by presolve: solving the LP",
                     this_solve_original_lp_time);
        this->run_data_.solve_time = this_solve_original_lp_time;
        return_status = interpretCallStatus(options_.log_options, call_status,
                                            return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromOptimizeModel(return_status, undo_mods);
        break;
      }
      case HighsPresolveStatus::kReduced: {
        HighsLp& reduced_lp = presolve_.getReducedProblem();
        reduced_lp.origin_name_ = "Reduced LP";
        reduced_lp.setMatrixDimensions();
        if (kAllowDeveloperAssert) {
          // Validate the reduced LP
          //
          // Although presolve can yield small values in the matrix,
          // they are only stripped out (by assessLp) in debug. This
          // suggests that they are no real danger to the simplex
          // solver. The only danger is pivoting on them, but that
          // implies that values of roughly that size have been chosen
          // in the ratio test. Even with the filter, values of 1e-9
          // could be in the matrix, and these would be bad
          // pivots. Hence, since the small values may play a
          // meaningful role in postsolve, then it's better to keep
          // them.
          //
          // ToDo. Analyse the extent of small value creation. See #1187
          assert(assessLp(reduced_lp, options_) == HighsStatus::kOk);
        } else {
          // Just analyse extreme values: none is changed
          reduced_lp.a_matrix_.assessExtremeValues("of presolved problem ",
                                                   options_);
        }
        call_status = cleanBounds(options_, reduced_lp);
        // Ignore any warning from clean bounds since the original LP
        // is still solved after presolve
        if (interpretCallStatus(options_.log_options, call_status,
                                return_status,
                                "cleanBounds") == HighsStatus::kError)
          return HighsStatus::kError;
        // Solving the presolved LP with strictly reduced dimensions
        // so ensure that the Ekk instance is cleared
        ekk_instance_.clear();
        ekk_instance_.lp_name_ = "Presolved LP";
        // Don't try dual cut-off when solving the presolved LP, as the
        // objective values aren't correct
        const double save_objective_bound = options_.objective_bound;
        options_.objective_bound = kHighsInf;
        timedSolveLp(reduced_lp, "Solving the presolved LP",
                     this_solve_presolved_lp_time);
        this->run_data_.solve_time = this_solve_presolved_lp_time;
        if (ekk_instance_.status_.initialised_for_solve) {
          // Record the pivot threshold resulting from solving the presolved LP
          // with simplex
          factor_pivot_threshold = ekk_instance_.info_.factor_pivot_threshold;
        }
        // Restore the dual objective cut-off
        options_.objective_bound = save_objective_bound;
        return_status = interpretCallStatus(options_.log_options, call_status,
                                            return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromOptimizeModel(return_status, undo_mods);
        presolved_lp_pdlp_iteration_count = info_.pdlp_iteration_count;
        have_optimal_solution = model_status_ == HighsModelStatus::kOptimal;
        no_incumbent_lp_solution_or_basis =
            model_status_ == HighsModelStatus::kInfeasible ||
            model_status_ == HighsModelStatus::kUnbounded ||
            model_status_ == HighsModelStatus::kUnboundedOrInfeasible ||
            model_status_ == HighsModelStatus::kTimeLimit ||
            model_status_ == HighsModelStatus::kIterationLimit ||
            model_status_ == HighsModelStatus::kInterrupt;
        if (no_incumbent_lp_solution_or_basis) {
          // Postsolve won't be performed, so clear the HEkk data
          // corresponding to the (strictly reduced) presolved LP here
          ekk_instance_.clear();
          setHighsModelStatusAndClearSolutionAndBasis(model_status_);
        }
        break;
      }
      case HighsPresolveStatus::kReducedToEmpty: {
        // Create a trivial optimal solution for postsolve to use
        solution_.clear();
        basis_.clear();
        basis_.debug_origin_name = "Presolve to empty";
        basis_.valid = true;
        basis_.alien = false;
        basis_.useful = true;
        basis_.was_alien = false;
        solution_.value_valid = true;
        solution_.dual_valid = true;
        have_optimal_solution = true;
        // Optimality will not be spotted if there's no basis
        // postsolve
        model_status_ = HighsModelStatus::kOptimal;
        this->run_data_.solve_time = 0;
        break;
      }
      case HighsPresolveStatus::kInfeasible: {
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kInfeasible);
        highsLogUser(log_options, HighsLogType::kInfo,
                     "Problem status detected on presolve: %s\n",
                     modelStatusToString(model_status_).c_str());
        return returnFromOptimizeModel(return_status, undo_mods);
      }
      case HighsPresolveStatus::kUnboundedOrInfeasible: {
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Problem status detected on presolve: %s\n",
            modelStatusToString(HighsModelStatus::kUnboundedOrInfeasible)
                .c_str());
        if (options_.allow_unbounded_or_infeasible) {
          setHighsModelStatusAndClearSolutionAndBasis(
              HighsModelStatus::kUnboundedOrInfeasible);
          return returnFromOptimizeModel(return_status, undo_mods);
        }
        // Presolve has returned kUnboundedOrInfeasible, but HiGHS
        // can't return this. Use primal simplex solver on the original
        // LP
        HighsStatus save_return_status = return_status;
        HighsOptions save_options = options_;
        options_.solver = "simplex";
        options_.simplex_strategy = kSimplexStrategyPrimal;
        timedSolveLp(incumbent_lp,
                     "Solving the original LP with primal simplex "
                     "to determine infeasible or unbounded",
                     this_solve_original_lp_time);
        // Recover the options
        options_ = save_options;
        // Primal simplex doesn't necessarily resolve the
        // kUnboundedOrInfeasible case, as it may time out or reach
        // the iteration count (#3383), or even fail
        assert(model_status_ == HighsModelStatus::kInfeasible ||
               model_status_ == HighsModelStatus::kUnbounded ||
               model_status_ == HighsModelStatus::kTimeLimit ||
               model_status_ == HighsModelStatus::kIterationLimit ||
               model_status_ == HighsModelStatus::kInterrupt);
        const bool resolved = model_status_ == HighsModelStatus::kInfeasible ||
                              model_status_ == HighsModelStatus::kUnbounded;
        if (resolved) {
          info_.valid = true;
        } else {
          // Unresolved, so return with result from presolve
          highsLogUser(log_options, HighsLogType::kInfo,
                       "Unable to resolve whether LP is unbounded or "
                       "infeasible: primal simplex yields model status of %s\n",
                       modelStatusToString(model_status_).c_str());

          model_status_ = HighsModelStatus::kUnboundedOrInfeasible;
          return_status = save_return_status;
          setHighsModelStatusAndClearSolutionAndBasis(
              HighsModelStatus::kUnboundedOrInfeasible);
        }
        return returnFromOptimizeModel(return_status, undo_mods);
      }
      case HighsPresolveStatus::kTimeout: {
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kTimeLimit);
        highsLogDev(log_options, HighsLogType::kWarning,
                    "Presolve reached timeout\n");
        return returnFromOptimizeModel(HighsStatus::kWarning, undo_mods);
      }
      case HighsPresolveStatus::kOutOfMemory: {
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kMemoryLimit);
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
      case HighsPresolveStatus::kException: {
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kSolveError);
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
      case HighsPresolveStatus::kMatrixValueError: {
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kPresolveError);
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
      default: {
        assert(model_presolve_status_ == HighsPresolveStatus::kNullError);
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kPresolveError);
        highsLogDev(log_options, HighsLogType::kError,
                    "Presolve returned status %s\n",
                    presolveStatusToString(model_presolve_status_).c_str());
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
    }
    // End of presolve
    //
    // Cases of infeasibility/unboundedness timeout and memory errors
    // all handled, so just the successes remain
    assert(model_presolve_status_ == HighsPresolveStatus::kNotPresolved ||
           model_presolve_status_ == HighsPresolveStatus::kNotReduced ||
           model_presolve_status_ == HighsPresolveStatus::kReduced ||
           model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty);

    // Postsolve. Does nothing if there were no reductions during presolve.

    // If presolve has been run assuming that there's no basis
    // postsolve - allowing sparsify to be used in presolve -
    // invalidate any basis
    if (!may_require_basis_postsolve) this->invalidateBasis();
    const bool have_optimal_reduced_solution =
        model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty ||
        (model_presolve_status_ == HighsPresolveStatus::kReduced &&
         model_status_ == HighsModelStatus::kOptimal);
    const bool have_unknown_reduced_solution =
        model_presolve_status_ == HighsPresolveStatus::kReduced &&
        model_status_ == HighsModelStatus::kUnknown;
    // Take a copy of the HighsInfo and HighsModelStatus for the
    // presolved LP to determine whether a PDLP clean-up should be
    // carried out
    const HighsInfo presolved_lp_info = this->info_;
    const HighsModelStatus presolved_lp_model_status = this->model_status_;
    if (have_optimal_reduced_solution || have_unknown_reduced_solution) {
      // ToDo Put this in a separate method
      assert(model_status_ == HighsModelStatus::kOptimal ||
             model_status_ == HighsModelStatus::kUnknown ||
             model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty);
      assert(model_presolve_status_ == HighsPresolveStatus::kReduced ||
             model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty);
      if (have_unknown_reduced_solution)
        highsLogUser(
            options_.log_options, HighsLogType::kWarning,
            "Running postsolve on non-optimal solution of reduced LP\n\n");
      // If presolve is nontrivial, extract the optimal solution
      // and (any) basis for the presolved problem in order to generate
      // the solution and basis for postsolve to use to generate a
      // solution and basis that is, hopefully, optimal.
      //
      // If there is a basis, optimality is confirmed or corrected by
      // hot-starting the simplex solver
      //
      // If there is no basis, then a KKT check is carried out, and
      // used to determine whether a PDLP clean-up should be done
      presolve_.data_.recovered_solution_ = solution_;
      presolve_.data_.recovered_basis_ = basis_;

      if (model_presolve_status_ == HighsPresolveStatus::kReduced)
        this->callLpKktCheck(presolve_.getReducedProblem(), "Before postsolve");

      this_postsolve_time = -timer_.read(timer_.postsolve_clock);
      timer_.start(timer_.postsolve_clock);
      HighsPostsolveStatus postsolve_status = runPostsolve();
      timer_.stop(timer_.postsolve_clock);
      this_postsolve_time += timer_.read(timer_.postsolve_clock);
      presolve_.info_.postsolve_time = this_postsolve_time;
      this->run_data_.postsolve_time = this_postsolve_time;

      if (postsolve_status == HighsPostsolveStatus::kSolutionRecovered) {
        // Indicate that nontrivial postsolve has been performed
        if (model_presolve_status_ == HighsPresolveStatus::kReduced ||
            model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty)
          highsLogUser(log_options, HighsLogType::kInfo,
                       "Performed postsolve\n");
        // Set solution and its status
        solution_.clear();
        solution_ = presolve_.data_.recovered_solution_;
        solution_.value_valid = true;
        if (!basis_.valid) {
          // Have a primal-dual solution, but no basis, since PDLP or
          // IPX without crossover were used. In the case of IPX, this
          // was because either run_crossover was "off" or "choose"
          // and IPX determined optimality
          solution_.dual_valid = true;
          this->invalidateBasis();
        } else {
          //
          // Hot-start the simplex solver for the incumbent LP
          //
          solution_.dual_valid = true;
          // Set basis and its status
          basis_.valid = true;
          basis_.useful = true;
          basis_.col_status = presolve_.data_.recovered_basis_.col_status;
          basis_.row_status = presolve_.data_.recovered_basis_.row_status;
          basis_.debug_origin_name += ": after postsolve";
          // Basic primal activities are wrong after postsolve, so
          // possibly skip KKT check
          const bool perform_kkt_check = true;
          if (perform_kkt_check) {
            // Possibly force debug to perform KKT check on what's
            // returned from postsolve
            const bool force_debug = false;
            HighsInt save_highs_debug_level = options_.highs_debug_level;
            if (force_debug)
              options_.highs_debug_level = kHighsDebugLevelCostly;
            if (debugHighsSolution("After returning from postsolve", options_,
                                   model_, solution_,
                                   basis_) == HighsDebugStatus::kLogicalError)
              return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
            options_.highs_debug_level = save_highs_debug_level;
          }
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions save_options = options_;
          const bool full_logging = false;
          if (full_logging) options_.log_dev_level = kHighsLogDevLevelVerbose;
          // Force the use of simplex to clean up
          options_.solver = kSimplexString;
          options_.simplex_strategy = kSimplexStrategyChoose;
          // Ensure that the parallel solver isn't used
          options_.simplex_min_concurrency = 1;
          options_.simplex_max_concurrency = 1;
          // Use any pivot threshold resulting from solving the presolved LP
          if (factor_pivot_threshold > 0)
            options_.factor_pivot_threshold = factor_pivot_threshold;
          // The basis returned from postsolve is just basic/nonbasic
          // and EKK expects a refined basis, so set it up now
          refineBasis(incumbent_lp, solution_, basis_);
          // Scrap the EKK data from solving the presolved LP
          ekk_instance_.invalidate();
          ekk_instance_.lp_name_ = "Postsolve LP";
          // Set up the iteration count and timing records so that
          // adding the corresponding values after callSolveLp gives
          // difference
          postsolve_iteration_count = -info_.simplex_iteration_count;
          timedSolveLp(
              incumbent_lp,
              "Solving the original LP from the solution after postsolve",
              this_solve_original_lp_time);
          this->run_data_.postsolve_time += this_solve_original_lp_time;
          // Determine the iteration count
          postsolve_iteration_count += info_.simplex_iteration_count;
          //
          return_status = HighsStatus::kOk;
          return_status = interpretCallStatus(log_options, call_status,
                                              return_status, "callSolveLp");
          // Recover the options
          options_ = save_options;
          if (return_status == HighsStatus::kError)
            return returnFromOptimizeModel(return_status, undo_mods);
          this->run_data_.num_simplex_iterations_after_postsolve =
              postsolve_iteration_count;
          if (postsolve_iteration_count > 0)
            highsLogUser(log_options, HighsLogType::kInfo,
                         "Required %d simplex iterations after postsolve\n",
                         int(postsolve_iteration_count));
        }
      } else {
        // Postsolve has failed, so any simplex data for the presolved
        // problem is useless
        ekk_instance_.clear();
        highsLogUser(log_options, HighsLogType::kError,
                     "Postsolve return status: %s\n",
                     utilPostsolveStatusToString(postsolve_status).c_str());
        setHighsModelStatusAndClearSolutionAndBasis(
            HighsModelStatus::kPostsolveError);
        return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
      }
    } else {
      this->run_data_.postsolve_time = 0;
    }
    // Final chance to resolve model_status_ ==
    // HighsModelStatus::kUnknown after solver that doesn't yield a
    // basis
    //
    // Possibly clean up using PDLP
    const bool consider_pdlp_cleanup = false;
    /*
      !basis_.valid && model_status_ == HighsModelStatus::kUnknown &&
      options_.allow_pdlp_cleanup && !options_.run_centring;
    */
    if (consider_pdlp_cleanup) {
      HighsInt pdlp_cleanup_iteration_limit = 0;
      if (tryPdlpCleanup(pdlp_cleanup_iteration_limit, presolved_lp_info)) {
        // Primal/dual infeasibilities/residuals can be magnified in
        // postsolve after PDLP, and IPX without crossover can fail,
        // both leading to model_status_ == HighsModelStatus::kUnknown,
        // so run PDLP on the original problem, starting from the
        // incumbent solution
        const std::string solver = options_.solver;
        const HighsInt pdlp_iteration_limit = options_.pdlp_iteration_limit;
        assert(useIpm(solver) || solver == kPdlpString);
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Unknown model status and no basis after initial solve "
            "so use PDLP with KKT tolerance = %g and iteration limit = %d to "
            "solve "
            "the original LP from the incumbent solution after postsolve\n",
            options_.kkt_tolerance, int(pdlp_cleanup_iteration_limit));
        options_.solver = kPdlpString;
        options_.pdlp_iteration_limit = pdlp_cleanup_iteration_limit;
        timedSolveLp(
            incumbent_lp,
            "Using PDLP to solve the original LP from the solution after "
            "postsolve",
            this_solve_original_lp_time);
        // Recover solver and time limit option values
        options_.solver = solver;
        options_.pdlp_iteration_limit = pdlp_iteration_limit;
        return_status = HighsStatus::kOk;
        return_status = interpretCallStatus(log_options, call_status,
                                            return_status, "callSolveLp");
        if (return_status == HighsStatus::kError)
          return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
        // Update the number of PDLP, iterations, since the
        // iteration count is reset to zero if PDLP is used to
        // clean up after postsolve
        info_.pdlp_iteration_count += presolved_lp_pdlp_iteration_count;
      }
    }
  }
  // Cycling can yield model_status_ == HighsModelStatus::kNotset,
  //  assert(model_status_ != HighsModelStatus::kNotset);

  // Unless the model status was determined using the strictly reduced LP, the
  // HiGHS info is valid
  if (!no_incumbent_lp_solution_or_basis) {
    this->callLpKktCheck(this->model_.lp_);
    info_.valid = true;
  }

  double lp_solve_final_time = timer_.read();
  double this_solve_time = lp_solve_final_time - initial_time;
  if (postsolve_iteration_count < 0) {
    highsLogDev(log_options, HighsLogType::kInfo, "Postsolve  : \n");
  } else {
    highsLogDev(log_options, HighsLogType::kInfo,
                "Postsolve  : %" HIGHSINT_FORMAT "\n",
                postsolve_iteration_count);
  }
  if (this_solve_time > 0)
    highsLogDev(log_options, HighsLogType::kInfo, "Time           : %8.2f\n",
                this_solve_time);
  if (this_presolve_time > 0)
    highsLogDev(log_options, HighsLogType::kInfo, "Time Pre       : %8.2f\n",
                this_presolve_time);
  if (this_solve_presolved_lp_time > 0)
    highsLogDev(log_options, HighsLogType::kInfo, "Time PreLP     : %8.2f\n",
                this_solve_presolved_lp_time);
  if (this_solve_original_lp_time > 0)
    highsLogDev(log_options, HighsLogType::kInfo, "Time OriginalLP: %8.2f\n",
                this_solve_original_lp_time);
  if (this_solve_time > 0) {
    highsLogDev(log_options, HighsLogType::kInfo, "For LP %16s",
                incumbent_lp.model_name_.c_str());
    double sum_time = 0;
    if (this_presolve_time > 0) {
      sum_time += this_presolve_time;
      HighsInt pct = (100 * this_presolve_time) / this_solve_time;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "    : Presolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_presolve_time, pct);
    }
    if (this_solve_presolved_lp_time > 0) {
      sum_time += this_solve_presolved_lp_time;
      HighsInt pct = (100 * this_solve_presolved_lp_time) / this_solve_time;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "    : Solve presolved LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_presolved_lp_time, pct);
    }
    if (this_postsolve_time > 0) {
      sum_time += this_postsolve_time;
      HighsInt pct = (100 * this_postsolve_time) / this_solve_time;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "    : Postsolve %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_postsolve_time, pct);
    }
    if (this_solve_original_lp_time > 0) {
      sum_time += this_solve_original_lp_time;
      HighsInt pct = (100 * this_solve_original_lp_time) / this_solve_time;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "    : Solve original LP %8.2f (%3" HIGHSINT_FORMAT "%%)",
                  this_solve_original_lp_time, pct);
    }
    highsLogDev(log_options, HighsLogType::kInfo, "\n");
    double rlv_time_difference =
        fabs(sum_time - this_solve_time) / this_solve_time;
    if (rlv_time_difference > 0.1) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Strange: Solve time = %g; Sum times = %g: relative "
                  "difference = %g\n",
                  this_solve_time, sum_time, rlv_time_difference);
    }
  }
  // Assess success according to the model status, regardless of
  // whether anything worse has happened earlier
  return_status = highsStatusFromHighsModelStatus(model_status_);
  return returnFromOptimizeModel(return_status, undo_mods);
}

HighsStatus Highs::getStandardFormLp(HighsInt& num_col, HighsInt& num_row,
                                     HighsInt& num_nz, double& offset,
                                     double* cost, double* rhs, HighsInt* start,
                                     HighsInt* index, double* value) {
  this->logHeader();
  if (!this->standard_form_valid_) {
    HighsStatus status = formStandardFormLp();
    assert(status == HighsStatus::kOk);
    if (status != HighsStatus::kOk) return status;
  }
  num_col = this->standard_form_cost_.size();
  num_row = this->standard_form_rhs_.size();
  num_nz = this->standard_form_matrix_.start_[num_col];
  offset = this->standard_form_offset_;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (cost) cost[iCol] = this->standard_form_cost_[iCol];
    if (start) start[iCol] = this->standard_form_matrix_.start_[iCol];
    if (index || value) {
      for (HighsInt iEl = this->standard_form_matrix_.start_[iCol];
           iEl < this->standard_form_matrix_.start_[iCol + 1]; iEl++) {
        if (index) index[iEl] = this->standard_form_matrix_.index_[iEl];
        if (value) value[iEl] = this->standard_form_matrix_.value_[iEl];
      }
    }
  }
  if (start) start[num_col] = this->standard_form_matrix_.start_[num_col];
  if (rhs) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      rhs[iRow] = this->standard_form_rhs_[iRow];
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getFixedLp(HighsLp& lp) const {
  this->logHeader();
  if (!this->model_.lp_.isMip()) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Incumbent model is not a MIP, so cannot form fixed LP\n");
    return HighsStatus::kError;
  }
  if (!this->solution_.value_valid) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Incumbent model does not have a valid solution, so cannot "
                 "form fixed LP\n");
    return HighsStatus::kError;
  }
  lp = this->model_.lp_;
  const std::vector<HighsVarType> integrality = this->model_.lp_.integrality_;
  lp.integrality_.clear();
  HighsInt num_non_conts_fractional = 0;
  double max_fractional = 0;
  for (HighsInt iCol = 0; iCol < this->model_.lp_.num_col_; iCol++) {
    double value = this->solution_.col_value[iCol];
    // Fix integer and semi-integer variables at their
    // value. Semi-continuous variables are fixed at zero if they are
    // closer to zero than their lower bound
    if (integrality[iCol] == HighsVarType::kInteger ||
        integrality[iCol] == HighsVarType::kSemiInteger ||
        (integrality[iCol] == HighsVarType::kSemiContinuous &&
         value < lp.col_lower_[iCol] - value)) {
      double fractional = fractionality(value);
      if (fractional > this->options_.mip_feasibility_tolerance) {
        num_non_conts_fractional++;
        max_fractional = std::max(fractional, max_fractional);
      }
      lp.col_lower_[iCol] = value;
      lp.col_upper_[iCol] = value;
    }
  }
  if (num_non_conts_fractional) {
    highsLogUser(
        options_.log_options, HighsLogType::kWarning,
        "Fixed LP has %d variables fixed at max fractional value of %g\n",
        int(num_non_conts_fractional), max_fractional);
    return HighsStatus::kWarning;
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getDualRay(bool& has_dual_ray, double* dual_ray_value) {
  this->logHeader();
  has_dual_ray = false;
  return getDualRayInterface(has_dual_ray, dual_ray_value);
}

HighsStatus Highs::getDualRaySparse(bool& has_dual_ray,
                                    HVector& row_ep_buffer) {
  this->logHeader();
  has_dual_ray = ekk_instance_.dual_ray_record_.index != kNoRayIndex;
  if (has_dual_ray) {
    ekk_instance_.setNlaPointersForLpAndScale(model_.lp_);
    row_ep_buffer.clear();
    row_ep_buffer.count = 1;
    row_ep_buffer.packFlag = true;
    HighsInt iRow = ekk_instance_.dual_ray_record_.index;
    row_ep_buffer.index[0] = iRow;
    row_ep_buffer.array[iRow] = ekk_instance_.dual_ray_record_.sign;

    ekk_instance_.btran(row_ep_buffer, ekk_instance_.info_.row_ep_density);
  }

  return HighsStatus::kOk;
}

HighsStatus Highs::getDualUnboundednessDirection(
    bool& has_dual_unboundedness_direction,
    double* dual_unboundedness_direction_value) {
  this->logHeader();
  if (dual_unboundedness_direction_value) {
    std::vector<double> dual_ray_value(this->model_.lp_.num_row_);
    HighsStatus status =
        getDualRay(has_dual_unboundedness_direction, dual_ray_value.data());
    if (status != HighsStatus::kOk || !has_dual_unboundedness_direction)
      return HighsStatus::kError;
    std::vector<double> dual_unboundedness_direction;
    this->model_.lp_.a_matrix_.productTransposeQuad(
        dual_unboundedness_direction, dual_ray_value);
    for (HighsInt iCol = 0; iCol < this->model_.lp_.num_col_; iCol++)
      dual_unboundedness_direction_value[iCol] =
          dual_unboundedness_direction[iCol];
  } else {
    return getDualRay(has_dual_unboundedness_direction, nullptr);
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getPrimalRay(bool& has_primal_ray,
                                double* primal_ray_value) {
  this->logHeader();
  has_primal_ray = false;
  return getPrimalRayInterface(has_primal_ray, primal_ray_value);
}

HighsStatus Highs::getRanging(HighsRanging& ranging) {
  this->logHeader();
  HighsProfiling profiling;
  const bool already_profiling = this->profiling_;
  if (!already_profiling) this->initializeSingleThreadedProfiling(&profiling);
  HighsStatus return_status = getRangingInterface();
  if (!already_profiling) this->clearProfiling();
  ranging = this->ranging_;
  return return_status;
}

HighsStatus Highs::feasibilityRelaxation(const double global_lower_penalty,
                                         const double global_upper_penalty,
                                         const double global_rhs_penalty,
                                         const double* local_lower_penalty,
                                         const double* local_upper_penalty,
                                         const double* local_rhs_penalty) {
  this->logHeader();
  std::vector<HighsInt> infeasible_row_subset;
  return elasticityFilter(global_lower_penalty, global_upper_penalty,
                          global_rhs_penalty, local_lower_penalty,
                          local_upper_penalty, local_rhs_penalty);
}

HighsStatus Highs::getIllConditioning(HighsIllConditioning& ill_conditioning,
                                      const bool constraint,
                                      const HighsInt method,
                                      const double ill_conditioning_bound) {
  // Cannot be const, as calculations ensure that the incumbent LP is
  // column-wise
  this->logHeader();
  if (!basis_.valid) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot get ill-conditioning without a valid basis\n");
    return HighsStatus::kError;
  }
  return computeIllConditioning(ill_conditioning, constraint, method,
                                ill_conditioning_bound);
}

HighsStatus Highs::getObjectiveBoundScaling(HighsInt& suggested_objective_scale,
                                            HighsInt& suggested_bound_scale) {
  this->logHeader();
  HighsUserScaleData data;
  initialiseUserScaleData(this->options_, data);
  this->model_.lp_.a_matrix_.ensureColwise();
  assessExcessiveObjectiveBoundScaling(this->options_, this->model_, data);
  suggested_objective_scale = data.suggested_user_objective_scale;
  suggested_bound_scale = data.suggested_user_bound_scale;
  return HighsStatus::kOk;
}

HighsStatus Highs::getIis(HighsIis& iis) {
  this->logHeader();
  HighsStatus return_status = this->getIisInterface();
  iis = this->iis_;
  return return_status;
}

HighsStatus Highs::getDualObjectiveValue(
    double& dual_objective_function_value) const {
  this->logHeader();
  bool have_dual_objective_value = computeDualObjectiveValue(
      model_, solution_, dual_objective_function_value);
  return have_dual_objective_value ? HighsStatus::kOk : HighsStatus::kError;
}

HighsStatus Highs::getBasicVariables(HighsInt* basic_variables) {
  this->logHeader();
  if (basic_variables == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasicVariables: basic_variables is NULL\n");
    return HighsStatus::kError;
  }
  HighsProfiling profiling;
  const bool already_profiling = this->profiling_;
  if (!already_profiling) this->initializeSingleThreadedProfiling(&profiling);
  HighsStatus status = getBasicVariablesInterface(basic_variables);
  if (!already_profiling) this->clearProfiling();
  return status;
}

HighsStatus Highs::getBasisInverseRowSparse(const HighsInt row,
                                            HVector& row_ep_buffer) {
  this->logHeader();
  ekk_instance_.setNlaPointersForLpAndScale(model_.lp_);
  row_ep_buffer.clear();
  row_ep_buffer.count = 1;
  row_ep_buffer.index[0] = row;
  row_ep_buffer.array[row] = 1;
  row_ep_buffer.packFlag = true;

  ekk_instance_.btran(row_ep_buffer, ekk_instance_.info_.row_ep_density);

  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisInverseRow(const HighsInt row, double* row_vector,
                                      HighsInt* row_num_nz,
                                      HighsInt* row_indices) {
  this->logHeader();
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisInverseRow: row_vector is NULL\n");
    return HighsStatus::kError;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt num_row = model_.lp_.num_row_;
  if (row < 0 || row >= num_row) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseRow\n",
                 row, num_row - 1);
    return HighsStatus::kError;
  }
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getBasisInverseRow");
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(num_row, 0);
  rhs[row] = 1;
  basisSolveInterface(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisInverseCol(const HighsInt col, double* col_vector,
                                      HighsInt* col_num_nz,
                                      HighsInt* col_indices) {
  this->logHeader();
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisInverseCol: col_vector is NULL\n");
    return HighsStatus::kError;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  HighsInt num_row = model_.lp_.num_row_;
  if (col < 0 || col >= num_row) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getBasisInverseCol\n",
                 col, num_row - 1);
    return HighsStatus::kError;
  }
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getBasisInverseCol");
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(num_row, 0);
  rhs[col] = 1;
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 HighsInt* solution_num_nz,
                                 HighsInt* solution_indices) {
  this->logHeader();
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisSolve: Xrhs is NULL\n");
    return HighsStatus::kError;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisSolve: solution_vector is NULL\n");
    return HighsStatus::kError;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getBasisSolve");
  HighsInt num_row = model_.lp_.num_row_;
  vector<double> rhs;
  rhs.assign(num_row, 0);
  for (HighsInt row = 0; row < num_row; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      false);
  return HighsStatus::kOk;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          HighsInt* solution_num_nz,
                                          HighsInt* solution_indices) {
  this->logHeader();
  if (Xrhs == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisTransposeSolve: Xrhs is NULL\n");
    return HighsStatus::kError;
  }
  if (solution_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getBasisTransposeSolve: solution_vector is NULL\n");
    return HighsStatus::kError;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getBasisTransposeSolve");
  HighsInt num_row = model_.lp_.num_row_;
  vector<double> rhs;
  rhs.assign(num_row, 0);
  for (HighsInt row = 0; row < num_row; row++) rhs[row] = Xrhs[row];
  basisSolveInterface(rhs, solution_vector, solution_num_nz, solution_indices,
                      true);
  return HighsStatus::kOk;
}

HighsStatus Highs::getReducedRow(const HighsInt row, double* row_vector,
                                 HighsInt* row_num_nz, HighsInt* row_indices,
                                 const double* pass_basis_inverse_row_vector) {
  this->logHeader();
  HighsLp& lp = model_.lp_;
  // Ensure that the LP is column-wise
  lp.ensureColwise();
  if (row_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getReducedRow: row_vector is NULL\n");
    return HighsStatus::kError;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not pass_basis_inverse_row_vector
  // NULL - it's the trigger to determine whether it's computed or not
  if (row < 0 || row >= lp.num_row_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Row index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT "] in getReducedRow\n",
                 row, lp.num_row_ - 1);
    return HighsStatus::kError;
  }
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getReducedRow");
  HighsInt num_row = lp.num_row_;
  vector<double> basis_inverse_row;
  double* basis_inverse_row_vector = (double*)pass_basis_inverse_row_vector;
  if (basis_inverse_row_vector == NULL) {
    vector<double> rhs;
    vector<HighsInt> col_indices;
    rhs.assign(num_row, 0);
    rhs[row] = 1;
    basis_inverse_row.resize(num_row, 0);
    // Form B^{-T}e_{row}
    basisSolveInterface(rhs, basis_inverse_row.data(), NULL, NULL, true);
    basis_inverse_row_vector = basis_inverse_row.data();
  }
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    double value = 0;
    for (HighsInt el = lp.a_matrix_.start_[col];
         el < lp.a_matrix_.start_[col + 1]; el++) {
      HighsInt row = lp.a_matrix_.index_[el];
      value += lp.a_matrix_.value_[el] * basis_inverse_row_vector[row];
    }
    row_vector[col] = 0;
    if (fabs(value) > kHighsTiny) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getReducedColumn(const HighsInt col, double* col_vector,
                                    HighsInt* col_num_nz,
                                    HighsInt* col_indices) {
  this->logHeader();
  HighsLp& lp = model_.lp_;
  // Ensure that the LP is column-wise
  lp.ensureColwise();
  if (col_vector == NULL) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getReducedColumn: col_vector is NULL\n");
    return HighsStatus::kError;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
  if (col < 0 || col >= lp.num_col_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Column index %" HIGHSINT_FORMAT
                 " out of range [0, %" HIGHSINT_FORMAT
                 "] in getReducedColumn\n",
                 col, lp.num_col_ - 1);
    return HighsStatus::kError;
  }
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getReducedColumn");
  HighsInt num_row = lp.num_row_;
  vector<double> rhs;
  rhs.assign(num_row, 0);
  for (HighsInt el = lp.a_matrix_.start_[col];
       el < lp.a_matrix_.start_[col + 1]; el++)
    rhs[lp.a_matrix_.index_[el]] = lp.a_matrix_.value_[el];
  basisSolveInterface(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::kOk;
}

HighsStatus Highs::getKappa(double& kappa, const bool exact,
                            const bool report) const {
  this->logHeader();
  if (!ekk_instance_.status_.has_invert)
    return invertRequirementError("getKappa");
  kappa = ekk_instance_.computeBasisCondition(this->model_.lp_, exact, report);
  return HighsStatus::kOk;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = this->model_.lp_;
  // Determine whether a new solution will be defined. If so,
  // the old solution and any basis are cleared
  const bool new_primal_solution =
      lp.num_col_ > 0 &&
      solution.col_value.size() >= static_cast<size_t>(lp.num_col_);
  const bool new_dual_solution =
      lp.num_row_ > 0 &&
      solution.row_dual.size() >= static_cast<size_t>(lp.num_row_);
  const bool new_solution = new_primal_solution || new_dual_solution;

  if (new_solution) {
    invalidateSolverData();
  } else {
    // Solution is rejected, so give a logging message and error
    // return
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "setSolution: User solution is rejected due to mismatch between "
        "size of col_value and row_dual vectors (%d, %d) and number "
        "of columns and rows in the model (%d, %d)\n",
        int(solution.col_value.size()), int(solution.row_dual.size()),
        int(lp.num_col_), int(lp.num_row_));
    return HighsStatus::kError;
  }

  HighsSolution modifiable_solution = solution;
  if (new_primal_solution) {
    // There should be no error return from assessUserColValues, since
    // this corresponds to index illegalities, but use it to warn about
    // infeasibilities in value
    return_status = assessUserColValues(lp, options_, lp.num_col_, nullptr,
                                        solution.col_value.data());
    assert(return_status != HighsStatus::kError);
    if (lp.num_row_ > 0) {
      // Worth computing the row values
      modifiable_solution.row_value.resize(lp.num_row_);
      // Matrix must be column-wise
      lp.a_matrix_.ensureColwise();
      return_status = interpretCallStatus(
          options_.log_options, calculateRowValuesQuad(lp, modifiable_solution),
          return_status, "calculateRowValuesQuad");
      if (return_status == HighsStatus::kError) return return_status;
    }
    solution_.col_value = modifiable_solution.col_value;
    solution_.row_value = modifiable_solution.row_value;
    solution_.value_valid = true;
  }
  if (new_dual_solution) {
    if (lp.num_col_ > 0) {
      // Worth computing the column duals
      modifiable_solution.col_dual.resize(lp.num_col_);
      // Matrix must be column-wise
      lp.a_matrix_.ensureColwise();
      return_status = interpretCallStatus(
          options_.log_options, calculateColDualsQuad(lp, modifiable_solution),
          return_status, "calculateColDuals");
      if (return_status == HighsStatus::kError) return return_status;
    }
    solution_.col_dual = modifiable_solution.col_dual;
    solution_.row_dual = modifiable_solution.row_dual;
    solution_.dual_valid = true;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::getColOrRowName(const HighsLp& lp, const bool is_col,
                                   const HighsInt index,
                                   std::string& name) const {
  this->logHeader();
  HighsInt num_index = is_col ? lp.num_col_ : lp.num_row_;
  if (index < 0 || index >= num_index) {
    highsLogUser(this->options_.log_options, HighsLogType::kError,
                 "Index %d for %s name is outside the range [0, "
                 "num_%s = %d)\n",
                 int(index), is_col ? "column" : "row", is_col ? "col" : "row",
                 int(num_index));
    return HighsStatus::kError;
  }
  const HighsInt num_index_name =
      is_col ? lp.col_names_.size() : lp.row_names_.size();
  if (index >= num_index_name) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Index %d for %s name is outside the range [0, "
                 "num_%s_name = %d)\n",
                 int(index), is_col ? "column" : "row", is_col ? "col" : "row",
                 int(num_index_name));
    return HighsStatus::kError;
  }
  name = is_col ? lp.col_names_[index] : lp.row_names_[index];
  return HighsStatus::kOk;
}

HighsStatus Highs::setSolution(const HighsInt num_entries,
                               const HighsInt* index, const double* value) {
  this->logHeader();
  HighsLp& lp = this->model_.lp_;
  if (lp.num_col_ == 0) return HighsStatus::kOk;
  // Check for illegalities, and warn about duplicates in index, but
  // don't check for feasibility of values, as this will be done when
  // setSolution is called below.
  const bool perform_feasibility_check = false;
  HighsStatus return_status = assessUserColValues(
      lp, options_, num_entries, index, value, perform_feasibility_check);
  if (return_status == HighsStatus::kError) return return_status;

  // Clear the solution, indicate the values not determined by the
  // user, and insert the values determined by the user
  HighsSolution new_solution;
  assignUserColValues(lp.num_col_, new_solution.col_value, num_entries, index,
                      value);
  return interpretCallStatus(options_.log_options, setSolution(new_solution),
                             return_status, "setSolution");
}

HighsStatus Highs::setCallback(HighsCallbackFunctionType user_callback,
                               void* user_callback_data) {
  this->logHeader();
  this->callback_.clear();
  this->callback_.user_callback = std::move(user_callback);
  this->callback_.user_callback_data = user_callback_data;

  options_.log_options.user_callback = this->callback_.user_callback;
  options_.log_options.user_callback_data = this->callback_.user_callback_data;
  options_.log_options.user_callback_active = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::setCallback(HighsCCallbackType c_callback,
                               void* user_callback_data) {
  this->logHeader();
  this->callback_.clear();
  this->callback_.user_callback = [c_callback](
                                      int a, const std::string& b,
                                      const HighsCallbackOutput* cb_out,
                                      HighsCallbackInput* cb_in, void* e) {
    HighsCallbackDataOut cc_out = static_cast<HighsCallbackDataOut>(*cb_out);
    HighsCallbackDataIn cc_in;
    if (cb_in) cc_in = static_cast<HighsCallbackDataIn>(*cb_in);
    c_callback(a, b.c_str(), &cc_out, &cc_in, e);
    if (cb_in) *cb_in = cc_in;  // copy the data in
  };
  this->callback_.user_callback_data = user_callback_data;

  options_.log_options.user_callback = this->callback_.user_callback;
  options_.log_options.user_callback_data = this->callback_.user_callback_data;
  options_.log_options.user_callback_active = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::startCallback(const int callback_type) {
  this->logHeader();
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return HighsStatus::kError;
  if (!this->callback_.user_callback) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot start callback when user_callback not defined\n");
    return HighsStatus::kError;
  }
  assert(int(this->callback_.active.size()) == kNumCallbackType);
  this->callback_.active[callback_type] = true;
  // Possibly modify the logging callback activity
  if (callback_type == kCallbackLogging)
    options_.log_options.user_callback_active = true;
  return HighsStatus::kOk;
}

HighsStatus Highs::startCallback(const HighsCallbackType callback_type) {
  this->logHeader();
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return HighsStatus::kError;
  if (!this->callback_.user_callback) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Cannot start callback when user_callback not defined\n");
    return HighsStatus::kError;
  }
  assert(int(this->callback_.active.size()) == kNumCallbackType);
  this->callback_.active[callback_type] = true;
  // Possibly modify the logging callback activity
  if (callback_type == kCallbackLogging)
    options_.log_options.user_callback_active = true;
  return HighsStatus::kOk;
}

HighsStatus Highs::stopCallback(const int callback_type) {
  this->logHeader();
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return HighsStatus::kError;
  if (!this->callback_.user_callback) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Cannot stop callback when user_callback not defined\n");
    return HighsStatus::kWarning;
  }
  assert(int(this->callback_.active.size()) == kNumCallbackType);
  this->callback_.active[callback_type] = false;
  // Possibly modify the logging callback activity
  if (callback_type == kCallbackLogging)
    options_.log_options.user_callback_active = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::stopCallback(const HighsCallbackType callback_type) {
  this->logHeader();
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return HighsStatus::kError;
  if (!this->callback_.user_callback) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Cannot stop callback when user_callback not defined\n");
    return HighsStatus::kWarning;
  }
  assert(int(this->callback_.active.size()) == kNumCallbackType);
  this->callback_.active[callback_type] = false;
  // Possibly modify the logging callback activity
  if (callback_type == kCallbackLogging)
    options_.log_options.user_callback_active = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::setBasis(const HighsBasis& basis,
                            const std::string& origin) {
  this->logHeader();
  HighsLp& lp = this->model_.lp_;
  if (basis.alien) {
    // An alien basis needs to be checked properly, since it may be
    // singular, or even incomplete.
    if (lp.num_row_ == 0) {
      // Special case where there are no rows, so no singularity
      // issues. All columns with basic status must be set nonbasic
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        basis_.col_status[iCol] =
            basis.col_status[iCol] == HighsBasisStatus::kBasic
                ? HighsBasisStatus::kNonbasic
                : basis.col_status[iCol];
      basis_.alien = false;
    } else {
      // Check whether a new basis can be defined
      if (!isBasisRightSize(lp, basis)) {
        highsLogUser(
            options_.log_options, HighsLogType::kError,
            "setBasis: User basis is rejected due to mismatch between "
            "size of column and row status vectors (%d, %d) and number "
            "of columns and rows in the model (%d, %d)\n",
            int(basis_.col_status.size()), int(basis_.row_status.size()),
            int(lp.num_col_), int(lp.num_row_));
        return HighsStatus::kError;
      }
      HighsBasis modifiable_basis = basis;
      modifiable_basis.was_alien = true;
      // May not have profiling set up - if Highs::setBasis is called
      // before Highs::run, so have to do it here, and has to be
      // single-threaded
      HighsProfiling profiling;
      const bool already_profiling = this->profiling_;
      if (!already_profiling)
        this->initializeSingleThreadedProfiling(&profiling);
      HighsLpSolverObject solver_object(lp, modifiable_basis, solution_, info_,
                                        ekk_instance_, callback_, options_,
                                        timer_);
      solver_object.setProfiling(this->profiling_);
      HighsStatus return_status = formSimplexLpBasisAndFactor(solver_object);
      if (!already_profiling) this->clearProfiling();
      if (return_status != HighsStatus::kOk) return HighsStatus::kError;
      // Update the HiGHS basis
      basis_ = std::move(modifiable_basis);
    }
  } else {
    // Check the user-supplied basis
    if (!isBasisConsistent(lp, basis)) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "setBasis: invalid basis\n");
      return HighsStatus::kError;
    }
    // Update the HiGHS basis
    basis_ = basis;
  }
  basis_.valid = true;
  basis_.useful = true;
  if (origin != "") basis_.debug_origin_name = origin;
  assert(basis_.debug_origin_name != "");
  assert(!basis_.alien);
  if (basis_.was_alien) {
    highsLogDev(
        options_.log_options, HighsLogType::kInfo,
        "Highs::setBasis Was alien = %-5s; Id = %9d; UpdateCount = %4d; Origin "
        "(%s)\n",
        highsBoolToString(basis_.was_alien).c_str(), (int)basis_.debug_id,
        (int)basis_.debug_update_count, basis_.debug_origin_name.c_str());
  }

  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
}

HighsStatus Highs::setBasis() {
  this->logHeader();
  // Invalidate the basis for HiGHS
  //
  // Don't set to logical basis since that causes presolve to be
  // skipped
  this->invalidateBasis();
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::kOk;
}

HighsStatus Highs::putIterate() {
  // Check that there is a simplex iterate to put
  if (!ekk_instance_.status_.has_invert) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "putIterate: no simplex iterate to put\n");
    return HighsStatus::kError;
  }
  ekk_instance_.putIterate();
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::getIterate() {
  // Check that there is a simplex iterate to get
  if (!ekk_instance_.status_.initialised_for_new_lp) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "getIterate: no simplex iterate to get\n");
    return HighsStatus::kError;
  }
  HighsStatus call_status = ekk_instance_.getIterate();
  if (call_status != HighsStatus::kOk) return call_status;
  // Get the corresponding HiGHS basis
  basis_ = ekk_instance_.getHighsBasis(model_.lp_);
  // Clear everything else
  invalidateModelStatusSolutionAndInfo();
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::addCol(const double cost, const double lower_bound,
                          const double upper_bound, const HighsInt num_new_nz,
                          const HighsInt* indices, const double* values) {
  this->logHeader();
  HighsInt starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

HighsStatus Highs::addCols(const HighsInt num_new_col, const double* costs,
                           const double* lower_bounds,
                           const double* upper_bounds,
                           const HighsInt num_new_nz, const HighsInt* starts,
                           const HighsInt* indices, const double* values) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  clearDerivedModelProperties();
  return_status = interpretCallStatus(
      options_.log_options,
      addColsInterface(num_new_col, costs, lower_bounds, upper_bounds,
                       num_new_nz, starts, indices, values),
      return_status, "addCols");
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::addVars(const HighsInt num_new_var, const double* lower,
                           const double* upper) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  // Avoid touching entry [0] of a vector of size 0
  if (num_new_var <= 0) return returnFromHighs(return_status);
  std::vector<double> cost;
  cost.assign(num_new_var, 0);
  return addCols(num_new_var, cost.data(), lower, upper, 0, nullptr, nullptr,
                 nullptr);
}

HighsStatus Highs::addRow(const double lower_bound, const double upper_bound,
                          const HighsInt num_new_nz, const HighsInt* indices,
                          const double* values) {
  this->logHeader();
  HighsInt starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

HighsStatus Highs::addRows(const HighsInt num_new_row,
                           const double* lower_bounds,
                           const double* upper_bounds,
                           const HighsInt num_new_nz, const HighsInt* starts,
                           const HighsInt* indices, const double* values) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  clearDerivedModelProperties();
  return_status = interpretCallStatus(
      options_.log_options,
      addRowsInterface(num_new_row, lower_bounds, upper_bounds, num_new_nz,
                       starts, indices, values),
      return_status, "addRows");
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeObjectiveSense(const ObjSense sense) {
  this->logHeader();
  if ((sense == ObjSense::kMinimize) !=
      (model_.lp_.sense_ == ObjSense::kMinimize)) {
    model_.lp_.sense_ = sense;
    // Nontrivial change
    clearDerivedModelProperties();
    invalidateModelStatusSolutionAndInfo();
  }
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::changeObjectiveOffset(const double offset) {
  this->logHeader();
  // Update the objective value
  info_.objective_function_value += (offset - model_.lp_.offset_);
  model_.lp_.offset_ = offset;
  presolved_model_.lp_.offset_ += offset;
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::changeColIntegrality(const HighsInt col,
                                        const HighsVarType integrality) {
  this->logHeader();
  return changeColsIntegrality(1, &col, &integrality);
}

HighsStatus Highs::changeColsIntegrality(const HighsInt from_col,
                                         const HighsInt to_col,
                                         const HighsVarType* integrality) {
  this->logHeader();
  clearPresolve();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_col, to_col, model_.lp_.num_col_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::changeColsIntegrality "
                 "is out of range [0, %d)\n",
                 int(from_col), int(to_col), int(model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  HighsStatus call_status =
      changeIntegralityInterface(index_collection, integrality);
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "changeIntegrality");
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

static HighsStatus analyseSetCreateError(const HighsLogOptions& log_options,
                                         const std::string& method,
                                         const HighsInt create_error,
                                         const bool ordered,
                                         const HighsInt num_set_entries,
                                         const HighsInt* set,
                                         const HighsInt dimension) {
  if (create_error == kIndexCollectionCreateIllegalSetSize) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Set supplied to Highs::%s has illegal size of %d\n",
                 method.c_str(), int(num_set_entries));
  } else if (create_error == kIndexCollectionCreateIllegalSetOrder) {
    if (ordered) {
      // Creating an index_collection data structure for the set
      // includes a test that the indices increase strictly. If this
      // is not the case then, since an increasing set was created
      // locally, it must contain duplicate entries. return with an
      // error
      highsLogUser(log_options, HighsLogType::kError,
                   "Set supplied to Highs::%s contains duplicate entries\n",
                   method.c_str());
    } else {
      highsLogUser(log_options, HighsLogType::kError,
                   "Set supplied to Highs::%s not ordered\n", method.c_str());
    }
  } else if (create_error < 0) {
    HighsInt illegal_set_index = -1 - create_error;
    HighsInt illegal_set_entry = set[illegal_set_index];
    highsLogUser(
        log_options, HighsLogType::kError,
        "Set supplied to Highs::%s has entry %d of %d out of range [0, %d)\n",
        method.c_str(), int(illegal_set_index), int(illegal_set_entry),
        int(dimension));
  }
  assert(create_error != kIndexCollectionCreateIllegalSetDimension);
  return HighsStatus::kError;
}

HighsStatus Highs::changeColsIntegrality(const HighsInt num_set_entries,
                                         const HighsInt* set,
                                         const HighsVarType* integrality) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  clearPresolve();
  // Ensure that the set and data are in ascending order
  std::vector<HighsVarType> local_integrality{integrality,
                                              integrality + num_set_entries};
  std::vector<HighsInt> local_set{set, set + num_set_entries};
  sortSetData(num_set_entries, local_set, integrality,
              local_integrality.data());
  HighsIndexCollection index_collection;
  const HighsInt create_error = create(index_collection, num_set_entries,
                                       local_set.data(), model_.lp_.num_col_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "changeColsIntegrality",
                                 create_error, true, num_set_entries,
                                 local_set.data(), model_.lp_.num_col_);
  HighsStatus return_status =
      changeIntegralityInterface(index_collection, local_integrality.data());
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColsIntegrality(const HighsInt* mask,
                                         const HighsVarType* integrality) {
  this->logHeader();
  clearPresolve();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_col_);
  assert(!create_error);
  (void)create_error;
  HighsStatus return_status =
      changeIntegralityInterface(index_collection, integrality);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColCost(const HighsInt col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

HighsStatus Highs::changeColsCost(const HighsInt from_col,
                                  const HighsInt to_col, const double* cost) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_col, to_col, model_.lp_.num_col_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::changeColsCost is out "
                 "of range [0, %d)\n",
                 int(from_col), int(to_col), int(model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  HighsStatus return_status = changeCostsInterface(index_collection, cost);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColsCost(const HighsInt num_set_entries,
                                  const HighsInt* set, const double* cost) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  // Check for NULL data in "set" version of changeColsCost since
  // values are sorted with set
  if (doubleUserDataNotNull(options_.log_options, cost, "column costs"))
    return HighsStatus::kError;
  clearDerivedModelProperties();
  // Ensure that the set and data are in ascending order
  std::vector<double> local_cost{cost, cost + num_set_entries};
  std::vector<HighsInt> local_set{set, set + num_set_entries};
  sortSetData(num_set_entries, local_set, cost, NULL, NULL, local_cost.data(),
              NULL, NULL);
  HighsIndexCollection index_collection;
  const HighsInt create_error = create(index_collection, num_set_entries,
                                       local_set.data(), model_.lp_.num_col_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "changeColsCost",
                                 create_error, true, num_set_entries,
                                 local_set.data(), model_.lp_.num_col_);
  HighsStatus return_status =
      changeCostsInterface(index_collection, local_cost.data());
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColsCost(const HighsInt* mask, const double* cost) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_col_);
  assert(!create_error);
  (void)create_error;
  HighsStatus return_status = changeCostsInterface(index_collection, cost);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColBounds(const HighsInt col, const double lower,
                                   const double upper) {
  this->logHeader();
  return changeColsBounds(1, &col, &lower, &upper);
}

HighsStatus Highs::changeColsBounds(const HighsInt from_col,
                                    const HighsInt to_col, const double* lower,
                                    const double* upper) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_col, to_col, model_.lp_.num_col_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::changeColsBounds is out "
                 "of range [0, %d)\n",
                 int(from_col), int(to_col), int(model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  HighsStatus return_status =
      changeColBoundsInterface(index_collection, lower, upper);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColsBounds(const HighsInt num_set_entries,
                                    const HighsInt* set, const double* lower,
                                    const double* upper) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  // Check for NULL data in "set" version of changeColsBounds since
  // values are sorted with set
  bool null_data = false;
  null_data = doubleUserDataNotNull(options_.log_options, lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options_.log_options, upper,
                                    "column upper bounds") ||
              null_data;
  if (null_data) return HighsStatus::kError;
  clearDerivedModelProperties();
  // Ensure that the set and data are in ascending order
  std::vector<double> local_lower{lower, lower + num_set_entries};
  std::vector<double> local_upper{upper, upper + num_set_entries};
  std::vector<HighsInt> local_set{set, set + num_set_entries};
  sortSetData(num_set_entries, local_set, lower, upper, NULL,
              local_lower.data(), local_upper.data(), NULL);
  HighsIndexCollection index_collection;
  const HighsInt create_error = create(index_collection, num_set_entries,
                                       local_set.data(), model_.lp_.num_col_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "changeColsBounds",
                                 create_error, true, num_set_entries,
                                 local_set.data(), model_.lp_.num_col_);
  HighsStatus return_status = changeColBoundsInterface(
      index_collection, local_lower.data(), local_upper.data());
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeColsBounds(const HighsInt* mask, const double* lower,
                                    const double* upper) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_col_);
  assert(!create_error);
  (void)create_error;
  HighsStatus return_status =
      changeColBoundsInterface(index_collection, lower, upper);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeRowBounds(const HighsInt row, const double lower,
                                   const double upper) {
  this->logHeader();
  return changeRowsBounds(1, &row, &lower, &upper);
}

HighsStatus Highs::changeRowsBounds(const HighsInt from_row,
                                    const HighsInt to_row, const double* lower,
                                    const double* upper) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_row, to_row, model_.lp_.num_row_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::changeRowsBounds is out "
                 "of range [0, %d)\n",
                 int(from_row), int(to_row), int(model_.lp_.num_row_));
    return HighsStatus::kError;
  }
  HighsStatus return_status =
      changeRowBoundsInterface(index_collection, lower, upper);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeRowsBounds(const HighsInt num_set_entries,
                                    const HighsInt* set, const double* lower,
                                    const double* upper) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  // Check for NULL data in "set" version of changeRowsBounds since
  // values are sorted with set
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(options_.log_options, lower, "row lower bounds") ||
      null_data;
  null_data =
      doubleUserDataNotNull(options_.log_options, upper, "row upper bounds") ||
      null_data;
  if (null_data) return HighsStatus::kError;
  clearDerivedModelProperties();
  // Ensure that the set and data are in ascending order
  std::vector<double> local_lower{lower, lower + num_set_entries};
  std::vector<double> local_upper{upper, upper + num_set_entries};
  std::vector<HighsInt> local_set{set, set + num_set_entries};
  sortSetData(num_set_entries, local_set, lower, upper, NULL,
              local_lower.data(), local_upper.data(), NULL);
  HighsIndexCollection index_collection;
  const HighsInt create_error = create(index_collection, num_set_entries,
                                       local_set.data(), model_.lp_.num_row_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "changeRowsBounds",
                                 create_error, true, num_set_entries,
                                 local_set.data(), model_.lp_.num_row_);
  HighsStatus return_status = changeRowBoundsInterface(
      index_collection, local_lower.data(), local_upper.data());
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeRowsBounds(const HighsInt* mask, const double* lower,
                                    const double* upper) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_row_);
  assert(!create_error);
  (void)create_error;
  HighsStatus return_status =
      changeRowBoundsInterface(index_collection, lower, upper);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::changeCoeff(const HighsInt row, const HighsInt col,
                               const double value) {
  this->logHeader();
  if (row < 0 || row >= model_.lp_.num_row_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Row %" HIGHSINT_FORMAT
                 " supplied to Highs::changeCoeff is not in the range [0, "
                 "%" HIGHSINT_FORMAT "]\n",
                 row, model_.lp_.num_row_);
    return HighsStatus::kError;
  }
  if (col < 0 || col >= model_.lp_.num_col_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Col %" HIGHSINT_FORMAT
                 " supplied to Highs::changeCoeff is not in the range [0, "
                 "%" HIGHSINT_FORMAT "]\n",
                 col, model_.lp_.num_col_);
    return HighsStatus::kError;
  }
  const double abs_value = std::fabs(value);
  if (abs_value >= options_.large_matrix_value) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "|Value| of %g supplied to Highs::changeCoeff is greater than %g\n",
        abs_value, options_.large_matrix_value);
    return HighsStatus::kError;
  }
  return returnFromHighs(changeCoeffInterface(row, col, value));
}

HighsStatus Highs::getObjectiveSense(ObjSense& sense) const {
  sense = model_.lp_.sense_;
  return HighsStatus::kOk;
}

HighsStatus Highs::getObjectiveOffset(double& offset) const {
  offset = model_.lp_.offset_;
  return HighsStatus::kOk;
}

HighsStatus Highs::getCols(const HighsInt from_col, const HighsInt to_col,
                           HighsInt& num_col, double* costs, double* lower,
                           double* upper, HighsInt& num_nz, HighsInt* start,
                           HighsInt* index, double* value) const {
  this->logHeader();
  if (from_col > to_col) {
    // Empty interval
    num_col = 0;
    num_nz = 0;
    return HighsStatus::kOk;
  }
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_col, to_col, model_.lp_.num_col_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::getCols is out of range "
                 "[0, %d)\n",
                 int(from_col), int(to_col), int(model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  getColsInterface(index_collection, num_col, costs, lower, upper, num_nz,
                   start, index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getCols(const HighsInt num_set_entries, const HighsInt* set,
                           HighsInt& num_col, double* costs, double* lower,
                           double* upper, HighsInt& num_nz, HighsInt* start,
                           HighsInt* index, double* value) const {
  this->logHeader();
  if (num_set_entries == 0) {
    // Empty interval
    num_col = 0;
    num_nz = 0;
    return HighsStatus::kOk;
  }
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, num_set_entries, set, model_.lp_.num_col_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "getCols", create_error,
                                 false, num_set_entries, set,
                                 model_.lp_.num_col_);
  getColsInterface(index_collection, num_col, costs, lower, upper, num_nz,
                   start, index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getCols(const HighsInt* mask, HighsInt& num_col,
                           double* costs, double* lower, double* upper,
                           HighsInt& num_nz, HighsInt* start, HighsInt* index,
                           double* value) const {
  this->logHeader();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_col_);
  assert(!create_error);
  (void)create_error;
  getColsInterface(index_collection, num_col, costs, lower, upper, num_nz,
                   start, index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getColName(const HighsInt col, std::string& name) const {
  this->logHeader();
  return getColOrRowName(this->model_.lp_, true, col, name);
}

HighsStatus Highs::getColByName(const std::string& name, HighsInt& col) {
  this->logHeader();
  HighsLp& lp = model_.lp_;
  if (!lp.col_names_.size()) return HighsStatus::kError;
  if (!lp.col_hash_.size()) lp.col_hash_.form(lp.col_names_);
  std::string from_method = "Highs::getColByName";
  const bool is_column = true;
  return getIndexFromName(options_.log_options, from_method, is_column, name,
                          lp.col_hash_, col, lp.col_names_);
}

HighsStatus Highs::getColIntegrality(const HighsInt col,
                                     HighsVarType& integrality) const {
  this->logHeader();
  if (col < 0 || col >= this->model_.lp_.num_col_) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Index %d for column integrality is outside the range [0, "
                 "num_col = %d)\n",
                 int(col), int(this->model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  integrality = static_cast<size_t>(col) < this->model_.lp_.integrality_.size()
                    ? this->model_.lp_.integrality_[col]
                    : HighsVarType::kContinuous;
  return HighsStatus::kOk;
}

HighsStatus Highs::getRows(const HighsInt from_row, const HighsInt to_row,
                           HighsInt& num_row, double* lower, double* upper,
                           HighsInt& num_nz, HighsInt* start, HighsInt* index,
                           double* value) const {
  this->logHeader();
  if (from_row > to_row) {
    // Empty interval
    num_row = 0;
    num_nz = 0;
    return HighsStatus::kOk;
  }
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_row, to_row, model_.lp_.num_row_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::getRows is out of range "
                 "[0, %d)\n",
                 int(from_row), int(to_row), int(model_.lp_.num_row_));
    return HighsStatus::kError;
  }
  getRowsInterface(index_collection, num_row, lower, upper, num_nz, start,
                   index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getRows(const HighsInt num_set_entries, const HighsInt* set,
                           HighsInt& num_row, double* lower, double* upper,
                           HighsInt& num_nz, HighsInt* start, HighsInt* index,
                           double* value) const {
  this->logHeader();
  if (num_set_entries == 0) {
    num_row = 0;
    num_nz = 0;
    return HighsStatus::kOk;
  }
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, num_set_entries, set, model_.lp_.num_row_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "getRows", create_error,
                                 false, num_set_entries, set,
                                 model_.lp_.num_row_);
  getRowsInterface(index_collection, num_row, lower, upper, num_nz, start,
                   index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getRows(const HighsInt* mask, HighsInt& num_row,
                           double* lower, double* upper, HighsInt& num_nz,
                           HighsInt* start, HighsInt* index,
                           double* value) const {
  this->logHeader();
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, model_.lp_.num_row_);
  assert(!create_error);
  (void)create_error;
  getRowsInterface(index_collection, num_row, lower, upper, num_nz, start,
                   index, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::getRowName(const HighsInt row, std::string& name) const {
  this->logHeader();
  return getColOrRowName(this->model_.lp_, false, row, name);
}

HighsStatus Highs::getRowByName(const std::string& name, HighsInt& row) {
  this->logHeader();
  HighsLp& lp = model_.lp_;
  if (!lp.row_names_.size()) return HighsStatus::kError;
  if (!lp.row_hash_.size()) lp.row_hash_.form(lp.row_names_);
  std::string from_method = "Highs::getRowByName";
  const bool is_column = false;
  return getIndexFromName(options_.log_options, from_method, is_column, name,
                          lp.row_hash_, row, lp.row_names_);
}

HighsStatus Highs::getCoeff(const HighsInt row, const HighsInt col,
                            double& value) const {
  this->logHeader();
  if (row < 0 || row >= model_.lp_.num_row_) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Row %" HIGHSINT_FORMAT
        " supplied to Highs::getCoeff is not in the range [0, %" HIGHSINT_FORMAT
        "]\n",
        row, model_.lp_.num_row_);
    return HighsStatus::kError;
  }
  if (col < 0 || col >= model_.lp_.num_col_) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Col %" HIGHSINT_FORMAT
        " supplied to Highs::getCoeff is not in the range [0, %" HIGHSINT_FORMAT
        "]\n",
        col, model_.lp_.num_col_);
    return HighsStatus::kError;
  }
  getCoeffInterface(row, col, value);
  return HighsStatus::kOk;
}

HighsStatus Highs::deleteCols(const HighsInt from_col, const HighsInt to_col) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_col, to_col, model_.lp_.num_col_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::deleteCols is out of "
                 "range [0, %d)\n",
                 int(from_col), int(to_col), int(model_.lp_.num_col_));
    return HighsStatus::kError;
  }
  deleteColsInterface(index_collection);
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::deleteCols(const HighsInt num_set_entries,
                              const HighsInt* set) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, num_set_entries, set, model_.lp_.num_col_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "deleteCols",
                                 create_error, false, num_set_entries, set,
                                 model_.lp_.num_col_);
  deleteColsInterface(index_collection);
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::deleteCols(HighsInt* mask) {
  this->logHeader();
  clearDerivedModelProperties();
  const HighsInt original_num_col = model_.lp_.num_col_;
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, original_num_col);
  assert(!create_error);
  (void)create_error;
  deleteColsInterface(index_collection);
  for (HighsInt iCol = 0; iCol < original_num_col; iCol++)
    mask[iCol] = index_collection.mask_[iCol];
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::deleteRows(const HighsInt from_row, const HighsInt to_row) {
  this->logHeader();
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, from_row, to_row, model_.lp_.num_row_);
  if (create_error) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Interval [%d, %d] supplied to Highs::deleteRows is out of "
                 "range [0, %d)\n",
                 int(from_row), int(to_row), int(model_.lp_.num_row_));
    return HighsStatus::kError;
  }
  deleteRowsInterface(index_collection);
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::deleteRows(const HighsInt num_set_entries,
                              const HighsInt* set) {
  this->logHeader();
  if (num_set_entries == 0) return HighsStatus::kOk;
  clearDerivedModelProperties();
  HighsIndexCollection index_collection;
  const HighsInt create_error =
      create(index_collection, num_set_entries, set, model_.lp_.num_row_);
  if (create_error)
    return analyseSetCreateError(options_.log_options, "deleteRows",
                                 create_error, false, num_set_entries, set,
                                 model_.lp_.num_row_);
  deleteRowsInterface(index_collection);
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::deleteRows(HighsInt* mask) {
  this->logHeader();
  clearDerivedModelProperties();
  const HighsInt original_num_row = model_.lp_.num_row_;
  HighsIndexCollection index_collection;
  const bool create_error = create(index_collection, mask, original_num_row);
  assert(!create_error);
  (void)create_error;
  deleteRowsInterface(index_collection);
  for (HighsInt iRow = 0; iRow < original_num_row; iRow++)
    mask[iRow] = index_collection.mask_[iRow];
  return returnFromHighs(HighsStatus::kOk);
}

HighsStatus Highs::scaleCol(const HighsInt col, const double scale_value) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  clearDerivedModelProperties();
  HighsStatus call_status = scaleColInterface(col, scale_value);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "scaleCol");
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::scaleRow(const HighsInt row, const double scale_value) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  clearDerivedModelProperties();
  HighsStatus call_status = scaleRowInterface(row, scale_value);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "scaleRow");
  if (return_status == HighsStatus::kError) return HighsStatus::kError;
  return returnFromHighs(return_status);
}

HighsStatus Highs::postsolve(const HighsSolution& solution) {
  this->logHeader();
  HighsBasis basis;
  return this->postsolve(solution, basis);
}

HighsStatus Highs::postsolve(const HighsSolution& solution,
                             const HighsBasis& basis) {
  this->logHeader();
  const bool can_run_postsolve =
      model_presolve_status_ == HighsPresolveStatus::kNotPresolved ||
      model_presolve_status_ == HighsPresolveStatus::kNotReduced ||
      model_presolve_status_ == HighsPresolveStatus::kReduced ||
      model_presolve_status_ == HighsPresolveStatus::kReducedToEmpty ||
      model_presolve_status_ == HighsPresolveStatus::kTimeout;
  if (!can_run_postsolve) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Cannot run postsolve with presolve status: %s\n",
                 presolveStatusToString(model_presolve_status_).c_str());
    return HighsStatus::kWarning;
  }
  HighsProfiling profiling;
  const bool already_profiling = this->profiling_;
  if (!already_profiling) this->initializeSingleThreadedProfiling(&profiling);
  HighsStatus return_status = callRunPostsolve(solution, basis);
  if (!already_profiling) this->clearProfiling();
  return returnFromHighs(return_status);
}

HighsStatus Highs::writeSolution(const std::string& filename,
                                 const HighsInt style) {
  this->logHeader();
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  FILE* file;
  HighsFileType file_type;
  call_status = openWriteFile(filename, "writeSolution", file, file_type);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  call_status = normaliseNames(this->options_.log_options, this->model_.lp_);
  return_status = interpretCallStatus(options_.log_options, call_status,
                                      return_status, "normaliseNames");
  assert(call_status != HighsStatus::kError);

  // Report to user that solution is being written
  if (filename != "")
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Writing the solution to %s\n", filename.c_str());
  writeSolutionFile(file, options_, model_, basis_, solution_, info_,
                    model_status_, style);
  if (style == kSolutionStyleSparse)
    return returnFromWriteSolution(file, return_status);
  if (style == kSolutionStyleRaw) {
    fprintf(file, "\n# Basis\n");
    writeBasisFile(file, options_, model_.lp_, basis_);
  }
  if (options_.ranging == kHighsOnString) {
    if (model_.isMip() || model_.isQp()) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Cannot determine ranging information for MIP or QP\n");
      return_status = HighsStatus::kError;
      return returnFromWriteSolution(file, return_status);
    }
    return_status =
        interpretCallStatus(options_.log_options, this->getRangingInterface(),
                            return_status, "getRangingInterface");
    if (return_status == HighsStatus::kError)
      return returnFromWriteSolution(file, return_status);
    fprintf(file, "\n# Ranging\n");
    writeRangingFile(file, model_.lp_, info_.objective_function_value, basis_,
                     solution_, ranging_, style);
  }
  return returnFromWriteSolution(file, return_status);
}

HighsStatus Highs::readSolution(const std::string& filename,
                                const HighsInt style) {
  this->logHeader();
  return readSolutionFile(filename, options_, model_.lp_, basis_, solution_,
                          style);
}

HighsStatus Highs::assessPrimalSolution(bool& valid, bool& integral,
                                        bool& feasible) const {
  this->logHeader();
  return assessLpPrimalSolution("", options_, model_.lp_, solution_, valid,
                                integral, feasible);
}

std::string Highs::highsStatusToString(const HighsStatus status) const {
  switch (status) {
    case HighsStatus::kOk:
      return "Ok";
    case HighsStatus::kWarning:
      return "Warning";
    case HighsStatus::kError:
      return "Error";
    default:
      assert(1 == 0);
      return "Unrecognised Highs status";
  }
}

std::string Highs::presolveStatusToString(
    const HighsPresolveStatus presolve_status) const {
  switch (presolve_status) {
    case HighsPresolveStatus::kNotPresolved:
      return "Not presolved";
    case HighsPresolveStatus::kNotReduced:
      return "Not reduced";
    case HighsPresolveStatus::kInfeasible:
      return "Infeasible";
    case HighsPresolveStatus::kUnboundedOrInfeasible:
      return "Unbounded or infeasible";
    case HighsPresolveStatus::kReduced:
      return "Reduced";
    case HighsPresolveStatus::kReducedToEmpty:
      return "Reduced to empty";
    case HighsPresolveStatus::kTimeout:
      return "Timeout";
    case HighsPresolveStatus::kOutOfMemory:
      return "Memory allocation error";
    case HighsPresolveStatus::kException:
      return "Exception error";
    case HighsPresolveStatus::kMatrixValueError:
      return "Matrix value error";
    default:
      assert(1 == 0);
      return "Unrecognised presolve status";
  }
}

std::string Highs::modelStatusToString(
    const HighsModelStatus model_status) const {
  return utilModelStatusToString(model_status);
}

std::string Highs::solutionStatusToString(
    const HighsInt solution_status) const {
  return utilSolutionStatusToString(solution_status);
}

std::string Highs::basisStatusToString(
    const HighsBasisStatus basis_status) const {
  return utilBasisStatusToString(basis_status).full_;
}

std::string Highs::basisValidityToString(const HighsInt basis_validity) const {
  return utilBasisValidityToString(basis_validity);
}

std::string Highs::presolveRuleTypeToString(
    const HighsInt presolve_rule) const {
  return utilPresolveRuleTypeToString(presolve_rule);
}

// Private methods
void Highs::deprecationMessage(const std::string& method_name,
                               const std::string& alt_method_name) const {
  if (alt_method_name.compare("None") == 0) {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Method %s is deprecated: no alternative method\n",
                 method_name.c_str());
  } else {
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Method %s is deprecated: alternative method is %s\n",
                 method_name.c_str(), alt_method_name.c_str());
  }
}

HighsPresolveStatus Highs::runPresolve(const bool force_lp_presolve,
                                       const bool force_presolve) {
  presolve_.clear();
  // Exit if presolve is set to off (unless presolve is forced)
  if (options_.presolve == kHighsOffString && !force_presolve)
    return HighsPresolveStatus::kNotPresolved;

  if (model_.isEmpty()) {
    // Empty models shouldn't reach here, but this status would cause
    // no harm if one did
    assert(1 == 0);
    return HighsPresolveStatus::kNotReduced;
  }

  // Ensure that the LP is column-wise
  HighsLp& original_lp = model_.lp_;
  original_lp.ensureColwise();

  if (original_lp.num_col_ == 0 && original_lp.num_row_ == 0)
    return HighsPresolveStatus::kNullError;

  // Ensure that the timer is running
  if (!timer_.running()) timer_.start();
  double start_presolve = timer_.read();

  // Set time limit.
  if (options_.time_limit > 0 && options_.time_limit < kHighsInf) {
    double left = options_.time_limit - start_presolve;
    if (left <= 0) {
      highsLogDev(options_.log_options, HighsLogType::kError,
                  "Time limit reached while reading in matrix\n");
      return HighsPresolveStatus::kTimeout;
    }

    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Time limit set: reading matrix took %.2g, presolve "
                "time left: %.2g\n",
                start_presolve, left);
  }

  // Presolve.
  HighsPresolveStatus presolve_return_status =
      HighsPresolveStatus::kNotPresolved;
  if (model_.isMip() && !force_lp_presolve) {
    // Use presolve for MIP
    //
    // Presolved model is extracted now since it's part of solver,
    // which is lost on return
    HighsMipSolver solver(callback_, options_, original_lp, solution_);
    solver.setProfiling(this->profiling_);
    // Start the MIP solver's timer so that timeout in presolve can be
    // identified
    solver.timer_.start();
    // Only place that HighsMipSolver::runMipPresolve is called
    solver.runMipPresolve(options_.presolve_reduction_limit);
    presolve_return_status = solver.getPresolveStatus();
    // Assign values to data members of presolve_
    presolve_.data_.reduced_lp_ = solver.getPresolvedModel();
    presolve_.data_.postSolveStack = solver.getPostsolveStack();
    presolve_.presolve_status_ = presolve_return_status;
    //    presolve_.data_.presolve_log_ =
  } else {
    // Use presolve for LP
    presolve_.init(original_lp, timer_);
    presolve_.options_ = &options_;
    if (options_.time_limit > 0 && options_.time_limit < kHighsInf) {
      double current = timer_.read();
      double time_init = current - start_presolve;
      double left = presolve_.options_->time_limit - time_init;
      if (left <= 0) {
        highsLogDev(options_.log_options, HighsLogType::kError,
                    "Time limit reached while copying matrix into presolve.\n");
        return HighsPresolveStatus::kTimeout;
      }
      highsLogDev(options_.log_options, HighsLogType::kVerbose,
                  "Time limit set: copying matrix took %.2g, presolve "
                  "time left: %.2g\n",
                  time_init, left);
    }

    presolve_return_status = presolve_.run();
  }

  highsLogDev(options_.log_options, HighsLogType::kVerbose,
              "presolve_.run() returns status: %s\n",
              presolveStatusToString(presolve_return_status).c_str());

  // Update reduction counts.
  assert(presolve_return_status == presolve_.presolve_status_);
  presolve_log_ = presolve_.getPresolveLog();
  switch (presolve_.presolve_status_) {
    case HighsPresolveStatus::kReduced: {
      HighsLp& reduced_lp = presolve_.getReducedProblem();
      presolve_.info_.n_cols_removed =
          original_lp.num_col_ - reduced_lp.num_col_;
      presolve_.info_.n_rows_removed =
          original_lp.num_row_ - reduced_lp.num_row_;
      presolve_.info_.n_nnz_removed =
          (HighsInt)original_lp.numNz() - (HighsInt)reduced_lp.numNz();
      // Clear any scaling information inherited by the reduced LP
      reduced_lp.clearScale();
      assert(reduced_lp.legalDimensions("RunPresolve: reduced_lp",
                                        options_.log_options));
      break;
    }
    case HighsPresolveStatus::kReducedToEmpty: {
      presolve_.info_.n_cols_removed = original_lp.num_col_;
      presolve_.info_.n_rows_removed = original_lp.num_row_;
      presolve_.info_.n_nnz_removed = (HighsInt)original_lp.numNz();
      break;
    }
    default:
      break;
  }
  // Presolve creates integrality vector for an LP, so clear it
  if (!model_.isMip()) presolve_.data_.reduced_lp_.integrality_.clear();

  return presolve_return_status;
}

HighsPostsolveStatus Highs::runPostsolve() {
  // assert(presolve_.has_run_);
  const bool have_primal_solution =
      presolve_.data_.recovered_solution_.value_valid;
  // Need at least a primal solution
  if (!have_primal_solution)
    return HighsPostsolveStatus::kNoPrimalSolutionError;
  const bool have_dual_solution =
      presolve_.data_.recovered_solution_.dual_valid;
  const HighsInt report_3040_col = -578;
  const bool have_basis = presolve_.data_.recovered_basis_.valid;
  presolve_.data_.postSolveStack.undo(
      options_, presolve_.data_.recovered_solution_,
      presolve_.data_.recovered_basis_, 0, report_3040_col);
  if (have_basis &&
      !isBasisConsistent(this->model_.lp_, presolve_.data_.recovered_basis_)) {
    // Recovered basis was assumed to be consistent, but #3323
    // exposes that in exceptional circumstances it may not be,
    // so check here and return solver error if it is
    // inconsistent
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Highs::runPostsolve: Error in basis after postsolve\n");
    return HighsPostsolveStatus::kBasisError;
  }
  // Compute the row activities
  assert(model_.lp_.a_matrix_.isColwise());
  calculateRowValuesQuad(model_.lp_, presolve_.data_.recovered_solution_);

  if (have_dual_solution && model_.lp_.sense_ == ObjSense::kMaximize) {
    // Negate the dual values since the incumbent LP is a maximization
    for (HighsInt iCol = 0; iCol < model_.lp_.num_col_; iCol++)
      presolve_.data_.recovered_solution_.col_dual[iCol] *= -1;
    for (HighsInt iRow = 0; iRow < model_.lp_.num_row_; iRow++)
      presolve_.data_.recovered_solution_.row_dual[iRow] *= -1;
  }

  // Ensure that the postsolve status is used to set
  // presolve_.postsolve_status_, as well as being returned
  HighsPostsolveStatus postsolve_status =
      HighsPostsolveStatus::kSolutionRecovered;
  presolve_.postsolve_status_ = postsolve_status;
  return postsolve_status;
}

void Highs::clearDerivedModelProperties() {
  this->clearPresolve();
  this->clearStandardFormLp();
  this->clearRayRecords();
}

void Highs::clearPresolve() {
  model_presolve_status_ = HighsPresolveStatus::kNotPresolved;
  presolved_model_.clear();
  presolve_.clear();
}

void Highs::clearStandardFormLp() {
  standard_form_valid_ = false;
  standard_form_offset_ = 0;
  standard_form_cost_.clear();
  standard_form_rhs_.clear();
  standard_form_matrix_.clear();
}

void Highs::invalidateSolverData() {
  invalidateSolverDualData();
  invalidateSolution();
  invalidateBasis();
  invalidateEkk();
  clearIis();
}

void Highs::invalidateSolverDualData() {
  invalidateModelStatus();
  invalidateRanging();
  invalidateInfo();
  invalidateRunData();
}

void Highs::invalidateModelStatusSolutionAndInfo() {
  invalidateModelStatusAndInfo();
  invalidateSolution();
}

void Highs::invalidateModelStatusAndInfo() {
  invalidateModelStatus();
  invalidateRanging();
  invalidateInfo();
  invalidateRunData();
  clearIis();
}

void Highs::invalidateModelStatus() {
  model_status_ = HighsModelStatus::kNotset;
}

void Highs::invalidateSolution() {
  info_.primal_solution_status = kSolutionStatusNone;
  info_.dual_solution_status = kSolutionStatusNone;
  info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  this->solution_.invalidate();
}

void Highs::invalidateBasis() {
  info_.basis_validity = kBasisValidityInvalid;
  this->basis_.invalidate();
}

void Highs::invalidateInfo() { info_.invalidate(); }

void Highs::invalidateRunData() { run_data_.invalidate(); }

void Highs::invalidateRanging() { ranging_.invalidate(); }

void Highs::invalidateEkk() { ekk_instance_.invalidate(); }

void Highs::clearIis() { iis_.clear(); }

HighsStatus Highs::completeSolutionFromDiscreteAssignment() {
  // Determine whether the current solution of a MIP is feasible and,
  // if not, try to assign values to continuous variables and discrete
  // variables not at integer values to achieve a feasible
  // solution. Valuable in the case where users make a heuristic
  // (partial) assignment of discrete variables
  assert(model_.isMip() && solution_.value_valid);
  HighsLp& lp = model_.lp_;
  // Determine whether the solution contains undefined values, in
  // order to decide whether to check its feasibility
  const bool contains_undefined_values = solution_.hasUndefined();
  if (!contains_undefined_values) {
    bool valid, integral, feasible;
    // Determine whether this solution is integer feasible
    HighsStatus return_status = assessLpPrimalSolution(
        "", options_, lp, solution_, valid, integral, feasible);
    assert(return_status != HighsStatus::kError);
    if (return_status == HighsStatus::kError) return HighsStatus::kError;
    assert(valid);
    // If the current solution is integer feasible, then it can be
    // used by MIP solver to get a primal bound
    if (feasible) return HighsStatus::kOk;
  }
  // Save the column bounds and integrality in preparation for fixing
  // the discrete variables when user-supplied values are
  // integer
  std::vector<double> save_col_lower = lp.col_lower_;
  std::vector<double> save_col_upper = lp.col_upper_;
  std::vector<HighsVarType> save_integrality = lp.integrality_;
  const bool have_integrality = (lp.integrality_.size() != 0);
  assert(have_integrality);
  // Count the number of fixed and unfixed discrete variables
  HighsInt num_fixed_discrete_variable = 0;
  HighsInt num_unfixed_discrete_variable = 0;
  const bool assess_as_mip = true;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const double primal = solution_.col_value[iCol];
    // Default value is lower bound, unless primal is integer for a
    // discrete variable
    solution_.col_value[iCol] = lp.col_lower_[iCol];
    if (lp.integrality_[iCol] == HighsVarType::kContinuous) continue;
    // Fix discrete variable if its value is defined and integer
    if (primal == kHighsUndefined) {
      num_unfixed_discrete_variable++;
    } else {
      const double lower = lp.col_lower_[iCol];
      const double upper = lp.col_upper_[iCol];
      const HighsVarType type = lp.integrality_[iCol];
      double bound_infeasibility = 0;
      double integer_infeasibility = 0;
      assessColPrimalSolution(assess_as_mip, options_, primal, lower, upper,
                              type, bound_infeasibility, integer_infeasibility);
      if (integer_infeasibility > options_.mip_feasibility_tolerance ||
          bound_infeasibility > options_.mip_feasibility_tolerance) {
        num_unfixed_discrete_variable++;
      } else {
        // Variable is integer and bound feasible, so fix it at this
        // value and remove its integrality
        num_fixed_discrete_variable++;
        lp.col_lower_[iCol] = primal;
        lp.col_upper_[iCol] = primal;
        lp.integrality_[iCol] = HighsVarType::kContinuous;
      }
    }
  }
  assert(!solution_.hasUndefined());
  const HighsInt num_discrete_variable =
      num_unfixed_discrete_variable + num_fixed_discrete_variable;
  const HighsInt num_continuous_variable = lp.num_col_ - num_discrete_variable;
  assert(num_continuous_variable >= 0);
  bool call_run = true;
  const bool few_fixed_discrete_variables =
      10 * num_fixed_discrete_variable < num_discrete_variable;
  if (num_unfixed_discrete_variable == 0) {
    // Solution is integer valued
    if (num_continuous_variable == 0) {
      // There are no continuous variables, so no feasible solution can be
      // deduced
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "No continuous variables, so user-supplied values of "
                   "discrete variables cannot yield "
                   "feasible solution\n");
      call_run = false;
    } else {
      // Solve an LP, so clear all integrality
      lp.integrality_.clear();
      highsLogUser(
          options_.log_options, HighsLogType::kInfo,
          "Attempting to find feasible solution "
          "by solving LP for user-supplied values of discrete variables\n");
    }
  } else {
    // There are unfixed discrete variables
    if (few_fixed_discrete_variables) {
      // Too few discrete variables are fixed so warn, but still
      // attempt to complete a feasible solution
      highsLogUser(
          options_.log_options, HighsLogType::kWarning,
          "User-supplied values fix only %d / %d discrete variables, "
          "so attempt to complete a feasible solution may be expensive\n",
          int(num_fixed_discrete_variable), int(num_discrete_variable));
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Attempting to find feasible solution "
                   "by solving MIP for user-supplied values of %d / %d "
                   "discrete variables\n",
                   int(num_fixed_discrete_variable),
                   int(num_discrete_variable));
    }
  }
  HighsStatus return_status = HighsStatus::kOk;
  // Clear the current solution since either the user solution has
  // been used to fix (a subset of) discrete variables - so a valid
  // solution will be obtained from optimizeModel() if the local model is
  // feasible - or it's not worth using the user solution
  solution_.clear();
  if (call_run) {
    // Solve the model, using mip_max_start_nodes for
    // mip_max_nodes...
    const HighsInt mip_max_nodes = options_.mip_max_nodes;
    options_.mip_max_nodes = options_.mip_max_start_nodes;
    // Solve the model
    basis_.clear();
    if (this->profiling_) {
      // Should not already be a sub-MIP, as handling a user-supplied
      // solution
      assert(!this->profiling_->isSubMip());
    }
    return_status = this->optimizeModel();
    // Reset any profiling data
    if (this->profiling_) this->resetProfiling();
    // ... remembering to recover the original value of mip_max_nodes
    options_.mip_max_nodes = mip_max_nodes;
  }
  // Recover the column bounds and integrality
  lp.col_lower_ = save_col_lower;
  lp.col_upper_ = save_col_upper;
  lp.integrality_ = save_integrality;
  // Handle the error return
  if (return_status == HighsStatus::kError) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Highs::optimizeModel() error trying to find feasible solution\n");
    return HighsStatus::kError;
  }
  return HighsStatus::kOk;
}

// The method below runs calls solveLp for the given LP
HighsStatus Highs::callSolveLp(HighsLp& lp, const std::string& message) {
  HighsStatus return_status = HighsStatus::kOk;

  // Check that the model is column-wise
  assert(lp.a_matrix_.isColwise());

  HighsLpSolverObject solver_object(lp, basis_, solution_, info_, ekk_instance_,
                                    callback_, options_, timer_);
  solver_object.setProfiling(this->profiling_);

  // Solve the LP
  return_status = solveLp(solver_object, message);
  // Extract the model status
  model_status_ = solver_object.model_status_;

  return return_status;
}

// The method below runs calls solveMip for the given MIP
HighsStatus Highs::callSolveMip(HighsLp& lp, const std::string& message) {
  // Need to call invalidateSolverData() - why? - so save any valid
  // previous solution
  HighsSolution previous_solution;
  if (this->solution_.value_valid) {
    previous_solution.value_valid = true;
    previous_solution.col_value = std::move(this->solution_.col_value);
    previous_solution.row_value = std::move(this->solution_.row_value);
  }
  invalidateSolverData();
  // Recover any valid previous solution
  if (previous_solution.value_valid) {
    solution_.value_valid = true;
    solution_.col_value = std::move(previous_solution.col_value);
    solution_.row_value = std::move(previous_solution.row_value);
  }

  HighsStatus return_status = HighsStatus::kOk;

  // Check that the model is column-wise
  assert(lp.a_matrix_.isColwise());

  HighsMipSolverObject solver_object(lp, solution_,
                                     saved_objective_and_solution_, info_,
                                     callback_, options_, timer_);
  solver_object.setProfiling(this->profiling_);

  // Solve the MIP
  return_status = solveMip(solver_object, message);
  // Extract the model status
  model_status_ = solver_object.model_status_;

  return return_status;
}

// The method below runs calls solveQp for the given QP
HighsStatus Highs::callSolveQp(HighsModel& model, const std::string& message) {
  HighsStatus return_status = HighsStatus::kOk;

  // Check that the model is column-wise
  assert(model.lp_.a_matrix_.isColwise());

  HighsQpSolverObject solver_object(model, basis_, solution_, info_, callback_,
                                    options_, timer_);
  solver_object.setProfiling(this->profiling_);

  // Solve the QP
  return_status = solveQp(solver_object, message);
  // Extract the model status
  model_status_ = solver_object.model_status_;
  return return_status;
}

// Only called from Highs::postsolve
HighsStatus Highs::callRunPostsolve(const HighsSolution& solution,
                                    const HighsBasis& basis) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  const HighsLp& presolved_lp = presolve_.getReducedProblem();

  // Must at least have a primal column solution of the right size
  if (HighsInt(solution.col_value.size()) != presolved_lp.num_col_) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Primal solution provided to postsolve is of size %d rather than %d\n",
        int(solution.col_value.size()), int(presolved_lp.num_col_));
    return HighsStatus::kError;
  }
  // Check any basis that is supplied
  const bool basis_supplied =
      basis.col_status.size() > 0 || basis.row_status.size() > 0 || basis.valid;
  if (basis_supplied) {
    if (!isBasisConsistent(presolved_lp, basis)) {
      highsLogUser(
          options_.log_options, HighsLogType::kError,
          "Basis provided to postsolve is incorrect size or inconsistent\n");
      return HighsStatus::kError;
    }
  }
  // Copy in the solution provided
  presolve_.data_.recovered_solution_ = solution;
  // Ignore any row values
  presolve_.data_.recovered_solution_.row_value.assign(presolved_lp.num_row_,
                                                       0);
  presolve_.data_.recovered_solution_.value_valid = true;

  if (this->model_.isMip() && !basis.valid) {
    // Postsolving a MIP without a valid basis - which, if valid,
    // would imply that the relaxation had been solved, a case handled
    // below
    //
    // Ignore any dual values
    presolve_.data_.recovered_solution_.dual_valid = false;
    presolve_.data_.recovered_solution_.col_dual.clear();
    presolve_.data_.recovered_solution_.row_dual.clear();
    // Ignore any basis
    presolve_.data_.recovered_basis_.valid = false;

    HighsPostsolveStatus postsolve_status = runPostsolve();

    if (postsolve_status == HighsPostsolveStatus::kSolutionRecovered) {
      this->solution_ = presolve_.data_.recovered_solution_;
      this->model_status_ = HighsModelStatus::kUnknown;
      invalidateInfo();
      HighsLp& lp = this->model_.lp_;
      this->info_.objective_function_value =
          computeObjectiveValue(lp, this->solution_);
      const bool is_qp = this->model_.isQp();
      assert(!is_qp);
      const bool get_residuals = true;
      getKktFailures(this->options_, is_qp, this->model_.lp_,
                     this->model_.lp_.col_cost_, this->solution_, this->info_,
                     get_residuals);
      double& max_integrality_violation = this->info_.max_integrality_violation;
      max_integrality_violation = 0;
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (lp.integrality_[iCol] == HighsVarType::kInteger) {
          max_integrality_violation =
              std::max(fractionality(this->solution_.col_value[iCol]),
                       max_integrality_violation);
        }
      }
      highsLogUser(
          options_.log_options, HighsLogType::kWarning,
          "Postsolve performed for MIP, but model status cannot be known\n");
    } else {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Postsolve return status: %s\n",
                   utilPostsolveStatusToString(postsolve_status).c_str());
      setHighsModelStatusAndClearSolutionAndBasis(
          HighsModelStatus::kPostsolveError);
    }
  } else {
    // Postsolving an LP, or a MIP after solving the relaxation
    // (identified by passing a valid basis).
    //
    // If there are dual values, make sure that both vectors are the
    // right size
    const bool dual_supplied =
        presolve_.data_.recovered_solution_.col_dual.size() > 0 ||
        presolve_.data_.recovered_solution_.row_dual.size() > 0 ||
        presolve_.data_.recovered_solution_.dual_valid;
    if (dual_supplied) {
      if (!isRowDualSolutionRightSize(presolved_lp,
                                      presolve_.data_.recovered_solution_)) {
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Row dual solution provided to postsolve is of size %d "
                     "rather than %d\n",
                     int(solution.row_dual.size()), int(presolved_lp.num_row_));
        return HighsStatus::kError;
      }
      if (!isColDualSolutionRightSize(presolved_lp,
                                      presolve_.data_.recovered_solution_)) {
        highsLogUser(options_.log_options, HighsLogType::kError,
                     "Column dual solution provided to postsolve is of size %d "
                     "rather than %d\n",
                     int(solution.col_dual.size()), int(presolved_lp.num_col_));
        return HighsStatus::kError;
      }
      presolve_.data_.recovered_solution_.dual_valid = true;
    } else {
      presolve_.data_.recovered_solution_.dual_valid = false;
    }
    // Copy in the basis provided. It's already been checked for
    // consistency, so the basis is valid iff it was supplied
    presolve_.data_.recovered_basis_ = basis;
    presolve_.data_.recovered_basis_.valid = basis_supplied;

    HighsPostsolveStatus postsolve_status = runPostsolve();

    if (postsolve_status == HighsPostsolveStatus::kSolutionRecovered) {
      highsLogDev(options_.log_options, HighsLogType::kVerbose,
                  "Postsolve finished\n");
      // Set solution and its status
      solution_.clear();
      solution_ = presolve_.data_.recovered_solution_;
      assert(solution_.value_valid);
      if (!solution_.dual_valid) {
        solution_.col_dual.assign(model_.lp_.num_col_, 0);
        solution_.row_dual.assign(model_.lp_.num_row_, 0);
      }
      basis_ = presolve_.data_.recovered_basis_;
      // Validity of the solution and basis should be inherited
      //
      // solution_.value_valid = true;
      // solution_.dual_valid = true;
      //
      // Set basis and its status
      //
      // basis_.valid = true;
      // basis_.useful = true;
      // basis_.col_status = presolve_.data_.recovered_basis_.col_status;
      // basis_.row_status = presolve_.data_.recovered_basis_.row_status;
      basis_.debug_origin_name += ": after postsolve";
      if (basis_.valid) {
        // Save the options to allow the best simplex strategy to be
        // used
        HighsOptions save_options = options_;
        options_.simplex_strategy = kSimplexStrategyChoose;
        // Ensure that the parallel solver isn't used
        options_.simplex_min_concurrency = 1;
        options_.simplex_max_concurrency = 1;
        // Use any pivot threshold resulting from solving the presolved LP
        // if (factor_pivot_threshold > 0)
        //    options_.factor_pivot_threshold = factor_pivot_threshold;
        // The basis returned from postsolve is just basic/nonbasic
        // and EKK expects a refined basis, so set it up now
        HighsLp& incumbent_lp = model_.lp_;
        refineBasis(incumbent_lp, solution_, basis_);
        // Scrap the EKK data from solving the presolved LP
        ekk_instance_.invalidate();
        ekk_instance_.lp_name_ = "Postsolve LP";
        timer_.start(timer_.solve_clock);
        call_status = callSolveLp(
            incumbent_lp,
            "Solving the original LP from the solution after postsolve");
        // Determine the timing record
        timer_.stop(timer_.solve_clock);
        return_status = interpretCallStatus(options_.log_options, call_status,
                                            return_status, "callSolveLp");
        // Recover the options
        options_ = save_options;
        const bool is_qp = this->model_.isQp();
        assert(!is_qp);
        const bool get_residuals = true;
        getKktFailures(this->options_, is_qp, this->model_.lp_,
                       this->model_.lp_.col_cost_, this->solution_, this->info_,
                       get_residuals);
        if (return_status == HighsStatus::kError) {
          // Set undo_mods = false, since passing models requiring
          // modification to Highs::presolve is illegal
          const bool undo_mods = false;
          return returnFromOptimizeModel(return_status, undo_mods);
        }
      } else {
        this->basis_.clear();
        // callLpKktCheck will set LPs with model status
        // HighsModelStatus::kUnknown; to HighsModelStatus::kOptimal
        // if relative primal and dual infeasibilities, and the
        // relative primal-dual objective error, are acceptable
        //
        // If there are no dual values, then dual infeasibility
        // measures and the primal-dual objective error are undefined
        // so set to infinity
        //
        // To prevent this from being done with MIPs - if spurious
        // dual values have been set by a user! - set their model
        // status to HighsModelStatus::kNotset
        this->model_status_ = this->model_.lp_.isMip()
                                  ? HighsModelStatus::kNotset
                                  : HighsModelStatus::kUnknown;
        this->callLpKktCheck(this->model_.lp_);
        info_.valid = true;
        highsLogUser(options_.log_options, HighsLogType::kInfo,
                     "\nPure postsolve yields primal %s basis: model "
                     "status is %s\n",
                     solution_.dual_valid ? "and dual solution, but no"
                                          : "but no dual solution or",
                     modelStatusToString(model_status_).c_str());
      }
    } else {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Postsolve return status: %s\n",
                   utilPostsolveStatusToString(postsolve_status).c_str());
      setHighsModelStatusAndClearSolutionAndBasis(
          HighsModelStatus::kPostsolveError);
      // Set undo_mods = false, since passing models requiring
      // modification to Highs::presolve is illegal
      const bool undo_mods = false;
      return returnFromOptimizeModel(HighsStatus::kError, undo_mods);
    }
  }
  call_status = highsStatusFromHighsModelStatus(model_status_);
  return_status =
      interpretCallStatus(options_.log_options, call_status, return_status,
                          "highsStatusFromHighsModelStatus");
  return return_status;
}

// End of public methods
void Highs::logHeader() const {
  if (!options_.output_flag) return;
  //  if (!*options_.log_options.output_flag) return;
  static std::once_flag flag;

  std::call_once(flag, [this]() {
    highsLogHeader(options_.log_options, options_.log_githash);
  });
  return;
}

void Highs::reportModel(const HighsModel& model) {
  reportLp(options_.log_options, model.lp_, HighsLogType::kVerbose);
  if (model.hessian_.dim_) {
    assert(model.hessian_.format_ == HessianFormat::kTriangular);
    const HighsInt dim = model.hessian_.dim_;
    reportHessian(options_.log_options, dim, model.hessian_.start_[dim],
                  model.hessian_.start_.data(), model.hessian_.index_.data(),
                  model.hessian_.value_.data());
  } else if (model.hessian_.isOracle()) {
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Hessian is represented via an oracle\n");
  }
}

// Actions to take if there is a new Highs basis
void Highs::newHighsBasis() {
  // Clear any simplex basis
  ekk_instance_.updateStatus(LpAction::kNewBasis);
}

// Ensure that the HiGHS solution and basis have the same size as the
// model, and that the HiGHS basis is kept up-to-date with any solved
// basis
void Highs::forceHighsSolutionBasisSize() {
  // Ensure that the HiGHS solution and basis vectors are the right size
  //
  // If the column or row, values or duals, are increased in size,
  // ensure that values are assigned to the new entries, and that
  // validity and status values are over-written in HighsSolution,
  // HighsBasis and HighsInfo
  const HighsInt num_col = this->model_.lp_.num_col_;
  const HighsInt num_row = this->model_.lp_.num_row_;
  // Values
  if (solution_.col_value.size() < static_cast<size_t>(num_col) ||
      solution_.row_value.size() < static_cast<size_t>(num_row)) {
    solution_.value_valid = false;
    info_.primal_solution_status = kSolutionStatusNone;
  }
  solution_.col_value.resize(num_col, 0);
  solution_.row_value.resize(num_row, 0);
  // Duals
  if (solution_.col_dual.size() < static_cast<size_t>(num_col) ||
      solution_.row_dual.size() < static_cast<size_t>(num_row)) {
    solution_.dual_valid = false;
    info_.dual_solution_status = kSolutionStatusNone;
  }
  solution_.col_dual.resize(num_col, 0);
  solution_.row_dual.resize(num_row, 0);
  // Basis
  if (basis_.col_status.size() != static_cast<size_t>(num_col) ||
      basis_.row_status.size() != static_cast<size_t>(num_row)) {
    basis_.valid = false;
    basis_.useful = false;
    info_.basis_validity = kBasisValidityInvalid;
  }
  basis_.col_status.resize(num_col, HighsBasisStatus::kNonbasic);
  basis_.row_status.resize(num_row, HighsBasisStatus::kBasic);
}

void Highs::setHighsModelStatusAndClearSolutionAndBasis(
    const HighsModelStatus model_status) {
  model_status_ = model_status;
  invalidateSolution();
  invalidateBasis();
  info_.valid = true;
}

HighsStatus Highs::openWriteFile(const std::string& filename,
                                 const std::string& method_name, FILE*& file,
                                 HighsFileType& file_type) const {
  file_type = HighsFileType::kFull;
  if (filename == "") {
    // Empty file name: use stdout
    file = stdout;
  } else {
    file = fopen(filename.c_str(), "w");
    if (file == 0) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Cannot open writable file \"%s\" in %s\n", filename.c_str(),
                   method_name.c_str());
      return HighsStatus::kError;
    }
    const char* dot = strrchr(filename.c_str(), '.');
    if (dot && dot != filename) {
      if (strcmp(dot + 1, "mps") == 0) {
        file_type = HighsFileType::kMps;
      } else if (strcmp(dot + 1, "lp") == 0) {
        file_type = HighsFileType::kLp;
      } else if (strcmp(dot + 1, "md") == 0) {
        file_type = HighsFileType::kMd;
      }
    }
  }
  return HighsStatus::kOk;
}

// Always called when returning from Highs::writeSolution
HighsStatus Highs::returnFromWriteSolution(FILE* file,
                                           const HighsStatus return_status) {
  if (file != stdout) fclose(file);
  return return_status;
}

// Applies checks before returning from optimizeModel()
HighsStatus Highs::returnFromOptimizeModel(const HighsStatus run_return_status,
                                           const bool undo_mods) {
  assert(!called_return_from_optimize_model);
  HighsStatus return_status = highsStatusFromHighsModelStatus(model_status_);
  if (return_status != run_return_status) {
    highsLogDev(
        options_.log_options, HighsLogType::kError,
        "Highs::returnFromOptimizeModel: run_return_status = %d != %d = "
        "return_status = highsStatusFromHighsModelStatus(model_status_ = %s)\n",
        int(run_return_status), int(return_status),
        modelStatusToString(model_status_).c_str());
  }
  assert(return_status == run_return_status);
  //  return_status = run_return_status;
  switch (model_status_) {
      // First consider the error returns
    case HighsModelStatus::kNotset:
    case HighsModelStatus::kLoadError:
    case HighsModelStatus::kModelError:
    case HighsModelStatus::kPresolveError:
    case HighsModelStatus::kSolveError:
    case HighsModelStatus::kPostsolveError:
    case HighsModelStatus::kMemoryLimit:
      // Don't clear the model status!
      //      invalidateSolverData();
      invalidateInfo();
      invalidateRunData();
      invalidateSolution();
      invalidateBasis();
      assert(return_status == HighsStatus::kError);
      break;

      // Then consider the OK returns
    case HighsModelStatus::kModelEmpty:
      invalidateInfo();
      invalidateRunData();
      invalidateSolution();
      invalidateBasis();
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kOptimal:
      // The following is an aspiration
      //
      // assert(info_.primal_solution_status == kSolutionStatusFeasible);
      //
      // assert(info_.dual_solution_status == kSolutionStatusFeasible);
      assert(model_status_ == HighsModelStatus::kNotset ||
             model_status_ == HighsModelStatus::kOptimal);
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
      // For kInfeasible, will not have a basis, if infeasibility was
      // detected in presolve or by IPX without crossover
      assert(return_status == HighsStatus::kOk);
      break;

    case HighsModelStatus::kUnboundedOrInfeasible:
      if (options_.allow_unbounded_or_infeasible ||
          (useIpm(options_.solver) &&
           options_.run_crossover == kHighsOnString) ||
          (options_.solver == kPdlpString) || model_.isMip()) {
        assert(return_status == HighsStatus::kOk);
      } else {
        // This model status is not permitted unless IPM is run
        // without crossover, if PDLP is used, or if the problem is a
        // MIP. This latter qualification was introduced in the case
        // where MIP presolve identifies kUnboundedOrInfeasible and
        // the problem is not then solved to resolve the model
        // status. However, as exposed by #3383, it's possible that a
        // LP relaxation yields kUnboundedOrInfeasible because the
        // simplex iteration limit is zero.
        highsLogUser(
            options_.log_options, HighsLogType::kError,
            "returnFromHighs: HighsModelStatus::kUnboundedOrInfeasible is not "
            "permitted\n");
        assert(options_.allow_unbounded_or_infeasible);
        return_status = HighsStatus::kError;
      }
      break;

      // Finally consider the warning returns
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kSolutionLimit:
    case HighsModelStatus::kInterrupt:
    case HighsModelStatus::kHighsInterrupt:
    case HighsModelStatus::kUnknown:
      assert(return_status == HighsStatus::kWarning);
      break;
    default:
      // All cases should have been considered so assert on reaching here
      assert(1 == 0);
  }
  // Now to check what's available with each model status
  //
  const bool have_info = info_.valid;
  const bool have_primal_solution = solution_.value_valid;
  const bool have_dual_solution = solution_.dual_valid;
  // Can't have a dual solution without a primal solution
  assert(have_primal_solution || !have_dual_solution);
  //  const bool have_solution = have_primal_solution && have_dual_solution;
  const bool have_basis = basis_.valid;
  switch (model_status_) {
    case HighsModelStatus::kNotset:
    case HighsModelStatus::kLoadError:
    case HighsModelStatus::kModelError:
    case HighsModelStatus::kPresolveError:
    case HighsModelStatus::kSolveError:
    case HighsModelStatus::kPostsolveError:
    case HighsModelStatus::kModelEmpty:
    case HighsModelStatus::kMemoryLimit:
      // No info, primal solution or basis
      assert(have_info == false);
      assert(have_primal_solution == false);
      assert(have_basis == false);
      break;
    case HighsModelStatus::kOptimal:
    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kUnboundedOrInfeasible:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kSolutionLimit:
    case HighsModelStatus::kInterrupt:
    case HighsModelStatus::kHighsInterrupt:
    case HighsModelStatus::kUnknown:
      // Have info and primal solution (unless infeasible). No primal solution
      // in some other case, too!
      assert(have_info == true);
      break;
    default:
      // All cases should have been considered so assert on reaching here
      assert(1 == 0);
  }
  if (have_primal_solution) {
    if (debugPrimalSolutionRightSize(options_, model_.lp_, solution_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (have_dual_solution) {
    if (debugDualSolutionRightSize(options_, model_.lp_, solution_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (have_basis) {
    if (debugBasisRightSize(options_, model_.lp_, basis_) ==
        HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (have_primal_solution) {
    // Debug the Highs solution - needs primal values at least
    if (debugHighsSolution("Return from optimizeModel()", options_, model_,
                           solution_, basis_, model_status_,
                           info_) == HighsDebugStatus::kLogicalError)
      return_status = HighsStatus::kError;
  }
  if (debugInfo(options_, model_.lp_, basis_, solution_, info_,
                model_status_) == HighsDebugStatus::kLogicalError)
    return_status = HighsStatus::kError;

  // Record that returnFromOptimizeModel() has been called
  called_return_from_optimize_model = true;

  if (undo_mods) {
    // Restore any infinite costs
    this->restoreInfCost(return_status);

    // Unapply any modifications that have not yet been unapplied
    this->model_.lp_.unapplyMods();
    //    undo_mods = false;
  }

  // Unless solved as a MIP, report on the solution
  const bool solved_as_mip = model_.isMip() && !options_.solve_relaxation;
  if (!solved_as_mip) reportSolvedLpQpStats();
  return returnFromHighs(return_status);
}

HighsStatus Highs::returnFromHighs(HighsStatus highs_return_status) {
  // Applies checks before returning from HiGHS
  HighsStatus return_status = highs_return_status;

  forceHighsSolutionBasisSize();

  const bool consistent =
      debugHighsBasisConsistent(options_, model_.lp_, basis_) !=
      HighsDebugStatus::kLogicalError;
  if (!consistent) {
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "returnFromHighs: Supposed to be a HiGHS basis, but not consistent\n");
    assert(consistent);
    return_status = HighsStatus::kError;
  }
  // Check that any retained Ekk data - basis and NLA - are OK
  bool retained_ekk_data_ok = ekk_instance_.debugRetainedDataOk(model_.lp_) !=
                              HighsDebugStatus::kLogicalError;
  if (!retained_ekk_data_ok) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "returnFromHighs: Retained Ekk data not OK\n");
    assert(retained_ekk_data_ok);
    return_status = HighsStatus::kError;
  }
  // Check that returnFromOptimizeModel() has been called
  if (!called_return_from_optimize_model) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "Highs::returnFromHighs() called with "
                "called_return_from_optimize_model false\n");
    assert(called_return_from_optimize_model);
  }
  // Stop the HiGHS run clock if it is running
  if (timer_.running()) timer_.stop();
  const bool dimensions_ok =
      model_.lp_.legalDimensions("returnFromHighs", options_.log_options);
  if (!dimensions_ok) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "LP Dimension error in returnFromHighs()\n");
    return_status = HighsStatus::kError;
  }
  assert(dimensions_ok);
  if (ekk_instance_.status_.has_nla) {
    const bool lp_factor_row_compatible =
        ekk_instance_.lpFactorRowCompatible(model_.lp_.num_row_);
    if (!lp_factor_row_compatible) {
      highsLogDev(options_.log_options, HighsLogType::kWarning,
                  "Highs::returnFromHighs(): LP and HFactor have inconsistent "
                  "numbers of rows\n");
      assert(lp_factor_row_compatible);
      // Clear Ekk entirely
      ekk_instance_.clear();
    }
  }
  return return_status;
}

void Highs::reportSolvedLpQpStats() {
  if (!options_.output_flag) return;
  HighsLogOptions& log_options = options_.log_options;
  if (this->model_.lp_.model_name_.length())
    highsLogUser(log_options, HighsLogType::kInfo, "Model name          : %s\n",
                 model_.lp_.model_name_.c_str());
  highsLogUser(log_options, HighsLogType::kInfo, "Model status        : %s\n",
               modelStatusToString(model_status_).c_str());
  if (info_.valid) {
    if (info_.simplex_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Simplex   iterations: %" HIGHSINT_FORMAT "\n",
                   info_.simplex_iteration_count);
    if (info_.ipm_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "IPM       iterations: %" HIGHSINT_FORMAT "\n",
                   info_.ipm_iteration_count);
    if (info_.crossover_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Crossover iterations: %" HIGHSINT_FORMAT "\n",
                   info_.crossover_iteration_count);
    if (info_.pdlp_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "PDLP      iterations: %" HIGHSINT_FORMAT "\n",
                   info_.pdlp_iteration_count);
    if (info_.qp_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "QP ASM    iterations: %" HIGHSINT_FORMAT "\n",
                   info_.qp_iteration_count);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Objective value     : %17.10e\n",
                 info_.objective_function_value);
  }
  if (solution_.dual_valid && info_.primal_dual_objective_error < kHighsInf)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "P-D objective error : %17.10e\n",
                 info_.primal_dual_objective_error);
  if (!options_.timeless_log) {
    double run_time = timer_.read();
    highsLogUser(log_options, HighsLogType::kInfo,
                 "HiGHS run time      : %13.2f\n", run_time);
  }
}

HighsStatus Highs::crossover(const HighsSolution& user_solution) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLogOptions& log_options = options_.log_options;
  HighsLp& lp = model_.lp_;
  if (lp.isMip()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Cannot apply crossover to solve MIP\n");
    return_status = HighsStatus::kError;
  } else if (model_.isQp()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Cannot apply crossover to solve QP\n");
    return_status = HighsStatus::kError;
  } else {
    clearSolver();
    solution_ = user_solution;
    // Use IPX crossover to try to form a basic solution
    return_status = callCrossover(options_, model_.lp_, basis_, solution_,
                                  model_status_, info_, callback_);
    if (return_status == HighsStatus::kError) return return_status;
    // Get the objective and any KKT failures
    info_.objective_function_value =
        model_.lp_.objectiveValue(solution_.col_value);
    getLpKktFailures(options_, model_.lp_, solution_, basis_, info_);
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::openLogFile(const std::string& log_file) {
  highsOpenLogFile(options_.log_options, options_.records, log_file);
  return HighsStatus::kOk;
}

HighsStatus Highs::closeLogFile() {
  // Work with a copy of the pointer for brevity
  FILE* log_stream = this->options_.log_options.log_stream;
  if (log_stream != nullptr) {
    assert(log_stream != stdout);
    fflush(log_stream);
    fclose(log_stream);
    // Set the true log_stream to nullptr to give a test whether it
    // has been closed (and avoid trying to close it again which
    // causes an error)
    this->options_.log_options.log_stream = nullptr;
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::initializeMultiThreading() {
  highs::parallel::initialize_scheduler(this->options_.threads);
  this->max_threads_ = highs::parallel::num_threads();
  HighsLogOptions& log_options = this->options_.log_options;
  unsigned int available_cores = highs::parallel::available_core_count();
  highsLogDev(log_options, HighsLogType::kInfo,
              "Scheduler initialized: %d threads (physical cores: %d)\n",
              static_cast<int>(this->max_threads_),
              static_cast<int>(available_cores));
  if (this->max_threads_ > static_cast<int>(available_cores))
    highsLogDev(
        log_options, HighsLogType::kWarning,
        "Thread count (%d) exceeds physical cores (%d); consider reducing "
        "threads or calling Highs::resetThreadScheduler()\n",
        static_cast<int>(this->max_threads_),
        static_cast<int>(available_cores));
  if (this->options_.threads != 0 &&
      this->max_threads_ != this->options_.threads) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Option 'threads' is set to %d but the scheduler has already been "
        "initialized to use %d threads. The previous scheduler instance can "
        "be destroyed by calling Highs::resetThreadScheduler().\n",
        int(options_.threads), int(this->max_threads_));
    return HighsStatus::kError;
  }
  if (this->max_threads_ <= 0) {
    highsLogDev(log_options, HighsLogType::kError, "max_threads() returns %d\n",
                int(this->max_threads_));
    assert(this->max_threads_ > 0);
    return HighsStatus::kError;
  }
  highsLogDev(log_options, HighsLogType::kDetailed,
              "Running with %d thread(s)\n", int(max_threads_));

  return HighsStatus::kOk;
}

void Highs::resetThreadScheduler(bool blocking) {
  HighsTaskExecutor::shutdown(blocking);
}

void Highs::resetGlobalScheduler(bool blocking) {
  // Deprecated in favour of resetThreadScheduler().
  // Cannot call deprecationMessage here.
  resetThreadScheduler(blocking);
}

HighsStatus Highs::initializeExtras(const std::string& directory,
                                    bool search_if_empty) {
  return HighsExtras::initialize(directory, search_if_empty);
}

void Highs::initializeSingleThreadedProfiling(HighsProfiling* profiling) {
  if (!profiling) return;
  profiling->multi_threaded = false;
  this->initializeProfiling(profiling);
}

void Highs::initializeProfiling(HighsProfiling* profiling) {
  if (!profiling) return;
  assert(!this->profiling_);
  if (this->profiling_) return;
  // Only initialize profiling if this->profiling_ is nullptr
  const bool sub_solver = this->options_.log_dev_level > 0;  // true; //
  // Cannot perform MIP profiling without sub-solver profiling
  const bool mip = sub_solver && kHighsAnalysisLevelMipTime &
                                     this->options_.highs_analysis_level;
  profiling->initialize(this->timer_, sub_solver, mip);
  profiling->model_name_ = this->model_.lp_.model_name_;
  this->setProfiling(profiling);
}

void Highs::resetProfiling() {
  if (!this->profiling_) return;
  HighsProfiling* profiling = this->profiling_;
  this->clearProfiling();
  this->initializeProfiling(profiling);
}

void Highs::clearProfiling() {
  if (!this->profiling_) return;
  this->profiling_->clear();
  this->profiling_ = nullptr;
}

void Highs::setProfiling(HighsProfiling* profiling) {
  assert(profiling);
  this->profiling_ = profiling;
}
