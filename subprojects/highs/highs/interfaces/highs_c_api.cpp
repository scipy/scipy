/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "highs_c_api.h"

#include "Highs.h"

HighsInt Highs_lpCall(const HighsInt num_col, const HighsInt num_row,
                      const HighsInt num_nz, const HighsInt a_format,
                      const HighsInt sense, const double offset,
                      const double* col_cost, const double* col_lower,
                      const double* col_upper, const double* row_lower,
                      const double* row_upper, const HighsInt* a_start,
                      const HighsInt* a_index, const double* a_value,
                      double* col_value, double* col_dual, double* row_value,
                      double* row_dual, HighsInt* col_basis_status,
                      HighsInt* row_basis_status, HighsInt* model_status) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  *model_status = kHighsModelStatusNotset;
  HighsStatus status = highs.passModel(
      num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower,
      col_upper, row_lower, row_upper, a_start, a_index, a_value);
  if (status == HighsStatus::kError) return (HighsInt)status;

  status = highs.run();

  if (status == HighsStatus::kOk) {
    const HighsSolution& solution = highs.getSolution();
    const HighsBasis& basis = highs.getBasis();
    *model_status = (HighsInt)highs.getModelStatus();
    const HighsInfo& info = highs.getInfo();

    const bool copy_col_value =
        col_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_col_dual =
        col_dual != nullptr &&
        info.dual_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_col_basis = col_basis_status != nullptr && basis.valid;
    for (HighsInt i = 0; i < num_col; i++) {
      if (copy_col_value) col_value[i] = solution.col_value[i];
      if (copy_col_dual) col_dual[i] = solution.col_dual[i];
      if (copy_col_basis) col_basis_status[i] = (HighsInt)basis.col_status[i];
    }

    const bool copy_row_value =
        row_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_row_dual =
        row_dual != nullptr &&
        info.dual_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_row_basis = row_basis_status != nullptr && basis.valid;
    for (HighsInt i = 0; i < num_row; i++) {
      if (copy_row_value) row_value[i] = solution.row_value[i];
      if (copy_row_dual) row_dual[i] = solution.row_dual[i];
      if (copy_row_basis) row_basis_status[i] = (HighsInt)basis.row_status[i];
    }
  }
  return (HighsInt)status;
}

HighsInt Highs_mipCall(const HighsInt num_col, const HighsInt num_row,
                       const HighsInt num_nz, const HighsInt a_format,
                       const HighsInt sense, const double offset,
                       const double* col_cost, const double* col_lower,
                       const double* col_upper, const double* row_lower,
                       const double* row_upper, const HighsInt* a_start,
                       const HighsInt* a_index, const double* a_value,
                       const HighsInt* integrality, double* col_value,
                       double* row_value, HighsInt* model_status) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  *model_status = kHighsModelStatusNotset;
  HighsStatus status = highs.passModel(
      num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower,
      col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality);
  if (status == HighsStatus::kError) return (HighsInt)status;

  status = highs.run();

  if (status == HighsStatus::kOk) {
    const HighsSolution& solution = highs.getSolution();
    *model_status = (HighsInt)highs.getModelStatus();
    const HighsInfo& info = highs.getInfo();
    const bool copy_col_value =
        col_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;

    if (copy_col_value) {
      for (HighsInt i = 0; i < num_col; i++)
        col_value[i] = solution.col_value[i];
    }
    const bool copy_row_value =
        row_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    if (copy_row_value) {
      for (HighsInt i = 0; i < num_row; i++)
        row_value[i] = solution.row_value[i];
    }
  }
  return (HighsInt)status;
}

HighsInt Highs_qpCall(
    const HighsInt num_col, const HighsInt num_row, const HighsInt num_nz,
    const HighsInt q_num_nz, const HighsInt a_format, const HighsInt q_format,
    const HighsInt sense, const double offset, const double* col_cost,
    const double* col_lower, const double* col_upper, const double* row_lower,
    const double* row_upper, const HighsInt* a_start, const HighsInt* a_index,
    const double* a_value, const HighsInt* q_start, const HighsInt* q_index,
    const double* q_value, double* col_value, double* col_dual,
    double* row_value, double* row_dual, HighsInt* col_basis_status,
    HighsInt* row_basis_status, HighsInt* model_status) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  *model_status = kHighsModelStatusNotset;
  HighsStatus status = highs.passModel(
      num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset,
      col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index,
      a_value, q_start, q_index, q_value);
  if (status == HighsStatus::kError) return (HighsInt)status;

  status = highs.run();

  if (status == HighsStatus::kOk) {
    const HighsSolution& solution = highs.getSolution();
    const HighsBasis& basis = highs.getBasis();
    *model_status = (HighsInt)highs.getModelStatus();
    const HighsInfo& info = highs.getInfo();

    const bool copy_col_value =
        col_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_col_dual =
        col_dual != nullptr &&
        info.dual_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_col_basis = col_basis_status != nullptr && basis.valid;
    for (HighsInt i = 0; i < num_col; i++) {
      if (copy_col_value) col_value[i] = solution.col_value[i];
      if (copy_col_dual) col_dual[i] = solution.col_dual[i];
      if (copy_col_basis) col_basis_status[i] = (HighsInt)basis.col_status[i];
    }

    const bool copy_row_value =
        row_value != nullptr &&
        info.primal_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_row_dual =
        row_dual != nullptr &&
        info.dual_solution_status != SolutionStatus::kSolutionStatusNone;
    const bool copy_row_basis = row_basis_status != nullptr && basis.valid;
    for (HighsInt i = 0; i < num_row; i++) {
      if (copy_row_value) row_value[i] = solution.row_value[i];
      if (copy_row_dual) row_dual[i] = solution.row_dual[i];
      if (copy_row_basis) row_basis_status[i] = (HighsInt)basis.row_status[i];
    }
  }
  return (HighsInt)status;
}

void* Highs_create(void) { return new Highs(); }

void Highs_destroy(void* highs) { delete (Highs*)highs; }

const char* Highs_version(void) { return highsVersion(); }
HighsInt Highs_versionMajor(void) { return highsVersionMajor(); }
HighsInt Highs_versionMinor(void) { return highsVersionMinor(); }
HighsInt Highs_versionPatch(void) { return highsVersionPatch(); }
const char* Highs_githash(void) { return highsGithash(); }

HighsInt Highs_presolve(void* highs) {
  return (HighsInt)((Highs*)highs)->presolve();
}

HighsInt Highs_run(void* highs) { return (HighsInt)((Highs*)highs)->run(); }

HighsInt Highs_postsolve(void* highs, const double* col_value,
                         const double* col_dual, const double* row_dual) {
  const HighsLp& presolved_lp = ((Highs*)highs)->getPresolvedLp();
  HighsInt num_col = presolved_lp.num_col_;
  HighsInt num_row = presolved_lp.num_row_;
  // Create a HighsSolution from what's been passed
  HighsSolution solution;
  if (col_value) {
    solution.value_valid = true;
    solution.col_value.resize(num_col);
    // No need for primal row values, but resize the vector for later
    // use
    solution.row_value.resize(num_row);
  }
  if (col_dual || row_dual) {
    // If column or row duals are passed, assume that they are
    // valid. If either is a null pointer, then the corresponding
    // vector will have no data, and the size check will fail
    solution.dual_valid = true;
    if (col_dual) solution.col_dual.resize(num_col);
    if (row_dual) solution.row_dual.resize(num_row);
  }
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (col_value) solution.col_value[iCol] = col_value[iCol];
    if (col_dual) solution.col_dual[iCol] = col_dual[iCol];
  }
  if (row_dual) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      solution.row_dual[iRow] = row_dual[iRow];
  }
  return (HighsInt)((Highs*)highs)->postsolve(solution);
}

HighsInt Highs_readModel(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->readModel(std::string(filename));
}

HighsInt Highs_writeModel(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeModel(std::string(filename));
}

HighsInt Highs_writePresolvedModel(void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writePresolvedModel(std::string(filename));
}

HighsInt Highs_writeSolution(const void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)
      ->writeSolution(std::string(filename), kSolutionStyleRaw);
}

HighsInt Highs_writeSolutionPretty(const void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)
      ->writeSolution(std::string(filename), kSolutionStylePretty);
}

HighsInt Highs_passLp(void* highs, const HighsInt num_col,
                      const HighsInt num_row, const HighsInt num_nz,
                      const HighsInt a_format, const HighsInt sense,
                      const double offset, const double* col_cost,
                      const double* col_lower, const double* col_upper,
                      const double* row_lower, const double* row_upper,
                      const HighsInt* a_start, const HighsInt* a_index,
                      const double* a_value) {
  return (HighsInt)((Highs*)highs)
      ->passModel(num_col, num_row, num_nz, a_format, sense, offset, col_cost,
                  col_lower, col_upper, row_lower, row_upper, a_start, a_index,
                  a_value);
}

HighsInt Highs_passMip(void* highs, const HighsInt num_col,
                       const HighsInt num_row, const HighsInt num_nz,
                       const HighsInt a_format, const HighsInt sense,
                       const double offset, const double* col_cost,
                       const double* col_lower, const double* col_upper,
                       const double* row_lower, const double* row_upper,
                       const HighsInt* a_start, const HighsInt* a_index,
                       const double* a_value, const HighsInt* integrality) {
  return (HighsInt)((Highs*)highs)
      ->passModel(num_col, num_row, num_nz, a_format, sense, offset, col_cost,
                  col_lower, col_upper, row_lower, row_upper, a_start, a_index,
                  a_value, integrality);
}

HighsInt Highs_passModel(void* highs, const HighsInt num_col,
                         const HighsInt num_row, const HighsInt num_nz,
                         const HighsInt q_num_nz, const HighsInt a_format,
                         const HighsInt q_format, const HighsInt sense,
                         const double offset, const double* col_cost,
                         const double* col_lower, const double* col_upper,
                         const double* row_lower, const double* row_upper,
                         const HighsInt* a_start, const HighsInt* a_index,
                         const double* a_value, const HighsInt* q_start,
                         const HighsInt* q_index, const double* q_value,
                         const HighsInt* integrality) {
  return (HighsInt)((Highs*)highs)
      ->passModel(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense,
                  offset, col_cost, col_lower, col_upper, row_lower, row_upper,
                  a_start, a_index, a_value, q_start, q_index, q_value,
                  integrality);
}

HighsInt Highs_passHessian(void* highs, const HighsInt dim,
                           const HighsInt num_nz, const HighsInt format,
                           const HighsInt* start, const HighsInt* index,
                           const double* value) {
  return (HighsInt)((Highs*)highs)
      ->passHessian(dim, num_nz, format, start, index, value);
}

HighsInt Highs_passHessianOracle(void* highs, const HighsInt dim,
                                 HighsCHessianFunctionType oracleCall,
                                 void* oracle_data) {
  auto status =
      static_cast<Highs*>(highs)->passHessian(dim, oracleCall, oracle_data);
  return static_cast<HighsInt>(status);
}

HighsInt Highs_checkHessianOracle(void* highs) {
  auto status = static_cast<Highs*>(highs)->checkHessianOracle();
  return static_cast<HighsInt>(status);
}

HighsInt Highs_passLinearObjectives(const void* highs,
                                    const HighsInt num_linear_objective,
                                    const double* weight, const double* offset,
                                    const double* coefficients,
                                    const double* abs_tolerance,
                                    const double* rel_tolerance,
                                    const HighsInt* priority) {
  HighsInt status = Highs_clearLinearObjectives(highs);
  if (status != kHighsStatusOk) return status;
  HighsLinearObjective linear_objective;
  for (HighsInt iObj = 0; iObj < num_linear_objective; iObj++) {
    HighsInt num_col = Highs_getNumCol(highs);
    linear_objective.weight = weight[iObj];
    linear_objective.offset = offset[iObj];
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      linear_objective.coefficients.push_back(
          coefficients[iObj * num_col + iCol]);
    linear_objective.abs_tolerance = abs_tolerance[iObj];
    linear_objective.rel_tolerance = rel_tolerance[iObj];
    linear_objective.priority = priority[iObj];
    linear_objective.weight = weight[iObj];
    status =
        HighsInt(((Highs*)highs)->addLinearObjective(linear_objective, iObj));
    if (status != kHighsStatusOk) return status;
    linear_objective.coefficients.clear();
  }
  return kHighsStatusOk;
}

HighsInt Highs_addLinearObjective(const void* highs, const double weight,
                                  const double offset,
                                  const double* coefficients,
                                  const double abs_tolerance,
                                  const double rel_tolerance,
                                  const HighsInt priority) {
  HighsLinearObjective linear_objective;
  HighsInt num_col = Highs_getNumCol(highs);
  linear_objective.weight = weight;
  linear_objective.offset = offset;
  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    linear_objective.coefficients.push_back(coefficients[iCol]);
  linear_objective.abs_tolerance = abs_tolerance;
  linear_objective.rel_tolerance = rel_tolerance;
  linear_objective.priority = priority;
  linear_objective.weight = weight;
  return HighsInt(((Highs*)highs)->addLinearObjective(linear_objective));
}

HighsInt Highs_clearLinearObjectives(const void* highs) {
  return HighsInt(((Highs*)highs)->clearLinearObjectives());
}

HighsInt Highs_passRowName(const void* highs, const HighsInt row,
                           const char* name) {
  return (HighsInt)((Highs*)highs)->passRowName(row, std::string(name));
}

HighsInt Highs_passColName(const void* highs, const HighsInt col,
                           const char* name) {
  return (HighsInt)((Highs*)highs)->passColName(col, std::string(name));
}

HighsInt Highs_passModelName(const void* highs, const char* name) {
  return (HighsInt)((Highs*)highs)->passModelName(std::string(name));
}

HighsInt Highs_readOptions(const void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->readOptions(filename);
}

HighsInt Highs_clear(void* highs) { return (HighsInt)((Highs*)highs)->clear(); }

HighsInt Highs_clearModel(void* highs) {
  return (HighsInt)((Highs*)highs)->clearModel();
}

HighsInt Highs_clearSolver(void* highs) {
  return (HighsInt)((Highs*)highs)->clearSolver();
}

HighsInt Highs_releaseMemory(void* highs) {
  return (HighsInt)((Highs*)highs)->releaseMemory();
}

HighsInt Highs_setBoolOptionValue(void* highs, const char* option,
                                  const HighsInt value) {
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), (bool)value);
}

HighsInt Highs_setIntOptionValue(void* highs, const char* option,
                                 const HighsInt value) {
  return (HighsInt)((Highs*)highs)->setOptionValue(std::string(option), value);
}

HighsInt Highs_setDoubleOptionValue(void* highs, const char* option,
                                    const double value) {
  return (HighsInt)((Highs*)highs)->setOptionValue(std::string(option), value);
}

HighsInt Highs_setStringOptionValue(void* highs, const char* option,
                                    const char* value) {
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_getNumOptions(const void* highs) {
  return ((Highs*)highs)->getNumOptions();
}

HighsInt Highs_getOptionName(const void* highs, const HighsInt index,
                             char** name) {
  std::string name_v;
  HighsInt retcode = (HighsInt)((Highs*)highs)->getOptionName(index, &name_v);
  // Guess we have to add one (for Windows, lol!) because char* is
  // null-terminated
  const HighsInt malloc_size = sizeof(char) * (name_v.length() + 1);
  *name = (char*)malloc(malloc_size);
  strcpy(*name, name_v.c_str());
  return retcode;
}

HighsInt Highs_getBoolOptionValue(const void* highs, const char* option,
                                  HighsInt* value) {
  return Highs_getBoolOptionValues(highs, option, value, nullptr);
}

HighsInt Highs_getIntOptionValue(const void* highs, const char* option,
                                 HighsInt* value) {
  return Highs_getIntOptionValues(highs, option, value, nullptr, nullptr,
                                  nullptr);
}

HighsInt Highs_getDoubleOptionValue(const void* highs, const char* option,
                                    double* value) {
  return Highs_getDoubleOptionValues(highs, option, value, nullptr, nullptr,
                                     nullptr);
}

HighsInt Highs_getStringOptionValue(const void* highs, const char* option,
                                    char* value) {
  return Highs_getStringOptionValues(highs, option, value, nullptr);
}

HighsInt Highs_getOptionType(const void* highs, const char* option,
                             HighsInt* type) {
  HighsOptionType t;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getOptionType(std::string(option), t);
  *type = (HighsInt)t;
  return retcode;
}

HighsInt Highs_getBoolOptionValues(const void* highs, const char* option,
                                   HighsInt* current_value,
                                   HighsInt* default_value) {
  bool current_v;
  bool default_v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)
          ->getBoolOptionValues(std::string(option), &current_v, &default_v);
  if (current_value) *current_value = current_v;
  if (default_value) *default_value = default_v;
  return retcode;
}

HighsInt Highs_getIntOptionValues(const void* highs, const char* option,
                                  HighsInt* current_value, HighsInt* min_value,
                                  HighsInt* max_value,
                                  HighsInt* default_value) {
  return (HighsInt)((Highs*)highs)
      ->getIntOptionValues(std::string(option), current_value, min_value,
                           max_value, default_value);
}

HighsInt Highs_getDoubleOptionValues(const void* highs, const char* option,
                                     double* current_value, double* min_value,
                                     double* max_value, double* default_value) {
  return (HighsInt)((Highs*)highs)
      ->getDoubleOptionValues(std::string(option), current_value, min_value,
                              max_value, default_value);
}

HighsInt Highs_getStringOptionValues(const void* highs, const char* option,
                                     char* current_value, char* default_value) {
  std::string current_v;
  std::string default_v;
  // Inherited from Highs_getStringOptionValue: cannot see why this
  // was ever useful. Must assume that current_value is of length at
  // least 7, which isn't necessarily true
  //
  //  if (current_value) memset(current_value, 0, 7);
  //  if (default_value) memset(default_value, 0, 7);
  HighsInt retcode =
      (HighsInt)((Highs*)highs)
          ->getStringOptionValues(std::string(option), &current_v, &default_v);
  // current_value and default_value are nullptr by default
  if (current_value) strcpy(current_value, current_v.c_str());
  if (default_value) strcpy(default_value, default_v.c_str());
  return retcode;
}

HighsInt Highs_resetOptions(void* highs) {
  return (HighsInt)((Highs*)highs)->resetOptions();
}

HighsInt Highs_writeOptions(const void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeOptions(filename);
}

HighsInt Highs_writeOptionsDeviations(const void* highs, const char* filename) {
  return (HighsInt)((Highs*)highs)->writeOptions(filename, true);
}

HighsInt Highs_getIntInfoValue(const void* highs, const char* info,
                               HighsInt* value) {
  return (HighsInt)((Highs*)highs)->getInfoValue(info, *value);
}

HighsInt Highs_getDoubleInfoValue(const void* highs, const char* info,
                                  double* value) {
  return (HighsInt)((Highs*)highs)->getInfoValue(info, *value);
}

HighsInt Highs_getInt64InfoValue(const void* highs, const char* info,
                                 int64_t* value) {
  return (HighsInt)((Highs*)highs)->getInfoValue(info, *value);
}

HighsInt Highs_getInfoType(const void* highs, const char* info,
                           HighsInt* type) {
  HighsInfoType t;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getInfoType(std::string(info), t);
  *type = (HighsInt)t;
  return retcode;
}

HighsInt Highs_getSolution(const void* highs, double* col_value,
                           double* col_dual, double* row_value,
                           double* row_dual) {
  if (((Highs*)highs)->getInfo().primal_solution_status == kSolutionStatusNone)
    return kHighsStatusError;

  const HighsSolution& solution = ((Highs*)highs)->getSolution();

  if (col_value != nullptr) {
    for (size_t i = 0; i < solution.col_value.size(); i++) {
      col_value[i] = solution.col_value[i];
    }
  }

  if (col_dual != nullptr) {
    for (size_t i = 0; i < solution.col_dual.size(); i++) {
      col_dual[i] = solution.col_dual[i];
    }
  }

  if (row_value != nullptr) {
    for (size_t i = 0; i < solution.row_value.size(); i++) {
      row_value[i] = solution.row_value[i];
    }
  }

  if (row_dual != nullptr) {
    for (size_t i = 0; i < solution.row_dual.size(); i++) {
      row_dual[i] = solution.row_dual[i];
    }
  }
  return kHighsStatusOk;
}

HighsInt Highs_getBasis(const void* highs, HighsInt* col_status,
                        HighsInt* row_status) {
  const HighsBasis& basis = ((Highs*)highs)->getBasis();
  if (!basis.valid) return kHighsStatusError;
  for (size_t i = 0; i < basis.col_status.size(); i++) {
    col_status[i] = static_cast<HighsInt>(basis.col_status[i]);
  }

  for (size_t i = 0; i < basis.row_status.size(); i++) {
    row_status[i] = static_cast<HighsInt>(basis.row_status[i]);
  }
  return kHighsStatusOk;
}

HighsInt Highs_getModelStatus(const void* highs) {
  return (HighsInt)((Highs*)highs)->getModelStatus();
}

HighsInt Highs_getDualRay(const void* highs, HighsInt* has_dual_ray,
                          double* dual_ray_value) {
  bool v;
  HighsInt retcode = (HighsInt)((Highs*)highs)->getDualRay(v, dual_ray_value);
  *has_dual_ray = (HighsInt)v;
  return retcode;
}

HighsInt Highs_getDualUnboundednessDirection(
    const void* highs, HighsInt* has_dual_unboundedness_direction,
    double* dual_unboundedness_direction_value) {
  bool v;
  HighsInt retcode = (HighsInt)((Highs*)highs)
                         ->getDualUnboundednessDirection(
                             v, dual_unboundedness_direction_value);
  *has_dual_unboundedness_direction = (HighsInt)v;
  return retcode;
}

HighsInt Highs_getPrimalRay(const void* highs, HighsInt* has_primal_ray,
                            double* primal_ray_value) {
  bool v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getPrimalRay(v, primal_ray_value);
  *has_primal_ray = (HighsInt)v;
  return retcode;
}

double Highs_getObjectiveValue(const void* highs) {
  return ((Highs*)highs)->getObjectiveValue();
}

HighsInt Highs_getBasicVariables(const void* highs, HighsInt* basic_variables) {
  return (HighsInt)((Highs*)highs)->getBasicVariables(basic_variables);
}

HighsInt Highs_getBasisInverseRow(const void* highs, const HighsInt row,
                                  double* row_vector, HighsInt* row_num_nz,
                                  HighsInt* row_index) {
  return (HighsInt)((Highs*)highs)
      ->getBasisInverseRow(row, row_vector, row_num_nz, row_index);
}

HighsInt Highs_getBasisInverseCol(const void* highs, const HighsInt col,
                                  double* col_vector, HighsInt* col_num_nz,
                                  HighsInt* col_index) {
  return (HighsInt)((Highs*)highs)
      ->getBasisInverseCol(col, col_vector, col_num_nz, col_index);
}

HighsInt Highs_getBasisSolve(const void* highs, const double* rhs,
                             double* solution_vector, HighsInt* solution_num_nz,
                             HighsInt* solution_index) {
  return (HighsInt)((Highs*)highs)
      ->getBasisSolve(rhs, solution_vector, solution_num_nz, solution_index);
}

HighsInt Highs_getBasisTransposeSolve(const void* highs, const double* rhs,
                                      double* solution_vector,
                                      HighsInt* solution_nz,
                                      HighsInt* solution_index) {
  return (HighsInt)((Highs*)highs)
      ->getBasisTransposeSolve(rhs, solution_vector, solution_nz,
                               solution_index);
}

HighsInt Highs_getReducedRow(const void* highs, const HighsInt row,
                             double* row_vector, HighsInt* row_num_nz,
                             HighsInt* row_index) {
  return (HighsInt)((Highs*)highs)
      ->getReducedRow(row, row_vector, row_num_nz, row_index);
}

HighsInt Highs_getReducedColumn(const void* highs, const HighsInt col,
                                double* col_vector, HighsInt* col_num_nz,
                                HighsInt* col_index) {
  return (HighsInt)((Highs*)highs)
      ->getReducedColumn(col, col_vector, col_num_nz, col_index);
}

HighsInt Highs_setBasis(void* highs, const HighsInt* col_status,
                        const HighsInt* row_status) {
  HighsBasis basis;
  const HighsInt num__col = Highs_getNumCol(highs);
  if (num__col > 0) {
    basis.col_status.resize(num__col);
    for (HighsInt i = 0; i < num__col; i++) {
      if (col_status[i] == (HighsInt)HighsBasisStatus::kLower) {
        basis.col_status[i] = HighsBasisStatus::kLower;
      } else if (col_status[i] == (HighsInt)HighsBasisStatus::kBasic) {
        basis.col_status[i] = HighsBasisStatus::kBasic;
      } else if (col_status[i] == (HighsInt)HighsBasisStatus::kUpper) {
        basis.col_status[i] = HighsBasisStatus::kUpper;
      } else if (col_status[i] == (HighsInt)HighsBasisStatus::kZero) {
        basis.col_status[i] = HighsBasisStatus::kZero;
      } else if (col_status[i] == (HighsInt)HighsBasisStatus::kNonbasic) {
        basis.col_status[i] = HighsBasisStatus::kNonbasic;
      } else {
        return (HighsInt)HighsStatus::kError;
      }
    }
  }
  const HighsInt num__row = Highs_getNumRow(highs);
  if (num__row > 0) {
    basis.row_status.resize(num__row);
    for (HighsInt i = 0; i < num__row; i++) {
      if (row_status[i] == (HighsInt)HighsBasisStatus::kLower) {
        basis.row_status[i] = HighsBasisStatus::kLower;
      } else if (row_status[i] == (HighsInt)HighsBasisStatus::kBasic) {
        basis.row_status[i] = HighsBasisStatus::kBasic;
      } else if (row_status[i] == (HighsInt)HighsBasisStatus::kUpper) {
        basis.row_status[i] = HighsBasisStatus::kUpper;
      } else if (row_status[i] == (HighsInt)HighsBasisStatus::kZero) {
        basis.row_status[i] = HighsBasisStatus::kZero;
      } else if (row_status[i] == (HighsInt)HighsBasisStatus::kNonbasic) {
        basis.row_status[i] = HighsBasisStatus::kNonbasic;
      } else {
        return (HighsInt)HighsStatus::kError;
      }
    }
  }
  return (HighsInt)((Highs*)highs)->setBasis(basis);
}

HighsInt Highs_setLogicalBasis(void* highs) {
  return (HighsInt)((Highs*)highs)->setBasis();
}

HighsInt Highs_setSolution(void* highs, const double* col_value,
                           const double* row_value, const double* col_dual,
                           const double* row_dual) {
  HighsSolution solution;
  const HighsInt num__col = Highs_getNumCol(highs);
  if (num__col > 0) {
    if (col_value) {
      solution.col_value.resize(num__col);
      for (HighsInt i = 0; i < num__col; i++)
        solution.col_value[i] = col_value[i];
    }
    if (col_dual) {
      solution.col_dual.resize(num__col);
      for (HighsInt i = 0; i < num__col; i++)
        solution.col_dual[i] = col_dual[i];
    }
  }
  const HighsInt num__row = Highs_getNumRow(highs);
  if (num__row > 0) {
    if (row_value) {
      solution.row_value.resize(num__row);
      for (HighsInt i = 0; i < num__row; i++)
        solution.row_value[i] = row_value[i];
    }
    if (row_dual) {
      solution.row_dual.resize(num__row);
      for (HighsInt i = 0; i < num__row; i++)
        solution.row_dual[i] = row_dual[i];
    }
  }

  return (HighsInt)((Highs*)highs)->setSolution(solution);
}

HighsInt Highs_setSparseSolution(void* highs, const HighsInt num_entries,
                                 const HighsInt* index, const double* value) {
  return (HighsInt)((Highs*)highs)->setSolution(num_entries, index, value);
}

HighsInt Highs_setCallback(void* highs, HighsCCallbackType user_callback,
                           void* user_callback_data) {
  auto status = static_cast<Highs*>(highs)->setCallback(user_callback,
                                                        user_callback_data);
  return static_cast<HighsInt>(status);
}

HighsInt Highs_startCallback(void* highs, const HighsInt callback_type) {
  return (HighsInt)((Highs*)highs)->startCallback(callback_type);
}

HighsInt Highs_stopCallback(void* highs, const HighsInt callback_type) {
  return (HighsInt)((Highs*)highs)->stopCallback(callback_type);
}

double Highs_getRunTime(const void* highs) {
  return (double)((Highs*)highs)->getRunTime();
}

HighsInt Highs_zeroAllClocks(const void* highs) {
  ((Highs*)highs)->zeroAllClocks();
  return (HighsInt)HighsStatus::kOk;
}

HighsInt Highs_addCol(void* highs, const double cost, const double lower,
                      const double upper, const HighsInt num_new_nz,
                      const HighsInt* index, const double* value) {
  return (HighsInt)((Highs*)highs)
      ->addCol(cost, lower, upper, num_new_nz, index, value);
}

HighsInt Highs_addCols(void* highs, const HighsInt num_new_col,
                       const double* costs, const double* lower,
                       const double* upper, const HighsInt num_new_nz,
                       const HighsInt* starts, const HighsInt* index,
                       const double* value) {
  return (HighsInt)((Highs*)highs)
      ->addCols(num_new_col, costs, lower, upper, num_new_nz, starts, index,
                value);
}

HighsInt Highs_addVar(void* highs, const double lower, const double upper) {
  return (HighsInt)((Highs*)highs)->addVar(lower, upper);
}

HighsInt Highs_addVars(void* highs, const HighsInt num_new_var,
                       const double* lower, const double* upper) {
  return (HighsInt)((Highs*)highs)->addVars(num_new_var, lower, upper);
}

HighsInt Highs_addRow(void* highs, const double lower, const double upper,
                      const HighsInt num_new_nz, const HighsInt* index,
                      const double* value) {
  return (HighsInt)((Highs*)highs)
      ->addRow(lower, upper, num_new_nz, index, value);
}

HighsInt Highs_addRows(void* highs, const HighsInt num_new_row,
                       const double* lower, const double* upper,
                       const HighsInt num_new_nz, const HighsInt* starts,
                       const HighsInt* index, const double* value) {
  return (HighsInt)((Highs*)highs)
      ->addRows(num_new_row, lower, upper, num_new_nz, starts, index, value);
}

HighsInt Highs_ensureColwise(void* highs) {
  return (HighsInt)((Highs*)highs)->ensureColwise();
}

HighsInt Highs_ensureRowwise(void* highs) {
  return (HighsInt)((Highs*)highs)->ensureRowwise();
}

HighsInt Highs_changeObjectiveSense(void* highs, const HighsInt sense) {
  ObjSense pass_sense = ObjSense::kMinimize;
  if (sense == (HighsInt)ObjSense::kMaximize) pass_sense = ObjSense::kMaximize;
  return (HighsInt)((Highs*)highs)->changeObjectiveSense(pass_sense);
}

HighsInt Highs_changeObjectiveOffset(void* highs, const double offset) {
  return (HighsInt)((Highs*)highs)->changeObjectiveOffset(offset);
}

HighsInt Highs_changeColIntegrality(void* highs, const HighsInt col,
                                    const HighsInt integrality) {
  return (HighsInt)((Highs*)highs)
      ->changeColIntegrality(col, (HighsVarType)integrality);
}

HighsInt Highs_changeColsIntegralityByRange(void* highs,
                                            const HighsInt from_col,
                                            const HighsInt to_col,
                                            const HighsInt* integrality) {
  vector<HighsVarType> pass_integrality;
  HighsInt num_ix = to_col - from_col + 1;
  if (num_ix > 0) {
    pass_integrality.resize(num_ix);
    for (HighsInt ix = 0; ix < num_ix; ix++) {
      pass_integrality[ix] = (HighsVarType)integrality[ix];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(from_col, to_col, pass_integrality.data());
}

HighsInt Highs_changeColsIntegralityBySet(void* highs,
                                          const HighsInt num_set_entries,
                                          const HighsInt* set,
                                          const HighsInt* integrality) {
  vector<HighsVarType> pass_integrality;
  if (num_set_entries > 0) {
    pass_integrality.resize(num_set_entries);
    for (HighsInt ix = 0; ix < num_set_entries; ix++) {
      pass_integrality[ix] = (HighsVarType)integrality[ix];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(num_set_entries, set, pass_integrality.data());
}

HighsInt Highs_changeColsIntegralityByMask(void* highs, const HighsInt* mask,
                                           const HighsInt* integrality) {
  const HighsInt num__col = Highs_getNumCol(highs);
  vector<HighsVarType> pass_integrality;
  if (num__col > 0) {
    pass_integrality.resize(num__col);
    for (HighsInt iCol = 0; iCol < num__col; iCol++) {
      pass_integrality[iCol] = (HighsVarType)integrality[iCol];
    }
  }
  return (HighsInt)((Highs*)highs)
      ->changeColsIntegrality(mask, pass_integrality.data());
}

HighsInt Highs_clearIntegrality(void* highs) {
  return (HighsInt)((Highs*)highs)->clearIntegrality();
}

HighsInt Highs_changeColCost(void* highs, const HighsInt col,
                             const double cost) {
  return (HighsInt)((Highs*)highs)->changeColCost(col, cost);
}

HighsInt Highs_changeColsCostByRange(void* highs, const HighsInt from_col,
                                     const HighsInt to_col,
                                     const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(from_col, to_col, cost);
}

HighsInt Highs_changeColsCostBySet(void* highs, const HighsInt num_set_entries,
                                   const HighsInt* set, const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(num_set_entries, set, cost);
}

HighsInt Highs_changeColsCostByMask(void* highs, const HighsInt* mask,
                                    const double* cost) {
  return (HighsInt)((Highs*)highs)->changeColsCost(mask, cost);
}

HighsInt Highs_changeColBounds(void* highs, const HighsInt col,
                               const double lower, const double upper) {
  return (HighsInt)((Highs*)highs)->changeColBounds(col, lower, upper);
}

HighsInt Highs_changeColsBoundsByRange(void* highs, const HighsInt from_col,
                                       const HighsInt to_col,
                                       const double* lower,
                                       const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeColsBounds(from_col, to_col, lower, upper);
}

HighsInt Highs_changeColsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeColsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeColsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return (HighsInt)((Highs*)highs)->changeColsBounds(mask, lower, upper);
}

HighsInt Highs_changeRowBounds(void* highs, const HighsInt row,
                               const double lower, const double upper) {
  return (HighsInt)((Highs*)highs)->changeRowBounds(row, lower, upper);
}

HighsInt Highs_changeRowsBoundsByRange(void* highs, const HighsInt from_row,
                                       const HighsInt to_row,
                                       const double* lower,
                                       const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeRowsBounds(from_row, to_row, lower, upper);
}

HighsInt Highs_changeRowsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper) {
  return (HighsInt)((Highs*)highs)
      ->changeRowsBounds(num_set_entries, set, lower, upper);
}

HighsInt Highs_changeRowsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower,
                                      const double* upper) {
  return (HighsInt)((Highs*)highs)->changeRowsBounds(mask, lower, upper);
}

HighsInt Highs_changeCoeff(void* highs, const HighsInt row, const HighsInt col,
                           const double value) {
  return (HighsInt)((Highs*)highs)->changeCoeff(row, col, value);
}

HighsInt Highs_getObjectiveSense(const void* highs, HighsInt* sense) {
  ObjSense get_sense;
  HighsStatus status = ((Highs*)highs)->getObjectiveSense(get_sense);
  *sense = (HighsInt)get_sense;
  return (HighsInt)status;
}

HighsInt Highs_getObjectiveOffset(const void* highs, double* offset) {
  return (HighsInt)((Highs*)highs)->getObjectiveOffset(*offset);
}

HighsInt Highs_getColsByRange(const void* highs, const HighsInt from_col,
                              const HighsInt to_col, HighsInt* num_col,
                              double* costs, double* lower, double* upper,
                              HighsInt* num_nz, HighsInt* matrix_start,
                              HighsInt* matrix_index, double* matrix_value) {
  HighsInt local_num_col, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getCols(from_col, to_col, local_num_col, costs, lower, upper,
                    local_num_nz, matrix_start, matrix_index, matrix_value);
  *num_col = local_num_col;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

HighsInt Highs_getColsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_col,
                            double* costs, double* lower, double* upper,
                            HighsInt* num_nz, HighsInt* matrix_start,
                            HighsInt* matrix_index, double* matrix_value) {
  HighsInt local_num_col, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getCols(num_set_entries, set, local_num_col, costs, lower, upper,
                    local_num_nz, matrix_start, matrix_index, matrix_value);
  *num_col = local_num_col;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

HighsInt Highs_getColsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_col, double* costs, double* lower,
                             double* upper, HighsInt* num_nz,
                             HighsInt* matrix_start, HighsInt* matrix_index,
                             double* matrix_value) {
  HighsInt local_num_col, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getCols(mask, local_num_col, costs, lower, upper, local_num_nz,
                    matrix_start, matrix_index, matrix_value);
  *num_col = local_num_col;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsByRange(const void* highs, const HighsInt from_row,
                              const HighsInt to_row, HighsInt* num_row,
                              double* lower, double* upper, HighsInt* num_nz,
                              HighsInt* matrix_start, HighsInt* matrix_index,
                              double* matrix_value) {
  HighsInt local_num_row, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getRows(from_row, to_row, local_num_row, lower, upper, local_num_nz,
                    matrix_start, matrix_index, matrix_value);
  *num_row = local_num_row;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_row,
                            double* lower, double* upper, HighsInt* num_nz,
                            HighsInt* matrix_start, HighsInt* matrix_index,
                            double* matrix_value) {
  HighsInt local_num_row, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getRows(num_set_entries, set, local_num_row, lower, upper,
                    local_num_nz, matrix_start, matrix_index, matrix_value);
  *num_row = local_num_row;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

HighsInt Highs_getRowsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_row, double* lower, double* upper,
                             HighsInt* num_nz, HighsInt* matrix_start,
                             HighsInt* matrix_index, double* matrix_value) {
  HighsInt local_num_row, local_num_nz;
  HighsStatus status =
      ((Highs*)highs)
          ->getRows(mask, local_num_row, lower, upper, local_num_nz,
                    matrix_start, matrix_index, matrix_value);
  *num_row = local_num_row;
  *num_nz = local_num_nz;
  return (HighsInt)status;
}

static HighsInt Highs_getHighsLpColOrRowName(const void* highs,
                                             const HighsLp& lp,
                                             const bool is_col,
                                             const HighsInt index, char* name) {
  std::string name_v;
  HighsStatus status =
      ((Highs*)highs)->getColOrRowName(lp, is_col, index, name_v);
  if (status == HighsStatus::kError) return kHighsStatusError;
  strcpy(name, name_v.c_str());
  return kHighsStatusOk;
}

HighsInt Highs_getRowName(const void* highs, const HighsInt row, char* name) {
  return Highs_getHighsLpColOrRowName(highs, ((Highs*)highs)->getLp(), false,
                                      row, name);
}

HighsInt Highs_getRowByName(const void* highs, const char* name,
                            HighsInt* row) {
  HighsInt local_row;
  HighsInt retcode = (HighsInt)((Highs*)highs)->getRowByName(name, local_row);
  *row = local_row;
  return retcode;
}

HighsInt Highs_getColName(const void* highs, const HighsInt col, char* name) {
  return Highs_getHighsLpColOrRowName(highs, ((Highs*)highs)->getLp(), true,
                                      col, name);
}

HighsInt Highs_getColByName(const void* highs, const char* name,
                            HighsInt* col) {
  HighsInt local_col;
  HighsInt retcode = (HighsInt)((Highs*)highs)->getColByName(name, local_col);
  *col = local_col;
  return retcode;
}

HighsInt Highs_getColIntegrality(const void* highs, const HighsInt col,
                                 HighsInt* integrality) {
  HighsVarType integrality_v;
  HighsInt retcode =
      (HighsInt)((Highs*)highs)->getColIntegrality(col, integrality_v);
  *integrality = HighsInt(integrality_v);
  return retcode;
}

HighsInt Highs_deleteColsByRange(void* highs, const HighsInt from_col,
                                 const HighsInt to_col) {
  return (HighsInt)((Highs*)highs)->deleteCols(from_col, to_col);
}

HighsInt Highs_deleteColsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return (HighsInt)((Highs*)highs)->deleteCols(num_set_entries, set);
}

HighsInt Highs_deleteColsByMask(void* highs, HighsInt* mask) {
  return (HighsInt)((Highs*)highs)->deleteCols(mask);
}

HighsInt Highs_deleteRowsByRange(void* highs, const HighsInt from_row,
                                 const HighsInt to_row) {
  return (HighsInt)((Highs*)highs)->deleteRows(from_row, to_row);
}

HighsInt Highs_deleteRowsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set) {
  return (HighsInt)((Highs*)highs)->deleteRows(num_set_entries, set);
}

HighsInt Highs_deleteRowsByMask(void* highs, HighsInt* mask) {
  return (HighsInt)((Highs*)highs)->deleteRows(mask);
}

HighsInt Highs_scaleCol(void* highs, const HighsInt col,
                        const double scaleval) {
  return (HighsInt)((Highs*)highs)->scaleCol(col, scaleval);
}

HighsInt Highs_scaleRow(void* highs, const HighsInt row,
                        const double scaleval) {
  return (HighsInt)((Highs*)highs)->scaleRow(row, scaleval);
}

double Highs_getInfinity(const void* highs) {
  return ((Highs*)highs)->getInfinity();
}

HighsInt Highs_getSizeofHighsInt(const void* highs) {
  return ((Highs*)highs)->getSizeofHighsInt();
}

HighsInt Highs_getNumCol(const void* highs) {
  return ((Highs*)highs)->getNumCol();
}

HighsInt Highs_getNumRow(const void* highs) {
  return ((Highs*)highs)->getNumRow();
}

HighsInt Highs_getNumNz(const void* highs) {
  return ((Highs*)highs)->getNumNz();
}

HighsInt Highs_getHessianNumNz(const void* highs) {
  return ((Highs*)highs)->getHessianNumNz();
}

HighsInt Highs_getPresolvedNumCol(const void* highs) {
  return ((Highs*)highs)->getPresolvedLp().num_col_;
}

HighsInt Highs_getPresolvedNumRow(const void* highs) {
  return ((Highs*)highs)->getPresolvedLp().num_row_;
}

HighsInt Highs_getPresolvedNumNz(const void* highs) {
  return ((Highs*)highs)->getPresolvedLp().numNz();
}

// Gets pointers to all the public data members of HighsLp: avoids
// duplicate code in Highs_getModel, Highs_getPresolvedLp,
static HighsInt Highs_getHighsLpData(const HighsLp& lp, const HighsInt a_format,
                                     HighsInt* num_col, HighsInt* num_row,
                                     HighsInt* num_nz, HighsInt* sense,
                                     double* offset, double* col_cost,
                                     double* col_lower, double* col_upper,
                                     double* row_lower, double* row_upper,
                                     HighsInt* a_start, HighsInt* a_index,
                                     double* a_value, HighsInt* integrality) {
  const MatrixFormat desired_a_format =
      a_format == HighsInt(MatrixFormat::kColwise) ? MatrixFormat::kColwise
                                                   : MatrixFormat::kRowwise;
  *sense = (HighsInt)lp.sense_;
  *offset = lp.offset_;
  *num_col = lp.num_col_;
  *num_row = lp.num_row_;
  *num_nz = 0;  // In case one of the matrix dimensions is zero
  if (*num_col > 0) {
    if (col_cost)
      memcpy(col_cost, lp.col_cost_.data(), *num_col * sizeof(double));
    if (col_lower)
      memcpy(col_lower, lp.col_lower_.data(), *num_col * sizeof(double));
    if (col_upper)
      memcpy(col_upper, lp.col_upper_.data(), *num_col * sizeof(double));
  }
  if (*num_row > 0) {
    if (row_lower)
      memcpy(row_lower, lp.row_lower_.data(), *num_row * sizeof(double));
    if (row_upper)
      memcpy(row_upper, lp.row_upper_.data(), *num_row * sizeof(double));
  }

  // Nothing to do if one of the matrix dimensions is zero
  if (*num_col > 0 && *num_row > 0) {
    // Determine the desired orientation and number of start entries to
    // be copied
    const HighsInt num_start_entries =
        desired_a_format == MatrixFormat::kColwise ? *num_col : *num_row;
    if ((desired_a_format == MatrixFormat::kColwise &&
         lp.a_matrix_.isColwise()) ||
        (desired_a_format == MatrixFormat::kRowwise &&
         lp.a_matrix_.isRowwise())) {
      // Incumbent format is OK
      *num_nz = lp.numNz();
      if (a_start)
        memcpy(a_start, lp.a_matrix_.start_.data(),
               num_start_entries * sizeof(HighsInt));
      if (a_index)
        memcpy(a_index, lp.a_matrix_.index_.data(), *num_nz * sizeof(HighsInt));
      if (a_value)
        memcpy(a_value, lp.a_matrix_.value_.data(), *num_nz * sizeof(double));
    } else {
      // Take a copy and transpose it
      HighsSparseMatrix local_matrix = lp.a_matrix_;
      if (desired_a_format == MatrixFormat::kColwise) {
        assert(local_matrix.isRowwise());
        local_matrix.ensureColwise();
      } else {
        assert(local_matrix.isColwise());
        local_matrix.ensureRowwise();
      }
      *num_nz = local_matrix.numNz();
      if (a_start)
        memcpy(a_start, local_matrix.start_.data(),
               num_start_entries * sizeof(HighsInt));
      if (a_index)
        memcpy(a_index, local_matrix.index_.data(), *num_nz * sizeof(HighsInt));
      if (a_value)
        memcpy(a_value, local_matrix.value_.data(), *num_nz * sizeof(double));
    }
  }
  if (HighsInt(lp.integrality_.size()) && integrality) {
    for (int iCol = 0; iCol < *num_col; iCol++)
      integrality[iCol] = HighsInt(lp.integrality_[iCol]);
  }
  return kHighsStatusOk;
}

HighsInt Highs_getModel(const void* highs, const HighsInt a_format,
                        const HighsInt q_format, HighsInt* num_col,
                        HighsInt* num_row, HighsInt* num_nz, HighsInt* q_num_nz,
                        HighsInt* sense, double* offset, double* col_cost,
                        double* col_lower, double* col_upper, double* row_lower,
                        double* row_upper, HighsInt* a_start, HighsInt* a_index,
                        double* a_value, HighsInt* q_start, HighsInt* q_index,
                        double* q_value, HighsInt* integrality) {
  HighsInt return_status = Highs_getHighsLpData(
      ((Highs*)highs)->getLp(), a_format, num_col, num_row, num_nz, sense,
      offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start,
      a_index, a_value, integrality);
  if (return_status != kHighsStatusOk) return return_status;
  const HighsHessian& hessian = ((Highs*)highs)->getModel().hessian_;
  if (hessian.dim_ > 0) {
    *q_num_nz = hessian.start_[*num_col];
    if (q_start)
      memcpy(q_start, hessian.start_.data(), *num_col * sizeof(HighsInt));
    if (q_index)
      memcpy(q_index, hessian.index_.data(), *q_num_nz * sizeof(HighsInt));
    if (q_value)
      memcpy(q_value, hessian.value_.data(), *q_num_nz * sizeof(double));
  }
  return kHighsStatusOk;
}

HighsInt Highs_getLp(const void* highs, const HighsInt a_format,
                     HighsInt* num_col, HighsInt* num_row, HighsInt* num_nz,
                     HighsInt* sense, double* offset, double* col_cost,
                     double* col_lower, double* col_upper, double* row_lower,
                     double* row_upper, HighsInt* a_start, HighsInt* a_index,
                     double* a_value, HighsInt* integrality) {
  return Highs_getHighsLpData(((Highs*)highs)->getLp(), a_format, num_col,
                              num_row, num_nz, sense, offset, col_cost,
                              col_lower, col_upper, row_lower, row_upper,
                              a_start, a_index, a_value, integrality);
}

HighsInt Highs_getFixedLp(const void* highs, const HighsInt a_format,
                          HighsInt* num_col, HighsInt* num_row,
                          HighsInt* num_nz, HighsInt* sense, double* offset,
                          double* col_cost, double* col_lower,
                          double* col_upper, double* row_lower,
                          double* row_upper, HighsInt* a_start,
                          HighsInt* a_index, double* a_value) {
  HighsLp lp;
  HighsInt status = HighsInt(((Highs*)highs)->getFixedLp(lp));
  if (status == kHighsStatusError) return status;
  HighsInt* integrality = nullptr;
  Highs_getHighsLpData(lp, a_format, num_col, num_row, num_nz, sense, offset,
                       col_cost, col_lower, col_upper, row_lower, row_upper,
                       a_start, a_index, a_value, integrality);
  return status;
}

HighsInt Highs_getPresolvedLp(const void* highs, const HighsInt a_format,
                              HighsInt* num_col, HighsInt* num_row,
                              HighsInt* num_nz, HighsInt* sense, double* offset,
                              double* col_cost, double* col_lower,
                              double* col_upper, double* row_lower,
                              double* row_upper, HighsInt* a_start,
                              HighsInt* a_index, double* a_value,
                              HighsInt* integrality) {
  return Highs_getHighsLpData(((Highs*)highs)->getPresolvedLp(), a_format,
                              num_col, num_row, num_nz, sense, offset, col_cost,
                              col_lower, col_upper, row_lower, row_upper,
                              a_start, a_index, a_value, integrality);
}

HighsInt Highs_getPresolvedColName(const void* highs, const HighsInt col,
                                   char* name) {
  return Highs_getHighsLpColOrRowName(highs, ((Highs*)highs)->getPresolvedLp(),
                                      true, col, name);
}

HighsInt Highs_getPresolvedRowName(const void* highs, const HighsInt row,
                                   char* name) {
  return Highs_getHighsLpColOrRowName(highs, ((Highs*)highs)->getPresolvedLp(),
                                      false, row, name);
}

HighsInt Highs_getIis(void* highs, HighsInt* iis_num_col, HighsInt* iis_num_row,
                      HighsInt* col_index, HighsInt* row_index,
                      HighsInt* col_bound, HighsInt* row_bound,
                      HighsInt* col_status, HighsInt* row_status) {
  HighsIis iis;
  HighsInt status = (HighsInt)((Highs*)highs)->getIis(iis);
  if (status == (HighsInt)HighsStatus::kError) return status;
  *iis_num_col = iis.col_index_.size();
  *iis_num_row = iis.row_index_.size();
  if (col_index != nullptr) {
    for (size_t i = 0; i < static_cast<size_t>(*iis_num_col); i++) {
      col_index[i] = iis.col_index_[i];
    }
  }
  if (row_index != nullptr) {
    for (size_t i = 0; i < static_cast<size_t>(*iis_num_row); i++) {
      row_index[i] = iis.row_index_[i];
    }
  }
  if (col_bound != nullptr) {
    for (size_t i = 0; i < static_cast<size_t>(*iis_num_col); i++) {
      col_bound[i] = iis.col_bound_[i];
    }
  }
  if (row_bound != nullptr) {
    for (size_t i = 0; i < static_cast<size_t>(*iis_num_row); i++) {
      row_bound[i] = iis.row_bound_[i];
    }
  }
  if (col_status != nullptr) {
    for (size_t i = 0;
         i < static_cast<size_t>(((Highs*)highs)->getLp().num_col_); i++) {
      col_status[i] = iis.col_status_[i];
    }
  }
  if (row_status != nullptr) {
    for (size_t i = 0;
         i < static_cast<size_t>(((Highs*)highs)->getLp().num_row_); i++) {
      row_status[i] = iis.row_status_[i];
    }
  }
  return status;
}

HighsInt Highs_getIisLp(const void* highs, const HighsInt a_format,
                        HighsInt* num_col, HighsInt* num_row, HighsInt* num_nz,
                        HighsInt* sense, double* offset, double* col_cost,
                        double* col_lower, double* col_upper, double* row_lower,
                        double* row_upper, HighsInt* a_start, HighsInt* a_index,
                        double* a_value, HighsInt* integrality) {
  return Highs_getHighsLpData(((Highs*)highs)->getIisLp(), a_format, num_col,
                              num_row, num_nz, sense, offset, col_cost,
                              col_lower, col_upper, row_lower, row_upper,
                              a_start, a_index, a_value, integrality);
}

HighsInt Highs_crossover(void* highs, const HighsInt num_col,
                         const HighsInt num_row, const double* col_value,
                         const double* col_dual, const double* row_dual) {
  HighsSolution solution;
  if (col_value) {
    solution.value_valid = true;
    solution.col_value.resize(num_col);
    for (int col = 0; col < num_col; col++)
      solution.col_value[col] = col_value[col];
  }

  if (col_dual && row_dual) {
    solution.dual_valid = true;
    solution.col_dual.resize(num_col);
    solution.row_dual.resize(num_row);
    for (int col = 0; col < num_col; col++)
      solution.col_dual[col] = col_dual[col];
    for (int row = 0; row < num_row; row++)
      solution.row_dual[row] = row_dual[row];
  }

  return (HighsInt)((Highs*)highs)->crossover(solution);
}

HighsInt Highs_getRanging(
    void* highs,
    //
    double* col_cost_up_value, double* col_cost_up_objective,
    HighsInt* col_cost_up_in_var, HighsInt* col_cost_up_ou_var,
    double* col_cost_dn_value, double* col_cost_dn_objective,
    HighsInt* col_cost_dn_in_var, HighsInt* col_cost_dn_ou_var,
    double* col_bound_up_value, double* col_bound_up_objective,
    HighsInt* col_bound_up_in_var, HighsInt* col_bound_up_ou_var,
    double* col_bound_dn_value, double* col_bound_dn_objective,
    HighsInt* col_bound_dn_in_var, HighsInt* col_bound_dn_ou_var,
    double* row_bound_up_value, double* row_bound_up_objective,
    HighsInt* row_bound_up_in_var, HighsInt* row_bound_up_ou_var,
    double* row_bound_dn_value, double* row_bound_dn_objective,
    HighsInt* row_bound_dn_in_var, HighsInt* row_bound_dn_ou_var) {
  HighsRanging ranging;
  HighsInt status = (HighsInt)((Highs*)highs)->getRanging(ranging);
  if (status == (HighsInt)HighsStatus::kError) return status;
  HighsInt num_col = ((Highs*)highs)->getNumCol();
  HighsInt num_row = ((Highs*)highs)->getNumRow();
  if (col_cost_up_value)
    memcpy(col_cost_up_value, ranging.col_cost_up.value_.data(),
           num_col * sizeof(double));
  if (col_cost_up_objective)
    memcpy(col_cost_up_objective, ranging.col_cost_up.objective_.data(),
           num_col * sizeof(double));
  if (col_cost_up_in_var)
    memcpy(col_cost_up_in_var, ranging.col_cost_up.in_var_.data(),
           num_col * sizeof(HighsInt));
  if (col_cost_up_ou_var)
    memcpy(col_cost_up_ou_var, ranging.col_cost_up.ou_var_.data(),
           num_col * sizeof(HighsInt));

  if (col_cost_dn_value)
    memcpy(col_cost_dn_value, ranging.col_cost_dn.value_.data(),
           num_col * sizeof(double));
  if (col_cost_dn_objective)
    memcpy(col_cost_dn_objective, ranging.col_cost_dn.objective_.data(),
           num_col * sizeof(double));
  if (col_cost_dn_in_var)
    memcpy(col_cost_dn_in_var, ranging.col_cost_dn.in_var_.data(),
           num_col * sizeof(HighsInt));
  if (col_cost_dn_ou_var)
    memcpy(col_cost_dn_ou_var, ranging.col_cost_dn.ou_var_.data(),
           num_col * sizeof(HighsInt));

  if (col_bound_up_value)
    memcpy(col_bound_up_value, ranging.col_bound_up.value_.data(),
           num_col * sizeof(double));
  if (col_bound_up_objective)
    memcpy(col_bound_up_objective, ranging.col_bound_up.objective_.data(),
           num_col * sizeof(double));
  if (col_bound_up_in_var)
    memcpy(col_bound_up_in_var, ranging.col_bound_up.in_var_.data(),
           num_col * sizeof(HighsInt));
  if (col_bound_up_ou_var)
    memcpy(col_bound_up_ou_var, ranging.col_bound_up.ou_var_.data(),
           num_col * sizeof(HighsInt));

  if (col_bound_dn_value)
    memcpy(col_bound_dn_value, ranging.col_bound_dn.value_.data(),
           num_col * sizeof(double));
  if (col_bound_dn_objective)
    memcpy(col_bound_dn_objective, ranging.col_bound_dn.objective_.data(),
           num_col * sizeof(double));
  if (col_bound_dn_in_var)
    memcpy(col_bound_dn_in_var, ranging.col_bound_dn.in_var_.data(),
           num_col * sizeof(HighsInt));
  if (col_bound_dn_ou_var)
    memcpy(col_bound_dn_ou_var, ranging.col_bound_dn.ou_var_.data(),
           num_col * sizeof(HighsInt));

  if (row_bound_up_value)
    memcpy(row_bound_up_value, ranging.row_bound_up.value_.data(),
           num_row * sizeof(double));
  if (row_bound_up_objective)
    memcpy(row_bound_up_objective, ranging.row_bound_up.objective_.data(),
           num_row * sizeof(double));
  if (row_bound_up_in_var)
    memcpy(row_bound_up_in_var, ranging.row_bound_up.in_var_.data(),
           num_row * sizeof(HighsInt));
  if (row_bound_up_ou_var)
    memcpy(row_bound_up_ou_var, ranging.row_bound_up.ou_var_.data(),
           num_row * sizeof(HighsInt));

  if (row_bound_dn_value)
    memcpy(row_bound_dn_value, ranging.row_bound_dn.value_.data(),
           num_row * sizeof(double));
  if (row_bound_dn_objective)
    memcpy(row_bound_dn_objective, ranging.row_bound_dn.objective_.data(),
           num_row * sizeof(double));
  if (row_bound_dn_in_var)
    memcpy(row_bound_dn_in_var, ranging.row_bound_dn.in_var_.data(),
           num_row * sizeof(HighsInt));
  if (row_bound_dn_ou_var)
    memcpy(row_bound_dn_ou_var, ranging.row_bound_dn.ou_var_.data(),
           num_row * sizeof(HighsInt));

  return status;
}

HighsInt Highs_feasibilityRelaxation(void* highs,
                                     const double global_lower_penalty,
                                     const double global_upper_penalty,
                                     const double global_rhs_penalty,
                                     const double* local_lower_penalty,
                                     const double* local_upper_penalty,
                                     const double* local_rhs_penalty) {
  return (HighsInt)((Highs*)highs)
      ->feasibilityRelaxation(global_lower_penalty, global_upper_penalty,
                              global_rhs_penalty, local_lower_penalty,
                              local_upper_penalty, local_rhs_penalty);
}

void Highs_resetThreadScheduler(HighsInt blocking) {
  Highs::resetThreadScheduler(blocking != 0);
}

const void* Highs_getCallbackDataOutItem(const HighsCallbackDataOut* data_out,
                                         const char* item_name) {
  // Accessor function for HighsCallbackDataOut
  //
  // Remember that pointers in HighsCallbackDataOut don't need to be referenced!
  if (!strcmp(item_name, kHighsCallbackDataOutLogTypeName)) {
    return (void*)(&data_out->log_type);
  } else if (!strcmp(item_name, kHighsCallbackDataOutRunningTimeName)) {
    return (void*)(&data_out->running_time);
  } else if (!strcmp(item_name,
                     kHighsCallbackDataOutSimplexIterationCountName)) {
    return (void*)(&data_out->simplex_iteration_count);
  } else if (!strcmp(item_name, kHighsCallbackDataOutIpmIterationCountName)) {
    return (void*)(&data_out->ipm_iteration_count);
  } else if (!strcmp(item_name, kHighsCallbackDataOutPdlpIterationCountName)) {
    return (void*)(&data_out->pdlp_iteration_count);
  } else if (!strcmp(item_name, kHighsCallbackDataOutQpasmIterationCountName)) {
    return (void*)(&data_out->qpasm_iteration_count);
  } else if (!strcmp(item_name,
                     kHighsCallbackDataOutObjectiveFunctionValueName)) {
    return (void*)(&data_out->objective_function_value);
  } else if (!strcmp(item_name, kHighsCallbackDataOutMipNodeCountName)) {
    return (void*)(&data_out->mip_node_count);
  } else if (!strcmp(item_name,
                     kHighsCallbackDataOutMipTotalLpIterationsName)) {
    return (void*)(&data_out->mip_total_lp_iterations);
  } else if (!strcmp(item_name, kHighsCallbackDataOutMipPrimalBoundName)) {
    return (void*)(&data_out->mip_primal_bound);
  } else if (!strcmp(item_name, kHighsCallbackDataOutMipDualBoundName)) {
    return (void*)(&data_out->mip_dual_bound);
  } else if (!strcmp(item_name, kHighsCallbackDataOutMipGapName)) {
    return (void*)(&data_out->mip_gap);
  } else if (!strcmp(item_name, kHighsCallbackDataOutMipSolutionName)) {
    return (void*)(data_out->mip_solution);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolNumColName)) {
    return (void*)(&data_out->cutpool_num_col);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolNumCutName)) {
    return (void*)(&data_out->cutpool_num_cut);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolNumNzName)) {
    return (void*)(&data_out->cutpool_num_nz);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolStartName)) {
    return (void*)(data_out->cutpool_start);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolIndexName)) {
    return (void*)(data_out->cutpool_index);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolValueName)) {
    return (void*)(data_out->cutpool_value);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolLowerName)) {
    return (void*)(data_out->cutpool_lower);
  } else if (!strcmp(item_name, kHighsCallbackDataOutCutpoolUpperName)) {
    return (void*)(data_out->cutpool_upper);
  } else if (!strcmp(item_name, kHighsCallbackDataOutQpSolutionName)) {
    return (void*)(data_out->qp_solution);
  }
  return nullptr;
}

HighsInt Highs_setCallbackSolution(HighsCallbackDataIn* data_in,
                                   HighsInt num_entries, const double* value) {
  if (data_in != nullptr && data_in->cbdata != nullptr) {
    HighsCallbackInput* obj = static_cast<HighsCallbackInput*>(data_in->cbdata);
    return static_cast<int>(obj->setSolution(num_entries, value));
  } else
    return static_cast<int>(HighsStatus::kError);
}

HighsInt Highs_setCallbackSparseSolution(HighsCallbackDataIn* data_in,
                                         HighsInt num_entries,
                                         const HighsInt* index,
                                         const double* value) {
  if (data_in != nullptr && data_in->cbdata != nullptr) {
    HighsCallbackInput* obj = static_cast<HighsCallbackInput*>(data_in->cbdata);
    return static_cast<int>(obj->setSolution(num_entries, index, value));
  } else
    return static_cast<int>(HighsStatus::kError);
}

HighsInt Highs_repairCallbackSolution(HighsCallbackDataIn* data_in) {
  if (data_in != nullptr && data_in->cbdata != nullptr) {
    HighsCallbackInput* obj = static_cast<HighsCallbackInput*>(data_in->cbdata);
    return static_cast<int>(obj->repairSolution());
  } else
    return static_cast<int>(HighsStatus::kError);
}

// *********************
// * Deprecated methods*
// *********************

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
const char* Highs_compilationDate(void) { return "Deprecated"; }
#pragma GCC diagnostic pop

HighsInt Highs_call(const HighsInt num_col, const HighsInt num_row,
                    const HighsInt num_nz, const double* col_cost,
                    const double* col_lower, const double* col_upper,
                    const double* row_lower, const double* row_upper,
                    const HighsInt* a_start, const HighsInt* a_index,
                    const double* a_value, double* col_value, double* col_dual,
                    double* row_value, double* row_dual,
                    HighsInt* col_basis_status, HighsInt* row_basis_status,
                    HighsInt* model_status) {
  printf(
      "Method Highs_call is deprecated: alternative method is Highs_lpCall\n");
  const HighsInt aformat_columnwise = 1;
  const HighsInt sense = 1;
  const double offset = 0;
  return Highs_lpCall(num_col, num_row, num_nz, aformat_columnwise, sense,
                      offset, col_cost, col_lower, col_upper, row_lower,
                      row_upper, a_start, a_index, a_value, col_value, col_dual,
                      row_value, row_dual, col_basis_status, row_basis_status,
                      model_status);
}

HighsInt Highs_setOptionValue(void* highs, const char* option,
                              const char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setOptionValue",
                           "Highs_setStringOptionValue");
  return (HighsInt)((Highs*)highs)
      ->setOptionValue(std::string(option), std::string(value));
}

HighsInt Highs_runQuiet(void* highs) {
  ((Highs*)highs)->deprecationMessage("Highs_runQuiet", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_setHighsLogfile(void* highs, const void* logfile) {
  ((Highs*)highs)->deprecationMessage("Highs_setHighsLogfile", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_setHighsOutput(void* highs, const void* outputfile) {
  ((Highs*)highs)->deprecationMessage("Highs_setHighsOutput", "None");
  return (HighsInt)((Highs*)highs)->setOptionValue("output_flag", false);
}

HighsInt Highs_getIterationCount(const void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getIterationCount", "Highs_getIntInfoValue");
  return (HighsInt)((Highs*)highs)->getInfo().simplex_iteration_count;
}

HighsInt Highs_getSimplexIterationCount(const void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getSimplexIterationCount",
                           "Highs_getIntInfoValue");
  return (HighsInt)((Highs*)highs)->getInfo().simplex_iteration_count;
}

HighsInt Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                       const HighsInt value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsBoolOptionValue",
                           "Highs_setBoolOptionValue");
  return Highs_setBoolOptionValue(highs, option, value);
}

HighsInt Highs_setHighsIntOptionValue(void* highs, const char* option,
                                      const HighsInt value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsIntOptionValue",
                           "Highs_setIntOptionValue");
  return Highs_setIntOptionValue(highs, option, value);
}

HighsInt Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                         const double value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsDoubleOptionValue",
                           "Highs_setDoubleOptionValue");
  return Highs_setDoubleOptionValue(highs, option, value);
}

HighsInt Highs_setHighsStringOptionValue(void* highs, const char* option,
                                         const char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsStringOptionValue",
                           "Highs_setStringOptionValue");
  return Highs_setStringOptionValue(highs, option, value);
}

HighsInt Highs_setHighsOptionValue(void* highs, const char* option,
                                   const char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_setHighsOptionValue", "Highs_setOptionValue");
  return Highs_setOptionValue(highs, option, value);
}

HighsInt Highs_getHighsBoolOptionValue(const void* highs, const char* option,
                                       HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsBoolOptionValue",
                           "Highs_getBoolOptionValue");
  return Highs_getBoolOptionValue(highs, option, value);
}

HighsInt Highs_getHighsIntOptionValue(const void* highs, const char* option,
                                      HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsIntOptionValue",
                           "Highs_getIntOptionValue");
  return Highs_getIntOptionValue(highs, option, value);
}

HighsInt Highs_getHighsDoubleOptionValue(const void* highs, const char* option,
                                         double* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsDoubleOptionValue",
                           "Highs_getDoubleOptionValue");
  return Highs_getDoubleOptionValue(highs, option, value);
}

HighsInt Highs_getHighsStringOptionValue(const void* highs, const char* option,
                                         char* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsStringOptionValue",
                           "Highs_getStringOptionValue");
  return Highs_getStringOptionValue(highs, option, value);
}

HighsInt Highs_getHighsOptionType(const void* highs, const char* option,
                                  HighsInt* type) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsOptionType", "Highs_getOptionType");
  return Highs_getOptionType(highs, option, type);
}

HighsInt Highs_resetHighsOptions(void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_resetHighsOptions", "Highs_resetOptions");
  return Highs_resetOptions(highs);
}

HighsInt Highs_getHighsIntInfoValue(const void* highs, const char* info,
                                    HighsInt* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsIntInfoValue",
                           "Highs_getIntInfoValue");
  return Highs_getIntInfoValue(highs, info, value);
}

HighsInt Highs_getHighsDoubleInfoValue(const void* highs, const char* info,
                                       double* value) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsDoubleInfoValue",
                           "Highs_getDoubleInfoValue");
  return Highs_getDoubleInfoValue(highs, info, value);
}

HighsInt Highs_getNumCols(const void* highs) { return Highs_getNumCol(highs); }
HighsInt Highs_getNumRows(const void* highs) { return Highs_getNumRow(highs); }

double Highs_getHighsRunTime(const void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsRunTime", "Highs_getRunTime");
  return Highs_getRunTime(highs);
}

double Highs_getHighsInfinity(const void* highs) {
  ((Highs*)highs)
      ->deprecationMessage("Highs_getHighsInfinity", "Highs_getInfinity");
  return Highs_getInfinity(highs);
}

HighsInt Highs_getScaledModelStatus(const void* highs) {
  return (HighsInt)((Highs*)highs)->getModelStatus();
}

// Deprecated in favour of resetThreadScheduler().
// Cannot call deprecationMessage here.
void Highs_resetGlobalScheduler(HighsInt blocking) {
  Highs::resetGlobalScheduler(blocking != 0);
}

HighsInt Highs_initializeExtras(const char* directory,
                                HighsInt search_if_empty) {
  std::string dir((directory != nullptr) ? directory : "");
  return (HighsInt)Highs::initializeExtras(dir, search_if_empty != 0);
}