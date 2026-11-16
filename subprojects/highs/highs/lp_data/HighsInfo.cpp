/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInfo.cpp
 * @brief
 */
#include "lp_data/HighsInfo.h"

#include <cassert>

#include "lp_data/HighsOptions.h"

void HighsInfo::invalidate() {
  valid = false;
  mip_node_count = -1;
  simplex_iteration_count = -1;
  ipm_iteration_count = -1;
  crossover_iteration_count = -1;
  pdlp_iteration_count = -1;
  qp_iteration_count = -1;
  basis_validity = kBasisValidityInvalid;
  centring_status = kCentringStatusNone;
  objective_function_value = 0;
  mip_dual_bound = 0;
  mip_gap = kHighsInf;
  max_integrality_violation = kHighsIllegalInfeasibilityMeasure;
  this->invalidateKkt();
  primal_dual_integral = -kHighsInf;
}

void HighsInfo::invalidateKkt() {
  this->invalidatePrimalKkt();
  this->invalidateDualKkt();
}

void HighsInfo::invalidatePrimalKkt() {
  primal_solution_status = kSolutionStatusNone;
  num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  num_semi_infeasibilities = kHighsIllegalInfeasibilityCount;
  max_semi_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_semi_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  num_relative_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  max_relative_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  num_primal_residual_errors = kHighsIllegalResidualCount;
  max_primal_residual_error = kHighsIllegalResidualMeasure;
  num_relative_primal_residual_errors = kHighsIllegalResidualCount;
  max_relative_primal_residual_error = kHighsIllegalResidualMeasure;
  num_complementarity_violations = kHighsIllegalComplementarityCount;
  max_complementarity_violation = kHighsIllegalComplementarityViolation;
  primal_dual_objective_error = kHighsIllegalComplementarityViolation;
  active_bound_norm = kHighsIllegalDoubleMeasure;
}

void HighsInfo::invalidateDualKkt() {
  dual_solution_status = kSolutionStatusNone;
  num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  num_relative_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  max_relative_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  num_dual_residual_errors = kHighsIllegalResidualCount;
  max_dual_residual_error = kHighsIllegalResidualMeasure;
  num_relative_dual_residual_errors = kHighsIllegalResidualCount;
  max_relative_dual_residual_error = kHighsIllegalResidualMeasure;
  num_complementarity_violations = kHighsIllegalComplementarityCount;
  max_complementarity_violation = kHighsIllegalComplementarityViolation;
  primal_dual_objective_error = kHighsIllegalComplementarityViolation;
  active_cost_norm = kHighsIllegalDoubleMeasure;
}

bool HighsInfo::equal(const HighsInfo& info_) const {
  if (info_.valid != this->valid) return false;
  if (info_.mip_node_count != this->mip_node_count) return false;
  if (info_.simplex_iteration_count != this->simplex_iteration_count)
    return false;
  if (info_.ipm_iteration_count != this->ipm_iteration_count) return false;
  if (info_.crossover_iteration_count != this->crossover_iteration_count)
    return false;
  if (info_.pdlp_iteration_count != this->pdlp_iteration_count) return false;
  if (info_.qp_iteration_count != this->qp_iteration_count) return false;
  if (info_.primal_solution_status != this->primal_solution_status)
    return false;
  if (info_.dual_solution_status != this->dual_solution_status) return false;
  if (info_.basis_validity != this->basis_validity) return false;
  if (info_.centring_status != this->centring_status) return false;
  if (info_.objective_function_value != this->objective_function_value)
    return false;
  if (info_.mip_dual_bound != this->mip_dual_bound) return false;
  if (info_.mip_gap != this->mip_gap) return false;
  if (info_.max_integrality_violation != this->max_integrality_violation)
    return false;
  if (info_.num_primal_infeasibilities != this->num_primal_infeasibilities)
    return false;
  if (info_.max_primal_infeasibility != this->max_primal_infeasibility)
    return false;
  if (info_.sum_primal_infeasibilities != this->sum_primal_infeasibilities)
    return false;
  if (info_.num_dual_infeasibilities != this->num_dual_infeasibilities)
    return false;
  if (info_.max_dual_infeasibility != this->max_dual_infeasibility)
    return false;
  if (info_.sum_dual_infeasibilities != this->sum_dual_infeasibilities)
    return false;
  if (info_.num_semi_infeasibilities != this->num_semi_infeasibilities)
    return false;
  if (info_.max_semi_infeasibility != this->max_semi_infeasibility)
    return false;
  if (info_.sum_semi_infeasibilities != this->sum_semi_infeasibilities)
    return false;
  if (info_.num_relative_primal_infeasibilities !=
      this->num_relative_primal_infeasibilities)
    return false;
  if (info_.max_relative_primal_infeasibility !=
      this->max_relative_primal_infeasibility)
    return false;
  if (info_.num_relative_dual_infeasibilities !=
      this->num_relative_dual_infeasibilities)
    return false;
  if (info_.max_relative_dual_infeasibility !=
      this->max_relative_dual_infeasibility)
    return false;
  if (info_.num_primal_residual_errors != this->num_primal_residual_errors)
    return false;
  if (info_.max_primal_residual_error != this->max_primal_residual_error)
    return false;
  if (info_.num_dual_residual_errors != this->num_dual_residual_errors)
    return false;
  if (info_.max_dual_residual_error != this->max_dual_residual_error)
    return false;
  if (info_.num_relative_primal_residual_errors !=
      this->num_relative_primal_residual_errors)
    return false;
  if (info_.max_relative_primal_residual_error !=
      this->max_relative_primal_residual_error)
    return false;
  if (info_.num_relative_dual_residual_errors !=
      this->num_relative_dual_residual_errors)
    return false;
  if (info_.max_relative_dual_residual_error !=
      this->max_relative_dual_residual_error)
    return false;
  if (info_.num_complementarity_violations !=
      this->num_complementarity_violations)
    return false;
  if (info_.max_complementarity_violation !=
      this->max_complementarity_violation)
    return false;
  if (info_.primal_dual_objective_error != this->primal_dual_objective_error)
    return false;
  if (info_.active_cost_norm != this->active_cost_norm) return false;
  if (info_.active_bound_norm != this->active_bound_norm) return false;
  return true;
}

static std::string infoEntryTypeToString(const HighsInfoType type) {
  if (type == HighsInfoType::kInt64) {
    return "int64_t";
  } else if (type == HighsInfoType::kInt) {
    return "HighsInt";
  } else {
    return "double";
  }
}

InfoStatus getInfoIndex(const HighsLogOptions& report_log_options,
                        const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        HighsInt& index) {
  HighsInt num_info = info_records.size();
  for (index = 0; index < num_info; index++)
    if (info_records[index]->name == name) return InfoStatus::kOk;
  highsLogUser(report_log_options, HighsLogType::kError,
               "getInfoIndex: Info \"%s\" is unknown\n", name.c_str());
  return InfoStatus::kUnknownInfo;
}

InfoStatus checkInfo(const HighsLogOptions& report_log_options,
                     const std::vector<InfoRecord*>& info_records) {
  bool error_found = false;
  HighsInt num_info = info_records.size();
  for (HighsInt index = 0; index < num_info; index++) {
    std::string name = info_records[index]->name;
    HighsInfoType type = info_records[index]->type;
    // Check that there are no other info with the same name
    for (HighsInt check_index = 0; check_index < num_info; check_index++) {
      if (check_index == index) continue;
      std::string check_name = info_records[check_index]->name;
      if (check_name == name) {
        highsLogUser(report_log_options, HighsLogType::kError,
                     "checkInfo: Info %" HIGHSINT_FORMAT
                     " (\"%s\") has the same name as info %" HIGHSINT_FORMAT
                     " \"%s\"\n",
                     index, name.c_str(), check_index, check_name.c_str());
        error_found = true;
      }
    }
    if (type == HighsInfoType::kInt64) {
      // Check int64_t info
      InfoRecordInt64& info = ((InfoRecordInt64*)info_records[index])[0];
      // Check that there are no other info with the same value pointers
      int64_t* value_pointer = info.value;
      for (HighsInt check_index = 0; check_index < num_info; check_index++) {
        if (check_index == index) continue;
        InfoRecordInt64& check_info =
            ((InfoRecordInt64*)info_records[check_index])[0];
        if (check_info.type == HighsInfoType::kInt64) {
          if (check_info.value == value_pointer) {
            highsLogUser(report_log_options, HighsLogType::kError,
                         "checkInfo: Info %" HIGHSINT_FORMAT
                         " (\"%s\") has the same value "
                         "pointer as info %" HIGHSINT_FORMAT " (\"%s\")\n",
                         index, info.name.c_str(), check_index,
                         check_info.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsInfoType::kInt) {
      // Check HighsInt info
      InfoRecordInt& info = ((InfoRecordInt*)info_records[index])[0];
      // Check that there are no other info with the same value pointers
      HighsInt* value_pointer = info.value;
      for (HighsInt check_index = 0; check_index < num_info; check_index++) {
        if (check_index == index) continue;
        InfoRecordInt& check_info =
            ((InfoRecordInt*)info_records[check_index])[0];
        if (check_info.type == HighsInfoType::kInt) {
          if (check_info.value == value_pointer) {
            highsLogUser(report_log_options, HighsLogType::kError,
                         "checkInfo: Info %" HIGHSINT_FORMAT
                         " (\"%s\") has the same value "
                         "pointer as info %" HIGHSINT_FORMAT " (\"%s\")\n",
                         index, info.name.c_str(), check_index,
                         check_info.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsInfoType::kDouble) {
      // Check double info
      InfoRecordDouble& info = ((InfoRecordDouble*)info_records[index])[0];
      // Check that there are no other info with the same value pointers
      double* value_pointer = info.value;
      for (HighsInt check_index = 0; check_index < num_info; check_index++) {
        if (check_index == index) continue;
        InfoRecordDouble& check_info =
            ((InfoRecordDouble*)info_records[check_index])[0];
        if (check_info.type == HighsInfoType::kDouble) {
          if (check_info.value == value_pointer) {
            highsLogUser(report_log_options, HighsLogType::kError,
                         "checkInfo: Info %" HIGHSINT_FORMAT
                         " (\"%s\") has the same value "
                         "pointer as info %" HIGHSINT_FORMAT " (\"%s\")\n",
                         index, info.name.c_str(), check_index,
                         check_info.name.c_str());
            error_found = true;
          }
        }
      }
    }
  }
  if (error_found) return InfoStatus::kIllegalValue;
  highsLogUser(report_log_options, HighsLogType::kInfo,
               "checkInfo: Info are OK\n");
  return InfoStatus::kOk;
}

#ifndef HIGHSINT64
InfoStatus getLocalInfoValue(const HighsLogOptions& report_log_options,
                             const std::string& name, const bool valid,
                             const std::vector<InfoRecord*>& info_records,
                             int64_t& value) {
  HighsInt index;
  InfoStatus status =
      getInfoIndex(report_log_options, name, info_records, index);
  if (status != InfoStatus::kOk) return status;
  if (!valid) return InfoStatus::kUnavailable;
  HighsInfoType type = info_records[index]->type;
  if (type != HighsInfoType::kInt64) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "getInfoValue: Info \"%s\" requires value of type %s, not int64_t\n",
        name.c_str(), infoEntryTypeToString(type).c_str());
    return InfoStatus::kIllegalValue;
  }
  InfoRecordInt64 info = ((InfoRecordInt64*)info_records[index])[0];
  value = *info.value;
  return InfoStatus::kOk;
}
#endif

InfoStatus getLocalInfoValue(const HighsLogOptions& report_log_options,
                             const std::string& name, const bool valid,
                             const std::vector<InfoRecord*>& info_records,
                             HighsInt& value) {
  HighsInt index;
  InfoStatus status =
      getInfoIndex(report_log_options, name, info_records, index);
  if (status != InfoStatus::kOk) return status;
  if (!valid) return InfoStatus::kUnavailable;
  HighsInfoType type = info_records[index]->type;
  bool type_ok = type == HighsInfoType::kInt;
  // When HIGHSINT64 is "on", value is int64_t and this method is used
  // get HighsInfo values of type HighsInfoType::kInt64
#ifdef HIGHSINT64
  type_ok = type_ok || type == HighsInfoType::kInt64;
#endif
  if (!type_ok) {
    std::string illegal_type = "HighsInt";
#ifdef HIGHSINT64
    illegal_type += " or int64_t";
#endif
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "getInfoValue: Info \"%s\" requires value of type %s, not %s\n",
        name.c_str(), infoEntryTypeToString(type).c_str(),
        illegal_type.c_str());
    return InfoStatus::kIllegalValue;
  }
  if (type == HighsInfoType::kInt) {
    InfoRecordInt info = ((InfoRecordInt*)info_records[index])[0];
    value = *info.value;
  } else {
    assert(type == HighsInfoType::kInt64);
    InfoRecordInt64 info = ((InfoRecordInt64*)info_records[index])[0];
    value = *info.value;
  }
  return InfoStatus::kOk;
}

InfoStatus getLocalInfoValue(const HighsLogOptions& report_log_options,
                             const std::string& name, const bool valid,
                             const std::vector<InfoRecord*>& info_records,
                             double& value) {
  HighsInt index;
  InfoStatus status =
      getInfoIndex(report_log_options, name, info_records, index);
  if (status != InfoStatus::kOk) return status;
  if (!valid) return InfoStatus::kUnavailable;
  HighsInfoType type = info_records[index]->type;
  if (type != HighsInfoType::kDouble) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "getInfoValue: Info \"%s\" requires value of type %s, not double\n",
        name.c_str(), infoEntryTypeToString(type).c_str());
    return InfoStatus::kIllegalValue;
  }
  InfoRecordDouble info = ((InfoRecordDouble*)info_records[index])[0];
  value = *info.value;
  return InfoStatus::kOk;
}

InfoStatus getLocalInfoType(const HighsLogOptions& report_log_options,
                            const std::string& name,
                            const std::vector<InfoRecord*>& info_records,
                            HighsInfoType& type) {
  HighsInt index;
  InfoStatus status =
      getInfoIndex(report_log_options, name, info_records, index);
  if (status != InfoStatus::kOk) return status;
  type = info_records[index]->type;
  return InfoStatus::kOk;
}

HighsStatus writeInfoToFile(FILE* file, const bool valid, const HighsInfo& info,
                            const HighsFileType file_type) {
  return writeInfoToFile(file, valid, info.records, file_type);
}

HighsStatus writeInfoToFile(FILE* file, const bool valid,
                            const std::vector<InfoRecord*>& info_records,
                            const HighsFileType file_type) {
  const bool documentation_file = file_type == HighsFileType::kMd;
  if (!documentation_file && !valid) return HighsStatus::kWarning;
  if (documentation_file || valid) reportInfo(file, info_records, file_type);
  return HighsStatus::kOk;
}

void reportInfo(FILE* file, const std::vector<InfoRecord*>& info_records,
                const HighsFileType file_type) {
  HighsInt num_info = info_records.size();
  for (HighsInt index = 0; index < num_info; index++) {
    HighsInfoType type = info_records[index]->type;
    if (type == HighsInfoType::kInt64) {
      reportInfo(file, ((InfoRecordInt64*)info_records[index])[0], file_type);
    } else if (type == HighsInfoType::kInt) {
      reportInfo(file, ((InfoRecordInt*)info_records[index])[0], file_type);
    } else {
      reportInfo(file, ((InfoRecordDouble*)info_records[index])[0], file_type);
    }
  }
}

void reportInfo(FILE* file, const InfoRecordInt64& info,
                const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: long integer\n\n",
            highsInsertMdEscapes(info.name).c_str(),
            highsInsertMdEscapes(info.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: int64_t]\n%s = %" PRId64 "\n",
            info.description.c_str(), info.name.c_str(), *info.value);
  } else {
    fprintf(file, "%-30s = %" PRId64 "\n", info.name.c_str(), *info.value);
  }
}

void reportInfo(FILE* file, const InfoRecordInt& info,
                const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: integer\n\n",
            highsInsertMdEscapes(info.name).c_str(),
            highsInsertMdEscapes(info.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: HighsInt]\n%s = %" HIGHSINT_FORMAT "\n",
            info.description.c_str(), info.name.c_str(), *info.value);
  } else {
    fprintf(file, "%-30s = %" HIGHSINT_FORMAT "\n", info.name.c_str(),
            *info.value);
  }
}

void reportInfo(FILE* file, const InfoRecordDouble& info,
                const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: double\n\n",
            highsInsertMdEscapes(info.name).c_str(),
            highsInsertMdEscapes(info.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: double]\n%s = %g\n",
            info.description.c_str(), info.name.c_str(), *info.value);
  } else {
    fprintf(file, "%-30s = %g\n", info.name.c_str(), *info.value);
  }
}
