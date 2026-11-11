/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRunData.cpp
 * @brief
 */
#include "lp_data/HighsRunData.h"

#include <cassert>

#include "lp_data/HighsOptions.h"

void HighsRunData::invalidate() {
  valid = false;
  presolved_model_num_col = kHighsIllegalIntMeasure;
  presolved_model_num_row = kHighsIllegalIntMeasure;
  presolved_model_num_nz = kHighsIllegalIntMeasure;
  num_simplex_iterations_after_postsolve = kHighsIllegalIntMeasure;
  presolve_time = kHighsIllegalDoubleMeasure;
  solve_time = kHighsIllegalDoubleMeasure;
  postsolve_time = kHighsIllegalDoubleMeasure;
}

bool HighsRunData::equal(const HighsRunData& run_data_) const {
  if (run_data_.valid != this->valid) return false;

  if (run_data_.presolved_model_num_col != this->presolved_model_num_col)
    return false;
  if (run_data_.presolved_model_num_row != this->presolved_model_num_row)
    return false;
  if (run_data_.presolved_model_num_nz != this->presolved_model_num_nz)
    return false;
  if (run_data_.num_simplex_iterations_after_postsolve !=
      this->num_simplex_iterations_after_postsolve)
    return false;
  if (run_data_.presolve_time != this->presolve_time) return false;
  if (run_data_.solve_time != this->solve_time) return false;
  if (run_data_.postsolve_time != this->postsolve_time) return false;
  return true;
}

static std::string run_dataEntryTypeToString(const HighsRunDataType type) {
  if (type == HighsRunDataType::kInt64) {
    return "int64_t";
  } else if (type == HighsRunDataType::kInt) {
    return "HighsInt";
  } else {
    return "double";
  }
}

RunDataStatus getRunDataIndex(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<RunDataRecord*>& run_data_records, HighsInt& index) {
  HighsInt num_run_data = run_data_records.size();
  for (index = 0; index < num_run_data; index++)
    if (run_data_records[index]->name == name) return RunDataStatus::kOk;
  highsLogUser(report_log_options, HighsLogType::kError,
               "getRunDataIndex: Run data \"%s\" is unknown\n", name.c_str());
  return RunDataStatus::kUnknownRunData;
}

#ifndef HIGHSINT64
RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    int64_t& value) {
  HighsInt index;
  RunDataStatus status =
      getRunDataIndex(report_log_options, name, run_data_records, index);
  if (status != RunDataStatus::kOk) return status;
  if (!valid) return RunDataStatus::kUnavailable;
  HighsRunDataType type = run_data_records[index]->type;
  if (type != HighsRunDataType::kInt64) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getRunDataValue: RunData \"%s\" requires value of type %s, "
                 "not int64_t\n",
                 name.c_str(), run_dataEntryTypeToString(type).c_str());
    return RunDataStatus::kIllegalValue;
  }
  RunDataRecordInt64 run_data =
      ((RunDataRecordInt64*)run_data_records[index])[0];
  value = *run_data.value;
  return RunDataStatus::kOk;
}
#endif

RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    HighsInt& value) {
  HighsInt index;
  RunDataStatus status =
      getRunDataIndex(report_log_options, name, run_data_records, index);
  if (status != RunDataStatus::kOk) return status;
  if (!valid) return RunDataStatus::kUnavailable;
  HighsRunDataType type = run_data_records[index]->type;
  bool type_ok = type == HighsRunDataType::kInt;
  // When HIGHSINT64 is "on", value is int64_t and this method is used
  // get HighsRunData values of type HighsRunDataType::kInt64
#ifdef HIGHSINT64
  type_ok = type_ok || type == HighsRunDataType::kInt64;
#endif
  if (!type_ok) {
    std::string illegal_type = "HighsInt";
#ifdef HIGHSINT64
    illegal_type += " or int64_t";
#endif
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "getRunDataValue: RunData \"%s\" requires value of type %s, not %s\n",
        name.c_str(), run_dataEntryTypeToString(type).c_str(),
        illegal_type.c_str());
    return RunDataStatus::kIllegalValue;
  }
  if (type == HighsRunDataType::kInt) {
    RunDataRecordInt run_data = ((RunDataRecordInt*)run_data_records[index])[0];
    value = *run_data.value;
  } else {
    assert(type == HighsRunDataType::kInt64);
    RunDataRecordInt64 run_data =
        ((RunDataRecordInt64*)run_data_records[index])[0];
    value = *run_data.value;
  }
  return RunDataStatus::kOk;
}

RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    double& value) {
  HighsInt index;
  RunDataStatus status =
      getRunDataIndex(report_log_options, name, run_data_records, index);
  if (status != RunDataStatus::kOk) return status;
  if (!valid) return RunDataStatus::kUnavailable;
  HighsRunDataType type = run_data_records[index]->type;
  if (type != HighsRunDataType::kDouble) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getRunDataValue: RunData \"%s\" requires value of type %s, "
                 "not double\n",
                 name.c_str(), run_dataEntryTypeToString(type).c_str());
    return RunDataStatus::kIllegalValue;
  }
  RunDataRecordDouble run_data =
      ((RunDataRecordDouble*)run_data_records[index])[0];
  value = *run_data.value;
  return RunDataStatus::kOk;
}

RunDataStatus getLocalRunDataType(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<RunDataRecord*>& run_data_records,
    HighsRunDataType& type) {
  HighsInt index;
  RunDataStatus status =
      getRunDataIndex(report_log_options, name, run_data_records, index);
  if (status != RunDataStatus::kOk) return status;
  type = run_data_records[index]->type;
  return RunDataStatus::kOk;
}

HighsStatus writeRunDataToFile(FILE* file, const bool valid,
                               const HighsRunData& run_data,
                               const HighsFileType file_type) {
  return writeRunDataToFile(file, valid, run_data.records, file_type);
}

HighsStatus writeRunDataToFile(
    FILE* file, const bool valid,
    const std::vector<RunDataRecord*>& run_data_records,
    const HighsFileType file_type) {
  const bool documentation_file = file_type == HighsFileType::kMd;
  if (!documentation_file && !valid) return HighsStatus::kWarning;
  if (documentation_file || valid)
    reportRunData(file, run_data_records, file_type);
  return HighsStatus::kOk;
}

void reportRunData(FILE* file,
                   const std::vector<RunDataRecord*>& run_data_records,
                   const HighsFileType file_type) {
  HighsInt num_run_data = run_data_records.size();
  for (HighsInt index = 0; index < num_run_data; index++) {
    HighsRunDataType type = run_data_records[index]->type;
    if (type == HighsRunDataType::kInt64) {
      reportRunData(file, ((RunDataRecordInt64*)run_data_records[index])[0],
                    file_type);
    } else if (type == HighsRunDataType::kInt) {
      reportRunData(file, ((RunDataRecordInt*)run_data_records[index])[0],
                    file_type);
    } else {
      reportRunData(file, ((RunDataRecordDouble*)run_data_records[index])[0],
                    file_type);
    }
  }
}

void reportRunData(FILE* file, const RunDataRecordInt64& run_data,
                   const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: long integer\n\n",
            highsInsertMdEscapes(run_data.name).c_str(),
            highsInsertMdEscapes(run_data.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: int64_t]\n%s = %" PRId64 "\n",
            run_data.description.c_str(), run_data.name.c_str(),
            *run_data.value);
  } else {
    fprintf(file, "%-30s = %" PRId64 "\n", run_data.name.c_str(),
            *run_data.value);
  }
}

void reportRunData(FILE* file, const RunDataRecordInt& run_data,
                   const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: integer\n\n",
            highsInsertMdEscapes(run_data.name).c_str(),
            highsInsertMdEscapes(run_data.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: HighsInt]\n%s = %" HIGHSINT_FORMAT "\n",
            run_data.description.c_str(), run_data.name.c_str(),
            *run_data.value);
  } else {
    fprintf(file, "%-30s = %" HIGHSINT_FORMAT "\n", run_data.name.c_str(),
            *run_data.value);
  }
}

void reportRunData(FILE* file, const RunDataRecordDouble& run_data,
                   const HighsFileType file_type) {
  if (file_type == HighsFileType::kMd) {
    fprintf(file, "## %s\n- %s\n- Type: double\n\n",
            highsInsertMdEscapes(run_data.name).c_str(),
            highsInsertMdEscapes(run_data.description).c_str());
  } else if (file_type == HighsFileType::kFull) {
    fprintf(file, "\n# %s\n# [type: double]\n%s = %g\n",
            run_data.description.c_str(), run_data.name.c_str(),
            *run_data.value);
  } else {
    fprintf(file, "%-30s = %g\n", run_data.name.c_str(), *run_data.value);
  }
}
