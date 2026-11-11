/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRunData.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_RUN_DATA_H_
#define LP_DATA_HIGHS_RUN_DATA_H_

#include <cstring>  // For strchr
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"

class HighsOptions;

enum class RunDataStatus {
  kOk = 0,
  kUnknownRunData,
  kIllegalValue,
  kUnavailable
};

class RunDataRecord {
 public:
  HighsRunDataType type;
  std::string name;
  std::string description;
  bool advanced;

  RunDataRecord(HighsRunDataType Xtype, const std::string& Xname,
                const std::string& Xdescription, bool Xadvanced) {
    this->type = Xtype;
    this->name = Xname;
    this->description = Xdescription;
    this->advanced = Xadvanced;
  }

  virtual ~RunDataRecord() {}
};

class RunDataRecordInt64 : public RunDataRecord {
 public:
  int64_t* value;
  int64_t default_value;
  RunDataRecordInt64(const std::string& Xname, const std::string& Xdescription,
                     bool Xadvanced, int64_t* Xvalue_pointer,
                     int64_t Xdefault_value)
      : RunDataRecord(HighsRunDataType::kInt64, Xname, Xdescription,
                      Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  virtual ~RunDataRecordInt64() {}
};

class RunDataRecordInt : public RunDataRecord {
 public:
  HighsInt* value;
  HighsInt default_value;
  RunDataRecordInt(const std::string& Xname, const std::string& Xdescription,
                   bool Xadvanced, HighsInt* Xvalue_pointer,
                   HighsInt Xdefault_value)
      : RunDataRecord(HighsRunDataType::kInt, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  virtual ~RunDataRecordInt() {}
};

class RunDataRecordDouble : public RunDataRecord {
 public:
  double* value;
  double default_value;
  RunDataRecordDouble(const std::string& Xname, const std::string& Xdescription,
                      bool Xadvanced, double* Xvalue_pointer,
                      double Xdefault_value)
      : RunDataRecord(HighsRunDataType::kDouble, Xname, Xdescription,
                      Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  virtual ~RunDataRecordDouble() {}
};

struct HighsRunDataStruct {
  bool valid;
  HighsInt presolved_model_num_col;
  HighsInt presolved_model_num_row;
  HighsInt presolved_model_num_nz;
  HighsInt num_simplex_iterations_after_postsolve;
  double presolve_time;
  double solve_time;
  double postsolve_time;
};

class HighsRunData : public HighsRunDataStruct {
 public:
  HighsRunData() { initRecords(); }

  HighsRunData(const HighsRunData& run_data) {
    initRecords();
    HighsRunDataStruct::operator=(run_data);
  }

  HighsRunData(HighsRunData&& run_data) {
    records = std::move(run_data.records);
    HighsRunDataStruct::operator=(std::move(run_data));
  }

  const HighsRunData& operator=(const HighsRunData& other) {
    if (&other != this) {
      if ((HighsInt)records.size() == 0) initRecords();
      HighsRunDataStruct::operator=(other);
    }
    return *this;
  }

  const HighsRunData& operator=(HighsRunData&& other) {
    if (&other != this) {
      if ((HighsInt)records.size() == 0) initRecords();
      HighsRunDataStruct::operator=(other);
    }
    return *this;
  }

  virtual ~HighsRunData() {
    if (records.size() > 0) deleteRecords();
  }

  void invalidate();
  bool equal(const HighsRunData& run_data_) const;

 private:
  void deleteRecords() {
    for (auto record : records) delete record;
  }

  void initRecords() {
    RunDataRecordInt* record_int;
    RunDataRecordDouble* record_double;
    const bool advanced = false;  // Not used

    record_int = new RunDataRecordInt("presolved_model_num_col",
                                      "Number of columns in presolved model",
                                      advanced, &presolved_model_num_col, 0);
    records.push_back(record_int);

    record_int = new RunDataRecordInt("presolved_model_num_row",
                                      "Number of rows in presolved model",
                                      advanced, &presolved_model_num_row, 0);
    records.push_back(record_int);

    record_int = new RunDataRecordInt("presolved_model_num_nz",
                                      "Number of nonzeros in presolved model",
                                      advanced, &presolved_model_num_nz, 0);
    records.push_back(record_int);

    record_int = new RunDataRecordInt(
        "num_simplex_iterations_after_postsolve",
        "Number of simplex iterations after postsolve", advanced,
        &num_simplex_iterations_after_postsolve, 0);
    records.push_back(record_int);

    record_double = new RunDataRecordDouble("presolve_time", "Presolve time",
                                            advanced, &presolve_time, 0);
    records.push_back(record_double);

    record_double = new RunDataRecordDouble("solve_time", "Solve time",
                                            advanced, &solve_time, 0);
    records.push_back(record_double);

    record_double = new RunDataRecordDouble("postsolve_time", "Postsolve time",
                                            advanced, &postsolve_time, 0);
    records.push_back(record_double);
  }

 public:
  std::vector<RunDataRecord*> records;
};

HighsStatus writeRunDataToFile(
    FILE* file, const bool valid, const HighsRunData& run_data,
    const HighsFileType file_type = HighsFileType::kFull);

RunDataStatus getRunDataIndex(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<RunDataRecord*>& run_data_records, HighsInt& index);

RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    int64_t& value);
RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    HighsInt& value);
RunDataStatus getLocalRunDataValue(
    const HighsLogOptions& report_log_options, const std::string& name,
    const bool valid, const std::vector<RunDataRecord*>& run_data_records,
    double& value);

RunDataStatus getLocalRunDataType(
    const HighsLogOptions& report_log_options, const std::string& name,
    const std::vector<RunDataRecord*>& run_data_records,
    HighsRunDataType& type);

HighsStatus writeRunDataToFile(
    FILE* file, const bool valid,
    const std::vector<RunDataRecord*>& run_data_records,
    const HighsFileType file_type = HighsFileType::kFull);
void reportRunData(FILE* file,
                   const std::vector<RunDataRecord*>& run_data_records,
                   const HighsFileType file_type = HighsFileType::kFull);
void reportRunData(FILE* file, const RunDataRecordInt64& run_data,
                   const HighsFileType file_type = HighsFileType::kFull);
void reportRunData(FILE* file, const RunDataRecordInt& run_data,
                   const HighsFileType file_type = HighsFileType::kFull);
void reportRunData(FILE* file, const RunDataRecordDouble& run_data,
                   const HighsFileType file_type = HighsFileType::kFull);

#endif
