/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInfo.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_INFO_H_
#define LP_DATA_HIGHS_INFO_H_

#include <cstring>  // For strchr
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"

class HighsOptions;

enum class InfoStatus { OK = 0, NO_FILE, UNKNOWN_INFO, ILLEGAL_VALUE };

class InfoRecord {
 public:
  HighsInfoType type;
  std::string name;
  std::string description;
  bool advanced;

  InfoRecord(HighsInfoType Xtype, std::string Xname, std::string Xdescription,
             bool Xadvanced) {
    this->type = Xtype;
    this->name = Xname;
    this->description = Xdescription;
    this->advanced = Xadvanced;
  }

  virtual ~InfoRecord() {}
};

class InfoRecordInt : public InfoRecord {
 public:
  int* value;
  int default_value;
  InfoRecordInt(std::string Xname, std::string Xdescription, bool Xadvanced,
                int* Xvalue_pointer, int Xdefault_value)
      : InfoRecord(HighsInfoType::INT, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  virtual ~InfoRecordInt() {}
};

class InfoRecordDouble : public InfoRecord {
 public:
  double* value;
  double default_value;
  InfoRecordDouble(std::string Xname, std::string Xdescription, bool Xadvanced,
                   double* Xvalue_pointer, double Xdefault_value)
      : InfoRecord(HighsInfoType::DOUBLE, Xname, Xdescription, Xadvanced) {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

  virtual ~InfoRecordDouble() {}
};

InfoStatus getInfoIndex(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        int& index);

InfoStatus checkInfo(const HighsOptions& options,
                     const std::vector<InfoRecord*>& info_records);
InfoStatus checkInfo(const InfoRecordInt& info);
InfoStatus checkInfo(const InfoRecordDouble& info);

InfoStatus getInfoValue(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        int& value);
InfoStatus getInfoValue(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        double& value);

HighsStatus writeInfoToFile(FILE* file,
                            const std::vector<InfoRecord*>& info_records,
                            const bool html = false);
void reportInfo(FILE* file, const std::vector<InfoRecord*>& info_records,
                const bool html = false);
void reportInfo(FILE* file, const InfoRecordInt& info, const bool html = false);
void reportInfo(FILE* file, const InfoRecordDouble& info,
                const bool html = false);

// For now, but later change so HiGHS properties are string based so that new
// info (for debug and testing too) can be added easily. The info below
// are just what has been used to parse info from argv.
// todo: when creating the new info don't forget underscores for class
// variables but no underscores for struct
struct HighsInfoStruct {
  int simplex_iteration_count;
  int ipm_iteration_count;
  int crossover_iteration_count;
  int primal_status;
  int dual_status;
  double objective_function_value;
  int num_primal_infeasibilities;
  double max_primal_infeasibility;
  double sum_primal_infeasibilities;
  int num_dual_infeasibilities;
  double max_dual_infeasibility;
  double sum_dual_infeasibilities;
};

class HighsInfo : public HighsInfoStruct {
 public:
  HighsInfo() { initRecords(); }

  HighsInfo(const HighsInfo& info) {
    initRecords();
    HighsInfoStruct::operator=(info);
  }

  HighsInfo(HighsInfo&& info) {
    records = std::move(info.records);
    HighsInfoStruct::operator=(std::move(info));
  }

  const HighsInfo& operator=(const HighsInfo& other) {
    if (&other != this) {
      if ((int)records.size() == 0) initRecords();
      HighsInfoStruct::operator=(other);
    }
    return *this;
  }

  const HighsInfo& operator=(HighsInfo&& other) {
    if (&other != this) {
      if ((int)records.size() == 0) initRecords();
      HighsInfoStruct::operator=(other);
    }
    return *this;
  }

  virtual ~HighsInfo() {
    if (records.size() > 0) deleteRecords();
  }

  void clear();

 private:
  void deleteRecords() {
    for (unsigned int i = 0; i < records.size(); i++) delete records[i];
  }

  void initRecords() {
    InfoRecordInt* record_int;
    InfoRecordDouble* record_double;
    bool advanced;
    advanced = false;

    record_int = new InfoRecordInt("simplex_iteration_count",
                                   "Iteration count for simplex solver",
                                   advanced, &simplex_iteration_count, 0);
    records.push_back(record_int);

    record_int = new InfoRecordInt("ipm_iteration_count",
                                   "Iteration count for IPM solver", advanced,
                                   &ipm_iteration_count, 0);
    records.push_back(record_int);

    record_int = new InfoRecordInt("crossover_iteration_count",
                                   "Iteration count for crossover", advanced,
                                   &crossover_iteration_count, 0);
    records.push_back(record_int);

    record_int = new InfoRecordInt(
        "primal_status",
        "Primal status of the model: -1 => Not set; 0 => No solution; 1 => "
        "Unknown; 2 => Infeasible point; 3 => Feasible point",
        advanced, &primal_status, (int)PrimalDualStatus::STATUS_NOTSET);
    records.push_back(record_int);

    record_int = new InfoRecordInt(
        "dual_status",
        "Dual status of the model: -1 => Not set; 0 => No solution; 1 => "
        "Unknown; 2 => Infeasible point; 3 => Feasible point",
        advanced, &dual_status, (int)PrimalDualStatus::STATUS_NOTSET);
    records.push_back(record_int);

    record_double = new InfoRecordDouble("objective_function_value",
                                         "Objective function value", advanced,
                                         &objective_function_value, 0);
    records.push_back(record_double);

    record_int = new InfoRecordInt("num_primal_infeasibilities",
                                   "Number of primal infeasibilities", advanced,
                                   &num_primal_infeasibilities, -1);
    records.push_back(record_int);

    record_double = new InfoRecordDouble(
        "max_primal_infeasibility", "Maximum primal infeasibility", advanced,
        &max_primal_infeasibility, 0);
    records.push_back(record_double);

    record_double = new InfoRecordDouble(
        "sum_primal_infeasibilities", "Sum of primal infeasibilities", advanced,
        &sum_primal_infeasibilities, 0);
    records.push_back(record_double);

    record_int = new InfoRecordInt("num_dual_infeasibilities",
                                   "Number of dual infeasibilities", advanced,
                                   &num_dual_infeasibilities, -1);
    records.push_back(record_int);

    record_double = new InfoRecordDouble("max_dual_infeasibility",
                                         "Maximum dual infeasibility", advanced,
                                         &max_dual_infeasibility, 0);
    records.push_back(record_double);

    record_double = new InfoRecordDouble(
        "sum_dual_infeasibilities", "Sum of dual infeasibilities", advanced,
        &sum_dual_infeasibilities, 0);
    records.push_back(record_double);
  }

 public:
  std::vector<InfoRecord*> records;
};

#endif
