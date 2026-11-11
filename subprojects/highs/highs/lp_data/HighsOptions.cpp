/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.cpp
 * @brief
 */
#include "lp_data/HighsOptions.h"

#include <algorithm>
#include <cassert>
#include <cctype>

#include "HighsExternalApi.h"
#include "util/stringutil.h"

template <typename RecordType>
static RecordType& getOptionRecord(OptionRecord* record) {
  return *static_cast<RecordType*>(record);
}

void highsOpenLogFile(HighsLogOptions& log_options,
                      std::vector<OptionRecord*>& option_records,
                      const std::string& log_file) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(log_options, "log_file", option_records, index);
  assert(status == OptionStatus::kOk);
  if (log_options.log_stream != nullptr) {
    // Current log file stream is not null, so flush and close it
    fflush(log_options.log_stream);
    fclose(log_options.log_stream);
    // Set the stream to null to give a test whether it has been
    // closed (and avoid trying to close it again which causes an
    // error)
    log_options.log_stream = nullptr;
  }
  assert(!log_options.log_stream);
  // If new log file name is not empty, open it, appending if possible
  //
  // If fopen fails then it returns nullptr, so log_file is open for
  // writing or nullptr
  if (log_file.compare(""))
    log_options.log_stream = fopen(log_file.c_str(), "a");
  OptionRecordString& option =
      getOptionRecord<OptionRecordString>(option_records[index]);
  option.assignvalue(log_file);
}

static std::string optionEntryTypeToString(const HighsOptionType type) {
  if (type == HighsOptionType::kBool) {
    return "bool";
  } else if (type == HighsOptionType::kInt) {
    return "HighsInt";
  } else if (type == HighsOptionType::kDouble) {
    return "double";
  } else {
    return "string";
  }
}

static std::string joinOptionValues(const std::vector<std::string>& values) {
  std::string result;
  for (size_t i = 0; i < values.size(); i++) {
    if (i > 0) result += (i == values.size() - 1) ? " or " : ", ";
    result += "\"" + values[i] + "\"";
  }
  return result;
}

bool optionOffChooseOnOk(const HighsLogOptions& report_log_options,
                         const string& name, const string& value) {
  if (value == kHighsOffString || value == kHighsChooseString ||
      value == kHighsOnString)
    return true;
  highsLogUser(
      report_log_options, HighsLogType::kError,
      "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
      name.c_str(),
      joinOptionValues({kHighsOffString, kHighsChooseString, kHighsOnString})
          .c_str());
  return false;
}

bool optionOffOnOk(const HighsLogOptions& report_log_options,
                   const string& name, const string& value) {
  if (value == kHighsOffString || value == kHighsOnString) return true;
  highsLogUser(report_log_options, HighsLogType::kError,
               "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
               name.c_str(),
               joinOptionValues({kHighsOffString, kHighsOnString}).c_str());
  return false;
}

static bool optionSolverValueOk(const HighsLogOptions& report_log_options,
                                const string& value,
                                const std::vector<std::string>& valid_solvers,
                                const string& option_name,
                                const string& option_description,
                                HighsLogType log_type) {
  for (const auto& solver : valid_solvers)
    if (value == solver) return true;
  if (value == kHipoString &&
      !HighsExternalApi::isAvailable<HighsExtras::hipo>()) {
    HighsExternalApi::logUnavailable<HighsExtras::hipo>(
        report_log_options, HighsLogType::kError,
        "The HiPO solver was requested via the \"%s\" option.",
        option_name.c_str());
    return false;
  }
  highsLogUser(report_log_options, log_type,
               "Value \"%s\" for %s option (\"%s\") is not one of %s\n",
               value.c_str(), option_description.c_str(), option_name.c_str(),
               joinOptionValues(valid_solvers).c_str());
  return false;
}

bool optionSolverOk(const HighsLogOptions& report_log_options,
                    const string& value) {
  std::vector<std::string> valid = {
      kHighsChooseString, kSimplexString, kIpmString,  kIpxString,
      kPdlpString,        kHiPdlpString,  kQpAsmString};
  if (HighsExternalApi::isAvailable<HighsExtras::hipo>())
    valid.insert(valid.begin() + 3, kHipoString);
  return optionSolverValueOk(report_log_options, value, valid, kSolverString,
                             "LP/QP solver", HighsLogType::kError);
}

bool optionMipLpSolverOk(const HighsLogOptions& report_log_options,
                         const string& value) {
  std::vector<std::string> valid = {kHighsChooseString, kSimplexString,
                                    kIpmString, kIpxString};
  if (HighsExternalApi::isAvailable<HighsExtras::hipo>())
    valid.insert(valid.begin() + 3, kHipoString);
  return optionSolverValueOk(report_log_options, value, valid,
                             kMipLpSolverString, "MIP LP solver",
                             HighsLogType::kError);
}

bool optionMipIpmSolverOk(const HighsLogOptions& report_log_options,
                          const string& value) {
  std::vector<std::string> valid = {kHighsChooseString, kIpmString, kIpxString};
  if (HighsExternalApi::isAvailable<HighsExtras::hipo>())
    valid.insert(valid.begin() + 2, kHipoString);
  return optionSolverValueOk(report_log_options, value, valid,
                             kMipIpmSolverString, "MIP IPM solver",
                             HighsLogType::kError);
}

bool optionHipoParallelTypeOk(const HighsLogOptions& report_log_options,
                              const string& value) {
  if (value == kHipoNodeString || value == kHipoTreeString ||
      value == kHipoBothString || value == kHighsChooseString)
    return true;
  highsLogUser(report_log_options, HighsLogType::kError,
               "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
               kHipoParallelString.c_str(),
               joinOptionValues({kHipoTreeString, kHipoNodeString,
                                 kHipoBothString, kHighsChooseString})
                   .c_str());
  return false;
}

bool optionHipoSystemOk(const HighsLogOptions& report_log_options,
                        const string& value) {
  if (value == kHipoNormalEqString || value == kHipoAugmentedString ||
      value == kHighsChooseString)
    return true;
  highsLogUser(report_log_options, HighsLogType::kError,
               "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
               kHipoSystemString.c_str(),
               joinOptionValues({kHipoNormalEqString, kHipoAugmentedString,
                                 kHighsChooseString})
                   .c_str());
  return false;
}

bool optionHipoOrderingOk(const HighsLogOptions& report_log_options,
                          const string& value) {
  if (value == kHipoAmdString || value == kHipoMetisString ||
      value == kHipoRcmString || value == kHighsChooseString)
    return true;
  highsLogUser(report_log_options, HighsLogType::kError,
               "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
               kHipoOrderingString.c_str(),
               joinOptionValues({kHipoAmdString, kHipoMetisString,
                                 kHipoRcmString, kHighsChooseString})
                   .c_str());
  return false;
}

bool optionHipoFactorOk(const HighsLogOptions& report_log_options,
                        const string& value) {
  if (value == kHipoFactorMultifrontal || value == kHipoFactorUplooking ||
      value == kHighsChooseString)
    return true;
  highsLogUser(report_log_options, HighsLogType::kError,
               "Value \"%s\" for %s option is not one of %s\n", value.c_str(),
               kHipoFactorString.c_str(),
               joinOptionValues({kHipoFactorMultifrontal, kHipoFactorUplooking,
                                 kHighsChooseString})
                   .c_str());
  return false;
}

bool boolFromString(std::string value, bool& bool_value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (value == "t" || value == "true" || value == "1" ||
      value == kHighsOnString) {
    bool_value = true;
  } else if (value == "f" || value == "false" || value == "0" ||
             value == kHighsOffString) {
    bool_value = false;
  } else {
    return false;
  }
  return true;
}

OptionStatus getOptionIndex(const HighsLogOptions& report_log_options,
                            const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            HighsInt& index) {
  HighsInt num_options = option_records.size();
  for (index = 0; index < num_options; index++)
    if (option_records[index]->name == name) return OptionStatus::kOk;
  highsLogUser(report_log_options, HighsLogType::kError,
               "getOptionIndex: Option \"%s\" is unknown\n", name.c_str());
  return OptionStatus::kUnknownOption;
}

template <typename RecordType>
static bool hasDuplicateValuePointers(
    const HighsLogOptions& report_log_options,
    const std::vector<OptionRecord*>& option_records, HighsInt index) {
  RecordType& option = getOptionRecord<RecordType>(option_records[index]);
  for (HighsInt check_index = 0;
       check_index < static_cast<HighsInt>(option_records.size());
       check_index++) {
    if (check_index == index) continue;
    if (option_records[check_index]->type != option.type) continue;
    RecordType& check_option =
        getOptionRecord<RecordType>(option_records[check_index]);
    if (check_option.value == option.value) {
      highsLogUser(report_log_options, HighsLogType::kError,
                   "checkOptions: Option %" HIGHSINT_FORMAT
                   " (\"%s\") has the same "
                   "value pointer as option %" HIGHSINT_FORMAT " (\"%s\")\n",
                   index, option.name.c_str(), check_index,
                   check_option.name.c_str());
      return true;
    }
  }
  return false;
}

OptionStatus checkOptions(const HighsLogOptions& report_log_options,
                          const std::vector<OptionRecord*>& option_records) {
  bool error_found = false;
  HighsInt num_options = option_records.size();
  for (HighsInt index = 0; index < num_options; index++) {
    std::string name = option_records[index]->name;
    HighsOptionType type = option_records[index]->type;
    // Check that there are no other options with the same name
    for (HighsInt check_index = 0; check_index < num_options; check_index++) {
      if (check_index == index) continue;
      std::string check_name = option_records[check_index]->name;
      if (check_name == name) {
        highsLogUser(report_log_options, HighsLogType::kError,
                     "checkOptions: Option %" HIGHSINT_FORMAT
                     " (\"%s\") has the same name as "
                     "option %" HIGHSINT_FORMAT " \"%s\"\n",
                     index, name.c_str(), check_index, check_name.c_str());
        error_found = true;
      }
    }
    if (type == HighsOptionType::kBool) {
      // Check bool option
      if (hasDuplicateValuePointers<OptionRecordBool>(report_log_options,
                                                      option_records, index))
        error_found = true;
    } else if (type == HighsOptionType::kInt) {
      // Check HighsInt option
      OptionRecordInt& option =
          getOptionRecord<OptionRecordInt>(option_records[index]);
      if (checkOption(report_log_options, option) != OptionStatus::kOk)
        error_found = true;
      if (hasDuplicateValuePointers<OptionRecordInt>(report_log_options,
                                                     option_records, index))
        error_found = true;
    } else if (type == HighsOptionType::kDouble) {
      // Check double option
      OptionRecordDouble& option =
          getOptionRecord<OptionRecordDouble>(option_records[index]);
      if (checkOption(report_log_options, option) != OptionStatus::kOk)
        error_found = true;
      if (hasDuplicateValuePointers<OptionRecordDouble>(report_log_options,
                                                        option_records, index))
        error_found = true;
    } else if (type == HighsOptionType::kString) {
      // Check string option
      if (hasDuplicateValuePointers<OptionRecordString>(report_log_options,
                                                        option_records, index))
        error_found = true;
    }
  }
  if (error_found) return OptionStatus::kIllegalValue;
  highsLogUser(report_log_options, HighsLogType::kInfo,
               "checkOptions: Options are OK\n");
  return OptionStatus::kOk;
}

OptionStatus checkOption(const HighsLogOptions& report_log_options,
                         const OptionRecordInt& option) {
  if (option.lower_bound > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOption: Option \"%s\" has inconsistent bounds "
                 "[%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT "]\n",
                 option.name.c_str(), option.lower_bound, option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  if (option.default_value < option.lower_bound ||
      option.default_value > option.upper_bound) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "checkOption: Option \"%s\" has default value %" HIGHSINT_FORMAT
        " "
        "inconsistent with bounds [%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
        "]\n",
        option.name.c_str(), option.default_value, option.lower_bound,
        option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  HighsInt value = *option.value;
  if (value < option.lower_bound || value > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOption: Option \"%s\" has value %" HIGHSINT_FORMAT
                 " inconsistent with "
                 "bounds [%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT "]\n",
                 option.name.c_str(), value, option.lower_bound,
                 option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  return OptionStatus::kOk;
}

OptionStatus checkOption(const HighsLogOptions& report_log_options,
                         const OptionRecordDouble& option) {
  if (option.lower_bound > option.upper_bound) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "checkOption: Option \"%s\" has inconsistent bounds [%g, %g]\n",
        option.name.c_str(), option.lower_bound, option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  if (option.default_value < option.lower_bound ||
      option.default_value > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOption: Option \"%s\" has default value %g "
                 "inconsistent with bounds [%g, %g]\n",
                 option.name.c_str(), option.default_value, option.lower_bound,
                 option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  double value = *option.value;
  if (value < option.lower_bound || value > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOption: Option \"%s\" has value %g inconsistent with "
                 "bounds [%g, %g]\n",
                 option.name.c_str(), value, option.lower_bound,
                 option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  return OptionStatus::kOk;
}

OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordInt& option, const HighsInt value) {
  if (value < option.lower_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOptionValue: Value %" HIGHSINT_FORMAT
                 " for option \"%s\" is below "
                 "lower bound of %" HIGHSINT_FORMAT "\n",
                 value, option.name.c_str(), option.lower_bound);
    return OptionStatus::kIllegalValue;
  } else if (value > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOptionValue: Value %" HIGHSINT_FORMAT
                 " for option \"%s\" is above "
                 "upper bound of %" HIGHSINT_FORMAT "\n",
                 value, option.name.c_str(), option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  return OptionStatus::kOk;
}

OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordDouble& option, const double value) {
  if (value < option.lower_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOptionValue: Value %g for option \"%s\" is below "
                 "lower bound of %g\n",
                 value, option.name.c_str(), option.lower_bound);
    return OptionStatus::kIllegalValue;
  } else if (value > option.upper_bound) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "checkOptionValue: Value %g for option \"%s\" is above "
                 "upper bound of %g\n",
                 value, option.name.c_str(), option.upper_bound);
    return OptionStatus::kIllegalValue;
  }
  return OptionStatus::kOk;
}

OptionStatus checkOptionValue(const HighsLogOptions& report_log_options,
                              OptionRecordString& option,
                              const std::string& value) {
  // Setting a string option. For some options only particular values
  // are permitted, so check them
  if (option.name == kPresolveString) {
    if (!optionOffChooseOnOk(report_log_options, option.name, value) &&
        value != "mip")
      return OptionStatus::kIllegalValue;
  } else if (option.name == kSolverString) {
    if (!optionSolverOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kMipLpSolverString) {
    if (!optionMipLpSolverOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kMipIpmSolverString) {
    if (!optionMipIpmSolverOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kParallelString) {
    if (!optionOffChooseOnOk(report_log_options, option.name, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kRunCrossoverString) {
    if (!optionOffChooseOnOk(report_log_options, option.name, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kRangingString) {
    if (!optionOffOnOk(report_log_options, option.name, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kHipoParallelString) {
    if (!optionHipoParallelTypeOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kHipoSystemString) {
    if (!optionHipoSystemOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kHipoOrderingString) {
    if (!optionHipoOrderingOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  } else if (option.name == kHipoFactorString) {
    if (!optionHipoFactorOk(report_log_options, value))
      return OptionStatus::kIllegalValue;
  }
  return OptionStatus::kOk;
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const bool value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, name, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kBool) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "setLocalOptionValue: Option \"%s\" cannot be assigned a bool\n",
        name.c_str());
    return OptionStatus::kIllegalValue;
  }
  return setLocalOptionValue(
      getOptionRecord<OptionRecordBool>(option_records[index]), value);
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const HighsInt value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, name, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kInt) {
    if (type == HighsOptionType::kDouble) {
      // Interpret integer as double
      double use_value = value;
      return setLocalOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordDouble>(option_records[index]),
          use_value);
    }
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "setLocalOptionValue: Option \"%s\" cannot be assigned an int\n",
        name.c_str());
    return OptionStatus::kIllegalValue;
  }
  return setLocalOptionValue(
      report_log_options,
      getOptionRecord<OptionRecordInt>(option_records[index]), value);
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 std::vector<OptionRecord*>& option_records,
                                 const double value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, name, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kDouble) {
    highsLogUser(
        report_log_options, HighsLogType::kError,
        "setLocalOptionValue: Option \"%s\" cannot be assigned a double\n",
        name.c_str());
    return OptionStatus::kIllegalValue;
  }
  return setLocalOptionValue(
      report_log_options,
      getOptionRecord<OptionRecordDouble>(option_records[index]), value);
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 HighsLogOptions& log_options,
                                 std::vector<OptionRecord*>& option_records,
                                 const std::string& value_passed) {
  // Trim any leading and trailing spaces
  std::string value_trim = value_passed;
  trim(value_trim, " ");
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, name, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type == HighsOptionType::kBool) {
    bool value_bool;
    bool return_status = boolFromString(value_trim, value_bool);
    if (!return_status) {
      highsLogUser(report_log_options, HighsLogType::kError,
                   "setLocalOptionValue: Value \"%s\" cannot be interpreted "
                   "as a bool\n",
                   value_trim.c_str());
      return OptionStatus::kIllegalValue;
    }
    return setLocalOptionValue(
        getOptionRecord<OptionRecordBool>(option_records[index]), value_bool);
  } else if (type == HighsOptionType::kInt) {
    // Check that the string only contains legitimate characters
    if (value_trim.find_first_not_of(kLegalIntegerCharacters) !=
        std::string::npos)
      return OptionStatus::kIllegalValue;
    HighsInt value_int;
    int scanned_num_char;
    const char* value_char = value_trim.c_str();
    sscanf(value_char, "%" HIGHSINT_FORMAT "%n", &value_int, &scanned_num_char);
    const int value_num_char = strlen(value_char);
    // Check that all characters in value_char are scanned by the format
    const bool converted_ok = scanned_num_char == value_num_char;
    if (!converted_ok) {
      highsLogDev(report_log_options, HighsLogType::kError,
                  "setLocalOptionValue: Value = \"%s\" converts via sscanf as "
                  "%" HIGHSINT_FORMAT
                  " "
                  "by scanning %" HIGHSINT_FORMAT " of %" HIGHSINT_FORMAT
                  " characters\n",
                  value_trim.c_str(), value_int, scanned_num_char,
                  value_num_char);
      return OptionStatus::kIllegalValue;
    }
    return setLocalOptionValue(
        report_log_options,
        getOptionRecord<OptionRecordInt>(option_records[index]), value_int);
  } else if (type == HighsOptionType::kDouble) {
    // Check that the string only contains legitimate characters -
    // after handling +/- inf
    double value_double = 0;
    tolower(value_trim);
    if (value_trim == "inf" || value_trim == "+inf") {
      value_double = kHighsInf;
    } else if (value_trim == "-inf") {
      value_double = -kHighsInf;
    } else {
      if (value_trim.find_first_not_of(kLegalDoubleCharacters) !=
          std::string::npos)
        return OptionStatus::kIllegalValue;
      value_double = atof(value_trim.c_str());
    }
    return setLocalOptionValue(
        report_log_options,
        getOptionRecord<OptionRecordDouble>(option_records[index]),
        value_double);
  } else {
    // Setting a string option value
    if (!name.compare(kLogFileString)) {
      OptionRecordString& option =
          getOptionRecord<OptionRecordString>(option_records[index]);
      std::string original_log_file = *(option.value);
      if (value_passed.compare(original_log_file)) {
        // Changing the name of the log file
        highsOpenLogFile(log_options, option_records, value_passed);
      }
    }
    if (!name.compare(kModelFileString)) {
      // Don't allow model filename to be changed - it's only an
      // option so that reading of run-time options works
      highsLogUser(report_log_options, HighsLogType::kError,
                   "setLocalOptionValue: model filename cannot be set\n");
      return OptionStatus::kUnknownOption;
    } else {
      return setLocalOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordString>(option_records[index]),
          value_passed);
    }
  }
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 const std::string& name,
                                 HighsLogOptions& log_options,
                                 std::vector<OptionRecord*>& option_records,
                                 const char* value) {
  // Handles values passed as explicit values in quotes
  std::string value_as_string(value);
  return setLocalOptionValue(report_log_options, name, log_options,
                             option_records, value_as_string);
}

OptionStatus setLocalOptionValue(OptionRecordBool& option, const bool value) {
  option.assignvalue(value);
  return OptionStatus::kOk;
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordInt& option,
                                 const HighsInt value) {
  OptionStatus return_status =
      checkOptionValue(report_log_options, option, value);
  if (return_status != OptionStatus::kOk) return return_status;
  option.assignvalue(value);
  return OptionStatus::kOk;
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordDouble& option,
                                 const double value) {
  OptionStatus return_status =
      checkOptionValue(report_log_options, option, value);
  if (return_status != OptionStatus::kOk) return return_status;
  option.assignvalue(value);
  return OptionStatus::kOk;
}

OptionStatus setLocalOptionValue(const HighsLogOptions& report_log_options,
                                 OptionRecordString& option,
                                 const std::string& value) {
  std::string possible_lower_case_value = value;
  // Trim any leading and trailing spaces
  trim(possible_lower_case_value, " ");
  // Possibly convert to lower case - not done for file names
  possibleLowerCaseOptionValue(option.name, possible_lower_case_value);
  OptionStatus return_status =
      checkOptionValue(report_log_options, option, possible_lower_case_value);
  if (return_status != OptionStatus::kOk) return return_status;
  option.assignvalue(possible_lower_case_value);
  return OptionStatus::kOk;
}

void possibleLowerCaseOptionValue(const std::string& name, std::string& value) {
  // Don't convert values for file name options to lower case
  if (name == kModelFileString || name == kReadBasisFileString ||
      name == kWriteBasisFileString || name == kOptionsFileString ||
      name == kWriteSolutionFileString || name == kWriteModelFileString ||
      name == kWritePresolvedModelFileString ||
      name == kWriteIisModelFileString || name == kReadSolutionFileString ||
      name == kLogFileString ||
#ifdef HIGHS_DEBUGSOL
      name == kMipDebugSolutionFileString ||
#endif
      name == kMipImprovingSolutionFileString)
    return;
  // Transform other options to lower case
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

OptionStatus passLocalOptions(const HighsLogOptions& report_log_options,
                              const HighsOptions& from_options,
                              HighsOptions& to_options) {
  // (Attempt to) set option value from the HighsOptions passed in
  OptionStatus return_status;
  //  std::string empty_file = "";
  //  std::string from_log_file = from_options.log_file;
  //  std::string original_to_log_file = to_options.log_file;
  //  FILE* original_to_log_stream =
  //  to_options.log_options.log_stream;
  HighsInt num_options = to_options.records.size();
  // Check all the option values before setting any of them - in case
  // to_options are the main Highs options. Checks are only needed for
  // HighsInt, double and string since bool values can't be illegal
  for (HighsInt index = 0; index < num_options; index++) {
    HighsOptionType type = to_options.records[index]->type;
    if (type == HighsOptionType::kInt) {
      HighsInt value = *(
          getOptionRecord<OptionRecordInt>(from_options.records[index]).value);
      return_status = checkOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordInt>(to_options.records[index]), value);
      if (return_status != OptionStatus::kOk) return return_status;
    } else if (type == HighsOptionType::kDouble) {
      double value =
          *(getOptionRecord<OptionRecordDouble>(from_options.records[index])
                .value);
      return_status = checkOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordDouble>(to_options.records[index]),
          value);
      if (return_status != OptionStatus::kOk) return return_status;
    } else if (type == HighsOptionType::kString) {
      std::string value =
          *(getOptionRecord<OptionRecordString>(from_options.records[index])
                .value);
      return_status = checkOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordString>(to_options.records[index]),
          value);
      if (return_status != OptionStatus::kOk) return return_status;
    }
  }
  // Checked from_options and found it to be OK, so set all the values
  for (HighsInt index = 0; index < num_options; index++) {
    HighsOptionType type = to_options.records[index]->type;
    if (type == HighsOptionType::kBool) {
      bool value = *(
          getOptionRecord<OptionRecordBool>(from_options.records[index]).value);
      return_status = setLocalOptionValue(
          getOptionRecord<OptionRecordBool>(to_options.records[index]), value);
      if (return_status != OptionStatus::kOk) return return_status;
    } else if (type == HighsOptionType::kInt) {
      HighsInt value = *(
          getOptionRecord<OptionRecordInt>(from_options.records[index]).value);
      return_status = setLocalOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordInt>(to_options.records[index]), value);
      if (return_status != OptionStatus::kOk) return return_status;
    } else if (type == HighsOptionType::kDouble) {
      double value =
          *(getOptionRecord<OptionRecordDouble>(from_options.records[index])
                .value);
      return_status = setLocalOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordDouble>(to_options.records[index]),
          value);
      if (return_status != OptionStatus::kOk) return return_status;
    } else {
      std::string value =
          *(getOptionRecord<OptionRecordString>(from_options.records[index])
                .value);
      return_status = setLocalOptionValue(
          report_log_options,
          getOptionRecord<OptionRecordString>(to_options.records[index]),
          value);
      if (return_status != OptionStatus::kOk) return return_status;
    }
  }
  /*
  if (from_log_file.compare(original_to_log_file)) {
    // The log file name has changed
    if (from_options.log_options.log_stream &&
        !original_to_log_file.compare(empty_file)) {
      // The stream corresponding to from_log_file is non-null and the
      // original log file name was empty, so to_options inherits the
      // stream, but associated with the (necessarily) non-empty name
      // from_log_file
      //
      // This ensures that the stream to Highs.log opened in
      // RunHighs.cpp is retained unless the log file name is changed.
      assert(from_log_file.compare(empty_file));
      assert(!original_to_log_stream);
      to_options.log_options.log_stream =
          from_options.log_options.log_stream;
    } else {
      highsOpenLogFile(to_options, to_options.log_file);
    }
  }
  */
  return OptionStatus::kOk;
}

OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& option,
    const std::vector<OptionRecord*>& option_records, bool* current_value,
    bool* default_value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, option, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kBool) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getLocalOptionValue: Option \"%s\" requires value of type "
                 "%s, not bool\n",
                 option.c_str(), optionEntryTypeToString(type).c_str());
    return OptionStatus::kIllegalValue;
  }
  OptionRecordBool& option_record =
      getOptionRecord<OptionRecordBool>(option_records[index]);
  if (current_value) *current_value = *(option_record.value);
  if (default_value) *default_value = option_record.default_value;
  return OptionStatus::kOk;
}

OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& option,
    const std::vector<OptionRecord*>& option_records, HighsInt* current_value,
    HighsInt* min_value, HighsInt* max_value, HighsInt* default_value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, option, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kInt) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getLocalOptionValue: Option \"%s\" requires value of type "
                 "%s, not HighsInt\n",
                 option.c_str(), optionEntryTypeToString(type).c_str());
    return OptionStatus::kIllegalValue;
  }
  OptionRecordInt& option_record =
      getOptionRecord<OptionRecordInt>(option_records[index]);
  if (current_value) *current_value = *(option_record.value);
  if (min_value) *min_value = option_record.lower_bound;
  if (max_value) *max_value = option_record.upper_bound;
  if (default_value) *default_value = option_record.default_value;
  return OptionStatus::kOk;
}

OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& option,
    const std::vector<OptionRecord*>& option_records, double* current_value,
    double* min_value, double* max_value, double* default_value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, option, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kDouble) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getLocalOptionValue: Option \"%s\" requires value of type "
                 "%s, not double\n",
                 option.c_str(), optionEntryTypeToString(type).c_str());
    return OptionStatus::kIllegalValue;
  }
  OptionRecordDouble& option_record =
      getOptionRecord<OptionRecordDouble>(option_records[index]);
  if (current_value) *current_value = *(option_record.value);
  if (min_value) *min_value = option_record.lower_bound;
  if (max_value) *max_value = option_record.upper_bound;
  if (default_value) *default_value = option_record.default_value;
  return OptionStatus::kOk;
}

OptionStatus getLocalOptionValues(
    const HighsLogOptions& report_log_options, const std::string& option,
    const std::vector<OptionRecord*>& option_records,
    std::string* current_value, std::string* default_value) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, option, option_records, index);
  if (status != OptionStatus::kOk) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::kString) {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "getLocalOptionValue: Option \"%s\" requires value of type "
                 "%s, not string\n",
                 option.c_str(), optionEntryTypeToString(type).c_str());
    return OptionStatus::kIllegalValue;
  }
  OptionRecordString& option_record =
      getOptionRecord<OptionRecordString>(option_records[index]);
  if (current_value) *current_value = *(option_record.value);
  if (default_value) *default_value = option_record.default_value;
  return OptionStatus::kOk;
}

OptionStatus getLocalOptionType(
    const HighsLogOptions& report_log_options, const std::string& option,
    const std::vector<OptionRecord*>& option_records, HighsOptionType* type) {
  HighsInt index;
  OptionStatus status =
      getOptionIndex(report_log_options, option, option_records, index);
  if (status != OptionStatus::kOk) return status;
  if (type) *type = option_records[index]->type;
  return OptionStatus::kOk;
}

void resetLocalOptions(std::vector<OptionRecord*>& option_records) {
  HighsInt num_options = option_records.size();
  for (HighsInt index = 0; index < num_options; index++) {
    HighsOptionType type = option_records[index]->type;
    if (type == HighsOptionType::kBool) {
      OptionRecordBool& option =
          getOptionRecord<OptionRecordBool>(option_records[index]);
      *(option.value) = option.default_value;
    } else if (type == HighsOptionType::kInt) {
      OptionRecordInt& option =
          getOptionRecord<OptionRecordInt>(option_records[index]);
      *(option.value) = option.default_value;
    } else if (type == HighsOptionType::kDouble) {
      OptionRecordDouble& option =
          getOptionRecord<OptionRecordDouble>(option_records[index]);
      *(option.value) = option.default_value;
    } else {
      OptionRecordString& option =
          getOptionRecord<OptionRecordString>(option_records[index]);
      *(option.value) = option.default_value;
    }
  }
}

HighsStatus writeOptionsToFile(FILE* file, const HighsLogOptions& log_options,
                               const std::vector<OptionRecord*>& option_records,
                               const bool report_only_deviations,
                               const HighsFileType file_type) {
  reportOptions(file, log_options, option_records, report_only_deviations,
                file_type);
  return HighsStatus::kOk;
}

void reportOptions(FILE* file, const HighsLogOptions& log_options,
                   const std::vector<OptionRecord*>& option_records,
                   const bool report_only_deviations,
                   const HighsFileType file_type) {
  HighsInt num_options = option_records.size();
  const bool not_md_or_full =
      file_type != HighsFileType::kMd && file_type != HighsFileType::kFull;
  if (file_type == HighsFileType::kMd)
    fprintf(file, "# [List of options](@id option-definitions)\n\n");
  for (HighsInt index = 0; index < num_options; index++) {
    HighsOptionType type = option_records[index]->type;
    if (not_md_or_full) {
      if (option_records[index]->name == kLogFileString) {
        // Default HiGHS log file name is "" so that deviations from it
        // trigger opening the log file. However, it's unnecessary to
        // report the deviation to kLogFileString, which is the default
        // non-empty log file name in HighsRun.cpp
        if (*getOptionRecord<OptionRecordString>(option_records[index]).value ==
            kHighsRunLogFile)
          continue;
      }
    }
    // Only report non-advanced options
    if (option_records[index]->advanced) {
      // Possibly skip the advanced options when creating Md file
      if (!kAdvancedInDocumentation) continue;
    }
    if (type == HighsOptionType::kBool) {
      reportOption(file, log_options,
                   getOptionRecord<OptionRecordBool>(option_records[index]),
                   report_only_deviations, file_type);
    } else if (type == HighsOptionType::kInt) {
      reportOption(file, log_options,
                   getOptionRecord<OptionRecordInt>(option_records[index]),
                   report_only_deviations, file_type);
    } else if (type == HighsOptionType::kDouble) {
      reportOption(file, log_options,
                   getOptionRecord<OptionRecordDouble>(option_records[index]),
                   report_only_deviations, file_type);
    } else {
      reportOption(file, log_options,
                   getOptionRecord<OptionRecordString>(option_records[index]),
                   report_only_deviations, file_type);
    }
  }
}

void reportOption(FILE* file, const HighsLogOptions& log_options,
                  const OptionRecordBool& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type) {
  if (!report_only_deviations || option.default_value != *option.value) {
    if (file_type == HighsFileType::kMd) {
      fprintf(file,
              "## [%s](@id option-%s)\n- %s\n- Type: boolean\n- Default: "
              "\"%s\"\n\n",
              highsInsertMdEscapes(option.name).c_str(),
              highsInsertMdId(option.name).c_str(),
              highsInsertMdEscapes(option.description).c_str(),
              highsBoolToString(option.default_value).c_str());
    } else if (file_type == HighsFileType::kFull) {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(
          file,
          "# [type: bool, advanced: %s, range: {false, true}, default: %s]\n",
          highsBoolToString(option.advanced).c_str(),
          highsBoolToString(option.default_value).c_str());
      fprintf(file, "%s = %s\n", option.name.c_str(),
              highsBoolToString(*option.value).c_str());
    } else {
      std::string line =
          highsFormatToString("Set option %s to %s\n", option.name.c_str(),
                              highsBoolToString(*option.value).c_str());
      if (file == stdout) {
        highsLogUser(log_options, HighsLogType::kInfo, "%s", line.c_str());
      } else {
        fprintf(file, "%s", line.c_str());
      }
    }
  }
}

void reportOption(FILE* file, const HighsLogOptions& log_options,
                  const OptionRecordInt& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type) {
  if (!report_only_deviations || option.default_value != *option.value) {
    if (file_type == HighsFileType::kMd) {
      fprintf(file,
              "## [%s](@id option-%s)\n- %s\n- Type: integer\n- Range: {%d, "
              "%d}\n- Default: %d\n\n",
              highsInsertMdEscapes(option.name).c_str(),
              highsInsertMdId(option.name).c_str(),
              highsInsertMdEscapes(option.description).c_str(),
              int(option.lower_bound), int(option.upper_bound),
              int(option.default_value));
    } else if (file_type == HighsFileType::kFull) {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file,
              "# [type: integer, advanced: %s, range: {%" HIGHSINT_FORMAT
              ", %d}, default: %d]\n",
              highsBoolToString(option.advanced).c_str(), option.lower_bound,
              int(option.upper_bound), int(option.default_value));
      fprintf(file, "%s = %d\n", option.name.c_str(), int(*option.value));
    } else {
      std::string line = highsFormatToString(
          "Set option %s to %d\n", option.name.c_str(), int(*option.value));
      if (file == stdout) {
        highsLogUser(log_options, HighsLogType::kInfo, "%s", line.c_str());
      } else {
        fprintf(file, "%s", line.c_str());
      }
    }
  }
}

void reportOption(FILE* file, const HighsLogOptions& log_options,
                  const OptionRecordDouble& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type) {
  if (!report_only_deviations || option.default_value != *option.value) {
    if (file_type == HighsFileType::kMd) {
      fprintf(file,
              "## [%s](@id option-%s)\n- %s\n- Type: double\n- Range: [%g, "
              "%g]\n- Default: %g\n\n",
              highsInsertMdEscapes(option.name).c_str(),
              highsInsertMdId(option.name).c_str(),
              highsInsertMdEscapes(option.description).c_str(),
              option.lower_bound, option.upper_bound, option.default_value);
    } else if (file_type == HighsFileType::kFull) {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file,
              "# [type: double, advanced: %s, range: [%g, %g], default: %g]\n",
              highsBoolToString(option.advanced).c_str(), option.lower_bound,
              option.upper_bound, option.default_value);
      fprintf(file, "%s = %g\n", option.name.c_str(), *option.value);
    } else {
      std::string line = highsFormatToString(
          "Set option %s to %g\n", option.name.c_str(), *option.value);
      if (file == stdout) {
        highsLogUser(log_options, HighsLogType::kInfo, "%s", line.c_str());
      } else {
        fprintf(file, "%s", line.c_str());
      }
    }
  }
}

void reportOption(FILE* file, const HighsLogOptions& log_options,
                  const OptionRecordString& option,
                  const bool report_only_deviations,
                  const HighsFileType file_type) {
  // Don't report for the options file if writing to an options file
  // Don't report options that can only be passed via the command line
  if (option.name == kOptionsFileString) return;
  // ToDo: are there others?

  if (!report_only_deviations || option.default_value != *option.value) {
    if (file_type == HighsFileType::kMd) {
      fprintf(file,
              "## [%s](@id option-%s)\n- %s\n- Type: string\n- Default: "
              "\"%s\"\n\n",
              highsInsertMdEscapes(option.name).c_str(),
              highsInsertMdId(option.name).c_str(),
              highsInsertMdEscapes(option.description).c_str(),
              option.default_value.c_str());
    } else if (file_type == HighsFileType::kFull) {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file, "# [type: string, advanced: %s, default: \"%s\"]\n",
              highsBoolToString(option.advanced).c_str(),
              option.default_value.c_str());
      fprintf(file, "%s = %s\n", option.name.c_str(), (*option.value).c_str());
    } else {
      std::string line =
          highsFormatToString("Set option %s to \"%s\"\n", option.name.c_str(),
                              (*option.value).c_str());
      if (file == stdout) {
        highsLogUser(log_options, HighsLogType::kInfo, "%s", line.c_str());
      } else {
        fprintf(file, "%s", line.c_str());
      }
    }
  }
}

void HighsOptions::setLogOptions() {
  this->log_options.output_flag = &this->output_flag;
  this->log_options.log_to_console = &this->log_to_console;
  this->log_options.log_dev_level = &this->log_dev_level;
}

void warnSolverInvalid(const HighsOptions& options,
                       const std::string& problem_type) {
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "Solver \"%s\" is not available for %s. Option \"%s\" ignored\n",
               options.solver.c_str(), problem_type.c_str(),
               kSolverString.c_str());
}
bool solverValidForLp(const std::string& solver) {
  return solver == kHighsChooseString || solver == kSimplexString ||
         solver == kIpmString || solver == kIpxString ||
         solver == kHipoString || solver == kPdlpString ||
         solver == kHiPdlpString;
}
bool solverValidForMip(const std::string& solver) {
  return solver == kHighsChooseString;
}
bool solverValidForQp(const std::string& solver) {
  return solver == kHighsChooseString || solver == kQpAsmString ||
         solver == kIpmString || solver == kHipoString;
}

bool useIpm(const std::string& solver) {
  return solver == kIpmString || solver == kHipoString || solver == kIpxString;
}

bool usePdlp(const std::string& solver) {
  return solver == kPdlpString || solver == kHiPdlpString;
}

bool mayRequireDualPostsolve(const HighsOptions* options) {
  if (usePdlp(options->solver)) return true;
  if (useIpm(options->solver)) {
    if (options->run_crossover == kHighsOffString) return true;
    if (options->run_crossover == kHighsChooseString) return true;
    // Currently if run_crossover is on, and crossover fails, then
    // simplex is started, so basis postsolve will be used. If this
    // policy changes and a primal-dual solution is found, then the
    // following will be required
    //
    //    if (options->run_crossover == kHighsOnString) return true;
  }
  return false;
}

bool mayRequireBasisPostsolve(const HighsOptions* options) {
  if (usePdlp(options->solver)) return false;
  if (useIpm(options->solver) && options->run_crossover == kHighsOffString)
    return false;
  return true;
}
