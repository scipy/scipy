/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsOptions.h"

inline const char* bool2string(bool b) { return b ? "true" : "false"; }

std::string optionEntryType2string(const HighsOptionType type) {
  if (type == HighsOptionType::BOOL) {
    return "bool";
  } else if (type == HighsOptionType::INT) {
    return "int";
  } else if (type == HighsOptionType::DOUBLE) {
    return "double";
  } else {
    return "string";
  }
}

bool commandLineOffChooseOnOk(FILE* logfile, const string& value) {
  if (value == off_string || value == choose_string || value == on_string)
    return true;
  HighsLogMessage(logfile, HighsMessageType::WARNING,
                  "Value \"%s\" is not one of \"%s\", \"%s\" or \"%s\"\n",
                  value.c_str(), off_string.c_str(), choose_string.c_str(),
                  on_string.c_str());
  return false;
}

bool commandLineSolverOk(FILE* logfile, const string& value) {
  if (value == simplex_string || value == choose_string || value == ipm_string)
    return true;
  HighsLogMessage(logfile, HighsMessageType::WARNING,
                  "Value \"%s\" is not one of \"%s\", \"%s\" or \"%s\"\n",
                  value.c_str(), simplex_string.c_str(), choose_string.c_str(),
                  ipm_string.c_str());
  return false;
}

bool boolFromString(const std::string value, bool& bool_value) {
  if (value == "t" || value == "true" || value == "T" || value == "True" ||
      value == "TRUE") {
    bool_value = true;
  } else if (value == "f" || value == "false" || value == "F" ||
             value == "False" || value == "FALSE") {
    bool_value = false;
  } else {
    return false;
  }
  return true;
}

OptionStatus getOptionIndex(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            int& index) {
  int num_options = option_records.size();
  for (index = 0; index < num_options; index++)
    if (option_records[index]->name == name) return OptionStatus::OK;
  HighsLogMessage(logfile, HighsMessageType::ERROR,
                  "getOptionIndex: Option \"%s\" is unknown", name.c_str());
  return OptionStatus::UNKNOWN_OPTION;
}

OptionStatus checkOptions(FILE* logfile,
                          const std::vector<OptionRecord*>& option_records) {
  bool error_found = false;
  int num_options = option_records.size();
  for (int index = 0; index < num_options; index++) {
    std::string name = option_records[index]->name;
    HighsOptionType type = option_records[index]->type;
    // Check that there are no other options with the same name
    for (int check_index = 0; check_index < num_options; check_index++) {
      if (check_index == index) continue;
      std::string check_name = option_records[check_index]->name;
      if (check_name == name) {
        HighsLogMessage(logfile, HighsMessageType::ERROR,
                        "checkOptions: Option %d (\"%s\") has the same name as "
                        "option %d \"%s\"",
                        index, name.c_str(), check_index, check_name.c_str());
        error_found = true;
      }
    }
    if (type == HighsOptionType::BOOL) {
      // Check bool option
      OptionRecordBool& option = ((OptionRecordBool*)option_records[index])[0];
      // Check that there are no other options with the same value pointers
      bool* value_pointer = option.value;
      for (int check_index = 0; check_index < num_options; check_index++) {
        if (check_index == index) continue;
        OptionRecordBool& check_option =
            ((OptionRecordBool*)option_records[check_index])[0];
        if (check_option.type == HighsOptionType::BOOL) {
          if (check_option.value == value_pointer) {
            HighsLogMessage(logfile, HighsMessageType::ERROR,
                            "checkOptions: Option %d (\"%s\") has the same "
                            "value pointer as option %d (\"%s\")",
                            index, option.name.c_str(), check_index,
                            check_option.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsOptionType::INT) {
      // Check int option
      OptionRecordInt& option = ((OptionRecordInt*)option_records[index])[0];
      if (checkOption(logfile, option) != OptionStatus::OK) error_found = true;
      // Check that there are no other options with the same value pointers
      int* value_pointer = option.value;
      for (int check_index = 0; check_index < num_options; check_index++) {
        if (check_index == index) continue;
        OptionRecordInt& check_option =
            ((OptionRecordInt*)option_records[check_index])[0];
        if (check_option.type == HighsOptionType::INT) {
          if (check_option.value == value_pointer) {
            HighsLogMessage(logfile, HighsMessageType::ERROR,
                            "checkOptions: Option %d (\"%s\") has the same "
                            "value pointer as option %d (\"%s\")",
                            index, option.name.c_str(), check_index,
                            check_option.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsOptionType::DOUBLE) {
      // Check double option
      OptionRecordDouble& option =
          ((OptionRecordDouble*)option_records[index])[0];
      if (checkOption(logfile, option) != OptionStatus::OK) error_found = true;
      // Check that there are no other options with the same value pointers
      double* value_pointer = option.value;
      for (int check_index = 0; check_index < num_options; check_index++) {
        if (check_index == index) continue;
        OptionRecordDouble& check_option =
            ((OptionRecordDouble*)option_records[check_index])[0];
        if (check_option.type == HighsOptionType::DOUBLE) {
          if (check_option.value == value_pointer) {
            HighsLogMessage(logfile, HighsMessageType::ERROR,
                            "checkOptions: Option %d (\"%s\") has the same "
                            "value pointer as option %d (\"%s\")",
                            index, option.name.c_str(), check_index,
                            check_option.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsOptionType::STRING) {
      // Check string option
      OptionRecordString& option =
          ((OptionRecordString*)option_records[index])[0];
      // Check that there are no other options with the same value pointers
      std::string* value_pointer = option.value;
      for (int check_index = 0; check_index < num_options; check_index++) {
        if (check_index == index) continue;
        OptionRecordString& check_option =
            ((OptionRecordString*)option_records[check_index])[0];
        if (check_option.type == HighsOptionType::STRING) {
          if (check_option.value == value_pointer) {
            HighsLogMessage(logfile, HighsMessageType::ERROR,
                            "checkOptions: Option %d (\"%s\") has the same "
                            "value pointer as option %d (\"%s\")",
                            index, option.name.c_str(), check_index,
                            check_option.name.c_str());
            error_found = true;
          }
        }
      }
    }
  }
  if (error_found) return OptionStatus::ILLEGAL_VALUE;
  HighsLogMessage(logfile, HighsMessageType::INFO,
                  "checkOptions: Options are OK");
  return OptionStatus::OK;
}

OptionStatus checkOption(FILE* logfile, const OptionRecordInt& option) {
  if (option.lower_bound > option.upper_bound) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "checkOption: Option \"%s\" has inconsistent bounds [%d, %d]",
        option.name.c_str(), option.lower_bound, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  if (option.default_value < option.lower_bound ||
      option.default_value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "checkOption: Option \"%s\" has default value %d "
                    "inconsistent with bounds [%d, %d]",
                    option.name.c_str(), option.default_value,
                    option.lower_bound, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  int value = *option.value;
  if (value < option.lower_bound || value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "checkOption: Option \"%s\" has value %d inconsistent with "
                    "bounds [%d, %d]",
                    option.name.c_str(), value, option.lower_bound,
                    option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus checkOption(FILE* logfile, const OptionRecordDouble& option) {
  if (option.lower_bound > option.upper_bound) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "checkOption: Option \"%s\" has inconsistent bounds [%g, %g]",
        option.name.c_str(), option.lower_bound, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  if (option.default_value < option.lower_bound ||
      option.default_value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "checkOption: Option \"%s\" has default value %g "
                    "inconsistent with bounds [%g, %g]",
                    option.name.c_str(), option.default_value,
                    option.lower_bound, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  double value = *option.value;
  if (value < option.lower_bound || value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "checkOption: Option \"%s\" has value %g inconsistent with "
                    "bounds [%g, %g]",
                    option.name.c_str(), value, option.lower_bound,
                    option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus checkOptionValue(FILE* logfile, OptionRecordInt& option,
                              const int value) {
  if (value < option.lower_bound) {
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "checkOptionValue: Value %d for option \"%s\" is below "
                    "lower bound of %d",
                    value, option.name.c_str(), option.lower_bound);
    return OptionStatus::ILLEGAL_VALUE;
  } else if (value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "checkOptionValue: Value %d for option \"%s\" is above "
                    "upper bound of %d",
                    value, option.name.c_str(), option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus checkOptionValue(FILE* logfile, OptionRecordDouble& option,
                              const double value) {
  if (value < option.lower_bound) {
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "checkOptionValue: Value %g for option \"%s\" is below "
                    "lower bound of %g",
                    value, option.name.c_str(), option.lower_bound);
    return OptionStatus::ILLEGAL_VALUE;
  } else if (value > option.upper_bound) {
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "checkOptionValue: Value %g for option \"%s\" is above "
                    "upper bound of %g",
                    value, option.name.c_str(), option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus checkOptionValue(FILE* logfile, OptionRecordString& option,
                              const std::string value) {
  // Setting a string option. For some options only particular values
  // are permitted, so check them
  if (option.name == presolve_string) {
    if (!commandLineOffChooseOnOk(logfile, value))
      return OptionStatus::ILLEGAL_VALUE;
  } else if (option.name == solver_string) {
    if (!commandLineSolverOk(logfile, value))
      return OptionStatus::ILLEGAL_VALUE;
  } else if (option.name == parallel_string) {
    if (!commandLineOffChooseOnOk(logfile, value))
      return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const bool value) {
  int index;
  //  printf("setOptionValue: \"%s\" with bool %d\n", name.c_str(), value);
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::BOOL) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "setOptionValue: Option \"%s\" cannot be assigned a bool",
                    name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(((OptionRecordBool*)option_records[index])[0], value);
}

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const int value) {
  int index;
  //  printf("setOptionValue: \"%s\" with int %d\n", name.c_str(), value);
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::INT) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "setOptionValue: Option \"%s\" cannot be assigned an int",
                    name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(logfile, ((OptionRecordInt*)option_records[index])[0],
                        value);
}

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const double value) {
  int index;
  //  printf("setOptionValue: \"%s\" with double %g\n", name.c_str(), value);
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::DOUBLE) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "setOptionValue: Option \"%s\" cannot be assigned a double",
                    name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(logfile,
                        ((OptionRecordDouble*)option_records[index])[0], value);
}

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const std::string value) {
  int index;
  //  printf("setOptionValue: \"%s\" with value string %s\n", name.c_str(),
  //  value.c_str());
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type == HighsOptionType::BOOL) {
    bool value_bool;
    bool return_status = boolFromString(value, value_bool);
    //    printf("boolFromString for \"%s\" returns %d from \"%s\" with status
    //    %d\n", name.c_str(), value_bool, value.c_str(), return_status);
    if (!return_status) {
      HighsLogMessage(
          logfile, HighsMessageType::ERROR,
          "setOptionValue: Value \"%s\" cannot be interpreted as a bool",
          value.c_str());
      return OptionStatus::ILLEGAL_VALUE;
    }
    return setOptionValue(((OptionRecordBool*)option_records[index])[0],
                          value_bool);
  } else if (type == HighsOptionType::INT) {
    int value_int;
    int scanned_num_char;
    const char* value_char = value.c_str();
    sscanf(value_char, "%d%n", &value_int, &scanned_num_char);
    const int value_num_char = strlen(value_char);
    const bool converted_ok = scanned_num_char == value_num_char;
    /*
      value_int = atoi(value_char);
      double value_double = atof(value_char);
      double value_int_double = value_int;
      converted_ok = value_double == value_int_double;
    */
    if (!converted_ok) {
      HighsLogMessage(
          logfile, HighsMessageType::ERROR,
          "setOptionValue: Value = \"%s\" converts via sscanf as %d "
          "by scanning %d of %d characters",
          value.c_str(), value_int, scanned_num_char, value_num_char);
      /*
        HighsLogMessage(logfile, HighsMessageType::ERROR,
                        "setOptionValue: Value = \"%s\" converts via atoi as %d
        " "so is %g as double, but as %g via atof", value.c_str(), value_int,
        value_int_double, value_double);
      */
      return OptionStatus::ILLEGAL_VALUE;
    }
    return setOptionValue(logfile, ((OptionRecordInt*)option_records[index])[0],
                          value_int);
  } else if (type == HighsOptionType::DOUBLE) {
    int value_int = atoi(value.c_str());
    double value_double = atof(value.c_str());
    double value_int_double = value_int;
    if (value_double == value_int_double) {
      HighsLogMessage(logfile, HighsMessageType::INFO,
                      "setOptionValue: Value = \"%s\" converts via atoi as %d "
                      "so is %g as double, and %g via atof\n",
                      value.c_str(), value_int, value_int_double, value_double);
    }
    return setOptionValue(logfile,
                          ((OptionRecordDouble*)option_records[index])[0],
                          atof(value.c_str()));
  } else {
    return setOptionValue(
        logfile, ((OptionRecordString*)option_records[index])[0], value);
  }
}

OptionStatus setOptionValue(FILE* logfile, const std::string& name,
                            std::vector<OptionRecord*>& option_records,
                            const char* value) {
  // Handles values passed as explicit values in quotes
  std::string value_as_string(value);
  return setOptionValue(logfile, name, option_records, value_as_string);
}

OptionStatus setOptionValue(OptionRecordBool& option, const bool value) {
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(FILE* logfile, OptionRecordInt& option,
                            const int value) {
  OptionStatus return_status = checkOptionValue(logfile, option, value);
  if (return_status != OptionStatus::OK) return return_status;
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(FILE* logfile, OptionRecordDouble& option,
                            const double value) {
  OptionStatus return_status = checkOptionValue(logfile, option, value);
  if (return_status != OptionStatus::OK) return return_status;
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(FILE* logfile, OptionRecordString& option,
                            const std::string value) {
  OptionStatus return_status = checkOptionValue(logfile, option, value);
  if (return_status != OptionStatus::OK) return return_status;
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus passOptions(FILE* logfile, const HighsOptions& from_options,
                         HighsOptions& to_options) {
  // (Attempt to) set option value from the HighsOptions passed in
  OptionStatus return_status;
  int num_options = to_options.records.size();
  // Check all the option values before setting any of them - in case
  // to_options are the main Highs options. Checks are only needed for
  // int, double and string since bool values can't be illegal
  for (int index = 0; index < num_options; index++) {
    HighsOptionType type = to_options.records[index]->type;
    if (type == HighsOptionType::INT) {
      int value = *(((OptionRecordInt*)from_options.records[index])[0].value);
      return_status = checkOptionValue(
          logfile, ((OptionRecordInt*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    } else if (type == HighsOptionType::DOUBLE) {
      double value =
          *(((OptionRecordDouble*)from_options.records[index])[0].value);
      return_status = checkOptionValue(
          logfile, ((OptionRecordDouble*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    } else if (type == HighsOptionType::STRING) {
      std::string value =
          *(((OptionRecordString*)from_options.records[index])[0].value);
      return_status = checkOptionValue(
          logfile, ((OptionRecordString*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    }
  }
  // Checked from_options and found it to be OK, so set all the values
  for (int index = 0; index < num_options; index++) {
    HighsOptionType type = to_options.records[index]->type;
    if (type == HighsOptionType::BOOL) {
      bool value = *(((OptionRecordBool*)from_options.records[index])[0].value);
      return_status = setOptionValue(
          ((OptionRecordBool*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    } else if (type == HighsOptionType::INT) {
      int value = *(((OptionRecordInt*)from_options.records[index])[0].value);
      return_status = setOptionValue(
          logfile, ((OptionRecordInt*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    } else if (type == HighsOptionType::DOUBLE) {
      double value =
          *(((OptionRecordDouble*)from_options.records[index])[0].value);
      return_status = setOptionValue(
          logfile, ((OptionRecordDouble*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    } else {
      std::string value =
          *(((OptionRecordString*)from_options.records[index])[0].value);
      return_status = setOptionValue(
          logfile, ((OptionRecordString*)to_options.records[index])[0], value);
      if (return_status != OptionStatus::OK) return return_status;
    }
  }
  return OptionStatus::OK;
}

OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            bool& value) {
  int index;
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::BOOL) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "getOptionValue: Option \"%s\" requires value of type %s, not bool",
        name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordBool option = ((OptionRecordBool*)option_records[index])[0];
  value = *option.value;
  return OptionStatus::OK;
}

OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            int& value) {
  int index;
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::INT) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "getOptionValue: Option \"%s\" requires value of type %s, not int",
        name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordInt option = ((OptionRecordInt*)option_records[index])[0];
  value = *option.value;
  return OptionStatus::OK;
}

OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            double& value) {
  int index;
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::DOUBLE) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "getOptionValue: Option \"%s\" requires value of type %s, not double",
        name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordDouble option = ((OptionRecordDouble*)option_records[index])[0];
  value = *option.value;
  return OptionStatus::OK;
}

OptionStatus getOptionValue(FILE* logfile, const std::string& name,
                            const std::vector<OptionRecord*>& option_records,
                            std::string& value) {
  int index;
  OptionStatus status = getOptionIndex(logfile, name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::STRING) {
    HighsLogMessage(
        logfile, HighsMessageType::ERROR,
        "getOptionValue: Option \"%s\" requires value of type %s, not string",
        name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordString option = ((OptionRecordString*)option_records[index])[0];
  value = *option.value;
  return OptionStatus::OK;
}

void resetOptions(std::vector<OptionRecord*>& option_records) {
  int num_options = option_records.size();
  for (int index = 0; index < num_options; index++) {
    HighsOptionType type = option_records[index]->type;
    if (type == HighsOptionType::BOOL) {
      OptionRecordBool& option = ((OptionRecordBool*)option_records[index])[0];
      *(option.value) = option.default_value;
    } else if (type == HighsOptionType::INT) {
      OptionRecordInt& option = ((OptionRecordInt*)option_records[index])[0];
      *(option.value) = option.default_value;
    } else if (type == HighsOptionType::DOUBLE) {
      OptionRecordDouble& option =
          ((OptionRecordDouble*)option_records[index])[0];
      *(option.value) = option.default_value;
    } else {
      OptionRecordString& option =
          ((OptionRecordString*)option_records[index])[0];
      *(option.value) = option.default_value;
    }
  }
}

HighsStatus writeOptionsToFile(FILE* file,
                               const std::vector<OptionRecord*>& option_records,
                               const bool report_only_non_default_values,
                               const bool html) {
  if (html) {
    fprintf(file, "<!DOCTYPE HTML>\n<html>\n\n<head>\n");
    fprintf(file, "  <title>HiGHS Options</title>\n");
    fprintf(file, "	<meta charset=\"utf-8\" />\n");
    fprintf(file,
            "	<meta name=\"viewport\" content=\"width=device-width, "
            "initial-scale=1, user-scalable=no\" />\n");
    fprintf(file,
            "	<link rel=\"stylesheet\" href=\"assets/css/main.css\" />\n");
    fprintf(file, "</head>\n");
    fprintf(file, "<body style=\"background-color:f5fafa;\"></body>\n\n");
    fprintf(file, "<h3>HiGHS Options</h3>\n\n");
    fprintf(file, "<ul>\n");
  }
  reportOptions(file, option_records, report_only_non_default_values, html);
  if (html) {
    fprintf(file, "</ul>\n");
    fprintf(file, "</body>\n\n</html>\n");
  }
  return HighsStatus::OK;
}

void reportOptions(FILE* file, const std::vector<OptionRecord*>& option_records,
                   const bool report_only_non_default_values, const bool html) {
  int num_options = option_records.size();
  for (int index = 0; index < num_options; index++) {
    HighsOptionType type = option_records[index]->type;
    //    fprintf(file, "\n# Option %1d\n", index);
    // Skip the advanced options when creating HTML
    if (html && option_records[index]->advanced) continue;
    if (type == HighsOptionType::BOOL) {
      reportOption(file, ((OptionRecordBool*)option_records[index])[0],
                   report_only_non_default_values, html);
    } else if (type == HighsOptionType::INT) {
      reportOption(file, ((OptionRecordInt*)option_records[index])[0],
                   report_only_non_default_values, html);
    } else if (type == HighsOptionType::DOUBLE) {
      reportOption(file, ((OptionRecordDouble*)option_records[index])[0],
                   report_only_non_default_values, html);
    } else {
      reportOption(file, ((OptionRecordString*)option_records[index])[0],
                   report_only_non_default_values, html);
    }
  }
}

void reportOption(FILE* file, const OptionRecordBool& option,
                  const bool report_only_non_default_values, const bool html) {
  if (!report_only_non_default_values ||
      option.default_value != *option.value) {
    if (html) {
      fprintf(file,
              "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
              option.name.c_str());
      fprintf(file, "%s<br>\n", option.description.c_str());
      fprintf(file,
              "type: bool, advanced: %s, range: {false, true}, default: %s\n",
              bool2string(option.advanced), bool2string(option.default_value));
      fprintf(file, "</li>\n");
    } else {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(
          file,
          "# [type: bool, advanced: %s, range: {false, true}, default: %s]\n",
          bool2string(option.advanced), bool2string(option.default_value));
      fprintf(file, "%s = %s\n", option.name.c_str(),
              bool2string(*option.value));
    }
  }
}

void reportOption(FILE* file, const OptionRecordInt& option,
                  const bool report_only_non_default_values, const bool html) {
  if (!report_only_non_default_values ||
      option.default_value != *option.value) {
    if (html) {
      fprintf(file,
              "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
              option.name.c_str());
      fprintf(file, "%s<br>\n", option.description.c_str());
      fprintf(file, "type: int, advanced: %s, range: {%d, %d}, default: %d\n",
              bool2string(option.advanced), option.lower_bound,
              option.upper_bound, option.default_value);
      fprintf(file, "</li>\n");
    } else {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file,
              "# [type: int, advanced: %s, range: {%d, %d}, default: %d]\n",
              bool2string(option.advanced), option.lower_bound,
              option.upper_bound, option.default_value);
      fprintf(file, "%s = %d\n", option.name.c_str(), *option.value);
    }
  }
}

void reportOption(FILE* file, const OptionRecordDouble& option,
                  const bool report_only_non_default_values, const bool html) {
  if (!report_only_non_default_values ||
      option.default_value != *option.value) {
    if (html) {
      fprintf(file,
              "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
              option.name.c_str());
      fprintf(file, "%s<br>\n", option.description.c_str());
      fprintf(file,
              "type: double, advanced: %s, range: [%g, %g], default: %g\n",
              bool2string(option.advanced), option.lower_bound,
              option.upper_bound, option.default_value);
      fprintf(file, "</li>\n");
    } else {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file,
              "# [type: double, advanced: %s, range: [%g, %g], default: %g]\n",
              bool2string(option.advanced), option.lower_bound,
              option.upper_bound, option.default_value);
      fprintf(file, "%s = %g\n", option.name.c_str(), *option.value);
    }
  }
}

void reportOption(FILE* file, const OptionRecordString& option,
                  const bool report_only_non_default_values, const bool html) {
  // Don't report for the options file if writing to an options file
  if (option.name == options_file_string) return;
  if (!report_only_non_default_values ||
      option.default_value != *option.value) {
    if (html) {
      fprintf(file,
              "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
              option.name.c_str());
      fprintf(file, "%s<br>\n", option.description.c_str());
      fprintf(file, "type: string, advanced: %s, default: \"%s\"\n",
              bool2string(option.advanced), option.default_value.c_str());
      fprintf(file, "</li>\n");
    } else {
      fprintf(file, "\n# %s\n", option.description.c_str());
      fprintf(file, "# [type: string, advanced: %s, default: \"%s\"]\n",
              bool2string(option.advanced), option.default_value.c_str());
      fprintf(file, "%s = %s\n", option.name.c_str(), (*option.value).c_str());
    }
  }
}
