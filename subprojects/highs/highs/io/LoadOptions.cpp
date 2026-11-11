/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "io/LoadOptions.h"

#include <fstream>

#include "util/stringutil.h"

// For extended options to be parsed from a file. Assuming options file is
// specified.
HighsLoadOptionsStatus loadOptionsFromFile(
    const HighsLogOptions& report_log_options, HighsOptions& options,
    const std::string& filename) {
  if (filename.size() == 0) return HighsLoadOptionsStatus::kEmpty;

  string line, option, value;
  HighsInt line_count = 0;
  // loadOptionsFromFile needs its own non-chars string since the
  // default setting in io/stringutil.h excludes \" and \' that can
  // appear in an MPS name - the only other place where trim() is used
  const std::string non_chars = "\t\n\v\f\r\"\' ";
  std::ifstream file(filename);
  if (file.is_open()) {
    while (file.good()) {
      getline(file, line);
      line_count++;
      size_t non_space = line.find_first_not_of(' ');
      if (line.size() == 0 || line[0] == '#' || non_space == std::string::npos)
        continue;
      size_t equals = line.find_first_of('=');
      if (equals == std::string::npos || equals + 1 >= line.size()) {
        highsLogUser(report_log_options, HighsLogType::kError,
                     "Error on line %" HIGHSINT_FORMAT
                     " (\"%s\") of options file\n",
                     line_count, line.c_str());
        return HighsLoadOptionsStatus::kError;
      }
      option = line.substr(0, equals);
      value = line.substr(equals + 1, line.size() - equals);
      trim(option, non_chars);
      trim(value, non_chars);
      if (setLocalOptionValue(report_log_options, option, options.log_options,
                              options.records, value) != OptionStatus::kOk) {
        highsLogUser(report_log_options, HighsLogType::kError,
                     "Cannot read value \"%s\" for option \"%s\"\n",
                     value.c_str(), option.c_str());
        return HighsLoadOptionsStatus::kError;
      }
    }
  } else {
    highsLogUser(report_log_options, HighsLogType::kError,
                 "Options file not found\n");
    return HighsLoadOptionsStatus::kError;
  }

  return HighsLoadOptionsStatus::kOk;
}
