/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "io/LoadOptions.h"

#include <fstream>

#include "util/stringutil.h"

// For extended options to be parsed from a file. Assuming options file is
// specified.
bool loadOptionsFromFile(HighsOptions& options) {
  if (options.options_file.size() == 0) return false;

  string line, option, value;
  int line_count = 0;
  std::ifstream file(options.options_file);
  if (file.is_open()) {
    while (file.good()) {
      getline(file, line);
      line_count++;
      if (line.size() == 0 || line[0] == '#') continue;

      int equals = line.find_first_of("=");
      if (equals < 0 || equals >= (int)line.size() - 1) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Error on line %d of options file.", line_count);
        return false;
      }
      option = line.substr(0, equals);
      value = line.substr(equals + 1, line.size() - equals);
      trim(option);
      trim(value);
      if (setOptionValue(options.logfile, option, options.records, value) !=
          OptionStatus::OK)
        return false;
    }
  } else {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Options file not found.");
    return false;
  }

  return true;
}
