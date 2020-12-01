/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LP_DATA_HIGHSRUNTIMEOPTIONS_H_
#define LP_DATA_HIGHSRUNTIMEOPTIONS_H_

#include "cxxopts.hpp"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "util/stringutil.h"

bool loadOptions(int argc, char** argv, HighsOptions& options) {
  try {
    cxxopts::Options cxx_options(argv[0], "HiGHS options");
    cxx_options.positional_help("[file]").show_positional_help();

    std::string presolve, solver, parallel;

    cxx_options.add_options()(model_file_string, "File of model to solve.",
                              cxxopts::value<std::vector<std::string>>())(
        presolve_string,
        "Presolve: \"choose\" by default - \"on\"/\"off\" are alternatives.",
        cxxopts::value<std::string>(presolve))(
        solver_string,
        "Solver: \"choose\" by default - \"simplex\"/\"ipm\" are alternatives.",
        cxxopts::value<std::string>(solver))(
        parallel_string,
        "Parallel solve: \"choose\" by default - \"on\"/\"off\" are "
        "alternatives.",
        cxxopts::value<std::string>(parallel))(time_limit_string,
                                               "Run time limit (double).",
                                               cxxopts::value<double>())(
        options_file_string, "File containing HiGHS options.",
        cxxopts::value<std::vector<std::string>>())("h, help", "Print help.");
    cxx_options.parse_positional("model_file");

    auto result = cxx_options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << cxx_options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count(model_file_string)) {
      auto& v = result[model_file_string].as<std::vector<std::string>>();
      if (v.size() > 1) {
        int nonEmpty = 0;
        for (int i = 0; i < (int)v.size(); i++) {
          std::string arg = v[i];
          if (trim(arg).size() > 0) {
            nonEmpty++;
            options.model_file = arg;
          }
        }
        if (nonEmpty > 1) {
          std::cout << "Multiple files not implemented.\n";
          return false;
        }
      } else {
        options.model_file = v[0];
      }
    }

    if (result.count(presolve_string)) {
      std::string value = result[presolve_string].as<std::string>();
      if (setOptionValue(options.logfile, presolve_string, options.records,
                         value) != OptionStatus::OK)
        return false;
    }

    if (result.count(solver_string)) {
      std::string value = result[solver_string].as<std::string>();
      if (setOptionValue(options.logfile, solver_string, options.records,
                         value) != OptionStatus::OK)
        return false;
    }

    if (result.count(parallel_string)) {
      std::string value = result[parallel_string].as<std::string>();
      if (setOptionValue(options.logfile, parallel_string, options.records,
                         value) != OptionStatus::OK)
        return false;
    }

    if (result.count(time_limit_string)) {
      double value = result[time_limit_string].as<double>();
      if (setOptionValue(options.logfile, time_limit_string, options.records,
                         value) != OptionStatus::OK)
        return false;
    }

    if (result.count(options_file_string)) {
      auto& v = result[options_file_string].as<std::vector<std::string>>();
      if (v.size() > 1) {
        std::cout << "Multiple options files not implemented.\n";
        return false;
      }
      options.options_file = v[0];
      if (!loadOptionsFromFile(options)) return false;
    }

  } catch (const cxxopts::OptionException& e) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Error parsing options: %s", e.what());
    return false;
  }

  if (options.model_file.size() == 0) {
    std::cout << "Please specify filename in .mps|.lp|.ems format.\n";
    return false;
  }

  return true;
}

#endif
