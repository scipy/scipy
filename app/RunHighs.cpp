/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HConfig.h"
#include "Highs.h"
#include "HighsIO.h"
#include "HighsMipSolver.h"
#include "HighsOptions.h"
#include "HighsRuntimeOptions.h"
#include "HighsTimer.h"
#include "LoadProblem.h"

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message = nullptr);
void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp);
void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, const Highs& highs);
HighsStatus callLpSolver(const HighsOptions& options, const HighsLp& lp,
                         FILE* output, int message_level, bool run_quiet);
HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp,
                          FILE* output, int message_level, bool run_quiet);

int main(int argc, char** argv) {
  printHighsVersionCopyright(stdout, ML_ALWAYS);

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  // Set message level.
  FILE* output;
  int message_level;
  output = options.output;
  message_level = options.message_level;

  bool run_quiet = false;  // true;//
  if (run_quiet) {
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "In main: running highs.run() quietly\n");
  }

  output = options.output;
  message_level = options.message_level;

  // Load problem.
  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  reportLpStatsOrError(output, message_level, read_status, lp);
  if (read_status == HighsStatus::Error) return (int)HighsStatus::Error;

  // Run LP or MIP solver.
  HighsStatus run_status = HighsStatus::Error;
  if (!options.mip) {
    run_status = callLpSolver(options, lp, output, message_level, run_quiet);
  } else {
    run_status = callMipSolver(options, lp, output, message_level, run_quiet);
  }

  return (int)run_status;
}

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message) {
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Running HiGHS %d.%d.%d [date: %s, git hash: %s]\n",
                    HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR,
                    HIGHS_VERSION_PATCH, HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Copyright (c) 2020 ERGO-Code under MIT licence terms\n\n");
#ifdef HiGHSDEV
  // Report on preprocessing macros
  if (message != nullptr) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "In %s\n", message);
  }
#ifdef OPENMP
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is not defined\n");
#endif

#ifdef SCIP_DEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is not defined\n");
#endif

#ifdef HiGHSDEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is not defined\n");
#endif
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Built with CMAKE_BUILD_TYPE=%s\n", CMAKE_BUILD_TYPE);
#endif
}

void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp) {
  if (read_status == HighsStatus::Error) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n");
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "LP       : %s\n",
                      lp.model_name_.c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Rows     : %d\n",
                      lp.numRow_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Cols     : %d\n",
                      lp.numCol_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Nonzeros : %d\n",
                      lp.Avalue_.size());
    if (lp.numInt_)
      HighsPrintMessage(output, message_level, ML_ALWAYS, "Integer  : %d\n",
                        lp.numInt_);
  }
}

void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, const Highs& highs) {
  if (run_status == HighsStatus::Error) {
    std::string statusname = HighsStatusToString(run_status);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n",
                      statusname.c_str());
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getHighsInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::OPTIMAL) {
        // The scaled model has been solved to optimality, but not the
        // unscaled model, flag this up, but report the scaled model
        // status
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Primal infeasibility: %10.3e (%d)\n",
                          highs_info.max_primal_infeasibility,
                          highs_info.num_primal_infeasibilities);
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Dual   infeasibility: %10.3e (%d)\n",
                          highs_info.max_dual_infeasibility,
                          highs_info.num_dual_infeasibilities);
        model_status = scaled_model_status;
      }
    }
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Model   status      : %s\n",
                      highs.highsModelStatusToString(model_status).c_str());
    /*
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Primal  status      : %s\n",
                      highs.highsPrimalDualStatusToString(highs_info.primal_status).c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Dual    status      : %s\n",
                      highs.highsPrimalDualStatusToString(highs_info.dual_status).c_str());
    */
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Simplex   iterations: %d\n",
                      highs_info.simplex_iteration_count);
    if (highs_info.ipm_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "IPM       iterations: %d\n",
                        highs_info.ipm_iteration_count);
    if (highs_info.crossover_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Crossover iterations: %d\n",
                        highs_info.crossover_iteration_count);
    if (model_status == HighsModelStatus::OPTIMAL) {
      double objective_function_value;
      highs.getHighsInfoValue("objective_function_value",
                              objective_function_value);
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Objective value     : %13.6e\n",
                        objective_function_value);
    }

    // Possibly write the solution to a file
    const HighsOptions& options = highs.getHighsOptions();
    if (options.write_solution_to_file)
      highs.writeSolution(options.solution_file, options.write_solution_pretty);
  }

  /*
  highs.writeSolution("", true);
  highs.writeSolution("", false);
  highs.writeHighsInfo("");
  highs.writeHighsInfo("HighsInfo.html");
  */
}

HighsStatus callLpSolver(const HighsOptions& options, const HighsLp& lp,
                         FILE* output, int message_level, bool run_quiet) {
  // Solve LP case.
  Highs highs;
  HighsStatus return_status = highs.passHighsOptions(options);
  if (return_status != HighsStatus::OK) {
    if (return_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from passHighsOptions\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "In main: fail return from passHighsOptions\n");
      return return_status;
    }
  }

  if (run_quiet) {
    highs.setHighsLogfile(NULL);
    highs.setHighsOutput(NULL);
  }

  HighsStatus init_status = highs.passModel(lp);
  if (init_status != HighsStatus::OK) {
    if (init_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return setting HighsLp\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Error setting HighsLp\n");
      return HighsStatus::Error;
    }
  }

  /*
  HighsStatus write_status;
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Writing model as MPS\n");
  write_status = highs.writeModel("write.mps");
  if (write_status != HighsStatus::OK) {
    if (write_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from highs.writeModel\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Error return from highs.writeModel\n");
    }
  }
  */

  // Write all the options to an options file
  // highs.writeHighsOptions("HighsOptions.set", false);
  // Write all the options as HTML
  // highs.writeHighsOptions("HighsOptions.html", false);
  // Possibly report options settings
  highs.writeHighsOptions("");  //, false);

  if (run_quiet)
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Before calling highs.run()\n");

  // Run HiGHS.
  HighsStatus run_status = highs.run();

  if (run_quiet)
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "After calling highs.run()\n");

  reportSolvedLpStats(output, message_level, run_status, highs);
  return run_status;
}

HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp,
                          FILE* output, int message_level, bool run_quiet) {
  HighsMipSolver solver(options, lp);
  HighsMipStatus status = solver.runMipSolver();
  switch (status) {
    case HighsMipStatus::kOptimal:
      return HighsStatus::OK;
    default:
      break;
  }
  return HighsStatus::Error;
}
