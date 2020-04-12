/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"

#ifdef OPENMP
#include "omp.h"
#endif

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"

#ifdef IPX_ON
#include "ipm/IpxWrapper.h"
#else
#include "ipm/IpxWrapperEmpty.h"
#endif

Highs::Highs() {
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const int value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsLogfile(FILE* logfile) {
  options_.logfile = logfile;
  return HighsStatus::OK;
}

HighsStatus Highs::setHighsOutput(FILE* output) {
  options_.output = output;
  return HighsStatus::OK;
}

HighsStatus Highs::readHighsOptions(const std::string filename) {
  if (filename.size() <= 0) {
    HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
                    "Empty file name so not reading options");
    return HighsStatus::Warning;
  }
  options_.options_file = filename;
  if (!loadOptionsFromFile(options_)) return HighsStatus::Error;
  return HighsStatus::OK;
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {
  if (passOptions(options_.logfile, options, options_) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

const HighsOptions& Highs::getHighsOptions() { return options_; }

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, int& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) ==
      OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsOptions(
    const std::string filename, const bool report_only_non_default_values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeHighsOptions", file, html);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  call_status = writeOptionsToFile(file, options_.records,
                                   report_only_non_default_values, html);
  return_status =
      interpretCallStatus(call_status, return_status, "writeOptionsToFile");
  return return_status;
}

const HighsOptions& Highs::getHighsOptions() const { return options_; }

const HighsInfo& Highs::getHighsInfo() const { return info_; }

HighsStatus Highs::getHighsInfoValue(const std::string& info, int& value) {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsInfoValue(const std::string& info,
                                     double& value) const {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsInfo(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeHighsInfo", file, html);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  call_status = writeInfoToFile(file, info_.records, html);
  return_status =
      interpretCallStatus(call_status, return_status, "writeInfoToFile");
  return return_status;
}

HighsStatus Highs::passModel(const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Copy the LP to the internal LP
  lp_ = lp;
  // Check validity of the LP, normalising its values (by default).
  call_status = assessLp(lp_, options_);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  return return_status;
}

HighsStatus Highs::readModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  Filereader* reader = Filereader::getFilereader(filename.c_str());
  HighsLp model;
  this->options_.model_file = filename;

  FilereaderRetcode call_code =
      reader->readModelFromFile(this->options_, model);
  delete reader;
  if (call_code != FilereaderRetcode::OK) {
    call_status = HighsStatus::Error;
    return_status =
        interpretCallStatus(call_status, return_status, "readModelFromFile");
    if (return_status == HighsStatus::Error) return return_status;
  }
  call_status = this->passModel(model);
  return_status = interpretCallStatus(call_status, return_status, "passModel");
  return return_status;
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp model = this->lp_;

  if (filename == "") {
    // Empty file name: report model on stdout
    reportLp(options_, model, 2);
    return_status = HighsStatus::OK;
  } else {
    Filereader* writer = Filereader::getFilereader(filename.c_str());
    call_status = writer->writeModelToFile(options_, filename.c_str(), model);
    delete writer;
    return_status =
        interpretCallStatus(call_status, return_status, "writeModelToFile");
  }
  return return_status;
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runLpSolver(..)
HighsStatus Highs::run() {
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
  assert(omp_max_threads > 0);
#ifdef HiGHSDEV
  if (omp_max_threads <= 0)
    printf("WARNING: omp_get_max_threads() returns %d\n", omp_max_threads);
  printf("Running with %d OMP thread(s)\n", omp_max_threads);
#endif
#endif
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  /*
if (options_.message_level >= 0) {
  printf("\n!! Actually solving an LP with %d cols, %d rows", lp_.numCol_,
lp_.numRow_); if (lp_.numCol_) printf(" and %d nonzeros",
lp_.Astart_[lp_.numCol_]); printf(":basis.valid_ = %d: basis_.valid_ = %d:
simplex_lp_status_.has_basis = %d!!\n\n", basis_.valid_, hmos_[0].basis_.valid_,
         hmos_[0].simplex_lp_status_.has_basis);
  if (basis_.valid_ != hmos_[0].basis_.valid_) {
    printf("NB %d = basis_.valid_ != hmos_[0].basis_.valid_ = %d\n",
basis_.valid_, hmos_[0].basis_.valid_);
  }
}
  */
  // If running as hsol, reset any changed options
  if (options_.run_as_hsol) setHsolOptions(options_);
  // Initialise the HiGHS model status values
  hmos_[0].scaled_model_status_ = HighsModelStatus::NOTSET;
  hmos_[0].unscaled_model_status_ = HighsModelStatus::NOTSET;
  model_status_ = hmos_[0].scaled_model_status_;
  scaled_model_status_ = hmos_[0].unscaled_model_status_;

#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  call_status = assessLp(lp_, options_);  //, normalise);
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
#endif

  // Return immediately if the LP has no columns
  if (!lp_.numCol_) {
    hmos_[0].unscaled_model_status_ = HighsModelStatus::MODEL_EMPTY;
    model_status_ = hmos_[0].unscaled_model_status_;
    return highsStatusFromHighsModelStatus(hmos_[0].unscaled_model_status_);
  }

  HighsSetIO(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.logfile, options_.records) != OptionStatus::OK)
    return HighsStatus::Error;
#endif
  HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                    "Solving %s\n", lp_.model_name_.c_str());

  double this_presolve_time = -1;
  double this_solve_presolved_lp_time = -1;
  double this_postsolve_time = -1;
  double this_solve_original_lp_time = -1;

  // Running as LP solver: start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  // Record the initial time and zero the overall iteration count
  double initial_time = timer_.readRunHighsClock();
  int postsolve_iteration_count = 0;
  // Define identifiers to refer to the HMO of the original LP (0) and
  // the HMO created when using presolve. The index of this HMO is 1
  // when solving a one-off LP, but greater than one if presolve has
  // been called multiple times. It's equal to the size of HMO
  const int original_hmo = 0;
  const int presolve_hmo = hmos_.size();
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  int solved_hmo = original_hmo;
  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  //  printf("\nHighs::run() 1: basis_.valid_ = %d\n", basis_.valid_);
  //  fflush(stdout);
  if (!basis_.valid_ && options_.presolve != off_string) {
    // No HiGHS basis so consider presolve
    hmos_[original_hmo].scaled_model_status_ = HighsModelStatus::NOTSET;
    // Presolve. runPresolve handles the level of presolving (0 = don't
    // presolve).

    //    printf("Writing before_presolve.mps\n");
    //    writeModel("before_presolve.mps");

    this_presolve_time = -timer_.read(timer_.presolve_clock);
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve, lp_, timer_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);
    this_presolve_time += timer_.read(timer_.presolve_clock);
    //    printf("\nHighs::run() 2: presolve status = %d\n",
    //    (int)presolve_status);fflush(stdout);
    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
        hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status =
            runLpSolver(hmos_[solved_hmo], "Not presolved: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::NotReduced: {
        hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        // Log the presolve reductions
        logPresolveReductions(hmos_[original_hmo].options_,
                              hmos_[original_hmo].lp_, false);
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status =
            runLpSolver(hmos_[solved_hmo],
                        "Problem not reduced by presolve: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_info.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::OK);
        call_status = cleanBounds(options_, reduced_lp);
        return_status =
            interpretCallStatus(call_status, return_status, "cleanBounds");
        if (return_status == HighsStatus::Error) return return_status;
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        logPresolveReductions(hmos_[original_hmo].options_,
                              hmos_[original_hmo].lp_, hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
        this_solve_presolved_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status =
            runLpSolver(hmos_[solved_hmo], "Solving the presolved LP");
        timer_.stop(timer_.solve_clock);
        this_solve_presolved_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        logPresolveReductions(hmos_[original_hmo].options_,
                              hmos_[original_hmo].lp_, true);
        hmos_[original_hmo].scaled_model_status_ = HighsModelStatus::OPTIMAL;
        hmos_[original_hmo].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
        // Proceed to postsolve.
        break;
      }
        //	printf("\nHighs::run() 3: presolve status = %d\n",
        //(int)presolve_status);fflush(stdout);
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        if (presolve_status == HighsPresolveStatus::Infeasible) {
          hmos_[original_hmo].scaled_model_status_ =
              HighsModelStatus::PRIMAL_INFEASIBLE;
          hmos_[original_hmo].unscaled_model_status_ =
              HighsModelStatus::PRIMAL_INFEASIBLE;
        } else {
          hmos_[original_hmo].scaled_model_status_ =
              HighsModelStatus::PRIMAL_UNBOUNDED;
          hmos_[original_hmo].unscaled_model_status_ =
              HighsModelStatus::PRIMAL_UNBOUNDED;
        }
        HighsLogMessage(
            options_.logfile, HighsMessageType::INFO,
            "Problem status detected on presolve: %s",
            highsModelStatusToString(hmos_[original_hmo].unscaled_model_status_)
                .c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

        model_status_ = hmos_[original_hmo].unscaled_model_status_;
        return HighsStatus::OK;
      }
      default: {
        // case HighsPresolveStatus::Error
        HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
                          "Presolve failed.");
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
        hmos_[original_hmo].unscaled_model_status_ =
            HighsModelStatus::PRESOLVE_ERROR;
        model_status_ = hmos_[original_hmo].unscaled_model_status_;
        return HighsStatus::Error;
      }
    }
    // Postsolve. Does nothing if there were no reductions during presolve.
    if (hmos_[solved_hmo].scaled_model_status_ == HighsModelStatus::OPTIMAL) {
      if (presolve_status == HighsPresolveStatus::Reduced ||
          presolve_status == HighsPresolveStatus::ReducedToEmpty) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
        presolve_info.reduced_solution_ = hmos_[solved_hmo].solution_;
        presolve_info.presolve_[0].setBasisInfo(
            hmos_[solved_hmo].basis_.col_status,
            hmos_[solved_hmo].basis_.row_status);
        // Run postsolve
        this_postsolve_time = -timer_.read(timer_.postsolve_clock);
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
        timer_.stop(timer_.postsolve_clock);
        this_postsolve_time += -timer_.read(timer_.postsolve_clock);
        if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
          HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                            "Postsolve finished\n");
          //
          // Now hot-start the simplex solver for the original_hmo:
          //
          // The original model hasn't been solved, so set up its solution
          // parameters
          resetModelStatusAndSolutionParams(hmos_[original_hmo]);
          // Set solution and its status
          hmos_[original_hmo].solution_ = presolve_info.recovered_solution_;
          //
          // Set basis and its status
          hmos_[original_hmo].basis_.col_status =
              presolve_info.presolve_[0].getColStatus();
          hmos_[original_hmo].basis_.row_status =
              presolve_info.presolve_[0].getRowStatus();
          hmos_[original_hmo].basis_.valid_ = true;
          analyseHighsBasicSolution(options_.logfile, hmos_[original_hmo],
                                    "after returning from postsolve");
          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions& options = hmos_[solved_hmo].options_;
          HighsOptions save_options = options;
          const bool full_logging = false;
          if (full_logging) options.message_level = ML_ALWAYS;
          // Force the use of simplex to clean up if IPM has been used
          // to solve the presolved problem
          if (options.solver == ipm_string) options.solver = simplex_string;
          options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
          // Ensure that the parallel solver isn't used
          options.highs_min_threads = 1;
          options.highs_max_threads = 1;
          hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
          int iteration_count0 =
              hmos_[solved_hmo]
                  .unscaled_solution_params_.simplex_iteration_count;
          this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
          timer_.start(timer_.solve_clock);
          call_status = runLpSolver(
              hmos_[solved_hmo],
              "Solving the original LP from the solution after postsolve");
          timer_.stop(timer_.solve_clock);
          this_solve_original_lp_time += timer_.read(timer_.solve_clock);
          return_status =
              interpretCallStatus(call_status, return_status, "runLpSolver");
          // Recover the options
          options = save_options;
          if (return_status == HighsStatus::Error) return return_status;
          int iteration_count1 =
              hmos_[solved_hmo]
                  .unscaled_solution_params_.simplex_iteration_count;
          postsolve_iteration_count = iteration_count1 - iteration_count0;
        }
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status
      hmos_[original_hmo].unscaled_model_status_ =
          hmos_[solved_hmo].unscaled_model_status_;
    }
  } else {
    // There is a valid basis for the problem or presolve is off
    solved_hmo = original_hmo;
    hmos_[solved_hmo].lp_.lp_name_ = "LP without presolve or with basis";
    // There is a valid HiGHS basis, so use it to initialise the basis
    // in the HMO to be solved
    if (basis_.valid_) hmos_[solved_hmo].basis_ = basis_;
    this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
    timer_.start(timer_.solve_clock);
    call_status = runLpSolver(hmos_[solved_hmo],
                              "Solving LP without presolve or with basis");
    timer_.stop(timer_.solve_clock);
    this_solve_original_lp_time += timer_.read(timer_.solve_clock);
    return_status =
        interpretCallStatus(call_status, return_status, "runLpSolver");
    if (return_status == HighsStatus::Error) return return_status;
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  //   assert(solved_hmo == original_hmo);
  // solved_hmo will be original_hmo unless the presolved LP is found to be
  // infeasible or unbounded

  if (!getHighsModelStatusAndInfo(solved_hmo)) return HighsStatus::Error;

  // Copy HMO solution/basis to HiGHS solution/basis: this resizes solution_ and
  // basis_ The HiGHS solution and basis have to come from the original_hmo for
  // them to have the right dimension.
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  // Report times
  if (hmos_[original_hmo].report_model_operations_clock) {
    std::vector<int> clockList{timer_.presolve_clock, timer_.solve_clock,
                               timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

  double lp_solve_final_time = timer_.readRunHighsClock();
  double this_solve_time = lp_solve_final_time - initial_time;
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Postsolve  : %d\n", postsolve_iteration_count);
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Time       : %0.3g\n", this_solve_time);
  if (this_solve_time > 0) {
    HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                      "For LP %16s",
                      hmos_[original_hmo].lp_.model_name_.c_str());
    double sum_time = 0;
    if (this_presolve_time > 0) {
      sum_time += this_presolve_time;
      int pct = (100 * this_presolve_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Presolve %0.3g (%3d%%)", this_presolve_time, pct);
    }
    if (this_solve_presolved_lp_time > 0) {
      sum_time += this_solve_presolved_lp_time;
      int pct = (100 * this_solve_presolved_lp_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Solve presolved LP %0.3g (%3d%%)",
                        this_solve_presolved_lp_time, pct);
    }
    if (this_postsolve_time > 0) {
      sum_time += this_postsolve_time;
      int pct = (100 * this_postsolve_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Postsolve %0.3g (%3d%%)", this_postsolve_time, pct);
    }
    if (this_solve_original_lp_time > 0) {
      sum_time += this_solve_original_lp_time;
      int pct = (100 * this_solve_original_lp_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Solve original LP %0.3g (%3d%%)",
                        this_solve_original_lp_time, pct);
    }
    HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                      "\n");
    double rlv_time_difference =
        fabs(sum_time - this_solve_time) / this_solve_time;
    if (rlv_time_difference > 0.1)
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        "Strange: Solve time = %g; Sum times = %g: relative "
                        "difference = %g\n",
                        this_solve_time, sum_time, rlv_time_difference);
  }
  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}

const HighsLp& Highs::getLp() const { return lp_; }

const HighsSolution& Highs::getSolution() const { return solution_; }

const HighsBasis& Highs::getBasis() const { return basis_; }

const HighsModelStatus& Highs::getModelStatus(const bool scaled_model) const {
  if (scaled_model) {
    return scaled_model_status_;
  } else {
    return model_status_;
  }
}

HighsStatus Highs::getBasicVariables(int* basic_variables) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_basis) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No basis available in getBasicVariables");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  int numCol = hmos_[0].lp_.numCol_;
  if (numRow != hmos_[0].simplex_lp_.numRow_) {
    HighsLogMessage(
        options_.logfile, HighsMessageType::ERROR,
        "Model LP and simplex LP row dimension difference (%d-%d=%d", numRow,
        hmos_[0].simplex_lp_.numRow_, numRow - hmos_[0].simplex_lp_.numRow_);
    return HighsStatus::Error;
  }
  for (int row = 0; row < numRow; row++) {
    int var = hmos_[0].simplex_basis_.basicIndex_[row];
    if (var < numCol) {
      basic_variables[row] = var;
    } else {
      basic_variables[row] = -(1 + var - numCol);
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseRow(const int row, double* row_vector,
                                      int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (row < 0 || row >= numRow) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getBasisInverseRow",
                    row, numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseRow");
    return HighsStatus::Error;
  }
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseCol(const int col, double* col_vector,
                                      int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (col < 0 || col >= numRow) {
    HighsLogMessage(
        options_.logfile, HighsMessageType::ERROR,
        "Column index %d out of range [0, %d] in getBasisInverseCol", col,
        numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseCol");
    return HighsStatus::Error;
  }
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[col] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 int* solution_num_nz, int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          int* solution_num_nz,
                                          int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisTransposeSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedRow(const int row, double* row_vector,
                                 int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (row < 0 || row >= hmos_[0].lp_.numRow_) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getReducedRow", row,
                    hmos_[0].lp_.numRow_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getReducedRow");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  vector<double> col_vector;
  vector<int> col_indices;
  int col_num_nz;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  col_vector.resize(numRow, 0);
  col_indices.resize(numRow, 0);
  HighsSimplexInterface simplex_interface(hmos_[0]);
  // Form B^{-T}e_{row}
  simplex_interface.basisSolve(rhs, &col_vector[0], &col_num_nz,
                               &col_indices[0], true);
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    double value = 0;
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      value += lp.Avalue_[el] * col_vector[row];
    }
    row_vector[col] = 0;
    if (fabs(value) > HIGHS_CONST_TINY) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedColumn(const int col, double* col_vector,
                                    int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (col < 0 || col >= hmos_[0].lp_.numCol_) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Column index %d out of range [0, %d] in getReducedColumn",
                    col, hmos_[0].lp_.numCol_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getReducedColumn");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    rhs[lp.Aindex_[el]] = lp.Avalue_[el];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("setSolution");
  // Check if solution is valid.
  assert((int)solution_.col_value.size() != 0 ||
         (int)solution_.col_value.size() != lp_.numCol_);
  assert((int)solution.col_dual.size() == 0 ||
         (int)solution.col_dual.size() == lp_.numCol_);
  assert((int)solution.row_dual.size() == 0 ||
         (int)solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  if (solution.col_value.size() > 0) {
    call_status = calculateRowValues(lp_, solution_);
    return_status =
        interpretCallStatus(call_status, return_status, "calculateRowValues");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (solution.row_dual.size() > 0) {
    call_status = calculateColDuals(lp_, solution_);
    return_status =
        interpretCallStatus(call_status, return_status, "calculateColDuals");
    if (return_status == HighsStatus::Error) return return_status;
  }
  return return_status;
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  underDevelopmentLogMessage("setBasis");
  if (!basisOk(options_.logfile, lp_, basis)) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "setBasis: invalid basis");
    return HighsStatus::Error;
  }
  basis_ = basis;
  basis_.valid_ = true;
  return HighsStatus::OK;
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int* indices,
                   const double* values) {
  int starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

bool Highs::addRows(const int num_new_row, const double* lower_bounds,
                    const double* upper_bounds, const int num_new_nz,
                    const int* starts, const int* indices,
                    const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("addRows");
  // Check that there is a HighsModelObject
  if (!haveHmo("addRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.addRows(num_new_row, lower_bounds, upper_bounds,
                                  num_new_nz, starts, indices, values);
  return_status = interpretCallStatus(call_status, return_status, "addRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const int num_new_nz,
                   const int* indices, const double* values) {
  int starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

bool Highs::addCols(const int num_new_col, const double* costs,
                    const double* lower_bounds, const double* upper_bounds,
                    const int num_new_nz, const int* starts, const int* indices,
                    const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("addCols");
  if (!haveHmo("addCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status =
      interface.addCols(num_new_col, costs, lower_bounds, upper_bounds,
                        num_new_nz, starts, indices, values);
  return_status = interpretCallStatus(call_status, return_status, "addCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeObjectiveSense(const ObjSense sense) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeObjectiveSense");
  if (!haveHmo("changeObjectiveSense")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeObjectiveSense(sense);
  return_status =
      interpretCallStatus(call_status, return_status, "changeObjectiveSense");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set,
                           const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsCost");
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(num_set_entries, set, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsCost");
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(mask, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColBounds(const int col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(num_set_entries, set, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int from_col, const int to_col,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(from_col, to_col, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(mask, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeRowBounds(const int row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeRowsBounds");
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(num_set_entries, set, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeRowsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeRowsBounds");
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(mask, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeCoeff(const int row, const int col, const double value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeCoeff");
  if (!haveHmo("changeCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCoefficient(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getObjectiveSense(ObjSense& sense) {
  underDevelopmentLogMessage("getObjectiveSense");
  if (!haveHmo("getObjectiveSense")) return false;
  sense = hmos_[0].lp_.sense_;
  return true;
}

bool Highs::getCols(const int from_col, const int to_col, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(from_col, to_col, num_col, costs, lower,
                                  upper, num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCols(const int n, const int* set, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(n, set, num_col, costs, lower, upper, num_nz,
                                  start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCols(const int* col_mask, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(col_mask, num_col, costs, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int from_row, const int to_row, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(from_row, to_row, num_row, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int num_set_entries, const int* set, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(num_set_entries, set, num_row, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int* mask, int& num_row, double* lower, double* upper,
                    int& num_nz, int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(mask, num_row, lower, upper, num_nz, start,
                                  index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCoeff(const int row, const int col, double& value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCoeff");
  if (!haveHmo("getCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);

  call_status = interface.getCoefficient(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "getCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(from_col, to_col);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(const int num_set_entries, const int* set) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(num_set_entries, set);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(mask);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(from_row, to_row);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(const int num_set_entries, const int* set) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(num_set_entries, set);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(mask);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

HighsStatus Highs::clearSolver() {
  underDevelopmentLogMessage("clearSolver");
  basis_.valid_ = false;
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void Highs::reportModelStatusSolutionBasis(const std::string message,
                                           const HighsModelStatus model_status,
                                           const HighsLp& lp,
                                           const HighsSolution& solution,
                                           const HighsBasis& basis) {
  printf(
      "\n%s\nModelStatus = %s; LP(%d, %d); solution (%d, %d; %d, %d); basis %d "
      "(%d, %d)\n\n",
      message.c_str(), utilHighsModelStatusToString(model_status).c_str(),
      lp.numCol_, lp.numRow_, (int)solution.col_value.size(),
      (int)solution.row_value.size(), (int)solution.col_dual.size(),
      (int)solution.row_dual.size(), basis.valid_, (int)basis.col_status.size(),
      (int)basis.row_status.size());
}
#endif

std::string Highs::highsModelStatusToString(
    const HighsModelStatus model_status) const {
  return utilHighsModelStatusToString(model_status);
}

std::string Highs::highsPrimalDualStatusToString(const int primal_dual_status) {
  return utilPrimalDualStatusToString(primal_dual_status);
}

// Private methods
HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolve == off_string) return HighsPresolveStatus::NotPresolved;

  if (info.lp_ == nullptr) return HighsPresolveStatus::NullError;

  if (info.presolve_.size() == 0) return HighsPresolveStatus::NotReduced;

  info.presolve_[0].load(*(info.lp_));

  // Initialize a new presolve class instance for the LP given in presolve info
  HighsPresolveStatus presolve_return_status = info.presolve_[0].presolve();

  if (presolve_return_status == HighsPresolveStatus::Reduced &&
      info.lp_->sense_ == ObjSense::MAXIMIZE)
    info.negateReducedCosts();

  return presolve_return_status;
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo& info) {
  if (info.presolve_.size() != 0) {
    // Handle max case.
    if (info.lp_->sense_ == ObjSense::MAXIMIZE) info.negateColDuals(true);

    bool solution_ok =
        isSolutionConsistent(info.getReducedProblem(), info.reduced_solution_);
    if (!solution_ok)
      return HighsPostsolveStatus::ReducedSolutionDimenionsError;

    // todo: error handling + see todo in run()
    info.presolve_[0].postsolve(info.reduced_solution_,
                                info.recovered_solution_);

    if (info.lp_->sense_ == ObjSense::MAXIMIZE) info.negateColDuals(false);

    return HighsPostsolveStatus::SolutionRecovered;
  } else {
    return HighsPostsolveStatus::NoPostsolve;
  }
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runLpSolver(HighsModelObject& model, const string message) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Reset unscaled and scaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndSolutionParams(model);
  HighsLogMessage(options_.logfile, HighsMessageType::INFO, message.c_str());
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  call_status = assessLp(lp_, options_);
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
#endif
  if (!model.lp_.numRow_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::Error) return return_status;
  } else if (options_.solver == ipm_string) {
    // Use IPM
#ifdef IPX_ON
    bool imprecise_solution;
    call_status =
        solveLpIpx(options_, timer_, model.lp_, imprecise_solution,
                   model.basis_, model.solution_, model.unscaled_model_status_,
                   model.unscaled_solution_params_);
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpIpx");
    if (return_status == HighsStatus::Error) return return_status;
    if (imprecise_solution) {
      // IPX+crossover has not obtained a solution satisfying the tolerances.
      // Use the simplex method to clean up
      call_status = solveLpSimplex(model);
      return_status =
          interpretCallStatus(call_status, return_status, "solveLpSimplex");
      if (return_status == HighsStatus::Error) return return_status;

      if (!isSolutionConsistent(model.lp_, model.solution_)) {
        HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                        "Inconsistent solution returned from solver");
        return HighsStatus::Error;
      }
    } else {
      // Set the scaled model status and solution params for completeness
      model.scaled_model_status_ = model.unscaled_model_status_;
      model.scaled_solution_params_ = model.unscaled_solution_params_;
    }
#else
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Model cannot be solved with IPM");
    return HighsStatus::Error;
#endif
  } else {
    // Use Simplex
    call_status = solveLpSimplex(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpSimplex");
    if (return_status == HighsStatus::Error) return return_status;

    if (!isSolutionConsistent(model.lp_, model.solution_)) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "Inconsistent solution returned from solver");
      return HighsStatus::Error;
    }
  }
  call_status = analyseHighsBasicSolution(
      options_.logfile, model.lp_, model.basis_, model.solution_,
      model.unscaled_model_status_, model.unscaled_solution_params_, message);
  return_status = interpretCallStatus(call_status, return_status,
                                      "analyseHighsBasicSolution");
  return return_status;
}

HighsStatus Highs::writeSolution(const std::string filename,
                                 const bool pretty) const {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  HighsBasis basis = this->basis_;
  HighsSolution solution = this->solution_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeSolution", file, html);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  std::cout << "warning: Feature under development" << std::endl;
  return HighsStatus::Warning;

  // writeSolutionToFile(file, lp, basis, solution, pretty);
  // return HighsStatus::OK;
}

bool Highs::updateHighsSolutionBasis() {
  if (!haveHmo("updateHighsSolutionBasis")) return false;
  solution_.col_value.resize(lp_.numCol_);
  solution_.row_value.resize(lp_.numRow_);
  solution_.col_dual.resize(lp_.numCol_);
  solution_.row_dual.resize(lp_.numRow_);
  hmos_[0].solution_.col_value.resize(lp_.numCol_);
  hmos_[0].solution_.row_value.resize(lp_.numRow_);
  hmos_[0].solution_.col_dual.resize(lp_.numCol_);
  hmos_[0].solution_.row_dual.resize(lp_.numRow_);

  if (hmos_[0].basis_.valid_) {
    basis_ = hmos_[0].basis_;
  } else {
    basis_.valid_ = false;
    basis_.col_status.resize(lp_.numCol_);
    basis_.row_status.resize(lp_.numRow_);
  }
  return true;
}

bool Highs::getHighsModelStatusAndInfo(const int solved_hmo) {
  if (!haveHmo("getHighsModelStatusAndInfo")) return false;

  model_status_ = hmos_[solved_hmo].unscaled_model_status_;
  scaled_model_status_ = hmos_[solved_hmo].scaled_model_status_;

  HighsSolutionParams& solution_params =
      hmos_[solved_hmo].unscaled_solution_params_;

  // Get the total simplex IPM and crossover iteration counts over all HMO
  info_.simplex_iteration_count = 0;
  info_.ipm_iteration_count = 0;
  info_.crossover_iteration_count = 0;
  int hmos_size = hmos_.size();
  for (int k = 0; k < hmos_size; k++) {
    info_.simplex_iteration_count +=
        hmos_[k].unscaled_solution_params_.simplex_iteration_count;
    info_.ipm_iteration_count +=
        hmos_[k].unscaled_solution_params_.ipm_iteration_count;
    info_.crossover_iteration_count +=
        hmos_[k].unscaled_solution_params_.crossover_iteration_count;
  }
  info_.primal_status = solution_params.primal_status;
  info_.dual_status = solution_params.dual_status;
  info_.objective_function_value = solution_params.objective_function_value;
  info_.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  info_.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  info_.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
  info_.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  info_.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  info_.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
  return true;
}

HighsStatus Highs::openWriteFile(const string filename,
                                 const string method_name, FILE*& file,
                                 bool& html) const {
  html = false;
  if (filename == "") {
    // Empty file name: use stdout
    file = stdout;
  } else {
    file = fopen(filename.c_str(), "w");
    if (file == 0) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "Cannot open writeable file \"%s\" in %s",
                      filename.c_str(), method_name.c_str());
      return HighsStatus::Error;
    }
    const char* dot = strrchr(filename.c_str(), '.');
    if (dot && dot != filename) html = strcmp(dot + 1, "html") == 0;
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getUseModelStatus(
    HighsModelStatus& use_model_status,
    const double unscaled_primal_feasibility_tolerance,
    const double unscaled_dual_feasibility_tolerance,
    const bool rerun_from_logical_basis) {
  if (model_status_ != HighsModelStatus::NOTSET) {
    use_model_status = model_status_;
  } else {
    // Handle the case where the status of the unscaled model is not set
    HighsStatus return_status = HighsStatus::OK;
    HighsStatus call_status;
    const double report = false;  // true;//
    if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                        unscaled_dual_feasibility_tolerance, report)) {
      use_model_status = HighsModelStatus::OPTIMAL;
    } else if (rerun_from_logical_basis) {
      std::string save_presolve = options_.presolve;
      basis_.valid_ = false;
      options_.presolve = on_string;
      call_status = run();
      return_status = interpretCallStatus(call_status, return_status, "run()");
      options_.presolve = save_presolve;
      if (return_status == HighsStatus::Error) return return_status;

      if (report)
        printf(
            "Unscaled model status was NOTSET: after running from logical "
            "basis it is %s\n",
            highsModelStatusToString(model_status_).c_str());

      if (model_status_ != HighsModelStatus::NOTSET) {
        use_model_status = model_status_;
      } else if (unscaledOptimal(unscaled_primal_feasibility_tolerance,
                                 unscaled_dual_feasibility_tolerance, report)) {
        use_model_status = HighsModelStatus::OPTIMAL;
      }
    } else {
      // Nothing to be done: use original unscaled model status
      use_model_status = model_status_;
    }
  }
  return HighsStatus::OK;
}

bool Highs::unscaledOptimal(const double unscaled_primal_feasibility_tolerance,
                            const double unscaled_dual_feasibility_tolerance,
                            const bool report) {
  if (scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    const double max_primal_infeasibility = info_.max_primal_infeasibility;
    const double max_dual_infeasibility = info_.max_dual_infeasibility;
    if (report)
      printf(
          "Scaled model status is OPTIMAL: max unscaled (primal / dual) "
          "infeasibilities are (%g / %g)\n",
          max_primal_infeasibility, max_dual_infeasibility);
    if ((max_primal_infeasibility > unscaled_primal_feasibility_tolerance) ||
        (max_dual_infeasibility > unscaled_dual_feasibility_tolerance)) {
      printf(
          "Use model status of NOTSET since max unscaled (primal / dual) "
          "infeasibilities are (%g / %g)\n",
          max_primal_infeasibility, max_dual_infeasibility);
    } else {
      if (report)
        printf(
            "Set unscaled model status to OPTIMAL since unscaled "
            "infeasibilities are tolerable\n");
      return true;
    }
  }
  return false;
}

bool Highs::haveHmo(const string method_name) {
  bool have_hmo = hmos_.size() > 0;
  assert(have_hmo);
#ifdef HiGHSDEV
  if (!have_hmo)
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Method %s called without any HighsModelObject",
                    method_name.c_str());
#endif
  return have_hmo;
}

void Highs::underDevelopmentLogMessage(const string method_name) {
  HighsLogMessage(
      options_.logfile, HighsMessageType::WARNING,
      "Method %s is still under development and behaviour may be unpredictable",
      method_name.c_str());
}
