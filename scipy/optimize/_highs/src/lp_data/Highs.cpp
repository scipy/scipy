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
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HighsSimplexInterface.h"
#include "util/HighsMatrixPic.h"

#ifdef OPENMP
#include "omp.h"
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

HighsStatus Highs::resetHighsOptions() {
  resetOptions(options_.records);
  return HighsStatus::OK;
}

HighsStatus Highs::writeHighsOptions(
    const std::string filename, const bool report_only_non_default_values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  return_status = interpretCallStatus(
      openWriteFile(filename, "writeHighsOptions", file, html), return_status,
      "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  return_status = interpretCallStatus(
      writeOptionsToFile(file, options_.records, report_only_non_default_values,
                         html),
      return_status, "writeOptionsToFile");
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
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  return_status =
      interpretCallStatus(openWriteFile(filename, "writeHighsInfo", file, html),
                          return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  return_status =
      interpretCallStatus(writeInfoToFile(file, info_.records, html),
                          return_status, "writeInfoToFile");
  return return_status;
}

// Methods below change the incumbent model or solver infomation
// associated with it. Hence returnFromHighs is called at the end of
// each
HighsStatus Highs::reset() {
  HighsStatus return_status = HighsStatus::OK;
  // Clear the status, solution, basis and info associated with any previous
  // model
  return_status =
      interpretCallStatus(clearSolver(), return_status, "clearSolver");
  if (return_status == HighsStatus::Error) return return_status;
  // Clear any HiGHS model object
  hmos_.clear();
  // Create a HiGHS model object for this LP
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));

  presolve_.clear();

  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::OK;
  // Copy the LP to the internal LP
  lp_ = lp;
  // Check validity of the LP, normalising its values
  return_status =
      interpretCallStatus(assessLp(lp_, options_), return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
  // Clear solver status, solution, basis and info associated with any
  // previous model; clear any HiGHS model object; create a HiGHS
  // model object for this LP
  return_status = interpretCallStatus(reset(), return_status, "reset");
  return returnFromHighs(return_status);
}

HighsStatus Highs::passModel(const int num_col, const int num_row,
                             const int num_nz, const double* costs,
                             const double* col_lower, const double* col_upper,
                             const double* row_lower, const double* row_upper,
                             const int* astart, const int* aindex,
                             const double* avalue) {
  HighsLp lp;
  lp.numCol_ = num_col;
  lp.numRow_ = num_row;
  if (num_col > 0) {
    assert(costs != NULL);
    assert(col_lower != NULL);
    assert(col_upper != NULL);
    lp.colCost_.assign(costs, costs + num_col);
    lp.colLower_.assign(col_lower, col_lower + num_col);
    lp.colUpper_.assign(col_upper, col_upper + num_col);
  }
  if (num_row > 0) {
    assert(row_lower != NULL);
    assert(row_upper != NULL);
    lp.rowLower_.assign(row_lower, row_lower + num_row);
    lp.rowUpper_.assign(row_upper, row_upper + num_row);
  }
  if (num_nz > 0) {
    assert(num_col > 0);
    assert(num_row > 0);
    assert(astart != NULL);
    assert(aindex != NULL);
    assert(avalue != NULL);
    lp.Astart_.assign(astart, astart + num_col);
    lp.Aindex_.assign(aindex, aindex + num_nz);
    lp.Avalue_.assign(avalue, avalue + num_nz);
  }
  lp.Astart_.resize(num_col + 1);
  lp.Astart_[num_col] = num_nz;
  return this->passModel(lp);
}

HighsStatus Highs::readModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  Filereader* reader = Filereader::getFilereader(filename);
  if (reader == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Model file %s not supported", filename.c_str());
    return HighsStatus::Error;
  }

  HighsLp model;
  this->options_.model_file = filename;

  FilereaderRetcode call_code =
      reader->readModelFromFile(this->options_, model);
  delete reader;
  if (call_code != FilereaderRetcode::OK) {
    interpretFilereaderRetcode(this->options_.logfile, filename.c_str(),
                               call_code);
    return_status = interpretCallStatus(HighsStatus::Error, return_status,
                                        "readModelFromFile");
    if (return_status == HighsStatus::Error) return return_status;
  }
  model.model_name_ = extractModelName(filename);
  return_status =
      interpretCallStatus(this->passModel(model), return_status, "passModel");
  return returnFromHighs(return_status);
}

HighsStatus Highs::clearModel() {
  HighsStatus return_status = HighsStatus::OK;
  // Remove all HighsModelObject entries
  hmos_.clear();
  // Set up with an empty LP so that addrows/cols can be used to build
  //  HighsLp empty_lp; lp_ = empty_lp;
  lp_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  return_status =
      interpretCallStatus(this->clearSolver(), return_status, "clearSolver");
  if (return_status == HighsStatus::Error) return return_status;
  return returnFromHighs(return_status);
}

HighsStatus Highs::readBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  // Try to read basis file into read_basis
  HighsBasis read_basis = this->basis_;
  return_status =
      interpretCallStatus(readBasisFile(options_, read_basis, filename),
                          return_status, "readBasis");
  if (return_status != HighsStatus::OK) return return_status;
  // Basis read OK: check whether it's consistent with the LP
  if (!isBasisConsistent(lp_, read_basis)) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "readBasis: invalid basis");
    return HighsStatus::Error;
  }
  // Update the HiGHS basis and invalidate any simplex basis for the model
  this->basis_ = read_basis;
  this->basis_.valid_ = true;
  if (hmos_.size() > 0) {
    HighsSimplexInterface interface(hmos_[0]);
    interface.clearBasis();
  }
  // Can't use returnFromHighs since...
  return HighsStatus::OK;
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsLp model = this->lp_;

  if (filename == "") {
    // Empty file name: report model on stdout
    reportLp(options_, model, 2);
    return_status = HighsStatus::OK;
  } else {
    Filereader* writer = Filereader::getFilereader(filename);
    if (writer == NULL) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "Model file %s not supported", filename.c_str());
      return HighsStatus::Error;
    }
    return_status =
        interpretCallStatus(writer->writeModelToFile(options_, filename, model),
                            return_status, "writeModelToFile");
    delete writer;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::writeBasis(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  return_status =
      interpretCallStatus(writeBasisFile(options_, this->basis_, filename),
                          return_status, "writeBasis");
  return returnFromHighs(return_status);
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runLpSolver(..)
HighsStatus Highs::run() {
#ifdef HiGHSDEV
  const int min_highs_debug_level = HIGHS_DEBUG_LEVEL_MIN;
  //      HIGHS_DEBUG_LEVEL_CHEAP;
  // HIGHS_DEBUG_LEVEL_COSTLY;
  // HIGHS_DEBUG_LEVEL_MAX;
  if (options_.highs_debug_level < min_highs_debug_level) {
    printf(
        "Highs::run() HiGHSDEV define so switching options_.highs_debug_level "
        "from %d to %d\n",
        options_.highs_debug_level, min_highs_debug_level);
    options_.highs_debug_level = min_highs_debug_level;
  }
  writeModel("HighsRunModel.mps");
  //  if (lp_.numRow_>0 && lp_.numCol_>0) writeLpMatrixPicToFile(options_,
  //  "LpMatrix", lp_);
#endif

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
  // Zero the HiGHS iteration counts
  zeroHighsIterationCounts(info_);
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
  // Determine whether a model has been loaded.
  assert((int)hmos_.size() <= 1);
  if (hmos_.size() == 0) {
    // No Highs model object, so load model according to value of
    // model_file
    if (options_.model_file.compare(FILENAME_DEFAULT) == 0) {
      // model_file is still default value, so return with error
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "No model can be loaded in run()");
      return_status = HighsStatus::Error;
      return returnFromRun(return_status);
    } else {
      std::string model_file = options_.model_file;
      call_status = readModel(model_file);
      return_status =
          interpretCallStatus(call_status, return_status, "readModel");
      if (return_status == HighsStatus::Error)
        return returnFromRun(return_status);
    }
  }
  // Ensure that there is exactly one Highs model object
  assert((int)hmos_.size() == 1);

  // Initialise the HiGHS model status values
  hmos_[0].scaled_model_status_ = HighsModelStatus::NOTSET;
  hmos_[0].unscaled_model_status_ = HighsModelStatus::NOTSET;
  model_status_ = hmos_[0].scaled_model_status_;
  scaled_model_status_ = hmos_[0].unscaled_model_status_;

#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  call_status = assessLp(lp_, options_);
  // If any errors have been found or normalisation carried out,
  // call_status will be ERROR or WARNING, so only valid return is OK.
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return returnFromRun(return_status);
#endif

  // Return immediately if the LP has no columns
  if (!lp_.numCol_) {
    model_status_ = HighsModelStatus::MODEL_EMPTY;
    scaled_model_status_ = model_status_;
    hmos_[0].unscaled_model_status_ = model_status_;
    hmos_[0].scaled_model_status_ = model_status_;
    return_status = highsStatusFromHighsModelStatus(model_status_);
    return returnFromRun(return_status);
  }

  HighsSetIO(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.logfile, options_.records) != OptionStatus::OK) {
    return_status = HighsStatus::Error;
    return returnFromRun(return_status);
  }
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
  // Record the initial time and set the postsolve iteration count to
  // -1 to identify whether it's not required
  double initial_time = timer_.readRunHighsClock();
  int postsolve_iteration_count = -1;
  // Define identifiers to refer to the HMO of the original LP (0) and
  // the HMO created when using presolve. The index of this HMO is 1
  // when solving a one-off LP, but greater than one if presolve has
  // been called multiple times. It's equal to the size of HMO
  const int original_hmo = 0;
  const int presolve_hmo = hmos_.size();
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  int solved_hmo = original_hmo;

  if (!basis_.valid_ && options_.presolve != off_string) {
    // No HiGHS basis so consider presolve
    //
    // If using IPX to solve the reduced LP, crossover must be run
    // since a basic solution is required by postsolve
    if (options_.solver == ipm_string && !options_.run_crossover) {
      HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
                      "Forcing IPX to use crossover after presolve");
      options_.run_crossover = true;
    }

    hmos_[original_hmo].scaled_model_status_ = HighsModelStatus::NOTSET;
    // Presolve. runPresolve handles the level of presolving (0 = don't
    // presolve).

    //    printf("Writing before_presolve.mps\n");
    //    writeModel("before_presolve.mps");

    // Run and time presolve.
    const double from_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time = -from_presolve_time;
    timer_.start(timer_.presolve_clock);

    HighsPresolveStatus presolve_status = runPresolve();
    timer_.stop(timer_.presolve_clock);
    const double to_presolve_time = timer_.read(timer_.presolve_clock);
    this_presolve_time += to_presolve_time;
    presolve_.info_.presolve_time = this_presolve_time;

    // Set an illegal local pivot threshold value that's updated after
    // solving the presolved LP - if simplex is used
    double factor_pivot_threshold = -1;

    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
        hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = runLpSolver(solved_hmo, "Not presolved: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::NotReduced: {
        hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        // Log the presolve reductions
        reportPresolveReductions(hmos_[original_hmo].options_,
                                 hmos_[original_hmo].lp_, false);
        this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = runLpSolver(
            solved_hmo, "Problem not reduced by presolve: solving the LP");
        timer_.stop(timer_.solve_clock);
        this_solve_original_lp_time += timer_.read(timer_.solve_clock);
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::OK);
        call_status = cleanBounds(options_, reduced_lp);
        // Ignore any warning from clean bounds since the original LP
        // is still solved after presolve
        if (interpretCallStatus(call_status, return_status, "cleanBounds") ==
            HighsStatus::Error)
          return HighsStatus::Error;
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        reportPresolveReductions(hmos_[original_hmo].options_,
                                 hmos_[original_hmo].lp_,
                                 hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
        // Don't try dual cut-off when solving the presolved LP, as the
        // objective values aren't correct
        //	HighsOptions& options = hmos_[solved_hmo].options_;
        //	HighsOptions save_options = options;
        const double save_dual_objective_value_upper_bound =
            options_.dual_objective_value_upper_bound;
        options_.dual_objective_value_upper_bound = HIGHS_CONST_INF;
        this_solve_presolved_lp_time = -timer_.read(timer_.solve_clock);
        timer_.start(timer_.solve_clock);
        call_status = runLpSolver(solved_hmo, "Solving the presolved LP");
        timer_.stop(timer_.solve_clock);
        this_solve_presolved_lp_time += timer_.read(timer_.solve_clock);
        if (hmos_[solved_hmo].simplex_lp_status_.valid) {
          // Record the pivot threshold resulting from solving the presolved LP
          // with simplex
          factor_pivot_threshold =
              hmos_[solved_hmo].simplex_info_.factor_pivot_threshold;
        }
        // Restore the dual objective cut-off
        options_.dual_objective_value_upper_bound =
            save_dual_objective_value_upper_bound;
        return_status =
            interpretCallStatus(call_status, return_status, "runLpSolver");
        if (return_status == HighsStatus::Error)
          return returnFromRun(return_status);

        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        reportPresolveReductions(hmos_[original_hmo].options_,
                                 hmos_[original_hmo].lp_, true);
        hmos_[original_hmo].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
        hmos_[original_hmo].scaled_model_status_ =
            hmos_[original_hmo].unscaled_model_status_;
        // Proceed to postsolve.
        break;
      }
        //	printf("\nHighs::run() 3: presolve status = %d\n",
        //(int)presolve_status);fflush(stdout);
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        if (presolve_status == HighsPresolveStatus::Infeasible) {
          model_status_ = HighsModelStatus::PRIMAL_INFEASIBLE;
        } else {
          model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
        }
        HighsLogMessage(options_.logfile, HighsMessageType::INFO,
                        "Problem status detected on presolve: %s",
                        highsModelStatusToString(model_status_).c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

        // Transfer the model status to the scaled model status and orriginal
        // HMO statuses;
        scaled_model_status_ = model_status_;
        hmos_[original_hmo].unscaled_model_status_ = model_status_;
        hmos_[original_hmo].scaled_model_status_ = model_status_;
        return_status = HighsStatus::OK;
        return returnFromRun(return_status);
      }
      case HighsPresolveStatus::Timeout: {
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
                          "Presolve reached timeout\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        return HighsStatus::Warning;
      }
      case HighsPresolveStatus::OptionsError: {
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
                          "Presolve options error.\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        return HighsStatus::Warning;
      }
      default: {
        // case HighsPresolveStatus::Error
        model_status_ = HighsModelStatus::PRESOLVE_ERROR;
        HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
                          "Presolve failed.\n");
        if (run_highs_clock_already_running) timer_.stopRunHighsClock();
        // Transfer the model status to the scaled model status and orriginal
        // HMO statuses;
        scaled_model_status_ = model_status_;
        hmos_[original_hmo].unscaled_model_status_ = model_status_;
        hmos_[original_hmo].scaled_model_status_ = model_status_;
        return_status = HighsStatus::Error;
        return returnFromRun(return_status);
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
        if (presolve_status == HighsPresolveStatus::ReducedToEmpty) {
          clearSolutionUtil(hmos_[solved_hmo].solution_);
          clearBasisUtil(hmos_[solved_hmo].basis_);
        }

        presolve_.data_.reduced_solution_ = hmos_[solved_hmo].solution_;
        presolve_.data_.reduced_basis_.col_status =
            hmos_[solved_hmo].basis_.col_status;
        presolve_.data_.reduced_basis_.row_status =
            hmos_[solved_hmo].basis_.row_status;

        this_postsolve_time = -timer_.read(timer_.postsolve_clock);
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve();
        timer_.stop(timer_.postsolve_clock);
        this_postsolve_time += -timer_.read(timer_.postsolve_clock);
        presolve_.info_.postsolve_time = this_postsolve_time;

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
          hmos_[original_hmo].solution_ = presolve_.data_.recovered_solution_;

          // Set basis and its status
          hmos_[original_hmo].basis_.valid_ = true;
          hmos_[original_hmo].basis_.col_status =
              presolve_.data_.recovered_basis_.col_status;
          hmos_[original_hmo].basis_.row_status =
              presolve_.data_.recovered_basis_.row_status;

          // Possibly force debug to perform KKT check on what's
          // returned from postsolve
          const bool force_debug = false;
          int save_highs_debug_level = options_.highs_debug_level;
          if (force_debug)
            options_.highs_debug_level = HIGHS_DEBUG_LEVEL_COSTLY;
          debugHighsBasicSolution("After returning from postsolve", options_,
                                  lp_, hmos_[original_hmo].basis_,
                                  hmos_[original_hmo].solution_);
          options_.highs_debug_level = save_highs_debug_level;

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
          // Use any pivot threshold resulting from solving the presolved LP
          if (factor_pivot_threshold > 0)
            options.factor_pivot_threshold = factor_pivot_threshold;

          hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
          int iteration_count0 = info_.simplex_iteration_count;
          this_solve_original_lp_time = -timer_.read(timer_.solve_clock);
          timer_.start(timer_.solve_clock);
          call_status = runLpSolver(
              solved_hmo,
              "Solving the original LP from the solution after postsolve");
          timer_.stop(timer_.solve_clock);
          this_solve_original_lp_time += timer_.read(timer_.solve_clock);
          return_status =
              interpretCallStatus(call_status, return_status, "runLpSolver");
          // Recover the options
          options = save_options;
          if (return_status == HighsStatus::Error)
            return returnFromRun(return_status);
          postsolve_iteration_count =
              info_.simplex_iteration_count - iteration_count0;
        } else {
          HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                          "Postsolve return status is %d\n",
                          (int)postsolve_status);
          model_status_ = HighsModelStatus::POSTSOLVE_ERROR;
          scaled_model_status_ = model_status_;
          hmos_[0].unscaled_model_status_ = model_status_;
          hmos_[0].scaled_model_status_ = model_status_;
          return_status = HighsStatus::Error;
          return returnFromRun(return_status);
        }
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status
      hmos_[original_hmo].unscaled_model_status_ =
          hmos_[solved_hmo].unscaled_model_status_;
      hmos_[original_hmo].scaled_model_status_ =
          hmos_[solved_hmo].scaled_model_status_;
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
    call_status =
        runLpSolver(solved_hmo, "Solving LP without presolve or with basis");
    timer_.stop(timer_.solve_clock);
    this_solve_original_lp_time += timer_.read(timer_.solve_clock);
    return_status =
        interpretCallStatus(call_status, return_status, "runLpSolver");
    if (return_status == HighsStatus::Error)
      return returnFromRun(return_status);
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  //   assert(solved_hmo == original_hmo);
  // solved_hmo will be original_hmo unless the presolved LP is found to be
  // infeasible or unbounded

  if (!getHighsModelStatusAndInfo(solved_hmo)) {
    return_status = HighsStatus::Error;
    return returnFromRun(return_status);
  }

  // Copy HMO solution/basis to HiGHS solution/basis: this resizes solution_ and
  // basis_ The HiGHS solution and basis have to come from the original_hmo for
  // them to have the right dimension.
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

  double lp_solve_final_time = timer_.readRunHighsClock();
  double this_solve_time = lp_solve_final_time - initial_time;
  if (postsolve_iteration_count < 0) {
    HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                      "Postsolve  : \n");
  } else {
    HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                      "Postsolve  : %d\n", postsolve_iteration_count);
  }
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Time       : %8.2f\n", this_solve_time);
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Time Pre   : %8.2f\n", this_presolve_time);
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Time PreLP : %8.2f\n", this_solve_presolved_lp_time);
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                    "Time PostLP: %8.2f\n", this_solve_original_lp_time);
  if (this_solve_time > 0) {
    HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                      "For LP %16s",
                      hmos_[original_hmo].lp_.model_name_.c_str());
    double sum_time = 0;
    if (this_presolve_time > 0) {
      sum_time += this_presolve_time;
      int pct = (100 * this_presolve_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Presolve %8.2f (%3d%%)", this_presolve_time, pct);
    }
    if (this_solve_presolved_lp_time > 0) {
      sum_time += this_solve_presolved_lp_time;
      int pct = (100 * this_solve_presolved_lp_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Solve presolved LP %8.2f (%3d%%)",
                        this_solve_presolved_lp_time, pct);
    }
    if (this_postsolve_time > 0) {
      sum_time += this_postsolve_time;
      int pct = (100 * this_postsolve_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Postsolve %8.2f (%3d%%)", this_postsolve_time, pct);
    }
    if (this_solve_original_lp_time > 0) {
      sum_time += this_solve_original_lp_time;
      int pct = (100 * this_solve_original_lp_time) / this_solve_time;
      HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
                        ": Solve original LP %8.2f (%3d%%)",
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
  return returnFromRun(return_status);
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

HighsStatus Highs::getDualRay(bool& has_dual_ray, double* dual_ray_value) {
  if (!haveHmo("getDualRay")) return HighsStatus::Error;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  return simplex_interface.getDualRay(has_dual_ray, dual_ray_value);
}

HighsStatus Highs::getPrimalRay(bool& has_primal_ray,
                                double* primal_ray_value) {
  underDevelopmentLogMessage("getPrimalRay");
  if (!haveHmo("getPrimalRay")) return HighsStatus::Error;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  return simplex_interface.getPrimalRay(has_primal_ray, primal_ray_value);
}

HighsStatus Highs::getRanging(HighsRanging& ranging) {
  underDevelopmentLogMessage("getRanging");
  if (!haveHmo("getRanging")) return HighsStatus::Error;
  return getHighsRanging(ranging, hmos_[0]);
}

HighsStatus Highs::getBasicVariables(int* basic_variables) {
  if (!haveHmo("getBasicVariables")) return HighsStatus::Error;
  if (basic_variables == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasicVariables: basic_variables is NULL");
    return HighsStatus::Error;
  }
  HighsSimplexInterface simplex_interface(hmos_[0]);
  return simplex_interface.getBasicVariables(basic_variables);
}

HighsStatus Highs::getBasisInverseRow(const int row, double* row_vector,
                                      int* row_num_nz, int* row_indices) {
  if (!haveHmo("getBasisInverseRow")) return HighsStatus::Error;
  if (row_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisInverseRow: row_vector is NULL");
    return HighsStatus::Error;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
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
  if (!haveHmo("getBasisInverseCol")) return HighsStatus::Error;
  if (col_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisInverseCol: col_vector is NULL");
    return HighsStatus::Error;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
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
  if (!haveHmo("getBasisSolve")) return HighsStatus::Error;
  if (Xrhs == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisSolve: Xrhs is NULL");
    return HighsStatus::Error;
  }
  if (solution_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisSolve: solution_vector is NULL");
    return HighsStatus::Error;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
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
  if (!haveHmo("getBasisTransposeSolve")) return HighsStatus::Error;
  if (Xrhs == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisTransposeSolve: Xrhs is NULL");
    return HighsStatus::Error;
  }
  if (solution_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getBasisTransposeSolve: solution_vector is NULL");
    return HighsStatus::Error;
  }
  // solution_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
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
                                 int* row_num_nz, int* row_indices,
                                 const double* pass_basis_inverse_row_vector) {
  if (!haveHmo("getReducedRow")) return HighsStatus::Error;
  if (row_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getReducedRow: row_vector is NULL");
    return HighsStatus::Error;
  }
  // row_indices can be NULL - it's the trigger that determines
  // whether they are identified or not pass_basis_inverse_row_vector
  // NULL - it's the trigger to determine whether it's computed or not
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
  vector<double> basis_inverse_row;
  double* basis_inverse_row_vector = (double*)pass_basis_inverse_row_vector;
  if (basis_inverse_row_vector == NULL) {
    vector<double> rhs;
    vector<int> col_indices;
    rhs.assign(numRow, 0);
    rhs[row] = 1;
    basis_inverse_row.resize(numRow, 0);
    HighsSimplexInterface simplex_interface(hmos_[0]);
    // Form B^{-T}e_{row}
    simplex_interface.basisSolve(rhs, &basis_inverse_row[0], NULL, NULL, true);
    basis_inverse_row_vector = &basis_inverse_row[0];
  }
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    double value = 0;
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      value += lp.Avalue_[el] * basis_inverse_row_vector[row];
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
  if (!haveHmo("getReducedColumn")) return HighsStatus::Error;
  if (col_vector == NULL) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "getReducedColumn: col_vector is NULL");
    return HighsStatus::Error;
  }
  // col_indices can be NULL - it's the trigger that determines
  // whether they are identified or not
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
    return_status = interpretCallStatus(calculateRowValues(lp_, solution_),
                                        return_status, "calculateRowValues");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (solution.row_dual.size() > 0) {
    return_status = interpretCallStatus(calculateColDuals(lp_, solution_),
                                        return_status, "calculateColDuals");
    if (return_status == HighsStatus::Error) return return_status;
  }
  return returnFromHighs(return_status);
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  // Check the user-supplied basis
  if (!isBasisConsistent(lp_, basis)) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "setBasis: invalid basis");
    return HighsStatus::Error;
  }
  // Update the HiGHS basis
  this->basis_ = basis;
  this->basis_.valid_ = true;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
  return HighsStatus::OK;
}

HighsStatus Highs::setBasis() {
  // Invalidate the basis for HiGHS Don't set to logical basis since
  // that causes presolve to be skipped
  this->basis_.valid_ = false;
  // Follow implications of a new HiGHS basis
  newHighsBasis();
  // Can't use returnFromHighs since...
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
  // Check that there is a HighsModelObject
  if (!haveHmo("addRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interpretCallStatus(
      interface.addRows(num_new_row, lower_bounds, upper_bounds, num_new_nz,
                        starts, indices, values),
      return_status, "addRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
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
  if (!haveHmo("addCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interpretCallStatus(
      interface.addCols(num_new_col, costs, lower_bounds, upper_bounds,
                        num_new_nz, starts, indices, values),
      return_status, "addCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeObjectiveSense(const ObjSense sense) {
  HighsStatus return_status = HighsStatus::OK;
  if (!haveHmo("changeObjectiveSense")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  return_status = interpretCallStatus(interface.changeObjectiveSense(sense),
                                      return_status, "changeObjectiveSense");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set,
                           const double* cost) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<int> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(index_collection, cost);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColBounds(const int col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int from_col, const int to_col,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<int> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeRowBounds(const int row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeRowsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<int> local_mask{mask, mask + lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(index_collection, lower, upper);
  return_status =
      interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::changeCoeff(const int row, const int col, const double value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("changeCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCoefficient(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "changeCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getObjectiveSense(ObjSense& sense) {
  if (!haveHmo("getObjectiveSense")) return false;
  sense = hmos_[0].lp_.sense_;
  return true;
}

bool Highs::getCols(const int from_col, const int to_col, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(index_collection, num_col, costs, lower,
                                  upper, num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCols(const int num_set_entries, const int* set, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(index_collection, num_col, costs, lower,
                                  upper, num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCols(const int* mask, int& num_col, double* costs, double* lower,
                    double* upper, int& num_nz, int* start, int* index,
                    double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<int> local_mask{mask, mask + lp_.numCol_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(index_collection, num_col, costs, lower,
                                  upper, num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const int from_row, const int to_row, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(index_collection, num_row, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const int num_set_entries, const int* set, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(index_collection, num_row, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getRows(const int* mask, int& num_row, double* lower, double* upper,
                    int& num_nz, int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local mask that is not const since
  // index_collection.mask_ cannot be const as it changes when
  // deleting rows/columns
  vector<int> local_mask{mask, mask + lp_.numRow_};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &local_mask[0];
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(index_collection, num_row, lower, upper,
                                  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::getCoeff(const int row, const int col, double& value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("getCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);

  call_status = interface.getCoefficient(row, col, value);
  return_status =
      interpretCallStatus(call_status, return_status, "getCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(const int num_set_entries, const int* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteCols(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numCol_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_row;
  index_collection.to_ = to_row;
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(const int num_set_entries, const int* set) {
  if (num_set_entries <= 0) return true;
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Create a local set that is not const since index_collection.set_
  // cannot be const as it may change if the set is not ordered
  vector<int> local_set{set, set + num_set_entries};
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_set_ = true;
  index_collection.set_ = &local_set[0];
  index_collection.set_num_entries_ = num_set_entries;
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::deleteRows(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp_.numRow_;
  index_collection.is_mask_ = true;
  index_collection.mask_ = &mask[0];
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(index_collection);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::scaleCol(const int col, const double scaleval) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("scaleCol")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.scaleCol(col, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleCol");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

bool Highs::scaleRow(const int row, const double scaleval) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (!haveHmo("scaleRow")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.scaleRow(row, scaleval);
  return_status = interpretCallStatus(call_status, return_status, "scaleRow");
  if (return_status == HighsStatus::Error) return false;
  return returnFromHighs(return_status) != HighsStatus::Error;
}

double Highs::getHighsInfinity() { return HIGHS_CONST_INF; }

double Highs::getHighsRunTime() { return timer_.readRunHighsClock(); }

HighsStatus Highs::clearSolver() {
  clearModelStatus();
  clearSolution();
  clearBasis();
  clearInfo();
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void Highs::reportModelStatusSolutionBasis(const std::string message,
                                           const int hmo_ix) {
  HighsModelStatus& model_status = model_status_;
  HighsModelStatus& scaled_model_status = scaled_model_status_;
  HighsSolution& solution = solution_;
  HighsBasis& basis = basis_;
  int unscaled_primal_status = info_.primal_status;
  int scaled_primal_status = unscaled_primal_status;
  int unscaled_dual_status = info_.dual_status;
  int scaled_dual_status = unscaled_dual_status;
  HighsLp& lp = lp_;
  if (hmo_ix >= 0) {
    assert(hmo_ix < (int)hmos_.size());
    model_status = hmos_[hmo_ix].unscaled_model_status_;
    scaled_model_status = hmos_[hmo_ix].scaled_model_status_;
    solution = hmos_[hmo_ix].solution_;
    basis = hmos_[hmo_ix].basis_;
    unscaled_primal_status =
        hmos_[hmo_ix].unscaled_solution_params_.primal_status;
    scaled_primal_status = hmos_[hmo_ix].scaled_solution_params_.primal_status;
    unscaled_dual_status = hmos_[hmo_ix].unscaled_solution_params_.dual_status;
    scaled_dual_status = hmos_[hmo_ix].scaled_solution_params_.dual_status;
    lp = hmos_[hmo_ix].lp_;
  }
  printf(
      "\n%s\nModel status = %s; Scaled model status = %s; LP(%d, %d); solution "
      "([%d:%d] %d, %d; [%d:%d] %d, %d); basis %d "
      "(%d, %d)\n\n",
      message.c_str(), utilHighsModelStatusToString(model_status).c_str(),
      utilHighsModelStatusToString(scaled_model_status).c_str(), lp.numCol_,
      lp.numRow_, unscaled_primal_status, scaled_primal_status,
      (int)solution.col_value.size(), (int)solution.row_value.size(),
      unscaled_dual_status, scaled_dual_status, (int)solution.col_dual.size(),
      (int)solution.row_dual.size(), basis.valid_, (int)basis.col_status.size(),
      (int)basis.row_status.size());
}
#endif

std::string Highs::highsModelStatusToString(
    const HighsModelStatus model_status) const {
  return utilHighsModelStatusToString(model_status);
}

std::string Highs::primalDualStatusToString(const int primal_dual_status) {
  return utilPrimalDualStatusToString(primal_dual_status);
}

// Private methods
HighsPresolveStatus Highs::runPresolve() {
  // Exit if the problem is empty or if presolve is set to off.
  if (options_.presolve == off_string) return HighsPresolveStatus::NotPresolved;
  if (lp_.numCol_ == 0 && lp_.numRow_ == 0)
    return HighsPresolveStatus::NullError;

  // Clear info from previous runs if lp_ has been modified.
  if (presolve_.has_run_) presolve_.clear();
  double start_presolve = timer_.readRunHighsClock();

  // Set time limit.
  if (options_.time_limit > 0 && options_.time_limit < HIGHS_CONST_INF) {
    double left = options_.time_limit - start_presolve;
    if (left <= 0) {
      HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                        "Time limit reached while reading in matrix\n");
      return HighsPresolveStatus::Timeout;
    }

    HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                      "Time limit set: reading matrix took %.2g, presolve "
                      "time left: %.2g\n",
                      start_presolve, left);
    presolve_.options_.time_limit = left;
  }

  // Presolve.
  presolve_.init(lp_, timer_);
  if (options_.time_limit > 0 && options_.time_limit < HIGHS_CONST_INF) {
    double current = timer_.readRunHighsClock();
    double time_init = current - start_presolve;
    double left = presolve_.options_.time_limit - time_init;
    if (left <= 0) {
      HighsPrintMessage(
          options_.output, options_.message_level, ML_VERBOSE,
          "Time limit reached while copying matrix into presolve.\n");
      return HighsPresolveStatus::Timeout;
    }

    HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                      "Time limit set: copying matrix took %.2g, presolve "
                      "time left: %.2g\n",
                      time_init, left);
    presolve_.options_.time_limit = options_.time_limit;
  }

  presolve_.data_.presolve_[0].message_level = options_.message_level;
  presolve_.data_.presolve_[0].output = options_.output;

  HighsPresolveStatus presolve_return_status = presolve_.run();

  // Handle max case.
  if (presolve_return_status == HighsPresolveStatus::Reduced &&
      lp_.sense_ == ObjSense::MAXIMIZE) {
    presolve_.negateReducedLpCost();
    presolve_.data_.reduced_lp_.sense_ = ObjSense::MAXIMIZE;
  }

  // Update reduction counts.
  switch (presolve_.presolve_status_) {
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_.getReducedProblem();
      presolve_.info_.n_cols_removed = lp_.numCol_ - reduced_lp.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_ - reduced_lp.numRow_;
      presolve_.info_.n_nnz_removed =
          (int)lp_.Avalue_.size() - (int)reduced_lp.Avalue_.size();
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      presolve_.info_.n_cols_removed = lp_.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_;
      presolve_.info_.n_nnz_removed = (int)lp_.Avalue_.size();
      break;
    }
    default:
      break;
  }
  return presolve_return_status;
}

HighsPostsolveStatus Highs::runPostsolve() {
  assert(presolve_.has_run_);
  bool solution_ok = isSolutionRightSize(presolve_.getReducedProblem(),
                                         presolve_.data_.reduced_solution_);
  if (!solution_ok) return HighsPostsolveStatus::ReducedSolutionDimenionsError;

  // Run postsolve
  if (presolve_.presolve_status_ != HighsPresolveStatus::Reduced &&
      presolve_.presolve_status_ != HighsPresolveStatus::ReducedToEmpty)
    return HighsPostsolveStatus::NoPostsolve;

  // Handle max case.
  if (lp_.sense_ == ObjSense::MAXIMIZE) presolve_.negateReducedLpColDuals(true);

  HighsPostsolveStatus postsolve_status =
      presolve_.data_.presolve_[0].postsolve(
          presolve_.data_.reduced_solution_, presolve_.data_.reduced_basis_,
          presolve_.data_.recovered_solution_,
          presolve_.data_.recovered_basis_);

  if (postsolve_status != HighsPostsolveStatus::SolutionRecovered)
    return postsolve_status;

  if (lp_.sense_ == ObjSense::MAXIMIZE)
    presolve_.negateReducedLpColDuals(false);

  return HighsPostsolveStatus::SolutionRecovered;
}

// The method below runs calls solveLp to solve the LP associated with
// a particular model, integrating the iteration counts into the
// overall values in HighsInfo
HighsStatus Highs::runLpSolver(const int model_index, const string message) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;

  // Check that the model index is OK
  bool model_index_ok = model_index >= 0 && model_index < (int)hmos_.size();
  assert(model_index_ok);
  if (!model_index_ok) return HighsStatus::Error;

  HighsModelObject& model = hmos_[model_index];

  // Transfer the LP solver iteration counts to this model
  HighsIterationCounts& iteration_counts = hmos_[model_index].iteration_counts_;
  copyHighsIterationCounts(info_, iteration_counts);

  // Solve the LP
  call_status = solveLp(model, message);
  return_status = interpretCallStatus(call_status, return_status, "solveLp");
  if (return_status == HighsStatus::Error) return return_status;

  // Transfer this model's LP solver iteration counts to HiGHS
  copyHighsIterationCounts(iteration_counts, info_);

  return returnFromHighs(return_status);
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

  //  std::cout << "warning: Feature under development" << std::endl;

  writeSolutionToFile(file, lp, basis, solution, pretty);

  //  return HighsStatus::Warning;
  return HighsStatus::OK;
}

// Actions to take if there is a new Highs basis
void Highs::newHighsBasis() {
  if (hmos_.size() > 0) {
    // Copy this basis to the HMO basis and clear any simplex basis
    hmos_[0].basis_ = this->basis_;
    HighsSimplexInterface interface(hmos_[0]);
    interface.clearBasis();
  }
}

// Ensure that the HiGHS solution and basis have the same size as the
// model, and that the HiGHS basis is kept up-to-date with any solved
// basis
void Highs::forceHighsSolutionBasisSize() {
  // Ensure that the HiGHS solution vectors are the right size
  solution_.col_value.resize(lp_.numCol_);
  solution_.row_value.resize(lp_.numRow_);
  solution_.col_dual.resize(lp_.numCol_);
  solution_.row_dual.resize(lp_.numRow_);
  // Ensure that the HiGHS basis vectors are the right size,
  // invalidating the basis if they aren't
  if ((int)basis_.col_status.size() != lp_.numCol_) {
    basis_.col_status.resize(lp_.numCol_);
    basis_.valid_ = false;
  }
  if ((int)basis_.row_status.size() != lp_.numRow_) {
    basis_.row_status.resize(lp_.numRow_);
    basis_.valid_ = false;
  }
}

bool Highs::getHighsModelStatusAndInfo(const int solved_hmo) {
  if (!haveHmo("getHighsModelStatusAndInfo")) return false;

  model_status_ = hmos_[solved_hmo].unscaled_model_status_;
  scaled_model_status_ = hmos_[solved_hmo].scaled_model_status_;

  HighsSolutionParams& solution_params =
      hmos_[solved_hmo].unscaled_solution_params_;

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

bool Highs::haveHmo(const string method_name) const {
  bool have_hmo = hmos_.size() > 0;
#ifdef HiGHSDEV
  if (!have_hmo)
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Method %s called without any HighsModelObject",
                    method_name.c_str());
#endif
  assert(have_hmo);
  return have_hmo;
}

void Highs::clearModelStatus() {
  model_status_ = HighsModelStatus::NOTSET;
  scaled_model_status_ = HighsModelStatus::NOTSET;
}

void Highs::clearSolution() {
  info_.primal_status = (int)PrimalDualStatus::STATUS_NOTSET;
  info_.dual_status = (int)PrimalDualStatus::STATUS_NOTSET;
  clearSolutionUtil(solution_);
}

void Highs::clearBasis() { clearBasisUtil(basis_); }

void Highs::clearInfo() { info_.clear(); }

// Applies checks before returning from run()
HighsStatus Highs::returnFromRun(const HighsStatus run_return_status) {
  bool have_solution = false;
  HighsStatus return_status = run_return_status;
  if (hmos_.size() == 0) {
    // No model has been loaded: ensure that the status, solution,
    // basis and info associated with any previous model are cleared
    clearSolver();
    return returnFromHighs(return_status);
  } else {
    // A model has been loaded: remove any additional HMO created when solving
    if (hmos_.size() > 1) hmos_.pop_back();
    // There should be only one entry in hmos_
    assert((int)hmos_.size() == 1);
    // Make sure that the unscaled status, solution, basis and info
    // are consistent with the scaled status
#ifdef HiGHSDEV
    reportModelStatusSolutionBasis("returnFromRun(HiGHS)");
    reportModelStatusSolutionBasis("returnFromRun(HMO_0)", 0);
#endif
    switch (scaled_model_status_) {
      // First consider the error returns
      case HighsModelStatus::NOTSET:
      case HighsModelStatus::LOAD_ERROR:
      case HighsModelStatus::MODEL_ERROR:
      case HighsModelStatus::PRESOLVE_ERROR:
      case HighsModelStatus::SOLVE_ERROR:
      case HighsModelStatus::POSTSOLVE_ERROR:
        clearSolver();
        assert(return_status == HighsStatus::Error);
        break;

      // Then consider the OK returns
      case HighsModelStatus::MODEL_EMPTY:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_UNBOUNDED:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::OPTIMAL:
        have_solution = true;
        // The following is an aspiration
        //        assert(info_.primal_status ==
        //                   (int)PrimalDualStatus::STATUS_FEASIBLE_POINT);
        //        assert(info_.dual_status ==
        //                   (int)PrimalDualStatus::STATUS_FEASIBLE_POINT);
        assert(model_status_ == HighsModelStatus::NOTSET ||
               model_status_ == HighsModelStatus::OPTIMAL);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::OK);
        break;

      // Finally consider the warning returns
      case HighsModelStatus::REACHED_TIME_LIMIT:
      case HighsModelStatus::REACHED_ITERATION_LIMIT:
        clearSolution();
        clearBasis();
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::Warning);
        break;
      case HighsModelStatus::DUAL_INFEASIBLE:
        clearSolution();
        // May have a basis, according to whether infeasibility was
        // detected in presolve or solve
        clearInfo();
        assert(model_status_ == scaled_model_status_);
        assert(return_status == HighsStatus::Warning);
        break;
    }
  }
  if (have_solution) debugSolutionRightSize(options_, lp_, solution_);
  bool have_basis = basis_.valid_;
  if (have_basis) {
    if (debugBasisRightSize(options_, lp_, basis_) ==
        HighsDebugStatus::LOGICAL_ERROR)
      return_status = HighsStatus::Error;
  }

  if (have_solution && have_basis) {
    if (debugHighsBasicSolution("Return from run()", options_, lp_, basis_,
                                solution_, info_, model_status_) ==
        HighsDebugStatus::LOGICAL_ERROR)
      return_status = HighsStatus::Error;
  }
  return returnFromHighs(return_status);
}

// Applies checks before returning from HiGHS
HighsStatus Highs::returnFromHighs(HighsStatus highs_return_status) {
  HighsStatus return_status = highs_return_status;

  forceHighsSolutionBasisSize();

  const bool consistent = debugBasisConsistent(options_, lp_, basis_) !=
                          HighsDebugStatus::LOGICAL_ERROR;
  if (!consistent) {
    HighsLogMessage(
        options_.logfile, HighsMessageType::ERROR,
        "returnFromHighs: Supposed to be a HiGHS basis, but not consistent");
    assert(consistent);
    return_status = HighsStatus::Error;
  }

  if (hmos_.size()) {
    const bool simplex_lp_ok =
        debugSimplexLp(hmos_[0]) != HighsDebugStatus::LOGICAL_ERROR;
    if (!simplex_lp_ok) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "returnFromHighs: Simplex LP not OK");
      assert(simplex_lp_ok);
      return_status = HighsStatus::Error;
    }
  }

  return return_status;
}
void Highs::underDevelopmentLogMessage(const string method_name) {
  HighsLogMessage(
      options_.logfile, HighsMessageType::WARNING,
      "Method %s is still under development and behaviour may be unpredictable",
      method_name.c_str());
}

void Highs::getPresolveReductionCounts(int& rows, int& cols, int& nnz) const {
  rows = presolve_.info_.n_rows_removed;
  cols = presolve_.info_.n_cols_removed;
  nnz = presolve_.info_.n_nnz_removed;
}
