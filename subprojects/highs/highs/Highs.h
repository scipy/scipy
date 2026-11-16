/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file Highs.h
 * @brief The HiGHS class
 */
#ifndef HIGHS_H_
#define HIGHS_H_

#include <sstream>

#include "lp_data/HighsCallback.h"
#include "lp_data/HighsIis.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsRanging.h"
#include "lp_data/HighsRunData.h"
#include "lp_data/HighsSolutionDebug.h"
#include "model/HighsModel.h"
#include "presolve/ICrash.h"
#include "presolve/PresolveComponent.h"

/**
 * @brief Return the version
 */
const char* highsVersion();

/**
 * @brief Return detailed version information, githash and compilation
 * date
 */
HighsInt highsVersionMajor();
HighsInt highsVersionMinor();
HighsInt highsVersionPatch();
const char* highsGithash();

/**
 * @brief Class to set parameters and run HiGHS
 */
class Highs {
 public:
  Highs();
  virtual ~Highs() { this->closeLogFile(); }

  /**
   * @brief Return the version as a string
   */
  std::string version() const { return highsVersion(); }

  /**
   * @brief Return major version
   */
  HighsInt versionMajor() const { return highsVersionMajor(); }

  /**
   * @brief Return minor version
   */
  HighsInt versionMinor() const { return highsVersionMinor(); }

  /**
   * @brief Return patch version
   */
  HighsInt versionPatch() const { return highsVersionPatch(); }

  /**
   * @brief Return githash
   */
  std::string githash() const { return highsGithash(); }

  /**
   * @brief Return third-party notice and attribution details
   */
  std::string getThirdPartyNotice() const;

  /**
   * @brief Reset the options and then call clearModel()
   */
  HighsStatus clear();

  /**
   * @brief Clear the incumbent model and then call clearSolver()
   */
  HighsStatus clearModel();

  /**
   * @brief Clear all solution data associated with the model
   */
  HighsStatus clearSolver();

  /**
   * @brief Clear all dual data associated with the model
   */
  HighsStatus clearSolverDualData();

  /**
   * @brief Release all retained memory back to the allocator
   *
   * Clears all solver state and shrinks internal vectors to free
   * unused capacity. Useful in long-running services that reuse a
   * Highs instance across multiple solves to prevent unbounded RSS
   * growth from heap fragmentation.
   */
  HighsStatus releaseMemory();

  /**
   * Methods for model input
   */

  /**
   * Every model loading module eventually uses
   * passModel(HighsModel model) to communicate the model to HiGHS.
   */

  /**
   * @brief Pass a HighsModel instance to Highs
   */
  HighsStatus passModel(HighsModel model);

  /**
   * @brief Pass a HighsLp instance to Highs
   */
  HighsStatus passModel(HighsLp lp);

  /**
   * @brief Pass a QP (possibly with integrality data) via pointers to vectors
   * of data
   */
  HighsStatus passModel(
      const HighsInt num_col, const HighsInt num_row, const HighsInt num_nz,
      const HighsInt q_num_nz, const HighsInt a_format, const HighsInt q_format,
      const HighsInt sense, const double offset, const double* col_cost,
      const double* col_lower, const double* col_upper, const double* row_lower,
      const double* row_upper, const HighsInt* a_start, const HighsInt* a_index,
      const double* a_value, const HighsInt* q_start, const HighsInt* q_index,
      const double* q_value, const HighsInt* integrality = nullptr);

  /**
   * @brief Pass an LP (possibly with integrality data) via pointers to vectors
   * of data
   */
  HighsStatus passModel(const HighsInt num_col, const HighsInt num_row,
                        const HighsInt num_nz, const HighsInt a_format,
                        const HighsInt sense, const double offset,
                        const double* col_cost, const double* col_lower,
                        const double* col_upper, const double* row_lower,
                        const double* row_upper, const HighsInt* a_start,
                        const HighsInt* a_index, const double* a_value,
                        const HighsInt* integrality = nullptr);

  /**
   * @brief Pass a HighsHessian instance for the incumbent model
   */
  HighsStatus passHessian(HighsHessian hessian_);

  /**
   * @brief Pass the Hessian for the incumbent model via pointers to vectors of
   * data
   */
  HighsStatus passHessian(const HighsInt dim, const HighsInt num_nz,
                          const HighsInt format, const HighsInt* start,
                          const HighsInt* index, const double* value);
  /**
   * @brief Pass a Hessian oracle for the incumbent model
   */
  HighsStatus passHessian(const HighsInt dim,
                          HighsHessianFunctionType oracleCall,
                          void* oracle_data,
                          HighsCHessianFunctionType c_oracleCall = nullptr);

  /**
   * @brief Checks any incumbent Hessian oracle
   */
  HighsStatus checkHessianOracle(const bool exit_on_first_error = false) const;

  /**
   * @brief Pass multiple linear objectives for the incumbent model
   */
  HighsStatus passLinearObjectives(
      const HighsInt num_linear_objective,
      const HighsLinearObjective* linear_objective);

  /**
   * @brief Add a linear objective for the incumbent model
   */
  HighsStatus addLinearObjective(const HighsLinearObjective& linear_objective,
                                 const HighsInt iObj = -1);

  /**
   * @brief Get number of linear objectives from the incumbent model
   */
  HighsInt getNumLinearObjectives() const {
    return static_cast<HighsInt>(multi_linear_objective_.size());
  }

  /**
   * @brief Get a linear objective from the incumbent model
   */
  const HighsLinearObjective& getLinearObjective(const HighsInt idx) const {
    assert(idx >= 0 &&
           idx < static_cast<HighsInt>(multi_linear_objective_.size()));
    return multi_linear_objective_[idx];
  }

  /**
   * @brief Clear the multiple linear objective data
   */
  HighsStatus clearLinearObjectives();

  /**
   * @brief Pass a column name to the incumbent model
   */
  HighsStatus passColName(const HighsInt col, const std::string& name);

  /**
   * @brief Pass a row name to the incumbent model
   */
  HighsStatus passRowName(const HighsInt row, const std::string& name);

  /**
   * @brief Pass a model name to the incumbent model
   */
  HighsStatus passModelName(const std::string& name);

  /**
   * @brief Read in a model
   */
  HighsStatus readModel(const std::string& filename);

  /**
   * @brief Read in a basis
   */
  HighsStatus readBasis(const std::string& filename);

  /**
   * @brief Generate a PBM image of the matrix nonzeros, and possibly
   * Hessian nonzeros. Note that the .pbm extnesion will be added to
   * the fime names passed
   */
  HighsStatus matrixImage(const std::string& matrix_image_filename,
                          const std::string& hessian_image_filename = "") const;

  /**
   * @brief Presolve the incumbent model, allowing the presolved model
   * to be extracted. Subsequent solution of the incumbent model will
   * only use presolve if there is no valid basis
   */
  HighsStatus presolve();

  /**
   * @brief Run the solver, applying file-based options and user
   * scaling before optimization
   */
  HighsStatus run();

  /**
   * @brief Postsolve the incumbent model using a solution
   */
  HighsStatus postsolve(const HighsSolution& solution);

  /**
   * @brief Postsolve the incumbent model using a solution and basis
   */
  HighsStatus postsolve(const HighsSolution& solution, const HighsBasis& basis);

  /**
   * @brief Write the current solution to a file in a given style
   */
  HighsStatus writeSolution(const std::string& filename,
                            const HighsInt style = kSolutionStyleRaw);

  /**
   * @brief Read a HiGHS solution file in a given style
   */
  HighsStatus readSolution(const std::string& filename,
                           const HighsInt style = kSolutionStyleRaw);

  /**
   * @brief Assess the validity, integrality and feasibility of the
   * current primal solution. Row values are computed and checked
   * against what's in Highs::solution_.row_value and, if the
   * differences exceed a tolerance, valid returns false.  If any of
   * valid, integral or feasible is false, then assessPrimalSolution
   * returns HighsStatus::kWarning.
   */
  HighsStatus assessPrimalSolution(bool& valid, bool& integral,
                                   bool& feasible) const;

  /**
   * Methods for HiGHS option input/output
   */

  /**
   * @brief Set an option to the bool/HighsInt/double/string value if it's
   * legal and, for bool/HighsInt/double, only if it's of the correct type
   */

  HighsStatus setOptionValue(const std::string& option, const bool value);

  HighsStatus setOptionValue(const std::string& option, const HighsInt value);

#ifdef HIGHSINT64
  HighsStatus setOptionValue(const std::string& option, const int value) {
    return setOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setOptionValue(const std::string& option, const double value);

  HighsStatus setOptionValue(const std::string& option,
                             const std::string& value);

  HighsStatus setOptionValue(const std::string& option, const char* value);

  /**
   * @brief Read option values from a file
   */
  HighsStatus readOptions(const std::string& filename);

  /**
   * @brief Pass a HighsOptions instance to Highs
   */
  HighsStatus passOptions(const HighsOptions& options);

  /**
   * @brief Get a const reference to the internal option values
   */
  const HighsOptions& getOptions() const { return options_; }

  /**
   * @brief Gets an option value as bool/HighsInt/double/string and, for
   * bool/int/double, only if it's of the correct type.
   *
   * NB Deprecate in v2.0, in order to replace with more general
   * get*OptionValues
   */
  HighsStatus getOptionValue(const std::string& option, bool& value) const {
    return this->getBoolOptionValues(option, &value);
  }

  HighsStatus getOptionValue(const std::string& option, HighsInt& value) const {
    return this->getIntOptionValues(option, &value);
  }

  HighsStatus getOptionValue(const std::string& option, double& value) const {
    return this->getDoubleOptionValues(option, &value);
  }

  HighsStatus getOptionValue(const std::string& option,
                             std::string& value) const {
    return this->getStringOptionValues(option, &value);
  }

  /**
   * @brief Get the type expected by an option
   */
  HighsStatus getOptionType(const std::string& option,
                            HighsOptionType& type) const {
    return this->getOptionType(option, &type);
  }

  /**
   * @brief Reset the options to the default values
   */
  HighsStatus resetOptions();

  /**
   * @brief Write (deviations from default values of) the options to a
   * file, using the standard format used to read options from a file.
   * Possible to write only deviations from default values.
   */
  HighsStatus writeOptions(const std::string& filename,  //!< The filename
                           const bool report_only_deviations = false) const;

  /**
   * @brief Returns the number of user-settable options
   */
  HighsInt getNumOptions() const {
    return this->options_.num_user_settable_options_;
  }

  /**
   * @brief Get the number of user-settable options
   */
  HighsStatus getOptionName(const HighsInt index, std::string* name) const;

  /**
   * @brief Get the type of an option
   */
  HighsStatus getOptionType(const std::string& option,
                            HighsOptionType* type) const;

  /**
   * @brief Get the current and default values of a bool option
   */
  HighsStatus getBoolOptionValues(const std::string& option,
                                  bool* current_value = nullptr,
                                  bool* default_value = nullptr) const;

  /**
   * @brief Get the current, min, max and default values of an int option
   */
  HighsStatus getIntOptionValues(const std::string& option,
                                 HighsInt* current_value = nullptr,
                                 HighsInt* min_value = nullptr,
                                 HighsInt* max_value = nullptr,
                                 HighsInt* default_value = nullptr) const;

  /**
   * @brief Get the current, min, max and default values of a double option
   */
  HighsStatus getDoubleOptionValues(const std::string& option,
                                    double* current_value = nullptr,
                                    double* min_value = nullptr,
                                    double* max_value = nullptr,
                                    double* default_value = nullptr) const;

  /**
   * @brief Get the current and default values of a string option
   */
  HighsStatus getStringOptionValues(const std::string& option,
                                    std::string* current_value = nullptr,
                                    std::string* default_value = nullptr) const;

  /**
   * @brief Get a const reference to the internal run data values.
   */
  const HighsRunData& getRunData() const { return run_data_; }

  /**
   * @brief Get an run_data value as HighsInt/int64_t/double, and only if
   * it's of the correct type.
   */

  HighsStatus getRunDataValue(const std::string& run_data,
                              HighsInt& value) const;

#ifndef HIGHSINT64
  HighsStatus getRunDataValue(const std::string& run_data,
                              int64_t& value) const;
#endif

  HighsStatus getRunDataValue(const std::string& run_data, double& value) const;

  HighsStatus getRunDataType(const std::string& run_data,
                             HighsRunDataType& type) const;

  /**
   * @brief Write run data values to a file, with the extension ".html"
   * producing HTML, otherwise using the standard format used to read
   * options from a file.
   */
  HighsStatus writeRunData(const std::string& filename = "") const;

  /**
   * @brief Get a const reference to the internal info values.
   */
  const HighsInfo& getInfo() const { return info_; }

  /**
   * @brief Get an info value as HighsInt/int64_t/double, and only if
   * it's of the correct type.
   */

  HighsStatus getInfoValue(const std::string& info, HighsInt& value) const;

#ifndef HIGHSINT64
  HighsStatus getInfoValue(const std::string& info, int64_t& value) const;
#endif

  HighsStatus getInfoValue(const std::string& info, double& value) const;

  HighsStatus getInfoType(const std::string& info, HighsInfoType& type) const;

  /**
   * @brief Write info values to a file, with the extension ".html"
   * producing HTML, otherwise using the standard format used to read
   * options from a file.
   */
  HighsStatus writeInfo(const std::string& filename = "") const;

  /**
   * @brief Get the value of infinity used by HiGHS
   */
  double getInfinity() const { return kHighsInf; }

  /**
   * @brief Get the size of HighsInt
   */
  HighsInt getSizeofHighsInt() const {
    return sizeof(options_.num_user_settable_options_);
  }

  /**
   * @brief Get the run time of HiGHS
   */
  double getRunTime() const { return timer_.read(); }

  /**
   * Methods for model output
   */

  /**
   * @brief Identify the standard form of the HighsLp instance in
   * HiGHS
   */
  HighsStatus getStandardFormLp(HighsInt& num_col, HighsInt& num_row,
                                HighsInt& num_nz, double& offset,
                                double* cost = nullptr, double* rhs = nullptr,
                                HighsInt* start = nullptr,
                                HighsInt* index = nullptr,
                                double* value = nullptr);

  /**
   * @brief Return a const reference to the internal presolved HighsLp
   * instance
   */
  const HighsLp& getPresolvedLp() const { return presolved_model_.lp_; }

  /**
   * @brief Return a const reference to the internal presolved
   * HighsModel instance
   */
  const HighsModel& getPresolvedModel() const { return presolved_model_; }

  /**
   * @brief Return a const reference to the logging data for presolve
   */
  const HighsPresolveLog& getPresolveLog() const { return presolve_log_; }

  /**
   * @brief Return a const pointer to the original column indices for
   * the presolved model
   */
  const HighsInt* getPresolveOrigColsIndex() const {
    return presolve_.data_.postSolveStack.getOrigColIndex().data();
  }

  /**
   * @brief Return a const pointer to the original row indices for the
   * presolved model
   */
  const HighsInt* getPresolveOrigRowsIndex() const {
    return presolve_.data_.postSolveStack.getOrigRowIndex().data();
  }

  /**
   * @brief Return an LP associated with a MIP and its solution, with
   * each integer variable fixed to the value it takes in the MIP
   * solution. If no solution is available, an error is returned.
   */
  HighsStatus getFixedLp(HighsLp& lp) const;

  /**
   * @brief Return a const reference to the incumbent LP
   */
  const HighsLp& getLp() const { return model_.lp_; }

  /**
   * @brief Return a const reference to the incumbent model
   */
  const HighsModel& getModel() const { return model_; }

  /**
   * @brief Return a const reference to the internal HighsSolution
   * instance
   */
  const HighsSolution& getSolution() const { return solution_; }

  /**
   * @brief Return a const reference to the internal IIS LP instance
   */
  const HighsLp& getIisLp() const { return iis_.model_.lp_; }

  /**
   * @brief Zero all clocks in the internal HighsTimer instance
   */
  void zeroAllClocks() { timer_.zeroAllClocks(); }

  /**
   * @brief Return a const reference to the internal HighsSolution
   * instance
   */
  const std::vector<HighsObjectiveSolution>& getSavedMipSolutions() const {
    return saved_objective_and_solution_;
  }

  /**
   * @brief Return a const reference to the internal ICrash info instance
   */
  const ICrashInfo& getICrashInfo() const { return icrash_info_; }

  /**
   * @brief Return a const reference to the internal HighsBasis instance
   */
  const HighsBasis& getBasis() const { return basis_; }

  /**
   * @brief Return the status for the incumbent model.
   */
  const HighsModelStatus& getModelStatus() const { return model_status_; }

  /**
   * @brief Returns the current model's presolve status
   */
  const HighsPresolveStatus& getModelPresolveStatus() const {
    return model_presolve_status_;
  }

  /**
   * @brief Indicate whether a dual unbounded ray exists, and (at the
   * expense of solving an LP) gets it if it does not and
   * dual_ray_value is not nullptr
   */
  HighsStatus getDualRay(bool& has_dual_ray, double* dual_ray_value = nullptr);

  /**
   * @brief Indicate whether a dual unbounded ray exists, and gets it
   * if it does
   */
  HighsStatus getDualRaySparse(bool& has_dual_ray, HVector& row_ep_buffer);

  /**
   * @brief Indicate whether a dual unboundedness direction exists,
   * and (at the expense of solving an LP) gets it if
   * dual_unboundedness_direction is not nullptr
   */
  HighsStatus getDualUnboundednessDirection(
      bool& has_dual_unboundedness_direction,
      double* dual_unboundedness_direction_value = nullptr);

  /**
   * @brief Indicate whether a primal unbounded ray exists, and (at
   * the expense of solving an LP) gets it if primal_ray is not
   * nullptr
   */
  HighsStatus getPrimalRay(bool& has_primal_ray,
                           double* primal_ray_value = nullptr);

  /**
   * @brief Get the ranging information for the current LP
   */
  HighsStatus getRanging(HighsRanging& ranging);

  /**
   * @brief Solve the feasibility relaxation problem
   */
  HighsStatus feasibilityRelaxation(const double global_lower_penalty,
                                    const double global_upper_penalty,
                                    const double global_rhs_penalty,
                                    const double* local_lower_penalty = nullptr,
                                    const double* local_upper_penalty = nullptr,
                                    const double* local_rhs_penalty = nullptr);

  /**
   * @brief Get the ill-conditioning information for the current basis
   */
  HighsStatus getIllConditioning(HighsIllConditioning& ill_conditioning,
                                 const bool constraint,
                                 const HighsInt method = 0,
                                 const double ill_conditioning_bound = 1e-4);

  /**
   * @brief Get the suggested objective and bound scaling for the incumbent
   * model
   */
  HighsStatus getObjectiveBoundScaling(HighsInt& suggested_objective_scale,
                                       HighsInt& suggested_bound_scale);

  /**
   * @brief Get (any) irreducible infeasible subsystem (IIS)
   * information for the incumbent model
   */
  HighsStatus getIis(HighsIis& iis);

  /**
   * @brief Get the current model objective function value
   */
  double getObjectiveValue() const { return info_.objective_function_value; }

  /**
   * @brief Try to get the current dual objective function value
   */
  HighsStatus getDualObjectiveValue(double& dual_objective_value) const;

  /**
   * Methods for operations with the invertible representation of the
   * current basis matrix
   */

  /**
   * @brief Returns true if an invertible representation of the
   * current basis matrix is available
   */
  bool hasInvert() const { return ekk_instance_.status_.has_invert; }

  /**
   * @brief Gets the basic variables in the order corresponding to
   * calls to getBasisInverseRow, getBasisInverseCol, getBasisSolve,
   * getBasisTransposeSolve, getReducedRow and
   * getReducedColumn. Non-negative entries are indices of columns,
   * and negative entries are -(row_index+1).
   */
  HighsStatus getBasicVariables(HighsInt* basic_variables);

  /**
   * @brief Form a row of \f$B^{-1}\f$ for basis matrix \f$B\f$,
   * returning the indices of the nonzeros unless row_num_nz is
   * nullptr
   */
  HighsStatus getBasisInverseRow(const HighsInt row, double* row_vector,
                                 HighsInt* row_num_nz = nullptr,
                                 HighsInt* row_indices = nullptr);

  /**
   * @brief Form a column of \f$B^{-1}\f$ for basis matrix \f$B\f$,
   * returning the indices of the nonzeros unless col_num_nz is
   * nullptr
   */
  HighsStatus getBasisInverseCol(const HighsInt col, double* col_vector,
                                 HighsInt* col_num_nz = nullptr,
                                 HighsInt* col_indices = nullptr);

  /**
   * @brief Form \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$, returning the indices of the nonzeros unless
   * solution_num_nz is nullptr
   */
  HighsStatus getBasisSolve(const double* rhs, double* solution_vector,
                            HighsInt* solution_num_nz = nullptr,
                            HighsInt* solution_indices = nullptr);

  /**
   * @brief Form \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$, returning the indices of the nonzeros unless
   * solution_num_nz is nullptr
   */
  HighsStatus getBasisTransposeSolve(const double* rhs, double* solution_vector,
                                     HighsInt* solution_num_nz = nullptr,
                                     HighsInt* solution_indices = nullptr);

  /**
   * @brief Form a row of \f$B^{-1}A\f$, returning the indices of the
   * nonzeros unless row_num_nz is nullptr, computing the row using
   * pass_basis_inverse_row_vector unless it is nullptr
   */
  HighsStatus getReducedRow(
      const HighsInt row, double* row_vector, HighsInt* row_num_nz = nullptr,
      HighsInt* row_indices = nullptr,
      const double* pass_basis_inverse_row_vector = nullptr);

  /**
   * @brief Form a column of \f$B^{-1}A\f$, returning the indices of
   * the nonzeros unless col_num_nz is nullptr
   */
  HighsStatus getReducedColumn(const HighsInt col, double* col_vector,
                               HighsInt* col_num_nz = nullptr,
                               HighsInt* col_indices = nullptr);

  /**
   * @brief Get the condition number of the current basis matrix,
   * possibly computing it exactly and reporting the error in the
   * approximate condition number
   */
  HighsStatus getKappa(double& kappa, const bool exact = false,
                       const bool report = false) const;

  /**
   * @brief Get the number of columns in the incumbent model
   */
  HighsInt getNumCol() const { return model_.lp_.num_col_; }

  /**
   * @brief Get the number of rows in the incumbent model
   */
  HighsInt getNumRow() const { return model_.lp_.num_row_; }

  /**
   * @brief Get the number of (constraint matrix) nonzeros in the incumbent
   * model
   */
  HighsInt getNumNz() const { return model_.lp_.numNz(); }

  /**
   * @brief Get the number of Hessian matrix nonzeros in the incumbent model
   */
  HighsInt getHessianNumNz() const { return model_.hessian_.numNz(); }

  /**
   * @brief Get the objective sense of the incumbent model
   */
  HighsStatus getObjectiveSense(ObjSense& sense) const;

  /**
   * @brief Get the objective offset of the incumbent model
   */
  HighsStatus getObjectiveOffset(double& offset) const;

  /**
   * @brief Get multiple columns from the model given by an interval [from_col,
   * to_col]
   */
  HighsStatus getCols(
      const HighsInt
          from_col,  //!< The index of the first column to get from the model
      const HighsInt
          to_col,  //!< The index of the last column to get from the model
      HighsInt& num_col,  //!< Number of columns got from the model
      double* cost,       //!< Array of size num_col with costs
      double* lower,      //!< Array of size num_col with lower bounds
      double* upper,      //!< Array of size num_col with upper bounds
      HighsInt& num_nz,   //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!< Array of size num_col with start indices of the columns
      HighsInt*
          index,     //!< Array of size num_nz with row indices for the columns
      double* value  //!< Array of size num_nz with row values for the columns
  ) const;

  /**
   * @brief Get multiple columns from the model given by a set
   */
  HighsStatus getCols(
      const HighsInt num_set_entries,  //!< The number of indices in the set
      const HighsInt* set,  //!< Array of size num_set_entries with indices of
                            //!< columns to get
      HighsInt& num_col,    //!< Number of columns got from the model
      double* cost,         //!< Array of size num_col with costs
      double* lower,        //!< Array of size num_col with lower bounds
      double* upper,        //!< Array of size num_col with upper bounds
      HighsInt& num_nz,     //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!< Array of size num_col with start indices of the columns
      HighsInt*
          index,     //!< Array of size num_nz with row indices for the columns
      double* value  //!< Array of size num_nz with row values for the columns
  ) const;

  /**
   * @brief Get multiple columns from the model given by a mask
   */
  HighsStatus getCols(
      const HighsInt* mask,  //!< Full length array with 1 => get; 0 => not
      HighsInt& num_col,     //!< Number of columns got from the model
      double* cost,          //!< Array of size num_col with cost
      double* lower,         //!< Array of size num_col with lower bounds
      double* upper,         //!< Array of size num_col with upper bounds
      HighsInt& num_nz,      //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!<  Array of size num_col with start indices of the columns
      HighsInt*
          index,     //!<  Array of size num_nz with row indices for the columns
      double* value  //!<  Array of size num_nz with row values for the columns
  ) const;

  /**
   * @brief Get a column name from the incumbent model
   */
  HighsStatus getColName(const HighsInt col, std::string& name) const;

  /**
   * @brief Get column index corresponding to name
   */
  HighsStatus getColByName(const std::string& name, HighsInt& col);

  /**
   * @brief Get a column integrality from the incumbent model
   */
  HighsStatus getColIntegrality(const HighsInt col,
                                HighsVarType& integrality) const;

  /**
   * @brief Get multiple rows from the model given by an interval [from_row,
   * to_row]
   */
  HighsStatus getRows(
      const HighsInt
          from_row,  //!< The index of the first row to get from the model
      const HighsInt
          to_row,         //!< The index of the last row to get from the model
      HighsInt& num_row,  //!< Number of rows got from the model
      double* lower,      //!< Array of size num_row with lower bounds
      double* upper,      //!< Array of size num_row with upper bounds
      HighsInt& num_nz,   //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!< Array of size num_row with start indices of the rows
      HighsInt*
          index,     //!< Array of size num_nz with column indices for the rows
      double* value  //!< Array of size num_nz with column values for the rows
  ) const;

  /**
   * @brief Get multiple rows from the model given by a set
   */
  HighsStatus getRows(
      const HighsInt num_set_entries,  //!< The number of indices in the set
      const HighsInt*
          set,  //!< Array of size num_set_entries with indices of rows to get
      HighsInt& num_row,  //!< Number of rows got from the model
      double* lower,      //!< Array of size num_row with lower bounds
      double* upper,      //!< Array of size num_row with upper bounds
      HighsInt& num_nz,   //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!< Array of size num_row with start indices of the rows
      HighsInt*
          index,     //!< Array of size num_nz with column indices for the rows
      double* value  //!< Array of size num_nz with column values for the rows
  ) const;

  /**
   * @brief Get multiple rows from the model given by a mask
   */
  HighsStatus getRows(
      const HighsInt* mask,  //!< Full length array with 1 => get; 0 => not
      HighsInt& num_row,     //!< Number of rows got from the model
      double* lower,         //!< Array of size num_row with lower bounds
      double* upper,         //!< Array of size num_row with upper bounds
      HighsInt& num_nz,      //!< Number of nonzeros got from the model
      HighsInt*
          start,  //!< Array of size num_row with start indices of the rows
      HighsInt*
          index,     //!< Array of size num_nz with column indices for the rows
      double* value  //!< Array of size num_nz with column values for the rows
  ) const;

  /**
   * @brief Get a row name from the incumbent model
   */
  HighsStatus getRowName(const HighsInt row, std::string& name) const;

  /**
   * @brief Get row index corresponding to name
   */
  HighsStatus getRowByName(const std::string& name, HighsInt& row);

  /**
   * @brief Get a matrix coefficient
   */
  HighsStatus getCoeff(const HighsInt row, const HighsInt col,
                       double& value) const;

  /**
   * @brief Write out the incumbent model to a file
   */
  HighsStatus writeModel(const std::string& filename = "");

  /**
   * @brief Write out the incumbent presolved model to a file
   */
  HighsStatus writePresolvedModel(const std::string& filename = "");

  /**
   * @brief Write out the internal IIS LP instance to a file
   */
  HighsStatus writeIisModel(const std::string& filename = "");

  /**
   * @brief Write out the given model to a file
   */
  HighsStatus writeLocalModel(HighsModel& model,
                              const std::string& filename = "");

  /**
   * @brief Write out the internal HighsBasis instance to a file
   */
  HighsStatus writeBasis(const std::string& filename = "");

  /**
   * Methods for incumbent model modification
   */

  /**
   * @brief Change the objective sense of the incumbent model
   */
  HighsStatus changeObjectiveSense(const ObjSense sense);

  /**
   * @brief Change the objective offset of the incumbent model
   */
  HighsStatus changeObjectiveOffset(const double offset);

  /**
   * @brief Change the integrality of a column
   */
  HighsStatus changeColIntegrality(const HighsInt col,
                                   const HighsVarType integrality);

  /**
   * @brief Change the integrality of multiple columns given by an
   * interval [from_col, to_col]
   */
  HighsStatus changeColsIntegrality(const HighsInt from_col,
                                    const HighsInt to_col,
                                    const HighsVarType* integrality);

  /**
   * @brief Change the integrality of multiple columns given by a set of indices
   */
  HighsStatus changeColsIntegrality(const HighsInt num_set_entries,
                                    const HighsInt* set,
                                    const HighsVarType* integrality);

  /**
   * @brief Change the integrality of multiple columns given by a mask
   * (full length array with 1 => change; 0 => not)
   */
  HighsStatus changeColsIntegrality(const HighsInt* mask,
                                    const HighsVarType* integrality);

  /**
   * @brief Clear the integrality of all columns
   */
  HighsStatus clearIntegrality() {
    this->model_.lp_.integrality_.clear();
    return HighsStatus::kOk;
  }

  /**
   * @brief Change the cost of a column
   */
  HighsStatus changeColCost(const HighsInt col, const double cost);

  /**
   * @brief Change the cost of multiple columns given by an interval [from_col,
   * to_col]
   */
  HighsStatus changeColsCost(const HighsInt from_col, const HighsInt to_col,
                             const double* cost);

  /**
   * @brief Change the cost of multiple columns given by a set of indices
   */
  HighsStatus changeColsCost(const HighsInt num_set_entries,
                             const HighsInt* set, const double* cost);

  /**
   * @brief Change the cost of multiple columns given by a mask
   * (full length array with 1 => change; 0 => not)
   */
  HighsStatus changeColsCost(const HighsInt* mask, const double* cost);

  /**
   * @brief Change the bounds of a column
   */
  HighsStatus changeColBounds(const HighsInt col, const double lower,
                              const double upper);

  /**
   * @brief Change the bounds of multiple columns given by an interval
   * [from_col, to_col]
   */
  HighsStatus changeColsBounds(const HighsInt from_col, const HighsInt to_col,
                               const double* lower, const double* upper);

  /**
   * @brief Change the bounds of multiple columns given by a set of indices
   */
  HighsStatus changeColsBounds(const HighsInt num_set_entries,
                               const HighsInt* set, const double* lower,
                               const double* upper);

  /**
   * @brief Change the bounds of multiple columns given by a mask (full
   * length array with 1 => change; 0 => not)
   */
  HighsStatus changeColsBounds(const HighsInt* mask, const double* lower,
                               const double* upper);

  /**
   * @brief Change the bounds of a row
   */
  HighsStatus changeRowBounds(const HighsInt row, const double lower,
                              const double upper);

  /**
   * @brief Change the bounds of multiple rows given by an interval [from_row,
   * to_row]
   */
  HighsStatus changeRowsBounds(const HighsInt from_row, const HighsInt to_row,
                               const double* lower, const double* upper);

  /**
   * @brief Change the bounds of multiple rows given by a set of indices
   */
  HighsStatus changeRowsBounds(const HighsInt num_set_entries,
                               const HighsInt* set, const double* lower,
                               const double* upper);

  /**
   * @brief Change the bounds of multiple rows given by a mask (full
   * length array with 1 => change; 0 => not)
   */
  HighsStatus changeRowsBounds(const HighsInt* mask, const double* lower,
                               const double* upper);

  /**
   * @brief Change a matrix coefficient
   */
  HighsStatus changeCoeff(const HighsInt row, const HighsInt col,
                          const double value);
  /**
   * @brief Sets the constraint matrix format of the incumbent model
   */
  HighsStatus setMatrixFormat(const MatrixFormat desired_format) {
    this->model_.lp_.setFormat(desired_format);
    return HighsStatus::kOk;
  }

  /**
   * @brief Adds a variable to the incumbent model, without the matrix
   * coefficients if num_new_nz = 0, in which case indices and values
   * arrays can be nullptr
   */
  HighsStatus addCol(const double cost, const double lower, const double upper,
                     const HighsInt num_new_nz, const HighsInt* indices,
                     const double* values);

  /**
   * @brief Adds multiple columns to the incumbent model, without the matrix
   * coefficients if num_new_nz = 0, in which case column-wise starts,
   * indices and values arrays can be nullptr
   */
  HighsStatus addCols(const HighsInt num_new_col, const double* cost,
                      const double* lower, const double* upper,
                      const HighsInt num_new_nz, const HighsInt* starts,
                      const HighsInt* indices, const double* values);

  /**
   * @brief Adds a variable to the incumbent model, without the cost or matrix
   * coefficients
   */
  HighsStatus addVar(const double lower = 0, const double upper = kHighsInf) {
    return this->addVars(1, &lower, &upper);
  }

  /**
   * @brief Adds multiple variables to the incumbent model, without the costs or
   * matrix coefficients
   */
  HighsStatus addVars(const HighsInt num_new_var, const double* lower,
                      const double* upper);

  /**
   * @brief Add a row to the incumbent model, without the matrix coefficients if
   * num_new_nz = 0, in which case indices and values arrays can be
   * nullptr
   */
  HighsStatus addRow(const double lower, const double upper,
                     const HighsInt num_new_nz, const HighsInt* indices,
                     const double* values);

  /**
   * @brief Adds multiple rows to the incumbent model, without the matrix
   * coefficients if num_new_nz = 0, in which case row-wise starts,
   * indices and values arrays can be nullptr
   */
  HighsStatus addRows(const HighsInt num_new_row, const double* lower,
                      const double* upper, const HighsInt num_new_nz,
                      const HighsInt* starts, const HighsInt* indices,
                      const double* values);

  HighsStatus ensureColwise() {
    this->model_.lp_.ensureColwise();
    return HighsStatus::kOk;
  }

  HighsStatus ensureRowwise() {
    this->model_.lp_.ensureRowwise();
    return HighsStatus::kOk;
  }

  /**
   * @brief Delete multiple columns from the incumbent model given by an
   * interval [from_col, to_col]
   */
  HighsStatus deleteCols(const HighsInt from_col, const HighsInt to_col);

  /**
   * @brief Delete multiple columns from the incumbent model given by a set
   */
  HighsStatus deleteCols(const HighsInt num_set_entries, const HighsInt* set);

  /**
   * @brief Delete multiple columns from the incumbent model given by
   * a mask (full length array with 1 => delete; 0 => keep). New index
   * of any column kept is returned in place of the value 0.  For
   * deleted columns, a value of -1 is returned.
   */
  HighsStatus deleteCols(HighsInt* mask);

  /**
   * @brief Delete multiple variables from the incumbent model given by an
   * interval [from_var, to_var]
   */
  HighsStatus deleteVars(const HighsInt from_var, const HighsInt to_var) {
    return deleteCols(from_var, to_var);
  }

  /**
   * @brief Delete multiple variables from the incumbent model given by a set
   */
  HighsStatus deleteVars(const HighsInt num_set_entries, const HighsInt* set) {
    return deleteCols(num_set_entries, set);
  }

  /**
   * @brief Delete multiple variables from the incumbent model given
   * by a mask (full length array with 1 => delete; 0 => keep). New
   * index of any variable not deleted is returned in place of the
   * value 0. For deleted variables, a value of -1 is returned.
   */
  HighsStatus deleteVars(HighsInt* mask) { return deleteCols(mask); }

  /**
   * @brief Delete multiple rows from the incumbent model given by an interval
   * [from_row, to_row]
   */
  HighsStatus deleteRows(const HighsInt from_row, const HighsInt to_row);

  /**
   * @brief Delete multiple rows from the incumbent model given by a set
   */
  HighsStatus deleteRows(const HighsInt num_set_entries, const HighsInt* set);

  /**
   * @brief Delete multiple rows from the incumbent model given by a
   * mask (full length array with 1 => delete; 0 => keep). New index
   * of any row not deleted is returned in place of the value 0. For
   * deleted rows, a value of -1 is returned.
   */
  HighsStatus deleteRows(HighsInt* mask);

  /**
   * @brief Scale a matrix column (and cost) by a constant - flipping bounds if
   * the constant is negative
   */
  HighsStatus scaleCol(const HighsInt col, const double scale_value);

  /**
   * @brief Scale a matrix row by a constant - flipping bounds if the constant
   * is negative
   */
  HighsStatus scaleRow(const HighsInt row, const double scale_value);

  /**
   * Other methods for specialist applications
   */

  /**
   * Methods for setting basis_ and solution_
   */

  /**
   * @brief Pass a HighsSolution instance to set the internal
   * HighsSolution instance. If any of col_value, col_dual and
   * row_dual is not set, the internal HighsSolution instance is not
   * updated
   */
  HighsStatus setSolution(const HighsSolution& solution);

  /**
   * @brief Pass a primal solution. If index is not a null pointer,
   * then it is assumed that value contains num_entries of packed
   * values, with index defing the corresponding primal solution
   * components. If index is a null pointer, then value is assumed to
   * be a full primal solution.
   *
   * It allows a full primal solution to be passed (for MIPs) without
   * requiring either a HighsSolution that contains (empty vectors of)
   * spurious dual information, or a full list of indices
   */
  HighsStatus setSolution(const HighsInt num_entries, const HighsInt* index,
                          const double* value);

  /**
   * @brief Set the callback method to use for HiGHS
   */
  HighsStatus setCallback(HighsCallbackFunctionType user_callback,
                          void* user_callback_data = nullptr);
  HighsStatus setCallback(HighsCCallbackType c_callback,
                          void* user_callback_data = nullptr);

  /**
   * @brief Start callback of given type
   */
  HighsStatus startCallback(const int callback_type);
  HighsStatus startCallback(const HighsCallbackType callback_type);

  /**
   * @brief Stop callback of given type
   */
  HighsStatus stopCallback(const int callback_type);
  HighsStatus stopCallback(const HighsCallbackType callback_type);

  /**
   * @brief Use the HighsBasis passed to set the internal HighsBasis
   * instance. The origin string is used to identify the source of the
   * HighsBasis instance.
   */
  HighsStatus setBasis(const HighsBasis& basis, const std::string& origin = "");

  /**
   * @brief Clear the internal HighsBasis instance
   */
  HighsStatus setBasis();

  /**
   * @brief Report profiling
   */
  void reportProfiling() const;

  /**
   * @brief Run IPX crossover from a given HighsSolution instance and,
   * if successful, set the internal HighsBasis and HighsSolution
   * instance
   */
  HighsStatus crossover(const HighsSolution& user_solution);

  /**
   * @brief Open a named log file
   */
  HighsStatus openLogFile(const std::string& log_file = "");

  /**
   * @brief Close any open log file
   */
  HighsStatus closeLogFile();

  /**
   * @brief Interpret common qualifiers to string values
   */
  std::string highsStatusToString(const HighsStatus status) const;
  std::string presolveStatusToString(
      const HighsPresolveStatus presolve_status) const;
  std::string modelStatusToString(const HighsModelStatus model_status) const;
  std::string solutionStatusToString(const HighsInt solution_status) const;
  std::string basisStatusToString(const HighsBasisStatus basis_status) const;
  std::string basisValidityToString(const HighsInt basis_validity) const;
  std::string presolveRuleTypeToString(const HighsInt presolve_rule) const;

  /**
   * @brief Ensures that the scheduler is initialized,
   * returning HighsStatus::kError if it has already been initialized,
   * but the threads option is nonzero and not equal to
   * this->max_threads_.
   */
  HighsStatus initializeMultiThreading();

  /**
   * @brief Releases all resources held by the scheduler owned by the calling
   * thread. It is thread-safe to call this function only from the thread that
   * owns the scheduler (not from a worker thread) and when no run() or
   * presolve() is in progress on that thread. After this function has
   * terminated it is guaranteed that eventually all previously created
   * scheduler threads will terminate and allocated memory will be released.
   * After this function has returned the option value for the number of threads
   * may be altered to a new value before the next call to run() or presolve().
   * If the given bool parameter has value true, then the function will not
   * return until all memory is freed, which might be desirable when debugging
   * heap memory but requires the calling thread to wait for all scheduler
   * threads to wake-up which is usually not necessary.
   */
  static void resetThreadScheduler(bool blocking = false);

  /**
   * @brief Initialize optional extras and cache the result.
   *
   * An empty directory searches default library directories only when
   * search_if_empty is true. Statically linked extras ignore both arguments.
   *
   * @return HighsStatus::kOk on success, HighsStatus::kError on failure.
   * Subsequent calls return the cached result without retrying, unless
   * extras are explicitly unloaded.
   */
  static HighsStatus initializeExtras(const std::string& directory,
                                      bool search_if_empty = true);

  /**
   * @brief If profiling is not nullptr, sets up profiling and copies
   * its pointer to Highs
   */
  void initializeProfiling(HighsProfiling* profiling);
  void initializeSingleThreadedProfiling(HighsProfiling* profiling);

  /**
   * @brief Clears and then initializes profiling
   */
  void resetProfiling();

  /**
   * @brief If Highs::profiling_ is not nullptr, clears profiling and
   * sets Highs::profiling_ to nullptr
   */
  void clearProfiling();

  /**
   * @brief Checks that pointer is not nullptr, and copies it to Highs
   */
  void setProfiling(HighsProfiling* profiling);

  // Start of advanced methods: only for internal use!

  // Nested methods below Highs::run()
  //
  // See highs/HighsRun.md
  HighsStatus optimizeHighs();
  // Used in MIP solver as minimal LP solve
  HighsStatus optimizeLp();
  HighsStatus optimizeModel();
  HighsStatus optimizeModelTryCatch();
  HighsStatus calledOptimizeModel();

  const HighsSimplexStats& getSimplexStats() const {
    return ekk_instance_.getSimplexStats();
  }
  void reportSimplexStats(FILE* file) const {
    ekk_instance_.reportSimplexStats(file);
  }

  /**
   * @brief Put a copy of the current iterate - basis; invertible
   * representation and dual edge weights - into storage within
   * HSimplexNla. Advanced method: for HiGHS MIP solver
   */
  HighsStatus putIterate();

  /**
   * @brief Get a copy of the iterate stored within HSimplexNla and
   * overwrite the current iterate. Advanced method: for HiGHS MIP
   * solver
   */
  HighsStatus getIterate();

  /**
   * @brief Get the dual edge weights (steepest/devex) in the order of
   * the basic indices or nullptr when they are not available.
   */
  const double* getDualEdgeWeights() const {
    return ekk_instance_.status_.has_dual_steepest_edge_weights
               ? ekk_instance_.dual_edge_weight_.data()
               : nullptr;
  }

  /**
   * @brief Gets the internal basic variable index array in the order
   * corresponding to calls to getBasisInverseRow, getBasisInverseCol,
   * getBasisSolve, getBasisTransposeSolve, getReducedRow and getReducedColumn.
   * Entries are indices of columns if in [0,num_col), and entries in [num_col,
   * num_col+num_row) are (num_col+row_index).
   */
  const HighsInt* getBasicVariablesArray() const {
    assert(ekk_instance_.status_.has_invert);
    return ekk_instance_.basis_.basicIndex_.data();
  }

  /**
   * @brief Form a row of \f$B^{-1}\f$ for basis matrix \f$B\f$,
   * returning the result in the given HVector buffer which is
   * expected to be setup with dimension num_row. The buffers
   * previous contents will be overwritten.
   */
  HighsStatus getBasisInverseRowSparse(const HighsInt row,
                                       HVector& row_ep_buffer);

  /**
   * @brief Get the primal simplex phase 1 dual values. Advanced
   * method: for HiGHS IIS calculation
   */
  const std::vector<double>& getPrimalPhase1Dual() const {
    return ekk_instance_.primal_phase1_dual_;
  }

  /**
   * @brief Generalisation of getColName and getRowName. Advanced
   * method: for HiGHS C++ and C API
   */
  HighsStatus getColOrRowName(const HighsLp& lp, const bool is_col,
                              const HighsInt index, std::string& name) const;

  /**
   * @brief Development methods
   */
  HighsInt defineClock(const char* name) {
    return this->timer_.clock_def(name);
  }
  void writeAllClocks() { this->timer_.writeAllClocks(); }
  HighsStatus clearModelNames() {
    this->model_.lp_.col_names_.clear();
    this->model_.lp_.row_names_.clear();
    return HighsStatus::kOk;
  }

  // Start of deprecated methods

  std::string compilationDate() const { return "deprecated"; }

  HighsStatus setLogCallback(void (*user_log_callback)(HighsLogType,
                                                       const char*, void*),
                             void* user_log_callback_data = nullptr);

  HighsInt getNumCols() const {
    deprecationMessage("getNumCols", "getNumCol");
    return getNumCol();
  }
  HighsInt getNumRows() const {
    deprecationMessage("getNumRows", "getNumRow");
    return getNumRow();
  }
  HighsInt getNumEntries() {
    deprecationMessage("getNumEntries", "getNumNz");
    return getNumNz();
  }

  HighsStatus setHighsOptionValue(const std::string& option, const bool value);

  HighsStatus setHighsOptionValue(const std::string& option,
                                  const HighsInt value);

#ifdef HIGHSINT64
  HighsStatus setHighsOptionValue(const std::string& option,
                                  const int value  //!< The option value
  ) {
    deprecationMessage("setHighsOptionValue", "setOptionValue");
    return setOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setHighsOptionValue(const std::string& option,
                                  const double value);

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const std::string& value  //!< The option value
  );

  HighsStatus setHighsOptionValue(const std::string& option, const char* value);

  HighsStatus readHighsOptions(const std::string& filename  //!< The filename
  );

  HighsStatus passHighsOptions(const HighsOptions& options  //!< The options
  );

  HighsStatus getHighsOptionValue(const std::string& option, bool& value);

  HighsStatus getHighsOptionValue(const std::string& option, HighsInt& value);

  HighsStatus getHighsOptionValue(const std::string& option, double& value);

  HighsStatus getHighsOptionValue(const std::string& option,
                                  std::string& value);

  HighsStatus getHighsOptionType(const std::string& option,
                                 HighsOptionType& type  //!< The option type
  );

  const HighsOptions& getHighsOptions() const;

  HighsStatus resetHighsOptions();

  HighsStatus writeHighsOptions(const std::string& filename,  //!< The filename
                                const bool report_only_deviations = true);

  HighsInt getSimplexIterationCount() {
    deprecationMessage("getSimplexIterationCount", "None");
    return info_.simplex_iteration_count;
  }

  HighsStatus setHighsLogfile(FILE* logfile = nullptr);

  HighsStatus setHighsOutput(FILE* output = nullptr);

  const HighsInfo& getHighsInfo() const;

  HighsStatus getHighsInfoValue(const std::string& info, HighsInt& value);

  HighsStatus getHighsInfoValue(const std::string& info,
                                double& value) const;  //!< The info value

  HighsStatus writeHighsInfo(const std::string& filename  //!< The filename
  );

  double getHighsInfinity();

  double getHighsRunTime();

  const HighsModelStatus& getModelStatus(const bool scaled_model) const;

  void logHeader() const;

  void deprecationMessage(const std::string& method_name,
                          const std::string& alt_method_name) const;

  /**
   * @brief Get the hot start basis data from the most recent simplex
   * solve. Advanced method: for HiGHS MIP solver
   */
  const HotStart& getHotStart() const { return ekk_instance_.hot_start_; }

  /**
   * @brief Set up for simplex using the supplied hot start
   * data. Advanced method: for HiGHS MIP solver
   */
  HighsStatus setHotStart(const HotStart& hot_start) {
    this->deprecationMessage("setHotStart", "None");
    return HighsStatus::kError;
  }

  /**
   * @brief Freeze the current internal HighsBasis instance and
   * standard NLA, returning a value to be used to recover this basis
   * and standard NLA at minimal cost. Advanced method: for HiGHS MIP
   * solver
   */
  HighsStatus freezeBasis(HighsInt& frozen_basis_id) {
    this->deprecationMessage("freezeBasis", "None");
    return HighsStatus::kError;
  }

  /**
   * @brief Unfreeze a frozen HighsBasis instance and standard NLA (if
   * possible). Advanced method: for HiGHS MIP solver
   */
  HighsStatus unfreezeBasis(const HighsInt frozen_basis_id) {
    this->deprecationMessage("unfreezeBasis", "None");
    return HighsStatus::kError;
  }

  /**
   * @brief Check that all frozen basis data has been
   * cleared. Advanced method: for HiGHS MIP solver
   */
  HighsStatus frozenBasisAllDataClear() {
    this->deprecationMessage("frozenBasisAllDataClear", "None");
    return HighsStatus::kError;
  }

  // Deprecated in favour of resetThreadScheduler().
  // Cannot call deprecationMessage here.
  static void resetGlobalScheduler(bool blocking = false);

  // End of deprecated methods
 private:
  HighsSolution solution_;
  HighsBasis basis_;
  ICrashInfo icrash_info_;

  HighsModel model_;
  std::vector<HighsLinearObjective> multi_linear_objective_;

  HighsModel presolved_model_;
  HighsTimer timer_;

  HighsCallback callback_;
  HighsOptions options_;
  HighsInfo info_;
  HighsRunData run_data_;
  HighsRanging ranging_;
  HighsIis iis_;
  std::vector<HighsObjectiveSolution> saved_objective_and_solution_;

  HighsPresolveStatus model_presolve_status_ =
      HighsPresolveStatus::kNotPresolved;
  HighsModelStatus model_status_ = HighsModelStatus::kNotset;

  bool standard_form_valid_ = false;
  double standard_form_offset_;
  std::vector<double> standard_form_cost_;
  std::vector<double> standard_form_rhs_;
  HighsSparseMatrix standard_form_matrix_;

  HEkk ekk_instance_;

  HighsPresolveLog presolve_log_;

  HighsProfiling* profiling_ = nullptr;

  HighsInt max_threads_ = 0;
  // This is strictly for debugging. It's used to check whether
  // returnFromOptimizeModel() was called after the previous call to
  // Highs::optimizeModel() and, assuming that this is always done, it checks
  // whether Highs::optimizeModel() is called recursively.
  bool called_return_from_optimize_model = true;
  HighsInt debug_optimize_call_num_ = 0;

  void reportModelStats() const;

  void exactResizeModel() {
    this->model_.lp_.exactResize();
    this->model_.hessian_.exactResize();
  }

  HighsStatus completeSolutionFromDiscreteAssignment();

  HighsStatus callSolveLp(HighsLp& lp, const std::string& message);
  HighsStatus callSolveMip(HighsLp& lp, const std::string& message);
  HighsStatus callSolveQp(HighsModel& model, const std::string& message);

  HighsStatus callRunPostsolve(const HighsSolution& solution,
                               const HighsBasis& basis);

  PresolveComponent presolve_;
  HighsPresolveStatus runPresolve(const bool force_lp_presolve,
                                  const bool force_presolve = false);
  HighsPostsolveStatus runPostsolve();

  HighsStatus openWriteFile(const std::string& filename,
                            const std::string& method_name, FILE*& file,
                            HighsFileType& file_type) const;

  void reportModel(const HighsModel& model);
  void newHighsBasis();
  void forceHighsSolutionBasisSize();
  //
  // For cases where there is no solution data for the model, but its
  // status is proved otherwise. Sets the model status, then clears any solution
  // and basis data
  void setHighsModelStatusAndClearSolutionAndBasis(
      const HighsModelStatus model_status);

  // Clears derived model properties (like any presolved model,
  // standard form LP, and ray information) that (unlike solution and
  // basis) cannot be updated
  void clearDerivedModelProperties();

  // Clears the presolved model and its status
  void clearPresolve();

  // Clears the standard form LP
  void clearStandardFormLp();

  // Clears the ray records
  void clearRayRecords() { this->ekk_instance_.clearRayRecords(); }
  //
  // Methods to clear solver data for users in Highs class members
  // before (possibly) updating them with data from trying to solve
  // the incumbent model.
  //
  // Invalidates all solver data in Highs class members by calling
  // invalidateModelStatus(), invalidateSolution(), invalidateBasis(),
  // invalidateRanging(), invalidateInfo(), invalidateRunData(),
  // invalidateEkk() and clearIis()
  void invalidateSolverData();

  // Invalidates all solver dual data in Highs class members by calling
  // invalidateModelStatus(), invalidateRanging(), and invalidateInfo()
  //
  // Used when only the objective changes
  void invalidateSolverDualData();
  //
  // Invalidates the model status, solution_ and info_
  void invalidateModelStatusSolutionAndInfo();
  //
  // Invalidates the model status and info_
  void invalidateModelStatusAndInfo();
  //
  // Sets model status to HighsModelStatus::kNotset
  void invalidateModelStatus();
  //
  // Invalidates primal and dual solution
  void invalidateSolution();
  //
  // Invalidates basis
  void invalidateBasis();
  //
  // Invalidates info_ and resets the values of its members
  void invalidateInfo();
  //
  // Invalidates run_data_ and resets the values of its members
  void invalidateRunData();
  //
  // Invalidates ranging_
  void invalidateRanging();

  // Invalidates ekk_instance_
  void invalidateEkk();

  // Clears iis_
  void clearIis();

  HighsStatus returnFromWriteSolution(FILE* file,
                                      const HighsStatus return_status);
  HighsStatus returnFromOptimizeModel(const HighsStatus return_status,
                                      const bool undo_mods);
  HighsStatus returnFromHighs(const HighsStatus return_status);
  void reportSolvedLpQpStats();

  // Interface methods
  HighsStatus formStandardFormLp();
  HighsStatus basisForSolution();
  HighsStatus addColsInterface(
      HighsInt ext_num_new_col, const double* ext_col_cost,
      const double* ext_col_lower, const double* ext_col_upper,
      HighsInt ext_num_new_nz, const HighsInt* ext_a_start,
      const HighsInt* ext_a_index, const double* ext_a_value);

  HighsStatus addRowsInterface(HighsInt ext_num_new_row,
                               const double* ext_row_lower,
                               const double* ext_row_upper,
                               HighsInt ext_num_new_nz,
                               const HighsInt* ext_ar_start,
                               const HighsInt* ext_ar_index,
                               const double* ext_ar_value);

  void deleteColsInterface(HighsIndexCollection& index_collection);

  void deleteRowsInterface(HighsIndexCollection& index_collection);

  void getColsInterface(const HighsIndexCollection& index_collection,
                        HighsInt& num_col, double* cost, double* lower,
                        double* upper, HighsInt& num_nz, HighsInt* start,
                        HighsInt* index, double* value) const;

  void getRowsInterface(const HighsIndexCollection& index_collection,
                        HighsInt& num_row, double* lower, double* upper,
                        HighsInt& num_nz, HighsInt* start, HighsInt* index,
                        double* value) const;

  void getCoeffInterface(const HighsInt ext_row, const HighsInt ext_col,
                         double& value) const;

  HighsStatus changeObjectiveSenseInterface(const ObjSense ext_sense);
  HighsStatus changeObjectiveOffsetInterface(const double ext_offset);
  HighsStatus changeIntegralityInterface(HighsIndexCollection& index_collection,
                                         const HighsVarType* usr_integrality);
  HighsStatus changeCostsInterface(HighsIndexCollection& index_collection,
                                   const double* usr_col_cost);

  bool feasibleWrtBounds(const bool columns = true) const;
  HighsStatus changeColBoundsInterface(HighsIndexCollection& index_collection,
                                       const double* usr_col_lower,
                                       const double* usr_col_upper);
  HighsStatus changeRowBoundsInterface(HighsIndexCollection& index_collection,
                                       const double* usr_row_lower,
                                       const double* usr_row_upper);
  HighsStatus changeCoeffInterface(const HighsInt ext_row,
                                   const HighsInt ext_col,
                                   const double ext_new_value);
  HighsStatus scaleColInterface(const HighsInt col, const double scale_value);
  HighsStatus scaleRowInterface(const HighsInt row, const double scale_value);

  void setNonbasicStatusInterface(const HighsIndexCollection& index_collection,
                                  const bool columns);
  void appendNonbasicColsToBasisInterface(const HighsInt ext_num_new_col);
  void appendBasicRowsToBasisInterface(const HighsInt ext_num_new_row);

  HighsStatus getBasicVariablesInterface(HighsInt* basic_variables);
  HighsStatus basisSolveInterface(const vector<double>& rhs,
                                  double* solution_vector,
                                  HighsInt* solution_num_nz,
                                  HighsInt* solution_indices, bool transpose);

  void zeroIterationCounts();

  HighsStatus getDualRayInterface(bool& has_dual_ray, double* dual_ray_value);

  HighsStatus getPrimalRayInterface(bool& has_primal_ray,
                                    double* primal_ray_value);
  HighsStatus getRangingInterface();

  HighsStatus getIisInterface();
  HighsStatus getIisInterfaceReturn(
      const HighsStatus return_status, const HighsOptions& original_options,
      const std::vector<HighsBool>& original_callbacks);

  HighsStatus elasticityFilterReturn(
      const HighsStatus return_status, const std::string& original_model_name,
      const HighsModelStatus original_model_status,
      const HighsInt original_num_col, const HighsInt original_num_row,
      const ObjSense original_sense,
      const std::vector<double>& original_col_cost,
      const std::vector<double>& original_col_lower,
      const std::vector<double>& original_col_upper,
      const std::vector<HighsVarType>& original_integrality);
  HighsStatus elasticityFilter(const double global_lower_penalty,
                               const double global_upper_penalty,
                               const double global_rhs_penalty,
                               const double* local_lower_penalty,
                               const double* local_upper_penalty,
                               const double* local_rhs_penalty,
                               const bool get_iis = false);
  HighsStatus extractIis(HighsInt& num_iis_col, HighsInt& num_iis_row,
                         HighsInt* iis_col_index, HighsInt* iis_row_index,
                         HighsInt* iis_col_bound, HighsInt* iis_row_bound);

  HighsStatus returnFromLexicographicOptimization(
      const HighsStatus return_status, HighsInt original_lp_num_row);
  HighsStatus multiobjectiveSolve();

  bool aFormatOk(const HighsInt num_nz, const HighsInt format);
  bool qFormatOk(const HighsInt num_nz, const HighsInt format);
  void clearZeroHessian();
  void callLpKktCheck(const HighsLp& lp, const std::string& message = "");
  HighsStatus invertRequirementError(const std::string& method_name) const;

  HighsStatus handleInfCost();
  void restoreInfCost(HighsStatus& return_status);
  HighsStatus optionChangeAction();

  HighsStatus userScale(HighsUserScaleData& data);
  HighsStatus userUnscale(HighsUserScaleData& data);
  HighsStatus userScaleModel(HighsUserScaleData& data);
  HighsStatus userScaleSolution(HighsUserScaleData& data,
                                bool update_kkt = false);

  HighsStatus computeIllConditioning(HighsIllConditioning& ill_conditioning,
                                     const bool constraint,
                                     const HighsInt method,
                                     const double ill_conditioning_bound);
  void formIllConditioningLp0(HighsLp& ill_conditioning_lp,
                              std::vector<HighsInt>& basic_var,
                              const bool constraint);
  void formIllConditioningLp1(HighsLp& ill_conditioning_lp,
                              std::vector<HighsInt>& basic_var,
                              const bool constraint,
                              const double ill_conditioning_bound);
  bool processedBoundsFindsInfeasibility();
  bool validLinearObjective(const HighsLinearObjective& linear_objective,
                            const HighsInt iObj) const;
  bool hasRepeatedLinearObjectivePriorities(
      const HighsLinearObjective* linear_objective = nullptr) const;

  bool tryPdlpCleanup(HighsInt& pdlp_cleanup_iteration_limit,
                      const HighsInfo& presolved_lp_info) const;
};

// Start of deprecated methods not in the Highs class

const char* highsCompilationDate();

#endif
