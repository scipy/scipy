/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsCallback.h"
#include "lp_data/HighsOptions.h"
#include "parallel/HighsParallel.h"

struct HighsMipSolverData;
class HighsCutPool;
struct HighsPseudocostInitialization;
class HighsCliqueTable;
class HighsImplications;

struct HighsTerminator {
  HighsInt num_instance;
  HighsInt my_instance;
  HighsModelStatus* record;
  void clear();
  void initialise(HighsInt num_instance_, HighsInt my_instance_,
                  HighsModelStatus* record_);
  HighsInt concurrency() const;
  void terminate();
  bool terminated() const;
  HighsModelStatus terminationStatus() const;
  void report(const HighsLogOptions& log_options) const;
};

struct MipViolation {
  double bound_violation;
  double integrality_violation;
  double row_violation;
  HighsInt num_bound_violations;
  HighsInt num_integrality_violations;
  HighsInt num_row_violations;
  HighsInt col_of_max_bound_violation;
  HighsInt col_of_max_integrality_violation;
  HighsInt row_of_max_row_violation;
  void clear();
  void copy(double& bound_violation_, double& row_violation_,
            double& integrality_violation_) const;
  void log(const HighsLogOptions& log_options, const double objective_value,
           const std::string& source) const;
};

class HighsMipSolver {
 public:
  HighsCallback* callback_;
  const HighsOptions* options_mip_;
  const HighsLp* model_;
  const HighsLp* orig_model_;
  HighsModelStatus modelstatus_;
  std::vector<double> solution_;
  double solution_objective_;
  double bound_violation_;
  double integrality_violation_;
  double row_violation_;
  // The following are only to return data to HiGHS, and are set in
  // HighsMipSolver::cleanupSolve
  double dual_bound_;
  double primal_bound_;
  double gap_;
  int64_t node_count_;
  int64_t total_lp_iterations_;
  double primal_dual_integral_;

  FILE* improving_solution_file_;
  std::vector<HighsObjectiveSolution> saved_objective_and_solution_;

  bool submip;
  HighsInt submip_level;
  HighsInt max_submip_level;
  const HighsBasis* rootbasis;
  const HighsPseudocostInitialization* pscostinit;
  const HighsCliqueTable* clqtableinit;
  const HighsImplications* implicinit;

  std::unique_ptr<HighsMipSolverData> mipdata_;

  HighsModelStatus termination_status_;
  HighsTerminator terminator_;

  void run();

  HighsInt numCol() const { return model_->num_col_; }

  HighsInt numRow() const { return model_->num_row_; }

  HighsInt numNonzero() const { return model_->a_matrix_.numNz(); }

  const double* colCost() const { return model_->col_cost_.data(); }

  double colCost(HighsInt col) const { return model_->col_cost_[col]; }

  const double* rowLower() const { return model_->row_lower_.data(); }

  double rowLower(HighsInt col) const { return model_->row_lower_[col]; }

  const double* rowUpper() const { return model_->row_upper_.data(); }

  double rowUpper(HighsInt col) const { return model_->row_upper_[col]; }

  const HighsVarType* variableType() const {
    return model_->integrality_.data();
  }

  HighsVarType variableType(HighsInt col) const {
    return model_->integrality_[col];
  }

  bool isColIntegral(const HighsInt col) const {
    return variableType(col) != HighsVarType::kContinuous;
  }

  bool isColInteger(const HighsInt col) const {
    return variableType(col) == HighsVarType::kInteger;
  }

  bool isColContinuous(const HighsInt col) const {
    return variableType(col) == HighsVarType::kContinuous;
  }

  HighsMipSolver(HighsCallback& callback, const HighsOptions& options,
                 const HighsLp& lp, const HighsSolution& solution,
                 bool submip = false, HighsInt submip_level = 0);

  ~HighsMipSolver();

  template <class F>
  void runTask(F&& f, highs::parallel::TaskGroup& tg, bool parallel_lock,
               bool force_serial,
               const std::vector<HighsInt>& indices = std::vector<HighsInt>(1,
                                                                            0));

  void setModel(const HighsLp& model) {
    model_ = &model;
    solution_objective_ = kHighsInf;
  }

  mutable HighsTimer timer_;
  HighsProfiling* profiling_ = nullptr;

  void cleanupSolve();
  void solvingReport(const std::string& solutionstatus) const;

  void runMipPresolve(const HighsInt presolve_reduction_limit);
  const HighsLp& getPresolvedModel() const;
  HighsPresolveStatus getPresolveStatus() const;
  presolve::HighsPostsolveStack getPostsolveStack() const;

  void callbackGetCutPool() const;
  bool solutionFeasible(const HighsLp* lp, const std::vector<double>& col_value,
                        const std::vector<double>* pass_row_value,
                        MipViolation& violation, HighsCDouble& obj) const;

  std::vector<HighsModelStatus> initialiseTerminatorRecord(
      HighsInt num_instance) const;
  void initialiseTerminator(HighsInt num_instance_ = 0,
                            HighsInt my_instance_ = kNoThreadInstance,
                            HighsModelStatus* record_ = nullptr);
  void initialiseTerminator(const HighsMipSolver& mip_solver);
  bool terminate() const {
    return this->termination_status_ != HighsModelStatus::kNotset;
  }
  HighsModelStatus terminationStatus() const {
    return this->termination_status_;
  }
  void setParallelLock(bool lock) const;
  void setProfiling(HighsProfiling* profiling);
  HighsInt getMaxNumWorkers() const {
    if (options_mip_->parallel == kHighsOffString || submip) {
      return 1;
    }
    return highs::parallel::num_threads();
  }
};

std::array<char, 128> getGapString(const double gap_,
                                   const double primal_bound_,
                                   const HighsOptions* options_mip_);
#endif
