/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsCallback.cpp
 * @brief
 */
#include "lp_data/HighsCallback.h"

#include <algorithm>
#include <cassert>
#include <utility>

#include "Highs.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsIntegers.h"

void HighsCallback::clearHighsCallbackOutput() {
  this->data_out.log_type = HighsLogType::kInfo;
  this->data_out.running_time = -1;
  this->data_out.simplex_iteration_count = -1;
  this->data_out.ipm_iteration_count = -1;
  this->data_out.pdlp_iteration_count = -1;
  this->data_out.qpasm_iteration_count = -1;
  this->data_out.objective_function_value = -kHighsInf;
  this->data_out.mip_node_count = -1;
  this->data_out.mip_primal_bound = kHighsInf;
  this->data_out.mip_dual_bound = -kHighsInf;
  this->data_out.mip_gap = -1;
  this->data_out.mip_solution.clear();
  this->data_out.cutpool_start.clear();
  this->data_out.cutpool_index.clear();
  this->data_out.cutpool_value.clear();
  this->data_out.cutpool_lower.clear();
  this->data_out.cutpool_upper.clear();
  this->data_out.qp_solution.clear();
  this->data_out.external_solution_query_origin =
      ExternalMipSolutionQueryOrigin::kExternalMipSolutionQueryOriginAfterSetup;
}

void HighsCallback::clearHighsCallbackInput() {
  size_t num_col = highs != nullptr ? highs->getNumCol() : 0;

  // make sure buffer size is correct and reset the contents if previously used
  if (this->data_in.user_has_solution ||
      num_col != this->data_in.user_solution.size()) {
    this->data_in.user_solution.assign(num_col, kHighsUndefined);
  }

  this->data_in.user_interrupt = false;
  this->data_in.user_has_solution = false;
}

void HighsCallback::clear() {
  this->user_callback = nullptr;
  this->user_callback_data = nullptr;
  this->active.assign(kNumCallbackType, false);
  this->clearHighsCallbackOutput();
  this->clearHighsCallbackInput();
}

bool HighsCallback::callbackActive(const int callback_type) {
  // Check that callback function has been defined
  if (!this->user_callback) return false;
  // Check that callback_type is within range
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return false;
  // Don't call callback if it is not active
  assert(this->active.size() > 0);
  if (!this->active[callback_type]) return false;
  return true;
}

bool HighsCallback::callbackAction(const int callback_type,
                                   const std::string& message) {
  if (!callbackActive(callback_type)) return false;
  this->user_callback(callback_type, message.c_str(), &this->data_out,
                      &this->data_in, this->user_callback_data);
  // Assess any action
  bool action = this->data_in.user_interrupt;

  // Check for no action if case not handled internally
  if (callback_type == kCallbackMipImprovingSolution ||
      callback_type == kCallbackMipSolution ||
      callback_type == kCallbackMipLogging ||
      callback_type == kCallbackMipGetCutPool ||
      callback_type == kCallbackMipDefineLazyConstraints ||
      callback_type == kCallbackMipUserSolution ||
      callback_type == kCallbackQpFirstFeasiblePoint)
    assert(!action);
  return action;
}

// Conversions for C API

// Return HighsCallbackDataOut (ie pointer-based) data structure from
// HighsCallbackOutput
HighsCallbackOutput::operator HighsCallbackDataOut() const {
  HighsCallbackDataOut data;
  data.cbdata = static_cast<void*>(const_cast<HighsCallbackOutput*>(this));
  data.log_type = static_cast<int>(log_type);
  data.running_time = running_time;
  data.simplex_iteration_count = simplex_iteration_count;
  data.ipm_iteration_count = ipm_iteration_count;
  data.pdlp_iteration_count = pdlp_iteration_count;
  data.qpasm_iteration_count = qpasm_iteration_count;
  data.objective_function_value = objective_function_value;

  data.mip_node_count = mip_node_count;
  data.mip_total_lp_iterations = mip_total_lp_iterations;
  data.mip_primal_bound = mip_primal_bound;
  data.mip_dual_bound = mip_dual_bound;
  data.mip_gap = mip_gap;
  data.mip_solution_size = mip_solution.size();
  data.mip_solution =
      mip_solution.empty() ? nullptr : const_cast<double*>(mip_solution.data());

  data.cutpool_num_col = cutpool_num_col;
  data.cutpool_num_cut = cutpool_lower.size();
  data.cutpool_num_nz = cutpool_value.size();
  data.cutpool_start = cutpool_start.empty()
                           ? nullptr
                           : const_cast<HighsInt*>(cutpool_start.data());
  data.cutpool_index = cutpool_index.empty()
                           ? nullptr
                           : const_cast<HighsInt*>(cutpool_index.data());
  data.cutpool_value = cutpool_value.empty()
                           ? nullptr
                           : const_cast<double*>(cutpool_value.data());
  data.cutpool_lower = cutpool_lower.empty()
                           ? nullptr
                           : const_cast<double*>(cutpool_lower.data());
  data.cutpool_upper = cutpool_upper.empty()
                           ? nullptr
                           : const_cast<double*>(cutpool_upper.data());

  data.qp_solution_size = qp_solution.size();
  data.qp_solution =
      qp_solution.empty() ? nullptr : const_cast<double*>(qp_solution.data());

  data.external_solution_query_origin =
      static_cast<HighsInt>(external_solution_query_origin);

  return data;
}

HighsCallbackInput::operator HighsCallbackDataIn() const {
  HighsCallbackDataIn data;
  data.cbdata = static_cast<void*>(const_cast<HighsCallbackInput*>(this));
  data.user_interrupt = user_interrupt ? 1 : 0;
  data.user_has_solution = user_has_solution ? 1 : 0;
  data.user_solution_size = user_solution.size();
  data.user_solution = user_solution.empty()
                           ? nullptr
                           : const_cast<double*>(user_solution.data());
  return data;
}

// we assume that user_solution.data() == data_in.user_solution
// and that user_solution.size() == data_in.user_solution_size
HighsCallbackInput HighsCallbackInput::operator=(
    const HighsCallbackDataIn& data_in) {
  assert(user_solution.data() == data_in.user_solution);
  user_interrupt = data_in.user_interrupt != 0;
  user_has_solution = data_in.user_has_solution != 0;
  return *this;
}

// User provides a full solution
HighsStatus HighsCallbackInput::setSolution(HighsInt num_entries,
                                            const double* value) {
  return this->setSolution(num_entries, nullptr, value);
}

// User provides a partial solution which, in a special case could -
// via the call from setSolution(num_entries, value), or by a user -
// be a full solution with/without index being nullptr
HighsStatus HighsCallbackInput::setSolution(HighsInt num_entries,
                                            const HighsInt* index,
                                            const double* value) {
  const auto& lp = highs->getLp();
  const auto& options = highs->getOptions();

  assert(user_solution.size() == static_cast<size_t>(lp.num_col_));

  HighsStatus return_status =
      assessUserColValues(lp, options, num_entries, index, value);
  if (return_status != HighsStatus::kOk) return return_status;

  assignUserColValues(lp.num_col_, user_solution, num_entries, index, value);

  user_has_solution = true;
  return HighsStatus::kOk;
}

HighsStatus HighsCallbackInput::repairSolution() {
  if (!user_has_solution) {
    highsLogUser(highs->getOptions().log_options, HighsLogType::kError,
                 "repairSolution: No user solution has been set\n");
    return HighsStatus::kError;
  } else if (user_solution.size() != static_cast<size_t>(highs->getNumCol())) {
    highsLogUser(highs->getOptions().log_options, HighsLogType::kError,
                 "repairSolution: User solution size %d does not match model "
                 "number of columns %d\n",
                 int(user_solution.size()), int(highs->getNumCol()));
    return HighsStatus::kError;
  } else {
    // naive approach to get/check feasible solution
    // solve another MIP after fixing variables
    Highs clone;
    clone.setOptionValue("output_flag", false);
    clone.passModel(highs->getModel());

    HighsVarType vtype = HighsVarType::kContinuous;
    double tolerance = highs->getOptions().mip_feasibility_tolerance;

    // fix the variables
    for (HighsInt i = 0; i < static_cast<HighsInt>(user_solution.size()); i++) {
      if (user_solution[i] != kHighsUndefined) {
        double value = user_solution[i];
        highs->getColIntegrality(i, vtype);

        if (vtype == HighsVarType::kInteger) {
          // check if the provided solution value is integerfeasible
          if (!HighsIntegers::isIntegral(value, tolerance)) {
            highsLogUser(highs->getOptions().log_options, HighsLogType::kError,
                         "repairSolution: User solution (index %d) is outside "
                         "integrality tolerance (value %f)\n",
                         i, value);
            return HighsStatus::kError;
          } else {
            // since the variable is integral, round it to avoid numerical
            // issues
            value = std::round(value);
          }
        }

        clone.changeColBounds(i, value, value);
      }
    }

    // set callback to stop at first feasible solution
    bool user_interrupt = false;

    HighsCallbackFunctionType mip_callback =
        [&](int callback_type, const std::string& message,
            const HighsCallbackOutput* data_out, HighsCallbackInput* data_in,
            void* user_callback_data) {
          if (callback_type == kCallbackMipSolution) {
            user_interrupt = true;
          } else {
            data_in->user_interrupt = user_interrupt;
          }
        };

    clone.setCallback(mip_callback);
    clone.startCallback(kCallbackMipSolution);
    clone.startCallback(kCallbackMipInterrupt);
    clone.startCallback(kCallbackSimplexInterrupt);
    clone.startCallback(kCallbackIpmInterrupt);
    clone.run();

    auto status = clone.getModelStatus();

    // check if the solution is feasible
    if (status == HighsModelStatus::kOptimal ||
        status == HighsModelStatus::kInterrupt) {
      // copy the solution
      user_solution = clone.getSolution().col_value;
      user_has_solution = true;
      return HighsStatus::kOk;
    } else {
      highsLogUser(highs->getOptions().log_options, HighsLogType::kError,
                   "repairSolution: No feasible solution found\n");
      user_has_solution = false;
      return HighsStatus::kError;
    }
  }
}
