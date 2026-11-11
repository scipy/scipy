/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInfoDebug.cpp
 * @brief
 */
#include "lp_data/HighsInfoDebug.h"

HighsDebugStatus debugInfo(const HighsOptions& options, const HighsLp& lp,
                           const HighsBasis& basis,
                           const HighsSolution& solution, const HighsInfo& info,
                           const HighsModelStatus model_status) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  HighsDebugStatus call_status;

  const bool have_info = info.valid;
  const bool have_primal_solution = solution.value_valid;
  const bool have_dual_solution = solution.dual_valid;
  const bool have_basis = basis.valid;
  switch (model_status) {
    case HighsModelStatus::kNotset:
    case HighsModelStatus::kLoadError:
    case HighsModelStatus::kModelError:
    case HighsModelStatus::kPresolveError:
    case HighsModelStatus::kSolveError:
    case HighsModelStatus::kPostsolveError:
    case HighsModelStatus::kModelEmpty:
    case HighsModelStatus::kMemoryLimit:
      // Should have no info, so check this and return
      assert(!have_primal_solution);
      assert(!have_dual_solution);
      assert(!have_basis);
      call_status = debugNoInfo(info);
      if (call_status != HighsDebugStatus::kOk) return_status = call_status;
      return return_status;
    case HighsModelStatus::kOptimal:
    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kUnboundedOrInfeasible:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit:
    case HighsModelStatus::kSolutionLimit:
    case HighsModelStatus::kInterrupt:
    case HighsModelStatus::kHighsInterrupt:
    case HighsModelStatus::kUnknown:
      // Should have info
      assert(have_info == true);
      if (have_primal_solution) {
        if (info.num_primal_infeasibilities < 0) {
          highsLogDev(options.log_options, HighsLogType::kError,
                      "Have primal solution but num_primal_infeasibilities = "
                      "%" HIGHSINT_FORMAT "\n",
                      info.num_primal_infeasibilities);
          return HighsDebugStatus::kLogicalError;
        } else if (info.num_primal_infeasibilities == 0) {
          if (info.primal_solution_status != kSolutionStatusFeasible) {
            highsLogDev(options.log_options, HighsLogType::kError,
                        "Have primal solution and no infeasibilities but "
                        "primal status = %" HIGHSINT_FORMAT "\n",
                        info.primal_solution_status);
            return HighsDebugStatus::kLogicalError;
          }
        } else {
          if (info.primal_solution_status != kSolutionStatusInfeasible) {
            highsLogDev(options.log_options, HighsLogType::kError,
                        "Have primal solution and infeasibilities but primal "
                        "status = %" HIGHSINT_FORMAT "\n",
                        info.primal_solution_status);
            return HighsDebugStatus::kLogicalError;
          }
        }
      } else {
        if (info.primal_solution_status != kSolutionStatusNone) {
          highsLogDev(
              options.log_options, HighsLogType::kError,
              "Have no primal solution but primal status = %" HIGHSINT_FORMAT
              "\n",
              info.primal_solution_status);
          return HighsDebugStatus::kLogicalError;
        }
      }
      if (have_dual_solution) {
        if (info.num_dual_infeasibilities < 0) {
          highsLogDev(options.log_options, HighsLogType::kError,
                      "Have dual solution but num_dual_infeasibilities = "
                      "%" HIGHSINT_FORMAT "\n",
                      info.num_dual_infeasibilities);
          return HighsDebugStatus::kLogicalError;
        } else if (info.num_dual_infeasibilities == 0) {
          if (info.dual_solution_status != kSolutionStatusFeasible) {
            highsLogDev(options.log_options, HighsLogType::kError,
                        "Have dual solution and no infeasibilities but dual "
                        "status = %" HIGHSINT_FORMAT "\n",
                        info.dual_solution_status);
            return HighsDebugStatus::kLogicalError;
          }
        } else {
          if (info.dual_solution_status != kSolutionStatusInfeasible) {
            highsLogDev(options.log_options, HighsLogType::kError,
                        "Have dual solution and infeasibilities but dual "
                        "status = %" HIGHSINT_FORMAT "\n",
                        info.dual_solution_status);
            return HighsDebugStatus::kLogicalError;
          }
        }
      } else {
        if (info.dual_solution_status != kSolutionStatusNone) {
          highsLogDev(
              options.log_options, HighsLogType::kError,
              "Have no dual solution but dual status = %" HIGHSINT_FORMAT "\n",
              info.dual_solution_status);
          return HighsDebugStatus::kLogicalError;
        }
      }
      break;
    default:
      // All cases should have been considered so assert on reaching here
      assert(1 == 0);
  }

  return return_status;
}

HighsDebugStatus debugNoInfo(const HighsInfo& info) {
  HighsInfo no_info;
  no_info.invalidate();
  bool error_found = false;
  const std::vector<InfoRecord*>& info_records = info.records;
  const std::vector<InfoRecord*>& no_info_records = no_info.records;
  HighsInt num_info = info_records.size();
  for (HighsInt index = 0; index < num_info; index++) {
    HighsInfoType type = info_records[index]->type;
    if (type == HighsInfoType::kInt64) {
      int v0 = (int)*(((InfoRecordInt64*)info_records[index])[0].value);
      int v1 = (int)*(((InfoRecordInt64*)info_records[index])[0].value);
      if (v0 != v1)
        printf("debugNoInfo: Index %" HIGHSINT_FORMAT " has %d != %d \n", index,
               v0, v1);
      error_found = (*(((InfoRecordInt64*)info_records[index])[0].value) !=
                     *(((InfoRecordInt64*)no_info_records[index])[0].value)) ||
                    error_found;
    } else if (type == HighsInfoType::kInt) {
      int v0 = (int)*(((InfoRecordInt*)info_records[index])[0].value);
      int v1 = (int)*(((InfoRecordInt*)info_records[index])[0].value);
      if (v0 != v1)
        printf("debugNoInfo: Index %" HIGHSINT_FORMAT " has %d != %d \n", index,
               v0, v1);
      error_found = (*(((InfoRecordInt*)info_records[index])[0].value) !=
                     *(((InfoRecordInt*)no_info_records[index])[0].value)) ||
                    error_found;
    } else if (type == HighsInfoType::kDouble) {
      double v0 = (double)*(((InfoRecordDouble*)info_records[index])[0].value);
      double v1 = (double)*(((InfoRecordDouble*)info_records[index])[0].value);
      if (v0 != v1)
        printf("debugNoInfo: Index %" HIGHSINT_FORMAT " has %g != %g \n", index,
               v0, v1);
      error_found = (*(((InfoRecordDouble*)info_records[index])[0].value) !=
                     *(((InfoRecordDouble*)no_info_records[index])[0].value)) ||
                    error_found;
    } else {
      assert(1 == 0);
    }
  }
  error_found = (info.valid != no_info.valid) || error_found;
  if (error_found) return HighsDebugStatus::kLogicalError;
  return HighsDebugStatus::kOk;
}
