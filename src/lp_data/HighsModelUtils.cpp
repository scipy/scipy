/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "lp_data/HighsModelUtils.h"

#include <algorithm>
#include <vector>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "util/HighsUtils.h"

#ifdef HiGHSDEV
void analyseModelBounds(const char* message, int numBd,
                        const std::vector<double>& lower,
                        const std::vector<double>& upper) {
  if (numBd == 0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix = 0; ix < numBd; ix++) {
    if (highs_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Finite lower bound and infinite upper bound: Lb
        numLb++;
      } else {
        // Finite lower bound and   finite upper bound:
        if (lower[ix] < upper[ix]) {
          // Distinct finite bounds: Bx
          numBx++;
        } else {
          // Equal finite bounds: Fx
          numFx++;
        }
      }
    }
  }
  printf("Analysing %d %s bounds\n", numBd, message);
  if (numFr > 0)
    printf("   Free:  %7d (%3d%%)\n", numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    printf("   LB:    %7d (%3d%%)\n", numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    printf("   UB:    %7d (%3d%%)\n", numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    printf("   Boxed: %7d (%3d%%)\n", numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    printf("   Fixed: %7d (%3d%%)\n", numFx, (100 * numFx) / numBd);
  printf("grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  printf("grep_CharMl,%d,%d,%d,%d,%d,%d\n", numBd, numFr, numLb, numUb, numBx,
         numFx);
}

#endif
std::string ch4VarStatus(const HighsBasisStatus status, const double lower,
                         const double upper) {
  switch (status) {
    case HighsBasisStatus::LOWER:
      if (lower == upper) {
        return "FX";
      } else {
        return "LB";
      }
      break;
    case HighsBasisStatus::BASIC:
      return "BS";
      break;
    case HighsBasisStatus::UPPER:
      return "UB";
      break;
    case HighsBasisStatus::ZERO:
      return "FR";
      break;
    case HighsBasisStatus::SUPER:
      return "SU";
      break;
    case HighsBasisStatus::NONBASIC:
      return "NB";
      break;
  }
  return "";
}

void reportModelBoundSol(FILE* file, const bool columns, const int dim,
                         const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<std::string>& names,
                         const std::vector<double>& primal,
                         const std::vector<double>& dual,
                         const std::vector<HighsBasisStatus>& status) {
  const bool have_names = names.size() > 0;
  const bool have_basis = status.size() > 0;
  const bool have_primal = primal.size() > 0;
  const bool have_dual = dual.size() > 0;
  std::string ch4_var_status;
  if (columns) {
    fprintf(file, "Columns\n");
  } else {
    fprintf(file, "Rows\n");
  }
  fprintf(
      file,
      "    Index Status        Lower        Upper       Primal         Dual");
  if (have_names) {
    fprintf(file, "  Name\n");
  } else {
    fprintf(file, "\n");
  }
  for (int ix = 0; ix < dim; ix++) {
    if (have_basis) {
      ch4_var_status = ch4VarStatus(status[ix], lower[ix], upper[ix]);
    } else {
      ch4_var_status = "";
    }
    fprintf(file, "%9d   %4s %12g %12g", ix, ch4_var_status.c_str(), lower[ix],
            upper[ix]);
    if (have_primal) {
      fprintf(file, " %12g", primal[ix]);
    } else {
      fprintf(file, "             ");
    }
    if (have_dual) {
      fprintf(file, " %12g", dual[ix]);
    } else {
      fprintf(file, "             ");
    }
    if (have_names) {
      fprintf(file, "  %-s\n", names[ix].c_str());
    } else {
      fprintf(file, "\n");
    }
  }
}

bool namesWithSpaces(const int num_name, const std::vector<std::string>& names,
                     const bool report) {
  bool names_with_spaces = false;
  for (int ix = 0; ix < num_name; ix++) {
    int space_pos = names[ix].find(" ");
    if (space_pos >= 0) {
      if (report)
        printf("Name |%s| contains a space character in position %d\n",
               names[ix].c_str(), space_pos);
      names_with_spaces = true;
    }
  }
  return names_with_spaces;
}

int maxNameLength(const int num_name, const std::vector<std::string>& names) {
  int max_name_length = 0;
  for (int ix = 0; ix < num_name; ix++)
    max_name_length = std::max((int)names[ix].length(), max_name_length);
  return max_name_length;
}

HighsStatus normaliseNames(const HighsOptions& options,
                           const std::string name_type, const int num_name,
                           std::vector<std::string>& names,
                           int& max_name_length) {
  // Record the desired maximum name length
  int desired_max_name_length = max_name_length;
  // First look for empty names
  int num_empty_name = 0;
  std::string name_prefix = name_type.substr(0, 1);
  bool names_with_spaces = false;
  for (int ix = 0; ix < num_name; ix++) {
    if ((int)names[ix].length() == 0) num_empty_name++;
  }
  // If there are no empty names - in which case they will all be
  // replaced - find the maximum name length
  if (!num_empty_name) max_name_length = maxNameLength(num_name, names);
  bool construct_names =
      num_empty_name || max_name_length > desired_max_name_length;
  if (construct_names) {
    // Construct names, either because they are empty names, or
    // because the existing names are too long

    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "There are empty or excessively-long %s names: using "
                    "constructed names with prefix %s",
                    name_type.c_str(), name_prefix.c_str());
    for (int ix = 0; ix < num_name; ix++)
      names[ix] = name_prefix + std::to_string(ix);
  } else {
    // Using original names, so look to see whether there are names with spaces
    names_with_spaces = namesWithSpaces(num_name, names);
  }
  // Find the final maximum name length
  max_name_length = maxNameLength(num_name, names);
  // Can't have names with spaces and more than 8 characters
  if (max_name_length > 8 && names_with_spaces) return HighsStatus::Error;
  if (construct_names) return HighsStatus::Warning;
  return HighsStatus::OK;
}

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper) {
  HighsBasisStatus checked_status;
  if (ideal_status == HighsBasisStatus::LOWER ||
      ideal_status == HighsBasisStatus::ZERO) {
    // Looking to give status LOWER or ZERO
    if (highs_isInfinity(-lower)) {
      // Lower bound is infinite
      if (highs_isInfinity(upper)) {
        // Upper bound is infinite
        checked_status = HighsBasisStatus::ZERO;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::UPPER;
      }
    } else {
      checked_status = HighsBasisStatus::LOWER;
    }
  } else {
    // Looking to give status UPPER
    if (highs_isInfinity(upper)) {
      // Upper bound is infinite
      if (highs_isInfinity(-lower)) {
        // Lower bound is infinite
        checked_status = HighsBasisStatus::ZERO;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::LOWER;
      }
    } else {
      checked_status = HighsBasisStatus::UPPER;
    }
  }
  return checked_status;
}

// Return a string representation of PrimalDualStatus
std::string utilPrimalDualStatusToString(const int primal_dual_status) {
  switch (primal_dual_status) {
    case PrimalDualStatus::STATUS_NOTSET:
      return "Not set";
      break;
    case PrimalDualStatus::STATUS_NO_SOLUTION:
      return "No solution";
      break;
    case PrimalDualStatus::STATUS_UNKNOWN:
      return "Point of unknown feasibility";
      break;
    case PrimalDualStatus::STATUS_INFEASIBLE_POINT:
      return "Infeasible point";
      break;
    case PrimalDualStatus::STATUS_FEASIBLE_POINT:
      return "Feasible point";
      break;
    default:
#ifdef HiGHSDEV
      printf("Primal/dual status %d not recognised\n", primal_dual_status);
#endif
      return "Unrecognised primal/dual status";
      break;
  }
  return "";
}

// Return a string representation of HighsModelStatus.
std::string utilHighsModelStatusToString(const HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::NOTSET:
      return "Not Set";
      break;
    case HighsModelStatus::LOAD_ERROR:
      return "Load error";
      break;
    case HighsModelStatus::MODEL_ERROR:
      return "Model error";
      break;
    case HighsModelStatus::PRESOLVE_ERROR:
      return "Presolve error";
      break;
    case HighsModelStatus::SOLVE_ERROR:
      return "Solve error";
      break;
    case HighsModelStatus::POSTSOLVE_ERROR:
      return "Postsolve error";
      break;
    case HighsModelStatus::MODEL_EMPTY:
      return "Model empty";
      break;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return "Infeasible";  //"Primal infeasible";
      break;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return "Unbounded";  //"Primal unbounded";
      break;
    case HighsModelStatus::OPTIMAL:
      return "Optimal";
      break;
    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return "Reached dual objective upper bound";
      break;
    case HighsModelStatus::REACHED_TIME_LIMIT:
      return "Reached time limit";
      break;
    case HighsModelStatus::REACHED_ITERATION_LIMIT:
      return "Reached iteration limit";
      break;
    default:
#ifdef HiGHSDEV
      printf("HiGHS model status %d not recognised\n", (int)model_status);
#endif
      return "Unrecognised HiGHS model status";
      break;
  }
  return "";
}

void copyHighsIterationCounts(const HighsIterationCounts& iteration_counts,
                              HighsInfo& info) {
  info.simplex_iteration_count = iteration_counts.simplex;
  info.ipm_iteration_count = iteration_counts.ipm;
  info.crossover_iteration_count = iteration_counts.crossover;
}

void copyHighsIterationCounts(const HighsInfo& info,
                              HighsIterationCounts& iteration_counts) {
  iteration_counts.simplex = info.simplex_iteration_count;
  iteration_counts.ipm = info.ipm_iteration_count;
  iteration_counts.crossover = info.crossover_iteration_count;
}

// Deduce the HighsStatus value corresponding to a HighsModelStatus value.
HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::NOTSET:
      return HighsStatus::Error;
    case HighsModelStatus::LOAD_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::MODEL_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::PRESOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::SOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::POSTSOLVE_ERROR:
      return HighsStatus::Error;
    case HighsModelStatus::MODEL_EMPTY:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return HighsStatus::OK;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return HighsStatus::OK;
    case HighsModelStatus::OPTIMAL:
      return HighsStatus::OK;
    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return HighsStatus::OK;
    case HighsModelStatus::REACHED_TIME_LIMIT:
      return HighsStatus::Warning;
    case HighsModelStatus::REACHED_ITERATION_LIMIT:
      return HighsStatus::Warning;
    default:
      return HighsStatus::Error;
  }
}
