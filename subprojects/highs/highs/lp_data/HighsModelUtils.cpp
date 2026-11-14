/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "lp_data/HighsModelUtils.h"

#include <algorithm>
#include <cfloat>
#include <map>
#include <sstream>
#include <vector>

#include "lp_data/HighsSolution.h"
#include "util/stringutil.h"

void analyseModelBounds(const HighsLogOptions& log_options, const char* message,
                        HighsInt numBd, const std::vector<double>& lower,
                        const std::vector<double>& upper) {
  if (numBd == 0) return;
  HighsInt numFr = 0;
  HighsInt numLb = 0;
  HighsInt numUb = 0;
  HighsInt numBx = 0;
  HighsInt numFx = 0;
  for (HighsInt ix = 0; ix < numBd; ix++) {
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
  highsLogDev(log_options, HighsLogType::kInfo,
              "Analysing %" HIGHSINT_FORMAT " %s bounds\n", numBd, message);
  if (numFr > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Free:  %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   LB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   UB:    %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Boxed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    highsLogDev(log_options, HighsLogType::kInfo,
                "   Fixed: %7" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT "%%)\n",
                numFx, (100 * numFx) / numBd);
  highsLogDev(log_options, HighsLogType::kInfo,
              "grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  highsLogDev(log_options, HighsLogType::kInfo,
              "grep_CharMl,%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
              ",%" HIGHSINT_FORMAT "\n",
              numBd, numFr, numLb, numUb, numBx, numFx);
}

std::string statusToString(const HighsBasisStatus status, const double lower,
                           const double upper) {
  switch (status) {
    case HighsBasisStatus::kLower:
      if (lower == upper) {
        return "FX";
      } else {
        return "LB";
      }
      break;
    case HighsBasisStatus::kBasic:
      return "BS";
      break;
    case HighsBasisStatus::kUpper:
      if (lower == upper) {
        return "FX";
      } else {
        return "UB";
      }
      break;
    case HighsBasisStatus::kZero:
      return "FR";
      break;
    case HighsBasisStatus::kNonbasic:
      return "NB";
      break;
  }
  return "";
}

std::string typeToString(const HighsVarType type) {
  switch (type) {
    case HighsVarType::kContinuous:
      return "Continuous";
    case HighsVarType::kInteger:
      return "Integer   ";
    case HighsVarType::kSemiContinuous:
      return "Semi-conts";
    case HighsVarType::kSemiInteger:
      return "Semi-int  ";
    case HighsVarType::kImplicitInteger:
      return "ImpliedInt";
  }
  return "";
}

void writeModelBoundSolution(
    FILE* file, const HighsLogOptions& log_options, const bool columns,
    const HighsInt dim, const std::vector<double>& lower,
    const std::vector<double>& upper, const std::vector<std::string>& names,
    const bool have_primal, const std::vector<double>& primal,
    const bool have_dual, const std::vector<double>& dual,
    const bool have_basis, const std::vector<HighsBasisStatus>& status,
    const HighsVarType* integrality) {
  assert(names.size() == static_cast<size_t>(dim));
  if (have_primal) assert(primal.size() == static_cast<size_t>(dim));
  if (have_dual) assert(dual.size() == static_cast<size_t>(dim));
  if (have_basis) assert(status.size() == static_cast<size_t>(dim));
  const bool have_integrality = integrality != NULL;
  std::stringstream ss;
  std::string s = columns ? "Columns\n" : "Rows\n";
  highsFprintfString(file, log_options, s);
  ss.str(std::string());
  ss << "    Index Status        Lower        Upper       Primal         Dual";
  if (have_integrality) ss << "  Type      ";
  ss << "  Name\n";
  highsFprintfString(file, log_options, ss.str());
  for (HighsInt ix = 0; ix < dim; ix++) {
    ss.str(std::string());
    std::string var_status_string =
        have_basis ? statusToString(status[ix], lower[ix], upper[ix]) : "";
    ss << highsFormatToString("%9" HIGHSINT_FORMAT "   %4s %12g %12g", ix,
                              var_status_string.c_str(), lower[ix], upper[ix]);
    if (have_primal) {
      ss << highsFormatToString(" %12g", primal[ix]);
    } else {
      ss << "             ";
    }
    if (have_dual) {
      ss << highsFormatToString(" %12g", dual[ix]);
    } else {
      ss << "             ";
    }
    if (have_integrality)
      ss << highsFormatToString("  %s", typeToString(integrality[ix]).c_str());
    ss << highsFormatToString("  %-s\n", names[ix].c_str());
    highsFprintfString(file, log_options, ss.str());
  }
}

void writeModelObjective(FILE* file, const HighsLogOptions& log_options,
                         const HighsModel& model,
                         const std::vector<double>& primal_solution) {
  HighsCDouble objective_value =
      model.lp_.objectiveCDoubleValue(primal_solution);
  objective_value += model.hessian_.objectiveCDoubleValue(primal_solution);
  writeObjectiveValue(file, log_options, (double)objective_value);
}

void writeLpObjective(FILE* file, const HighsLogOptions& log_options,
                      const HighsLp& lp,
                      const std::vector<double>& primal_solution) {
  HighsCDouble objective_value = lp.objectiveCDoubleValue(primal_solution);
  writeObjectiveValue(file, log_options, (double)objective_value);
}

void writeObjectiveValue(FILE* file, const HighsLogOptions& log_options,
                         const double objective_value) {
  auto obj_string = highsDoubleToString(objective_value,
                                        kHighsSolutionValueToStringTolerance);
  std::string s = highsFormatToString("Objective %s\n", obj_string.data());
  highsFprintfString(file, log_options, s);
}

void writeMipBoundGap(FILE* file, const HighsLogOptions& log_options,
                      const double mip_dual_bound, const double mip_gap) {
  std::string s = highsFormatToString(
      "Objective_bound %s\n",
      highsDoubleToString(mip_dual_bound, kHighsSolutionValueToStringTolerance)
          .data());
  highsFprintfString(file, log_options, s);
  s = highsFormatToString(
      "MIP_rel_gap %s\n",
      highsDoubleToString(mip_gap, kHighsSolutionValueToStringTolerance)
          .data());
  highsFprintfString(file, log_options, s);
}

void writePrimalSolution(FILE* file, const HighsLogOptions& log_options,
                         const HighsLp& lp,
                         const std::vector<double>& primal_solution,
                         const bool sparse, const bool partial) {
  // Use when writing out the solution file (when names can be assumed
  // to exist) and the improving solution in the MIP solver (when
  // names cannot be assumed to exist)
  HighsInt num_nonzero_primal_value = 0;
  const bool have_col_names = lp.col_names_.size() > 0;
  if (have_col_names)
    assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  const bool sparse_or_partial = sparse || partial;
  assert(!(sparse && partial));
  if (sparse_or_partial) {
    // Determine the number of nonzero primal solution values
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      if (primal_solution[iCol]) num_nonzero_primal_value++;
  }
  // Indicate the number of column values to be written out, depending
  // on whether format is sparse: either lp.num_col_ if not sparse, or
  // the negation of the number of nonzero values, if sparse

  std::stringstream ss;
  ss.str(std::string());
  HighsInt num_col_field =
      sparse_or_partial ? -num_nonzero_primal_value : lp.num_col_;
  std::string sparse_or_partial_string =
      sparse_or_partial
          ? (" " + (sparse ? kHighsSparseString : kHighsPartialString))
          : "";
  ss << highsFormatToString("# Columns %d%s\n", int(num_col_field),
                            sparse_or_partial_string.c_str());
  highsFprintfString(file, log_options, ss.str());
  for (HighsInt ix = 0; ix < lp.num_col_; ix++) {
    if (sparse_or_partial && !primal_solution[ix]) continue;
    auto val_string = highsDoubleToString(primal_solution[ix],
                                          kHighsSolutionValueToStringTolerance);
    // Don't invent names locally: if none exist, then indicate this -
    // so that the (sparse/partial) solution line format remains "name
    // value (index)"
    const std::string name = have_col_names ? lp.col_names_[ix] : "NoName";
    ss.str(std::string());
    ss << highsFormatToString("%-s %s", name.c_str(), val_string.data());
    if (sparse_or_partial) ss << highsFormatToString(" %d", int(ix));
    ss << "\n";
    highsFprintfString(file, log_options, ss.str());
  }
  fflush(file);
}

void writeModelSolution(FILE* file, const HighsLogOptions& log_options,
                        const HighsModel& model,
                        const HighsModelStatus model_status,
                        const HighsSolution& solution, const HighsInfo& info,
                        const bool sparse, const bool partial) {
  const HighsLp& lp = model.lp_;
  const bool have_primal = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  assert(lp.row_names_.size() == static_cast<size_t>(lp.num_row_));
  if (have_primal) {
    assert((int)solution.col_value.size() >= lp.num_col_);
    assert((int)solution.row_value.size() >= lp.num_row_);
    assert(info.primal_solution_status != kSolutionStatusNone);
  }
  if (have_dual) {
    assert((int)solution.col_dual.size() >= lp.num_col_);
    assert((int)solution.row_dual.size() >= lp.num_row_);
    assert(info.dual_solution_status != kSolutionStatusNone);
  }
  std::stringstream ss;
  highsFprintfString(file, log_options, "Model status\n");
  ss.str(std::string());
  ss << highsFormatToString("%s\n",
                            utilModelStatusToString(model_status).c_str());
  highsFprintfString(file, log_options, ss.str());
  highsFprintfString(file, log_options, "\n# Primal solution values\n");
  if (!have_primal || info.primal_solution_status == kSolutionStatusNone) {
    highsFprintfString(file, log_options, "None\n");
  } else {
    if (info.primal_solution_status == kSolutionStatusFeasible) {
      highsFprintfString(file, log_options, "Feasible\n");
    } else {
      assert(info.primal_solution_status == kSolutionStatusInfeasible);
      highsFprintfString(file, log_options, "Infeasible\n");
    }
    writeModelObjective(file, log_options, model, solution.col_value);
    if (model.isMip() && model_status == HighsModelStatus::kOptimal &&
        info.mip_gap < kHighsInf)
      writeMipBoundGap(file, log_options, info.mip_dual_bound, info.mip_gap);
    writePrimalSolution(file, log_options, model.lp_, solution.col_value,
                        sparse, partial);
    if (sparse || partial) return;
    ss.str(std::string());
    ss << highsFormatToString("# Rows %" HIGHSINT_FORMAT "\n", lp.num_row_);
    highsFprintfString(file, log_options, ss.str());
    for (HighsInt ix = 0; ix < lp.num_row_; ix++) {
      auto val_string = highsDoubleToString(
          solution.row_value[ix], kHighsSolutionValueToStringTolerance);
      ss.str(std::string());
      ss << highsFormatToString("%-s %s\n", lp.row_names_[ix].c_str(),
                                val_string.data());
      highsFprintfString(file, log_options, ss.str());
    }
  }
  highsFprintfString(file, log_options, "\n# Dual solution values\n");
  if (!have_dual || info.dual_solution_status == kSolutionStatusNone) {
    highsFprintfString(file, log_options, "None\n");
  } else {
    if (info.dual_solution_status == kSolutionStatusFeasible) {
      highsFprintfString(file, log_options, "Feasible\n");
    } else {
      assert(info.dual_solution_status == kSolutionStatusInfeasible);
      highsFprintfString(file, log_options, "Infeasible\n");
    }
    ss.str(std::string());
    ss << highsFormatToString("# Columns %" HIGHSINT_FORMAT "\n", lp.num_col_);
    highsFprintfString(file, log_options, ss.str());
    for (HighsInt ix = 0; ix < lp.num_col_; ix++) {
      auto valStr = highsDoubleToString(solution.col_dual[ix],
                                        kHighsSolutionValueToStringTolerance);
      ss.str(std::string());
      ss << highsFormatToString("%-s %s\n", lp.col_names_[ix].c_str(),
                                valStr.data());
      highsFprintfString(file, log_options, ss.str());
    }
    ss.str(std::string());
    ss << highsFormatToString("# Rows %" HIGHSINT_FORMAT "\n", lp.num_row_);
    highsFprintfString(file, log_options, ss.str());
    for (HighsInt ix = 0; ix < lp.num_row_; ix++) {
      auto val_string = highsDoubleToString(
          solution.row_dual[ix], kHighsSolutionValueToStringTolerance);
      ss.str(std::string());
      ss << highsFormatToString("%-s %s\n", lp.row_names_[ix].c_str(),
                                val_string.data());
      highsFprintfString(file, log_options, ss.str());
    }
  }
}

bool replaceSpacesByUnderscores(std::string& name) {
  // Find the first occurrence of the substring
  const std::string replace_word = " ";
  const std::string replace_by = "_";
  size_t pos = name.find(replace_word);
  const bool found_space = pos != std::string::npos;
  const std::string from_name = name;

  // Iterate through the string and replace all
  // occurrences
  while (pos != string::npos) {
    // Replace the substring with the specified string
    name.replace(pos, replace_word.size(), replace_by);

    // Find the next occurrence of the substring
    pos = name.find(replace_word, pos + replace_by.size());
  }
  return found_space;
}

bool hasNamesWithSpaces(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (hasNamesWithSpaces(log_options, true, lp.col_names_)) return true;
  return hasNamesWithSpaces(log_options, false, lp.row_names_);
}

bool hasNamesWithSpaces(const HighsLogOptions& log_options, const bool col,
                        const std::vector<std::string>& names) {
  HighsInt num_names_with_spaces = 0;
  HighsInt num_name = names.size();
  for (HighsInt ix = 0; ix < num_name; ix++)
    if (names[ix].find(" ") != std::string::npos) num_names_with_spaces++;
  if (num_names_with_spaces)
    highsLogDev(log_options, HighsLogType::kInfo,
                "There are %d %s names with spaces\n",
                HighsInt(num_names_with_spaces), col ? "column" : "row");
  return num_names_with_spaces > 0;
}

bool hasIllegalNameForLpFile(const std::vector<std::string>& names) {
  HighsInt num_name = names.size();
  for (HighsInt ix = 0; ix < num_name; ix++) {
    const std::string& name = names[ix];
    const std::string first_character = name.substr(0, 1);
    if (name.find_first_not_of(kLegalLpFileColRowNameChar) != std::string::npos)
      return true;
    if (first_character.find_first_not_of(kLegalLpFileColRowNameFirstChar) !=
        std::string::npos)
      return true;
  }
  return false;
}

HighsInt maxNameLength(const HighsLp& lp) {
  return std::max(maxNameLength(lp.col_names_), maxNameLength(lp.row_names_));
}

HighsInt maxNameLength(const std::vector<std::string>& names) {
  HighsInt num_name = names.size();
  HighsInt max_name_length = 0;
  for (HighsInt ix = 0; ix < num_name; ix++)
    max_name_length = std::max(HighsInt(names[ix].length()), max_name_length);
  return max_name_length;
}

HighsStatus normaliseNames(const HighsLogOptions& log_options, HighsLp& lp,
                           HighsFileType type) {
  HighsStatus call_status =
      normaliseNames(log_options, true, lp.num_col_, lp.col_name_prefix_,
                     lp.col_name_suffix_, lp.col_names_, lp.col_hash_, type);
  assert(call_status != HighsStatus::kError);
  // Retain call_status in case normaliseNames for rows returns
  // HighsStatus::kOk
  HighsStatus return_status = call_status;
  call_status =
      normaliseNames(log_options, false, lp.num_row_, lp.row_name_prefix_,
                     lp.row_name_suffix_, lp.row_names_, lp.row_hash_, type);
  assert(call_status != HighsStatus::kError);
  if (call_status != HighsStatus::kOk) return call_status;
  return return_status;
}

HighsStatus normaliseNames(const HighsLogOptions& log_options, bool column,
                           HighsInt num_name_required, std::string& name_prefix,
                           HighsInt& name_suffix,
                           std::vector<std::string>& names,
                           HighsNameHash& name_hash, HighsFileType type) {
  // First look for no names being required
  if (num_name_required <= 0) return HighsStatus::kOk;
  // First look for there being no names
  HighsInt max_name_length = maxNameLength(names);
  // Keep track of whether the names are invalid
  bool invalid_names = max_name_length == 0;
  names.resize(num_name_required);
  HighsInt num_blank = 0;
  HighsInt num_names_with_spaces = 0;
  const HighsInt from_name_suffix = name_suffix;
  if (!invalid_names) {
    invalid_names =
        type == HighsFileType::kLp && hasIllegalNameForLpFile(names);
    if (!invalid_names) {
      for (HighsInt ix = 0; ix < num_name_required; ix++) {
        if (HighsInt(names[ix].length()) == 0) {
          // Name is blank, so create one
          num_blank++;
          name_prefix =
              column ? kHighsUniqueColNamePrefix : kHighsUniquerowNamePrefix;
          names[ix] = name_prefix + std::to_string(name_suffix++);
        } else if (type == HighsFileType::kMps) {
          // For MPS files, replace spaces in names by an underscore
          if (replaceSpacesByUnderscores(names[ix])) num_names_with_spaces++;
        }
      }
    }
  }
  // Possibly check for duplicates
  invalid_names = invalid_names || name_hash.hasDuplicate(names);
  name_hash.name2index.clear();
  if (!invalid_names) {
    // Names are valid
    if (num_blank || num_names_with_spaces) {
      // Report on any modifications
      if (num_blank) {
        highsLogUser(
            log_options, HighsLogType::kWarning,
            "Replaced %d blank %-6s name%s by one%s with prefix \"%s\", "
            "beginning with suffix %d\n",
            int(num_blank), column ? "column" : "row",
            num_blank == 1 ? "" : "s", num_blank == 1 ? "" : "s",
            name_prefix.c_str(), int(from_name_suffix));
      }
      if (num_names_with_spaces) {
        highsLogUser(log_options, HighsLogType::kWarning,
                     "Replaced spaces in %d %-6s name%s by underscores\n",
                     int(num_names_with_spaces), column ? "column" : "row",
                     num_names_with_spaces == 1 ? "" : "s");
      }
      return HighsStatus::kWarning;
    } else {
      // Return OK!
      return HighsStatus::kOk;
    }
  }
  assert(invalid_names);
  // Names are not valid, so create them using minimal prefix, and
  // suffix starting from zero
  name_prefix =
      column ? kHighsMinimalColNamePrefix : kHighsMinimalrowNamePrefix;
  name_suffix = 0;
  assert(names.size() == static_cast<size_t>(num_name_required));
  highsLogUser(log_options, HighsLogType::kWarning,
               "%s names are not present, or contain %sduplicates: using "
               "names with prefix \"%s\", beginning with suffix %d\n",
               column ? "Column" : "Row   ",
               type == HighsFileType::kLp ? "invalid characters or " : "",
               name_prefix.c_str(), int(name_suffix));
  for (HighsInt ix = 0; ix < num_name_required; ix++)
    names[ix] = name_prefix + std::to_string(name_suffix++);
  return HighsStatus::kWarning;
}

HighsFileType getFileType(const std::string& filename) {
  std::string lower_case_extension = getFilenameExt(filename);
  tolower(lower_case_extension);
  if (lower_case_extension.compare("mps") == 0) {
    return HighsFileType::kMps;
  } else if (lower_case_extension.compare("lp") == 0) {
    return HighsFileType::kLp;
  }
  return HighsFileType::kMinimal;
}

void writeSolutionFile(FILE* file, const HighsOptions& options,
                       const HighsModel& model, const HighsBasis& basis,
                       const HighsSolution& solution, const HighsInfo& info,
                       const HighsModelStatus model_status,
                       const HighsInt style) {
  const bool have_primal = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  const bool have_basis = basis.valid;
  const HighsLp& lp = model.lp_;
  const HighsLogOptions& log_options = options.log_options;
  assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  assert(lp.row_names_.size() == static_cast<size_t>(lp.num_row_));
  if (style == kSolutionStyleOldRaw) {
    writeOldRawSolution(file, log_options, lp, basis, solution);
  } else if (style == kSolutionStylePretty) {
    const HighsVarType* integrality =
        lp.integrality_.size() > 0 ? lp.integrality_.data() : nullptr;
    writeModelBoundSolution(file, log_options, true, lp.num_col_, lp.col_lower_,
                            lp.col_upper_, lp.col_names_, have_primal,
                            solution.col_value, have_dual, solution.col_dual,
                            have_basis, basis.col_status, integrality);
    writeModelBoundSolution(file, log_options, false, lp.num_row_,
                            lp.row_lower_, lp.row_upper_, lp.row_names_,
                            have_primal, solution.row_value, have_dual,
                            solution.row_dual, have_basis, basis.row_status);
    highsFprintfString(file, log_options, "\n");
    std::stringstream ss;
    ss.str(std::string());
    // Log the model status
    ss << highsFormatToString("Model status: %s\n\n",
                              utilModelStatusToString(model_status).c_str());
    highsFprintfString(file, log_options, ss.str());
    // Log the objective value
    ss.str(std::string());
    auto obj_string = highsDoubleToString((double)info.objective_function_value,
                                          kHighsSolutionValueToStringTolerance);
    ss << highsFormatToString("Objective value: %s\n", obj_string.data());
    highsFprintfString(file, log_options, ss.str());
    if (model.isMip() && model_status == HighsModelStatus::kOptimal &&
        info.mip_gap < kHighsInf) {
      // Possibly log the MIP bounds and gap
      auto mip_dual_bound_string = highsDoubleToString(
          info.mip_dual_bound, kHighsSolutionValueToStringTolerance);
      auto mip_gap_string = highsGapToString(info.mip_gap);
      ss.str(std::string());
      ss << highsFormatToString("Objective bound: %s\n",
                                mip_dual_bound_string.data());
      ss << highsFormatToString("MIP rel gap: %s%%\n", mip_gap_string.data());
      highsFprintfString(file, log_options, ss.str());
    }
  } else if (style == kSolutionStyleGlpsolRaw ||
             style == kSolutionStyleGlpsolPretty) {
    if (model.isQp()) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Cannot write QP solution in Glpsol format\n");
      return;
    }
    const bool raw = style == kSolutionStyleGlpsolRaw;
    writeGlpsolSolution(file, options, model.lp_, basis, solution, model_status,
                        info, raw);
  } else {
    // Standard raw solution file, possibly sparse or partial => only
    // nonzero primal values, but different flag
    const bool sparse = style == kSolutionStyleSparse;
    const bool partial = style == kSolutionStylePartial;
    assert(!(sparse && partial));
    assert(style == kSolutionStyleRaw || sparse || partial);
    writeModelSolution(file, log_options, model, model_status, solution, info,
                       sparse, partial);
  }
}

void writeGlpsolCostRow(FILE* file, const HighsLogOptions& log_options,
                        const bool raw, const bool is_mip,
                        const HighsInt row_id,
                        const std::string& objective_name,
                        const double objective_function_value) {
  std::stringstream ss;
  ss.str(std::string());
  if (raw) {
    double double_value = objective_function_value;
    auto double_string = highsDoubleToString(
        double_value, kGlpsolSolutionValueToStringTolerance);
    // Last term of 0 for dual should (also) be blank when not MIP
    ss << highsFormatToString("i %d %s%s%s\n", (int)row_id, is_mip ? "" : "b ",
                              double_string.data(), is_mip ? "" : " 0");
  } else {
    ss << highsFormatToString("%6d ", (int)row_id);
    if (objective_name.length() <= 12) {
      ss << highsFormatToString("%-12s ", objective_name.c_str());
    } else {
      ss << highsFormatToString("%s\n%20s", objective_name.c_str(), "");
    }
    if (is_mip) {
      ss << highsFormatToString("   ");
    } else {
      ss << highsFormatToString("B  ");
    }
    ss << highsFormatToString("%13.6g %13s %13s \n", objective_function_value,
                              "", "");
  }
  highsFprintfString(file, log_options, ss.str());
}

void writeGlpsolSolution(FILE* file, const HighsOptions& options,
                         const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution,
                         const HighsModelStatus model_status,
                         const HighsInfo& info, const bool raw) {
  const bool have_value = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  const bool have_basis = basis.valid;
  const double kGlpsolHighQuality = 1e-9;
  const double kGlpsolMediumQuality = 1e-6;
  const double kGlpsolLowQuality = 1e-3;
  const double kGlpsolPrintAsZero = 1e-9;
  const HighsLogOptions& log_options = options.log_options;
  assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  assert(lp.row_names_.size() == static_cast<size_t>(lp.num_row_));
  // Determine number of nonzeros including the objective function
  // and, hence, determine whether there is an objective function
  HighsInt num_nz = lp.numNz();
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    if (lp.col_cost_[iCol]) num_nz++;
  const bool empty_cost_row = num_nz == lp.numNz();
  const bool has_objective = !empty_cost_row;
  // Writes the solution using the GLPK raw style (defined in
  // api/wrsol.c) or pretty style (defined in api/prsol.c)
  //
  // When writing out the row information (and hence the number of
  // rows and nonzeros), the case of the cost row is tricky
  // (particularly if it's empty) if HiGHS is to be able to reproduce
  // the (inconsistent) behaviour of Glpsol.
  //
  // If Glpsol is run from a .mod file then the cost row is reported
  // unless there is no objective [minimize/maximize "objname"]
  // statement in the .mod file. In this case, the N-row in the MPS
  // file is called "R0000000" and referred to below as being artificial.
  //
  // However, the position of a defined cost row depends on where the
  // objective appears in the .mod file. If Glpsol is run from a .mod
  // file, and reads a .sol file, it must be in the right format.
  //
  // HiGHS can't read ..mod files, so works from an MPS or LP file
  // generated by glpsol.
  //
  // An MPS file generated by glpsol will have the cost row in the
  // same position as it was in the .mod file
  //
  // An LP file generated by glpsol will have the objective defined
  // first, so the desired position of the cost row in the .sol file
  // is unavailable. The only option with this route is to define the
  // cost row location "by hand" using glpsol_cost_row_location
  //
  // If Glpsol is run from an LP or MPS file then the cost row is not
  // reported. This behaviour is defined by setting
  // glpsol_cost_row_location = -1;
  //
  // This inconsistent behaviour means that it must be possible to
  // tell HiGHS to suppress the cost row
  //
  const HighsInt cost_row_option = options.glpsol_cost_row_location;
  // Define cost_row_location
  //
  // It is indexed from 1 so that it matches the index printed on that
  // row...
  //
  // ... hence a location of zero means that the cost row isn't
  // reported
  HighsInt cost_row_location = 0;
  std::string artificial_cost_row_name = "R0000000";
  const bool artificial_cost_row =
      lp.objective_name_ == artificial_cost_row_name;
  if (artificial_cost_row)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "The cost row name of \"%s\" is assumed to be artificial and "
                 "will not be reported in the Glpsol solution file\n",
                 lp.objective_name_.c_str());

  if (cost_row_option <= kGlpsolCostRowLocationLast ||
      cost_row_option > lp.num_row_) {
    // Place the cost row last
    cost_row_location = lp.num_row_ + 1;
  } else if (cost_row_option == kGlpsolCostRowLocationNone) {
    // Don't report the cost row
    assert(cost_row_location == 0);
  } else if (cost_row_option == kGlpsolCostRowLocationNoneIfEmpty) {
    // This option allows the cost row to be omitted if it's empty.
    if (empty_cost_row && artificial_cost_row) {
      // The cost row is empty and artificial, so don't report it
      assert(cost_row_location == 0);
    } else {
      // Place the cost row according to lp.cost_row_location_
      if (lp.cost_row_location_ >= 0) {
        // The cost row location is known from the MPS file. NB To
        // index from zero whenever possible, lp.cost_row_location_ =
        // 0 if the cost row came first
        assert(lp.cost_row_location_ <= lp.num_row_);
        cost_row_location = lp.cost_row_location_ + 1;
      } else {
        // The location isn't known from an MPS file, so place it
        // last, giving a warning
        cost_row_location = lp.num_row_ + 1;
        highsLogUser(
            options.log_options, HighsLogType::kWarning,
            "The cost row for the Glpsol solution file is reported last since "
            "there is no indication of where it should be\n");
      }
    }
  } else {
    // Place the cost row according to the option value
    cost_row_location = cost_row_option;
  }
  assert(0 <= cost_row_location && cost_row_location <= lp.num_row_ + 1);
  // Despite being written in C, GLPSOL indexes rows (columns) from
  // 1..m (1..n) with - bizarrely! - m being one more than the number
  // of constraints if the cost vector is reported.
  const HighsInt num_row = lp.num_row_;
  const HighsInt num_col = lp.num_col_;
  // There's one more row and more nonzeros if the cost row is
  // reported
  const HighsInt delta_num_row = cost_row_location > 0;
  const HighsInt glpsol_num_row = num_row + delta_num_row;
  // If the cost row isn't reported, then the number of nonzeros is
  // just the number in the constraint matrix
  if (cost_row_location <= 0) num_nz = lp.numNz();
  // Record the discrete nature of the model
  HighsInt num_integer = 0;
  HighsInt num_binary = 0;
  bool is_mip = false;
  if ((HighsInt)lp.integrality_.size() == lp.num_col_) {
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] != HighsVarType::kContinuous) {
        is_mip = true;
        num_integer++;
        if (lp.col_lower_[iCol] == 0 && lp.col_upper_[iCol] == 1) num_binary++;
      }
    }
  }
  // Raw and pretty are the initially the same, but for the "c "
  // prefix to raw lines
  std::string line_prefix = "";
  if (raw) line_prefix = "c ";
  highsFprintfString(file, log_options,
                     highsFormatToString("%s%-12s%s\n", line_prefix.c_str(),
                                         "Problem:", lp.model_name_.c_str()));
  highsFprintfString(file, log_options,
                     highsFormatToString("%s%-12s%d\n", line_prefix.c_str(),
                                         "Rows:", (int)glpsol_num_row));
  std::stringstream ss;
  ss.str(std::string());
  ss << highsFormatToString("%s%-12s%d", line_prefix.c_str(),
                            "Columns:", (int)num_col);
  if (!raw && is_mip)
    ss << highsFormatToString(" (%d integer, %d binary)", (int)num_integer,
                              (int)num_binary);
  ss << highsFormatToString("\n");
  highsFprintfString(file, log_options, ss.str());
  highsFprintfString(file, log_options,
                     highsFormatToString("%s%-12s%d\n", line_prefix.c_str(),
                                         "Non-zeros:", (int)num_nz));
  // Use model_status to define the GLPK model_status_text and
  // solution_status_char, where the former is used to specify the
  // model status. GLPK uses a single character to specify the
  // solution status, and for LPs this is deduced from the primal and
  // dual solution status. However, for MIPs, it is defined according
  // to the model status, so only set solution_status_char for MIPs
  std::string model_status_text = "???";
  std::string solution_status_char = "?";
  switch (model_status) {
    case HighsModelStatus::kOptimal:
      if (is_mip) {
        model_status_text = "INTEGER OPTIMAL";
        solution_status_char = "o";
      } else {
        model_status_text = "OPTIMAL";
      }
      break;
    case HighsModelStatus::kInfeasible:
      if (is_mip) {
        model_status_text = "INTEGER EMPTY";
        solution_status_char = "n";
      } else {
        model_status_text = "INFEASIBLE (FINAL)";
      }
      break;
    case HighsModelStatus::kUnbounded:
      // No apparent case in wrmip.c
      model_status_text = "UNBOUNDED";
      if (is_mip) solution_status_char = "u";
      break;
    default:
      if (info.primal_solution_status == kSolutionStatusFeasible) {
        if (is_mip) {
          model_status_text = "INTEGER NON-OPTIMAL";
          solution_status_char = "f";
        } else {
          model_status_text = "FEASIBLE";
        }
      } else {
        model_status_text = "UNDEFINED";
        if (is_mip) solution_status_char = "u";
      }
      break;
  }
  assert(model_status_text != "???");
  if (is_mip) assert(solution_status_char != "?");
  highsFprintfString(file, log_options,
                     highsFormatToString("%s%-12s%s\n", line_prefix.c_str(),
                                         "Status:", model_status_text.c_str()));
  // If info is not valid, then cannot write more
  if (!info.valid) return;
  // Now write out the numerical information
  //
  // Determine the objective name to write out
  std::string objective_name = lp.objective_name_;
  // There are row names to be written out, so there must be a
  // non-trivial objective name
  assert(lp.objective_name_ != "");
  const bool has_objective_name = lp.objective_name_ != "";
  highsFprintfString(
      file, log_options,
      highsFormatToString(
          "%s%-12s%s%.10g (%s)\n", line_prefix.c_str(), "Objective:",
          !(has_objective && has_objective_name)
              ? ""
              : (objective_name + " = ").c_str(),
          has_objective ? info.objective_function_value : 0,
          lp.sense_ == ObjSense::kMinimize ? "MINimum" : "MAXimum"));
  // No space after "c" on blank line!
  if (raw) line_prefix = "c";
  highsFprintfString(file, log_options,
                     highsFormatToString("%s\n", line_prefix.c_str()));
  // Detailed lines are rather different
  if (raw) {
    ss.str(std::string());
    ss << highsFormatToString("s %s %d %d ", is_mip ? "mip" : "bas",
                              (int)glpsol_num_row, (int)num_col);
    if (is_mip) {
      ss << highsFormatToString("%s", solution_status_char.c_str());
    } else {
      if (info.primal_solution_status == kSolutionStatusNone) {
        ss << highsFormatToString("u");
      } else if (info.primal_solution_status == kSolutionStatusInfeasible) {
        ss << highsFormatToString("i");
      } else if (info.primal_solution_status == kSolutionStatusFeasible) {
        ss << highsFormatToString("f");
      } else {
        ss << highsFormatToString("?");
      }
      ss << highsFormatToString(" ");
      if (info.dual_solution_status == kSolutionStatusNone) {
        ss << highsFormatToString("u");
      } else if (info.dual_solution_status == kSolutionStatusInfeasible) {
        ss << highsFormatToString("i");
      } else if (info.dual_solution_status == kSolutionStatusFeasible) {
        ss << highsFormatToString("f");
      } else {
        ss << highsFormatToString("?");
      }
    }
    double double_value = has_objective ? info.objective_function_value : 0;
    auto double_string =
        highsDoubleToString(double_value, kHighsSolutionValueToStringTolerance);
    ss << highsFormatToString(" %s\n", double_string.data());
    highsFprintfString(file, log_options, ss.str());
  }
  // GLPK puts out i 1 b 0 0 etc if there's no primal point, but
  // that's meaningless at best, so HiGHS returns in that case
  if (!have_value) return;
  if (!raw) {
    ss.str(std::string());
    ss << highsFormatToString(
        "   No.   Row name   %s   Activity     Lower bound  "
        " Upper bound",
        have_basis ? "St" : "  ");
    if (have_dual) ss << highsFormatToString("    Marginal");
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
    ss.str(std::string());
    ss << highsFormatToString(
        "------ ------------ %s ------------- ------------- "
        "-------------",
        have_basis ? "--" : "  ");
    if (have_dual) ss << highsFormatToString(" -------------");
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }

  HighsInt row_id = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    row_id++;
    if (row_id == cost_row_location) {
      writeGlpsolCostRow(file, log_options, raw, is_mip, row_id, objective_name,
                         info.objective_function_value);
      row_id++;
    }
    ss.str(std::string());
    if (raw) {
      ss << highsFormatToString("i %d ", (int)row_id);
      if (is_mip) {
        // Complete the line if for a MIP
        double double_value = have_value ? solution.row_value[iRow] : 0;
        auto double_string = highsDoubleToString(
            double_value, kHighsSolutionValueToStringTolerance);
        ss << highsFormatToString("%s\n", double_string.data());
        highsFprintfString(file, log_options, ss.str());
        continue;
      }
    } else {
      ss << highsFormatToString("%6d ", (int)row_id);
      std::string row_name = lp.row_names_[iRow];
      if (row_name.length() <= 12) {
        ss << highsFormatToString("%-12s ", row_name.c_str());
      } else {
        ss << highsFormatToString("%s\n", row_name.c_str());
        highsFprintfString(file, log_options, ss.str());
        ss.str(std::string());
        ss << "                    ";
      }
    }
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    const double value = have_value ? solution.row_value[iRow] : 0;
    const double dual = have_dual ? solution.row_dual[iRow] : 0;
    std::string status_text = "  ";
    std::string status_char = "";
    if (have_basis) {
      switch (basis.row_status[iRow]) {
        case HighsBasisStatus::kBasic:
          status_text = "B ";
          status_char = "b";
          break;
        case HighsBasisStatus::kLower:
          status_text = lower == upper ? "NS" : "NL";
          status_char = lower == upper ? "s" : "l";
          break;
        case HighsBasisStatus::kUpper:
          status_text = lower == upper ? "NS" : "NU";
          status_char = lower == upper ? "s" : "u";
          break;
        case HighsBasisStatus::kZero:
          status_text = "NF";
          status_char = "f";
          break;
        default:
          status_text = "??";
          status_char = "?";
          break;
      }
    }
    if (raw) {
      ss << highsFormatToString("%s ", status_char.c_str());
      double double_value = have_value ? solution.row_value[iRow] : 0;
      auto double_string = highsDoubleToString(
          double_value, kHighsSolutionValueToStringTolerance);
      ss << highsFormatToString("%s ", double_string.data());
    } else {
      ss << highsFormatToString("%s ", status_text.c_str());
      ss << highsFormatToString(
          "%13.6g ", fabs(value) <= kGlpsolPrintAsZero ? 0.0 : value);
      if (lower > -kHighsInf)
        ss << highsFormatToString("%13.6g ", lower);
      else
        ss << highsFormatToString("%13s ", "");
      if (lower != upper && upper < kHighsInf)
        ss << highsFormatToString("%13.6g ", upper);
      else
        ss << highsFormatToString("%13s ", lower == upper ? "=" : "");
    }
    if (have_dual) {
      if (raw) {
        double double_value = solution.row_dual[iRow];
        auto double_string = highsDoubleToString(
            double_value, kHighsSolutionValueToStringTolerance);
        ss << highsFormatToString("%s", double_string.data());
      } else {
        // If the row is known to be basic, don't print the dual
        // value. If there's no basis, row cannot be known to be basic
        bool not_basic = have_basis;
        if (have_basis)
          not_basic = basis.row_status[iRow] != HighsBasisStatus::kBasic;
        if (not_basic) {
          if (fabs(dual) <= kGlpsolPrintAsZero)
            ss << highsFormatToString("%13s", "< eps");
          else
            ss << highsFormatToString("%13.6g ", dual);
        }
      }
    }
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }

  if (cost_row_location == lp.num_row_ + 1) {
    row_id++;
    writeGlpsolCostRow(file, log_options, raw, is_mip, row_id, objective_name,
                       info.objective_function_value);
  }
  if (!raw) highsFprintfString(file, log_options, "\n");

  if (!raw) {
    ss.str(std::string());
    ss << highsFormatToString(
        "   No. Column name  %s   Activity     Lower bound  "
        " Upper bound",
        have_basis ? "St" : "  ");
    if (have_dual) ss << highsFormatToString("    Marginal");
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
    ss.str(std::string());
    ss << highsFormatToString(
        "------ ------------ %s ------------- ------------- "
        "-------------",
        have_basis ? "--" : "  ");
    if (have_dual) ss << highsFormatToString(" -------------");
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }

  if (raw) line_prefix = "j ";
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    ss.str(std::string());
    if (raw) {
      ss << highsFormatToString("%s%d ", line_prefix.c_str(), (int)(iCol + 1));
      if (is_mip) {
        double double_value = have_value ? solution.col_value[iCol] : 0;
        auto double_string = highsDoubleToString(
            double_value, kHighsSolutionValueToStringTolerance);
        ss << highsFormatToString("%s\n", double_string.data());
        highsFprintfString(file, log_options, ss.str());
        continue;
      }
    } else {
      ss << highsFormatToString("%6d ", (int)(iCol + 1));
      std::string col_name = lp.col_names_[iCol];
      if (col_name.length() <= 12) {
        ss << highsFormatToString("%-12s ", col_name.c_str());
      } else {
        ss << highsFormatToString("%s\n", col_name.c_str());
        highsFprintfString(file, log_options, ss.str());
        ss.str(std::string());
        ss << "                    ";
      }
    }
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
    const double value = have_value ? solution.col_value[iCol] : 0;
    const double dual = have_dual ? solution.col_dual[iCol] : 0;
    std::string status_text = "  ";
    std::string status_char = "";
    if (have_basis) {
      switch (basis.col_status[iCol]) {
        case HighsBasisStatus::kBasic:
          status_text = "B ";
          status_char = "b";
          break;
        case HighsBasisStatus::kLower:
          status_text = lower == upper ? "NS" : "NL";
          status_char = lower == upper ? "s" : "l";
          break;
        case HighsBasisStatus::kUpper:
          status_text = lower == upper ? "NS" : "NU";
          status_char = lower == upper ? "s" : "u";
          break;
        case HighsBasisStatus::kZero:
          status_text = "NF";
          status_char = "f";
          break;
        default:
          status_text = "??";
          status_char = "?";
          break;
      }
    } else if (is_mip) {
      if (lp.integrality_[iCol] != HighsVarType::kContinuous)
        status_text = "* ";
    }
    if (raw) {
      ss << highsFormatToString("%s ", status_char.c_str());
      double double_value = have_value ? solution.col_value[iCol] : 0;
      auto double_string = highsDoubleToString(
          double_value, kHighsSolutionValueToStringTolerance);
      ss << highsFormatToString("%s ", double_string.data());
    } else {
      ss << highsFormatToString("%s ", status_text.c_str());
      ss << highsFormatToString(
          "%13.6g ", fabs(value) <= kGlpsolPrintAsZero ? 0.0 : value);
      if (lower > -kHighsInf)
        ss << highsFormatToString("%13.6g ", lower);
      else
        ss << highsFormatToString("%13s ", "");
      if (lower != upper && upper < kHighsInf)
        ss << highsFormatToString("%13.6g ", upper);
      else
        ss << highsFormatToString("%13s ", lower == upper ? "=" : "");
    }
    if (have_dual) {
      if (raw) {
        double double_value = solution.col_dual[iCol];
        auto double_string = highsDoubleToString(
            double_value, kHighsSolutionValueToStringTolerance);
        ss << highsFormatToString("%s", double_string.data());
      } else {
        // If the column is known to be basic, don't print the dual
        // value. If there's no basis, column cannot be known to be
        // basic
        bool not_basic = have_basis;
        if (have_basis)
          not_basic = basis.col_status[iCol] != HighsBasisStatus::kBasic;
        if (not_basic) {
          if (fabs(dual) <= kGlpsolPrintAsZero)
            ss << highsFormatToString("%13s", "< eps");
          else
            ss << highsFormatToString("%13.6g ", dual);
        }
      }
    }
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }
  if (raw) {
    highsFprintfString(file, log_options, "e o f\n");
    return;
  }
  HighsPrimalDualErrors errors;
  HighsInfo local_info;
  HighsInt absolute_error_index;
  double absolute_error_value;
  HighsInt relative_error_index;
  double relative_error_value;
  getLpKktFailures(options, lp, solution, basis, local_info, &errors, true);
  highsFprintfString(file, log_options, "\n");
  if (is_mip) {
    highsFprintfString(file, log_options, "Integer feasibility conditions:\n");
  } else {
    highsFprintfString(file, log_options,
                       "Karush-Kuhn-Tucker optimality conditions:\n");
  }
  highsFprintfString(file, log_options, "\n");
  // Primal residual
  absolute_error_value = errors.glpsol_max_primal_residual.absolute_value;
  absolute_error_index = errors.glpsol_max_primal_residual.absolute_index + 1;
  relative_error_value = errors.glpsol_max_primal_residual.relative_value;
  relative_error_index = errors.glpsol_max_primal_residual.relative_index + 1;
  if (!absolute_error_value) absolute_error_index = 0;
  if (!relative_error_value) relative_error_index = 0;
  ss.str(std::string());
  ss << highsFormatToString(
      "KKT.PE: max.abs.err = %.2e on row %d\n", absolute_error_value,
      absolute_error_index == 0 ? 0 : (int)absolute_error_index);
  highsFprintfString(file, log_options, ss.str());
  ss.str(std::string());
  ss << highsFormatToString(
      "        max.rel.err = %.2e on row %d\n", relative_error_value,
      absolute_error_index == 0 ? 0 : (int)relative_error_index);
  highsFprintfString(file, log_options, ss.str());
  ss.str(std::string());
  ss << highsFormatToString(
      "%8s%s\n", "",
      relative_error_value <= kGlpsolHighQuality     ? "High quality"
      : relative_error_value <= kGlpsolMediumQuality ? "Medium quality"
      : relative_error_value <= kGlpsolLowQuality    ? "Low quality"
                                                  : "PRIMAL SOLUTION IS WRONG");
  ss << "\n";
  highsFprintfString(file, log_options, ss.str());

  // Primal infeasibility
  absolute_error_value = errors.glpsol_max_primal_infeasibility.absolute_value;
  absolute_error_index =
      errors.glpsol_max_primal_infeasibility.absolute_index + 1;
  relative_error_value = errors.glpsol_max_primal_infeasibility.relative_value;
  relative_error_index =
      errors.glpsol_max_primal_infeasibility.relative_index + 1;
  if (!absolute_error_value) absolute_error_index = 0;
  if (!relative_error_value) relative_error_index = 0;
  bool on_col = absolute_error_index > 0 && absolute_error_index <= lp.num_col_;
  ss.str(std::string());
  ss << highsFormatToString("KKT.PB: max.abs.err = %.2e on %s %d\n",
                            absolute_error_value, on_col ? "column" : "row",
                            absolute_error_index <= lp.num_col_
                                ? (int)absolute_error_index
                                : (int)(absolute_error_index - lp.num_col_));
  highsFprintfString(file, log_options, ss.str());
  on_col = relative_error_index > 0 && relative_error_index <= lp.num_col_;
  ss.str(std::string());
  ss << highsFormatToString("        max.rel.err = %.2e on %s %d\n",
                            relative_error_value, on_col ? "column" : "row",
                            relative_error_index <= lp.num_col_
                                ? (int)relative_error_index
                                : (int)(relative_error_index - lp.num_col_));
  highsFprintfString(file, log_options, ss.str());
  ss.str(std::string());
  ss << highsFormatToString(
      "%8s%s\n", "",
      relative_error_value <= kGlpsolHighQuality     ? "High quality"
      : relative_error_value <= kGlpsolMediumQuality ? "Medium quality"
      : relative_error_value <= kGlpsolLowQuality
          ? "Low quality"
          : "PRIMAL SOLUTION IS INFEASIBLE");
  ss << "\n";
  highsFprintfString(file, log_options, ss.str());

  if (have_dual) {
    // Dual residual
    absolute_error_value = errors.glpsol_max_dual_residual.absolute_value;
    absolute_error_index = errors.glpsol_max_dual_residual.absolute_index + 1;
    relative_error_value = errors.glpsol_max_dual_residual.relative_value;
    relative_error_index = errors.glpsol_max_dual_residual.relative_index + 1;
    if (!absolute_error_value) absolute_error_index = 0;
    if (!relative_error_value) relative_error_index = 0;
    ss.str(std::string());
    ss << highsFormatToString("KKT.DE: max.abs.err = %.2e on column %d\n",
                              absolute_error_value, (int)absolute_error_index);
    ss << highsFormatToString("        max.rel.err = %.2e on column %d\n",
                              relative_error_value, (int)relative_error_index);
    ss << highsFormatToString(
        "%8s%s\n", "",
        relative_error_value <= kGlpsolHighQuality     ? "High quality"
        : relative_error_value <= kGlpsolMediumQuality ? "Medium quality"
        : relative_error_value <= kGlpsolLowQuality    ? "Low quality"
                                                    : "DUAL SOLUTION IS WRONG");
    ss << "\n";
    highsFprintfString(file, log_options, ss.str());

    // Dual infeasibility
    absolute_error_value = errors.glpsol_max_dual_infeasibility.absolute_value;
    absolute_error_index =
        errors.glpsol_max_dual_infeasibility.absolute_index + 1;
    relative_error_value = errors.glpsol_max_dual_infeasibility.relative_value;
    relative_error_index =
        errors.glpsol_max_dual_infeasibility.relative_index + 1;
    if (!absolute_error_value) absolute_error_index = 0;
    if (!relative_error_value) relative_error_index = 0;
    bool on_col =
        absolute_error_index > 0 && absolute_error_index <= lp.num_col_;
    ss.str(std::string());
    ss << highsFormatToString("KKT.DB: max.abs.err = %.2e on %s %d\n",
                              absolute_error_value, on_col ? "column" : "row",
                              absolute_error_index <= lp.num_col_
                                  ? (int)absolute_error_index
                                  : (int)(absolute_error_index - lp.num_col_));
    highsFprintfString(file, log_options, ss.str());
    on_col = relative_error_index > 0 && relative_error_index <= lp.num_col_;
    ss.str(std::string());
    ss << highsFormatToString("        max.rel.err = %.2e on %s %d\n",
                              relative_error_value, on_col ? "column" : "row",
                              relative_error_index <= lp.num_col_
                                  ? (int)relative_error_index
                                  : (int)(relative_error_index - lp.num_col_));
    highsFprintfString(file, log_options, ss.str());
    ss.str(std::string());
    ss << highsFormatToString(
        "%8s%s\n", "",
        relative_error_value <= kGlpsolHighQuality     ? "High quality"
        : relative_error_value <= kGlpsolMediumQuality ? "Medium quality"
        : relative_error_value <= kGlpsolLowQuality
            ? "Low quality"
            : "DUAL SOLUTION IS INFEASIBLE");
    ss << "\n";
    highsFprintfString(file, log_options, ss.str());
  }
  highsFprintfString(file, log_options, "End of output\n");
}

void writeOldRawSolution(FILE* file, const HighsLogOptions& log_options,
                         const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution) {
  const bool have_value = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  const bool have_basis = basis.valid;
  vector<double> use_col_value;
  vector<double> use_row_value;
  vector<double> use_col_dual;
  vector<double> use_row_dual;
  vector<HighsBasisStatus> use_col_status;
  vector<HighsBasisStatus> use_row_status;
  if (have_value) {
    use_col_value = solution.col_value;
    use_row_value = solution.row_value;
  }
  if (have_dual) {
    use_col_dual = solution.col_dual;
    use_row_dual = solution.row_dual;
  }
  if (have_basis) {
    use_col_status = basis.col_status;
    use_row_status = basis.row_status;
  }
  if (!have_value && !have_dual && !have_basis) return;
  highsFprintfString(
      file, log_options,
      highsFormatToString(
          "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT
          " : Number of columns and rows for primal or dual solution "
          "or basis\n",
          lp.num_col_, lp.num_row_));
  std::stringstream ss;
  ss.str(std::string());
  if (have_value) {
    ss << highsFormatToString("T");
  } else {
    ss << highsFormatToString("F");
  }
  ss << highsFormatToString(" Primal solution\n");
  highsFprintfString(file, log_options, ss.str());
  ss.str(std::string());
  if (have_dual) {
    ss << highsFormatToString("T");
  } else {
    ss << highsFormatToString("F");
  }
  ss << highsFormatToString(" Dual solution\n");
  highsFprintfString(file, log_options, ss.str());
  ss.str(std::string());
  if (have_basis) {
    ss << highsFormatToString("T");
  } else {
    ss << highsFormatToString("F");
  }
  ss << highsFormatToString(" Basis\n");
  highsFprintfString(file, log_options, ss.str());
  highsFprintfString(file, log_options, "Columns\n");
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    ss.str(std::string());
    if (have_value) ss << highsFormatToString("%.15g ", use_col_value[iCol]);
    if (have_dual) ss << highsFormatToString("%.15g ", use_col_dual[iCol]);
    if (have_basis)
      ss << highsFormatToString("%" HIGHSINT_FORMAT "",
                                (HighsInt)use_col_status[iCol]);
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }
  highsFprintfString(file, log_options, "Rows\n");
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    ss.str(std::string());
    if (have_value) ss << highsFormatToString("%.15g ", use_row_value[iRow]);
    if (have_dual) ss << highsFormatToString("%.15g ", use_row_dual[iRow]);
    if (have_basis)
      ss << highsFormatToString("%" HIGHSINT_FORMAT "",
                                (HighsInt)use_row_status[iRow]);
    ss << highsFormatToString("\n");
    highsFprintfString(file, log_options, ss.str());
  }
}

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper) {
  HighsBasisStatus checked_status;
  if (ideal_status == HighsBasisStatus::kLower ||
      ideal_status == HighsBasisStatus::kZero) {
    // Looking to give status LOWER or ZERO
    if (highs_isInfinity(-lower)) {
      // Lower bound is infinite
      if (highs_isInfinity(upper)) {
        // Upper bound is infinite
        checked_status = HighsBasisStatus::kZero;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::kUpper;
      }
    } else {
      checked_status = HighsBasisStatus::kLower;
    }
  } else {
    // Looking to give status UPPER
    if (highs_isInfinity(upper)) {
      // Upper bound is infinite
      if (highs_isInfinity(-lower)) {
        // Lower bound is infinite
        checked_status = HighsBasisStatus::kZero;
      } else {
        // Upper bound is finite
        checked_status = HighsBasisStatus::kLower;
      }
    } else {
      checked_status = HighsBasisStatus::kUpper;
    }
  }
  return checked_status;
}

// Return a string representation of SolutionStatus
std::string utilSolutionStatusToString(const HighsInt solution_status) {
  switch (solution_status) {
    case kSolutionStatusNone:
      return "None";
      break;
    case kSolutionStatusInfeasible:
      return "Infeasible";
      break;
    case kSolutionStatusFeasible:
      return "Feasible";
      break;
    default:
      assert(1 == 0);
      return "Unrecognised solution status";
  }
}

// Return a string representation of HighsBasisStatus
StatusString utilBasisStatusToString(const HighsBasisStatus basis_status) {
  StatusString status;
  switch (basis_status) {
    case HighsBasisStatus::kLower:
      status.full_ = "At lower/fixed bound";
      status.s2_ = "LO";
      break;
    case HighsBasisStatus::kBasic:
      status.full_ = "Basic";
      status.s2_ = "BS";
      break;
    case HighsBasisStatus::kUpper:
      status.full_ = "At upper bound";
      status.s2_ = "UP";
      break;
    case HighsBasisStatus::kZero:
      status.full_ = "Free at zero";
      status.s2_ = "ZE";
      break;
    case HighsBasisStatus::kNonbasic:
      status.full_ = "Nonbasic";
      status.s2_ = "NB";
      break;
    default:
      assert(1 == 0);
      status.full_ = "Unrecognised solution status";
      status.s2_ = "??";
  }
  return status;
}

// Return a string representation of basis validity
std::string utilBasisValidityToString(const HighsInt basis_validity) {
  if (basis_validity) {
    return "Valid";
  } else {
    return "Not valid";
  }
}

// Return a string representation of HighsModelStatus.
std::string utilModelStatusToString(const HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::kNotset:
      return "Not Set";
    case HighsModelStatus::kLoadError:
      return "Load error";
    case HighsModelStatus::kModelError:
      return "Model error";
    case HighsModelStatus::kPresolveError:
      return "Presolve error";
    case HighsModelStatus::kSolveError:
      return "Solve error";
    case HighsModelStatus::kPostsolveError:
      return "Postsolve error";
    case HighsModelStatus::kModelEmpty:
      return "Empty";
    case HighsModelStatus::kMemoryLimit:
      return "Memory limit reached";
    case HighsModelStatus::kOptimal:
      return "Optimal";
    case HighsModelStatus::kInfeasible:
      return "Infeasible";
    case HighsModelStatus::kUnboundedOrInfeasible:
      return "Primal infeasible or unbounded";
    case HighsModelStatus::kUnbounded:
      return "Unbounded";
    case HighsModelStatus::kObjectiveBound:
      return "Bound on objective reached";
    case HighsModelStatus::kObjectiveTarget:
      return "Target for objective reached";
    case HighsModelStatus::kTimeLimit:
      return "Time limit reached";
    case HighsModelStatus::kIterationLimit:
      return "Iteration limit reached";
    case HighsModelStatus::kSolutionLimit:
      return "Solution limit reached";
    case HighsModelStatus::kInterrupt:
      return "Interrupted by user";
    case HighsModelStatus::kHighsInterrupt:
      return "Interrupted by HiGHS";
    case HighsModelStatus::kUnknown:
      return "Unknown";
    default:
      assert(1 == 0);
      return "Unrecognised HiGHS model status";
  }
}

std::string utilPresolveRuleTypeToString(const HighsInt rule_type) {
  if (rule_type == kPresolveRuleEmptyRow) {
    return "Empty row";
  } else if (rule_type == kPresolveRuleSingletonRow) {
    return "Singleton row";
  } else if (rule_type == kPresolveRuleRedundantRow) {
    return "Redundant row";
  } else if (rule_type == kPresolveRuleEmptyCol) {
    return "Empty column";
  } else if (rule_type == kPresolveRuleFixedCol) {
    return "Fixed column";
  } else if (rule_type == kPresolveRuleDominatedCol) {
    return "Dominated column";
  } else if (rule_type == kPresolveRuleForcingRow) {
    return "Forcing row";
  } else if (rule_type == kPresolveRuleForcingCol) {
    return "Forcing column";
  } else if (rule_type == kPresolveRuleFreeColSubstitution) {
    return "Free column substitution";
  } else if (rule_type == kPresolveRuleDoubletonEquation) {
    return "Doubleton equation";
  } else if (rule_type == kPresolveRuleDependentEquations) {
    return "Dependent equations";
  } else if (rule_type == kPresolveRuleDependentFreeCols) {
    return "Dependent free columns";
  } else if (rule_type == kPresolveRuleAggregator) {
    return "Aggregator";
  } else if (rule_type == kPresolveRuleParallelRowsAndCols) {
    return "Parallel rows and columns";
  } else if (rule_type == kPresolveRuleSparsify) {
    return "Sparsify";
  } else if (rule_type == kPresolveRuleProbing) {
    return "Probing";
  } else if (rule_type == kPresolveRuleEnumeration) {
    return "Enumeration";
  } else if (rule_type == kPresolveRuleDualFixing) {
    return "Dual fixing";
  } else if (rule_type == kPresolveRuleZeroCostSingleton) {
    return "Zero cost singleton";
  } else if (rule_type == kPresolveRuleColStuffing) {
    return "Column stuffing";
  } else if (rule_type == kPresolveRuleDualFixProbing) {
    return "Dual-fixing probing";
  } else if (rule_type == kPresolveRuleInitialSweep) {
    return "Initial sweep";
  } else if (rule_type == kPresolveRuleFourierMotzkin) {
    return "Fourier-Motzkin";
  } else if (rule_type == kPresolveRuleWeaklyDominatedCol) {
    return "Weakly dominated column";
  }
  assert(1 == 0);
  return "????";
}

std::string utilPostsolveStatusToString(
    const HighsPostsolveStatus postsolve_status) {
  switch (postsolve_status) {
    case HighsPostsolveStatus::kNotPresolved:
      return "Not presolved";
    case HighsPostsolveStatus::kNoPrimalSolutionError:
      return "Primal solution error";
    case HighsPostsolveStatus::kSolutionRecovered:
      return "Solution recovered";
    case HighsPostsolveStatus::kBasisError:
      return "Basis error";
    default:
      assert(1 == 0);
      return "Unrecognised HiGHS postsolve status";
  }
}

// Deduce the HighsStatus value corresponding to a HighsModelStatus value.
HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status) {
  switch (model_status) {
    case HighsModelStatus::kNotset:
      return HighsStatus::kError;
    case HighsModelStatus::kLoadError:
      return HighsStatus::kError;
    case HighsModelStatus::kModelError:
      return HighsStatus::kError;
    case HighsModelStatus::kPresolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kSolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kPostsolveError:
      return HighsStatus::kError;
    case HighsModelStatus::kMemoryLimit:
      return HighsStatus::kError;
    case HighsModelStatus::kModelEmpty:
      return HighsStatus::kOk;
    case HighsModelStatus::kOptimal:
      return HighsStatus::kOk;
    case HighsModelStatus::kInfeasible:
      return HighsStatus::kOk;
    case HighsModelStatus::kUnboundedOrInfeasible:
      return HighsStatus::kOk;
    case HighsModelStatus::kUnbounded:
      return HighsStatus::kOk;
    case HighsModelStatus::kObjectiveBound:
      return HighsStatus::kOk;
    case HighsModelStatus::kObjectiveTarget:
      return HighsStatus::kOk;
    case HighsModelStatus::kTimeLimit:
      return HighsStatus::kWarning;
    case HighsModelStatus::kIterationLimit:
      return HighsStatus::kWarning;
    case HighsModelStatus::kSolutionLimit:
      return HighsStatus::kWarning;
    case HighsModelStatus::kInterrupt:
      return HighsStatus::kWarning;
    case HighsModelStatus::kHighsInterrupt:
      return HighsStatus::kWarning;
    case HighsModelStatus::kUnknown:
      return HighsStatus::kWarning;
    default:
      return HighsStatus::kError;
  }
}

std::string findModelObjectiveName(const HighsLp* lp,
                                   const HighsHessian* hessian) {
  // Return any non-trivial current objective name
  if (lp->objective_name_ != "") return lp->objective_name_;

  std::string objective_name = "";
  // Determine whether there is a nonzero cost vector
  bool has_objective = false;
  for (HighsInt iCol = 0; iCol < lp->num_col_; iCol++) {
    if (lp->col_cost_[iCol]) {
      has_objective = true;
      break;
    }
  }
  if (!has_objective && hessian) {
    // Zero cost vector, so only chance of an objective comes from any
    // Hessian
    has_objective = (hessian->dim_ != 0);
  }
  HighsInt pass = 0;
  for (;;) {
    // Loop until a valid name is found. Vanishingly unlikely to have
    // to pass more than once, since check for objective name
    // duplicating a row name is very unlikely to fail
    //
    // So set up an appropriate name (stem)
    if (has_objective) {
      objective_name = "Obj";
    } else {
      objective_name = "NoObj";
    }
    // If there are no row names, then the objective name is certainly
    // OK
    if (lp->row_names_.size() == 0) break;
    if (pass != 0) objective_name += std::to_string(pass);
    // Ensure that the objective name doesn't clash with any row names
    bool ok_objective_name = true;
    for (HighsInt iRow = 0; iRow < lp->num_row_; iRow++) {
      std::string trimmed_name = lp->row_names_[iRow];
      trimmed_name = trim(trimmed_name);
      if (objective_name == trimmed_name) {
        ok_objective_name = false;
        break;
      }
    }
    if (ok_objective_name) break;
    pass++;
  }
  assert(objective_name != "");
  return objective_name;
}

/*
void print_map(std::string comment, const std::map<std::string, HighsInt>& m)
{
    std::cout << comment;

  for (const auto& n : m)
      std::cout << n.first << " = " << n.second << "; ";
    std::cout << '\n';
}
*/

/*
bool repeatedNames(const std::vector<std::string> name) {
  const HighsInt num_name = name.size();
  // With no names, cannot have any repeated
  if (num_name == 0) return false;
  std::map<std::string, HighsInt> name_map;
  for (HighsInt ix = 0; ix < num_name; ix++) {
    auto search = name_map.find(name[ix]);
    if (search != name_map.end()) return true;
    //    printf("Search for %s yields %d\n", name[ix].c_str(),
    //    int(search->second));
    name_map.insert({name[ix], ix});
    //    print_map("Map\n", name_map);
  }
  return false;
}
*/
