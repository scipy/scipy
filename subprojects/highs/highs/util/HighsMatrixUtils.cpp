/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "util/HighsMatrixUtils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "util/HighsHash.h"

HighsStatus assessMatrix(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt vec_dim,
    const HighsInt num_vec, vector<HighsInt>& matrix_start,
    vector<HighsInt>& matrix_index, vector<double>& matrix_value,
    const double small_matrix_value, const double large_matrix_value,
    const bool sum_duplicates, const std::string* col_names,
    const std::string* row_names) {
  vector<HighsInt> matrix_p_end;
  const bool partitioned = false;
  return assessMatrix(log_options, class_name, matrix_name, vec_dim, num_vec,
                      partitioned, matrix_start, matrix_p_end, matrix_index,
                      matrix_value, small_matrix_value, large_matrix_value,
                      sum_duplicates, col_names, row_names);
}

HighsStatus assessMatrix(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt vec_dim,
    const HighsInt num_vec, const bool partitioned,
    vector<HighsInt>& matrix_start, vector<HighsInt>& matrix_p_end,
    vector<HighsInt>& matrix_index, vector<double>& matrix_value,
    const double small_matrix_value, const double large_matrix_value,
    const bool sum_duplicates, const std::string* col_names,
    const std::string* row_names) {
  // Mechanism for summing duplicates and stripping out small values
  // isn't coded for partitioned row-wise matrix, but it is generated
  // internally from a constraint matrix, so shouldn't have duplicates
  // (and they aren't summed in constraint matrices) and any small
  // values it inherits (from presolve) are respected
  if (partitioned) {
    // Only checking structure for a partitioned (row-wise) matrix, as
    // its values are inherited and duplicates are never summed
    assert(small_matrix_value == 0);
    assert(large_matrix_value == kHighsInf);
    assert(sum_duplicates == false);
  }
  if (!legalMatrixDimensions(log_options, class_name, matrix_name, num_vec,
                             partitioned, matrix_start, matrix_p_end,
                             matrix_index, matrix_value))
    return HighsStatus::kError;

  if ((vec_dim == 0 || num_vec == 0) && matrix_start.size() > 1 &&
      matrix_index.size() > 0 && matrix_value.size() > 0) {
    // Most likely that vec_dim and/or num_vec have not been set so,
    // rather than trimming matrix_start_, matrix_.index_ and
    // matrix_.value_ to correspond to an empty matrix, log an error
    // and return
    highsLogUser(log_options, HighsLogType::kError,
                 "%s %s number of vectors (%d) and vector dimension (%d) "
                 "with (start, index, value) sizes of (%d, %d, %d)\n",
                 class_name.c_str(), matrix_name.c_str(), int(vec_dim),
                 int(num_vec), int(matrix_start.size()),
                 int(matrix_index.size()), int(matrix_value.size()));
    highsLogUser(log_options, HighsLogType::kError,
                 "   Suggests that the number of vectors and/or vector "
                 "dimension have not been set\n");
    return HighsStatus::kError;
  }

  // Lambda for quoting column name (if possible) when logging anomalies
  auto possible_col_name = [&](const HighsInt ix) {
    std::string name = "";
    if (col_names) name = " (col " + col_names[ix] + ")";
    return name;
  };

  // Lambda for quoting row name (if possible) when logging anomalies
  auto possible_row_name = [&](const HighsInt component) {
    std::string name = "";
    if (row_names) name = " (row " + row_names[component] + ")";
    return name;
  };

  bool error_found = false;
  bool warning_found = false;

  const HighsInt num_nz = matrix_start[num_vec];
  // Assess the starts
  //
  // Check whether the first start is zero
  if (matrix_start[0]) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s %s matrix start vector begins with %" HIGHSINT_FORMAT
                 " rather than 0\n",
                 class_name.c_str(), matrix_name.c_str(), matrix_start[0]);
    return HighsStatus::kError;
  }
  // Set up previous_start for a fictitious previous empty packed vector
  HighsInt previous_start = matrix_start[0];
  // Set up this_start to be the first start in case num_vec = 0
  HighsInt this_start = matrix_start[0];
  HighsInt this_p_end = 0;
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    this_start = matrix_start[ix];
    bool this_start_too_small = this_start < previous_start;
    if (this_start_too_small) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s %s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal start of %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT
                   " = "
                   "previous start\n",
                   class_name.c_str(), matrix_name.c_str(), ix, this_start,
                   previous_start);
      return HighsStatus::kError;
    }
    if (partitioned) {
      this_p_end = matrix_p_end[ix];
      bool this_p_end_too_small = this_p_end < this_start;
      if (this_p_end_too_small) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s %s matrix packed vector %" HIGHSINT_FORMAT
                     " has illegal partition end of %" HIGHSINT_FORMAT
                     " < %" HIGHSINT_FORMAT
                     " = "
                     " start\n",
                     class_name.c_str(), matrix_name.c_str(), ix, this_p_end,
                     this_start);
        return HighsStatus::kError;
      }
    }
    previous_start = this_start;
  }
  bool this_start_too_big = this_start > num_nz;
  if (this_start_too_big) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "%s %s matrix packed vector %" HIGHSINT_FORMAT
        " has illegal start of %" HIGHSINT_FORMAT " > %" HIGHSINT_FORMAT
        " = "
        "number of nonzeros\n",
        class_name.c_str(), matrix_name.c_str(), num_vec, this_start, num_nz);
    return HighsStatus::kError;
  }
  if (partitioned) {
    bool this_p_end_too_big = this_p_end > num_nz;
    if (this_p_end_too_big) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s %s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal partition end of %" HIGHSINT_FORMAT
                   " > %" HIGHSINT_FORMAT
                   " = "
                   "number of nonzeros\n",
                   class_name.c_str(), matrix_name.c_str(), num_vec, this_p_end,
                   num_nz);
      return HighsStatus::kError;
    }
  }
  // Assess the indices and values
  // Count the number of acceptable indices/values
  HighsInt new_num_nz = 0;
  HighsInt num_small_value = 0;
  double max_small_value = 0;
  double min_small_value = kHighsInf;
  HighsInt num_large_value = 0;
  double max_large_value = 0;
  double min_large_value = kHighsInf;
  HighsInt num_duplicate = 0;
  const HighsInt illegal_el = -1;
  // When duplicates only have to be identified (and not summed) use
  // highs_hash
  const bool use_highs_hash = !sum_duplicates;
  HighsHashTable<HighsInt> highs_hash;
  // Otherwise, use index-key map to identify entry number of
  // duplicates so they can be summed
  std::unordered_map<HighsInt, HighsInt> index_el_map;

  // Lambda to handle matrix entry magnitude
  auto okMatrixValueMagnitude = [&](const HighsInt el) {
    double abs_value = fabs(matrix_value[el]);
    // Check that the value is not too large
    if (abs_value >= large_matrix_value) {
      if (max_large_value < abs_value) max_large_value = abs_value;
      if (min_large_value > abs_value) min_large_value = abs_value;
      num_large_value++;
      return false;
    }
    if (abs_value <= small_matrix_value) {
      if (max_small_value < abs_value) max_small_value = abs_value;
      if (min_small_value > abs_value) min_small_value = abs_value;
      num_small_value++;
      return false;
    }
    return true;
  };

  // Lambda to identify whether an index has already occurred: when
  // summing duplicates, also return the position of its entry
  auto findIndex = [&](const HighsInt component, HighsInt& previous_el) {
    if (use_highs_hash) return highs_hash.find(component) != nullptr;
    auto found_component = index_el_map.find(component);
    if (found_component == index_el_map.end()) return false;
    previous_el = found_component->second;
    return true;
  };

  // Lambda to record the occurrence of an index: when summing
  // duplicates, also record the position of its entry so that any
  // duplicates can be summed into it
  auto recordIndex = [&](const HighsInt component, const HighsInt el) {
    if (use_highs_hash) {
      highs_hash.insert(component);
    } else {
      index_el_map.insert({component, el});
    }
  };

  for (HighsInt ix = 0; ix < num_vec; ix++) {
    HighsInt from_el = matrix_start[ix];
    HighsInt to_el = matrix_start[ix + 1];
    // Account for any index-value pairs removed so far
    matrix_start[ix] = new_num_nz;
    for (HighsInt el = from_el; el < to_el; el++) {
      // Check the index
      HighsInt component = matrix_index[el];
      // Check that the index is non-negative
      bool legal_component = component >= 0;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s %s matrix packed vector %" HIGHSINT_FORMAT
                     "%s, entry %" HIGHSINT_FORMAT
                     ", is illegal index %" HIGHSINT_FORMAT "\n",
                     class_name.c_str(), matrix_name.c_str(), ix,
                     possible_col_name(ix).c_str(), el, component);
        return HighsStatus::kError;
      }
      // Check that the index does not exceed the vector dimension
      legal_component = component < vec_dim;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s %s matrix packed vector %" HIGHSINT_FORMAT
                     "%s, entry %" HIGHSINT_FORMAT
                     ", is illegal index "
                     "%12" HIGHSINT_FORMAT " >= %" HIGHSINT_FORMAT
                     " = vector dimension\n",
                     class_name.c_str(), matrix_name.c_str(), ix,
                     possible_col_name(ix).c_str(), el, component, vec_dim);
        return HighsStatus::kError;
      }
      // If not summing duplicates, can check matrix value magnitude.
      // This is done before checking for duplicates, so that an entry
      // that is ignored because its value is small is not identified as
      // a duplicate, and is not recorded, whatever the order of the
      // entries.
      if (!sum_duplicates && !okMatrixValueMagnitude(el)) continue;
      // Check whether the index has already occurred.
      HighsInt previous_el = illegal_el;
      if (findIndex(component, previous_el)) {
        if (sum_duplicates) {
          num_duplicate++;
          // Sum the duplicate entry, making sure that it's been
          // assigned
          assert(previous_el != illegal_el);
          matrix_value[previous_el] += matrix_value[el];
          continue;
        }
        highsLogUser(log_options, HighsLogType::kError,
                     "%s %s matrix packed vector %" HIGHSINT_FORMAT
                     "%s, has duplicate entries for index %" HIGHSINT_FORMAT
                     "%s\n",
                     class_name.c_str(), matrix_name.c_str(), ix,
                     possible_col_name(ix).c_str(), component,
                     possible_row_name(component).c_str());
        return HighsStatus::kError;
      }
      // Not a duplicate, so record the index and shift the index and
      // value of the OK entry to the new position in the index and value
      // vectors, and increment the new number of nonzeros
      recordIndex(component, new_num_nz);
      matrix_index[new_num_nz] = matrix_index[el];
      matrix_value[new_num_nz] = matrix_value[el];
      new_num_nz++;
    }
    if (sum_duplicates) {
      // Now that any duplicates have been summed, can check matrix
      // value magnitude
      from_el = matrix_start[ix];
      to_el = new_num_nz;
      // Reset new_num_nz
      new_num_nz = matrix_start[ix];
      for (HighsInt el = from_el; el < to_el; el++) {
        if (!okMatrixValueMagnitude(el)) continue;
        // Shift the index and value of the OK entry to the new
        // position in the index and value vectors, and increment
        // the new number of nonzeros
        matrix_index[new_num_nz] = matrix_index[el];
        matrix_value[new_num_nz] = matrix_value[el];
        new_num_nz++;
      }
    }

    if (use_highs_hash) {
      highs_hash.clear();
    } else {
      index_el_map.clear();
    }
  }  // Loop 0; num_vec
  if (num_duplicate) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "%s %s matrix contains %" HIGHSINT_FORMAT
                 " duplicate entr%s: summed\n",
                 class_name.c_str(), matrix_name.c_str(), num_duplicate,
                 highsIntToPlural(num_duplicate, true).c_str());
  }
  if (num_large_value) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s %s matrix contains %" HIGHSINT_FORMAT
                 " |value| in [%.2g, %.2g] at/above tolerance of %g\n",
                 class_name.c_str(), matrix_name.c_str(), num_large_value,
                 min_large_value, max_large_value, large_matrix_value);
    error_found = true;
  }
  if (max_small_value > 0) {
    // If the small values are explicit zeros, then no model
    // information is lost if they are ignored, so only report and
    // return a warning if small nonzeros are ignored
    assert(num_small_value > 0);
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s %s matrix contains %" HIGHSINT_FORMAT
                 " |value| in [%.2g, %.2g] "
                 "at/below tolerance of %g: ignored\n",
                 class_name.c_str(), matrix_name.c_str(), num_small_value,
                 min_small_value, max_small_value, small_matrix_value);
    warning_found = true;
  }
  if (partitioned && num_nz > new_num_nz) {
    // No entries must be lost from partitioned matrices
    const HighsInt loss = num_nz - new_num_nz;
    assert(loss > 0);
    highsLogUser(log_options, HighsLogType::kError,
                 "Partitioned matrix with %d entries has lost %d entr%s\n",
                 int(num_nz), int(loss), highsIntToPlural(loss, true).c_str());
    error_found = true;
  }

  if (num_nz < new_num_nz) {
    // No entries must be gained!
    const HighsInt gain = new_num_nz - num_nz;
    assert(gain > 0);
    highsLogUser(log_options, HighsLogType::kError,
                 "Matrix with %d entries has gained %d entr%s\n", int(num_nz),
                 int(gain), highsIntToPlural(gain, true).c_str());
    error_found = true;
  }
  matrix_start[num_vec] = new_num_nz;

  // The start, index and value vectors may be larger than necessary,
  // even if no entries have been removed from the matrix, so trim
  // them to prevent a bug in presolve.
  HighsInt trimmed = 0;

  auto logTrimming = [&](const std::string& name) {
    assert(trimmed > 0);
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Trimmed %d excess entr%s from %s %s data member %s\n",
                 int(trimmed), highsIntToPlural(trimmed, true).c_str(),
                 matrix_name.c_str(), class_name.c_str(), name.c_str());
    warning_found = true;
  };

  // For start, warn about any trimming
  trimmed = matrix_start.size() - num_vec - 1;
  if (trimmed) {
    matrix_start.resize(num_vec + 1);
    logTrimming("start_");
  }
  // For index and value, warn about trimming to num_nz, but not
  // trimming to new_num_nz, since this is logged otherwise
  trimmed = matrix_index.size() - num_nz;
  if (trimmed) logTrimming("index_");
  matrix_index.resize(new_num_nz);

  trimmed = matrix_value.size() - num_nz;
  if (trimmed) logTrimming("value_");
  matrix_value.resize(new_num_nz);

  HighsStatus return_status = HighsStatus::kOk;
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  return return_status;
}

bool legalMatrixDimensions(const HighsLogOptions& log_options,
                           const std::string& class_name,
                           const std::string& matrix_name,
                           const HighsInt num_vec, const bool partitioned,
                           const vector<HighsInt>& start,
                           const vector<HighsInt>& p_end,
                           const vector<HighsInt>& index,
                           const vector<double>& value) {
  bool ok = true;
  // Assess main dimensions
  const bool legal_num_vec = num_vec >= 0;
  if (!legal_num_vec)
    highsLogUser(log_options, HighsLogType::kError,
                 "Matrix dimension validation for %s %s fails on number of "
                 "vectors = %d < 0\n",
                 class_name.c_str(), matrix_name.c_str(), (int)num_vec);
  ok = legal_num_vec && ok;
  const bool legal_start_size = (HighsInt)start.size() >= num_vec + 1;
  if (!legal_start_size) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Matrix dimension validation for %s %s fails on start size = %d < %d = "
        "num vectors + 1\n",
        class_name.c_str(), matrix_name.c_str(), (int)start.size(),
        (int)(num_vec + 1));
  }
  ok = legal_start_size && ok;
  if (partitioned) {
    const bool legal_p_end_size = (HighsInt)p_end.size() >= num_vec;
    if (!legal_p_end_size)
      highsLogUser(
          log_options, HighsLogType::kError,
          "Matrix dimension validation for %s %s fails on p_end size = %d < %d "
          "= num vectors\n",
          class_name.c_str(), matrix_name.c_str(), (int)p_end.size(),
          (int)(num_vec + 1));
    ok = (HighsInt)p_end.size() >= num_vec && ok;
  }
  // Possibly check the sizes of the index and value vectors. Can only
  // do this with the number of nonzeros, and this is only known if
  // the start vector has a legal size. Setting num_nz = 0 otherwise
  // means that all tests pass, as they just check that the sizes of
  // the index and value vectors are non-negative.
  const HighsInt num_nz = legal_start_size ? start[num_vec] : 0;
  if (num_nz >= 0) {
    const bool legal_index_size = (HighsInt)index.size() >= num_nz;
    if (!legal_index_size)
      highsLogUser(
          log_options, HighsLogType::kError,
          "Matrix dimension validation for %s %s fails on index size = %d < %d "
          "= number of nonzeros\n",
          class_name.c_str(), matrix_name.c_str(), (int)index.size(),
          (int)num_nz);
    ok = legal_index_size && ok;
    const bool legal_value_size = (HighsInt)value.size() >= num_nz;
    if (!legal_value_size)
      highsLogUser(
          log_options, HighsLogType::kError,
          "Matrix dimension validation for %s %s fails on value size = %d < %d "
          "= number of nonzeros\n",
          class_name.c_str(), matrix_name.c_str(), (int)value.size(),
          (int)num_nz);
    ok = legal_value_size && ok;
  } else {
    highsLogUser(log_options, HighsLogType::kError,
                 "Matrix dimension validation for %s %s fails on number of "
                 "nonzeros = %d < 0\n",
                 class_name.c_str(), matrix_name.c_str(), (int)num_nz);
    ok = false;
  }
  return ok;
}

HighsStatus trimMatrixDimensions(
    const HighsLogOptions& log_options, const std::string& class_name,
    const std::string& matrix_name, const HighsInt num_vec,
    const bool partitioned, vector<HighsInt>& start, vector<HighsInt>& p_end,
    vector<HighsInt>& index, vector<double>& value) {
  assert(num_vec >= 0);
  assert(start.size() >= static_cast<size_t>(num_vec + 1));
  HighsInt trimmed = 0;
  HighsStatus status = HighsStatus::kOk;

  auto logTrimming = [&](const std::string& name) {
    assert(trimmed > 0);
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Trimmed %d excess entr%s from %s %s data member %s\n",
                 int(trimmed), highsIntToPlural(trimmed, true).c_str(),
                 matrix_name.c_str(), class_name.c_str(), name.c_str());
    status = HighsStatus::kWarning;
  };

  trimmed = start.size() - num_vec - 1;
  if (trimmed) {
    start.resize(num_vec + 1);
    logTrimming("start_");
  }

  const HighsInt num_nz = start[num_vec];
  assert(num_nz >= 0);

  if (partitioned) {
    trimmed = p_end.size() - num_vec;
    if (trimmed) {
      p_end.resize(num_vec);
      logTrimming("p_end_");
    }
  }

  trimmed = index.size() - num_nz;
  if (trimmed) {
    index.resize(num_nz);
    logTrimming("index_");
  }

  trimmed = value.size() - num_nz;
  if (trimmed) {
    value.resize(num_nz);
    logTrimming("value_");
  }

  return status;
}
