/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "util/HighsUtils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "HConfig.h"
//#include "io/HighsIO.h"
//#include "lp_data/HConst.h"

/*
int getOmpNumThreads() {
  return omp_get_num_threads()
}
*/

bool assessIndexCollection(const HighsOptions& options,
                           const HighsIndexCollection& index_collection) {
  // Check parameter for each technique of defining an index collection
  if (index_collection.is_interval_) {
    // Changing by interval: check the parameters and that check set and mask
    // are false
    if (index_collection.is_set_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both interval and set");
      return false;
    }
    if (index_collection.is_mask_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both interval and mask");
      return false;
    }
    if (index_collection.from_ < 0) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index interval lower limit is %d < 0",
                      index_collection.from_);
      return false;
    }
    if (index_collection.to_ > index_collection.dimension_ - 1) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index interval upper limit is %d > %d",
                      index_collection.to_, index_collection.dimension_ - 1);
      return false;
    }
  } else if (index_collection.is_set_) {
    // Changing by set: check the parameters and check that interval and mask
    // are false
    if (index_collection.is_interval_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both set and interval");
      return false;
    }
    if (index_collection.is_mask_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both set and mask");
      return false;
    }
    if (index_collection.set_ == NULL) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index set is NULL");
      return false;
    }
    // Check that the values in the vector of integers are ascending
    const int* set = index_collection.set_;
    const int set_entry_upper = index_collection.dimension_ - 1;
    int prev_set_entry = -1;
    for (int k = 0; k < index_collection.set_num_entries_; k++) {
      if (set[k] < 0 || set[k] > set_entry_upper) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Index set entry set[%d] = %d is out of bounds [0, %d]",
                        k, set[k], set_entry_upper);
        return false;
      }
      if (set[k] <= prev_set_entry) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Index set entry set[%d] = %d is not greater than "
                        "previous entry %d",
                        k, set[k], prev_set_entry);
        return false;
      }
      prev_set_entry = set[k];
    }
  } else if (index_collection.is_mask_) {
    // Changing by mask: check the parameters and check that set and interval
    // are false
    if (index_collection.mask_ == NULL) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index mask is NULL");
      return false;
    }
    if (index_collection.is_interval_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both mask and interval");
      return false;
    }
    if (index_collection.is_set_) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Index collection is both mask and set");
      return false;
    }
  } else {
    // No method defined
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Undefined index collection");
    return false;
  }
  return true;
}

bool limitsForIndexCollection(const HighsOptions& options,
                              const HighsIndexCollection& index_collection,
                              int& from_k, int& to_k) {
  if (index_collection.is_interval_) {
    from_k = index_collection.from_;
    to_k = index_collection.to_;
  } else if (index_collection.is_set_) {
    from_k = 0;
    to_k = index_collection.set_num_entries_ - 1;
  } else if (index_collection.is_mask_) {
    from_k = 0;
    to_k = index_collection.dimension_ - 1;
  } else {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Undefined index collection");
    return false;
  }
  return true;
}

void updateIndexCollectionOutInIndex(
    const HighsIndexCollection& index_collection, int& out_from_ix,
    int& out_to_ix, int& in_from_ix, int& in_to_ix, int& current_set_entry) {
  if (index_collection.is_interval_) {
    out_from_ix = index_collection.from_;
    out_to_ix = index_collection.to_;
    in_from_ix = index_collection.to_ + 1;
    in_to_ix = index_collection.dimension_ - 1;
  } else if (index_collection.is_set_) {
    out_from_ix = index_collection.set_[current_set_entry];
    out_to_ix = out_from_ix;  //+1;
    current_set_entry++;
    int current_set_entry0 = current_set_entry;
    for (int set_entry = current_set_entry0;
         set_entry < index_collection.set_num_entries_; set_entry++) {
      int ix = index_collection.set_[set_entry];
      if (ix > out_to_ix + 1) break;
      out_to_ix = index_collection.set_[current_set_entry];
      current_set_entry++;
    }
    in_from_ix = out_to_ix + 1;
    if (current_set_entry < index_collection.set_num_entries_) {
      in_to_ix = index_collection.set_[current_set_entry] - 1;
    } else {
      // Account for getting to the end of the set
      in_to_ix = index_collection.dimension_ - 1;
    }
  } else {
    out_from_ix = in_to_ix + 1;
    out_to_ix = index_collection.dimension_ - 1;
    for (int ix = in_to_ix + 1; ix < index_collection.dimension_; ix++) {
      if (!index_collection.mask_[ix]) {
        out_to_ix = ix - 1;
        break;
      }
    }
    in_from_ix = out_to_ix + 1;
    in_to_ix = index_collection.dimension_ - 1;
    for (int ix = out_to_ix + 1; ix < index_collection.dimension_; ix++) {
      if (index_collection.mask_[ix]) {
        in_to_ix = ix - 1;
        break;
      }
    }
  }
}

int dataSizeOfIndexCollection(const HighsIndexCollection& index_collection) {
  if (index_collection.is_set_) {
    return index_collection.set_num_entries_;
  } else {
    if (index_collection.is_interval_) {
      return index_collection.to_ - index_collection.from_ + 1;
    } else {
      return index_collection.dimension_;
    }
  }
}

bool intUserDataNotNull(FILE* logfile, const int* user_data,
                        const std::string name) {
  bool null_data = false;
  if (user_data == NULL) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "User-supplied %s are NULL", name.c_str());
    null_data = true;
  }
  assert(!null_data);
  return null_data;
}

bool doubleUserDataNotNull(FILE* logfile, const double* user_data,
                           const std::string name) {
  bool null_data = false;
  if (user_data == NULL) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "User-supplied %s are NULL", name.c_str());
    null_data = true;
  }
  assert(!null_data);
  return null_data;
}

double getNorm2(const std::vector<double> values) {
  double sum = 0;
  int values_size = values.size();
  for (int i = 0; i < values_size; i++) sum += values[i] * values[i];
  return sum;
}

bool highs_isInfinity(double val) {
  if (val >= HIGHS_CONST_INF) return true;
  return false;
}

double highsRelativeDifference(const double v0, const double v1) {
  return fabs(v0 - v1) / std::max(v0, std::max(v1, 1.0));
}

#ifdef HiGHSDEV
void analyseVectorValues(const char* message, int vecDim,
                         const std::vector<double>& vec, bool analyseValueList,
                         std::string model_name) {
  if (vecDim == 0) return;
  double log10 = log(10.0);
  const int nVK = 20;
  int nNz = 0;
  int nPosInfV = 0;
  int nNegInfV = 0;
  std::vector<int> posVK;
  std::vector<int> negVK;
  posVK.resize(nVK + 1, 0);
  negVK.resize(nVK + 1, 0);

  const int VLsMxZ = 10;
  std::vector<int> VLsK;
  std::vector<double> VLsV;
  VLsK.resize(VLsMxZ, 0);
  VLsV.resize(VLsMxZ, 0);
  // Ensure that 1.0 and -1.0 are counted
  const int PlusOneIx = 0;
  const int MinusOneIx = 1;
  bool excessVLsV = false;
  int VLsZ = 2;
  VLsV[PlusOneIx] = 1.0;
  VLsV[MinusOneIx] = -1.0;

  for (int ix = 0; ix < vecDim; ix++) {
    double v = vec[ix];
    double absV = std::fabs(v);
    int log10V;
    if (absV > 0) {
      // Nonzero value
      nNz++;
      if (highs_isInfinity(-v)) {
        //-Inf value
        nNegInfV++;
      } else if (highs_isInfinity(v)) {
        //+Inf value
        nPosInfV++;
      } else {
        // Finite nonzero value
        if (absV == 1) {
          log10V = 0;
        } else if (absV == 10) {
          log10V = 1;
        } else if (absV == 100) {
          log10V = 2;
        } else if (absV == 1000) {
          log10V = 3;
        } else {
          log10V = log(absV) / log10;
        }
        if (log10V >= 0) {
          int k = std::min(log10V, nVK);
          posVK[k]++;
        } else {
          int k = std::min(-log10V, nVK);
          negVK[k]++;
        }
      }
    }
    if (analyseValueList) {
      if (v == 1.0) {
        VLsK[PlusOneIx]++;
      } else if (v == -1.0) {
        VLsK[MinusOneIx]++;
      } else {
        int fdIx = -1;
        for (int ix = 2; ix < VLsZ; ix++) {
          if (v == VLsV[ix]) {
            fdIx = ix;
            break;
          }
        }
        if (fdIx == -1) {
          // New value
          if (VLsZ < VLsMxZ) {
            fdIx = VLsZ;
            VLsV[fdIx] = v;
            VLsK[fdIx]++;
            VLsZ++;
          } else {
            excessVLsV = true;
          }
        } else {
          // Existing value
          VLsK[fdIx]++;
        }
      }
    }
  }
  printf("%s of dimension %d with %d nonzeros (%3d%%): Analysis\n", message,
         vecDim, nNz, 100 * nNz / vecDim);
  if (nNegInfV > 0) printf("%12d values are -Inf\n", nNegInfV);
  if (nPosInfV > 0) printf("%12d values are +Inf\n", nPosInfV);
  int k = nVK;
  int vK = posVK[k];
  if (vK > 0) printf("%12d values satisfy 10^(%3d) <= v < Inf\n", vK, k);
  for (int k = nVK - 1; k >= 0; k--) {
    int vK = posVK[k];
    if (vK > 0)
      printf("%12d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, k, k + 1);
  }
  for (int k = 1; k <= nVK; k++) {
    int vK = negVK[k];
    if (vK > 0)
      printf("%12d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, -k, 1 - k);
  }
  vK = vecDim - nNz;
  if (vK > 0) printf("%12d values are zero\n", vK);
  if (analyseValueList) {
    printf("           Value distribution:");
    if (excessVLsV) printf(" More than %d different values", VLsZ);
    printf("\n            Value        Count\n");
    for (int ix = 0; ix < VLsZ; ix++) {
      int pct = ((100.0 * VLsK[ix]) / vecDim) + 0.5;
      printf("     %12g %12d (%3d%%)\n", VLsV[ix], VLsK[ix], pct);
    }
    printf("grep_value_distrib,%s,%d", model_name.c_str(), VLsZ);
    printf(",");
    if (excessVLsV) printf("!");
    for (int ix = 0; ix < VLsZ; ix++) printf(",%g", VLsV[ix]);
    printf("\n");
  }
}

void analyseMatrixSparsity(const char* message, int numCol, int numRow,
                           const std::vector<int>& Astart,
                           const std::vector<int>& Aindex) {
  if (numCol == 0) return;
  std::vector<int> rowCount;
  std::vector<int> colCount;

  rowCount.assign(numRow, 0);
  colCount.resize(numCol);

  for (int col = 0; col < numCol; col++) {
    colCount[col] = Astart[col + 1] - Astart[col];
    for (int el = Astart[col]; el < Astart[col + 1]; el++)
      rowCount[Aindex[el]]++;
  }
  const int maxCat = 10;
  std::vector<int> CatV;
  std::vector<int> rowCatK;
  std::vector<int> colCatK;
  CatV.resize(maxCat + 1);
  rowCatK.assign(maxCat + 1, 0);
  colCatK.assign(maxCat + 1, 0);

  CatV[1] = 1;
  for (int cat = 2; cat < maxCat + 1; cat++) {
    CatV[cat] = 2 * CatV[cat - 1];
  }

  int maxRowCount = 0;
  int maxColCount = 0;
  for (int col = 0; col < numCol; col++) {
    maxColCount = std::max(colCount[col], maxColCount);
    int fdCat = maxCat;
    for (int cat = 0; cat < maxCat - 1; cat++) {
      if (colCount[col] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    colCatK[fdCat]++;
  }

  for (int row = 0; row < numRow; row++) {
    maxRowCount = std::max(rowCount[row], maxRowCount);
    int fdCat = maxCat;
    for (int cat = 0; cat < maxCat - 1; cat++) {
      if (rowCount[row] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    rowCatK[fdCat]++;
  }

  printf("\n%s\n\n", message);
  int lastRpCat;
  for (int cat = 0; cat < maxCat + 1; cat++) {
    if (colCatK[cat]) lastRpCat = cat;
  }
  int cat = maxCat;
  if (colCatK[cat]) lastRpCat = cat;
  int sumK = 0;
  int pct;
  double v;
  int sumPct = 0;
  for (int cat = 0; cat < lastRpCat; cat++) {
    sumK += colCatK[cat];
    v = 100 * colCatK[cat];
    v = v / numCol + 0.5;
    pct = v;
    sumPct += pct;
    printf("%12d (%3d%%) columns of count in [%3d, %3d]\n", colCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  sumK += colCatK[cat];
  v = 100 * colCatK[cat];
  v = v / numCol + 0.5;
  pct = v;
  sumPct += pct;
  if (cat == maxCat) {
    printf("%12d (%3d%%) columns of count in [%3d, inf]\n", colCatK[cat], pct,
           CatV[cat]);
  } else {
    printf("%12d (%3d%%) columns of count in [%3d, %3d]\n", colCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }
  printf("Max count is %d / %d\n\n", maxColCount, numRow);

  lastRpCat = -1;
  for (int cat = 0; cat < maxCat + 1; cat++) {
    if (rowCatK[cat]) lastRpCat = cat;
  }
  cat = maxCat;
  if (rowCatK[cat]) lastRpCat = cat;
  sumK = 0;
  pct = 0;
  v = 0;
  sumPct = 0;
  for (int cat = 0; cat < lastRpCat; cat++) {
    sumK += rowCatK[cat];
    v = 100 * rowCatK[cat];
    v = v / numRow + 0.5;
    pct = v;
    sumPct += pct;
    printf("%12d (%3d%%)    rows of count in [%3d, %3d]\n", rowCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  sumK += rowCatK[cat];
  v = 100 * rowCatK[cat];
  v = v / numRow + 0.5;
  pct = v;
  sumPct += pct;
  if (cat == maxCat) {
    printf("%12d (%3d%%)    rows of count in [%3d, inf]\n", rowCatK[cat], pct,
           CatV[cat]);
  } else {
    printf("%12d (%3d%%)    rows of count in [%3d, %3d]\n", rowCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }
  printf("Max count is %d / %d\n", maxRowCount, numCol);
}

bool initialiseValueDistribution(const std::string distribution_name,
                                 const std::string value_name,
                                 const double min_value_limit,
                                 const double max_value_limit,
                                 const double base_value_limit,
                                 HighsValueDistribution& value_distribution) {
  assert(min_value_limit > 0);
  assert(max_value_limit > 0);
  assert(base_value_limit > 1);
  value_distribution.distribution_name_ = distribution_name;
  value_distribution.value_name_ = value_name;
  if (min_value_limit <= 0) return false;
  if (max_value_limit < min_value_limit) return false;
  int num_count;
  if (min_value_limit == max_value_limit) {
    // For counting values below and above a value
    num_count = 1;
  } else {
    if (base_value_limit <= 0) return false;
    const double log_ratio = log(max_value_limit / min_value_limit);
    const double log_base_value_limit = log(base_value_limit);
    //    printf("initialiseValueDistribution: log_ratio = %g;
    //    log_base_value_limit = %g; log_ratio/log_base_value_limit = %g\n",
    //	   log_ratio, log_base_value_limit, log_ratio/log_base_value_limit);
    num_count = log_ratio / log_base_value_limit + 1;
  }
  //  printf("initialiseValueDistribution: num_count = %d\n", num_count);
  value_distribution.count_.assign(num_count + 1, 0);
  value_distribution.limit_.assign(num_count, 0);
  value_distribution.limit_[0] = min_value_limit;
  //  printf("Interval  0 is [%10.4g, %10.4g)\n", 0.0,
  //  value_distribution.limit_[0]);
  for (int i = 1; i < num_count; i++) {
    value_distribution.limit_[i] =
        base_value_limit * value_distribution.limit_[i - 1];
    //    printf("Interval %2d is [%10.4g, %10.4g)\n", i,
    //    value_distribution.limit_[i-1], value_distribution.limit_[i]);
  }
  //  printf("Interval %2d is [%10.4g, inf)\n", num_count,
  //  value_distribution.limit_[num_count-1]);
  value_distribution.num_count_ = num_count;
  value_distribution.num_zero_ = 0;
  value_distribution.num_one_ = 0;
  value_distribution.min_value_ = HIGHS_CONST_INF;
  value_distribution.max_value_ = 0;
  value_distribution.sum_count_ = 0;
  return true;
}

bool updateValueDistribution(const double value,
                             HighsValueDistribution& value_distribution) {
  if (value_distribution.num_count_ < 0) return false;
  value_distribution.sum_count_++;
  const double abs_value = fabs(value);
  value_distribution.min_value_ =
      std::min(abs_value, value_distribution.min_value_);
  value_distribution.max_value_ =
      std::max(abs_value, value_distribution.max_value_);
  if (!abs_value) {
    value_distribution.num_zero_++;
    return true;
  }
  if (abs_value == 1.0) {
    value_distribution.num_one_++;
    return true;
  }
  for (int i = 0; i < value_distribution.num_count_; i++) {
    if (abs_value < value_distribution.limit_[i]) {
      value_distribution.count_[i]++;
      return true;
    }
  }
  value_distribution.count_[value_distribution.num_count_]++;
  return true;
}

double doublePercentage(const int of, const int in) {
  return ((100.0 * of) / in);
}

int integerPercentage(const int of, const int in) {
  const double double_percentage = ((100.0 * of) / in) + 0.4999;
  return (int)double_percentage;
}

bool printValueDistribution(const HighsValueDistribution& value_distribution,
                            const int mu) {
  if (value_distribution.sum_count_ <= 0) return false;
  const int num_count = value_distribution.num_count_;
  if (num_count < 0) return false;
  if (value_distribution.distribution_name_ != "")
    printf("\n%s\n", value_distribution.distribution_name_.c_str());
  std::string value_name = value_distribution.value_name_;
  bool not_reported_ones = true;
  int sum_count = value_distribution.num_zero_ + value_distribution.num_one_;
  double sum_percentage = 0;
  const double min_value = value_distribution.min_value_;
  for (int i = 0; i < num_count + 1; i++)
    sum_count += value_distribution.count_[i];
  if (!sum_count) return false;
  printf("Min value = %g\n", min_value);
  printf("     Minimum %svalue is %10.4g", value_name.c_str(), min_value);
  if (mu > 0) {
    printf("  corresponding to  %10d / %10d\n", (int)(min_value * mu), mu);
  } else {
    printf("\n");
  }
  printf("     Maximum %svalue is %10.4g", value_name.c_str(),
         value_distribution.max_value_);
  if (mu > 0) {
    printf("  corresponding to  %10d / %10d\n",
           (int)(value_distribution.max_value_ * mu), mu);
  } else {
    printf("\n");
  }
  int sum_report_count = 0;
  double percentage;
  int int_percentage;
  int count = value_distribution.num_zero_;
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) are %10.4g\n", count, value_name.c_str(),
           int_percentage, 0.0);
    sum_report_count += count;
  }
  count = value_distribution.count_[0];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) in (%10.4g, %10.4g)", count,
           value_name.c_str(), int_percentage, 0.0,
           value_distribution.limit_[0]);
    sum_report_count += count;
    if (mu > 0) {
      printf(" corresponding to (%10d, %10d)\n", 0,
             (int)(value_distribution.limit_[0] * mu));
    } else {
      printf("\n");
    }
  }
  for (int i = 1; i < num_count; i++) {
    if (not_reported_ones && value_distribution.limit_[i - 1] >= 1.0) {
      count = value_distribution.num_one_;
      if (count) {
        percentage = doublePercentage(count, sum_count);
        sum_percentage += percentage;
        int_percentage = percentage;
        printf("%12d %svalues (%3d%%) are             %10.4g", count,
               value_name.c_str(), int_percentage, 1.0);
        sum_report_count += count;
        if (mu > 0) {
          printf(" corresponding to %10d\n", mu);
        } else {
          printf("\n");
        }
      }
      not_reported_ones = false;
    }
    count = value_distribution.count_[i];
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) in [%10.4g, %10.4g)", count,
             value_name.c_str(), int_percentage,
             value_distribution.limit_[i - 1], value_distribution.limit_[i]);
      sum_report_count += count;
      if (mu > 0) {
        printf(" corresponding to [%10d, %10d)\n",
               (int)(value_distribution.limit_[i - 1] * mu),
               (int)(value_distribution.limit_[i] * mu));
      } else {
        printf("\n");
      }
    }
  }
  if (not_reported_ones && value_distribution.limit_[num_count - 1] >= 1.0) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) are             %10.4g", count,
             value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        printf("  corresponding to  %10d\n", mu);
      } else {
        printf("\n");
      }
    }
    not_reported_ones = false;
  }
  count = value_distribution.count_[num_count];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) in [%10.4g,        inf)", count,
           value_name.c_str(), int_percentage,
           value_distribution.limit_[num_count - 1]);
    sum_report_count += count;
    if (mu > 0) {
      printf(" corresponding to [%10d,        inf)\n",
             (int)(value_distribution.limit_[num_count - 1] * mu));
    } else {
      printf("\n");
    }
  }
  if (not_reported_ones) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) are             %10.4g", count,
             value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        printf("  corresponding to  %10d\n", mu);
      } else {
        printf("\n");
      }
    }
  }
  printf("%12d %svalues\n", sum_count, value_name.c_str());
  if (sum_report_count != sum_count)
    printf("ERROR: %d = sum_report_count != sum_count = %d\n", sum_report_count,
           sum_count);
  return true;
}
#endif

bool initialiseScatterData(const int max_num_point,
                           HighsScatterData& scatter_data) {
  if (max_num_point < 1) return false;
  scatter_data.max_num_point_ = max_num_point;
  scatter_data.num_point_ = 0;
  scatter_data.last_point_ = -1;
  scatter_data.value0_.resize(max_num_point);
  scatter_data.value1_.resize(max_num_point);
  scatter_data.have_regression_coeff_ = false;
  scatter_data.num_error_comparison_ = 0;
  scatter_data.num_awful_linear_ = 0;
  scatter_data.num_awful_log_ = 0;
  scatter_data.num_bad_linear_ = 0;
  scatter_data.num_bad_log_ = 0;
  scatter_data.num_fair_linear_ = 0;
  scatter_data.num_fair_log_ = 0;
  scatter_data.num_better_linear_ = 0;
  scatter_data.num_better_log_ = 0;
  return true;
}

bool updateScatterData(const double value0, const double value1,
                       HighsScatterData& scatter_data) {
  if (value0 <= 0 || value0 <= 0) return false;
  scatter_data.num_point_++;
  scatter_data.last_point_++;
  if (scatter_data.last_point_ == scatter_data.max_num_point_)
    scatter_data.last_point_ = 0;
  scatter_data.value0_[scatter_data.last_point_] = value0;
  scatter_data.value1_[scatter_data.last_point_] = value1;
  return true;
}

bool regressScatterData(HighsScatterData& scatter_data) {
  if (scatter_data.num_point_ < 5) return true;
  double log_x;
  double log_y;
  double sum_log_x = 0;
  double sum_log_y = 0;
  double sum_log_xlog_x = 0;
  double sum_log_xlog_y = 0;
  double x;
  double y;
  double sum_x = 0;
  double sum_y = 0;
  double sum_xx = 0;
  double sum_xy = 0;
  int point_num = 0;
  for (int pass = 0; pass < 2; pass++) {
    int from_point;
    int to_point;
    if (pass == 0) {
      from_point = scatter_data.last_point_;
      to_point = std::min(scatter_data.num_point_, scatter_data.max_num_point_);
    } else {
      from_point = 0;
      to_point = scatter_data.last_point_;
    }
    for (int point = from_point; point < to_point; point++) {
      point_num++;
      x = scatter_data.value0_[point];
      y = scatter_data.value1_[point];
      sum_x += x;
      sum_y += y;
      sum_xx += x * x;
      sum_xy += x * y;
      log_x = log(x);
      log_y = log(y);
      sum_log_x += log_x;
      sum_log_y += log_y;
      sum_log_xlog_x += log_x * log_x;
      sum_log_xlog_y += log_x * log_y;
    }
  }
  double double_num = 1.0 * point_num;
  // Linear regression
  double det = double_num * sum_xx - sum_x * sum_x;
  if (fabs(det) < 1e-8) return true;
  scatter_data.linear_coeff0_ = (sum_xx * sum_y - sum_x * sum_xy) / det;
  scatter_data.linear_coeff1_ = (-sum_x * sum_y + double_num * sum_xy) / det;
  // Log regression
  det = double_num * sum_log_xlog_x - sum_log_x * sum_log_x;
  if (fabs(det) < 1e-8) return true;
  scatter_data.log_coeff0_ =
      (sum_log_xlog_x * sum_log_y - sum_log_x * sum_log_xlog_y) / det;
  scatter_data.log_coeff0_ = exp(scatter_data.log_coeff0_);
  scatter_data.log_coeff1_ =
      (-sum_log_x * sum_log_y + double_num * sum_log_xlog_y) / det;
  // Look at the errors in the two approaches
  scatter_data.have_regression_coeff_ = true;
  if (scatter_data.num_point_ < scatter_data.max_num_point_) return true;

  scatter_data.num_error_comparison_++;
  computeScatterDataRegressionError(scatter_data);
  const double linear_error = scatter_data.linear_regression_error_;
  const double log_error = scatter_data.log_regression_error_;

  const bool report_awful_error = false;
  if (linear_error > awful_regression_error ||
      log_error > awful_regression_error) {
    if (linear_error > awful_regression_error) {
      scatter_data.num_awful_linear_++;
      if (report_awful_error)
        printf("Awful linear regression error = %g\n", linear_error);
    }
    if (log_error > awful_regression_error) {
      scatter_data.num_awful_log_++;
      if (report_awful_error)
        printf("Awful log regression error = %g\n", log_error);
    }
    if (report_awful_error)
      computeScatterDataRegressionError(scatter_data, true);
  }
  if (linear_error > bad_regression_error) scatter_data.num_bad_linear_++;
  if (log_error > bad_regression_error) scatter_data.num_bad_log_++;
  if (linear_error > fair_regression_error) scatter_data.num_fair_linear_++;
  if (log_error > fair_regression_error) scatter_data.num_fair_log_++;
  if (linear_error < log_error) {
    scatter_data.num_better_linear_++;
  } else if (linear_error > log_error) {
    scatter_data.num_better_log_++;
  }
  //  printf("Linear regression error = %g\n", linear_error);
  //  printf("Log    regression error = %g\n", log_error);
  return true;
}

bool predictFromScatterData(const HighsScatterData& scatter_data,
                            const double value0, double& predicted_value1,
                            const bool log_regression) {
  if (!scatter_data.have_regression_coeff_) return false;
  if (log_regression) {
    predicted_value1 =
        scatter_data.log_coeff0_ * pow(value0, scatter_data.log_coeff1_);
    return true;
  } else {
    predicted_value1 =
        scatter_data.linear_coeff0_ + scatter_data.linear_coeff1_ * value0;
    return true;
  }
}

bool computeScatterDataRegressionError(HighsScatterData& scatter_data,
                                       const bool print) {
  if (!scatter_data.have_regression_coeff_) return false;
  if (scatter_data.num_point_ < scatter_data.max_num_point_) return false;
  double sum_log_error = 0;
  if (print)
    printf(
        "Log regression\nPoint     Value0     Value1 PredValue1      Error\n");
  for (int point = 0; point < scatter_data.max_num_point_; point++) {
    double value0 = scatter_data.value0_[point];
    double value1 = scatter_data.value1_[point];
    double predicted_value1;
    if (predictFromScatterData(scatter_data, value0, predicted_value1, true)) {
      double error = fabs(predicted_value1 - value1);  // / fabs(value1);
      if (
          //	10*error > awful_regression_error &&
          print)
        printf("%5d %10.4g %10.4g %10.4g %10.4g\n", point, value0, value1,
               predicted_value1, error);
      sum_log_error += error;
    }
  }
  if (print)
    printf("                                       %10.4g\n", sum_log_error);
  double sum_linear_error = 0;
  if (print)
    printf(
        "Linear regression\nPoint     Value0     Value1 PredValue1      "
        "Error\n");
  for (int point = 0; point < scatter_data.max_num_point_; point++) {
    double value0 = scatter_data.value0_[point];
    double value1 = scatter_data.value1_[point];
    double predicted_value1;
    if (predictFromScatterData(scatter_data, value0, predicted_value1)) {
      double error = fabs(predicted_value1 - value1);  //  / fabs(value1);
      if (
          //	10*error > awful_regression_error &&
          print)
        printf("%5d %10.4g %10.4g %10.4g %10.4g\n", point, value0, value1,
               predicted_value1, error);
      sum_linear_error += error;
    }
  }
  if (print)
    printf("                                       %10.4g\n", sum_linear_error);
  scatter_data.log_regression_error_ = sum_log_error;
  scatter_data.linear_regression_error_ = sum_linear_error;
  return true;
}

bool printScatterData(std::string name, const HighsScatterData& scatter_data) {
  if (!scatter_data.num_point_) return true;
  double x;
  double y;
  int point_num = 0;
  printf("%s scatter data\n", name.c_str());
  const int to_point =
      std::min(scatter_data.num_point_, scatter_data.max_num_point_);
  for (int point = scatter_data.last_point_ + 1; point < to_point; point++) {
    point_num++;
    x = scatter_data.value0_[point];
    y = scatter_data.value1_[point];
    printf("%d,%10.4g,%10.4g,%d\n", point, x, y, point_num);
  }
  for (int point = 0; point <= scatter_data.last_point_; point++) {
    point_num++;
    x = scatter_data.value0_[point];
    y = scatter_data.value1_[point];
    printf("%d,%10.4g,%10.4g,%d\n", point, x, y, point_num);
  }
  printf("Linear regression coefficients,%10.4g,%10.4g\n",
         scatter_data.linear_coeff0_, scatter_data.linear_coeff1_);
  printf("Log    regression coefficients,%10.4g,%10.4g\n",
         scatter_data.log_coeff0_, scatter_data.log_coeff1_);
  return true;
}

void printScatterDataRegressionComparison(
    std::string name, const HighsScatterData& scatter_data) {
  if (!scatter_data.num_error_comparison_) return;
  printf("\n%s scatter data regression\n", name.c_str());
  printf("%10d regression error comparisons\n",
         scatter_data.num_error_comparison_);
  printf("%10d regression awful  linear (>%10.4g)\n",
         scatter_data.num_awful_linear_, awful_regression_error);
  printf("%10d regression awful  log    (>%10.4g)\n",
         scatter_data.num_awful_log_, awful_regression_error);
  printf("%10d regression bad    linear (>%10.4g)\n",
         scatter_data.num_bad_linear_, bad_regression_error);
  printf("%10d regression bad    log    (>%10.4g)\n", scatter_data.num_bad_log_,
         bad_regression_error);
  printf("%10d regression fair   linear (>%10.4g)\n",
         scatter_data.num_fair_linear_, fair_regression_error);
  printf("%10d regression fair   log    (>%10.4g)\n",
         scatter_data.num_fair_log_, fair_regression_error);
  printf("%10d regression better linear\n", scatter_data.num_better_linear_);
  printf("%10d regression better log\n", scatter_data.num_better_log_);
}
