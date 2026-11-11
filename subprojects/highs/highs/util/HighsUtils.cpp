/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "util/HighsUtils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <vector>

#include "util/HighsSort.h"

HighsInt create(HighsIndexCollection& index_collection, const HighsInt from_col,
                const HighsInt to_col, const HighsInt dimension) {
  if (from_col < 0) return kIndexCollectionCreateIllegalInterval;
  if (to_col >= dimension) return kIndexCollectionCreateIllegalInterval;
  index_collection.dimension_ = dimension;
  index_collection.is_interval_ = true;
  index_collection.from_ = from_col;
  index_collection.to_ = to_col;
  return kIndexCollectionCreateOk;
}

HighsInt create(HighsIndexCollection& index_collection,
                const HighsInt num_set_entries, const HighsInt* set,
                const HighsInt dimension) {
  // Create an index collection for the given set - so long as it is strictly
  // ordered
  if (num_set_entries < 0) return kIndexCollectionCreateIllegalSetSize;
  if (dimension < 0) return kIndexCollectionCreateIllegalSetDimension;
  index_collection.dimension_ = dimension;
  index_collection.is_set_ = true;
  index_collection.set_ = {set, set + num_set_entries};
  index_collection.set_num_entries_ = num_set_entries;
  if (!increasingSetOk(index_collection.set_, 1, 0, true))
    return kIndexCollectionCreateIllegalSetOrder;
  for (HighsInt ix = 0; ix < num_set_entries; ix++)
    if (set[ix] < 0 || set[ix] >= dimension) return -(ix + 1);
  return kIndexCollectionCreateOk;
}

HighsInt create(HighsIndexCollection& index_collection, const HighsInt* mask,
                const HighsInt dimension) {
  // Create an index collection for the given mask
  if (dimension < 0) return kIndexCollectionCreateIllegalMaskSize;
  index_collection.dimension_ = dimension;
  index_collection.is_mask_ = true;
  index_collection.mask_ = {mask, mask + dimension};
  return kIndexCollectionCreateOk;
}

void highsSparseTranspose(HighsInt numRow, HighsInt numCol,
                          const std::vector<HighsInt>& Astart,
                          const std::vector<HighsInt>& Aindex,
                          const std::vector<double>& Avalue,
                          std::vector<HighsInt>& ARstart,
                          std::vector<HighsInt>& ARindex,
                          std::vector<double>& ARvalue) {
  // Make a AR copy
  std::vector<HighsInt> iwork(numRow, 0);
  ARstart.resize(numRow + 1, 0);
  HighsInt AcountX = Astart[numCol];
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (HighsInt k = 0; k < AcountX; k++) {
    assert(Aindex[k] < numRow);
    iwork[Aindex[k]]++;
  }
  for (HighsInt i = 1; i <= numRow; i++)
    ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (HighsInt i = 0; i < numRow; i++) iwork[i] = ARstart[i];
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      HighsInt iRow = Aindex[k];
      HighsInt iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }
}

bool ok(const HighsIndexCollection& index_collection) {
  // Check parameter for each technique of defining an index collection
  if (index_collection.is_interval_) {
    // Changing by interval: check the parameters and that check set and mask
    // are false
    if (index_collection.is_set_) {
      printf("Index collection is both interval and set\n");
      return false;
    }
    if (index_collection.is_mask_) {
      printf("Index collection is both interval and mask\n");
      return false;
    }
    if (index_collection.from_ < 0) {
      printf("Index interval lower limit is %" HIGHSINT_FORMAT " < 0\n",
             index_collection.from_);
      return false;
    }
    if (index_collection.to_ > index_collection.dimension_ - 1) {
      printf("Index interval upper limit is %" HIGHSINT_FORMAT
             " > %" HIGHSINT_FORMAT "\n",
             index_collection.to_, index_collection.dimension_ - 1);
      return false;
    }
  } else if (index_collection.is_set_) {
    // Changing by set: check the parameters and check that interval and mask
    // are false
    if (index_collection.is_interval_) {
      printf("Index collection is both set and interval\n");
      return false;
    }
    if (index_collection.is_mask_) {
      printf("Index collection is both set and mask\n");
      return false;
    }
    if (index_collection.set_.empty()) {
      printf("Index set is NULL\n");
      return false;
    }
    // Check that the values in the vector of integers are ascending
    const vector<HighsInt>& set = index_collection.set_;
    const HighsInt num_entries = index_collection.set_num_entries_;
    const HighsInt entry_upper = index_collection.dimension_ - 1;
    HighsInt prev_set_entry = -1;
    for (HighsInt k = 0; k < num_entries; k++) {
      if (set[k] < 0 || set[k] > entry_upper) {
        printf("Index set entry set[%" HIGHSINT_FORMAT "] = %" HIGHSINT_FORMAT
               " is out of bounds [0, %" HIGHSINT_FORMAT "]\n",
               k, set[k], entry_upper);
        return false;
      }
      if (set[k] <= prev_set_entry) {
        printf("Index set entry set[%" HIGHSINT_FORMAT "] = %" HIGHSINT_FORMAT
               " is not greater than "
               "previous entry %" HIGHSINT_FORMAT "\n",
               k, set[k], prev_set_entry);
        return false;
      }
      prev_set_entry = set[k];
    }
    // This was the old check done independently, and should be
    // equivalent.
    assert(increasingSetOk(set, 0, entry_upper, true));
  } else if (index_collection.is_mask_) {
    // Changing by mask: check the parameters and check that set and interval
    // are false
    if (index_collection.mask_.empty()) {
      printf("Index mask is NULL\n");
      return false;
    }
    if (index_collection.is_interval_) {
      printf("Index collection is both mask and interval\n");
      return false;
    }
    if (index_collection.is_set_) {
      printf("Index collection is both mask and set\n");
      return false;
    }
  } else {
    // No method defined
    printf("Undefined index collection\n");
    return false;
  }
  return true;
}

void limits(const HighsIndexCollection& index_collection, HighsInt& from_k,
            HighsInt& to_k) {
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
    assert(1 == 0);
  }
}

void updateOutInIndex(const HighsIndexCollection& index_collection,
                      HighsInt& out_from_ix, HighsInt& out_to_ix,
                      HighsInt& in_from_ix, HighsInt& in_to_ix,
                      HighsInt& current_set_entry) {
  if (index_collection.is_interval_) {
    out_from_ix = index_collection.from_;
    out_to_ix = index_collection.to_;
    in_from_ix = index_collection.to_ + 1;
    in_to_ix = index_collection.dimension_ - 1;
  } else if (index_collection.is_set_) {
    out_from_ix = index_collection.set_[current_set_entry];
    out_to_ix = out_from_ix;
    current_set_entry++;
    HighsInt current_set_entry0 = current_set_entry;
    for (HighsInt set_entry = current_set_entry0;
         set_entry < index_collection.set_num_entries_; set_entry++) {
      HighsInt ix = index_collection.set_[set_entry];
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
    for (HighsInt ix = in_to_ix + 1; ix < index_collection.dimension_; ix++) {
      if (!index_collection.mask_[ix]) {
        out_to_ix = ix - 1;
        break;
      }
    }
    in_from_ix = out_to_ix + 1;
    in_to_ix = index_collection.dimension_ - 1;
    for (HighsInt ix = out_to_ix + 1; ix < index_collection.dimension_; ix++) {
      if (index_collection.mask_[ix]) {
        in_to_ix = ix - 1;
        break;
      }
    }
  }
}

HighsInt dataSize(const HighsIndexCollection& index_collection) {
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

bool highsVarTypeUserDataNotNull(const HighsLogOptions& log_options,
                                 const HighsVarType* user_data,
                                 const std::string& name) {
  bool null_data = false;
  if (user_data == NULL) {
    highsLogUser(log_options, HighsLogType::kError,
                 "User-supplied %s are NULL\n", name.c_str());
    null_data = true;
  }
  assert(!null_data);
  return null_data;
}

bool intUserDataNotNull(const HighsLogOptions& log_options,
                        const HighsInt* user_data, const std::string& name) {
  bool null_data = false;
  if (user_data == NULL) {
    highsLogUser(log_options, HighsLogType::kError,
                 "User-supplied %s are NULL\n", name.c_str());
    null_data = true;
  }
  assert(!null_data);
  return null_data;
}

bool doubleUserDataNotNull(const HighsLogOptions& log_options,
                           const double* user_data, const std::string& name) {
  bool null_data = false;
  if (user_data == NULL) {
    highsLogUser(log_options, HighsLogType::kError,
                 "User-supplied %s are NULL\n", name.c_str());
    null_data = true;
  }
  assert(!null_data);
  return null_data;
}

double getNorm2(const std::vector<double>& values) {
  double sum = 0;
  HighsInt values_size = values.size();
  for (HighsInt i = 0; i < values_size; i++) sum += values[i] * values[i];
  return sum;
}

bool highs_isInfinity(double val) {
  if (val >= kHighsInf) return true;
  return false;
}

double highsRelativeDifference(const double v0, const double v1) {
  return fabs(v0 - v1) / std::max(v0, std::max(v1, 1.0));
}

void analyseVectorValues(const HighsLogOptions* log_options,
                         const std::string& message, HighsInt vecDim,
                         const std::vector<double>& vec, bool analyseValueList,
                         const std::string& model_name) {
  assert(vecDim == int(vec.size()));
  if (vecDim == 0) return;
  double log10 = log(10.0);
  const HighsInt nVK = 20;
  HighsInt nNz = 0;
  HighsInt nPosInfV = 0;
  HighsInt nNegInfV = 0;
  std::vector<HighsInt> posVK;
  std::vector<HighsInt> negVK;
  posVK.resize(nVK + 1, 0);
  negVK.resize(nVK + 1, 0);

  const HighsInt VLsMxZ = 10;
  std::vector<HighsInt> VLsK;
  std::vector<double> VLsV;
  VLsK.resize(VLsMxZ, 0);
  VLsV.resize(VLsMxZ, 0);
  // Ensure that 1.0 and -1.0 are counted
  const HighsInt PlusOneIx = 0;
  const HighsInt MinusOneIx = 1;
  bool excessVLsV = false;
  HighsInt VLsZ = 2;
  VLsV[PlusOneIx] = 1.0;
  VLsV[MinusOneIx] = -1.0;
  double min_abs_value = kHighsInf;
  double max_abs_value = 0;
  for (HighsInt ix = 0; ix < vecDim; ix++) {
    double v = vec[ix];
    double absV = std::fabs(v);
    if (absV) {
      min_abs_value = std::min(absV, min_abs_value);
      max_abs_value = std::max(absV, max_abs_value);
    }
    HighsInt log10V;
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
          HighsInt k = std::min(log10V, nVK);
          posVK[k]++;
        } else {
          HighsInt k = std::min(-log10V, nVK);
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
        HighsInt fdIx = -1;
        for (HighsInt ix = 2; ix < VLsZ; ix++) {
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
  }  // for (HighsInt ix = 0; ix < vecDim; ix++)
  // If there are no nonzeros, min_abs_value retains its starting
  // value of inf
  if (!nNz) min_abs_value = 0;
  std::stringstream ss;
  ss.str(std::string());
  ss << message << " of dimension " << vecDim;
  if (nNz < vecDim)
    ss << message << " with " << nNz << " nonzeros ("
       << int(1e2 * double(nNz) / double(vecDim)) << "%%)";
  highsReportDevInfo(
      log_options,
      highsFormatToString("%s: values in [%.2g, %.2g]\n", ss.str().c_str(),
                          min_abs_value, max_abs_value));
  if (nNegInfV > 0)
    highsReportDevInfo(
        log_options, highsFormatToString(
                         "%12" HIGHSINT_FORMAT " values are -Inf\n", nNegInfV));
  if (nPosInfV > 0)
    highsReportDevInfo(
        log_options, highsFormatToString(
                         "%12" HIGHSINT_FORMAT " values are +Inf\n", nPosInfV));
  HighsInt k = nVK;
  HighsInt vK = posVK[k];
  if (vK > 0)
    highsReportDevInfo(log_options, highsFormatToString(
                                        "%12" HIGHSINT_FORMAT
                                        " values satisfy 10^(%3" HIGHSINT_FORMAT
                                        ") <= v < Inf\n",
                                        vK, k));
  for (HighsInt k = nVK - 1; k >= 0; k--) {
    HighsInt vK = posVK[k];
    if (vK > 0)
      highsReportDevInfo(
          log_options,
          highsFormatToString("%12" HIGHSINT_FORMAT
                              " values satisfy 10^(%3" HIGHSINT_FORMAT
                              ") <= v < 10^(%3" HIGHSINT_FORMAT ")\n",
                              vK, k - 1, k));
  }
  for (HighsInt k = 1; k <= nVK; k++) {
    HighsInt vK = negVK[k];
    if (vK > 0)
      highsReportDevInfo(
          log_options,
          highsFormatToString("%12" HIGHSINT_FORMAT
                              " values satisfy 10^(%3" HIGHSINT_FORMAT
                              ") <= v < 10^(%3" HIGHSINT_FORMAT ")\n",
                              vK, -k - 1, -k));
  }
  vK = vecDim - nNz;
  if (vK > 0)
    highsReportDevInfo(
        log_options,
        highsFormatToString("%12" HIGHSINT_FORMAT " values are zero\n", vK));
  if (analyseValueList) {
    highsReportDevInfo(log_options,
                       highsFormatToString("           Value distribution:"));
    if (excessVLsV)
      highsReportDevInfo(
          log_options,
          highsFormatToString(
              " More than %" HIGHSINT_FORMAT " different values", VLsZ));
    highsReportDevInfo(
        log_options, highsFormatToString("\n            Value        Count\n"));
    for (HighsInt ix = 0; ix < VLsZ; ix++) {
      if (!VLsK[ix]) continue;
      HighsInt pct = ((100.0 * VLsK[ix]) / vecDim) + 0.5;
      highsReportDevInfo(log_options,
                         highsFormatToString("     %12g %12" HIGHSINT_FORMAT
                                             " (%3" HIGHSINT_FORMAT "%%)\n",
                                             VLsV[ix], VLsK[ix], pct));
    }
    highsReportDevInfo(
        log_options,
        highsFormatToString("grep_value_distrib,%s,%" HIGHSINT_FORMAT "",
                            model_name.c_str(), VLsZ));
    highsReportDevInfo(log_options, highsFormatToString(","));
    if (excessVLsV) highsReportDevInfo(log_options, highsFormatToString("!"));
    for (HighsInt ix = 0; ix < VLsZ; ix++)
      highsReportDevInfo(log_options, highsFormatToString(",%g", VLsV[ix]));
    highsReportDevInfo(log_options, highsFormatToString("\n"));
  }
}

void analyseVectorValues(const HighsLogOptions* log_options,
                         const std::string& message, HighsInt vecDim,
                         const std::vector<HighsInt>& vec,
                         const std::string& model_name) {
  if (vecDim == 0) return;
  const HighsInt VLsMxZ = 10;
  std::vector<std::pair<HighsInt, HighsInt>> VLs;
  bool excessVLsV = false;
  HighsInt VLsZ = 0;
  HighsInt min_value = kHighsIInf;
  HighsInt max_value = 0;
  HighsInt nNz = 0;

  for (HighsInt ix = 0; ix < vecDim; ix++) {
    HighsInt v = vec[ix];
    min_value = std::min(v, min_value);
    max_value = std::max(v, max_value);
    if (v != 0) nNz++;

    HighsInt fdIx = -1;
    for (HighsInt ix = 0; ix < VLsZ; ix++) {
      if (v == VLs[ix].first) {
        fdIx = ix;
        break;
      }
    }
    if (fdIx == -1) {
      // New value
      if (VLsZ < VLsMxZ) {
        fdIx = VLsZ;
        VLs.push_back(std::make_pair(v, 1));
        VLsZ++;
      } else {
        excessVLsV = true;
      }
    } else {
      // Existing value
      VLs[fdIx].second++;
    }
  }
  assert(VLsZ == int(VLs.size()));
  std::sort(VLs.begin(), VLs.end());
  highsReportDevInfo(
      log_options,
      highsFormatToString(
          "%s of dimension %d with %d nonzeros (%3d%%) in [%d, %d]\n",
          message.c_str(), int(vecDim), int(nNz), int(100 * nNz / vecDim),
          int(min_value), int(max_value)));
  highsReportDevInfo(log_options,
                     highsFormatToString("           Value distribution:"));
  if (excessVLsV)
    highsReportDevInfo(
        log_options,
        highsFormatToString(" More than %" HIGHSINT_FORMAT " different values",
                            VLsZ));
  highsReportDevInfo(log_options,
                     highsFormatToString("\n            Value        Count\n"));
  for (HighsInt ix = 0; ix < VLsZ; ix++) {
    if (!VLs[ix].second) continue;
    HighsInt pct = ((100.0 * VLs[ix].second) / vecDim) + 0.5;
    highsReportDevInfo(log_options,
                       highsFormatToString("     %12d %12d (%3d%%)\n",
                                           VLs[ix].first, VLs[ix].second, pct));
  }
  /*
  highsReportDevInfo(log_options,
                     highsFormatToString("grep_value_distrib,%s,%"
  HIGHSINT_FORMAT "", model_name.c_str(), VLsZ));
  highsReportDevInfo(log_options, highsFormatToString(","));
  if (excessVLsV) highsReportDevInfo(log_options, highsFormatToString("!"));
  for (HighsInt ix = 0; ix < VLsZ; ix++)
    highsReportDevInfo(log_options, highsFormatToString(",%g", VLs[ix].first));
  highsReportDevInfo(log_options, highsFormatToString("\n"));
  */
}

void analyseMatrixSparsity(const HighsLogOptions& log_options,
                           const char* message, HighsInt numCol,
                           HighsInt numRow, const std::vector<HighsInt>& Astart,
                           const std::vector<HighsInt>& Aindex) {
  if (numCol == 0) return;
  std::vector<HighsInt> rowCount;
  std::vector<HighsInt> colCount;

  rowCount.assign(numRow, 0);
  colCount.resize(numCol);

  for (HighsInt col = 0; col < numCol; col++) {
    colCount[col] = Astart[col + 1] - Astart[col];
    for (HighsInt el = Astart[col]; el < Astart[col + 1]; el++)
      rowCount[Aindex[el]]++;
  }
  analyseVectorValues(&log_options, "Column counts", numCol, colCount);
  analyseVectorValues(&log_options, "Row counts", numRow, rowCount);

  const HighsInt maxCat = 10;
  std::vector<HighsInt> CatV;
  std::vector<HighsInt> rowCatK;
  std::vector<HighsInt> colCatK;
  CatV.resize(maxCat + 1);
  rowCatK.assign(maxCat + 1, 0);
  colCatK.assign(maxCat + 1, 0);

  CatV[1] = 1;
  for (HighsInt cat = 2; cat < maxCat + 1; cat++) {
    CatV[cat] = 2 * CatV[cat - 1];
  }

  HighsInt maxRowCount = 0;
  HighsInt maxColCount = 0;
  for (HighsInt col = 0; col < numCol; col++) {
    maxColCount = std::max(colCount[col], maxColCount);
    HighsInt fdCat = maxCat;
    for (HighsInt cat = 0; cat < maxCat - 1; cat++) {
      if (colCount[col] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    colCatK[fdCat]++;
  }

  for (HighsInt row = 0; row < numRow; row++) {
    maxRowCount = std::max(rowCount[row], maxRowCount);
    HighsInt fdCat = maxCat;
    for (HighsInt cat = 0; cat < maxCat - 1; cat++) {
      if (rowCount[row] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    rowCatK[fdCat]++;
  }

  highsLogDev(log_options, HighsLogType::kInfo, "\n%s\n\n", message);
  HighsInt lastRpCat = -1;
  for (HighsInt cat = 0; cat < maxCat + 1; cat++) {
    if (colCatK[cat]) lastRpCat = cat;
  }
  HighsInt cat = maxCat;
  if (colCatK[cat]) lastRpCat = cat;
  HighsInt pct;
  double v;
  for (HighsInt cat = 0; cat < lastRpCat; cat++) {
    v = 100 * colCatK[cat];
    v = v / numCol + 0.5;
    pct = v;
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%) columns of count in [%3" HIGHSINT_FORMAT
                ", %3" HIGHSINT_FORMAT "]\n",
                colCatK[cat], pct, CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  v = 100 * colCatK[cat];
  v = v / numCol + 0.5;
  pct = v;
  if (cat == maxCat) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%) columns of count in [%3" HIGHSINT_FORMAT ", inf]\n",
                colCatK[cat], pct, CatV[cat]);
  } else {
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%) columns of count in [%3" HIGHSINT_FORMAT
                ", %3" HIGHSINT_FORMAT "]\n",
                colCatK[cat], pct, CatV[cat], CatV[cat + 1] - 1);
  }
  highsLogDev(log_options, HighsLogType::kInfo,
              "Max count is %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT "\n\n",
              maxColCount, numRow);

  lastRpCat = -1;
  for (HighsInt cat = 0; cat < maxCat + 1; cat++) {
    if (rowCatK[cat]) lastRpCat = cat;
  }
  cat = maxCat;
  if (rowCatK[cat]) lastRpCat = cat;
  pct = 0;
  v = 0;
  for (HighsInt cat = 0; cat < lastRpCat; cat++) {
    v = 100 * rowCatK[cat];
    v = v / numRow + 0.5;
    pct = v;
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%)    rows of count in [%3" HIGHSINT_FORMAT
                ", %3" HIGHSINT_FORMAT "]\n",
                rowCatK[cat], pct, CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  v = 100 * rowCatK[cat];
  v = v / numRow + 0.5;
  pct = v;
  if (cat == maxCat) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%)    rows of count in [%3" HIGHSINT_FORMAT ", inf]\n",
                rowCatK[cat], pct, CatV[cat]);
  } else {
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
                "%%)    rows of count in [%3" HIGHSINT_FORMAT
                ", %3" HIGHSINT_FORMAT "]\n",
                rowCatK[cat], pct, CatV[cat], CatV[cat + 1] - 1);
  }
  highsLogDev(log_options, HighsLogType::kInfo,
              "Max count is %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT "\n",
              maxRowCount, numCol);
}

bool initialiseValueDistribution(const std::string& distribution_name,
                                 const std::string& value_name,
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
  HighsInt num_count;
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
  //  printf("initialiseValueDistribution: num_count = %" HIGHSINT_FORMAT "\n",
  //  num_count);
  value_distribution.count_.assign(num_count + 1, 0);
  value_distribution.limit_.assign(num_count, 0);
  value_distribution.limit_[0] = min_value_limit;
  //  printf("Interval  0 is [%10.4g, %10.4g)\n", 0.0,
  //  value_distribution.limit_[0]);
  for (HighsInt i = 1; i < num_count; i++) {
    value_distribution.limit_[i] =
        base_value_limit * value_distribution.limit_[i - 1];
    //    printf("Interval %2" HIGHSINT_FORMAT " is [%10.4g, %10.4g)\n", i,
    //    value_distribution.limit_[i-1], value_distribution.limit_[i]);
  }
  //  printf("Interval %2" HIGHSINT_FORMAT " is [%10.4g, inf)\n", num_count,
  //  value_distribution.limit_[num_count-1]);
  value_distribution.num_count_ = num_count;
  value_distribution.num_zero_ = 0;
  value_distribution.num_one_ = 0;
  value_distribution.min_value_ = kHighsInf;
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
  for (HighsInt i = 0; i < value_distribution.num_count_; i++) {
    if (abs_value < value_distribution.limit_[i]) {
      value_distribution.count_[i]++;
      return true;
    }
  }
  value_distribution.count_[value_distribution.num_count_]++;
  return true;
}

double doublePercentage(const HighsInt of, const HighsInt in) {
  return ((100.0 * of) / in);
}

HighsInt integerPercentage(const HighsInt of, const HighsInt in) {
  const double double_percentage = ((100.0 * of) / in) + 0.4999;
  return (HighsInt)double_percentage;
}

bool logValueDistribution(const HighsLogOptions& log_options,
                          const HighsValueDistribution& value_distribution,
                          const HighsInt mu) {
  if (value_distribution.sum_count_ <= 0) return false;
  const HighsInt num_count = value_distribution.num_count_;
  if (num_count < 0) return false;
  if (value_distribution.distribution_name_ != "")
    highsLogDev(log_options, HighsLogType::kInfo, "\n%s\n",
                value_distribution.distribution_name_.c_str());
  std::string value_name = value_distribution.value_name_;
  bool not_reported_ones = true;
  HighsInt sum_count =
      value_distribution.num_zero_ + value_distribution.num_one_;
  const double min_value = value_distribution.min_value_;
  for (HighsInt i = 0; i < num_count + 1; i++)
    sum_count += value_distribution.count_[i];
  if (!sum_count) return false;
  highsLogDev(log_options, HighsLogType::kInfo, "Min value = %g\n", min_value);
  highsLogDev(log_options, HighsLogType::kInfo,
              "     Minimum %svalue is %10.4g", value_name.c_str(), min_value);
  if (mu > 0) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "  corresponding to  %10" HIGHSINT_FORMAT
                " / %10" HIGHSINT_FORMAT "\n",
                (HighsInt)(min_value * mu), mu);
  } else {
    highsLogDev(log_options, HighsLogType::kInfo, "\n");
  }
  highsLogDev(log_options, HighsLogType::kInfo,
              "     Maximum %svalue is %10.4g", value_name.c_str(),
              value_distribution.max_value_);
  if (mu > 0) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "  corresponding to  %10" HIGHSINT_FORMAT
                " / %10" HIGHSINT_FORMAT "\n",
                (HighsInt)(value_distribution.max_value_ * mu), mu);
  } else {
    highsLogDev(log_options, HighsLogType::kInfo, "\n");
  }
  HighsInt sum_report_count = 0;
  double percentage;
  HighsInt int_percentage;
  HighsInt count = value_distribution.num_zero_;
  if (count) {
    percentage = doublePercentage(count, sum_count);
    int_percentage = percentage;
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                "%%) are %10.4g\n",
                count, value_name.c_str(), int_percentage, 0.0);
    sum_report_count += count;
  }
  count = value_distribution.count_[0];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    int_percentage = percentage;
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                "%%) in (%10.4g, %10.4g)",
                count, value_name.c_str(), int_percentage, 0.0,
                value_distribution.limit_[0]);
    sum_report_count += count;
    if (mu > 0) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  " corresponding to (%10" HIGHSINT_FORMAT
                  ", %10" HIGHSINT_FORMAT ")\n",
                  0, (HighsInt)(value_distribution.limit_[0] * mu));
    } else {
      highsLogDev(log_options, HighsLogType::kInfo, "\n");
    }
  }
  for (HighsInt i = 1; i < num_count; i++) {
    if (not_reported_ones && value_distribution.limit_[i - 1] >= 1.0) {
      count = value_distribution.num_one_;
      if (count) {
        percentage = doublePercentage(count, sum_count);
        int_percentage = percentage;
        highsLogDev(log_options, HighsLogType::kInfo,
                    "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                    "%%) are             %10.4g",
                    count, value_name.c_str(), int_percentage, 1.0);
        sum_report_count += count;
        if (mu > 0) {
          highsLogDev(log_options, HighsLogType::kInfo,
                      " corresponding to %10" HIGHSINT_FORMAT "\n", mu);
        } else {
          highsLogDev(log_options, HighsLogType::kInfo, "\n");
        }
      }
      not_reported_ones = false;
    }
    count = value_distribution.count_[i];
    if (count) {
      percentage = doublePercentage(count, sum_count);
      int_percentage = percentage;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                  "%%) in [%10.4g, %10.4g)",
                  count, value_name.c_str(), int_percentage,
                  value_distribution.limit_[i - 1],
                  value_distribution.limit_[i]);
      sum_report_count += count;
      if (mu > 0) {
        highsLogDev(log_options, HighsLogType::kInfo,
                    " corresponding to [%10" HIGHSINT_FORMAT
                    ", %10" HIGHSINT_FORMAT ")\n",
                    (HighsInt)(value_distribution.limit_[i - 1] * mu),
                    (HighsInt)(value_distribution.limit_[i] * mu));
      } else {
        highsLogDev(log_options, HighsLogType::kInfo, "\n");
      }
    }
  }
  if (not_reported_ones && value_distribution.limit_[num_count - 1] >= 1.0) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      int_percentage = percentage;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                  "%%) are             %10.4g",
                  count, value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "  corresponding to  %10" HIGHSINT_FORMAT "\n", mu);
      } else {
        highsLogDev(log_options, HighsLogType::kInfo, "\n");
      }
    }
    not_reported_ones = false;
  }
  count = value_distribution.count_[num_count];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    int_percentage = percentage;
    highsLogDev(log_options, HighsLogType::kInfo,
                "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                "%%) in [%10.4g,        inf)",
                count, value_name.c_str(), int_percentage,
                value_distribution.limit_[num_count - 1]);
    sum_report_count += count;
    if (mu > 0) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  " corresponding to [%10" HIGHSINT_FORMAT ",        inf)\n",
                  (HighsInt)(value_distribution.limit_[num_count - 1] * mu));
    } else {
      highsLogDev(log_options, HighsLogType::kInfo, "\n");
    }
  }
  if (not_reported_ones) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      int_percentage = percentage;
      highsLogDev(log_options, HighsLogType::kInfo,
                  "%12" HIGHSINT_FORMAT " %svalues (%3" HIGHSINT_FORMAT
                  "%%) are             %10.4g",
                  count, value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "  corresponding to  %10" HIGHSINT_FORMAT "\n", mu);
      } else {
        highsLogDev(log_options, HighsLogType::kInfo, "\n");
      }
    }
  }
  highsLogDev(log_options, HighsLogType::kInfo,
              "%12" HIGHSINT_FORMAT " %svalues\n", sum_count,
              value_name.c_str());
  if (sum_report_count != sum_count)
    highsLogDev(log_options, HighsLogType::kInfo,
                "ERROR: %" HIGHSINT_FORMAT
                " = sum_report_count != sum_count = %" HIGHSINT_FORMAT "\n",
                sum_report_count, sum_count);
  return true;
}

bool initialiseScatterData(const HighsInt max_num_point,
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
  if (value0 <= 0 || value1 <= 0) return false;
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
  HighsInt point_num = 0;
  for (HighsInt pass = 0; pass < 2; pass++) {
    HighsInt from_point;
    HighsInt to_point;
    if (pass == 0) {
      from_point = scatter_data.last_point_;
      to_point = std::min(scatter_data.num_point_, scatter_data.max_num_point_);
    } else {
      from_point = 0;
      to_point = scatter_data.last_point_;
    }
    for (HighsInt point = from_point; point < to_point; point++) {
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
  for (HighsInt point = 0; point < scatter_data.max_num_point_; point++) {
    double value0 = scatter_data.value0_[point];
    double value1 = scatter_data.value1_[point];
    double predicted_value1;
    if (predictFromScatterData(scatter_data, value0, predicted_value1, true)) {
      double error = fabs(predicted_value1 - value1);  // / fabs(value1);
      if (
          //	10*error > awful_regression_error &&
          print)
        printf("%5" HIGHSINT_FORMAT " %10.4g %10.4g %10.4g %10.4g\n", point,
               value0, value1, predicted_value1, error);
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
  for (HighsInt point = 0; point < scatter_data.max_num_point_; point++) {
    double value0 = scatter_data.value0_[point];
    double value1 = scatter_data.value1_[point];
    double predicted_value1;
    if (predictFromScatterData(scatter_data, value0, predicted_value1)) {
      double error = fabs(predicted_value1 - value1);  //  / fabs(value1);
      if (
          //	10*error > awful_regression_error &&
          print)
        printf("%5" HIGHSINT_FORMAT " %10.4g %10.4g %10.4g %10.4g\n", point,
               value0, value1, predicted_value1, error);
      sum_linear_error += error;
    }
  }
  if (print)
    printf("                                       %10.4g\n", sum_linear_error);
  scatter_data.log_regression_error_ = sum_log_error;
  scatter_data.linear_regression_error_ = sum_linear_error;
  return true;
}

bool printScatterData(const std::string& name,
                      const HighsScatterData& scatter_data) {
  if (!scatter_data.num_point_) return true;
  double x;
  double y;
  HighsInt point_num = 0;
  printf("%s scatter data\n", name.c_str());
  const HighsInt to_point =
      std::min(scatter_data.num_point_, scatter_data.max_num_point_);
  for (HighsInt point = scatter_data.last_point_ + 1; point < to_point;
       point++) {
    point_num++;
    x = scatter_data.value0_[point];
    y = scatter_data.value1_[point];
    printf("%" HIGHSINT_FORMAT ",%10.4g,%10.4g,%" HIGHSINT_FORMAT "\n", point,
           x, y, point_num);
  }
  for (HighsInt point = 0; point <= scatter_data.last_point_; point++) {
    point_num++;
    x = scatter_data.value0_[point];
    y = scatter_data.value1_[point];
    printf("%" HIGHSINT_FORMAT ",%10.4g,%10.4g,%" HIGHSINT_FORMAT "\n", point,
           x, y, point_num);
  }
  printf("Linear regression coefficients,%10.4g,%10.4g\n",
         scatter_data.linear_coeff0_, scatter_data.linear_coeff1_);
  printf("Log    regression coefficients,%10.4g,%10.4g\n",
         scatter_data.log_coeff0_, scatter_data.log_coeff1_);
  return true;
}

void printScatterDataRegressionComparison(
    const std::string& name, const HighsScatterData& scatter_data) {
  if (!scatter_data.num_error_comparison_) return;
  printf("\n%s scatter data regression\n", name.c_str());
  printf("%10" HIGHSINT_FORMAT " regression error comparisons\n",
         scatter_data.num_error_comparison_);
  printf("%10" HIGHSINT_FORMAT " regression awful  linear (>%10.4g)\n",
         scatter_data.num_awful_linear_, awful_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression awful  log    (>%10.4g)\n",
         scatter_data.num_awful_log_, awful_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression bad    linear (>%10.4g)\n",
         scatter_data.num_bad_linear_, bad_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression bad    log    (>%10.4g)\n",
         scatter_data.num_bad_log_, bad_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression fair   linear (>%10.4g)\n",
         scatter_data.num_fair_linear_, fair_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression fair   log    (>%10.4g)\n",
         scatter_data.num_fair_log_, fair_regression_error);
  printf("%10" HIGHSINT_FORMAT " regression better linear\n",
         scatter_data.num_better_linear_);
  printf("%10" HIGHSINT_FORMAT " regression better log\n",
         scatter_data.num_better_log_);
}

double nearestPowerOfTwoScale(const double value) {
  int exp_scale;
  // Decompose value into a normalized fraction and an integral power
  // of two.
  //
  // If arg is zero, returns zero and stores zero in *exp. Otherwise
  // (if arg is not zero), if no errors occur, returns the value x in
  // the range (-1;-0.5], [0.5; 1) and stores an integer value in *exp
  // such that x×2(*exp)=arg
  double check_x = std::frexp(value, &exp_scale);
  if (std::fabs(check_x) == 0.5) {
    check_x *= 2;
    exp_scale--;
  }
  const double check_value = check_x * std::pow(2, exp_scale);
  assert(check_value == value);
  exp_scale = -exp_scale;
  // Multiply a floating point value x(=1) by the number 2 raised to
  // the exp power
  double scale = std::ldexp(1, exp_scale);
  return scale;
}

bool termIsZero(const double val, const double lb, const double ub,
                const double epsilon) {
  const double absVal = std::abs(val);
  if (absVal <= std::numeric_limits<double>::min()) return true;
  return absVal <= epsilon &&
         absVal * std::max(std::abs(lb), std::abs(ub)) <= epsilon;
}

void highsAssert(const bool assert_condition, const std::string& message) {
  if (assert_condition) return;
  printf("Failing highsAssert(\"%s\")\n", message.c_str());
#ifdef NDEBUG
  // Standard assert won't trigger abort, so do it explicitly
  printf("assert(%s) failed ...\n", message.c_str());
  fflush(stdout);
  abort();
#else
  // Standard assert
  assert(assert_condition);
#endif
}

bool highsPause(const bool pause_condition, const std::string& message) {
  if (!pause_condition) return pause_condition;
  printf("Satisfying highsPause(\"%s\")\n", message.c_str());
  char str[100];
  printf("Enter any value to continue:");
  fflush(stdout);
  if (fgets(str, 100, stdin) != nullptr) {
    printf("You entered: \"%s\"\n", str);
    fflush(stdout);
  }
  return pause_condition;
}

bool handleExceptionIsOom(const HighsLogOptions& log_options,
                          const std::string& source,
                          const std::exception& exception) {
  highsLogUser(log_options, HighsLogType::kError, "Exception %s in %s\n",
               exception.what(), source.c_str());
  return dynamic_cast<const std::bad_alloc*>(&exception);
}
