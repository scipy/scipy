#include "DataCollector.h"

#include "FactorHighsSettings.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

// Functions to manage DataCollector

DataCollector::DataCollector() {
  times_.resize(kTimeSize);
  blas_calls_.resize(kTimeBlasEnd - kTimeBlasStart + 1);
}
void DataCollector::append() {
  // add an empty IterData object to the record
  iter_data_record_.push_back(IterData());
}
IterData& DataCollector::back() {
  // access most recent record of data
  if (iter_data_record_.empty()) append();
  return iter_data_record_.back();
}

//
// Data-collecting functions
//

void DataCollector::sumTime(TimeItems i, double t) {
#if HIPO_TIMING_LEVEL > 0
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(mutex_);
  times_[i] += t;
#if HIPO_TIMING_LEVEL >= 3
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++blas_calls_[i - kTimeBlasStart];
#endif
#endif
}
void DataCollector::setExtremeEntries(double minD, double maxD, double minoffD,
                                      double maxoffD) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(mutex_);
  back().minD = std::min(back().minD, minD);
  back().maxD = std::max(back().maxD, maxD);
  back().minL = std::min(back().minL, minoffD);
  back().maxL = std::max(back().maxL, maxoffD);
#endif
}
void DataCollector::countRegPiv() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Increase the number of dynamically regularised pivots.
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_reg_piv;
#endif
}
void DataCollector::countSwap() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_swap;
#endif
}
void DataCollector::count2x2() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_2x2;
#endif
}
void DataCollector::setWrongSign(double p) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_wrong_sign;
  back().max_wrong_sign = std::max(back().max_wrong_sign, std::abs(p));
#endif
}
void DataCollector::setMaxReg(double new_reg) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Keep track of maximum regularisation used.
  std::lock_guard<std::mutex> lock(mutex_);
  back().max_reg = std::max(back().max_reg, new_reg);
#endif
}

void DataCollector::setNorms(double norm1, double maxdiag) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  back().M_norm1 = norm1;
  back().M_maxdiag = maxdiag;
#endif
}

//
// Printing
//

void DataCollector::printTimes(const Logger& logger) const {
#if HIPO_TIMING_LEVEL >= 1

  std::stringstream log_stream;

  log_stream << "----------------------------------------------------\n";
  log_stream << "Analyse time            \t" << fix(times_[kTimeAnalyse], 8, 4)
             << '\n';

#if HIPO_TIMING_LEVEL >= 2

  log_stream << "\tTree:                   "
             << fix(times_[kTimeAnalyseTree], 8, 4) << " ("
             << fix(times_[kTimeAnalyseTree] / times_[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tCounts:                 "
             << fix(times_[kTimeAnalyseCount], 8, 4) << " ("
             << fix(times_[kTimeAnalyseCount] / times_[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
  log_stream << "\tSupernodes:             "
             << fix(times_[kTimeAnalyseSn], 8, 4) << " ("
             << fix(times_[kTimeAnalyseSn] / times_[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tReorder:                "
             << fix(times_[kTimeAnalyseReorder], 8, 4) << " ("
             << fix(times_[kTimeAnalyseReorder] / times_[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
  log_stream << "\tSn sparsity pattern:    "
             << fix(times_[kTimeAnalysePattern], 8, 4) << " ("
             << fix(times_[kTimeAnalysePattern] / times_[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
  log_stream << "\tRelative indices:       "
             << fix(times_[kTimeAnalyseRelInd], 8, 4) << " ("
             << fix(times_[kTimeAnalyseRelInd] / times_[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
#endif

  log_stream << "----------------------------------------------------\n";
  log_stream << "Factorise time          \t"
             << fix(times_[kTimeFactorise], 8, 4) << "\n";

#if HIPO_TIMING_LEVEL >= 2
  log_stream << "\tPrepare fact:           "
             << fix(times_[kTimeFactorisePrepare], 8, 4) << " ("
             << fix(times_[kTimeFactorisePrepare] / times_[kTimeFactorise] *
                        100,
                    4, 1)
             << "%)\n";
  log_stream << "\tAssemble original:      "
             << fix(times_[kTimeFactoriseAssembleOriginal], 8, 4) << " ("
             << fix(times_[kTimeFactoriseAssembleOriginal] /
                        times_[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tAssemble children in F: "
             << fix(times_[kTimeFactoriseAssembleChildrenFrontal], 8, 4) << " ("
             << fix(times_[kTimeFactoriseAssembleChildrenFrontal] /
                        times_[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tAssemble children in C: "
             << fix(times_[kTimeFactoriseAssembleChildrenClique], 8, 4) << " ("
             << fix(times_[kTimeFactoriseAssembleChildrenClique] /
                        times_[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tDense factorisation:    "
             << fix(times_[kTimeFactoriseDenseFact], 8, 4) << " ("
             << fix(times_[kTimeFactoriseDenseFact] / times_[kTimeFactorise] *
                        100,
                    4, 1)
             << "%)\n";

  log_stream << "\tTerminate:              "
             << fix(times_[kTimeFactoriseTerminate], 8, 4) << " ("
             << fix(times_[kTimeFactoriseTerminate] / times_[kTimeFactorise] *
                        100,
                    4, 1)
             << "%)\n";

  log_stream << "\t\tmain:           " << fix(times_[kTimeDenseFact_main], 8, 4)
             << "\n";
  log_stream << "\t\tSchur:          "
             << fix(times_[kTimeDenseFact_schur], 8, 4) << "\n";
  log_stream << "\t\tkernel:         "
             << fix(times_[kTimeDenseFact_kernel], 8, 4) << "\n";
  log_stream << "\t\tconvert:        "
             << fix(times_[kTimeDenseFact_convert], 8, 4) << "\n";
  log_stream << "\t\tpivoting:       "
             << fix(times_[kTimeDenseFact_pivoting], 8, 4) << "\n";
#endif

  log_stream << "----------------------------------------------------\n";
  log_stream << "Solve time              \t" << fix(times_[kTimeSolve], 8, 4)
             << "\n";

#if HIPO_TIMING_LEVEL >= 2
  log_stream << "\tPrepare solve:          "
             << fix(times_[kTimeSolvePrepare], 8, 4) << " ("
             << fix(times_[kTimeSolvePrepare] / times_[kTimeSolve] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tSolve:                  "
             << fix(times_[kTimeSolveSolve], 8, 4) << " ("
             << fix(times_[kTimeSolveSolve] / times_[kTimeSolve] * 100, 4, 1)
             << "%)\n";
  log_stream << "\t\tdense:          "
             << fix(times_[kTimeSolveSolve_dense], 8, 4) << "\n";
  log_stream << "\t\tsparse:         "
             << fix(times_[kTimeSolveSolve_sparse], 8, 4) << "\n";
  log_stream << "\t\tswap:           "
             << fix(times_[kTimeSolveSolve_swap], 8, 4) << "\n";
#endif
  log_stream << "----------------------------------------------------\n";

#if HIPO_TIMING_LEVEL >= 3

  double total_blas_time =
      times_[kTimeBlas_copy] + times_[kTimeBlas_axpy] + times_[kTimeBlas_scal] +
      times_[kTimeBlas_swap] + times_[kTimeBlas_gemv] + times_[kTimeBlas_trsv] +
      times_[kTimeBlas_tpsv] + times_[kTimeBlas_ger] + times_[kTimeBlas_trsm] +
      times_[kTimeBlas_syrk] + times_[kTimeBlas_gemm];

  log_stream << "BLAS time               \t" << fix(total_blas_time, 8, 4)
             << '\n';
  log_stream << "\tcopy:           \t" << fix(times_[kTimeBlas_copy], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_copy] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_copy - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\taxpy:           \t" << fix(times_[kTimeBlas_axpy], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_axpy] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_axpy - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tscal:           \t" << fix(times_[kTimeBlas_scal], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_scal] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_scal - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tswap:           \t" << fix(times_[kTimeBlas_swap], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_swap] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_swap - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tgemv:           \t" << fix(times_[kTimeBlas_gemv], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_gemv] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_gemv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttrsv:           \t" << fix(times_[kTimeBlas_trsv], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_trsv] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_trsv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttpsv:           \t" << fix(times_[kTimeBlas_tpsv], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_tpsv] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_tpsv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tger:            \t" << fix(times_[kTimeBlas_ger], 8, 4)
             << " (" << fix(times_[kTimeBlas_ger] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_ger - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttrsm:           \t" << fix(times_[kTimeBlas_trsm], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_trsm] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_trsm - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tsyrk:           \t" << fix(times_[kTimeBlas_syrk], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_syrk] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_syrk - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tgemm:           \t" << fix(times_[kTimeBlas_gemm], 8, 4)
             << " ("
             << fix(times_[kTimeBlas_gemm] / total_blas_time * 100, 4, 1)
             << "%) in "
             << integer(blas_calls_[kTimeBlas_gemm - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "----------------------------------------------------\n";
#endif

  logger.print(log_stream.str().c_str());
#endif
}

void DataCollector::printIter(const Logger& logger) const {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA

  std::stringstream log_stream;

  log_stream << "\niter |   min D     max D     min L     max L  |"
                "    reg   swap    2x2    ws | "
                " max_reg   max_ws    norm1     maxdiag|\n";

  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];

    log_stream << integer(i, 3) << "  |";
    log_stream << sci(iter.minD, 9, 1) << " ";
    log_stream << sci(iter.maxD, 9, 1) << " ";
    log_stream << sci(iter.minL, 9, 1) << " ";
    log_stream << sci(iter.maxL, 9, 1) << " |";
    log_stream << integer(iter.n_reg_piv, 6) << " ";
    log_stream << integer(iter.n_swap, 6) << " ";
    log_stream << integer(iter.n_2x2, 6) << " ";
    log_stream << integer(iter.n_wrong_sign, 6) << " |";
    log_stream << sci(iter.max_reg, 9, 1) << " ";
    log_stream << sci(iter.max_wrong_sign, 9, 1) << " ";
    log_stream << sci(iter.M_norm1, 9, 1) << " ";
    log_stream << sci(iter.M_maxdiag, 9, 1) << "|\n";
  }
  logger.print(log_stream.str().c_str());
#endif
}

}  // namespace hipo