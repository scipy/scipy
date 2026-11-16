#include "Symbolic.h"

#include <iostream>

#include "FactorHighsSettings.h"

namespace hipo {

void TreeSchedule::clear() {
  sn_per_task.clear();
  task_parent.clear();
  valid = false;
}

Int64 Symbolic::nz() const { return nz_; }
double Symbolic::flops() const { return flops_; }
double Symbolic::spops() const { return spops_; }
double Symbolic::critops() const { return critops_; }
Int Symbolic::size() const { return n_; }
Int Symbolic::sn() const { return sn_; }
double Symbolic::fillin() const { return fillin_; }
Int Symbolic::rows(Int64 i) const { return rows_[i]; }
Int64 Symbolic::ptr(Int i) const { return ptr_[i]; }
Int Symbolic::snStart(Int i) const { return sn_start_[i]; }
Int Symbolic::snParent(Int i) const { return sn_parent_[i]; }
Int Symbolic::relindCols(Int i) const { return relind_cols_[i]; }
Int Symbolic::relindClique(Int i, Int j) const { return relind_clique_[i][j]; }
Int Symbolic::consecutiveSums(Int i, Int j) const {
  return consecutive_sums_[i][j];
}
Int64 Symbolic::cliqueBlockStart(Int sn, Int bl) const {
  return clique_block_start_[sn][bl];
}
Int64 Symbolic::cliqueSize(Int sn) const {
  return clique_block_start_[sn].back();
}
Int64 Symbolic::maxStackSize() const { return max_stack_size_; }
Int Symbolic::largestFront() const { return largest_front_; }
double Symbolic::storage() const { return serial_storage_; }
Int Symbolic::depth() const { return tree_depth_; }
double Symbolic::solveTreeSpeedup() const {
  return ops_solve_ / critops_solve_;
}

const std::vector<Int64>& Symbolic::ptr() const { return ptr_; }
const std::vector<Int>& Symbolic::iperm() const { return iperm_; }
const std::vector<Int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<Int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<Int>& Symbolic::pivotSign() const { return pivot_sign_; }
const TreeSchedule& Symbolic::schedule() const { return schedule_solve_; }

static std::string memoryString(double mem) {
  std::stringstream ss;

  if (mem < 1024)
    ss << sci(mem, 0, 1) << " B";
  else if (mem < 1024 * 1024)
    ss << sci(mem / 1024, 0, 1) << " KB";
  else if (mem < 1024 * 1024 * 1024)
    ss << sci(mem / 1024 / 1024, 0, 1) << " MB";
  else
    ss << sci(mem / 1024 / 1024 / 1024, 0, 1) << " GB";

  return ss.str();
}

void Symbolic::print(const Logger& logger, bool verbose) const {
  std::stringstream log_stream;
  log_stream << "\nFactorisation statistics\n";
  log_stream << textline("Size:") << sci(n_, 0, 2) << '\n';
  log_stream << textline("Nnz:") << sci(nz_, 0, 2) << '\n';
  log_stream << textline("Fill-in:") << fix(fillin_, 0, 2) << '\n';
  log_stream << textline("Serial memory:") << memoryString(serial_storage_)
             << '\n';
  log_stream << textline("Flops:") << sci(flops_, 0, 1) << '\n';
  if (verbose) {
    log_stream << textline("Sparse ops:") << sci(spops_, 0, 1) << '\n';
    log_stream << textline("Critical ops:") << sci(critops_, 0, 1) << '\n';
    log_stream << textline("Max tree speedup:") << fix(flops_ / critops_, 0, 2)
               << '\n';
    log_stream << textline("Tree depth:") << integer(tree_depth_, 0) << '\n';
    log_stream << textline("Number of solve tasks:")
               << integer(schedule_solve_.count(), 0) << '\n';
    log_stream << textline("Solve tree speedup:")
               << fix(ops_solve_ / critops_solve_, 0, 2) << '\n';
    log_stream << textline("Artificial nz:") << sci(artificial_nz_, 0, 1)
               << '\n';
    log_stream << textline("Artificial ops:") << sci(artificial_ops_, 0, 1)
               << '\n';
    log_stream << textline("Largest front:") << integer(largest_front_, 0)
               << '\n';
    log_stream << textline("Largest supernode:") << integer(largest_sn_, 0)
               << '\n';
    log_stream << textline("Supernodes:") << integer(sn_, 0) << '\n';
    log_stream << textline("Sn size <= 1:") << integer(sn_size_1_, 0) << '\n';
    log_stream << textline("Sn size <= 10:") << integer(sn_size_10_, 0) << '\n';
    log_stream << textline("Sn size <= 100:") << integer(sn_size_100_, 0)
               << '\n';
    log_stream << textline("Sn avg size:") << sci((double)n_ / sn_, 0, 1)
               << '\n';
  }
  logger.print(log_stream.str().c_str());
}

}  // namespace hipo