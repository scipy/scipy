#include "HybridHybridFormatHandler.h"

#include <cassert>
#include <cstring>

#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

HybridHybridFormatHandler::HybridHybridFormatHandler(
    const Symbolic& S, Int sn, DataCollector& data,
    std::vector<double>& frontal, double* clique_ptr, const FHoptions& FH_opt)
    : FormatHandler(S, sn, frontal, clique_ptr, FH_opt),
      data_{data},
      nb_{FH_opt_.nb} {
  // initialise frontal and clique
  initFrontal();

  // if CliqueStack is used, clique_ptr already points to a valid region of
  // memory for the clique. Otherwise, allocate it locally.
  if (!clique_ptr_) initClique();
}

void HybridHybridFormatHandler::initFrontal() {
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;
  diag_start_.resize(n_blocks);
  Int64 frontal_size =
      getDiagStart(ldf_, sn_size_, nb_, n_blocks, diag_start_) +
      extra_space_frontal;
  frontal_.resize(frontal_size);
  std::memset(frontal_.data(), 0, frontal_size * sizeof(double));

  // NB: extra_space_frontal is not strictly needed. However, it removes some
  // weird problem on windows in debug. Who knows what's happening...

  // frontal_ is actually allocated just the first time, then the memory is
  // reused from the previous factorisations and just initialised.
}

void HybridHybridFormatHandler::initClique() {
  clique_.resize(S_->cliqueSize(sn_));

  // If the clique size is zero, do not access the underlying pointer. This
  // causes strange issues on windows. It's not a problem if clique_ptr_ remains
  // null, because it will never be used in that case.
  if (!clique_.empty()) clique_ptr_ = clique_.data();
}

void HybridHybridFormatHandler::assembleFrontal(Int i, Int j, double val) {
  Int block_id = j / nb_;
  Int ldb = ldf_ - block_id * nb_;
  Int i_in_block = i - block_id * nb_;
  Int j_in_block = j - block_id * nb_;
  frontal_[diag_start_[block_id] + i_in_block + ldb * j_in_block] = val;
}

void HybridHybridFormatHandler::assembleFrontalMultiple(
    Int& num, const double* child_data, Int child_size, Int child_sn, Int row_c,
    Int col_c, Int row_f, Int col_f) {
  const Int block_child_id = col_c / nb_;
  const Int jb_child = std::min(nb_, child_size - nb_ * block_child_id);
  const Int row_c_local = row_c - block_child_id * nb_;
  const Int col_c_local = col_c - block_child_id * nb_;
  const Int64 start_block_child =
      S_->cliqueBlockStart(child_sn, block_child_id);

  Int block_frontal_id = col_f / nb_;
  Int ldb = ldf_ - block_frontal_id * nb_;
  Int row_f_local = row_f - block_frontal_id * nb_;
  Int col_f_local = col_f - block_frontal_id * nb_;

  if (num > kMinConsecutiveSums)
    callAndTime_daxpy(
        num, 1.0,
        &child_data[start_block_child + col_c_local + jb_child * row_c_local],
        jb_child,
        &frontal_[diag_start_[block_frontal_id] + row_f_local +
                  ldb * col_f_local],
        1, data_);
  else {
    frontal_[diag_start_[block_frontal_id] + row_f_local + ldb * col_f_local] +=
        child_data[start_block_child + col_c_local + jb_child * row_c_local];
    num = 1;
  }
}

Int HybridHybridFormatHandler::denseFactorise(double reg_thresh) {
  Int status;

  // either clique is valid, or clique is not needed
  assert(clique_ptr_ || ldf_ == sn_size_);

  status = denseFactFP2FH(frontal_.data(), ldf_, sn_size_, nb_, data_);
  if (status) return status;

  // find the position within pivot_sign corresponding to this supernode
  Int sn_start = S_->snStart(sn_);
  const Int* pivot_sign = &S_->pivotSign().data()[sn_start];

  status = denseFactFH('H', ldf_, sn_size_, frontal_.data(), clique_ptr_,
                       pivot_sign, reg_thresh, local_reg_.data(), swaps_.data(),
                       pivot_2x2_.data(), data_, FH_opt_);

  return status;
}

void HybridHybridFormatHandler::assembleClique(const double* child_data,
                                               Int child_size, Int child_sn) {
  // assemble the child clique into the current clique by blocks of columns.
  // within a block, assemble by rows.

  const Int blocks = (child_size - 1) / nb_ + 1;

  Int row_start{};

  // go through the blocks of columns of the child sn
  for (Int block = 0; block < blocks; ++block) {
    const Int64 block_start = S_->cliqueBlockStart(child_sn, block);

    const Int col_start = row_start;
    const Int col_end = std::min(col_start + nb_, child_size);

    // go through the rows within this block
    for (Int row = row_start; row < child_size; ++row) {
      const Int row_clique = S_->relindClique(child_sn, row) - sn_size_;

      // already assembled into frontal
      if (row_clique < 0) continue;

      // go through the columns of the block
      Int col = col_start;
      while (col < col_end) {
        Int col_clique = S_->relindClique(child_sn, col);
        if (col_clique < sn_size_) {
          ++col;
          continue;
        }
        col_clique -= sn_size_;

        // information and sizes of child sn
        const Int jb_child = std::min(nb_, child_size - nb_ * block);
        const Int row_local = row - block * nb_;
        const Int col_local = col - block * nb_;

        // sun consecutive entries in a row.
        // consecutive need to be reduced, to account for edge of the block
        const Int zeros_stored_row =
            std::max((Int)0, jb_child - (row - row_start) - 1);
        Int consecutive = S_->consecutiveSums(child_sn, col);
        const Int left_in_child = col_end - col - zeros_stored_row;
        consecutive = std::min(consecutive, left_in_child);

        // consecutive need to account also for edge of block in parent
        const Int block_in_parent = col_clique / nb_;
        const Int col_end_parent = std::min((block_in_parent + 1) * nb_, ldc_);
        const Int left_in_parent = col_end_parent - col_clique;
        consecutive = std::min(consecutive, left_in_parent);

        // needed to deal with zeros stored in upper right part of block
        if (consecutive == 0) break;

        // information and sizes of current sn
        const Int jb_clique = std::min(nb_, ldc_ - nb_ * block_in_parent);
        const Int row_clique_local = row_clique - block_in_parent * nb_;
        const Int col_clique_local = col_clique - block_in_parent * nb_;
        const Int64 block_start_clique =
            S_->cliqueBlockStart(sn_, block_in_parent);

        if (consecutive > kMinConsecutiveSums) {
          callAndTime_daxpy(
              consecutive, 1.0,
              &child_data[block_start + col_local + jb_child * row_local], 1,
              &clique_ptr_[block_start_clique + col_clique_local +
                           jb_clique * row_clique_local],
              1, data_);

          col += consecutive;
        } else {
          clique_ptr_[block_start_clique + col_clique_local +
                      jb_clique * row_clique_local] +=
              child_data[block_start + col_local + jb_child * row_local];

          col++;
        }
      }
    }

    row_start += nb_;
  }
}

void HybridHybridFormatHandler::assembleChild(Int child_sn,
                                              const double* child_data) {
  HIPO_CLOCK_CREATE;

  const Int child_begin = S_->snStart(child_sn);
  const Int child_end = S_->snStart(child_sn + 1);
  const Int child_sn_size = child_end - child_begin;
  const Int child_clique_size =
      S_->ptr(child_sn + 1) - S_->ptr(child_sn) - child_sn_size;

  // ASSEMBLE INTO FRONTAL
  HIPO_CLOCK_START(2);
  // go through the columns of the contribution of the child
  for (Int col = 0; col < child_clique_size; ++col) {
    // relative index of column in the frontal matrix
    Int col_f = S_->relindClique(child_sn, col);

    if (col_f < sn_size_) {
      // assemble into frontal

      // go through the rows of the contribution of the child
      Int row = col;
      while (row < child_clique_size) {
        // relative index of the entry in the matrix frontal
        const Int row_f = S_->relindClique(child_sn, row);

        // how many entries to sum
        Int consecutive = S_->consecutiveSums(child_sn, row);

        assembleFrontalMultiple(consecutive, child_data, child_clique_size,
                                child_sn, row, col, row_f, col_f);

        row += consecutive;
      }
    } else
      break;
  }
  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenFrontal);

  // ASSEMBLE INTO CLIQUE
  HIPO_CLOCK_START(2);
  assembleClique(child_data, child_clique_size, child_sn);
  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenClique);
}

void HybridHybridFormatHandler::extremeEntries() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  double minD = 1e100;
  double maxD = 0.0;
  double minoffD = 1e100;
  double maxoffD = 0.0;

  // number of blocks of columns
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;

  // index to access frontal
  Int64 index{};

  // go through blocks of columns for this supernode
  for (Int j = 0; j < n_blocks; ++j) {
    // number of columns in the block
    const Int jb = std::min(nb_, sn_size_ - nb_ * j);

    for (Int k = 0; k < jb; ++k) {
      // off diagonal entries
      for (Int i = 0; i < k; ++i) {
        if (frontal_[index] != 0.0) {
          minoffD = std::min(minoffD, std::abs(frontal_[index]));
          maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
        }
        index++;
      }

      // diagonal entry
      minD = std::min(minD, std::abs(1.0 / frontal_[index]));
      maxD = std::max(maxD, std::abs(1.0 / frontal_[index]));

      index += jb - k;
    }

    const Int64 entries_left = (Int64)(ldf_ - nb_ * j - jb) * jb;

    for (Int64 i = 0; i < entries_left; ++i) {
      if (frontal_[index] != 0.0) {
        minoffD = std::min(minoffD, std::abs(frontal_[index]));
        maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
      }
      index++;
    }
  }

  data_.setExtremeEntries(minD, maxD, minoffD, maxoffD);
#endif
}

}  // namespace hipo