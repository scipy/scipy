#include "FormatHandler.h"

#include "CallAndTimeBlas.h"
#include "DenseFact.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

FormatHandler::FormatHandler(const Symbolic& S, Int sn,
                             std::vector<double>& frontal, double* clique_ptr,
                             const FHoptions& FH_opt)
    : S_{&S},
      sn_{sn},
      sn_size_{S_->snStart(sn_ + 1) - S_->snStart(sn_)},
      ldf_{(Int)(S_->ptr(sn_ + 1) - S_->ptr(sn_))},
      ldc_{ldf_ - sn_size_},
      FH_opt_{FH_opt},
      frontal_{frontal},
      clique_ptr_{clique_ptr} {
  local_reg_.resize(sn_size_);
  swaps_.resize(sn_size_);
  pivot_2x2_.resize(sn_size_);
}

void FormatHandler::terminate(std::vector<double>& clique,
                              std::vector<double>& total_reg,
                              std::vector<Int>& swaps,
                              std::vector<double>& pivot_2x2) {
  // Move local copies of data into their final position.
  // In this way, the shared objects sn_columns_ and schur_contribution_ are
  // accessed only here, while a local copy is used for the assembly and dense
  // factorisation. This should avoid the problem of false sharing.

  clique = std::move(clique_);
  swaps = std::move(swaps_);
  pivot_2x2 = std::move(pivot_2x2_);

  // Move local regularisation into total regularisation.
  for (Int i = 0; i < sn_size_; ++i)
    total_reg[S_->snStart(sn_) + i] = local_reg_[i];

  // This function should not require a lock, since all threads access different
  // locations within total_reg, sn_columns and schur_contribution.
}

}  // namespace hipo
