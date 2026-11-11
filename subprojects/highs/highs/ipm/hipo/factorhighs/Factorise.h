#ifndef FACTORHIGHS_FACTORISE_H
#define FACTORHIGHS_FACTORISE_H

#include <cmath>

#include "CliqueStack.h"
#include "FactorHighsOptions.h"
#include "Numeric.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

class Factorise {
 public:
  // matrix to factorise
  std::vector<Int> rowsM_{};
  std::vector<Int> ptrM_{};
  std::vector<double> valM_{};
  Int n_{};
  Int nzM_{};

  // symbolic factorisation
  const Symbolic& S_;

  // children in supernodal elimination tree
  std::vector<Int> first_child_{};
  std::vector<Int> next_child_{};

  // reverse linked lists of chidlren
  std::vector<Int> first_child_reverse_{};
  std::vector<Int> next_child_reverse_{};

  // generated elements, aka Schur complements.
  std::vector<std::vector<double>> schur_contribution_{};

  // columns of L, stored as dense supernodes
  // This memory is managed outside of Factorise, so that it can be reused for
  // all ipm iterations.
  std::vector<std::vector<double>>& sn_columns_;

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<Int>> swaps_{};

  // Information about 2x2 pivots.
  // If pivot_2x2[sn][i] == 0, 1x1 pivot was used.
  // If pivot_2x2[sn][i] != 0, 2x2 pivot was used and pivot_2x2[sn][i] stores
  //  the off-diagonal pivot entry (of the 2x2 inverse).
  std::vector<std::vector<double>> pivot_2x2_{};

  // largest diagonal element in the original matrix and norms of columns
  double max_diag_{};
  double min_diag_{};
  double M_norm1_{};

  // regularisation
  std::vector<double> total_reg_{};

  // flag to stop computation
  std::atomic<bool> flag_stop_{false};

  const Logger* logger_;
  DataCollector& data_;
  const FHoptions& FH_opt_;

  CliqueStack* stack_;

  void permute(const std::vector<Int>& iperm);
  void processTreeSerial();
  void processTreeParallel();
  void processSupernode(Int sn);

 public:
  Factorise(const Symbolic& S, Int n, Int nz, const Int* rowsM, const Int* ptrM,
            const double* valM, const FHoptions& FH_opt, const Logger* logger,
            DataCollector& data, std::vector<std::vector<double>>& sn_columns,
            CliqueStack* stack);

  bool run(Numeric& num);
};

}  // namespace hipo

#endif
