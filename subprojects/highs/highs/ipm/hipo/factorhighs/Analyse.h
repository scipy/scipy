#ifndef FACTORHIGHS_ANALYSE_H
#define FACTORHIGHS_ANALYSE_H

#include <algorithm>
#include <vector>

#include "DataCollector.h"
#include "FactorHighsOptions.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

// Class to perform the analyse phase of the factorisation.
// The final symbolic factorisation is stored in an object of type Symbolic.
class Analyse {
  bool ready_ = false;

  // Matrix to be factorised, stored in upper and lower triangular format
  std::vector<Int> rows_upper_{};
  std::vector<Int> ptr_upper_{};
  std::vector<Int> rows_lower_{};
  std::vector<Int> ptr_lower_{};

  // info about matrix and factor
  Int n_{};
  Int nz_{};
  Int64 nz_factor_{};
  double dense_ops_{};
  double dense_ops_norelax_{};
  double sparse_ops_{};
  double critical_ops_{};
  std::vector<Int> signs_{};

  // Permutation and inverse permutation
  std::vector<Int> perm_{};
  std::vector<Int> iperm_{};

  // Elimination tree
  std::vector<Int> parent_{};

  // postorder of the elimination tree
  std::vector<Int> postorder_{};

  // number of entries in each column of L
  std::vector<Int> col_count_{};

  // sparsity pattern of supernodes of L
  std::vector<Int> rows_sn_{};
  std::vector<Int64> ptr_sn_{};

  std::vector<Int> sn_indices_{};

  // fundamental supernodes information
  Int sn_count_{};
  Int64 artificial_nz_{};
  std::vector<Int> sn_belong_{};
  std::vector<Int> sn_start_{};
  std::vector<Int> sn_parent_{};

  // temporary storage for relaxing supernodes
  std::vector<Int64> fake_nz_{};
  std::vector<Int> merged_into_{};
  Int merged_sn_{};

  // relative indices of original columns wrt L columns
  std::vector<Int> relind_cols_{};

  // relative indices of clique wrt parent
  std::vector<std::vector<Int>> relind_clique_{};

  // information about consecutive indices in relindClique
  std::vector<std::vector<Int>> consecutive_sums_{};

  // estimate of maximum storage
  double serial_storage_{};

  Int64 max_stack_size_{};

  std::vector<std::vector<Int64>> clique_block_start_{};

  const FHoptions& FH_opt_;

  const Logger* logger_;
  DataCollector& data_;

  TreeSchedule schedule_solve_;

  double ops_solve_;
  double critical_ops_solve_;
  std::vector<double> task_ops_solve_;

  // Functions to perform analyse phase
  void permute(const std::vector<Int>& iperm);
  void eTree();
  void postorder();
  void colCount();
  void fundamentalSupernodes();
  void relaxSupernodes();
  void relaxSnMaxNz(double& flops, double& spops);
  double doRelaxSnMaxNz(double& flops, double& spops, Int64 max_artificial_nz);
  void afterRelaxSn();
  void relaxSnNetOps(double& flops, double& spops);
  void snPattern();
  void relativeIndCols();
  void relativeIndClique();
  void reorderChildren();
  void computeStorage(Int fr, Int sz, Int64& fr_entries,
                      Int64& cl_entries) const;
  void computeCriticalPath();
  void computeCriticalPathSolve();
  void computeBlockStart();
  void computeStackSize();
  void computeTreeScheduleSolve();
  Int checkOverflow() const;

 public:
  // Constructor: matrix must be in lower triangular format
  Analyse(Int n, Int nz, const Int* rows, const Int* ptr, const Int* signs,
          const FHoptions& FH_opt, const Logger* logger, DataCollector& data,
          const Int* perm);

  // Run analyse phase and save the result in Symbolic object S
  Int run(Symbolic& S);
};

}  // namespace hipo

#endif