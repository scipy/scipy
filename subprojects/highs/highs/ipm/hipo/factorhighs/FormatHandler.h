#ifndef FACTORHIGHS_FORMAT_HANDLER_H
#define FACTORHIGHS_FORMAT_HANDLER_H

#include <vector>

#include "FactorHighsOptions.h"
#include "FactorHighsSettings.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Interface class to handle different formats of dense matrices during the
// factorise phase.
// Any implementation of a specific format needs to define:
// - initFrontal: to initialise the frontal matrix with the correct number of
//                elements; the entries should be set to zero.
// - initClique: to initialise the clique matrix with the correct number of
//                elements; the entries should be left uninitialised.
// - assembleFrontal: to set a specific entry of frontal (used to assemble the
//                original matrix)
// - assembleFrontalMultiple: to sum a given number of consecutive entries into
//                frontal (used to assemble the child supernodes)
// - assembleClique: to sum the contributions of a given child supernode into
//                clique.
// - denseFactorise: to perform the dense partial factorisation of frontal,
//                storing the Schur complement in clique.

class FormatHandler {
 protected:
  // symbolic object
  const Symbolic* S_;

  // supernode being processed
  const Int sn_{};

  // size of the supernode
  const Int sn_size_{};

  // size of the front
  const Int ldf_{};

  // size of the clique
  const Int ldc_{};

  const FHoptions& FH_opt_;

  // local copies to be moved at the end
  std::vector<double>& frontal_;
  std::vector<double> clique_{};
  double* clique_ptr_ = nullptr;
  std::vector<double> local_reg_{};
  std::vector<Int> swaps_{};
  std::vector<double> pivot_2x2_{};

 public:
  FormatHandler(const Symbolic& S, Int sn, std::vector<double>& frontal,
                double* clique_ptr, const FHoptions& FH_opt);
  void terminate(std::vector<double>& clique, std::vector<double>& total_reg,
                 std::vector<Int>& swaps, std::vector<double>& pivot_2x2);

  // avoid copies
  FormatHandler(const FormatHandler&) = delete;
  FormatHandler& operator=(const FormatHandler&) = delete;

  // virtual destructor
  virtual ~FormatHandler() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void initFrontal() = 0;
  virtual void initClique() = 0;
  virtual void assembleFrontal(Int i, Int j, double val) = 0;
  virtual void assembleChild(Int child_sn, const double* child) = 0;
  virtual Int denseFactorise(double reg_thresh) = 0;

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual void extremeEntries() {}
};

const Int extra_space_frontal = 10;

}  // namespace hipo

#endif