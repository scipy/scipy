#ifndef HIPO_LINEAR_SOLVER_H
#define HIPO_LINEAR_SOLVER_H

#include <vector>

#include "Model.h"
#include "Options.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "model/HighsHessian.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

// Interface class for solving augmented system or normal equations.
//
// Any linear solver needs to define the functions:
// - factorAS: factorise the augmented system
// - solveAS: solve a linear system with the augmented system
// - factorNE: factorise the normal equations
// - solveNE: solve a linear system with the normal equations
// - clear: reset the data structure for the next factorisation.
//
// The linear solver may also define functions:
// - setup: perform any preliminary calculation (e.g. symbolic factorisation)
// - flops: return number of flops needed for factorisation
// - spops: return number of sparse ops needed for factorisation
// - nz: return number of nonzeros in factorisation
// - getReg: extract regularisation
//
// NB: forming the normal equations or augmented system is delegated to the
// linear solver chosen, so that only the appropriate data (upper triangle,
// lower triangle, or else) is constructed.

enum LinearSolverType { kFactorHighsType = 1, kUpLookingType };

class LinearSolver {
 public:
  bool valid_ = false;

  // default constructor
  LinearSolver() = default;

  // avoid copies
  LinearSolver(const LinearSolver&) = delete;
  LinearSolver& operator=(const LinearSolver&) = delete;

  // virtual destructor
  virtual ~LinearSolver() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual Int factorAS(const std::vector<double>& scaling) = 0;

  virtual Int solveAS(const std::vector<double>& rhs_x,
                      const std::vector<double>& rhs_y,
                      std::vector<double>& lhs_x,
                      std::vector<double>& lhs_y) = 0;

  virtual Int factorNE(const std::vector<double>& scaling) = 0;

  virtual Int solveNE(const std::vector<double>& rhs,
                      std::vector<double>& lhs) = 0;

  virtual void clear() = 0;

  virtual LinearSolverType type() const = 0;

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual Int setup() { return 0; }
  virtual double flops() const { return 0; }
  virtual double spops() const { return 0; }
  virtual double nz() const { return 0; }
  virtual void getReg(std::vector<double>& reg){};
};

inline std::string LinearSolverTypeToString(LinearSolverType t) {
  static const std::map<LinearSolverType, std::string> type_map{
      {kFactorHighsType, "FactorHighs"},
      {kUpLookingType, "UpLooking"},
  };
  
  auto found = type_map.find(t);
  if (found != type_map.end()) return found->second;
  return "XUnknown";
}

}  // namespace hipo

#endif