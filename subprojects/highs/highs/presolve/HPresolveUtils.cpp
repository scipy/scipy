/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolveUtils.h"

#include "mip/HighsCliqueTable.h"

namespace presolve {

PresolveProbingScope::PresolveProbingScope(HighsCliqueTable& cliquetable)
    : cliquetable_(cliquetable) {
  cliquetable_.setinPresolveProbingFlag(true);
}

PresolveProbingScope::~PresolveProbingScope() {
  cliquetable_.setinPresolveProbingFlag(false);
}

SingletonRowResult computeSingletonRowBounds(
    double val, double rowLower, double rowUpper, double colLower,
    double colUpper, double primalFeastol, double maxAbsColVal, bool isIntegral,
    double& lb, double& ub, bool& lowerTightened, bool& upperTightened) {
  if (val > 0) {
    if (colUpper * val <= rowUpper + primalFeastol &&
        colLower * val >= rowLower - primalFeastol)
      return SingletonRowResult::kRedundant;
  } else {
    if (colLower * val <= rowUpper + primalFeastol &&
        colUpper * val >= rowLower - primalFeastol)
      return SingletonRowResult::kRedundant;
  }

  double newColUpper = kHighsInf;
  double newColLower = -kHighsInf;
  if (val > 0) {
    if (rowUpper != kHighsInf) newColUpper = rowUpper / val;
    if (rowLower != -kHighsInf) newColLower = rowLower / val;
  } else {
    if (rowUpper != kHighsInf) newColLower = rowUpper / val;
    if (rowLower != -kHighsInf) newColUpper = rowLower / val;
  }

  // use either the primal feasibility tolerance for the bound constraint or
  // for the singleton row including scaling, whichever is tighter.
  const double boundTol = std::max(primalFeastol / std::max(1.0, std::abs(val)),
                                   std::numeric_limits<double>::epsilon());

  lowerTightened = newColLower > colLower + boundTol;
  upperTightened = newColUpper < colUpper - boundTol;

  if (lowerTightened) {
    if (isIntegral) newColLower = std::ceil(newColLower - boundTol);
    lb = newColLower;
  } else
    lb = colLower;

  if (upperTightened) {
    if (isIntegral) newColUpper = std::floor(newColUpper + boundTol);
    ub = newColUpper;
  } else
    ub = colUpper;

  // check whether the bounds are equal in tolerances
  if (ub <= lb + primalFeastol) {
    // bounds could be infeasible or equal in tolerances, first check infeasible
    if (ub < lb - primalFeastol) return SingletonRowResult::kPrimalInfeasible;

    // bounds are equal in tolerances, if they have a slight infeasibility below
    // those tolerances or they have a slight numerical distance which changes
    // the largest contribution below feasibility tolerance then we can safely
    // set the bound to one of the values. To heuristically get rid of numerical
    // errors we choose the bound that was not tightened, or the midpoint if
    // both where tightened.
    if (ub < lb ||
        (ub > lb &&
         (ub - lb) * std::max(std::abs(val), maxAbsColVal) <= primalFeastol)) {
      if (lowerTightened && upperTightened) {
        ub = 0.5 * (ub + lb);
        lb = ub;
        lowerTightened = lb > colLower;
        upperTightened = ub < colUpper;
      } else if (lowerTightened) {
        lb = ub;
        lowerTightened = lb > colLower;
      } else {
        ub = lb;
        upperTightened = ub < colUpper;
      }
    }
  }

  return SingletonRowResult::kBoundsTightened;
}

}  // namespace presolve
