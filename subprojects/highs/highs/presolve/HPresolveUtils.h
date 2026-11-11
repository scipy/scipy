/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_HIGHS_PRESOLVE_UTILS_H_
#define PRESOLVE_HIGHS_PRESOLVE_UTILS_H_

#include "lp_data/HConst.h"

#define HPRESOLVE_CHECKED_CALL(presolveCall)            \
  do {                                                  \
    HPresolve::Result __result = presolveCall;          \
    if (__result != presolve::HPresolve::Result::kOk) { \
      return __result;                                  \
    }                                                   \
  } while (0)

// As HPRESOLVE_CHECKED_CALL, timing the call with the given presolve clock if
// active is true. The timer is stopped on every exit, including the early
// return
#define HPRESOLVE_TIMED_CHECKED_CALL_IF(active, presolveClock, presolveCall)  \
  do {                                                                        \
    PresolveTimerScope presolveTimerScope_(analysis_, presolveClock, active); \
    HPRESOLVE_CHECKED_CALL(presolveCall);                                     \
  } while (0)

// As HPRESOLVE_CHECKED_CALL, timing the call with the given presolve clock
#define HPRESOLVE_TIMED_CHECKED_CALL(presolveClock, presolveCall) \
  HPRESOLVE_TIMED_CHECKED_CALL_IF(true, presolveClock, presolveCall)

class HighsCliqueTable;

namespace presolve {

// Sets the clique table's presolve probing flag and resets it when the scope
// is left, on every exit
class PresolveProbingScope {
 public:
  explicit PresolveProbingScope(HighsCliqueTable& cliquetable);
  ~PresolveProbingScope();
  PresolveProbingScope(const PresolveProbingScope&) = delete;
  PresolveProbingScope& operator=(const PresolveProbingScope&) = delete;

 private:
  HighsCliqueTable& cliquetable_;
};

enum class SingletonRowResult {
  kRedundant,
  kPrimalInfeasible,
  kBoundsTightened,
};

SingletonRowResult computeSingletonRowBounds(
    double val, double rowLower, double rowUpper, double colLower,
    double colUpper, double primalFeastol, double maxAbsColVal, bool isIntegral,
    double& lb, double& ub, bool& lowerTightened, bool& upperTightened);

}  // namespace presolve
#endif
