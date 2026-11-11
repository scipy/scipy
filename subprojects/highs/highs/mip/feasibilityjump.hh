/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"

#define FJ_LOG_PREFIX "Feasibility Jump: "

// TIP: clang-format the reference feasibilityjump.hh before diffing with this
namespace external_feasibilityjump {

enum RowType {
  Equal,
  Lte,
  Gte,
};

enum VarType { Continuous, Integer };

enum CallbackControlFlow {
  Terminate,
  Continue,
};

struct FJStatus {
  size_t totalEffort;
  size_t effortSinceLastImprovement;
  int numVars;
  double solutionObjectiveValue;
  double* solution;
};

// Measures if two doubles are equal within a tolerance
static bool eq(double a, double b, double tol) { return fabs(a - b) < tol; }

struct IdxCoeff {
  uint32_t idx;
  double coeff;

  IdxCoeff(uint32_t idx, double coeff) : idx(idx), coeff(coeff) {}
};

struct Var {
  VarType vartype;
  double lb;
  double ub;
  double objectiveCoeff;
  std::vector<IdxCoeff> coeffs;
};

struct Constraint {
  RowType sense;
  double rhs;
  std::vector<IdxCoeff> coeffs;
  double weight;
  double incumbentLhs;
  int32_t violatedIdx;

  // Computes the constraint's contribution to the feasibility score:
  // If the constraint is satisfied by the given LHS value, returns 0.
  // If the constraint is violated by the given LHS value, returns -|lhs-rhs|.
  double score(double lhs) {
    if (sense == RowType::Equal)
      return -fabs(lhs - rhs);
    else if (sense == RowType::Lte)
      return -(std::max(0., lhs - rhs));
    else
      return -(std::max(0., rhs - lhs));
  }
};

// A potential new value for a variable, including its score.
struct Move {
  double value;
  double score;

  static Move undef() {
    Move move;
    move.value = NAN;
    move.score = -std::numeric_limits<double>::infinity();
    return move;
  }
};

// Represents a modification of the LHS in a constraint, for a specific
// variable/constraint combination.The `modifyMove` function below is used to
// update the score of a `Move` to reflect the LHS modification.
struct LhsModification {
  uint32_t varIdx;
  uint32_t constraintIdx;
  double coeff;
  double oldLhs;
  double newLhs;
};

// Stores the MIP problem, an incumbent assignment, and the set of constraints
// that are violated in the current incumbent assignment. This set is maintained
// when changes are given to the incumbent assignment using `setValue`.
class Problem {
  const double equalityTolerance;
  const double violationTolerance;

 public:
  Problem(double equalityTolerance, double violationTolerance)
      : equalityTolerance(equalityTolerance),
        violationTolerance(violationTolerance) {}
  std::vector<Var> vars;
  std::vector<Constraint> constraints;
  std::vector<double> incumbentAssignment;
  std::vector<uint32_t> violatedConstraints;
  bool usedRelaxContinuous = false;

  size_t nNonzeros = 0;
  double incumbentObjective = NAN;

  int addVar(VarType vartype, double lb, double ub, double objCoeff) {
    auto idx = vars.size();
    Var var;
    var.vartype = vartype;
    var.lb = lb;
    var.ub = ub;
    var.objectiveCoeff = objCoeff;
    vars.push_back(var);
    incumbentAssignment.push_back(lb);
    return idx;
  }

  int addConstraint(RowType sense, double rhs, int numCoeffs, int* rowVarIdxs,
                    double* rowCoeffs, int relax_continuous) {
    if (relax_continuous) usedRelaxContinuous = true;

    // If we are relaxing continuous variables, an equality needs to be split
    // into Gte and Lte.
    if (relax_continuous > 0 && sense == RowType::Equal)
      if (std::any_of(rowVarIdxs, rowVarIdxs + numCoeffs, [&](int varIdx) {
            return vars[varIdx].vartype == VarType::Continuous;
          })) {
        addConstraint(RowType::Gte, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                      relax_continuous);
        addConstraint(RowType::Lte, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                      relax_continuous);
        return INT_MAX;
      }

    std::vector<IdxCoeff> coeffs;
    for (int i = 0; i < numCoeffs; i += 1) {
      if (relax_continuous > 0 &&
          vars[rowVarIdxs[i]].vartype == VarType::Continuous) {
        if (sense == RowType::Lte) {
          if (rowCoeffs[i] >= 0.)
            rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
          else
            rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
        } else if (sense == RowType::Gte) {
          if (rowCoeffs[i] >= 0.)
            rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
          else
            rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
        } else
          return INT_MIN;
      } else
        coeffs.emplace_back(rowVarIdxs[i], rowCoeffs[i]);
    }

    if (coeffs.empty()) {
      bool ok;
      if (sense == RowType::Lte)
        ok = 0 <= rhs + equalityTolerance;
      else if (sense == RowType::Gte)
        ok = 0 + equalityTolerance >= rhs;
      else
        ok = eq(0, rhs, equalityTolerance);

      return ok ? INT_MAX : INT_MIN;
    }

    int newConstraintIdx = constraints.size();
    for (auto& c : coeffs) {
      vars[c.idx].coeffs.emplace_back(newConstraintIdx, c.coeff);
    }

    nNonzeros += coeffs.size();
    Constraint newConstraint;
    newConstraint.coeffs = coeffs;
    newConstraint.incumbentLhs = NAN;
    newConstraint.violatedIdx = -1;
    newConstraint.rhs = rhs;
    newConstraint.sense = sense;
    newConstraint.weight = 1.0;

    constraints.push_back(newConstraint);
    return newConstraintIdx;
  }

  void resetIncumbent(double* initialValues) {
    // Set the initial values, if given.

    if (initialValues)
      for (size_t i = 0; i < vars.size(); i += 1)
        incumbentAssignment[i] = initialValues[i];
    // std::copy(initialValues, initialValues + vars.size(),
    // incumbentAssignment);

    // Reset the incumbent objective.
    incumbentObjective = 0;
    for (size_t i = 0; i < vars.size(); i += 1)
      incumbentObjective += vars[i].objectiveCoeff * incumbentAssignment[i];

    // Reset the constraint LHSs and the violatedConstraints list.
    violatedConstraints.clear();
    for (size_t cIdx = 0; cIdx < constraints.size(); cIdx += 1) {
      Constraint& cstr = constraints[cIdx];

      cstr.incumbentLhs = 0.0;
      for (auto& vc : cstr.coeffs)
        cstr.incumbentLhs += vc.coeff * incumbentAssignment[vc.idx];

      if (cstr.score(cstr.incumbentLhs) < -violationTolerance) {
        cstr.violatedIdx = violatedConstraints.size();
        violatedConstraints.push_back(cIdx);
      } else
        cstr.violatedIdx = -1;
    }
  }

  // Updates a variable assignment for `varIdx` to `newValue`.
  // Takes a function parameter f that receives a LhsModification
  // for every variable/constraint combination (except for `varIdx` itself)
  // where the LHS of the constraint has changed.
  template <typename F>
  size_t setValue(uint32_t varIdx, double newValue, F f) {
    size_t dt = 0;
    double oldValue = incumbentAssignment[varIdx];
    double delta = (newValue - oldValue);
    incumbentAssignment[varIdx] = newValue;
    incumbentObjective += vars[varIdx].objectiveCoeff * delta;
    // printf("Setting v%d to from %g to value %g\n", varIdx, oldValue,
    // newValue);

    // Update the LHSs of all involved constraints.
    for (auto& cstrCoeff : vars[varIdx].coeffs) {
      double oldLhs = constraints[cstrCoeff.idx].incumbentLhs;
      double newLhs = oldLhs + cstrCoeff.coeff * delta;
      constraints[cstrCoeff.idx].incumbentLhs = newLhs;
      double newCost = constraints[cstrCoeff.idx].score(newLhs);

      // Add/remove from the violatedConstraints list.
      if (newCost < -violationTolerance &&
          constraints[cstrCoeff.idx].violatedIdx == -1) {
        // Became violated.
        constraints[cstrCoeff.idx].violatedIdx = violatedConstraints.size();
        violatedConstraints.push_back(cstrCoeff.idx);
      }
      if (newCost >= -violationTolerance &&
          constraints[cstrCoeff.idx].violatedIdx != -1) {
        // Became satisfied.
        auto lastViolatedIdx = violatedConstraints.size() - 1;
        auto lastConstraintIdx = violatedConstraints[lastViolatedIdx];
        auto thisViolatedIdx = constraints[cstrCoeff.idx].violatedIdx;
        std::swap(violatedConstraints[thisViolatedIdx],
                  violatedConstraints[lastViolatedIdx]);
        constraints[lastConstraintIdx].violatedIdx = thisViolatedIdx;
        constraints[cstrCoeff.idx].violatedIdx = -1;
        violatedConstraints.pop_back();
      }

      // Now, report the changes in LHS for other variables.
      dt += constraints[cstrCoeff.idx].coeffs.size();
      for (auto& varCoeff : constraints[cstrCoeff.idx].coeffs) {
        if (varCoeff.idx != varIdx) {
          LhsModification m;
          m.varIdx = varCoeff.idx;
          m.constraintIdx = cstrCoeff.idx;
          m.coeff = varCoeff.coeff;
          m.oldLhs = oldLhs;
          m.newLhs = newLhs;
          f(m);
        }
      }
    }

    return dt;
  }
};

static void modifyMove(LhsModification mod, Problem& problem, Move& move) {
  Constraint& c = problem.constraints[mod.constraintIdx];
  auto incumbent = problem.incumbentAssignment[mod.varIdx];
  double oldModifiedLhs = mod.oldLhs + mod.coeff * (move.value - incumbent);
  double oldScoreTerm =
      c.weight * (c.score(oldModifiedLhs) - c.score(mod.oldLhs));
  double newModifiedLhs = mod.newLhs + mod.coeff * (move.value - incumbent);
  double newScoreTerm =
      c.weight * (c.score(newModifiedLhs) - c.score(mod.newLhs));
  move.score += newScoreTerm - oldScoreTerm;
}

// Stores current moves and computes updated jump values for
// the "Jump" move type.
class JumpMove {
  const double equalityTolerance;
  std::vector<Move> moves;
  std::vector<std::pair<double, double>> bestShiftBuffer;

 public:
  JumpMove(double equalityTolerance) : equalityTolerance(equalityTolerance) {}
  void init(Problem& problem) { moves.resize(problem.vars.size()); }

  template <typename F>
  void forEachVarMove(int32_t varIdx, F f) {
    f(moves[varIdx]);
  }

  void updateValue(Problem& problem, uint32_t varIdx) {
    bestShiftBuffer.clear();
    auto varIncumbentValue = problem.incumbentAssignment[varIdx];
    double lower = problem.vars[varIdx].lb;
    double upper = problem.vars[varIdx].ub;
    double currentValue = lower;  // can be -inf!

    double currentSlope = 0.0;

    // printf(" updatevalue lb %g ub %g numcells %d\n",
    //        lower,
    //        upper, problem.vars[varIdx].coeffs.size());

    for (auto& cell : problem.vars[varIdx].coeffs) {
      auto& constraint = problem.constraints[cell.idx];

      std::vector<std::pair<double, double>> constraintBounds;
      if (constraint.sense == RowType::Lte)
        constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(),
                                      constraint.rhs);
      else if (constraint.sense == RowType::Gte)
        constraintBounds.emplace_back(constraint.rhs,
                                      std::numeric_limits<double>::infinity());
      else {
        constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(),
                                      constraint.rhs);
        constraintBounds.emplace_back(constraint.rhs,
                                      std::numeric_limits<double>::infinity());
      }

      for (auto& bound : constraintBounds) {
        double residualIncumbent =
            constraint.incumbentLhs - cell.coeff * varIncumbentValue;

        double boundForValidLower, boundForValidUpper;
        if (cell.coeff < 0.) {
          boundForValidLower = bound.second;
          boundForValidUpper = bound.first;
        } else {
          boundForValidLower = bound.first;
          boundForValidUpper = bound.second;
        }
        std::pair<double, double> validRange = {
            ((1.0 / cell.coeff) * (boundForValidLower - residualIncumbent)),
            ((1.0 / cell.coeff) * (boundForValidUpper - residualIncumbent)),
        };

        if (problem.vars[varIdx].vartype == VarType::Integer)
          validRange = {
              std::ceil(validRange.first - equalityTolerance),
              std::floor(validRange.second + equalityTolerance),
          };

        if (validRange.first > currentValue) {
          currentSlope -= constraint.weight;
          if (validRange.first < upper)
            bestShiftBuffer.emplace_back(validRange.first, constraint.weight);
        }

        if (validRange.second <= currentValue) {
          currentSlope += constraint.weight;
        } else if (validRange.second < upper)
          bestShiftBuffer.emplace_back(validRange.second, constraint.weight);
      }
    }

    if (std::isfinite(lower)) {
      bestShiftBuffer.emplace_back(lower, 0);
    }
    if (std::isfinite(upper)) {
      bestShiftBuffer.emplace_back(upper, 0);
    }
    std::sort(bestShiftBuffer.begin(), bestShiftBuffer.end());

    // We only care about relative scores, so can pick an arbitrary start level
    double currentScore = 0.0;
    currentValue =
        bestShiftBuffer.empty() ? varIncumbentValue : bestShiftBuffer[0].first;

    double bestScore = currentScore;
    double bestValue = currentValue;

    // printf("evaluating best shift buffer size %d \n",
    // bestShiftBuffer.size());
    for (auto& item : bestShiftBuffer) {
      currentScore += (item.first - currentValue) * currentSlope;
      currentSlope += item.second;
      currentValue = item.first;

      // printf("bestshift cscore %g cslope %g cval %g bestval %g bestscore
      // %g\n", currentScore,currentSlope, currentValue, bestScore, bestValue
      // );

      if (eq(bestValue, varIncumbentValue, equalityTolerance) ||
          (!eq(currentValue, varIncumbentValue, equalityTolerance) &&
           currentScore < bestScore)) {
        // Different to everywhere else (!)
        // Here the score represents infeasibility (not improvement)
        // Hence the minimization
        bestScore = currentScore;
        bestValue = currentValue;
      }

      // Slope is always increasing, so if we have a valid value, we can quit
      // as soon as the slope turns nonnegative, since we must already have
      // visited the minimum.
      if (!eq(bestValue, varIncumbentValue, equalityTolerance) &&
          currentSlope >= 0.)
        break;
    }

    // printf("Setting jump for %d to from %g to %g\n", varIdx,
    // varIncumbentValue, moves[varIdx].value);
    moves[varIdx].value = bestValue;
  }
};

class FeasibilityJumpSolver {
  Problem problem;
  JumpMove jumpMove;

  std::vector<uint32_t> goodVarsSet;
  std::vector<int32_t> goodVarsSetIdx;

  std::mt19937 rng;

  double bestObjective = std::numeric_limits<double>::infinity();
  double objectiveWeight = 0.0;
  size_t bestViolationScore = SIZE_MAX;
  size_t effortAtLastCallback = 0;
  size_t effortAtLastImprovement = 0;
  size_t effortAtLastLogging = 0;
  size_t totalEffort = 0;

  double weightUpdateDecay;
  double weightUpdateIncrement = 1.0;

  size_t nBumps = 0;

  bool need_logging_header = true;

  // The probability of choosing a random positive-score variable.
  const double randomVarProbability = 0.001;

  // The probability of choosing a variable using a random constraint's
  // non-zero coefficient after updating weights.

  const double randomCellProbability = 0.01;

  // The number of moves to evaluate, if there are many positive-score
  // variables available.
  const size_t maxMovesToEvaluate = 25;

  const int kLoggingFrequency = 100000;  // Originally 100000
  const int kLoggingHeaderFrequency = 50;
  const size_t kMinEffortToLogging = 500 * kLoggingFrequency;
  const double kMinRelativeObjectiveImprovement = 1e-4;

  const HighsLogOptions& logOptions;

 public:
  FeasibilityJumpSolver(const HighsLogOptions& _logOptions, int seed = 0,
                        double equalityTolerance = 1e-5,
                        double violationTolerance = 1e-5,
                        double _weightUpdateDecay = 1.0)
      : problem(equalityTolerance, violationTolerance),
        jumpMove(equalityTolerance),
        logOptions(_logOptions) {
    weightUpdateDecay = _weightUpdateDecay;
    rng = std::mt19937(seed);
  }

  int addVar(VarType vartype, double lb, double ub, double objCoeff) {
    goodVarsSetIdx.push_back(-1);
    return problem.addVar(vartype, lb, ub, objCoeff);
  }

  int addConstraint(RowType sense, double rhs, int numCoeffs, int* rowVarIdxs,
                    double* rowCoeffs, int relax_continuous) {
    return problem.addConstraint(sense, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                                 relax_continuous);
  }

  int solve(double* initialValues,
            const std::function<CallbackControlFlow(FJStatus)>& callback) {
    assert(callback);
    highsLogDev(logOptions, HighsLogType::kInfo,
                FJ_LOG_PREFIX
                "starting solve. weightUpdateDecay=%g, relaxContinuous=%d  \n",
                weightUpdateDecay, problem.usedRelaxContinuous);

    init(initialValues);

    effortAtLastLogging = -kMinEffortToLogging;  // Enabling step=0 logging
    int num_logging_lines_since_header = 0;
    for (int step = 0; step < INT_MAX; step += 1) {
      if (user_terminate(callback, nullptr)) break;

      if (step % kLoggingFrequency == 0 &&
          totalEffort > effortAtLastLogging + kMinEffortToLogging) {
        if (need_logging_header) {
          logging(0, true);
          num_logging_lines_since_header = 0;
        }
        logging(step);
        num_logging_lines_since_header++;
        need_logging_header =
            num_logging_lines_since_header == kLoggingHeaderFrequency;
      }

      if (problem.violatedConstraints.size() < bestViolationScore) {
        effortAtLastImprovement = totalEffort;
        bestViolationScore = problem.violatedConstraints.size();
      }

      if (problem.violatedConstraints.empty() &&
          problem.incumbentObjective < bestObjective) {
        double relative_objective_improvement =
            (bestObjective - problem.incumbentObjective) /
            std::max(1.0, std::fabs(problem.incumbentObjective));
        if (relative_objective_improvement > kMinRelativeObjectiveImprovement) {
          effortAtLastImprovement = totalEffort;
          bestObjective = problem.incumbentObjective;
          if (user_terminate(callback, problem.incumbentAssignment.data()))
            break;
          // Repeat the header in case user callback logs new solution
          need_logging_header = true;
        }
      }

      if (problem.vars.size() == 0) break;

      uint32_t var = selectVariable();
      if (var == UINT_MAX) {
        break;
      }
      doVariableMove(var);
    }

    return 0;
  }

 private:
  void logging(const int step, const bool header = false) {
    const HighsLogType logType = HighsLogType::kDetailed;
    if (header) {
      highsLogDev(logOptions, logType,
                  FJ_LOG_PREFIX
                  "       step  violations     good    bumps       effort (per "
                  "step)          Objective\n");
    } else {
      highsLogDev(logOptions, logType,
                  " %10d    %8zd   %6zd %8zd %12zd    %6zd          %10.4g\n",
                  step, problem.violatedConstraints.size(), goodVarsSet.size(),
                  nBumps, totalEffort, step > 0 ? totalEffort / step : 0,
                  problem.incumbentObjective);
      effortAtLastLogging = totalEffort;
    }
  }

  void init(double* initialValues) {
    problem.resetIncumbent(initialValues);
    jumpMove.init(problem);
    totalEffort += problem.nNonzeros;

    // Reset the variable scores.
    goodVarsSet.clear();
    for (size_t i = 0; i < problem.vars.size(); i += 1) resetMoves(i);
  }

  uint32_t selectVariable() {
    if (!goodVarsSet.empty()) {
      if (std::uniform_real_distribution<double>(0., 1.)(rng) <
          randomVarProbability)
        return goodVarsSet[rng() % goodVarsSet.size()];

      auto sampleSize = std::min(maxMovesToEvaluate, goodVarsSet.size());
      totalEffort += sampleSize;

      double bestScore = -std::numeric_limits<double>::infinity();
      uint32_t bestVar = UINT_MAX;
      for (size_t i = 0; i < sampleSize; i++) {
        auto setidx = rng() % goodVarsSet.size();
        auto varIdx = goodVarsSet[setidx];
        Move move = bestMove(varIdx);
        if (move.score > bestScore) {
          bestScore = move.score;
          bestVar = varIdx;
        }
      }
      return bestVar;
    }

    // Local minimum, update weights.
    updateWeights();

    if (!problem.violatedConstraints.empty()) {
      size_t cstrIdx =
          problem
              .violatedConstraints[rng() % problem.violatedConstraints.size()];
      auto& constraint = problem.constraints[cstrIdx];

      if (std::uniform_real_distribution<double>(0., 1.)(rng) <
          randomCellProbability)
        return constraint.coeffs[rng() % constraint.coeffs.size()].idx;

      double bestScore = -std::numeric_limits<double>::infinity();
      uint32_t bestVarIdx = UINT_MAX;

      for (auto& cell : constraint.coeffs) {
        Move move = bestMove(cell.idx);
        if (move.score > bestScore) {
          bestScore = move.score;
          bestVarIdx = cell.idx;
        }
      }
      return bestVarIdx;
    }

    // Fallback to random choice.
    return rng() % problem.vars.size();
  }

  void updateWeights() {
    highsLogDev(logOptions, HighsLogType::kVerbose,
                FJ_LOG_PREFIX "Reached a local minimum.\n");
    nBumps += 1;
    bool rescaleAllWeights = false;
    size_t dt = 0;

    if (problem.violatedConstraints.empty()) {
      objectiveWeight += weightUpdateIncrement;
      if (objectiveWeight > 1.0e20) rescaleAllWeights = true;

      dt += problem.vars.size();
      for (size_t varIdx = 0; varIdx < problem.vars.size(); varIdx += 1)
        forEachMove(varIdx, [&](Move& move) {
          // -= to align objective minimization with score maximization
          move.score -= weightUpdateIncrement *
                        problem.vars[varIdx].objectiveCoeff *
                        (move.value - problem.incumbentAssignment[varIdx]);
        });
    } else {
      for (auto& cIdx : problem.violatedConstraints) {
        auto& constraint = problem.constraints[cIdx];
        constraint.weight += weightUpdateIncrement;
        if (constraint.weight > 1.0e20) rescaleAllWeights = true;

        dt += constraint.coeffs.size();
        for (auto& cell : constraint.coeffs) {
          forEachMove(cell.idx, [&](Move& move) {
            double candidateLhs =
                constraint.incumbentLhs +
                cell.coeff *
                    (move.value - problem.incumbentAssignment[cell.idx]);
            double diff = weightUpdateIncrement *
                          (constraint.score(candidateLhs) -
                           constraint.score(constraint.incumbentLhs));
            move.score += diff;
          });

          updateGoodMoves(cell.idx);
        }
      }
    }

    weightUpdateIncrement /= weightUpdateDecay;
    if (rescaleAllWeights) {
      weightUpdateIncrement *= 1.0e-20;
      objectiveWeight *= 1.0e-20;

      for (auto& c : problem.constraints) c.weight *= 1.0e-20;
      dt += problem.constraints.size();

      for (size_t i = 0; i < problem.vars.size(); i += 1) resetMoves(i);
    }

    totalEffort += dt;
  }

  Move bestMove(uint32_t varIdx) {
    Move best = Move::undef();
    forEachMove(varIdx, [&](Move& move) {
      if (move.score > best.score) best = move;
    });

    return best;
  }

  void doVariableMove(uint32_t varIdx) {
    // First, we get the best move for the variable;
    auto m = bestMove(varIdx);
    auto newValue = m.value;
    // assert(!isnan(newValue));

    // Update the incumbent solution.
    // printf("Setting var %d from %g to %g for a score of %g\n", varIdx,
    // oldValue, newValue, m.score);

    totalEffort += problem.setValue(varIdx, newValue, [&](LhsModification mod) {
      forEachMove(mod.varIdx, [&](Move& m) { modifyMove(mod, problem, m); });
      updateGoodMoves(mod.varIdx);
    });

    resetMoves(varIdx);
  }

  void updateGoodMoves(int32_t varIdx) {
    bool anyGoodMoves = bestMove(varIdx).score > 0.;
    if (anyGoodMoves && goodVarsSetIdx[varIdx] == -1) {
      // Became good, add to good set.
      goodVarsSetIdx[varIdx] = goodVarsSet.size();
      goodVarsSet.push_back(varIdx);
    } else if (!anyGoodMoves && goodVarsSetIdx[varIdx] != -1) {
      // Became bad, remove from good set.
      auto lastSetIdx = goodVarsSet.size() - 1;
      auto lastVarIdx = goodVarsSet[lastSetIdx];
      auto thisSetIdx = goodVarsSetIdx[varIdx];
      std::swap(goodVarsSet[thisSetIdx], goodVarsSet[lastSetIdx]);
      goodVarsSetIdx[lastVarIdx] = thisSetIdx;
      goodVarsSetIdx[varIdx] = -1;
      goodVarsSet.pop_back();
    }
  }

  template <typename F>
  void forEachMove(int32_t varIdx, F f) {
    jumpMove.forEachVarMove(varIdx, f);

    // TODO: here, we can add more move types.
    // upDownMove.forEachVarMove(varIdx, f);
  }

  void resetMoves(uint32_t varIdx) {
    totalEffort += problem.vars[varIdx].coeffs.size();
    jumpMove.updateValue(problem, varIdx);

    forEachMove(varIdx, [&](Move& move) {
      move.score = 0.0;
      // -= to align objective minimization with score maximization
      move.score -= objectiveWeight * problem.vars[varIdx].objectiveCoeff *
                    (move.value - problem.incumbentAssignment[varIdx]);

      for (auto& cell : problem.vars[varIdx].coeffs) {
        auto& constraint = problem.constraints[cell.idx];
        auto candidateLhs =
            constraint.incumbentLhs +
            cell.coeff * (move.value - problem.incumbentAssignment[varIdx]);
        move.score +=
            constraint.weight * (constraint.score(candidateLhs) -
                                 constraint.score(constraint.incumbentLhs));
      }
    });

    updateGoodMoves(varIdx);
  }

  bool user_terminate(
      const std::function<CallbackControlFlow(FJStatus)>& callback,
      double* solution) {
    const int CALLBACK_EFFORT = 500000;  // Originally 500000
    if (solution != nullptr ||
        totalEffort - effortAtLastCallback > CALLBACK_EFFORT) {
      highsLogDev(logOptions, HighsLogType::kVerbose,
                  FJ_LOG_PREFIX "calling user termination.\n");
      effortAtLastCallback = totalEffort;

      FJStatus status;
      status.totalEffort = totalEffort;
      status.effortSinceLastImprovement = totalEffort - effortAtLastImprovement;

      status.solution = solution;
      status.numVars = problem.vars.size();
      status.solutionObjectiveValue = problem.incumbentObjective;

      auto result = callback(status);
      if (result == CallbackControlFlow::Terminate) {
        highsLogDev(logOptions, HighsLogType::kVerbose,
                    FJ_LOG_PREFIX "quitting.\n");
        return true;
      }
    }
    return false;
  }
};

}  // namespace external_feasibilityjump
