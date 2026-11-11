/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "util/HighsLinearSumBounds.h"

#include <algorithm>

void HighsLinearSumBounds::add(HighsInt sum, HighsInt var, double coefficient) {
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{1});
  handleImplVarUpper(sum, coefficient, getImplVarUpper(sum, var), HighsInt{1});
  handleImplVarLower(sum, coefficient, getImplVarLower(sum, var), HighsInt{1});
}

void HighsLinearSumBounds::remove(HighsInt sum, HighsInt var,
                                  double coefficient) {
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{-1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{-1});
  handleImplVarUpper(sum, coefficient, getImplVarUpper(sum, var), HighsInt{-1});
  handleImplVarLower(sum, coefficient, getImplVarLower(sum, var), HighsInt{-1});
}

void HighsLinearSumBounds::updatedVarUpper(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarUpper) {
  handleVarUpper(sum, coefficient, oldVarUpper, HighsInt{-1});
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{1});
  updatedImplVarUpper(sum, var, coefficient, oldVarUpper, implVarUpper[var],
                      implVarUpperSource[var]);
}

void HighsLinearSumBounds::updatedVarLower(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarLower) {
  handleVarLower(sum, coefficient, oldVarLower, HighsInt{-1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{1});
  updatedImplVarLower(sum, var, coefficient, oldVarLower, implVarLower[var],
                      implVarLowerSource[var]);
}

void HighsLinearSumBounds::updatedImplVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarUpper,
                                               HighsInt oldImplVarUpperSource) {
  updatedImplVarUpper(sum, var, coefficient, varUpper[var], oldImplVarUpper,
                      oldImplVarUpperSource);
}

void HighsLinearSumBounds::updatedImplVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarLower,
                                               HighsInt oldImplVarLowerSource) {
  updatedImplVarLower(sum, var, coefficient, varLower[var], oldImplVarLower,
                      oldImplVarLowerSource);
}

double HighsLinearSumBounds::getResidualSumLower(HighsInt sum, HighsInt var,
                                                 double coefficient,
                                                 HighsInt boundVar,
                                                 double boundVarCoefficient,
                                                 double boundVarValue) const {
  HighsCDouble activity = sumLower[sum];
  HighsInt numInfs = numInfSumLower[sum];
  residual(
      numInfs, activity,
      (coefficient > 0 ? getImplVarLower(sum, var) : getImplVarUpper(sum, var)),
      coefficient, boundVar, boundVarCoefficient,
      (boundVar != -1
           ? (boundVarCoefficient > 0 ? getImplVarLower(sum, boundVar)
                                      : getImplVarUpper(sum, boundVar))
           : kHighsInf),
      boundVarValue);
  return (numInfs == 0 ? static_cast<double>(activity) : -kHighsInf);
}

double HighsLinearSumBounds::getResidualSumUpper(HighsInt sum, HighsInt var,
                                                 double coefficient,
                                                 HighsInt boundVar,
                                                 double boundVarCoefficient,
                                                 double boundVarValue) const {
  HighsCDouble activity = sumUpper[sum];
  HighsInt numInfs = numInfSumUpper[sum];
  residual(
      numInfs, activity,
      (coefficient < 0 ? getImplVarLower(sum, var) : getImplVarUpper(sum, var)),
      coefficient, boundVar, boundVarCoefficient,
      (boundVar != -1
           ? (boundVarCoefficient < 0 ? getImplVarLower(sum, boundVar)
                                      : getImplVarUpper(sum, boundVar))
           : kHighsInf),
      boundVarValue);
  return (numInfs == 0 ? static_cast<double>(activity) : kHighsInf);
}

double HighsLinearSumBounds::getResidualSumLowerOrig(
    HighsInt sum, HighsInt var, double coefficient, HighsInt boundVar,
    double boundVarCoefficient, double boundVarValue) const {
  HighsCDouble activity = sumLowerOrig[sum];
  HighsInt numInfs = numInfSumLowerOrig[sum];
  residual(numInfs, activity, (coefficient > 0 ? varLower[var] : varUpper[var]),
           coefficient, boundVar, boundVarCoefficient,
           (boundVar != -1 ? (boundVarCoefficient > 0 ? varLower[boundVar]
                                                      : varUpper[boundVar])
                           : kHighsInf),
           boundVarValue);
  return (numInfs == 0 ? static_cast<double>(activity) : -kHighsInf);
}

double HighsLinearSumBounds::getResidualSumUpperOrig(
    HighsInt sum, HighsInt var, double coefficient, HighsInt boundVar,
    double boundVarCoefficient, double boundVarValue) const {
  HighsCDouble activity = sumUpperOrig[sum];
  HighsInt numInfs = numInfSumUpperOrig[sum];
  residual(numInfs, activity, (coefficient < 0 ? varLower[var] : varUpper[var]),
           coefficient, boundVar, boundVarCoefficient,
           (boundVar != -1 ? (boundVarCoefficient < 0 ? varLower[boundVar]
                                                      : varUpper[boundVar])
                           : kHighsInf),
           boundVarValue);
  return (numInfs == 0 ? static_cast<double>(activity) : kHighsInf);
}

void HighsLinearSumBounds::shrink(const std::vector<HighsInt>& newIndices,
                                  HighsInt newSize) {
  for (size_t i = 0; i != newIndices.size(); ++i) {
    if (newIndices[i] != -1) {
      sumLower[newIndices[i]] = sumLower[i];
      sumUpper[newIndices[i]] = sumUpper[i];
      numInfSumLower[newIndices[i]] = numInfSumLower[i];
      numInfSumUpper[newIndices[i]] = numInfSumUpper[i];
      sumLowerOrig[newIndices[i]] = sumLowerOrig[i];
      sumUpperOrig[newIndices[i]] = sumUpperOrig[i];
      numInfSumLowerOrig[newIndices[i]] = numInfSumLowerOrig[i];
      numInfSumUpperOrig[newIndices[i]] = numInfSumUpperOrig[i];
    }
  }

  sumLower.resize(newSize);
  sumUpper.resize(newSize);
  numInfSumLower.resize(newSize);
  numInfSumUpper.resize(newSize);
  sumLowerOrig.resize(newSize);
  sumUpperOrig.resize(newSize);
  numInfSumLowerOrig.resize(newSize);
  numInfSumUpperOrig.resize(newSize);
}

double HighsLinearSumBounds::getImplVarUpper(HighsInt sum, HighsInt var) const {
  return getImplVarUpper(sum, varUpper[var], implVarUpper[var],
                         implVarUpperSource[var]);
}

double HighsLinearSumBounds::getImplVarLower(HighsInt sum, HighsInt var) const {
  return getImplVarLower(sum, varLower[var], implVarLower[var],
                         implVarLowerSource[var]);
}

double HighsLinearSumBounds::getImplVarUpper(
    HighsInt sum, double myVarUpper, double myImplVarUpper,
    HighsInt myImplVarUpperSource) const {
  return (myImplVarUpperSource == sum ? myVarUpper
                                      : std::min(myImplVarUpper, myVarUpper));
}

double HighsLinearSumBounds::getImplVarLower(
    HighsInt sum, double myVarLower, double myImplVarLower,
    HighsInt myImplVarLowerSource) const {
  return (myImplVarLowerSource == sum ? myVarLower
                                      : std::max(myImplVarLower, myVarLower));
}

bool HighsLinearSumBounds::implBndsValidAfterCoefChange(HighsInt sum,
                                                        HighsInt var,
                                                        double oldCoef,
                                                        double newCoef) {
  HighsInt oldNumInfLower = numInfSumLower[sum];
  HighsCDouble oldSumLower = sumLower[sum];
  HighsInt oldNumInfUpper = numInfSumUpper[sum];
  HighsCDouble oldSumUpper = sumUpper[sum];

  if (oldCoef != 0.0) remove(sum, var, oldCoef);
  if (newCoef != 0.0) add(sum, var, newCoef);

  return numInfSumLower[sum] <= oldNumInfLower &&
         sumLower[sum] >= oldSumLower &&
         numInfSumUpper[sum] <= oldNumInfUpper && sumUpper[sum] <= oldSumUpper;
}

void HighsLinearSumBounds::update(HighsInt& numInfs, HighsCDouble& activity,
                                  HighsInt direction, double bound,
                                  double coefficient) const {
  if (std::abs(bound) != kHighsInf)
    activity += direction * static_cast<HighsCDouble>(bound) * coefficient;
  else
    numInfs += direction;
}

void HighsLinearSumBounds::update(HighsInt& numInfs, HighsCDouble& activity,
                                  double oldBound, double newBound,
                                  double coefficient) const {
  if (oldBound == newBound) return;
  update(numInfs, activity, HighsInt{-1}, oldBound, coefficient);
  update(numInfs, activity, HighsInt{1}, newBound, coefficient);
}

void HighsLinearSumBounds::residual(HighsInt& numInfs, HighsCDouble& activity,
                                    double varBound, double coefficient,
                                    HighsInt boundVar,
                                    double boundVarCoefficient,
                                    double oldBoundVarBound,
                                    double newBoundVarBound) const {
  // remove contribution of 'var'
  update(numInfs, activity, HighsInt{-1}, varBound, coefficient);
  // update contribution of 'boundVar'
  if (boundVar != -1)
    update(numInfs, activity, oldBoundVarBound, newBoundVarBound,
           boundVarCoefficient);
}

void HighsLinearSumBounds::handleVarUpper(HighsInt sum, double coefficient,
                                          double myVarUpper,
                                          HighsInt direction) {
  update(coefficient > 0 ? numInfSumUpperOrig[sum] : numInfSumLowerOrig[sum],
         coefficient > 0 ? sumUpperOrig[sum] : sumLowerOrig[sum], direction,
         myVarUpper, coefficient);
}

void HighsLinearSumBounds::handleVarLower(HighsInt sum, double coefficient,
                                          double myVarLower,
                                          HighsInt direction) {
  update(coefficient > 0 ? numInfSumLowerOrig[sum] : numInfSumUpperOrig[sum],
         coefficient > 0 ? sumLowerOrig[sum] : sumUpperOrig[sum], direction,
         myVarLower, coefficient);
}

void HighsLinearSumBounds::handleImplVarUpper(HighsInt sum, double coefficient,
                                              double myImplVarUpper,
                                              HighsInt direction) {
  update(coefficient > 0 ? numInfSumUpper[sum] : numInfSumLower[sum],
         coefficient > 0 ? sumUpper[sum] : sumLower[sum], direction,
         myImplVarUpper, coefficient);
}

void HighsLinearSumBounds::handleImplVarLower(HighsInt sum, double coefficient,
                                              double myImplVarLower,
                                              HighsInt direction) {
  update(coefficient > 0 ? numInfSumLower[sum] : numInfSumUpper[sum],
         coefficient > 0 ? sumLower[sum] : sumUpper[sum], direction,
         myImplVarLower, coefficient);
}

void HighsLinearSumBounds::updatedImplVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldVarUpper,
                                               double oldImplVarUpper,
                                               HighsInt oldImplVarUpperSource) {
  double oldVUpper =
      getImplVarUpper(sum, oldVarUpper, oldImplVarUpper, oldImplVarUpperSource);
  double vUpper = getImplVarUpper(sum, var);

  if (vUpper == oldVUpper) return;

  handleImplVarUpper(sum, coefficient, oldVUpper, HighsInt{-1});
  handleImplVarUpper(sum, coefficient, vUpper, HighsInt{1});
}

void HighsLinearSumBounds::updatedImplVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldVarLower,
                                               double oldImplVarLower,
                                               HighsInt oldImplVarLowerSource) {
  double oldVLower =
      getImplVarLower(sum, oldVarLower, oldImplVarLower, oldImplVarLowerSource);
  double vLower = getImplVarLower(sum, var);

  if (vLower == oldVLower) return;

  handleImplVarLower(sum, coefficient, oldVLower, HighsInt{-1});
  handleImplVarLower(sum, coefficient, vLower, HighsInt{1});
}
