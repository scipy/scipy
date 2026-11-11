/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsCutGeneration.h"

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "util/HighsIntegers.h"

HighsCutGeneration::HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                                       HighsCutPool& cutpool)
    : lpRelaxation(lpRelaxation),
      cutpool(cutpool),
      randgen(lpRelaxation.getMipSolver().options_mip_->random_seed +
              lpRelaxation.getNumLpIterations() + cutpool.getNumCuts()),
      feastol(lpRelaxation.getMipSolver().mipdata_->feastol),
      epsilon(lpRelaxation.getMipSolver().mipdata_->epsilon) {}

bool HighsCutGeneration::determineCover(bool lpSol) {
  if (rhs <= 10 * feastol) return false;

  cover.clear();
  cover.reserve(rowlen);

  for (HighsInt j = 0; j != rowlen; ++j) {
    if (!isintegral[j]) continue;

    if (lpSol && solval[j] <= feastol) continue;

    cover.push_back(j);
  }

  HighsInt maxCoverSize = cover.size();
  HighsInt coversize = 0;
  HighsInt r = randgen.integer();
  coverweight = 0.0;
  if (lpSol) {
    // take all variables that sit at their upper bound always into the cover
    coversize = std::partition(cover.begin(), cover.end(),
                               [&](HighsInt j) {
                                 return solval[j] >= upper[j] - feastol;
                               }) -
                cover.begin();

    for (HighsInt i = 0; i != coversize; ++i) {
      HighsInt j = cover[i];

      assert(solval[j] >= upper[j] - feastol);

      coverweight += vals[j] * upper[j];
    }

    // sort the remaining variables by the contribution to the rows activity in
    // the current solution
    pdqsort(cover.begin() + coversize, cover.begin() + maxCoverSize,
            [&](HighsInt i, HighsInt j) {
              if (upper[i] < 1.5 && upper[j] > 1.5) return true;
              if (upper[i] > 1.5 && upper[j] < 1.5) return false;

              double contributionA = solval[i] * vals[i];
              double contributionB = solval[j] * vals[j];

              if (contributionA > contributionB + feastol) return true;
              if (contributionA < contributionB - feastol) return false;
              // for equal contributions take the larger coefficients first
              // because this makes some of the lifting functions more likely
              // to generate a facet
              // if the value is equal too, choose a random tiebreaker based
              // on hashing the column index and the current number of pool
              // cuts
              if (std::abs(vals[i] - vals[j]) <= feastol)
                return HighsHashHelpers::hash(std::make_pair(inds[i], r)) >
                       HighsHashHelpers::hash(std::make_pair(inds[j], r));
              return vals[i] > vals[j];
            });
  } else {
    // the current solution
    const auto& nodequeue = lpRelaxation.getMipSolver().mipdata_->nodequeue;

    pdqsort(cover.begin() + coversize, cover.begin() + maxCoverSize,
            [&](HighsInt i, HighsInt j) {
              if (solval[i] > feastol && solval[j] <= feastol) return true;
              if (solval[i] <= feastol && solval[j] > feastol) return false;

              int64_t numNodesA;
              int64_t numNodesB;

              numNodesA = complementation[i] ? nodequeue.numNodesDown(inds[i])
                                             : nodequeue.numNodesUp(inds[i]);

              numNodesB = complementation[j] ? nodequeue.numNodesDown(inds[j])
                                             : nodequeue.numNodesUp(inds[j]);

              if (numNodesA > numNodesB) return true;
              if (numNodesA < numNodesB) return false;

              return HighsHashHelpers::hash(std::make_pair(inds[i], r)) >
                     HighsHashHelpers::hash(std::make_pair(inds[j], r));
            });
  }

  const double minlambda =
      std::max(10 * feastol, feastol * std::abs(double(rhs)));

  for (; coversize != maxCoverSize; ++coversize) {
    double lambda = double(coverweight - rhs);
    if (lambda > minlambda) break;

    HighsInt j = cover[coversize];
    coverweight += vals[j] * upper[j];
  }
  if (coversize == 0) return false;

  coverweight.renormalize();
  lambda = coverweight - rhs;

  if (lambda <= minlambda) return false;

  cover.resize(coversize);
  assert(lambda > feastol);

  return true;
}

void HighsCutGeneration::separateLiftedKnapsackCover() {
  const HighsInt coversize = cover.size();

  std::vector<double> S;
  S.resize(coversize);
  std::vector<int8_t> coverflag;
  coverflag.resize(rowlen);
  pdqsort_branchless(cover.begin(), cover.end(),
                     [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; });

  HighsCDouble abartmp = vals[cover[0]];
  HighsCDouble sigma = lambda;
  for (HighsInt i = 1; i != coversize; ++i) {
    HighsCDouble delta = abartmp - vals[cover[i]];
    HighsCDouble kdelta = double(i) * delta;
    if (double(kdelta) < double(sigma)) {
      abartmp = vals[cover[i]];
      sigma -= kdelta;
    } else {
      abartmp -= sigma * (1.0 / i);
      sigma = 0.0;
      break;
    }
  }

  if (double(sigma) > 0) abartmp = HighsCDouble(rhs) / double(coversize);

  double abar = double(abartmp);

  HighsCDouble sum = 0.0;
  HighsInt cplussize = 0;
  for (HighsInt i = 0; i != coversize; ++i) {
    sum += std::min(abar, vals[cover[i]]);
    S[i] = double(sum);

    if (vals[cover[i]] > abar + feastol) {
      ++cplussize;
      coverflag[cover[i]] = 1;
    } else
      coverflag[cover[i]] = -1;
  }
  assert(std::abs(double(sum - rhs) / double(rhs)) <= 1e-14);
  bool halfintegral = false;

  /* define the lifting function */
  auto g = [&](double z) {
    double hfrac = z / abar;
    double coef = 0.0;

    HighsInt h = std::floor(hfrac + 0.5);
    if (h != 0 && std::abs(hfrac - h) * std::max(1.0, abar) <= epsilon &&
        h <= cplussize - 1) {
      halfintegral = true;
      coef = 0.5;
    }

    h = std::max(h - 1, HighsInt{0});
    for (; h < coversize; ++h) {
      if (z <= S[h] + feastol) break;
    }

    return coef + h;
  };

  rhs = coversize - 1;

  for (HighsInt i = 0; i != rowlen; ++i) {
    if (vals[i] == 0.0) continue;
    if (coverflag[i] == -1) {
      vals[i] = 1;
    } else {
      vals[i] = g(vals[i]);
    }
  }

  if (halfintegral) {
    rhs *= 2;
    for (HighsInt i = 0; i != rowlen; ++i) vals[i] *= 2;
  }

  // resulting cut is always integral
  integralSupport = true;
  integralCoefficients = true;
}

bool HighsCutGeneration::separateLiftedMixedBinaryCover() {
  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  HighsInt coversize = cover.size();
  std::vector<double> S;
  S.resize(coversize);
  std::vector<HighsBool> coverflag;
  coverflag.resize(rowlen);

  if (coversize == 0) return false;

  for (HighsInt i = 0; i != coversize; ++i) coverflag[cover[i]] = true;

  pdqsort_branchless(cover.begin(), cover.end(),
                     [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; });
  HighsCDouble sum = 0;

  HighsInt p = coversize;
  for (HighsInt i = 0; i != coversize; ++i) {
    if (vals[cover[i]] - lambda <= epsilon) {
      p = i;
      break;
    }
    sum += vals[cover[i]];
    S[i] = double(sum);
  }
  if (p == 0) return false;
  /* define the lifting function */
  auto phi = [&](double a) {
    for (HighsInt i = 0; i < p; ++i) {
      if (a <= S[i] - lambda) return double(i * lambda);

      if (a <= S[i]) return double((i + 1) * lambda + (HighsCDouble(a) - S[i]));
    }

    return double(p * lambda + (HighsCDouble(a) - S[p - 1]));
  };

  rhs = -lambda;

  integralCoefficients = false;
  integralSupport = true;
  for (HighsInt i = 0; i != rowlen; ++i) {
    if (!isintegral[i]) {
      if (vals[i] < 0)
        integralSupport = false;
      else
        vals[i] = 0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = std::min(vals[i], double(lambda));
      rhs += vals[i];
    } else {
      vals[i] = phi(vals[i]);
    }
  }

  return true;
}

bool HighsCutGeneration::separateLiftedMixedIntegerCover() {
  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  HighsInt coversize = cover.size();

  HighsInt l = -1;

  std::vector<HighsBool> coverflag;
  coverflag.resize(rowlen);
  for (HighsInt i : cover) coverflag[i] = true;

  auto comp = [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; };
  pdqsort_branchless(cover.begin(), cover.end(), comp);

  std::vector<HighsCDouble> a;
  std::vector<HighsCDouble> u;
  std::vector<HighsCDouble> m;

  a.resize(coversize);
  u.resize(coversize + 1);
  m.resize(coversize + 1);

  HighsCDouble usum = 0.0;
  HighsCDouble msum = 0.0;
  // set up the partial sums of the upper bounds, and the contributions
  for (HighsInt c = 0; c != coversize; ++c) {
    HighsInt i = cover[c];

    u[c] = usum;
    m[c] = msum;
    a[c] = vals[i];
    double ub = upper[i];
    usum += ub;
    msum += ub * a[c];
  }

  u[coversize] = usum;
  m[coversize] = msum;

  // determine which variable in the cover we want to create the MIR inequality
  // from which we lift we try to select a variable to have the highest chance
  // of satisfying the facet conditions for the superadditive lifting function
  // gamma to be satisfied.
  HighsInt lpos = -1;
  HighsInt bestlCplusend = -1;
  double bestlVal = 0.0;
  bool bestlAtUpper = true;

  for (HighsInt i = 0; i != coversize; ++i) {
    HighsInt j = cover[i];
    double ub = upper[j];

    bool atUpper = solval[j] >= ub - feastol;
    if (atUpper && !bestlAtUpper) continue;

    double mju = ub * vals[j];
    HighsCDouble mu = mju - lambda;

    if (mu <= 10 * feastol) continue;
    if (std::abs(vals[j]) < 1000 * feastol) continue;

    double mudival = double(mu / vals[j]);
    if (HighsIntegers::isIntegral(mudival, feastol)) continue;
    double eta = ceil(mudival);

    HighsCDouble ulminusetaplusone = HighsCDouble(ub) - eta + 1.0;
    HighsCDouble cplusthreshold = ulminusetaplusone * vals[j];

    HighsInt cplusend =
        std::upper_bound(cover.begin(), cover.end(), double(cplusthreshold),
                         [&](double cplusthreshold, HighsInt i) {
                           return cplusthreshold > vals[i];
                         }) -
        cover.begin();

    HighsCDouble mcplus = m[cplusend];
    if (i < cplusend) mcplus -= mju;

    double jlVal = double(mcplus + eta * vals[j]);

    if (jlVal > bestlVal || (!atUpper && bestlAtUpper)) {
      lpos = i;
      bestlCplusend = cplusend;
      bestlVal = jlVal;
      bestlAtUpper = atUpper;
    }
  }

  if (lpos == -1) return false;

  l = cover[lpos];
  HighsCDouble al = vals[l];
  double upperl = upper[l];
  HighsCDouble mlu = upperl * al;
  HighsCDouble mu = mlu - lambda;

  a.resize(bestlCplusend);
  cover.resize(bestlCplusend);
  u.resize(bestlCplusend + 1);
  m.resize(bestlCplusend + 1);

  if (lpos < bestlCplusend) {
    a.erase(a.begin() + lpos);
    cover.erase(cover.begin() + lpos);
    u.erase(u.begin() + lpos + 1);
    m.erase(m.begin() + lpos + 1);
    for (HighsInt i = lpos + 1; i < bestlCplusend; ++i) {
      u[i] -= upperl;
      m[i] -= mlu;
    }
  }

  HighsInt cplussize = a.size();

  assert(mu > 10 * feastol);

  double mudival = double(mu / al);
  double eta = ceil(mudival);
  HighsCDouble r = mu - floor(mudival) * HighsCDouble(al);
  // we multiply with r and it is important that it does not flip the sign
  // so we safe guard against tiny numerical errors here
  if (r < 0) r = 0;

  HighsCDouble ulminusetaplusone = HighsCDouble(upperl) - eta + 1.0;
  HighsCDouble cplusthreshold = ulminusetaplusone * al;

  HighsInt kmin = floor(eta - upperl - 0.5);

  auto phi_l = [&](double a) {
    assert(a < 0);

    int64_t k = std::min(int64_t(a / double(al)), int64_t(-1));

    for (; k >= kmin; --k) {
      if (a >= k * al + r) {
        assert(a < (k + 1) * al);
        return double(a - (k + 1) * r);
      }

      if (a >= k * al) {
        assert(a < k * al + r);
        return double(k * (al - r));
      }
    }

    assert(a <= -lambda + epsilon);
    return double(kmin * (al - r));
  };

  int64_t kmax = floor(upperl - eta + 0.5);

  auto gamma_l = [&](double z) {
    assert(z > 0);
    for (HighsInt i = 0; i < cplussize; ++i) {
      HighsInt upperi = static_cast<HighsInt>(upper[cover[i]]);

      for (HighsInt h = 0; h <= upperi; ++h) {
        HighsCDouble mih = m[i] + h * a[i];
        HighsCDouble uih = u[i] + h;
        HighsCDouble mihplusdeltai = mih + a[i] - cplusthreshold;
        if (z <= mihplusdeltai) {
          assert(mih <= z + epsilon);
          return double(uih * ulminusetaplusone * (al - r));
        }

        int64_t k = ((int64_t)(double)((z - mihplusdeltai) / al)) - 1;
        for (; k <= kmax; ++k) {
          if (z <= mihplusdeltai + k * al + r) {
            assert(mihplusdeltai + k * al <= z);
            return double((uih * ulminusetaplusone + k) * (al - r));
          }

          if (z <= mihplusdeltai + (k + 1) * al) {
            assert(mihplusdeltai + k * al + r <= z);
            return double((uih * ulminusetaplusone) * (al - r) + z - mih -
                          a[i] + cplusthreshold - (k + 1) * r);
          }
        }
      }
    }

    int64_t p = ((int64_t)(double)((z - m[cplussize]) / al)) - 1;
    for (;; ++p) {
      if (z <= m[cplussize] + p * al + r) {
        assert(m[cplussize] + p * al <= z);
        return double((u[cplussize] * ulminusetaplusone + p) * (al - r));
      }

      if (z <= m[cplussize] + (p + 1) * al) {
        assert(m[cplussize] + p * al + r <= z);
        return double((u[cplussize] * ulminusetaplusone) * (al - r) + z -
                      m[cplussize] - (p + 1) * r);
      }
    }
  };

  rhs = (HighsCDouble(upperl) - eta) * r - lambda;
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt i = 0; i != rowlen; ++i) {
    if (vals[i] == 0.0) continue;
    if (!isintegral[i]) {
      if (vals[i] < 0.0)
        integralSupport = false;
      else
        vals[i] = 0.0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = -phi_l(-vals[i]);
      rhs += vals[i] * upper[i];
    } else {
      vals[i] = gamma_l(vals[i]);
    }
  }

  return true;
}

static double fast_floor(double x) { return (int64_t)x - (x < (int64_t)x); }

bool HighsCutGeneration::cmirCutGenerationHeuristic(double minEfficacy,
                                                    bool onlyInitialCMIRScale) {
  using std::abs;
  using std::floor;
  using std::max;
  using std::sqrt;

  // initialize indicators
  integralSupport = false;
  integralCoefficients = false;

  double continuouscontribution = 0.0;
  double continuoussqrnorm = 0.0;

  deltas.clear();
  deltas.reserve(rowlen + 3);
  integerinds.clear();
  integerinds.reserve(rowlen);

  double maxabsdelta = 0.0;
  constexpr double maxCMirScale = 1e6;
  constexpr double f0min = 0.005;
  constexpr double f0max = 0.995;

  complementation.resize(rowlen);

  for (HighsInt i = 0; i != rowlen; ++i) {
    if (isintegral[i]) {
      integerinds.push_back(i);

      if (upper[i] < 2 * solval[i]) flipComplementation(i);

      if (onlyInitialCMIRScale) continue;

      if (solval[i] > feastol) {
        double delta = abs(vals[i]);
        if (delta <= 1e-4 || delta == maxabsdelta) continue;
        maxabsdelta = max(maxabsdelta, delta);
        deltas.push_back(delta);
      }
    } else {
      updateViolationAndNorm(i, vals[i], continuouscontribution,
                             continuoussqrnorm);
    }
  }

  if (continuoussqrnorm == 0 && deltas.size() > 1) {
    double intScale = HighsIntegers::integralScale(deltas, feastol, kHighsTiny);

    if (intScale != 0.0 && intScale <= 1e4) {
      double scalrhs = double(rhs) * intScale;
      double downrhs = fast_floor(scalrhs);

      double f0 = scalrhs - downrhs;
      if (f0 >= f0min && f0 <= f0max) deltas.push_back(1.0 / intScale);
    }
  }

  deltas.push_back(std::min(1.0, initialScale));
  if (!onlyInitialCMIRScale)
    deltas.push_back(maxabsdelta + std::min(1.0, initialScale));

  if (deltas.empty()) return false;

  pdqsort(deltas.begin(), deltas.end());
  double curdelta = deltas[0];
  for (size_t i = 1; i < deltas.size(); ++i) {
    if (deltas[i] - curdelta <= 10 * feastol)
      deltas[i] = 0.0;
    else
      curdelta = deltas[i];
  }

  deltas.erase(std::remove(deltas.begin(), deltas.end(), 0.0), deltas.end());
  double bestdelta = -1;
  double bestefficacy = minEfficacy;

  for (double delta : deltas) {
    double scale = 1.0 / delta;
    double scalrhs = double(rhs) * scale;
    double downrhs = fast_floor(scalrhs);

    double f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) continue;
    double oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) continue;

    double contscale = scale * oneoveroneminusf0;
    double sqrnorm = contscale * contscale * continuoussqrnorm;
    double viol = contscale * continuouscontribution - downrhs;

    for (HighsInt j : integerinds) {
      double scalaj = vals[j] * scale;
      double downaj = fast_floor(scalaj + kHighsTiny);
      double fj = scalaj - downaj;
      double aj = downaj + std::max(0.0, (fj - f0) * oneoveroneminusf0);
      updateViolationAndNorm(j, aj, viol, sqrnorm);
    }

    // TODO: Add some minimum violation check, e.g., 0.001 * feastol?
    double efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  /* try if multiplying best delta by 2 4 or 8 gives a better efficacy */
  for (HighsInt k = 1; !onlyInitialCMIRScale && k <= 3; ++k) {
    double delta = bestdelta * (1 << k);
    double scale = 1.0 / delta;
    double scalrhs = double(rhs) * scale;
    double downrhs = fast_floor(scalrhs);
    double f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) continue;

    double oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) continue;

    double contscale = scale * oneoveroneminusf0;
    double sqrnorm = contscale * contscale * continuoussqrnorm;
    double viol = contscale * continuouscontribution - downrhs;

    for (HighsInt j : integerinds) {
      double scalaj = vals[j] * scale;
      double downaj = fast_floor(scalaj + kHighsTiny);
      double fj = scalaj - downaj;
      double aj = downaj + std::max(0.0, (fj - f0) * oneoveroneminusf0);
      updateViolationAndNorm(j, aj, viol, sqrnorm);
    }

    double efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  assert(bestdelta != -1);

  // try to flip complementation of integers to increase efficacy
  for (HighsInt k : integerinds) {
    if (upper[k] == kHighsInf) continue;
    if (solval[k] <= feastol) continue;

    flipComplementation(k);

    double delta = bestdelta;
    double scale = 1.0 / delta;
    double scalrhs = double(rhs) * scale;
    double downrhs = fast_floor(scalrhs);

    double f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) {
      flipComplementation(k);
      continue;
    }

    double oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) {
      flipComplementation(k);
      continue;
    }

    double contscale = scale * oneoveroneminusf0;
    double sqrnorm = contscale * contscale * continuoussqrnorm;
    double viol = contscale * continuouscontribution - downrhs;

    for (HighsInt j : integerinds) {
      double scalaj = vals[j] * scale;
      double downaj = fast_floor(scalaj + kHighsTiny);
      double fj = scalaj - downaj;
      double aj = downaj + std::max(0.0, (fj - f0) * oneoveroneminusf0);
      updateViolationAndNorm(j, aj, viol, sqrnorm);
    }

    double efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestefficacy = efficacy;
    } else {
      flipComplementation(k);
    }
  }

  HighsCDouble scale = 1.0 / HighsCDouble(bestdelta);
  HighsCDouble scalrhs = rhs * scale;
  double downrhs = floor(double(scalrhs));

  HighsCDouble f0 = scalrhs - downrhs;
  HighsCDouble oneoveroneminusf0 = 1.0 / (1.0 - f0);

  rhs = downrhs * bestdelta;
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt j = 0; j != rowlen; ++j) {
    if (vals[j] == 0.0) continue;
    if (!isintegral[j]) {
      if (vals[j] > 0.0)
        vals[j] = 0.0;
      else {
        vals[j] = double(vals[j] * oneoveroneminusf0);
        integralSupport = false;
      }
    } else {
      HighsCDouble scalaj = scale * vals[j];
      double downaj = floor(double(scalaj + kHighsTiny));
      HighsCDouble fj = scalaj - downaj;
      HighsCDouble aj = downaj;
      if (fj > f0) aj += (fj - f0) * oneoveroneminusf0;

      vals[j] = double(aj * bestdelta);
    }
  }

#ifndef NDEBUG
  // Check if the computed cut has the correct efficacy
  {
    double checkviol = -downrhs * bestdelta;
    double checknorm = 0.0;
    for (HighsInt j = 0; j != rowlen; ++j) {
      if (vals[j] == 0.0) continue;
      updateViolationAndNorm(j, vals[j], checkviol, checknorm);
    }
    if (checknorm != 0.0) {
      double checkefficacy = checkviol / sqrt(checknorm);
      // the efficacy can become infinite if the cut 0 <= -1 is derived
      assert(fabs(checkefficacy - bestefficacy) < 0.001 ||
             fabs(checkefficacy) >= 1e30);
    }
  }
#endif

  return true;
}

double HighsCutGeneration::scale(double val) {
  int expshift = 0;
  std::frexp(val, &expshift);
  expshift = -expshift;

  // Don't scale the coefficients by more than +1024 (violations can increase)
  expshift = std::min(10, expshift);

  // Scale rhs
  rhs = ldexp(rhs, expshift);

  // Scale row
  for (HighsInt i = 0; i != rowlen; ++i)
    vals[i] = std::ldexp(vals[i], expshift);

  // Return scaling factor
  return std::ldexp(1.0, expshift);
}

bool HighsCutGeneration::postprocessCut(const HighsDomain& globaldom) {
  // right hand sides slightly below zero are likely due to numerical errors and
  // can cause numerical troubles with scaling, so set them to zero
  if (rhs < 0 && rhs > -epsilon) rhs = 0;

  if (integralSupport && integralCoefficients) {
    // if the cut is known to be integral no postprocessing is needed and we
    // simply remove zero coefficients
    for (HighsInt i = rowlen - 1; i >= 0; --i) {
      if (vals[i] == 0.0) {
        --rowlen;
        inds[i] = inds[rowlen];
        vals[i] = vals[rowlen];
      }
    }
    return true;
  }

  // determine maximal absolute coefficient
  double maxAbsValue = 0.0;
  for (HighsInt i = 0; i != rowlen; ++i)
    maxAbsValue = std::max(std::abs(vals[i]), maxAbsValue);

  // determine minimal allowed coefficient
  double minCoefficientValue = 100 * feastol * std::max(maxAbsValue, 1e-3);

  // remove small coefficients and check whether the remaining support is
  // integral
  integralSupport = true;
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0) continue;
    if (std::abs(vals[i]) <= minCoefficientValue) {
      if (vals[i] < 0) {
        double ub = globaldom.col_upper_[inds[i]];
        if (ub == kHighsInf)
          return false;
        else
          rhs -= ub * vals[i];
      } else {
        double lb = globaldom.col_lower_[inds[i]];
        if (lb == -kHighsInf)
          return false;
        else
          rhs -= lb * vals[i];
      }

      vals[i] = 0.0;
      continue;
    }

    if (integralSupport && !lpRelaxation.isColIntegral(inds[i]))
      integralSupport = false;
  }

  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    }
  }

  if (rowlen == 0) return false;

  if (integralSupport) {
    // integral support -> determine scale to make all coefficients integral
    double intscale =
        HighsIntegers::integralScale(vals, rowlen, feastol, epsilon);

    bool scaleSmallestValToOne = true;

    if (intscale != 0.0 &&
        intscale * std::max(1.0, maxAbsValue) <= (double)(uint64_t{1} << 52)) {
      // A scale to make all value integral was found. The scale is only
      // rejected if it is in a range where not all integral values are
      // representable in double precision anymore or cannot be correctly
      // rounded by adding 0.5 and casting to an int. The latter starts at 2^52
      // + 1 since adding 0.5 will round upwards to the next even number for
      // that magnitude. So the largest acceptable value is 2^52 and when at
      // most that value we want to always use the scale to adjust the
      // coefficients and right hand side for numerical safety reasons. If the
      // resulting integral values are too large, however, we scale the cut down
      // by shifting the exponent.
      rhs.renormalize();
      rhs *= intscale;
      maxAbsValue = HighsIntegers::nearestInteger(maxAbsValue * intscale);
      for (HighsInt i = 0; i != rowlen; ++i) {
        HighsCDouble scaleval = intscale * HighsCDouble(vals[i]);
        double intval = HighsIntegers::nearestInteger(double(scaleval));
        double delta = double(scaleval - intval);

        vals[i] = intval;

        // if the coefficient would be strengthened by rounding, we add the
        // upperbound constraint to make it exactly integral instead and
        // therefore weaken the right hand side
        if (delta < 0.0) {
          double ub = globaldom.col_upper_[inds[i]];
          if (ub == kHighsInf) return false;

          rhs -= delta * ub;
        } else {
          double lb = globaldom.col_lower_[inds[i]];
          if (lb == -kHighsInf) return false;

          rhs -= delta * lb;
        }
      }

      // finally we can round down the right hand side. Therefore in most cases
      // small errors for which the upper bound constraints where used and the
      // right hand side was weakened, do not weaken the final cut.
      rhs = floor(rhs + feastol);

      if (intscale * maxAbsValue * feastol < 0.5) {
        // integral scale leads to small enough values, accept scale
        scaleSmallestValToOne = false;
        integralCoefficients = true;
      }
    }

    if (scaleSmallestValToOne) {
      // integral scale lead to very large coefficient values. We now shift the
      // exactly integral values down such that the smallest coefficient is
      // around 1
      double minAbsValue = kHighsInf;
      for (HighsInt i = 0; i != rowlen; ++i)
        minAbsValue = std::min(std::abs(vals[i]), minAbsValue);

      scale(minAbsValue - epsilon);
    }
  } else {
    // the support is not integral, scale cut to have the largest coefficient
    // around 1.0
    scale(maxAbsValue - epsilon);
  }

  return true;
}

bool HighsCutGeneration::preprocessBaseInequality(bool& hasUnboundedInts,
                                                  bool& hasGeneralInts,
                                                  bool& hasContinuous) {
  // preprocess the inequality before cut generation
  // 1. Determine the maximal activity to check for trivial redundancy and
  // tighten coefficients
  // 2. Check for presence of continuous variables and unbounded integers as not
  // all methods for cut generation are applicable in that case
  // 3. Remove coefficients that are below the feasibility tolerance to avoid
  // numerical troubles, use bound constraints to cancel them and
  // reject base inequalities where that is not possible due to unbounded
  // variables
  hasUnboundedInts = false;
  hasContinuous = false;
  hasGeneralInts = false;
  HighsInt numZeros = 0;

  double maxact = -feastol;
  double maxAbsVal = 0;
  for (HighsInt i = 0; i < rowlen; ++i)
    maxAbsVal = std::max(std::abs(vals[i]), maxAbsVal);

  initialScale = scale(maxAbsVal);

  isintegral.resize(rowlen);
  for (HighsInt i = 0; i != rowlen; ++i) {
    // we do not want to have integral variables with small coefficients as this
    // may lead to numerical instabilities during cut generation
    // Therefore we relax integral variables with small coefficients to
    // continuous ones because they still might have a non-negligible
    // contribution e.g. when they come from integral rows with coefficients on
    // the larger side. When we relax them to continuous variables they will be
    // complemented so that their solution value is closest to zero and then
    // will be relaxed if their value is positive or their maximal contribution
    // is below feasibility tolerance.
    isintegral[i] =
        lpRelaxation.isColIntegral(inds[i]) && std::abs(vals[i]) > 10 * feastol;

    if (!isintegral[i]) {
      // complement non-integer variable (cmir separation heuristic complements
      // integral variables in the same way)
      if (upper[i] < 2 * solval[i]) {
        if (complementation.empty()) complementation.resize(rowlen);

        flipComplementation(i);
      }

      // relax positive continuous variables and those with small contributions
      if (vals[i] > 0 || std::abs(vals[i]) * upper[i] <= 10 * feastol) {
        // printf("remove: vals[i] = %g  upper[i] = %g\n", vals[i], upper[i]);
        if (vals[i] < 0) {
          if (upper[i] == kHighsInf) return false;
          rhs -= vals[i] * upper[i];
        }

        ++numZeros;
        vals[i] = 0.0;
        continue;
      }

      hasContinuous = true;
    } else {
      if (upper[i] != kHighsInf &&
          std::abs(vals[i] * upper[i]) > kHighsInvMacheps)
        return false;
      if (upper[i] >= kHighsIInf) {
        hasUnboundedInts = true;
        hasGeneralInts = true;
      } else if (upper[i] != 1.0) {
        hasGeneralInts = true;
      }

      if (vals[i] > 0) maxact += vals[i] * upper[i];
    }
  }

  HighsInt maxLen = 100 + 0.15 * (lpRelaxation.numCols());

  if (rowlen - numZeros > maxLen) {
    HighsInt numCancel = rowlen - numZeros - maxLen;
    std::vector<HighsInt> cancelNzs;

    for (HighsInt i = 0; i != rowlen; ++i) {
      double cancelSlack = vals[i] > 0 ? solval[i] : upper[i] - solval[i];
      if (cancelSlack <= feastol) cancelNzs.push_back(i);
    }

    if ((HighsInt)cancelNzs.size() < numCancel) return false;
    if ((HighsInt)cancelNzs.size() > numCancel)
      std::partial_sort(cancelNzs.begin(), cancelNzs.begin() + numCancel,
                        cancelNzs.end(), [&](HighsInt a, HighsInt b) {
                          return std::abs(vals[a]) < std::abs(vals[b]);
                        });

    for (HighsInt i = 0; i < numCancel; ++i) {
      HighsInt j = cancelNzs[i];

      if (vals[j] < 0) {
        rhs -= vals[j] * upper[j];
      } else
        maxact -= vals[j] * upper[j];

      vals[j] = 0.0;
    }

    numZeros += numCancel;
  }

  if (numZeros != 0) {
    // remove zeros in place
    if (complementation.empty()) {
      for (HighsInt i = rowlen - 1; i >= 0; --i) {
        if (vals[i] == 0.0) {
          --rowlen;
          inds[i] = inds[rowlen];
          vals[i] = vals[rowlen];
          upper[i] = upper[rowlen];
          solval[i] = solval[rowlen];
          isintegral[i] = isintegral[rowlen];
          if (--numZeros == 0) break;
        }
      }
    } else {
      for (HighsInt i = rowlen - 1; i >= 0; --i) {
        if (vals[i] == 0.0) {
          --rowlen;
          inds[i] = inds[rowlen];
          vals[i] = vals[rowlen];
          upper[i] = upper[rowlen];
          solval[i] = solval[rowlen];
          isintegral[i] = isintegral[rowlen];
          complementation[i] = complementation[rowlen];
          if (--numZeros == 0) break;
        }
      }
    }
  }

  return maxact > rhs;
}

#if 0
static void checkNumerics(const double* vals, HighsInt len, double rhs) {
  double maxAbsCoef = 0.0;
  double minAbsCoef = kHighsInf;
  HighsCDouble sqrnorm = 0;
  for (HighsInt i = 0; i < len; ++i) {
    sqrnorm += vals[i] * vals[i];
    maxAbsCoef = std::max(std::abs(vals[i]), maxAbsCoef);
    minAbsCoef = std::min(std::abs(vals[i]), minAbsCoef);
  }

  double norm = double(sqrt(sqrnorm));

  // printf("length: %" HIGHSINT_FORMAT
  //       ", minCoef: %g, maxCoef, %g, norm %g, rhs: %g, dynamism=%g\n",
  //       len, minAbsCoef, maxAbsCoef, norm, rhs, maxAbsCoef / minAbsCoef);
}
#endif

bool HighsCutGeneration::generateCut(HighsTransformedLp& transLp,
                                     std::vector<HighsInt>& inds_,
                                     std::vector<double>& vals_, double& rhs_,
                                     bool onlyInitialCMIRScale) {
#if 0
  if (vals_.size() > 1) {
    std::vector<HighsInt> indsCheck_ = inds_;
    std::vector<double> valsCheck_ = vals_;
    double tmprhs_ = rhs_;
    bool intsPositive = true;
    if (!transLp.transform(valsCheck_, upper, solval, indsCheck_, tmprhs_,
                           intsPositive))
      return false;

    rowlen = indsCheck_.size();
    this->inds = indsCheck_.data();
    this->vals = valsCheck_.data();
    this->rhs = tmprhs_;
    complementation.clear();
    bool hasUnboundedInts = false;
    bool hasGeneralInts = false;
    bool hasContinuous = false;
    // printf("before preprocessing of base inequality:\n");
    checkNumerics(vals, rowlen, double(rhs));
    if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                  hasContinuous))
      return false;
    // printf("after preprocessing of base inequality:\n");
    checkNumerics(vals, rowlen, double(rhs));

    tmprhs_ = (double)rhs;
    valsCheck_.resize(rowlen);
    indsCheck_.resize(rowlen);
    if (!transLp.untransform(valsCheck_, indsCheck_, tmprhs_)) return false;

    // printf("after untransform of base inequality:\n");
    checkNumerics(vals, rowlen, double(rhs));

    // finally check whether the cut is violated
    rowlen = indsCheck_.size();
    inds = indsCheck_.data();
    vals = valsCheck_.data();
    lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
        inds, vals, rowlen, tmprhs_);
  }
#endif

  bool intsPositive = true;
  if (!transLp.transform(vals_, upper, solval, inds_, rhs_, intsPositive))
    return false;

  rowlen = inds_.size();
  this->inds = inds_.data();
  this->vals = vals_.data();
  this->rhs = rhs_;
  complementation.clear();
  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;
  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  // it can happen that there is an unbounded integer variable during the
  // transform call so that the integers are not transformed to positive values.
  // Now the call to preprocessBaseInequality may have removed the unbounded
  // integer, e.g. due to a small coefficient value, so that we can still use
  // the lifted inequalities instead of cmir. We need to make sure, however,
  // that the cut values are transformed to positive coefficients first, which
  // we do below.
  if (!hasUnboundedInts && !intsPositive) {
    complementation.resize(rowlen);

    for (HighsInt i = 0; i != rowlen; ++i) {
      if (vals[i] > 0 || !isintegral[i]) continue;

      flipComplementation(i);
    }
  }

  // try to generate a cut
  if (!tryGenerateCut(inds_, vals_, hasUnboundedInts, hasGeneralInts,
                      hasContinuous, 10 * feastol, onlyInitialCMIRScale))
    return false;

  // remove the complementation if exists
  removeComplementation();

  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    }
  }

  // transform the cut back into the original space, i.e. remove the bound
  // substitution and replace implicit slack variables
  rhs_ = (double)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);
  if (!transLp.untransform(vals_, inds_, rhs_)) return false;

  rowlen = inds_.size();
  inds = inds_.data();
  vals = vals_.data();
  rhs = rhs_;

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(inds, vals,
                                                               rowlen, rhs_);
  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(transLp.getGlobaldom())) return false;
  rhs_ = (double)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), rowlen, rhs_);

  // finally determine the violation of the cut in the original space
  HighsCDouble violation = -rhs_;
  const auto& sol = lpRelaxation.getSolution().col_value;
  for (HighsInt i = 0; i != rowlen; ++i) violation += sol[inds[i]] * vals_[i];

  if (violation <= 10 * feastol) return false;

  transLp.getGlobaldom().tightenCoefficients(inds, vals, rowlen, rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                     vals_.data(), inds_.size(), rhs_,
                                     integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::generateConflict(const HighsDomain& localdomain,
                                          const HighsDomain& globaldom,
                                          std::vector<HighsInt>& proofinds,
                                          std::vector<double>& proofvals,
                                          double& proofrhs) {
  this->inds = proofinds.data();
  this->vals = proofvals.data();
  this->rhs = proofrhs;
  rowlen = proofinds.size();

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds, vals, rowlen, proofrhs);

  complementation.assign(rowlen, false);

  upper.resize(rowlen);
  solval.resize(rowlen);

  double activity = 0.0;
  for (HighsInt i = 0; i != rowlen; ++i) {
    HighsInt col = inds[i];

    upper[i] = globaldom.col_upper_[col] - globaldom.col_lower_[col];

    solval[i] =
        vals[i] < 0
            ? std::min(globaldom.col_upper_[col], localdomain.col_upper_[col])
            : std::max(globaldom.col_lower_[col], localdomain.col_lower_[col]);
    if (vals[i] < 0 && globaldom.col_upper_[col] != kHighsInf) {
      rhs -= globaldom.col_upper_[col] * vals[i];
      vals[i] = -vals[i];
      complementation[i] = true;

      solval[i] = globaldom.col_upper_[col] - solval[i];
    } else {
      rhs -= globaldom.col_lower_[col] * vals[i];
      complementation[i] = false;
      solval[i] = solval[i] - globaldom.col_lower_[col];
    }

    activity += solval[i] * vals[i];
  }

  if (activity > rhs) {
    double solScale = double(rhs) / activity;
    for (HighsInt i = 0; i != rowlen; ++i) solval[i] *= solScale;
  }

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;

  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  // try to generate a cut
  if (!tryGenerateCut(proofinds, proofvals, hasUnboundedInts, hasGeneralInts,
                      hasContinuous, feastol, false, false, false))
    return false;

  // remove the complementation
  if (!complementation.empty()) {
    for (HighsInt i = 0; i != rowlen; ++i) {
      if (complementation[i]) {
        rhs -= globaldom.col_upper_[inds[i]] * vals[i];
        vals[i] = -vals[i];
      } else
        rhs += globaldom.col_lower_[inds[i]] * vals[i];
    }
  }

  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(globaldom)) return false;

  proofvals.resize(rowlen);
  proofinds.resize(rowlen);
  proofrhs = (double)rhs;

  bool cutintegral = integralSupport && integralCoefficients;

  globaldom.tightenCoefficients(proofinds.data(), proofvals.data(), rowlen,
                                proofrhs);

  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(),
                                     proofinds.data(), proofvals.data(), rowlen,
                                     proofrhs, cutintegral, true, true, true);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::finalizeAndAddCut(const HighsDomain& globaldom,
                                           std::vector<HighsInt>& inds_,
                                           std::vector<double>& vals_,
                                           double& rhs_) {
  complementation.clear();
  rowlen = inds_.size();
  inds = inds_.data();
  vals = vals_.data();
  rhs = rhs_;

  integralSupport = true;
  integralCoefficients = false;
  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    } else {
      integralSupport &= lpRelaxation.isColIntegral(inds[i]);
    }
  }

  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(inds, vals,
                                                               rowlen, rhs_);
  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut(globaldom)) return false;
  rhs_ = (double)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), rowlen, rhs_);

  // finally determine the violation of the cut in the original space
  HighsCDouble violation = -rhs_;
  const auto& sol = lpRelaxation.getSolution().col_value;
  for (HighsInt i = 0; i != rowlen; ++i) violation += sol[inds[i]] * vals_[i];

  if (violation <= 10 * feastol) return false;

  globaldom.tightenCoefficients(inds, vals, rowlen, rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                     vals_.data(), inds_.size(), rhs_,
                                     integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

void HighsCutGeneration::flipComplementation(HighsInt index) {
  // only variables with finite upper bounds can be complemented
  assert(upper[index] != kHighsInf);

  // flip complementation
  complementation[index] = !complementation[index];
  solval[index] = upper[index] - solval[index];
  rhs -= upper[index] * vals[index];
  vals[index] = -vals[index];
}

void HighsCutGeneration::removeComplementation() {
  // remove the complementation if it exists
  if (complementation.empty()) return;
  for (HighsInt i = 0; i != rowlen; ++i) {
    if (complementation[i]) flipComplementation(i);
  }
}

void HighsCutGeneration::updateViolationAndNorm(HighsInt index, double aj,
                                                double& violation,
                                                double& norm) const {
  // update violation
  violation += aj * solval[index];

  // update norm
  if (aj > 0 && solval[index] <= feastol) return;
  if (aj < 0 && solval[index] >= upper[index] - feastol) return;
  norm += aj * aj;
}

bool HighsCutGeneration::tryGenerateCut(std::vector<HighsInt>& inds_,
                                        std::vector<double>& vals_,
                                        bool hasUnboundedInts,
                                        bool hasGeneralInts, bool hasContinuous,
                                        double minEfficacy,
                                        bool onlyInitialCMIRScale,
                                        bool allowRejectCut, bool lpSol) {
  // use cmir if there are unbounded integer variables
  if (hasUnboundedInts)
    return cmirCutGenerationHeuristic(minEfficacy, onlyInitialCMIRScale);

  // 0. Save data before determining cover and applying lifting functions
  tmpVals.assign(vals, vals + rowlen);
  tmpInds.assign(inds, inds + rowlen);
  tmpComplementation = complementation;
  tmpSolval = solval;
  HighsCDouble tmpRhs = rhs;

  // 1. Determine a cover, cover does not need to be minimal as neither of
  //    the lifting functions have minimality of the cover as necessary facet
  //    condition
  bool success = false;
  bool saveIntegalSupport = false;
  bool saveIntegralCoefficients = false;
  do {
    if (!determineCover(lpSol)) break;

    // 2. use superadditive lifting function depending on structure of base
    //    inequality:
    //    We have 3 lifting functions available for pure binary knapsack sets,
    //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
    if (!hasContinuous && !hasGeneralInts) {
      separateLiftedKnapsackCover();
      success = true;
    } else if (hasGeneralInts) {
      success = separateLiftedMixedIntegerCover();
    } else {
      assert(hasContinuous);
      assert(!hasGeneralInts);
      success = separateLiftedMixedBinaryCover();
    }
  } while (false);

  double minMirEfficacy = minEfficacy;
  if (success) {
    // save data that might otherwise be overwritten when calling the cmir
    // separator
    saveIntegalSupport = integralSupport;
    saveIntegralCoefficients = integralCoefficients;

    // compute violation and squared norm
    double violation = -double(rhs);
    double sqrnorm = 0.0;
    for (HighsInt i = 0; i < rowlen; ++i) {
      updateViolationAndNorm(i, vals[i], violation, sqrnorm);
    }

    // compute efficacy (distance cut off)
    double efficacy = violation / std::sqrt(sqrnorm);
    if (allowRejectCut && efficacy <= minEfficacy) {
      // reject cut
      success = false;
      rhs = tmpRhs;
    } else {
      // accept cut and increase minimum efficiency requirement for cmir cut
      minMirEfficacy += efficacy;
      std::swap(tmpRhs, rhs);
    }
  }

  // restore indices and values; lifting methods do not modify complementation
  // and, thus, complementation-related data does not have to be restored here.
  inds = tmpInds.data();
  vals = tmpVals.data();

  if (cmirCutGenerationHeuristic(minMirEfficacy, onlyInitialCMIRScale)) {
    // take the cmir cut as it is better
    inds_.swap(tmpInds);
    vals_.swap(tmpVals);
    inds = inds_.data();
    vals = vals_.data();
    return true;
  } else if (success) {
    // take the previous lifted cut as cmir could not improve
    // we restore the old complementation vector, right hand side, and values
    rhs = tmpRhs;
    // note that the solution vector solval also needs to be restored because it
    // depends on the complementation. it would be OK not to restore solval, if
    // there would be a guarantee that it is not used from here on.
    complementation.swap(tmpComplementation);
    solval.swap(tmpSolval);
    inds = inds_.data();
    vals = vals_.data();
    // restore indicators
    integralSupport = saveIntegalSupport;
    integralCoefficients = saveIntegralCoefficients;
    return true;
  } else
    // neither cmir nor lifted cut successful
    return false;
}
