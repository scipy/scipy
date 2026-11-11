#include "Iterate.h"

#include "Parameters.h"
#include "ipm/IpxWrapper.h"
#include "ipm/hipo/factorhighs/FactorHighsSettings.h"
#include "model/HighsHessianUtils.h"

namespace hipo {

NewtonDir::NewtonDir(Int m, Int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

void NewtonDir::clear() {
  std::fill(x.begin(), x.end(), 0.0);
  std::fill(y.begin(), y.end(), 0.0);
  std::fill(xl.begin(), xl.end(), 0.0);
  std::fill(xu.begin(), xu.end(), 0.0);
  std::fill(zl.begin(), zl.end(), 0.0);
  std::fill(zu.begin(), zu.end(), 0.0);
}

void NewtonDir::add(const NewtonDir& d) {
  vectorAdd(x, 1.0, d.x, 1.0);
  vectorAdd(y, 1.0, d.y, 1.0);
  vectorAdd(xl, 1.0, d.xl, 1.0);
  vectorAdd(xu, 1.0, d.xu, 1.0);
  vectorAdd(zl, 1.0, d.zl, 1.0);
  vectorAdd(zu, 1.0, d.zu, 1.0);
}

Iterate::Iterate(Model& model_input, Info& info_input, Regularisation& r)
    : model{model_input},
      info{info_input},
      delta(model.m(), model.n()),
      regul{r} {
  clearIter();
  clearRes();
  clearIres();
  best_mu = 0;
}

bool Iterate::isNan() const {
  return (isNanVector(x) || isNanVector(xl) || isNanVector(xu) ||
          isNanVector(y) || isNanVector(zl) || isNanVector(zu));
}
bool Iterate::isInf() const {
  return (isInfVector(x) || isInfVector(xl) || isInfVector(xu) ||
          isInfVector(y) || isInfVector(zl) || isInfVector(zu));
}
bool Iterate::isResNan() const {
  return (isNanVector(res.r1) || isNanVector(res.r2) || isNanVector(res.r3) ||
          isNanVector(res.r4) || isNanVector(res.r5) || isNanVector(res.r6));
}
bool Iterate::isResInf() const {
  return (isInfVector(res.r1) || isInfVector(res.r2) || isInfVector(res.r3) ||
          isInfVector(res.r4) || isInfVector(res.r5) || isInfVector(res.r6));
}
bool Iterate::isDirNan(const NewtonDir& d) const {
  return (isNanVector(d.x) || isNanVector(d.xl) || isNanVector(d.xu) ||
          isNanVector(d.y) || isNanVector(d.zl) || isNanVector(d.zu));
}
bool Iterate::isDirInf(const NewtonDir& d) const {
  return (isInfVector(d.x) || isInfVector(d.xl) || isInfVector(d.xu) ||
          isInfVector(d.y) || isInfVector(d.zl) || isInfVector(d.zu));
}

void Iterate::computeMu() {
  mu = 0.0;
  Int number_finite_bounds{};
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) {
      mu += xl[i] * zl[i];
      ++number_finite_bounds;
    }
    if (model.hasUb(i)) {
      mu += xu[i] * zu[i];
      ++number_finite_bounds;
    }
  }
  mu /= number_finite_bounds;

  if (best_mu > 0.0)
    best_mu = std::min(best_mu, mu);
  else
    best_mu = mu;
}
void Iterate::computeScaling() {
  scaling.assign(model.n(), 0.0);

  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) scaling[i] += zl[i] / xl[i];
    if (model.hasUb(i)) scaling[i] += zu[i] / xu[i];

    // slow down the growth of theta
    if (scaling[i] < 1e-12) scaling[i] = sqrt(1e-12 * scaling[i]);
  }

  // compute extremes of theta
  data.back().min_theta = kHighsInf;
  data.back().max_theta = 0.0;
  for (double d : scaling) {
    if (d != 0.0) {
      data.back().min_theta = std::min(data.back().min_theta, 1.0 / d);
      data.back().max_theta = std::max(data.back().max_theta, 1.0 / d);
    }
  }
}
void Iterate::products() {
  double& min_prod = data.back().min_prod;
  double& max_prod = data.back().max_prod;
  Int& num_small = data.back().num_small_prod;
  Int& num_large = data.back().num_large_prod;

  min_prod = kHighsInf;
  max_prod = 0.0;
  num_small = 0;
  num_large = 0;

  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) {
      double prod = xl[i] * zl[i] / mu;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
    if (model.hasUb(i)) {
      double prod = xu[i] * zu[i] / mu;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
  }
}

void Iterate::indicators() {
  primalObj();
  dualObj();
  primalInfeasUnscaled();
  dualInfeasUnscaled();
  pdGap();
  products();
}

void Iterate::primalObj() {
  pobj = model.offset() + dotProd(x, model.c());
  if (model.qp()) pobj += model.sense() * model.Q().objectiveValue(x);
}
void Iterate::dualObj() {
  dobj = model.offset() + dotProd(y, model.b());
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) dobj += model.lb(i) * zl[i];
    if (model.hasUb(i)) dobj -= model.ub(i) * zu[i];
  }
  if (model.qp()) dobj -= model.sense() * model.Q().objectiveValue(x);
}
void Iterate::pdGap() {
  // relative primal-dual gap
  // Uses the offset-free objectives.
  // Assume pobj = offset + p, dobj = offset + d.
  // If |offset| >> |p|,|d| then |pobj+dobj| ~ 2*offset and the gap gets small
  // much earlier than it should.
  // If offset ~ -p, -d, then |pobj+dobj| ~ 0, so the gap gets small much later
  // than it should.
  pdgap = std::abs(pobj - dobj) /
          (1.0 + 0.5 * std::abs(pobj + dobj - 2 * model.offset()));
}

void Iterate::primalInfeas() {
  // relative infinity norm of scaled primal residuals
  pinf = infNorm(res.r1);
  pinf = std::max(pinf, infNorm(res.r2));
  pinf = std::max(pinf, infNorm(res.r3));
  pinf /= (1.0 + model.normScaledRhs());
}
void Iterate::dualInfeas() {
  // relative infinity norm of scaled dual residual
  dinf = infNorm(res.r4) / (1.0 + model.normScaledObj());
}
void Iterate::primalInfeasUnscaled() {
  // relative infinity norm of unscaled primal residuals
  pinf = 0.0;
  for (Int i = 0; i < model.m(); ++i) {
    double val = std::abs(res.r1[i]);
    if (model.scaled()) val /= model.rowScale(i);
    pinf = std::max(pinf, val);
  }
  for (Int i = 0; i < model.n(); ++i) {
    double val = std::abs(res.r2[i]);
    if (model.scaled()) val *= model.colScale(i);
    pinf = std::max(pinf, val);

    val = std::abs(res.r3[i]);
    if (model.scaled()) val *= model.colScale(i);
    pinf = std::max(pinf, val);
  }
  pinf /= (1.0 + model.normUnscaledRhs());
}
void Iterate::dualInfeasUnscaled() {
  // relative infinity norm of unscaled dual residual
  dinf = 0.0;
  for (Int i = 0; i < model.n(); ++i) {
    double val = std::abs(res.r4[i]);
    if (model.scaled()) val /= model.colScale(i);
    dinf = std::max(dinf, val);
  }
  dinf /= (1.0 + model.normUnscaledObj());
}

void Iterate::residual1234() {
  Clock clock;

  // res1
  res.r1 = model.b();
  model.A().alphaProductPlusY(-1.0, x, res.r1);

  // res2
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i))
      res.r2[i] = model.lb(i) - x[i] + xl[i];
    else
      res.r2[i] = 0.0;
  }

  // res3
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasUb(i))
      res.r3[i] = model.ub(i) - x[i] - xu[i];
    else
      res.r3[i] = 0.0;
  }

  // res4
  res.r4 = model.c();
  model.A().alphaProductPlusY(-1.0, y, res.r4, true);
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) res.r4[i] -= zl[i];
    if (model.hasUb(i)) res.r4[i] += zu[i];
  }
  if (model.qp()) model.Q().alphaProductPlusY(model.sense(), x, res.r4);

  info.times[kResidualsTime] += clock.stop();
}
void Iterate::residual56(double sigma) {
  Clock clock;

  for (Int i = 0; i < model.n(); ++i) {
    // res5
    if (model.hasLb(i))
      res.r5[i] = sigma * mu - xl[i] * zl[i];
    else
      res.r5[i] = 0.0;

    // res6
    if (model.hasUb(i))
      res.r6[i] = sigma * mu - xu[i] * zu[i];
    else
      res.r6[i] = 0.0;
  }

  info.times[kResidualsTime] += clock.stop();
}

std::vector<double> Iterate::residual7(const Residuals& r) const {
  std::vector<double> res7(r.r4);
  for (Int i = 0; i < model.n(); ++i) {
    if (model.hasLb(i)) res7[i] -= ((r.r5[i] + zl[i] * r.r2[i]) / xl[i]);
    if (model.hasUb(i)) res7[i] += ((r.r6[i] - zu[i] * r.r3[i]) / xu[i]);
  }
  return res7;
}
std::vector<double> Iterate::residual8(const Residuals& r,
                                       const std::vector<double>& res7) const {
  std::vector<double> res8(r.r1);
  std::vector<double> temp(res7);

  // temp = (Theta^-1+Rp+Q)^-1 * res7
  for (Int i = 0; i < model.n(); ++i) {
    double denom = scaling[i] + regul.primal;
    if (model.qp()) denom += model.sense() * model.Q().diag(i);
    temp[i] /= denom;
  }

  // res8 += A * temp
  model.A().alphaProductPlusY(1.0, temp, res8);

  return res8;
}

void Iterate::clearIter() {
  x.assign(model.n(), 0.0);
  xl.assign(model.n(), 0.0);
  xu.assign(model.n(), 0.0);
  y.assign(model.m(), 0.0);
  zl.assign(model.n(), 0.0);
  zu.assign(model.n(), 0.0);
}
void Iterate::clearRes() {
  res.r1.assign(model.m(), 0.0);
  res.r2.assign(model.n(), 0.0);
  res.r3.assign(model.n(), 0.0);
  res.r4.assign(model.n(), 0.0);
  res.r5.assign(model.n(), 0.0);
  res.r6.assign(model.n(), 0.0);
}
void Iterate::clearDir() { delta.clear(); }
void Iterate::clearIres() {
  ires.r1.assign(model.m(), 0.0);
  ires.r2.assign(model.n(), 0.0);
  ires.r3.assign(model.n(), 0.0);
  ires.r4.assign(model.n(), 0.0);
  ires.r5.assign(model.n(), 0.0);
  ires.r6.assign(model.n(), 0.0);
}

void Iterate::dropToComplementarity(std::vector<double>& x_cmp,
                                    std::vector<double>& y_cmp,
                                    std::vector<double>& z_cmp) const {
  x_cmp.assign(model.n(), 0.0);
  z_cmp.assign(model.n(), 0.0);
  y_cmp = y;

  for (Int j = 0; j < model.n(); ++j) {
    // value of x_[j] within bounds
    double xj = std::max(x[j], model.lb(j));
    xj = std::min(xj, model.ub(j));

    // FIXED VARIABLE
    if (model.lb(j) == model.ub(j)) {
      x_cmp[j] = model.lb(j);
      z_cmp[j] = zl[j] - zu[j];
    }

    // BOTH BOUNDS FINITE
    else if (model.hasLb(j) && model.hasUb(j)) {
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        // xlj/zlj <= xuj/zuj
        // Primal lower is smaller than primal upper, wrt respective duals
        if (zl[j] >= xl[j]) {
          // drop x to lower bound, set z positive
          x_cmp[j] = model.lb(j);
          z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      } else {
        // xuj/zuj < xlj/zlj
        // Primal upper is smaller than primal lower, wrt respective duals
        if (zu[j] >= xu[j]) {
          // drop x to upper bound, set z negative
          x_cmp[j] = model.ub(j);
          z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      }
    }

    // LOWER BOUND FINITE
    else if (model.hasLb(j)) {
      if (zl[j] >= xl[j]) {
        // drop x to lower bound, set z positive
        x_cmp[j] = model.lb(j);
        z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // UPPER BOUND FINITE
    else if (model.hasUb(j)) {
      if (zu[j] >= xu[j]) {
        // drop x to upper bound, set z negative
        x_cmp[j] = model.ub(j);
        z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // NO BOUNDS
    else {
      x_cmp[j] = xj;
      z_cmp[j] = 0.0;
    }
  }
}

double Iterate::infeasAfterDropping() const {
  // Compute estimate of residuals after dropping to complementarity (taken from
  // ipx).

  double pinf_max = 0.0;
  double dinf_max = 0.0;

  for (Int j = 0; j < model.n(); ++j) {
    double xdrop = 0.0;
    double zdrop = 0.0;

    if (model.hasLb(j) && model.hasUb(j)) {
      // BOTH BOUNDS FINITE
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        if (zl[j] >= xl[j])
          xdrop = x[j] - model.lb(j);
        else
          zdrop = zl[j] - zu[j];
      } else {
        if (zu[j] >= xu[j])
          xdrop = x[j] - model.ub(j);
        else
          zdrop = zl[j] - zu[j];
      }
    }

    // LOWER BOUND FINITE
    else if (model.hasLb(j)) {
      if (zl[j] >= xl[j])
        xdrop = x[j] - model.lb(j);
      else {
        // zu[j] is zero
        zdrop = zl[j];
      }
    }

    // UPPER BOUND FINITE
    else if (model.hasUb(j)) {
      if (zu[j] >= xu[j])
        xdrop = x[j] - model.ub(j);
      else {
        // zl[j] is zero
        zdrop = -zu[j];
      }
    }

    // largest entry in column j of A
    double Amax = 0.0;
    for (Int el = model.A().start_[j]; el < model.A().start_[j + 1]; ++el)
      Amax = std::max(Amax, std::abs(model.A().value_[el]));

    pinf_max = std::max(pinf_max, std::abs(xdrop) * Amax);
    dinf_max = std::max(dinf_max, std::abs(zdrop));
  }

  pinf_max /= (1.0 + model.normScaledRhs());
  dinf_max /= (1.0 + model.normScaledObj());

  return std::max(pinf_max, dinf_max);
}

void Iterate::finalResiduals() const {
  // If ipx has been used, the information is already available, otherwise,
  // compute it.

  if (!info.ipx_used) {
    info.p_obj = pobj;
    info.d_obj = dobj;
    info.pd_gap = pdgap;

    double& pinf = info.p_res_abs;
    for (Int i = 0; i < model.m(); ++i) {
      double val = std::abs(res.r1[i]);
      if (model.scaled()) val /= model.rowScale(i);
      pinf = std::max(pinf, val);
    }
    for (Int i = 0; i < model.n(); ++i) {
      double val = std::abs(res.r2[i]);
      if (model.scaled()) val *= model.colScale(i);
      pinf = std::max(pinf, val);

      val = std::abs(res.r3[i]);
      if (model.scaled()) val *= model.colScale(i);
      pinf = std::max(pinf, val);
    }
    info.p_res_rel = pinf / (1.0 + model.normUnscaledRhs());

    double& dinf = info.d_res_abs;
    for (Int i = 0; i < model.n(); ++i) {
      double val = std::abs(res.r4[i]);
      if (model.scaled()) val /= model.colScale(i);
      dinf = std::max(dinf, val);
    }
    info.d_res_rel = dinf / (1.0 + model.normUnscaledObj());

  } else {
    info.p_res_abs = info.ipx_info.abs_presidual;
    info.p_res_rel = info.ipx_info.rel_presidual;
    info.d_res_abs = info.ipx_info.abs_dresidual;
    info.d_res_rel = info.ipx_info.rel_dresidual;
    info.p_obj = info.ipx_info.pobjval;
    info.d_obj = info.ipx_info.dobjval;
    info.pd_gap = std::abs(info.ipx_info.rel_objgap);
  }
}

void Iterate::getReg(LinearSolver& LS, const std::string& nla) {
  if (nla == kHipoNormalEqString)
    total_reg.resize(model.m());
  else
    total_reg.resize(model.m() + model.n());

  // extract regularisation
  LS.getReg(total_reg);

  // easy access to primal/dual regularisation
  if (nla == kHipoNormalEqString) {
    Rp = nullptr;
    Rd = total_reg.data();
  } else {
    Rp = total_reg.data();
    Rd = model.m() > 0 ? &total_reg[model.n()] : nullptr;
  }
}

void Iterate::residuals6x6(const NewtonDir& d) {
  const std::vector<double>& dx = d.x;
  const std::vector<double>& dy = d.y;
  const std::vector<double>& dxl = d.xl;
  const std::vector<double>& dxu = d.xu;
  const std::vector<double>& dzl = d.zl;
  const std::vector<double>& dzu = d.zu;
  const Int m = model.m();
  const Int n = model.n();
  assert(Rd || m == 0);

  // res1,2,3,4,5,6 contain the rhs of the linear system

  // ires1 = res1 - A * dx - Rd * dy
  ires.r1 = res.r1;
  model.A().alphaProductPlusY(-1.0, dx, ires.r1);
  for (Int i = 0; i < m; ++i) {
    ires.r1[i] -= Rd[i] * dy[i];
  }

  // ires2 = res2 - dx + dxl
  for (Int i = 0; i < n; ++i)
    if (model.hasLb(i))
      ires.r2[i] = res.r2[i] - dx[i] + dxl[i];
    else
      ires.r2[i] = 0.0;

  // ires3 = res3 - dx - dxu
  for (Int i = 0; i < n; ++i)
    if (model.hasUb(i))
      ires.r3[i] = res.r3[i] - dx[i] - dxu[i];
    else
      ires.r3[i] = 0.0;

  // ires4 = res4 - A^T * dy - dzl + dzu + Q * dx + Rp * dx
  ires.r4 = res.r4;
  for (Int i = 0; i < n; ++i) {
    if (model.hasLb(i)) ires.r4[i] -= dzl[i];
    if (model.hasUb(i)) ires.r4[i] += dzu[i];
  }
  model.A().alphaProductPlusY(-1.0, dy, ires.r4, true);
  for (Int i = 0; i < n; ++i) {
    double reg_p = Rp ? Rp[i] : regul.primal;
    ires.r4[i] += reg_p * dx[i];
  }
  if (model.qp()) model.Q().alphaProductPlusY(model.sense(), dx, ires.r4);

  // ires5 = res5 - zl * dxl - xl * dzl
  for (Int i = 0; i < n; ++i) {
    if (model.hasLb(i))
      ires.r5[i] = res.r5[i] - zl[i] * dxl[i] - xl[i] * dzl[i];
    else
      ires.r5[i] = 0.0;
  }

  // ires6 = res6 - zu * dxu - xu * dzu
  for (Int i = 0; i < n; ++i) {
    if (model.hasUb(i))
      ires.r6[i] = res.r6[i] - zu[i] * dxu[i] - xu[i] * dzu[i];
    else
      ires.r6[i] = 0.0;
  }
}

void Iterate::assertConsistency(Int n, Int m) const {
  assert(static_cast<Int>(x.size()) == n);
  assert(static_cast<Int>(xl.size()) == n);
  assert(static_cast<Int>(xu.size()) == n);
  assert(static_cast<Int>(y.size()) == m);
  assert(static_cast<Int>(zl.size()) == n);
  assert(static_cast<Int>(zu.size()) == n);
}

void Iterate::makeStep(double alpha_primal, double alpha_dual) {
  if (std::min(alpha_primal, alpha_dual) < 0.05)
    ++bad_iter_;
  else
    bad_iter_ = 0;

  vectorAdd(x, 1.0, delta.x, alpha_primal);
  vectorAdd(xl, 1.0, delta.xl, alpha_primal);
  vectorAdd(xu, 1.0, delta.xu, alpha_primal);
  vectorAdd(y, 1.0, delta.y, alpha_dual);
  vectorAdd(zl, 1.0, delta.zl, alpha_dual);
  vectorAdd(zu, 1.0, delta.zu, alpha_dual);
}

void Iterate::saveBest(double feas_tol, double opt_tol, Int iter) {
  const double primal_violation = pinf / feas_tol;
  const double dual_violation = dinf / feas_tol;
  const double gap_violation = pdgap / opt_tol;

  double violation = std::max(primal_violation, dual_violation);
  violation = std::max(violation, gap_violation);

  if (violation < best_violation) {
    best_x = x;
    best_xl = xl;
    best_xu = xu;
    best_y = y;
    best_zl = zl;
    best_zu = zu;
    best_violation = violation;
    best_iter = iter;
    model.saveBounds();
  }
}

bool Iterate::resetBest(Int iter) {
  if (iter == best_iter) return false;

  x = best_x;
  xl = best_xl;
  xu = best_xu;
  y = best_y;
  zl = best_zl;
  zu = best_zu;
  model.resetBounds();

  computeMu();
  residual1234();
  indicators();

  return true;
}

bool Iterate::stagnation(std::stringstream& log_stream) {
  // too many iterations in a row with small stepsize
  bool stagnation = (bad_iter_ >= kMaxBadIter);

  // the next tests are aimed at problems where things are going
  // catastrophically bad, which may continue iterating forever otherwise.

  // dx / x or dy / y becoming way too small
  const double thresh_dxy_xy = 1e-30;

  std::vector<double> temp = delta.x;
  vectorDivide(temp, x);
  const double dx_x_max = infNorm(temp);

  temp = delta.y;
  vectorDivide(temp, y);
  const double dy_y_max = infNorm(temp);

  largest_dx_x_ = std::max(largest_dx_x_, dx_x_max);
  largest_dy_y_ = std::max(largest_dy_y_, dy_y_max);

  if (dx_x_max < largest_dx_x_ * thresh_dxy_xy ||
      dy_y_max < largest_dy_y_ * thresh_dxy_xy) {
    stagnation = true;
    log_stream << "Bad direction ratios, dx_x " << sci(dx_x_max, 0, 1)
               << " (max " << sci(largest_dx_x_, 0, 1) << "), dy_y "
               << sci(dy_y_max, 0, 1) << " (max " << sci(largest_dy_y_, 0, 1)
               << ")\n";
  }

  // infeasibilities jumping back up
  const double thresh_inf_to_best = 1e12;

  if (pinf > 0.0) best_pinf_ = std::min(best_pinf_, pinf);
  if (dinf > 0.0) best_dinf_ = std::min(best_dinf_, dinf);

  if (pinf > thresh_inf_to_best * best_pinf_ ||
      dinf > thresh_inf_to_best * best_dinf_) {
    stagnation = true;
    log_stream << "Bad infeasibility, pinf " << sci(pinf, 0, 1) << " (best "
               << sci(best_pinf_, 0, 1) << "), dinf " << sci(dinf, 0, 1)
               << " (best " << sci(best_dinf_, 0, 1) << ")\n";
  }

  return stagnation;
}

}  // namespace hipo