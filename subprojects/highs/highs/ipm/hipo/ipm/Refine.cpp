#include "Solver.h"

// Iterative refinement on the 6x6 system

namespace hipo {

void Solver::refine(NewtonDir& delta) {
  Clock clock;

  NewtonDir correction(m_, n_);
  NewtonDir temp(m_, n_);

  // compute the residuals of the linear system in it_->ires
  clock.start();
  it_->residuals6x6(delta);
  info_.times[kRefinementResTime] += clock.stop();

  clock.start();
  double omega = computeOmega(delta);
  info_.times[kRefinementOmegaTime] += clock.stop();

  double old_omega{};

  for (Int iter = 0; iter < kRefineMaxIter; ++iter) {
    if (omega < kRefineTol) break;

    correction.clear();
    solve6x6(correction, it_->ires);

    temp = delta;
    temp.add(correction);

    clock.start();
    it_->residuals6x6(temp);
    info_.times[kRefinementResTime] += clock.stop();

    old_omega = omega;

    clock.start();
    omega = computeOmega(temp);
    info_.times[kRefinementOmegaTime] += clock.stop();

    if (omega < old_omega) {
      delta = temp;
    } else {
      omega = old_omega;
      break;
    }
  }

  // save worst residual seen in this ipm iteration
  it_->data.back().omega = std::max(it_->data.back().omega, omega);
}

static void updateOmega(double tau, double& omega1, double& omega2, double num,
                        double den1, double den2, double den3) {
  if (den1 + den2 > tau) {
    double omega = num / (den1 + den2);
    omega1 = std::max(omega1, omega);
  } else {
    double omega = num / (den1 + den3);
    omega2 = std::max(omega2, omega);
  }
}

double Solver::computeOmega(const NewtonDir& delta) const {
  // Evaluate iterative refinement progress. Based on "Solving sparse linear
  // systems with sparse backward error", Arioli, Demmel, Duff.

  // Compute |A| * |dx| and |A^T| * |dy|
  std::vector<double> abs_prod_A(m_);
  std::vector<double> abs_prod_At(n_);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = model_.A().start_[col]; el < model_.A().start_[col + 1];
         ++el) {
      Int row = model_.A().index_[el];
      double val = model_.A().value_[el];
      abs_prod_A[row] += std::abs(val) * std::abs(delta.x[col]);
      abs_prod_At[col] += std::abs(val) * std::abs(delta.y[row]);
    }
  }

  double inf_norm_delta{};
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.x));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xu));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.y));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zu));

  // rhs is in it_->res
  // residual of linear system is in it_->ires

  // tau_i =
  // tau_const * (inf_norm_rows(big 6x6 matrix)_i * inf_norm_delta + |rhs_i|)

  const double tau_const = 1000.0 * (5 * n_ + m_) * 1e-16;
  double omega_1{}, omega_2{};

  assert(it_->Rd || m_ == 0);

  // First block
  for (Int i = 0; i < m_; ++i) {
    const double tau =
        tau_const * (std::max(model_.infNormRows(i), std::abs(it_->Rd[i])) *
                         inf_norm_delta +
                     std::abs(it_->res.r1[i]));
    updateOmega(
        tau, omega_1, omega_2,
        // numerator
        std::abs(it_->ires.r1[i]),
        // denominators
        abs_prod_A[i] + std::abs(it_->Rd[i]) * std::abs(delta.y[i]),
        std::abs(it_->res.r1[i]),
        (model_.oneNormRows(i) + std::abs(it_->Rd[i])) * inf_norm_delta);
  }

  // Second block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      const double tau =
          tau_const * (1.0 * inf_norm_delta + std::abs(it_->res.r2[i]));
      updateOmega(tau, omega_1, omega_2,
                  // numerator
                  std::abs(it_->ires.r2[i]),
                  // denominators
                  std::abs(delta.x[i]) + std::abs(delta.xl[i]),
                  std::abs(it_->res.r2[i]), 2 * inf_norm_delta

      );
    }
  }

  // Third block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      const double tau =
          tau_const * (1.0 * inf_norm_delta + std::abs(it_->res.r3[i]));
      updateOmega(tau, omega_1, omega_2,
                  // numerator
                  std::abs(it_->ires.r3[i]),
                  // denominators
                  std::abs(delta.x[i]) + std::abs(delta.xu[i]),
                  std::abs(it_->res.r3[i]), 2 * inf_norm_delta);
    }
  }

  // Fourth block
  for (Int i = 0; i < n_; ++i) {
    double reg_p = it_->Rp ? it_->Rp[i] : regul_.primal;

    const double tau =
        tau_const *
        (std::max(std::max(model_.infNormCols(i), 1.0), std::abs(reg_p)) *
             inf_norm_delta +
         std::abs(it_->res.r4[i]));

    double denom1 = abs_prod_At[i];
    if (model_.hasLb(i)) denom1 += std::abs(delta.zl[i]);
    if (model_.hasUb(i)) denom1 += std::abs(delta.zu[i]);
    denom1 += std::abs(reg_p) * std::abs(delta.x[i]);

    updateOmega(tau, omega_1, omega_2,
                // numerator
                std::abs(it_->ires.r4[i]),
                // denominators
                denom1, std::abs(it_->res.r4[i]),
                (model_.oneNormCols(i) + 2 + std::abs(reg_p)) * inf_norm_delta);
  }

  // Fifth block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      const double tau =
          tau_const * (std::max(std::abs(it_->xl[i]), std::abs(it_->zl[i])) *
                           inf_norm_delta +
                       std::abs(it_->res.r5[i]));

      updateOmega(
          tau, omega_1, omega_2,
          // numerator
          std::abs(it_->ires.r5[i]),
          // denominators
          std::abs(it_->zl[i]) * std::abs(delta.xl[i]) +
              std::abs(it_->xl[i]) * std::abs(delta.zl[i]),
          std::abs(it_->res.r5[i]),
          (std::abs(it_->zl[i]) + std::abs(it_->xl[i])) * inf_norm_delta);
    }
  }

  // Sixth block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      const double tau =
          tau_const * (std::max(std::abs(it_->xu[i]), std::abs(it_->zu[i])) *
                           inf_norm_delta +
                       std::abs(it_->res.r6[i]));

      updateOmega(
          tau, omega_1, omega_2,
          // numerator
          std::abs(it_->ires.r6[i]),
          // denominators
          std::abs(it_->zu[i]) * std::abs(delta.xu[i]) +
              std::abs(it_->xu[i]) * std::abs(delta.zu[i]),
          std::abs(it_->res.r6[i]),
          (std::abs(it_->zu[i]) + std::abs(it_->xu[i])) * inf_norm_delta);
    }
  }

  return omega_1 + omega_2;
}

}  // namespace hipo