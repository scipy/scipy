#include "ipm/ipx/ipm.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>
#include "ipm/ipx/timer.h"
#include "ipm/ipx/utils.h"

namespace ipx {

struct IPM::Step {
    Step(Int m, Int n) : x(n+m), xl(n+m), xu(n+m), y(m), zl(n+m), zu(n+m) {}
    Vector x, xl, xu, y, zl, zu;
    Step& operator+=(const Step& rhs) {
        x += rhs.x; xl += rhs.xl; xu += rhs.xu;
        y += rhs.y; zl += rhs.zl; zu += rhs.zu;
	return *this;
    }
};

IPM::IPM(const Control& control) : control_(control) {}

void IPM::StartingPoint(KKTSolver* kkt, Iterate* iterate, Info* info) {
    kkt_ = kkt;
    iterate_ = iterate;
    info_ = info;
    PrintHeader();
    ComputeStartingPoint();
    if (info->errflag == 0)
        PrintOutput();
    // Set status_ipm.
    if (info->errflag == IPX_ERROR_user_interrupt) {
        info->errflag = 0;
        info->status_ipm = IPX_STATUS_user_interrupt;
    } else if (info->errflag == IPX_ERROR_time_interrupt) {
        info->errflag = 0;
        info->status_ipm = IPX_STATUS_time_limit;
    } else if (info->errflag) {
        info->status_ipm = IPX_STATUS_failed;
    } else {
        info->status_ipm = IPX_STATUS_not_run;
    }
}

void IPM::Driver(KKTSolver* kkt, Iterate* iterate, Info* info) {
    const Model& model = iterate->model();
    const Int m = model.rows();
    const Int n = model.cols();
    Step step(m, n);
    kkt_ = kkt;
    iterate_ = iterate;
    info_ = info;
    num_bad_iter_ = 0;

    while (true) {
        if (iterate->term_crit_reached()) {
            info->status_ipm = IPX_STATUS_optimal;
            break;
        }
        if (num_bad_iter_ >= 5 ||
            (best_complementarity_ > 0.0 &&
            iterate_->complementarity() > kDivergeTol * best_complementarity_)) {
            // No progress in reducing the complementarity gap.
            // Check if model seems to be primal or dual infeasible.
            bool dualized = iterate_->model().dualized();
            double pobjective = iterate_->pobjective_after_postproc();
            double dobjective = iterate_->dobjective_after_postproc();
            if (iterate_->presidual() > iterate_->feasibility_tol() * (1.0 + model.norm_bounds()) &&
		dobjective > std::max(10.0 * std::abs(pobjective), 1.0)) {
                // Dual objective tends to positive infinity. Looks like the
                // model is dual unbounded, i.e. primal infeasible.
                info->status_ipm = dualized ?
                    IPX_STATUS_dual_infeas : IPX_STATUS_primal_infeas;
            } else if (iterate_->dresidual() > iterate_->feasibility_tol() * (1.0 + model.norm_c()) &&
		       pobjective < -std::max(10.0 * std::abs(dobjective), 1.0)) {
                // Primal objective tends to negative infinity. Looks like the
                // model is primal unbounded, i.e. dual infeasible.
                info->status_ipm = dualized ?
                    IPX_STATUS_primal_infeas : IPX_STATUS_dual_infeas;
            }
            else {
                info->status_ipm = IPX_STATUS_no_progress;
            }
            break;
        }
        if (info->iter >= maxiter_) {
            info->status_ipm = IPX_STATUS_iter_limit;
            break;
        }
        if ((info->errflag = control_.InterruptCheck(info->iter)) != 0)
            break;
        kkt->Factorize(iterate, info);
        if (info->errflag)
            break;
        Predictor(step);
        if (info->errflag)
            break;
        AddCorrector(step);
        if (info->errflag)
            break;
        MakeStep(step);
        info->iter++;
        PrintOutput();
    }

    // Set status_ipm if errflag terminated IPM.
    if (info->errflag) {
        if (info->errflag == IPX_ERROR_user_interrupt) {
	    info->errflag = 0;
	    info->status_ipm = IPX_STATUS_user_interrupt;
	} else if (info->errflag == IPX_ERROR_time_interrupt) {
            info->errflag = 0;
            info->status_ipm = IPX_STATUS_time_limit;
        } else {
            info->status_ipm = IPX_STATUS_failed;
        }
    }

    if (control_.runCentring() &&
	info->status_ipm == IPX_STATUS_optimal && !info->centring_tried) {
      // Centrality of a point is evaluated by the quantities
      //  min (xj * zj) / mu
      //  max (xj * zj) / mu
      // Ideally, they are in the interval [0.1,10.0].
      // As soon as the ratio
      //  max (xj * zj) / min (xj * zj)
      // is below centringRatioTolerance, the point is considered centred.
      // If the new point after centring has a ratio that is lower than the
      // previous ratio times centringRatioReduction, then the step is
      // accepted. Otherwise, the step is rejected and no more centring steps are
      // performed.
      //
      // If IPM is optimal and centring has not yet run, run centring
      // (to avoid running it twice during initial IPM and main IPM).
      control_.hLog("Performing centring steps...\n");

      // freeze mu to its current value
      const double mu_frozen = iterate_->mu();

      // assess and print centrality of current point
      AssessCentrality(iterate_->xl(), iterate_->xu(), iterate_->zl(),
		       iterate_->zu(), iterate_->mu());
      double prev_ratio = centring_ratio;
      Int prev_bad_products = bad_products;

      info->centring_success = false;
      // if ratio is below tolerance, point is centred
      if (prev_ratio < control_.centringRatioTolerance()) {
	control_.hLog("\tPoint is now centred\n");
	info->centring_success = true;
      } else {
	// perform centring steps
	bool centring_complete = false;
	for (int ii = 0; ii < control_.maxCentringSteps(); ++ii) {
	  // compute centring step
	  Centring(step, mu_frozen);
	  
	  // assess whether to take the step
	  bool accept = EvaluateCentringStep(step, prev_ratio, prev_bad_products);
	  if (!accept) {
	    control_.hLog("\tPoint cannot be centred further\n");
	    centring_complete = true;
	    break;
	  }

	  // take the step and print output
	  MakeStep(step, true);
	  info->iter++;
	  PrintOutput();
	  AssessCentrality(iterate_->xl(), iterate_->xu(), iterate_->zl(),
			   iterate_->zu(), iterate_->mu());
	  prev_ratio = centring_ratio;
	  prev_bad_products = bad_products;
	    
	  // if ratio is below tolerance, point is centred
	  if (prev_ratio < control_.centringRatioTolerance()) {
	    control_.hLog("\tPoint is now centred\n");
	    info->centring_success = true;
	    centring_complete = true;
	    break;
	  }
	}
	if (!centring_complete) {
	  std::stringstream h_logging_stream;
	  h_logging_stream.str(std::string());
	  h_logging_stream << "\tPoint could not be centred within "
			   << control_.maxCentringSteps() << " iterations\n";
	  control_.hLog(h_logging_stream);
	}
      }
      info->centring_tried = true;
    } // if (control_.runCentring() && info->status_ipm ==
      // IPX_STATUS_optimal && !info->centring_tried)
}

void IPM::ComputeStartingPoint() {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Vector& b = model.b();
    const Vector& c = model.c();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    Vector x(n+m), xl(n+m), xu(n+m), y(m), zl(n+m), zu(n+m);
    Vector rb(m);               // workspace

    // Factorize the KKT matrix with the identity matrix in the (1,1) block.
    kkt_->Factorize(nullptr, info_);
    if (info_->errflag)
        return;

    // Set x within its bounds and compute the minimum norm solution dx to
    // AI*dx = (b-AI*x). Then update x := x + dx to obtain a feasible point.
    rb = b;
    for (Int j = 0; j < n+m; j++) {
        double xj = 0.0;
        if (xj < lb[j])
            xj = lb[j];
        if (xj > ub[j])
            xj = ub[j];
        x[j] = xj;
        if (xj != 0.0)
            ScatterColumn(AI, j, -xj, rb);
    }
    double tol = 0.1 * Infnorm(rb);
    zl = 0.0;
    kkt_->Solve(zl, rb, tol, xl, y, info_);
    if (info_->errflag)
        return;
    x += xl;

    // Compute xl, xu and shift to become positive.
    double xinfeas = 0.0;
    for (Int j = 0; j < n+m; j++) {
        xl[j] = x[j]-lb[j];
        xinfeas = std::max(xinfeas, -xl[j]);
        xu[j] = ub[j]-x[j];
        xinfeas = std::max(xinfeas, -xu[j]);
    }
    double xshift1 = 1.0 + 1.5*xinfeas;
    xl += xshift1;
    xu += xshift1;

    const double cnorm = Twonorm(c);
    if (cnorm == 0.0) {
        // Special treatment for zero objective.
        for (Int j = 0; j < n+m; j++) {
            zl[j] = std::isfinite(lb[j]) ? 1.0 : 0.0;
            zu[j] = std::isfinite(ub[j]) ? 1.0 : 0.0;
        }
    } else {
        // Compute y as the least-squares solution to AI'*y=c.
        // Recompute zl = c-AI'*y because the KKT system is solved
        // approximately with a residual in the first block equation.
        rb = 0.0;
        double tol = 0.1 * Infnorm(c);
        kkt_->Solve(c, rb, tol, zl, y, info_);
        if (info_->errflag)
            return;
        zl = c;
        MultiplyAdd(AI, y, -1.0, zl, 'T');

        // When c lies in range(AI'), then the dual slack variables are (close
        // to) zero, and the initial point would be almost complementary but
        // usually not primal feasible. To prevent this from happening, add
        // a fraction of the objective to zl and adjust y. In exact computation
        // this does not affect dual feasibility.
        const double znorm = Twonorm(zl);
        const double rho = 0.05;
        if (znorm < rho*cnorm) {
            zl += rho * c;
            y *= (1.0-rho);
        }

        // Split dual slack solution into zl, zu and shift to become positive.
        double zinfeas = 0.0;
        for (Int j = 0; j < n+m; j++) {
            double zval = zl[j];
            zl[j] = 0.0;
            zu[j] = 0.0;
            if (std::isfinite(lb[j]) && std::isfinite(ub[j])) {
                zl[j] = 0.5*zval;
                zu[j] = -0.5*zval;
            }
            else if (std::isfinite(lb[j]))
                zl[j] = zval;
            else if (std::isfinite(ub[j]))
                zu[j] = -zval;
            zinfeas = std::max(zinfeas, -zl[j]);
            zinfeas = std::max(zinfeas, -zu[j]);
        }
        double zshift1 = 1.0 + 1.5*zinfeas;
        for (Int j = 0; j < n+m; j++) {
            if (std::isfinite(lb[j]))
                zl[j] += zshift1;
            if (std::isfinite(ub[j]))
                zu[j] += zshift1;
        }
    }

    // Level pairwise complementarity products.
    double xsum = 1.0;
    double zsum = 1.0;
    double mu = 1.0;
    for (Int j = 0; j < n+m; j++) {
        if (std::isfinite(lb[j])) {
            xsum += xl[j];
            zsum += zl[j];
            mu += xl[j]*zl[j];
        }
        if (std::isfinite(ub[j])) {
            xsum += xu[j];
            zsum += zu[j];
            mu += xu[j]*zu[j];
        }
    }
    double xshift2 = 0.5*mu/zsum;
    double zshift2 = 0.5*mu/xsum;
    xl += xshift2;
    xu += xshift2;
    for (Int j = 0; j < n+m; j++) {
        if (std::isfinite(lb[j]))
            zl[j] += zshift2;
        if (std::isfinite(ub[j]))
            zu[j] += zshift2;
    }
    iterate_->Initialize(x, xl, xu, y, zl, zu);
    best_complementarity_ = iterate_->complementarity();
}

// Computes maximum alpha such that x + alpha*dx >= 0.
// The blocking index is returned in blocking_index if not NULL.
static double StepToBoundary(const Vector& x, const Vector& dx,
                             Int* blocking_index, double alpha = 1.0) {
    const Int n = x.size();
    const double damp = 1.0 - std::numeric_limits<double>::epsilon();
    assert(damp < 1.0);

    Int iblock = -1;
    for (Int i = 0; i < n; i++) {
        assert(x[i] >= 0.0);
        if (x[i]+alpha*dx[i] < 0.0) {
            alpha = -(x[i]*damp) / dx[i];
            assert(x[i]+alpha*dx[i] >= 0.0);
            iblock = i;
        }
    }
    assert(alpha >= 0.0);
    if (blocking_index)
        *blocking_index = iblock;
    return alpha;
}

void IPM::Predictor(Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();

    // sl = -xl.*zl
    Vector sl(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_lb(j))
            sl[j] = -xl[j]*zl[j];
        else
            sl[j] = 0.0;
    assert(AllFinite(sl));

    // su = -xu.*zu
    Vector su(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_ub(j))
            su[j] = -xu[j]*zu[j];
        else
            su[j] = 0.0;
    assert(AllFinite(su));

    SolveNewtonSystem(&iterate_->rb()[0], &iterate_->rc()[0],
                      &iterate_->rl()[0], &iterate_->ru()[0], &sl[0], &su[0],
                      step);
}

void IPM::AddCorrector(Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    const Vector& dxl = step.xl;
    const Vector& dxu = step.xu;
    const Vector& dzl = step.zl;
    const Vector& dzu = step.zu;
    const double mu = iterate_->mu();

    // Choose centering parameter.
    double step_xl = StepToBoundary(xl, dxl, nullptr);
    double step_xu = StepToBoundary(xu, dxu, nullptr);
    double step_zl = StepToBoundary(zl, dzl, nullptr);
    double step_zu = StepToBoundary(zu, dzu, nullptr);
    double maxp = std::min(step_xl, step_xu);
    double maxd = std::min(step_zl, step_zu);
    double muaff = 0.0;
    Int num_finite = 0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->has_barrier_lb(j)) {
            assert(std::isfinite(xl[j]));
            assert(xl[j] != 0.0);
            muaff += (xl[j]+maxp*dxl[j]) * (zl[j]+maxd*dzl[j]);
            num_finite++;
        }
        if (iterate_->has_barrier_ub(j)) {
            assert(std::isfinite(xu[j]));
            assert(xu[j] != 0.0);
            muaff += (xu[j]+maxp*dxu[j]) * (zu[j]+maxd*dzu[j]);
            num_finite++;
        }
    }
    assert(std::isfinite(muaff));
    muaff /= num_finite;
    double ratio = muaff / mu;
    double sigma = ratio * ratio * ratio;

    // sl = -xl.*zl + sigma*mu - dxl.*dzl
    Vector sl(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_lb(j))
            sl[j] = -xl[j]*zl[j] + sigma*mu - dxl[j]*dzl[j];
        else
            sl[j] = 0.0;
    assert(AllFinite(sl));

    // su = -xu.*zu + sigma*mu - dxu.*dzu
    Vector su(n+m);
    for (Int j = 0; j < n+m; j++)
        if (iterate_->has_barrier_ub(j))
            su[j] = -xu[j]*zu[j] + sigma*mu - dxu[j]*dzu[j];
        else
            su[j] = 0.0;
    assert(AllFinite(su));

    SolveNewtonSystem(&iterate_->rb()[0], &iterate_->rc()[0],
                      &iterate_->rl()[0], &iterate_->ru()[0], &sl[0], &su[0],
                      step);
}

void IPM::Centring(Step& step, double mu) {
  const Model& model = iterate_->model();
  const Int m = model.rows();
  const Int n = model.cols();
  const Vector& xl = iterate_->xl();
  const Vector& xu = iterate_->xu();
  const Vector& zl = iterate_->zl();
  const Vector& zu = iterate_->zu();

  Vector sl(n + m);
  Vector su(n + m);

  // Set sigma to 1 for pure centring
  const double sigma = 1.0;

  // sl = -xl.*zl + sigma*mu
  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_lb(j)) {
      sl[j] = -xl[j] * zl[j] + sigma * mu;
    } else {
      sl[j] = 0.0;
    }
  }
  assert(AllFinite(sl));

  // su = -xu.*zu + sigma*mu
  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_ub(j)) {
      su[j] = -xu[j] * zu[j] + sigma * mu;
    } else {
      su[j] = 0.0;
    }
  }
  assert(AllFinite(su));

  SolveNewtonSystem(&iterate_->rb()[0], &iterate_->rc()[0], &iterate_->rl()[0],
                    &iterate_->ru()[0], &sl[0], &su[0], step);
}

void IPM::AssessCentrality(const Vector& xl, const Vector& xu,
                             const Vector& zl, const Vector& zu, double mu,
                             bool print) {
  // The function computes the ratio
  //  min(x_j * z_j) / max(x_j * z_j)
  // and the number of products x_j * z_j that are not in the interval
  //  [0.1 * mu, 10 * mu]
  // and prints information to screen if print is on.

  const Int m = iterate_->model().rows();
  const Int n = iterate_->model().cols();

  double minxz = kHighsInf;
  double maxxz = 0.0;

  const double gamma = 0.1;
  bad_products = 0;

  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_lb(j)) {
      const double product = xl[j] * zl[j];
      if (product < gamma * mu || product > mu / gamma){
        ++bad_products;
      }
      minxz = std::min(minxz, product);
      maxxz = std::max(maxxz, product);
    }
  }

  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_ub(j)) {
      const double product = xu[j] * zu[j];
      if (product < gamma * mu || product > mu / gamma){
        ++bad_products;
      }
      minxz = std::min(minxz, product);
      maxxz = std::max(maxxz, product);
    }
  }

  maxxz = std::max(maxxz, mu);
  minxz = std::min(minxz, mu);

  centring_ratio = maxxz / minxz;

  if (print) {
    std::stringstream h_logging_stream;
    h_logging_stream.str(std::string());
    h_logging_stream << "\txj*zj in [ "
		     << Scientific(minxz / mu, 8, 2) << ", "
		     << Scientific(maxxz / mu, 8, 2) << "]; Ratio = "
		     << Scientific(centring_ratio, 8, 2) << "; (xj*zj / mu) not_in [0.1, 10]: "
		     << bad_products << "\n"; 
    control_.hLog(h_logging_stream);
  }
}

bool IPM::EvaluateCentringStep(const Step& step, double prev_ratio, Int prev_bad) {
  // The function returns true is the step is to be accepted.
  // The step is accepted if the ratio of the new point is not worse 
  // than the previous one times centringRatioReduction or if the 
  // number of outliers products is reduced.

  StepSizes(step, true);

  const Int n = iterate_->model().cols();
  const Int m = iterate_->model().rows();

  Vector xl_temp = iterate_->xl();
  Vector xu_temp = iterate_->xu();
  Vector zl_temp = iterate_->zl();
  Vector zu_temp = iterate_->zu();

  // perform temporary step
  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_lb(j)) {
      xl_temp[j] += step_primal_ * step.xl[j];
    }
    if (iterate_->has_barrier_ub(j)) {
      xu_temp[j] += step_primal_ * step.xu[j];
    }
    if (iterate_->has_barrier_lb(j)) {
      zl_temp[j] += step_dual_ * step.zl[j];
    }
    if (iterate_->has_barrier_ub(j)) {
      zu_temp[j] += step_dual_ * step.zu[j];
    }
  }

  // compute temporary mu
  double mu_temp = 0.0;
  Int num_finite = 0;
  for (Int j = 0; j < n + m; j++) {
    if (iterate_->has_barrier_lb(j)) {
      mu_temp += xl_temp[j] * zl_temp[j];
      ++num_finite;
    }
    if (iterate_->has_barrier_ub(j)) {
      mu_temp += xu_temp[j] * zu_temp[j];
      ++num_finite;
    }
  }
  mu_temp /= num_finite;

  // assess quality of temporary point
  AssessCentrality(xl_temp, xu_temp, zl_temp, zu_temp, mu_temp, false);

  // accept the step if the new ratio is not more than centringRatioReduction
  // times the previous one, or if the new point has fewer outliers
  return (centring_ratio < control_.centringRatioReduction() * prev_ratio || 
          bad_products < prev_bad);
}

void IPM::StepSizes(const Step& step, bool isCentring) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    const Vector& dxl = step.xl;
    const Vector& dxu = step.xu;
    const Vector& dzl = step.zl;
    const Vector& dzu = step.zu;
    const double gammaf = 0.9;
    const double gammaa = 1.0 / (1.0-gammaf);

    Int block_xl, block_xu, block_zl, block_zu;
    double step_xl = StepToBoundary(xl, dxl, &block_xl);
    double step_xu = StepToBoundary(xu, dxu, &block_xu);
    double step_zl = StepToBoundary(zl, dzl, &block_zl);
    double step_zu = StepToBoundary(zu, dzu, &block_zu);
    double maxp = std::fmin(step_xl, step_xu);
    double maxd = std::fmin(step_zl, step_zu);
    double mufull = 0.0;
    Int num_finite = 0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->has_barrier_lb(j)) {
            assert(std::isfinite(xl[j]));
            assert(xl[j] != 0.0);
            mufull += (xl[j]+maxp*dxl[j]) * (zl[j]+maxd*dzl[j]);
            num_finite++;
        }
        if (iterate_->has_barrier_ub(j)) {
            assert(std::isfinite(xu[j]));
            assert(xu[j] != 0.0);
            mufull += (xu[j]+maxp*dxu[j]) * (zu[j]+maxd*dzu[j]);
            num_finite++;
        }
    }
    assert(std::isfinite(mufull));
    mufull /= num_finite;
    mufull /= gammaa;

    double alphap = 1.0;
    double alphad = 1.0;
    Int blockp = -1;
    Int blockd = -1;
    if (maxp < 1.0) {
        if (step_xl <= step_xu) {
            double buffer;
            blockp = block_xl;
            buffer = mufull / (zl[blockp] + maxd*dzl[blockp]);
            alphap = (xl[blockp]-buffer) / (-dxl[blockp]);
        } else {
            double buffer;
            blockp = block_xu;
            buffer = mufull / (zu[blockp] + maxd*dzu[blockp]);
            alphap = (xu[blockp]-buffer) / (-dxu[blockp]);
        }
        alphap = std::max(alphap, gammaf*maxp);
        alphap = std::min(alphap, 1.0);
        assert(blockp >= 0.0);
    }
    if (maxd < 1.0) {
        if (step_zl <= step_zu) {
            double buffer;
            blockd = block_zl;
            buffer = mufull / (xl[blockd] + maxp*dxl[blockd]);
            alphad = (zl[blockd]-buffer) / (-dzl[blockd]);
        } else {
            double buffer;
            blockd = block_zu;
            buffer = mufull / (xu[blockd] + maxp*dxu[blockd]);
            alphad = (zu[blockd]-buffer) / (-dzu[blockd]);
        }
        alphad = std::max(alphad, gammaf*maxd);
        alphad = std::min(alphad, 1.0);
        assert(blockd >= 0.0);
    }
    step_primal_ = std::min(alphap, 1.0-1e-6);
    step_dual_   = std::min(alphad, 1.0-1e-6);

    if (isCentring){
        // When computing stepsizes for a centring step, reduce them
        // by centringAlphaScaling. This ensures that the point is 
        // well centred and does not get too close to the boundary.
        step_primal_ = alphap * control_.centringAlphaScaling();
        step_dual_ = alphad * control_.centringAlphaScaling();
    }
}

void IPM::MakeStep(const Step& step, bool isCentring) {
    StepSizes(step, isCentring);

    const double presidual = iterate_->presidual();
    const double dresidual = iterate_->dresidual();
    const double complementarity = iterate_->complementarity();
    
    iterate_->Update(step_primal_, &step.x[0], &step.xl[0], &step.xu[0],
                     step_dual_,   &step.y[0], &step.zl[0], &step.zu[0]);
    if (!isCentring){
        // A step of length alpha should reduce the residuals by the factor
        // (1-alpha). If a residual that is not yet small hardly decreases
        // although the step is long, the Newton systems are inconsistent,
        // as for a primal or dual infeasible LP. The complementarity gap may
        // still decrease, so the iteration must be counted as bad for
        // Driver() to detect infeasibility. Similarly, if the
        // complementarity gap doesn't decrease although the step is long,
        // the iterate has stagnated (for example when rounding errors in
        // the objective prevent the optimality test from being satisfied).
        const double large_step_thresh = 0.9;
        const double reduction_thresh = 0.9;

        const Model& model = iterate_->model();
        const double tol = iterate_->feasibility_tol();
        const bool presidual_stalled =
            step_primal_ >= large_step_thresh 
            && iterate_->presidual() > reduction_thresh * presidual 
            && iterate_->presidual() > tol * (1.0+model.norm_bounds());
        const bool dresidual_stalled =
            step_dual_ >= large_step_thresh 
            && iterate_->dresidual() > reduction_thresh * dresidual 
            && iterate_->dresidual() > tol * (1.0+model.norm_c());
        const bool complementarity_stalled =
            std::min(step_primal_, step_dual_) >= large_step_thresh 
            && iterate_->complementarity() > reduction_thresh * complementarity;

        if (std::min(step_primal_, step_dual_) < 0.05 ||
            presidual_stalled || dresidual_stalled || complementarity_stalled)
            num_bad_iter_++;
        else
            num_bad_iter_ = 0;
        best_complementarity_ =
            std::min(best_complementarity_, iterate_->complementarity());
    }
}

void IPM::SolveNewtonSystem(const double* rb, const double* rc,
                                    const double* rl, const double* ru,
                                    const double* sl, const double* su,
                                    Step& step) {
    const Model& model = iterate_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Vector& xl = iterate_->xl();
    const Vector& xu = iterate_->xu();
    const Vector& zl = iterate_->zl();
    const Vector& zu = iterate_->zu();
    Vector& dx = step.x;
    Vector& dxl = step.xl;
    Vector& dxu = step.xu;
    Vector& dy = step.y;
    Vector& dzl = step.zl;
    Vector& dzu = step.zu;

    // Build RHS for KKT system.
    Vector rhs1(n+m), rhs2(m);
    if (rc) {
        for (Int j = 0; j < n+m; j++)
            rhs1[j] = -rc[j];
    }
    for (Int j = 0; j < n+m; j++) {
        double rlj = rl ? rl[j] : 0.0;
        double ruj = ru ? ru[j] : 0.0;
        if (iterate_->has_barrier_lb(j))
            rhs1[j] += (sl[j] + zl[j]*rlj) / xl[j];
        if (iterate_->has_barrier_ub(j))
            rhs1[j] -= (su[j] - zu[j]*ruj) / xu[j];
        if (iterate_->StateOf(j) == Iterate::State::fixed)
            rhs1[j] = 0.0;
    }
    assert(AllFinite(rhs1));
    if (rb)
        std::copy(rb, rb+m, std::begin(rhs2));

    // Solve KKT system.
    double tol =  control_.kkt_tol() * std::sqrt(iterate_->mu());
    kkt_->Solve(rhs1, rhs2, tol, dx, dy, info_);
    if (info_->errflag)
        return;

    // Recover solution to Newton system.
    dy *= -1.0;
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            dxl[j] = 0.0;
            dzl[j] = 0.0;
        } else if (iterate_->StateOf(j) == Iterate::State::free) {
            assert(!rl || rl[j] == 0.0);
            dxl[j] = 0.0;
            dzl[j] = 0.0;
        } else {
            double rlj = rl ? rl[j] : 0.0;
            dxl[j] = dx[j] - rlj;
            dzl[j] = (sl[j] - zl[j]*dxl[j]) / xl[j];
        }
    }
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            dxu[j] = 0.0;
            dzu[j] = 0.0;
        } else if (iterate_->StateOf(j) == Iterate::State::free) {
            assert(!ru || ru[j] == 0.0);
            dxu[j] = 0.0;
            dzu[j] = 0.0;
        } else {
            double ruj = ru ? ru[j] : 0.0;
            dxu[j] = ruj - dx[j];
            dzu[j] = (su[j] - zu[j]*dxu[j]) / xu[j];
        }
    }
    assert(AllFinite(dxl));
    assert(AllFinite(dxu));
    assert(AllFinite(dzl));
    assert(AllFinite(dzu));

    // Shift residual to the last two block equations.
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::barrier) {
            assert(std::isfinite(xl[j]) || std::isfinite(xu[j]));
            double atdy = DotColumn(AI, j, dy);
            double rcj = rc ? rc[j] : 0.0;
            if (std::isfinite(xl[j]) && std::isfinite(xu[j])) {
                if (zl[j]*xu[j] >= zu[j]*xl[j])
                    dzl[j] = rcj + dzu[j] - atdy;
                else
                    dzu[j] = -rcj + dzl[j] + atdy;
            } else if (std::isfinite(xl[j])) {
                dzl[j] = rcj + dzu[j] - atdy;
            } else {
                dzu[j] = -rcj + dzl[j] + atdy;
            }
        }
    }

    #ifndef NDEBUG
    // Check solution for free and fixed variables.
    for (Int j = 0; j < n+m; j++) {
        if (iterate_->StateOf(j) == Iterate::State::fixed) {
            assert(dx[j] == 0.0);
            assert(dxl[j] == 0.0 && dxu[j] == 0.0);
            assert(dzl[j] == 0.0 && dzu[j] == 0.0);
        }
        if (iterate_->StateOf(j) == Iterate::State::free)
            assert(dzl[j] == 0.0 && dzu[j] == 0.0);
    }
    #endif
}

void IPM::PrintHeader() {
  std::stringstream h_logging_stream;
  h_logging_stream.str(std::string());
  h_logging_stream
    << " " << Format("Iter", 4)
    << "  " << Format("primal obj", 15)
    << "  " << Format("dual obj", 15)
    << "  " << Format("pinf", 9)
    << "  " << Format("dinf", 9)
    << "  " << Format("gap", 8);
  if (!control_.timelessLog())
    h_logging_stream << "   " << Format("time", 7);
  control_.hLog(h_logging_stream);
  control_.Debug()
    << "  " << Format("stepsizes", 9)
    << "  " << Format("pivots", 7) << " " << Format("kktiter", 7)
    << "  " << Format("P.fixed", 7) << " " << Format("D.fixed", 7);
  control_.Debug(4) << "  " << Format("svdmin(B)", 9);
  control_.Debug(4) << "  " << Format("density", 8);
  control_.hLog("\n");
}

void IPM::PrintOutput() {
    const bool ipm_optimal = iterate_->feasible() && iterate_->optimal();

    double logging_pobj = iterate_->pobjective_after_postproc();
    double logging_dobj = iterate_->dobjective_after_postproc();
    double logging_presidual = iterate_->presidual();
    double logging_dresidual = iterate_->dresidual();

    // Now logging relative primal and dual infeasibility, and also
    // the relative primal dual objective gap
    logging_presidual /= iterate_->bounds_measure_;
    logging_dresidual /= iterate_->costs_measure_;
    double logging_gap = std::abs(logging_pobj - logging_dobj) /
      (1.0+0.5 *std::fabs(logging_pobj + logging_dobj));

    std::stringstream h_logging_stream;
    h_logging_stream.str(std::string());
    h_logging_stream
      << " "  << Format(info_->iter, 3)
      << (ipm_optimal ? "*" : " ")
      << "  " << Scientific(logging_pobj, 15, 8)
      << "  " << Scientific(logging_dobj, 15, 8)
      << "  " << Scientific(logging_presidual, 9, 2)
      << "  " << Scientific(logging_dresidual, 9, 2)
      << "  " << Scientific(logging_gap, 8, 2);
    //    h_logging_stream << "  " << Scientific(iterate_->mu(), 8, 2);
    if (!control_.timelessLog())
      h_logging_stream << "  " << time(control_.Elapsed());
    control_.hLog(h_logging_stream);
    control_.Debug()
      << "  " << Fixed(step_primal_, 4, 2) << " " << Fixed(step_dual_, 4, 2)
      << "  " << Format(kkt_->basis_changes(), 7)
      << " "  << Format(kkt_->iter(), 7);
    control_.Debug()
      << "  " << Format(info_->dual_dropped, 7)
      << " "  << Format(info_->primal_dropped, 7);

    const Basis* basis = kkt_->basis();
    if (basis) {
        if (control_.Debug(4)) {
            control_.Debug(4) << "  "
                              << Scientific(basis->MinSingularValue(), 9, 2);
            Timer timer;
            double density = basis->DensityInverse();
            info_->time_symb_invert += timer.Elapsed();
            control_.Debug(4) << "  " << Scientific(density, 8, 2);
        }
    } else {
        control_.Debug(4) << "  " << Format("-", 9);
        control_.Debug(4) << "  " << Format("-", 8);
    }
    control_.hLog("\n");
}

}  // namespace ipx
