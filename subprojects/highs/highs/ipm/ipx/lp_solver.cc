#include "ipm/ipx/lp_solver.h"
#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>
#include "ipm/ipx/crossover.h"
#include "ipm/ipx/info.h"
#include "ipm/ipx/kkt_solver_basis.h"
#include "ipm/ipx/kkt_solver_diag.h"
#include "ipm/ipx/starting_basis.h"
#include "ipm/ipx/utils.h"

namespace ipx {

Int LpSolver::LoadModel(Int num_var, const double offset,
                        const double* obj, const double* lb,
                        const double* ub, Int num_constr, const Int* Ap,
                        const Int* Ai, const double* Ax, const double* rhs,
                        const char* constr_type) {
    ClearModel();
    Int errflag = model_.Load(control_, num_constr, num_var, Ap, Ai, Ax, rhs,
                              constr_type, offset, obj, lb, ub);
    model_.GetInfo(&info_);
    return errflag;
}

Int LpSolver::LoadIPMStartingPoint(const double* x, const double* xl,
                                   const double* xu, const double* slack,
                                   const double* y, const double* zl,
                                   const double* zu) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    x_start_.resize(n+m);
    xl_start_.resize(n+m);
    xu_start_.resize(n+m);
    y_start_.resize(m);
    zl_start_.resize(n+m);
    zu_start_.resize(n+m);
    Int errflag = model_.PresolveIPMStartingPoint(x, xl, xu, slack, y, zl, zu,
                                                  x_start_, xl_start_,
                                                  xu_start_, y_start_,
                                                  zl_start_, zu_start_);
    if (errflag) {
        ClearIPMStartingPoint();
        return errflag;
    }
    MakeIPMStartingPointValid();
    return 0;
}

Int LpSolver::Solve() {
    if (model_.empty())
        return info_.status = IPX_STATUS_no_model;
    ClearSolution();
    control_.ResetTimer();
    control_.OpenLogfile();
    control_.hLog("IPX version 1.0\n");

    InteriorPointSolve();

    const bool run_crossover_on = control_.run_crossover() == 1;
    const bool run_crossover_choose = control_.run_crossover() == -1;
    const bool run_crossover_not_off = run_crossover_choose || run_crossover_on;
    const bool run_crossover =
      (info_.status_ipm == IPX_STATUS_optimal && run_crossover_on) ||
      (info_.status_ipm == IPX_STATUS_imprecise && run_crossover_not_off);
    //        if ((info_.status_ipm == IPX_STATUS_optimal ||
    //             info_.status_ipm == IPX_STATUS_imprecise) && run_crossover_on) {
    if (run_crossover) {
      if (run_crossover_on) {
	control_.hLog("Running crossover as requested\n");
      } else if (run_crossover_choose) {
	assert(info_.status_ipm == IPX_STATUS_imprecise);
	control_.hLog("Running crossover since IPX is imprecise\n");
      } else {
	assert(run_crossover_on || run_crossover_choose);
      }
      BuildCrossoverStartingPoint();
      RunCrossover();
    }
    if (basis_) {
      info_.ftran_sparse = basis_->frac_ftran_sparse();
      info_.btran_sparse = basis_->frac_btran_sparse();
      info_.time_lu_invert = basis_->time_factorize();
      info_.time_lu_update = basis_->time_update();
      info_.time_ftran = basis_->time_ftran();
      info_.time_btran = basis_->time_btran();
      info_.mean_fill = basis_->mean_fill();
      info_.max_fill = basis_->max_fill();
    }
    if (info_.status_ipm == IPX_STATUS_primal_infeas ||
	info_.status_ipm == IPX_STATUS_dual_infeas ||
	info_.status_crossover == IPX_STATUS_primal_infeas ||
	info_.status_crossover == IPX_STATUS_dual_infeas) {
      // When IPM or crossover detect the model to be infeasible
      // (currently only the former is implemented), then the problem is
      // solved.
      info_.status = IPX_STATUS_solved;
    } else {
      Int method_status = run_crossover ?
	info_.status_crossover : info_.status_ipm;
      if (method_status == IPX_STATUS_optimal ||
	  method_status == IPX_STATUS_imprecise)
	info_.status = IPX_STATUS_solved;
      else
	info_.status = IPX_STATUS_stopped;
    }

    PrintSummary();

    info_.time_total = control_.Elapsed();
    control_.Debug(2) << info_;
    control_.CloseLogfile();
    if (control_.reportBasisData())
      basis_->reportBasisData();
    return info_.status;
}

Info LpSolver::GetInfo() const {
    return info_;
}

Int LpSolver::GetInteriorSolution(double* x, double* xl, double* xu,
                                  double* slack, double* y, double* zl,
                                  double* zu) const {
    if (!iterate_)
        return -1;
    model_.PostsolveInteriorSolution(
        iterate_->x(), iterate_->xl(), iterate_->xu(),
        iterate_->y(), iterate_->zl(), iterate_->zu(),
        x, xl, xu, slack, y, zl, zu);
    return 0;
}

Int LpSolver::GetBasicSolution(double* x, double* slack, double* y, double* z,
                               Int* cbasis, Int* vbasis) const {
    if (basic_statuses_.empty())
        return -1;
    model_.PostsolveBasicSolution(x_crossover_, y_crossover_, z_crossover_,
                                  basic_statuses_, x, slack, y, z);
    model_.PostsolveBasis(basic_statuses_, cbasis, vbasis);
    return 0;
}

Parameters LpSolver::GetParameters() const {
    return control_.parameters();
}

void LpSolver::SetParameters(Parameters new_parameters) {
    control_.parameters(new_parameters);
}

void LpSolver::SetCallback(HighsCallback* callback) {
    control_.callback(callback);
}

void LpSolver::ClearModel() {
    model_.clear();
    ClearSolution();
    ClearIPMStartingPoint();
}

void LpSolver::ClearIPMStartingPoint() {
    x_start_.resize(0);
    xl_start_.resize(0);
    xu_start_.resize(0);
    y_start_.resize(0);
    zl_start_.resize(0);
    zu_start_.resize(0);
}

Int LpSolver::CrossoverFromStartingPoint(const double* x_start,
                                         const double* slack_start,
                                         const double* y_start,
                                         const double* z_start) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    const SparseMatrix& AI = model_.AI();

    ClearSolution();
    control_.hLog("Crossover from starting point\n");

    x_crossover_.resize(n+m);
    y_crossover_.resize(m);
    z_crossover_.resize(n+m);
    crossover_weights_.resize(0);
    model_.PresolveStartingPoint(x_start, slack_start, y_start, z_start,
                                 x_crossover_, y_crossover_, z_crossover_);

    // Check that starting point is complementary and satisfies bound and sign
    // conditions.
    for (Int j = 0; j < n+m; j++) {
        if (x_crossover_[j] < lb[j] || x_crossover_[j] > ub[j])
            return IPX_ERROR_invalid_vector;
        if (x_crossover_[j] != lb[j] && z_crossover_[j] > 0.0)
            return IPX_ERROR_invalid_vector;
        if (x_crossover_[j] != ub[j] && z_crossover_[j] < 0.0)
            return IPX_ERROR_invalid_vector;
    }

    // Construct starting basis.
    basis_.reset(new Basis(control_, model_));
    if (control_.crash_basis()) {
        // Take columns in the following order of priority:
        // - free columns
        // - columns between their bounds, in increasing number of nonzeros
        // - columns with zero dual, in increasing number of nonzeros
        // - Fixed columns and those with nonzero dual
        Timer timer;
        Vector colweight(n+m);
        for (Int j = 0; j < n+m; j++) {
            Int nz = AI.entries(j);
            if (lb[j] == ub[j])
                colweight[j] = 0.0;
            else if (std::isinf(lb[j]) && std::isinf(ub[j]))
                colweight[j] = INFINITY;
            else if (z_crossover_[j] != 0.0)
                colweight[j] = 0.0;
            else if (x_crossover_[j] != lb[j] && x_crossover_[j] != ub[j])
                colweight[j] = m + (m-nz+1);
            else
                colweight[j] = m-nz+1;
        }
        basis_->ConstructBasisFromWeights(&colweight[0], &info_);
        info_.time_starting_basis += timer.Elapsed();
        if (info_.errflag) {
            ClearSolution();
            return 0;
        }
    }

    RunCrossover();
    return 0;
}

Int LpSolver::GetIterate(double* x, double* y, double* zl, double* zu,
                         double* xl, double* xu) {
    if (!iterate_)
        return -1;
    if (x)
        std::copy(std::begin(iterate_->x()), std::end(iterate_->x()), x);
    if (y)
        std::copy(std::begin(iterate_->y()), std::end(iterate_->y()), y);
    if (zl)
        std::copy(std::begin(iterate_->zl()), std::end(iterate_->zl()), zl);
    if (zu)
        std::copy(std::begin(iterate_->zu()), std::end(iterate_->zu()), zu);
    if (xl)
        std::copy(std::begin(iterate_->xl()), std::end(iterate_->xl()), xl);
    if (xu)
        std::copy(std::begin(iterate_->xu()), std::end(iterate_->xu()), xu);
    return 0;
}

// Returns a vector of basic statuses that is consistent with the basis and
// the bounds from the model.
static std::vector<Int> BuildBasicStatuses(const Basis& basis) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    std::vector<Int> basic_statuses(n+m);
    for (Int j = 0; j < n+m; j++) {
        if (basis.IsBasic(j)) {
            basic_statuses[j] = IPX_basic;
        } else if (std::isfinite(lb[j])) {
            basic_statuses[j] = IPX_nonbasic_lb;
        } else if (std::isfinite(ub[j])) {
            basic_statuses[j] = IPX_nonbasic_ub;
        } else {
            basic_statuses[j] = IPX_superbasic;
        }
    }
    return basic_statuses;
}

Int LpSolver::GetBasis(Int* cbasis, Int* vbasis) {
    if (!basis_)
        return -1;
    if (!basic_statuses_.empty()) {
        // crossover provides basic statuses
        model_.PostsolveBasis(basic_statuses_, cbasis, vbasis);
    } else {
        model_.PostsolveBasis(BuildBasicStatuses(*basis_), cbasis, vbasis);
    }
    return 0;
}

Int LpSolver::GetKKTMatrix(Int* AIp, Int* AIi, double* AIx, double* g) {
    if (!iterate_)
        return -1;
    if (AIp && AIi && AIx) {
        const SparseMatrix& AI = model_.AI();
        std::copy_n(AI.colptr(), AI.cols()+1, AIp);
        std::copy_n(AI.rowidx(), AI.entries(), AIi);
        std::copy_n(AI.values(), AI.entries(), AIx);
    }
    if (g) {
        Int m = model_.rows();
        Int n = model_.cols();
        for (Int j = 0; j < n+m; j++) {
            switch (iterate_->StateOf(j)) {
            case Iterate::State::fixed:
                g[j] = INFINITY;
                break;
            case Iterate::State::free:
                g[j] = 0.0;
                break;
            case Iterate::State::barrier:
                g[j] = iterate_->zl(j)/iterate_->xl(j) +
                    iterate_->zu(j)/iterate_->xu(j);
                assert(std::isfinite(g[j]));
                assert(g[j] > 0.0);
                break;
            default:
                assert(0);
            }
        }
    }
    return 0;
}

Int LpSolver::SymbolicInvert(Int* rowcounts, Int* colcounts) {
    if (!basis_)
        return -1;
    basis_->SymbolicInvert(rowcounts, colcounts);
    return 0;
}

void LpSolver::ClearSolution() {
    iterate_.reset(nullptr);
    basis_.reset(nullptr);
    x_crossover_.resize(0);
    y_crossover_.resize(0);
    z_crossover_.resize(0);
    crossover_weights_.resize(0);
    basic_statuses_.clear();
    basic_statuses_.shrink_to_fit();
    info_ = Info();
    // Restore info entries that belong to model.
    model_.GetInfo(&info_);

}

void LpSolver::InteriorPointSolve() {
    if (control_.runCentring()) {
      control_.hLog("Interior point solve for analytic centre\n");
    } else {
      control_.hLog("Interior point solve\n");
    }

    // Allocate new iterate and set tolerances for IPM termination test.
    iterate_.reset(new Iterate(model_));
    iterate_->feasibility_tol(control_.ipm_feasibility_tol());
    iterate_->optimality_tol(control_.ipm_optimality_tol());
    if (control_.run_crossover())
        iterate_->start_crossover_tol(control_.start_crossover_tol());

    RunIPM();

    iterate_->Postprocess();
    iterate_->EvaluatePostsolved(&info_);

    // Declare status_ipm "imprecise" if the IPM terminated optimal but the
    // solution after postprocessing/postsolve does not satisfy tolerances.
    if (info_.status_ipm == IPX_STATUS_optimal) {
        if (std::abs(info_.rel_objgap) > control_.ipm_optimality_tol() ||
            info_.rel_presidual > control_.ipm_feasibility_tol() ||
            info_.rel_dresidual > control_.ipm_feasibility_tol())
            info_.status_ipm = IPX_STATUS_imprecise;
    }

    // Assess the success of analytic centre calculation
    if (info_.centring_tried)
      info_.status_ipm = info_.centring_success ? IPX_STATUS_optimal : IPX_STATUS_imprecise;
}

void LpSolver::RunIPM() {
    IPM ipm(control_);
    info_.centring_tried = false;
    info_.centring_success = false;

    if (x_start_.size() != 0) {
        control_.hLog(" Using starting point provided by user. Skipping initial iterations.\n");
        iterate_->Initialize(x_start_, xl_start_, xu_start_,
                             y_start_, zl_start_, zu_start_);
    }
    else {
        ComputeStartingPoint(ipm);
        if (info_.status_ipm != IPX_STATUS_not_run)
            return;
        RunInitialIPM(ipm);
        if (info_.status_ipm != IPX_STATUS_not_run)
            return;
    }
    BuildStartingBasis();
    if (info_.status_ipm != IPX_STATUS_not_run || 
        info_.centring_tried)
        return;
    RunMainIPM(ipm);
}

void LpSolver::MakeIPMStartingPointValid() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    Vector& xl = xl_start_;
    Vector& xu = xu_start_;
    Vector& zl = zl_start_;
    Vector& zu = zu_start_;

    Int numComplementarityProducts = 0;
    double sumComplementarityProducts = 0.0;
    for (Int j = 0; j < n+m; ++j) {
        if (xl[j] > 0.0 && zl[j] > 0.0) {
            sumComplementarityProducts += xl[j] * zl[j];
            numComplementarityProducts++;
        }
        if (xu[j] > 0.0 && zu[j] > 0.0) {
            sumComplementarityProducts += xu[j] * zu[j];
            numComplementarityProducts++;
        }
    }
    const double mu = numComplementarityProducts ?
        sumComplementarityProducts / numComplementarityProducts : 1.0;

    for (Int j = 0; j < n+m; ++j) {
        if (std::isfinite(lb[j])) {
            assert(std::isfinite(xl[j]) && xl[j] >= 0.0);
            assert(std::isfinite(zl[j]) && zl[j] >= 0.0);
            if (xl[j] == 0.0 && zl[j] == 0.0)
                xl[j] = zl[j] = std::sqrt(mu);
            else if (xl[j] == 0.0)
                xl[j] = mu / zl[j];
            else if (zl[j] == 0.0)
                zl[j] = mu / xl[j];
        } else {
            assert(xl[j] == INFINITY);
            assert(zl[j] == 0.0);
        }
        if (std::isfinite(ub[j])) {
            assert(std::isfinite(xu[j]) && xu[j] >= 0.0);
            assert(std::isfinite(zu[j]) && zu[j] >= 0.0);
            if (xu[j] == 0.0 && zu[j] == 0.0)
                xu[j] = zu[j] = std::sqrt(mu);
            else if (xu[j] == 0.0)
                xu[j] = mu / zu[j];
            else if (zu[j] == 0.0)
                zu[j] = mu / xu[j];
        } else {
            assert(xu[j] == INFINITY);
            assert(zu[j] == 0.0);
        }
    }
}

void LpSolver::ComputeStartingPoint(IPM& ipm) {
    Timer timer;
    KKTSolverDiag kkt(control_, model_);

    // If the starting point procedure fails, then iterate_ remains as
    // initialized by the constructor, which is a valid state for
    // postprocessing/postsolving.
    ipm.StartingPoint(&kkt, iterate_.get(), &info_);
    info_.time_ipm1 += timer.Elapsed();
}

void LpSolver::RunInitialIPM(IPM& ipm) {
    Timer timer;
    KKTSolverDiag kkt(control_, model_);

    Int switchiter = control_.switchiter();
    if (switchiter < 0) {
        // Switch iteration not specified by user. Run as long as KKT solver
        // converges within min(500,10+m/20) iterations.
        Int m = model_.rows();
        kkt.maxiter(std::min(500l, (long) (10+m/20) ));
        ipm.maxiter(control_.ipm_maxiter());
    } else {
        ipm.maxiter(std::min(switchiter, control_.ipm_maxiter()));
    }
    ipm.Driver(&kkt, iterate_.get(), &info_);
    switch (info_.status_ipm) {
    case IPX_STATUS_optimal:
        // If the IPM reached its termination criterion in the initial
        // iterations (happens rarely), we still call the IPM again with basis
        // preconditioning. In exact arithmetic it would terminate without an
        // additional iteration. A starting basis is then available for
        // crossover.
        info_.status_ipm = IPX_STATUS_not_run;
        break;
    case IPX_STATUS_no_progress:
        info_.status_ipm = IPX_STATUS_not_run;
        break;
    case IPX_STATUS_failed:
        info_.status_ipm = IPX_STATUS_not_run;
        info_.errflag = 0;
        break;
    case IPX_STATUS_iter_limit:
        if (info_.iter < control_.ipm_maxiter()) // stopped at switchiter
            info_.status_ipm = IPX_STATUS_not_run;
    }
    info_.time_ipm1 += timer.Elapsed();
}

void LpSolver::BuildStartingBasis() {
    if (control_.stop_at_switch() < 0) {
        info_.status_ipm = IPX_STATUS_debug;
        return;
    }
    control_.hLog(" Constructing starting basis...\n");
    basis_.reset(new Basis(control_, model_));
    StartingBasis(iterate_.get(), basis_.get(), &info_);
    if (info_.errflag == IPX_ERROR_user_interrupt) {
        info_.errflag = 0;
        info_.status_ipm = IPX_STATUS_user_interrupt;
        return;
    } else if (info_.errflag == IPX_ERROR_time_interrupt) {
        info_.errflag = 0;
        info_.status_ipm = IPX_STATUS_time_limit;
        return;
    } else if (info_.errflag) {
        info_.status_ipm = IPX_STATUS_failed;
        return;
    }
    if (model_.dualized()) {
        std::swap(info_.dependent_rows, info_.dependent_cols);
        std::swap(info_.rows_inconsistent, info_.cols_inconsistent);
    }
    if (control_.stop_at_switch() > 0) {
        info_.status_ipm = IPX_STATUS_debug;
        return;
    }
    if (info_.rows_inconsistent) {
        info_.status_ipm = IPX_STATUS_primal_infeas;
        return;
    }
    if (info_.cols_inconsistent) {
        info_.status_ipm = IPX_STATUS_dual_infeas;
        return;
    }
}

void LpSolver::RunMainIPM(IPM& ipm) {
    KKTSolverBasis kkt(control_, *basis_);
    Timer timer;
    ipm.PrintHeader();
    ipm.maxiter(control_.ipm_maxiter());
    ipm.Driver(&kkt, iterate_.get(), &info_);
    info_.time_ipm2 = timer.Elapsed();
}

void LpSolver::BuildCrossoverStartingPoint() {
    const Int m = model_.rows();
    const Int n = model_.cols();

    // Construct a complementary primal-dual point from the final IPM iterate.
    // This usually increases the residuals to Ax=b and A'y+z=c.
    x_crossover_.resize(n+m);
    y_crossover_.resize(m);
    z_crossover_.resize(n+m);
    iterate_->DropToComplementarity(x_crossover_, y_crossover_, z_crossover_);

    // Perform dual pushes in increasing order and primal pushes in decreasing
    // order of the scaling factors from the final IPM iterate.
    crossover_weights_.resize(n+m);
    for (Int j = 0; j < n+m; j++)
        crossover_weights_[j] = iterate_->ScalingFactor(j);
}

void LpSolver::RunCrossover() {
    assert(basis_);
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();
    const double *weights = NULL;
    basic_statuses_.clear();

    if (crossover_weights_.size()) {
        assert((Int)crossover_weights_.size() == n+m);
        weights = &crossover_weights_[0];
    }

    Crossover crossover(control_);
    crossover.PushAll(basis_.get(), x_crossover_, y_crossover_, z_crossover_,
                      weights, &info_);
    info_.time_crossover =
        crossover.time_primal() + crossover.time_dual();
    info_.updates_crossover =
        crossover.primal_pivots() + crossover.dual_pivots();
    if (info_.status_crossover != IPX_STATUS_optimal) {
        // Crossover failed. Discard solution.
        x_crossover_.resize(0);
        y_crossover_.resize(0);
        z_crossover_.resize(0);
        return;
    }

    // Recompute vertex solution and set basic statuses.
    basis_->ComputeBasicSolution(x_crossover_, y_crossover_, z_crossover_);
    basic_statuses_.resize(n+m);
    for (size_t j = 0; j < basic_statuses_.size(); j++) {
        if (basis_->IsBasic(j)) {
            basic_statuses_[j] = IPX_basic;
        } else {
            if (lb[j] == ub[j])
                basic_statuses_[j] = z_crossover_[j] >= 0.0 ?
                    IPX_nonbasic_lb : IPX_nonbasic_ub;
            else if (x_crossover_[j] == lb[j])
                basic_statuses_[j] = IPX_nonbasic_lb;
            else if (x_crossover_[j] == ub[j])
                basic_statuses_[j] = IPX_nonbasic_ub;
            else
                basic_statuses_[j] = IPX_superbasic;
        }
    }
    control_.Debug()
        << Textline("Bound violation of basic solution:")
        << sci2(PrimalInfeasibility(model_, x_crossover_)) << '\n'
        << Textline("Dual sign violation of basic solution:")
        << sci2(DualInfeasibility(model_, x_crossover_, z_crossover_)) << '\n';
    control_.Debug()
        << Textline("Minimum singular value of basis matrix:")
        << sci2(basis_->MinSingularValue()) << '\n';

    // Declare crossover status "imprecise" if the vertex solution defined by
    // the final basis does not satisfy tolerances.
    model_.EvaluateBasicSolution(x_crossover_, y_crossover_, z_crossover_,
                                 basic_statuses_, &info_);
    if (info_.primal_infeas > control_.pfeasibility_tol() ||
        info_.dual_infeas > control_.dfeasibility_tol())
        info_.status_crossover = IPX_STATUS_imprecise;
}

void LpSolver::PrintSummary() {
  std::stringstream h_logging_stream;
  h_logging_stream.str(std::string());
  h_logging_stream << "Summary\n";
  if (!control_.timelessLog())
    h_logging_stream << Textline("Runtime:") << fix2(control_.Elapsed()) << "s\n";
  h_logging_stream << Textline("Status interior point solve:")
                   << StatusString(info_.status_ipm) << '\n'
                   << Textline("Status crossover:")
                   << StatusString(info_.status_crossover) << '\n';
  control_.hLog(h_logging_stream);
  if (info_.status_ipm == IPX_STATUS_optimal ||
      info_.status_ipm == IPX_STATUS_imprecise) {
    h_logging_stream
      << Textline("objective value:") << sci8(info_.pobjval) << '\n'
      << Textline("interior solution primal residual (abs/rel):")
      << sci2(info_.abs_presidual) << " / " << sci2(info_.rel_presidual)
      << '\n'
      << Textline("interior solution dual residual (abs/rel):")
      << sci2(info_.abs_dresidual) << " / " << sci2(info_.rel_dresidual)
      << '\n'
      << Textline("interior solution objective gap (abs/rel):")
      << sci2(info_.pobjval-info_.dobjval) << " / "
      << sci2(info_.rel_objgap)  << '\n';
    control_.hLog(h_logging_stream);
  }
  if (info_.status_crossover == IPX_STATUS_optimal ||
      info_.status_crossover == IPX_STATUS_imprecise) {
    h_logging_stream
      << Textline("basic solution primal infeasibility:")
      << sci2(info_.primal_infeas) << '\n'
      << Textline("basic solution dual infeasibility:")
      << sci2(info_.dual_infeas) << '\n';
    control_.hLog(h_logging_stream);
  }
}

}  // namespace ipx
