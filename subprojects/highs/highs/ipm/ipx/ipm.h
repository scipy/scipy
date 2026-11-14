#ifndef IPX_IPM_H_
#define IPX_IPM_H_

#include "ipm/ipx/control.h"
#include "ipm/ipx/kkt_solver.h"
#include "ipm/ipx/iterate.h"

namespace ipx {

// IPM implements an interior point method based on KKTSolver and Iterate.
// The algorithm is a variant of Mehrotra's [1] predictor-corrector method
// that requires two linear system solves per iteration.
//
// [1] S. Mehrotra, "On the implementation of a primal-dual interior point
//     method", SIAM J. Optim., 2 (1992).

class IPM {
public:
    explicit IPM(const Control& control);

    // Initializes @iterate with a starting point for Driver(). The KKT solver
    // must allow Factorize(NULL, info) (see kkt_solver.h).
    // On return info->status_ipm is
    // IPX_STATUS_not_run    if successful,
    // IPX_STATUS_time_limit if the KKT solver was interrupted by time limit,
    // IPX_STATUS_failed     if the KKT solver failed with info->errflag.
    // If the method did not terminate successfully, @iterate is unchanged.
    void StartingPoint(KKTSolver* kkt, Iterate* iterate, Info* info);

    // Updates @iterate by interior point iterations. On return ipm_status is
    // IPX_STATUS_optimal       if iterate->term_crit_reached() is true,
    // IPX_STATUS_iter_limit    if info->iter >= maxiter(),
    // IPX_STATUS_no_progress   if no progress over a number of iterations,
    // IPX_STATUS_time_limit    if interrupted by time limit,
    // IPX_STATUS_failed        if the KKT solver failed with info->errflag.
    void Driver(KKTSolver* kkt, Iterate* iterate, Info* info);

    Int maxiter() const { return maxiter_; }
    void maxiter(Int i) { maxiter_ = i; }
    void PrintHeader();

private:
    struct Step;
    // IPM terminates when the complementarity gap of the current iterate
    // exceeds kDivergeTol times the smallest complementarity gap of all
    // iterates so far.
    static constexpr double kDivergeTol = 1e6;

    void ComputeStartingPoint();
    void Predictor(Step& step);
    void AddCorrector(Step& step);
    void Centring(Step& step, double mu_to_use);
    void AssessCentrality(const Vector& xl, const Vector& xu, const Vector& zl, 
                          const Vector& zu,double mu, bool print = true);
    bool EvaluateCentringStep(const Step& step, double prev_ratio, Int prev_bad);
    void StepSizes(const Step& step, bool isCentring = false);
    void MakeStep(const Step& step, bool isCentring = false);
    // Reduces the following linear system to KKT form:
    //  [ AI                 ] [dx ]    [rb]
    //  [ I  -I              ] [dxl] =  [rl]
    //  [ I      I           ] [dxu]    [ru]
    //  [          AI'  I -I ] [dy ]    [rc]
    //  [    Zl        Xl    ] [dzl]    [sl]
    //  [       Zu        Xu ] [dzu]    [su]
    // Each of @rb, @rc, @rl and @ru can be NULL, in which case its entries are
    // assumed to be 0.0. This is currently not used, but was implemented for
    // computing centrality correctors.
    void SolveNewtonSystem(const double* rb, const double* rc,
                           const double* rl, const double* ru,
                           const double* sl, const double* su, Step& lhs);
    void PrintOutput();

    const Control& control_;
    KKTSolver* kkt_{nullptr};
    Iterate* iterate_{nullptr};
    Info* info_{nullptr};

    double step_primal_{0.0}, step_dual_{0.0};
    // Counts the # bad iterations since the last good iteration. An iteration
    // is bad if the primal or dual step size is < 0.05.
    Int num_bad_iter_{0};
    // Smallest complementarity gap of all iterates so far.
    double best_complementarity_{0.0};

    Int maxiter_{-1};

    // indicators of centrality for centring steps
    double centring_ratio{0.0};
    Int bad_products{0};
};

}  // namespace ipx

#endif  // IPX_IPM_H_
