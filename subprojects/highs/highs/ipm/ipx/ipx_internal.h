#ifndef IPX_INTERNAL_H_
#define IPX_INTERNAL_H_

#include <cstring>
#include <valarray>

#include "ipm/ipx/ipx_config.h"
#include "ipm/ipx/ipx_info.h"
#include "ipm/ipx/ipx_parameters.h"
#include "ipm/ipx/ipx_status.h"

namespace ipx {

using Int = ipxint;

struct Info : public ipx_info {
  Info() { std::memset(this, 0, sizeof(Info)); }
};

struct Parameters : public ipx_parameters {
  Parameters() {
    display = 1;
    logfile = nullptr;
    print_interval = 60.0;
    analyse_basis_data = false;
    time_limit = -1.0;
    dualize = -1;
    scale = 1;
    ipm_maxiter = 300;
    ipm_feasibility_tol = 1e-6;
    ipm_optimality_tol = 1e-8;
    ipm_drop_primal = 1e-9;
    ipm_drop_dual = 1e-9;
    kkt_tol = 0.3;
    crash_basis = 1;
    dependency_tol = 1e-6;
    volume_tol = 2.0;
    rows_per_slice = 10000;
    maxskip_updates = 10;
    lu_kernel = 0;
    lu_pivottol = 0.0625;
    run_crossover = 1;
    start_crossover_tol = 1e-8;
    pfeasibility_tol = 1e-7;
    dfeasibility_tol = 1e-7;
    debug = 0;
    switchiter = -1;
    stop_at_switch = 0;
    update_heuristic = 1;
    maxpasses = -1;
    run_centring = 0;
    max_centring_steps = 5;
    centring_ratio_tolerance = 100.0;
    centring_ratio_reduction = 1.5;
    centring_alpha_scaling = 0.5;
    bad_products_tolerance = 3;
    highs_logging = false;
    timeless_log = false;
    log_options = nullptr;
  }

  Parameters(const ipx_parameters& p) : ipx_parameters(p) {}
};

using Vector = std::valarray<double>;

// A vector is treated sparse if it has no more than kHypersparseThreshold * dim
// nonzeros.
static constexpr double kHypersparseThreshold = 0.1;

// When LU factorization is used for rank detection, columns of the active
// submatrix whose maximum entry is <= kLuDependencyTol are removed immediately
// without choosing a pivot.
static constexpr double kLuDependencyTol = 1e-3;

// A fresh LU factorization is considered unstable if
//   ||b-Bx|| / (||b||+||B||*||x||) > kLuStabilityThreshold,
// where x=B\b is computed from the LU factors, b has components +/- 1 that are
// chosen to make x large, and ||.|| is the 1-norm. An unstable factorization
// triggers tightening of the pivot tolerance and refactorization.
static constexpr double kLuStabilityThreshold = 1e-12;

// A Forrest-Tomlin LU update is declared numerically unstable if the relative
// error in the new diagonal entry of U is larger than kFtDiagErrorTol.
static constexpr double kFtDiagErrorTol = 1e-8;

}  // namespace ipx

#endif // IPX_INTERNAL_H_
