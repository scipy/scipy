// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "kkt_solver.h"
#include "timer.h"

namespace ipx {

void KKTSolver::Factorize(Iterate* pt, Info* info) {
    Timer timer;
    _Factorize(pt, info);
    info->time_kkt_factorize += timer.Elapsed();
}

void KKTSolver::Solve(const Vector& a, const Vector& b, double tol,
                      Vector& x, Vector& y, Info* info) {
    Timer timer;
    _Solve(a, b, tol, x, y, info);
    info->time_kkt_solve += timer.Elapsed();
}

Int KKTSolver::iter() const { return _iter(); }
Int KKTSolver::basis_changes() const { return _basis_changes(); }
const Basis* KKTSolver::basis() const { return _basis(); }

}  // namespace ipx
