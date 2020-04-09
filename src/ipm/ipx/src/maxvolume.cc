// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "maxvolume.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include "timer.h"
#include "utils.h"

namespace ipx {

Maxvolume::Maxvolume(const Control& control) : control_(control) {}

Int Maxvolume::RunSequential(const double* colscale, Basis& basis) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    IndexedVector ftran(m);
    Timer timer;
    Int errflag = 0;

    const Int maxpasses = control_.maxpasses();
    const double volumetol = std::max(control_.volume_tol(), 1.0);

    // Maintain a copy of the inverse scaling factors of basic variables.
    // If a variable has status BASIC_FREE, its inverse scaling factor is
    // assumed zero, regardless of the entry in colscale. This ensures that
    // BASIC_FREE variables are never pivoted out of the basis.
    Vector invscale_basic(m);
    for (Int p = 0; p < m; p++) {
        Int j = basis[p];
        if (basis.StatusOf(j) == Basis::BASIC)
            invscale_basic[p] = colscale ? 1.0/colscale[j] : 1.0;
    }

    Reset();
    while (passes_ < maxpasses || maxpasses < 0) {
        tblnnz_ = 0;
        tblmax_ = 0.0;
        frobnorm_squared_ = 0.0;
        Int updates_last = 0;   // # basis updates in this pass
        std::vector<Int> candidates = Sortperm(n+m, colscale, false);
        while(!candidates.empty()) {
            Int j = candidates.back();
            const double dj = colscale ? colscale[j] : 1.0;
            if (dj == 0.0)      // all remaining columns have scaling factor 0
                break;
            if (basis.StatusOf(j) != Basis::NONBASIC) {
                candidates.pop_back();
                continue;
            }
            if ((errflag = control_.InterruptCheck()) != 0)
                break;
            basis.SolveForUpdate(j, ftran);
            Int pmax = -1;
            double vmax = 0.0;
            auto search_pivot = [&](Int p, double x) {
                double v = std::abs(x) * invscale_basic[p] * dj;
                if (v > vmax) {
                    vmax = v;
                    pmax = p;
                }
                tblnnz_ += v != 0;
                frobnorm_squared_ += v*v;
            };
            for_each_nonzero(ftran, search_pivot);
            tblmax_ = std::max(tblmax_, vmax);
            if (vmax <= volumetol) {
                skipped_++;
                candidates.pop_back();
                continue;
            }

            const Int jb = basis[pmax];
            assert(basis.StatusOf(jb) == Basis::BASIC);
            bool exchanged;
            errflag = basis.ExchangeIfStable(jb, j, ftran[pmax], -1, &exchanged);
            if (errflag)
                break;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            invscale_basic[pmax] = 1.0 / dj;
            updates_last++;
            volinc_ += std::log2(vmax);
            candidates.pop_back();
        }
        updates_ += updates_last;
        passes_++;
        if (updates_last == 0 || errflag != 0)
            break;
    }
    time_ = timer.Elapsed();
    return errflag;
}

struct Maxvolume::Slice {
    Slice(Int m, Int n) :
        colscale(n+m), invscale_basic(m), tblrow_used(m), colweights(n+m),
        lhs(m), row(n+m), work(m) {}
    Vector colscale;
    Vector invscale_basic;
    std::vector<bool> tblrow_used;
    Vector colweights;
    IndexedVector lhs, row;
    Vector work;
};

Int Maxvolume::RunHeuristic(const double* colscale, Basis& basis) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    Slice slice(m, n);
    Int errflag = 0;
    Timer timer;

    Reset();
    Int num_slices = 5 + std::max((long)(m/control_.rows_per_slice()), 0l);
    num_slices = std::min(num_slices, m);

    // Maintain a copy of the inverse scaling factors of basic variables.
    for (Int p = 0; p < m; p++) {
        Int j = basis[p];
        if (basis.StatusOf(j) == Basis::BASIC) {
            slice.invscale_basic[p] = colscale ? 1.0/colscale[j] : 1.0;
            assert(std::isfinite(slice.invscale_basic[p]));
        }
    }

    // We need to work with a copy of the column scaling factors, in which
    // skipped columns are set to zero (and not scanned again).
    for (Int j = 0; j < n+m; j++) {
        if (basis.StatusOf(j) == Basis::NONBASIC)
            slice.colscale[j] = colscale ? colscale[j] : 1.0;
    }

    // Split tableau matrix into num_slices row slices. In each call to
    // Driver() there is exactly one row from each slice that has
    // tblrow_used[p] == 1, while the other rows in the slice have
    // tblrow_used[p] == 0.
    std::vector<Int> perm = Sortperm(m, &slice.invscale_basic[0], false);
    for (Int s = 0; s < num_slices; s++) {
        for (Int i = 0; i < m; i++)
            slice.tblrow_used[perm[i]] = i%num_slices == s;
        errflag = Driver(basis, slice);
        if (errflag)
            break;
    }

    time_ = timer.Elapsed();
    passes_ = -1;
    slices_ = num_slices;
    return errflag;
}

Int Maxvolume::updates() const { return updates_; }
Int Maxvolume::skipped() const { return skipped_; }
Int Maxvolume::passes() const { return passes_; }
Int Maxvolume::slices() const { return slices_; }
double Maxvolume::volinc() const { return volinc_; }
double Maxvolume::time() const { return time_; }
Int Maxvolume::tblnnz() const { return tblnnz_; }
double Maxvolume::tblmax() const { return tblmax_; }
double Maxvolume::frobnorm_squared() const { return frobnorm_squared_; }

void Maxvolume::Reset() {
    updates_ = 0;
    skipped_ = 0;
    passes_ = 0;
    slices_ = 0;
    volinc_ = 0.0;;
    time_ = 0.0;
    tblnnz_ = 0;
    tblmax_ = 0.0;
    frobnorm_squared_ = 0.0;
}

// Returns a vector holding the indices of the largest and second largest
// entry in @weights. The index of the largest entry is at the back.
static std::vector<Int> FindLargest(const Vector& weights) {
    const Int n = weights.size();
    Int jmax = 0;               // index of largest element
    Int jmax2 = 0;              // index of second largest element
    double wmax = 0.0;
    double wmax2 = 0.0;

    for (Int j = 0; j < n; j++) {
        double w = std::abs(weights[j]);
        if (w > wmax) {
            wmax2 = wmax;
            wmax = w;
            jmax2 = jmax;
            jmax = j;
        }
        else if (w > wmax2) {
            wmax2 = w;
            jmax2 = j;
        }
    }
    return {jmax2, jmax};
}

Int Maxvolume::Driver(Basis& basis, Slice& slice) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    Int errflag = 0;

    const double volumetol = std::max(control_.volume_tol(), 1.0);
    const Int maxskip = control_.maxskip_updates();

    Vector& colscale = slice.colscale;
    Vector& invscale_basic = slice.invscale_basic;
    const std::vector<bool>& tblrow_used = slice.tblrow_used;
    Vector& colweights = slice.colweights;
    IndexedVector& lhs = slice.lhs;
    IndexedVector& row = slice.row;
    Vector& work = slice.work;

    // Compute column weights.
    for (Int p = 0; p < m; p++)
        work[p] = tblrow_used[p] ? invscale_basic[p] : 0.0;
    basis.SolveDense(work, work, 'T');
    for (Int j = 0; j < n+m; j++) {
        if (colscale[j]) {
            assert(basis.StatusOf(j) == Basis::NONBASIC);
            double sum = DotColumn(AI, j, work);
            colweights[j] = sum * colscale[j];
        }
        else
            colweights[j] = 0.0;
    }

    std::vector<Int> candidates;
    Int skipped = 0;
    while (true) {
        // Pick column with maximum weight.
        if (candidates.empty()) {
            candidates = FindLargest(colweights);
            assert(!candidates.empty());
        }
        Int jn = candidates.back();
        assert(jn >= 0 && jn < n+m);
        const double weight = colweights[jn];
        if (weight == 0.0)
            break;
        assert(basis.StatusOf(jn) == Basis::NONBASIC);
        assert(colscale[jn] > 0.0);

        if ((errflag = control_.InterruptCheck()) != 0)
            break;

        // Find maximum scaled FTRAN entry.
        basis.SolveForUpdate(jn, lhs);
        const Int pmax = ScaleFtran(colscale[jn], invscale_basic, lhs);
        const double scaled_pivot = lhs[pmax];
        const double vmax = std::abs(scaled_pivot);

        // Skip the column if exchange does not increase volume enough.
        if (vmax <= volumetol) {
            colweights[jn] = 0.0;
            colscale[jn] = 0.0;
            candidates.pop_back();
            if (++skipped > maxskip && maxskip >= 0)
                break;
            continue;
        }

        // Recompute column weight from FTRAN.
        double weight_recomp = 0.0;
        auto sum_used =  [&](Int p, double x) {
            if (tblrow_used[p])
                weight_recomp += x;
        };
        for_each_nonzero(lhs, sum_used);
        assert(std::isfinite(weight_recomp));

        // Update basis.
        const Int jb = basis[pmax];
        basis.TableauRow(jb, lhs, row, true);
        const double pivot = row[jn];
        if (std::abs(pivot) < 1e-3) {
            control_.Debug(3)
                << " |pivot| " << sci2(std::abs(pivot))
                << "(maxvolume)\n";
        }
        bool exchanged;
        errflag = basis.ExchangeIfStable(jb, jn, pivot, 0, &exchanged);
        if (errflag)
            break;
        if (!exchanged)         // factorization was unstable, try again
            continue;
        updates_++;
        volinc_ += std::log2(vmax);

        // Update colscale and invscale_basic.
        const double dn = colscale[jn];
        const double dbinv = invscale_basic[pmax];
        assert(colscale[jb] == 0.0);
        colscale[jb] = 1.0 / invscale_basic[pmax];
        invscale_basic[pmax] = 1.0 / colscale[jn];
        colscale[jn] = 0.0;
        assert(std::isfinite(colscale[jb]));
        assert(std::isfinite(invscale_basic[pmax]));

        // Update column weights.
        const double alpha = (tblrow_used[pmax] - weight_recomp) / (dn*pivot);
        assert(std::isfinite(alpha));
        auto add = [&](Int j, double x) {
            colweights[j] += alpha * x * colscale[j];
        };
        for_each_nonzero(row, add);
        colweights[jb] = tblrow_used[pmax] + alpha/dbinv;
        colweights[jn] = 0.0;
        candidates.clear();
    }

    skipped_ += skipped;
    return errflag;
}

Int Maxvolume::ScaleFtran(double colscale_jn, const Vector& invscale_basic,
                          IndexedVector& ftran) {
    double vmax = 0.0;
    Int pmax = 0;
    auto fcn = [&](Int p, double& pivot) {
        double scaled_pivot = pivot * colscale_jn * invscale_basic[p];
        double v = std::abs(scaled_pivot);
        if (v > vmax && std::abs(pivot) > kPivotZeroTol) {
            vmax = v;
            pmax = p;
        }
        pivot = scaled_pivot;
    };
    for_each_nonzero(ftran, fcn);
    return pmax;
}

}  // namespace ipx
