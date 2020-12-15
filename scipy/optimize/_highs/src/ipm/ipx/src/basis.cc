// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "basis.h"
#include <algorithm>
#include <cmath>
#include <tuple>
#include "basiclu_kernel.h"
#include "basiclu_wrapper.h"
#include "forrest_tomlin.h"
#include "guess_basis.h"
#include "power_method.h"
#include "symbolic_invert.h"
#include "timer.h"
#include "utils.h"

namespace ipx {

Basis::Basis(const Control& control, const Model& model)
    : control_(control), model_(model) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    basis_.resize(m);
    map2basis_.resize(n+m);
    if (control_.lu_kernel() <= 0) {
        lu_.reset(new BasicLu(control_, m));
    } else {
        std::unique_ptr<LuFactorization> lu(new BasicLuKernel);
        lu_.reset(new ForrestTomlin(control_, m, lu));
    }
    lu_->pivottol(control_.lu_pivottol());
    SetToSlackBasis();
}

void Basis::FixNonbasicVariable(Int j) {
    if (StatusOf(j) == NONBASIC_FIXED)
        return;
    assert(StatusOf(j) == NONBASIC);
    assert(map2basis_[j] == -1);
    map2basis_[j] = -2;
}

void Basis::FreeBasicVariable(Int j) {
    const Int m = model_.rows();
    if (StatusOf(j) == BASIC_FREE)
        return;
    assert(StatusOf(j) == BASIC);
    assert(map2basis_[j] >= 0 && map2basis_[j] < m);
    map2basis_[j] += m;
}

void Basis::UnfixVariables() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    for (Int j = 0; j < n+m; j++)
        if (map2basis_[j] == -2)
            map2basis_[j] = -1;
}

void Basis::UnfreeVariables() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    for (Int j = 0; j < n+m; j++)
        if (map2basis_[j] >= m)
            map2basis_[j] -= m;
}

void Basis::SetToSlackBasis() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    for (Int i = 0; i < m; i++)
        basis_[i] = n+i;
    for (Int j = 0; j < n; j++)
        map2basis_[j] = -1;
    for (Int i = 0; i < m; i++)
        map2basis_[n+i] = i;
    Int err = Factorize();
    // factorization of slack basis cannot fail other than out of memory
    assert(err == 0);
    (void)(err);
}

Int Basis::Load(const int* basic_status) {
    const Int m = model_.rows();
    const Int n = model_.cols();

    // Change member variables only when basis is valid.
    std::vector<Int> basis, map2basis(n+m);
    Int p = 0;
    for (Int j = 0; j < n+m; j++) {
        switch (basic_status[j]) {
        case NONBASIC_FIXED:
            map2basis[j] = -2;
            break;
        case NONBASIC:
            map2basis[j] = -1;
            break;
        case BASIC:
            basis.push_back(j);
            map2basis[j] = p++;
            break;
        case BASIC_FREE:
            basis.push_back(j);
            map2basis[j] = p++ + m;
            break;
        default:
            return IPX_ERROR_invalid_basis;
        }
    }
    if (p != m)
        return IPX_ERROR_invalid_basis;

    std::copy(basis.begin(), basis.end(), basis_.begin());
    std::copy(map2basis.begin(), map2basis.end(), map2basis_.begin());
    return Factorize();
}

Int Basis::Factorize() {
    const Int m = model_.rows();
    const SparseMatrix& AI = model_.AI();
    Timer timer;

    // Build column pointers for passing to LU factorization.
    std::vector<Int> begin(m), end(m);
    for (Int i = 0; i < m; i++) {
        assert(basis_[i] >= 0);
        begin[i] = AI.begin(basis_[i]);
        end[i] = AI.end(basis_[i]);
    }

    Int err = 0;                // return code
    while (true) {
        Int flag = lu_->Factorize(begin.data(), end.data(), AI.rowidx(),
                                  AI.values(), false);
        num_factorizations_++;
        fill_factors_.push_back(lu_->fill_factor());
        if (flag & 2) {
            AdaptToSingularFactorization();
            err = IPX_ERROR_basis_singular;
            break;
        }
        if ((flag & 1) && TightenLuPivotTol()) {
            // The factorization was numerically unstable and the pivot
            // tolerance could be tightened. Repeat.
            continue;
        }
        if (flag & 1) {
            control_.Debug(3)
                << " LU factorization unstable with pivot tolerance "
                << lu_->pivottol() << '\n';
            // ignore instability
        }
        break;
    }
    time_factorize_ += timer.Elapsed();
    factorization_is_fresh_ = true;
    return err;
}

bool Basis::FactorizationIsFresh() const {
    return factorization_is_fresh_;
}

void Basis::GetLuFactors(SparseMatrix *L, SparseMatrix *U, Int *rowperm,
                         Int *colperm) const {
    assert(FactorizationIsFresh());
    lu_->GetFactors(L, U, rowperm, colperm, nullptr);
}

void Basis::SolveDense(const Vector& rhs, Vector& lhs, char trans) const {
    lu_->SolveDense(rhs, lhs, trans);
}

void Basis::SolveForUpdate(Int j, IndexedVector& lhs) {
    const Int p = PositionOf(j);
    Timer timer;
    if (p < 0) {                // ftran
        const SparseMatrix& AI = model_.AI();
        Int begin = AI.begin(j);
        Int end = AI.end(j);
        Int nz = end-begin;
        const Int* bi = AI.rowidx() + begin;
        const double* bx = AI.values() + begin;
        lu_->FtranForUpdate(nz, bi, bx, lhs);
        num_ftran_++;
        if (lhs.sparse())
            num_ftran_sparse_++;
        time_ftran_ += timer.Elapsed();
    }
    else {                      // btran
        lu_->BtranForUpdate(p, lhs);
        num_btran_++;
        if (lhs.sparse())
            num_btran_sparse_++;
        time_btran_ += timer.Elapsed();
    }
}

void Basis::SolveForUpdate(Int j) {
    const Int p = PositionOf(j);
    Timer timer;
    if (p < 0) {
        const SparseMatrix& AI = model_.AI();
        Int begin = AI.begin(j);
        Int end = AI.end(j);
        Int nz = end-begin;
        const Int* bi = AI.rowidx() + begin;
        const double* bx = AI.values() + begin;
        lu_->FtranForUpdate(nz, bi, bx);
        time_ftran_ += timer.Elapsed();
    }
    else {
        lu_->BtranForUpdate(p);
        time_btran_ += timer.Elapsed();
    }
}

void Basis::TableauRow(Int jb, IndexedVector& btran, IndexedVector& row,
                       bool ignore_fixed) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    assert(IsBasic(jb));
    SolveForUpdate(jb, btran);

    // Estimate if tableau row is sparse.
    bool is_sparse = btran.sparse();
    if (is_sparse) {
        const SparseMatrix& AIt = model_.AIt();
        const Int* bi = btran.pattern();
        Int nz = 0;
        for (Int k = 0; k < btran.nnz(); k++) {
            Int i = bi[k];
            nz += AIt.end(i) - AIt.begin(i);
        }
        nz /= 2;                // guess for overlap
        if (nz > kHypersparseThreshold*n)
            is_sparse = false;
    }

    if (is_sparse) {
        // sparse vector * sparse matrix: accesses A rowwise
        const SparseMatrix& AIt = model_.AIt();
        const Int* Ati = AIt.rowidx();
        const double* Atx = AIt.values();
        const Int* bi = btran.pattern();
        row.set_to_zero();
        Int* row_pattern = row.pattern();
        Int nz = 0;
        for (Int k = 0; k < btran.nnz(); k++) {
            Int i = bi[k];
            double temp = btran[i];
            Int begin = AIt.begin(i);
            Int end = AIt.end(i);
            for (Int p = begin; p < end; p++) {
                Int j = Ati[p];
                if (map2basis_[j] == -1 || (map2basis_[j] == -2 && !ignore_fixed))
                {
                    map2basis_[j] -= 2; // mark column
                    row_pattern[nz++] = j;
                }
                if (map2basis_[j] < -2) // marked column
                    row[j] += temp * Atx[p];
            }
        }
        for (Int k = 0; k < nz; k++) // reset marked
            map2basis_[row_pattern[k]] += 2;
        row.set_nnz(nz);
    } else {
        // dense vector * sparse matrix: accesses A columnwise
        const SparseMatrix& AI = model_.AI();
        const Int* Ai = AI.rowidx();
        const double* Ax = AI.values();
        for (Int j = 0; j < n+m; j++) {
            double result = 0.0;
            if (map2basis_[j] == -1 || (map2basis_[j] == -2 && !ignore_fixed)) {
                Int begin = AI.begin(j);
                Int end = AI.end(j);
                for (Int p = begin; p < end; p++)
                    result += Ax[p] * btran[Ai[p]];
            }
            row[j] = result;
        }
        row.InvalidatePattern();
    }
}

Int Basis::ExchangeIfStable(Int jb, Int jn, double tableau_entry, int sys,
                       bool* exchanged) {
    assert(IsBasic(jb));
    assert(IsNonbasic(jn));
    if (sys > 0)                // forward system needs to be solved
        SolveForUpdate(jn);
    if (sys < 0)                // transposed system needs to be solved
        SolveForUpdate(jb);

    // Update factorization.
    *exchanged = false;
    Timer timer;
    Int err = lu_->Update(tableau_entry);
    time_update_ += timer.Elapsed();
    if (err != 0) {
        if (FactorizationIsFresh() && !TightenLuPivotTol())
            return IPX_ERROR_basis_too_ill_conditioned;
        control_.Debug(3)
            << " stability check forced refactorization after "
            << lu_->updates()-1 << " updates\n";
        return Factorize();     // refactorizes the old basis
    }

    // Update basis.
    Int ib = PositionOf(jb);
    assert(basis_[ib] == jb);
    basis_[ib] = jn;
    map2basis_[jn] = ib;        // status now BASIC
    map2basis_[jb] = -1;        // status now NONBASIC
    num_updates_++;
    factorization_is_fresh_ = false;
    *exchanged = true;

    if (lu_->NeedFreshFactorization())
        return Factorize();
    return 0;
}

void Basis::ComputeBasicSolution(Vector& x, Vector& y, Vector& z) const {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& b = model_.b();
    const Vector& c = model_.c();
    const SparseMatrix& AI = model_.AI();

    // Compute x[basic] so that Ax=b. Use y as workspace.
    y = b;
    for (Int j = 0; j < n+m; j++)
        if (IsNonbasic(j))
            ScatterColumn(AI, j, -x[j], y);
    SolveDense(y, y, 'N');
    for (Int p = 0; p < m; p++)
        x[basis_[p]] = y[p];

    // Compute y and z[nonbasic] so that AI'y+z=c.
    for (Int p = 0; p < m; p++)
        y[p] = c[basis_[p]] - z[basis_[p]];
    SolveDense(y, y, 'T');
    for (Int j = 0; j < n+m; j++) {
        if (IsNonbasic(j))
            z[j] = c[j] - DotColumn(AI, j, y);
    }
}

void Basis::ConstructBasisFromWeights(const double* colscale, Info* info) {
    // const Int m = model_.rows();
    // const Int n = model_.cols();
    assert(colscale);
    info->errflag = 0;
    info->dependent_rows = 0;
    info->dependent_cols = 0;

    if (control_.crash_basis()) {
        CrashBasis(colscale);
        double sigma = MinSingularValue();
        control_.Debug()
            << Textline("Minimum singular value of crash basis:") << sci2(sigma)
            << '\n';
        Repair(info);
        if (info->basis_repairs < 0) {
            control_.Log() << " discarding crash basis\n";
            SetToSlackBasis();
        }
        else if (info->basis_repairs > 0) {
            sigma = MinSingularValue();
            control_.Debug()
                << Textline("Minimum singular value of repaired crash basis:")
                << sci2(sigma) << '\n';
        }
    } else {
        SetToSlackBasis();
    }
    PivotFreeVariablesIntoBasis(colscale, info);
    if (info->errflag)
        return;
    PivotFixedVariablesOutOfBasis(colscale, info);
    if (info->errflag)
        return;
}

double Basis::MinSingularValue() const {
    const Int m = model_.rows();
    Vector v(m);

    // Computes maximum eigenvalue of inverse(B*B').
    double lambda = PowerMethod(
        [this](const Vector& x, Vector& fx) {
            SolveDense(x, fx, 'N');
            SolveDense(fx, fx, 'T'); }, v);
    return std::sqrt(1.0/lambda);
}

void Basis::SymbolicInvert(Int* rowcounts, Int* colcounts) const {
    ipx::SymbolicInvert(model_, basis_, rowcounts, colcounts);
}

double Basis::DensityInverse() const {
    Int m = model_.rows();
    std::vector<Int> rowcounts(m);
    SymbolicInvert(rowcounts.data(), nullptr);
    // Accumulating rowcounts would result in overflow for large LPs.
    double density = 0.0;
    for (Int i = 0; i < m; i++)
        density += 1.0*rowcounts[i]/m;
    return density/m;
}

const Model& Basis::model() const {
    return model_;
}

Int Basis::factorizations() const {
    return num_factorizations_;
}

Int Basis::updates_total() const {
    return num_updates_;
}

double Basis::frac_ftran_sparse() const {
    return 1.0 * num_ftran_sparse_ / num_ftran_;
}

double Basis::frac_btran_sparse() const {
    return 1.0 * num_btran_sparse_ / num_btran_;
}

double Basis::time_factorize() const {
    return time_factorize_;
}

double Basis::time_ftran() const {
    return time_ftran_;
}

double Basis::time_btran() const {
    return time_btran_;
}

double Basis::time_update() const {
    return time_update_;
}

double Basis::mean_fill() const {
    if (fill_factors_.empty())
        return 0.0;
    double mean = 1.0;
    Int num_factors = fill_factors_.size();
    for (double f : fill_factors_)
        mean *= std::pow(f, 1.0/num_factors);
    return mean;
}

double Basis::max_fill() const {
    if (fill_factors_.empty())
        return 0.0;
    return *std::max_element(fill_factors_.begin(), fill_factors_.end());
}

Int Basis::AdaptToSingularFactorization() {
    const Int m = model_.rows();
    const Int n = model_.cols();
    std::vector<Int> rowperm(m), colperm(m), dependent_cols;

    lu_->GetFactors(nullptr, nullptr, rowperm.data(), colperm.data(),
                    &dependent_cols);
    for (Int k : dependent_cols) {
        // Column p of the basis matrix was replaced by the i-th unit
        // column. Insert the corresponding slack variable jn into
        // position p of the basis.
        Int p = colperm[k];
        Int i = rowperm[k];
        Int jb = basis_[p];
        Int jn = n+i;
        assert(map2basis_[jn] < 0);
        basis_[p] = jn;
        map2basis_[jn] = p; // now BASIC at position p
        if (jb >= 0)
            map2basis_[jb] = -1; // now NONBASIC
    }
    return dependent_cols.size();
}

bool Basis::TightenLuPivotTol() {
    double tol = lu_->pivottol();
    if (tol <= 0.05)
        lu_->pivottol(0.1);
    else if (tol <= 0.25)
        lu_->pivottol(0.3);
    else if (tol <= 0.5)
        lu_->pivottol(0.9);
    else
        return false;
    control_.Log()
        << " LU pivot tolerance tightened to " << lu_->pivottol() << '\n';
    return true;
}

void Basis::CrashBasis(const double* colweights) {
    const Int m = model_.rows();

    // Make a guess for a basis. Then use LU factorization with a strict
    // absolute pivot tolerance to remove dependent columns. This is not a
    // rank revealing factorization, but it detects many dependencies in
    // practice.
    std::vector<Int> cols_guessed = GuessBasis(control_, model_, colweights);
    assert((int)cols_guessed.size() <= m);
    assert((int)cols_guessed.size() == m); // at the moment

    // Initialize the Basis object and factorize the (partial) basis. If
    // basis_[p] is negative, the p-th column of the basis matrix is zero,
    // and a slack column will be inserted by CrashFacorize().
    std::fill(basis_.begin(), basis_.end(), -1);
    std::fill(map2basis_.begin(), map2basis_.end(), -1);
    for (Int k = 0; k < (Int) cols_guessed.size(); k++) {
        basis_[k] = cols_guessed[k];
        assert(map2basis_[basis_[k]] == -1); // must not have duplicates
        map2basis_[basis_[k]] = k;
    }
    Int num_dropped = 0;
    CrashFactorize(&num_dropped);
    control_.Debug()
        << Textline("Number of columns dropped from guessed basis:")
        << num_dropped << '\n';
    (void)(m);
}

// Rook search for a large entry in inverse(B). If successful, returns [p,i,x]
// where x is the entry at position (p,i) of inverse(B). Notice that p
// corresponds to a column of B and i corresponds to a row of B. If overflow
// occurs, returns [-1,-1,INFINITY]. See [1] for a discussion of the algorithm.
//
// [1] N.J. Higham and S.D. Relton, "Estimating the Largest Elements of a
//     Matrix", MIMS EPrint 2015.116, University of Manchester (2015).
//
static std::tuple<Int,Int,double> InverseSearch(const Basis& basis,
                                                Vector& work) {
    const Int m = work.size();
    double inverse_max = 0.0;

    for (Int i = 0; i < m; i++)
        work[i] = 1.0/(i+1);

    while (true) {
        basis.SolveDense(work, work, 'N');
        if (!AllFinite(work))
            break;
        Int pmax = FindMaxAbs(work);
        work = 0.0;
        work[pmax] = 1.0;
        basis.SolveDense(work, work, 'T');
        if (!AllFinite(work))
            break;
        Int imax = FindMaxAbs(work);
        double inverse_entry = work[imax];
        if (std::abs(inverse_entry) <= 2.0*inverse_max)
            return std::make_tuple(pmax, imax, inverse_entry);
        inverse_max = std::abs(inverse_entry);
        work = 0.0;
        work[imax] = 1.0;
    }
    return std::make_tuple(-1,-1,INFINITY); // failure
}

void Basis::Repair(Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    Vector work(m);
    info->basis_repairs = 0;

    while (true) {
        std::tuple<Int,Int,double> entry = InverseSearch(*this, work);
        Int pmax = std::get<0>(entry);
        Int imax = std::get<1>(entry);
        double pivot = std::get<2>(entry);
        if (pmax < 0 || imax < 0 || !std::isfinite(pivot)) {
            info->basis_repairs = -1;
            break;
        }
        if (std::abs(pivot) < kBasisRepairThreshold)
            break;
        Int jb = basis_[pmax];
        Int jn = n + imax;
        if (!IsNonbasic(jn)) {
            info->basis_repairs = -2;
            break;
        }
        if (info->basis_repairs >= kMaxBasisRepair) {
            info->basis_repairs = -3;
            break;
        }
        SolveForUpdate(jb);
        SolveForUpdate(jn);
        CrashExchange(jb, jn, pivot, 0, nullptr);
        info->basis_repairs++;
        control_.Debug(3) << " basis repair: |pivot| = "
                      << sci2(std::abs(pivot)) << '\n';
    }
}

void Basis::CrashFactorize(Int* num_dropped) {
    const Int m = model_.rows();
    const SparseMatrix& AI = model_.AI();
    Timer timer;

    // Build column pointers for passing to LU factorization. A negative index
    // in basis_ means an empty slot. This option is kept for use by the crash
    // procedure, if an incomplete basis was constructed. The zero column in the
    // basis matrix causes a singularity in the LU factorization, so that the
    // empty slot will be replaced by a slack variable below.
    std::vector<Int> begin(m), end(m);
    for (Int i = 0; i < m; i++) {
        if (basis_[i] >= 0) {
            begin[i] = AI.begin(basis_[i]);
            end[i] = AI.end(basis_[i]);
        } else {
            begin[i] = 0;
            end[i] = 0;
        }
    }
    Int flag = lu_->Factorize(begin.data(), end.data(), AI.rowidx(),
                              AI.values(), true);
    num_factorizations_++;
    fill_factors_.push_back(lu_->fill_factor());
    Int ndropped = 0;
    if (flag & 2)
        ndropped = AdaptToSingularFactorization();
    if (num_dropped)
        *num_dropped = ndropped;

    time_factorize_ += timer.Elapsed();
    factorization_is_fresh_ = true;

    #ifndef NDEBUG
    // All empty slots must have been replaced by slack variables.
    for (Int i = 0; i < m; i++)
        assert(basis_[i] >= 0);
    #endif
}

void Basis::CrashExchange(Int jb, Int jn, double tableau_entry, int sys,
                          Int* num_dropped) {
    assert(IsBasic(jb));
    assert(IsNonbasic(jn));
    if (sys > 0)                // forward system needs to be solved
        SolveForUpdate(jn);
    if (sys < 0)                // transposed system needs to be solved
        SolveForUpdate(jb);

    // Update basis.
    Int ib = PositionOf(jb);
    assert(basis_[ib] == jb);
    basis_[ib] = jn;
    map2basis_[jn] = ib;        // status now BASIC
    map2basis_[jb] = -1;        // status now NONBASIC
    num_updates_++;
    factorization_is_fresh_ = false;

    // Update factorization.
    if (num_dropped)
        *num_dropped = 0;
    Timer timer;
    Int err = lu_->Update(tableau_entry);
    time_update_ += timer.Elapsed();
    if (err != 0 || lu_->NeedFreshFactorization()) {
        control_.Debug(3) << " refactorization required in CrashExchange()\n";
        CrashFactorize(num_dropped);
    }
}

void Basis::PivotFreeVariablesIntoBasis(const double* colweights, Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    IndexedVector ftran(m);
    const double dependency_tol = std::max(0.0, control_.dependency_tol());
    info->errflag = 0;
    info->dependent_cols = 0;
    Int stability_pivots = 0;

    // Maintain stack of free nonbasic variables.
    std::vector<Int> remaining;
    for (Int j = 0; j < n+m; j++) {
        if (std::isinf(colweights[j]) && map2basis_[j] < 0)
            remaining.push_back(j);
    }
    control_.Debug() << Textline("Number of free variables nonbasic:")
                  << remaining.size() << '\n';

    control_.ResetPrintInterval();
    while (!remaining.empty()) {
        Int jn = remaining.back();
        assert(std::isinf(colweights[jn]));
        assert(map2basis_[jn] < 0);
        if ((info->errflag = control_.InterruptCheck()) != 0)
            return;

        SolveForUpdate(jn, ftran);
        Int pmax = -1;
        Int pmax_nonfree = -1;
        double fmax = 0.0;
        double fmax_nonfree = 0.0;
        auto update_max = [&](Int p, double f) {
            f = std::abs(f);
            if (f > fmax) {
                fmax = f;
                pmax = p;
            }
            if (!std::isinf(colweights[basis_[p]])) {
                if (f > fmax_nonfree) {
                    fmax_nonfree = f;
                    pmax_nonfree = p;
                }
            }
        };
        for_each_nonzero(ftran, update_max);

        if (fmax > 4.0 && fmax_nonfree < 1.0) {
            Int jb = basis_[pmax];
            assert(std::isinf(colweights[jb]));
            bool exchanged;
            info->errflag = ExchangeIfStable(jb, jn, ftran[pmax], -1,
                                             &exchanged);
            if (info->errflag)
                return;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            remaining.pop_back();
            remaining.push_back(jb);
            info->updates_start++;
            stability_pivots++;
        } else {
            if (fmax_nonfree <= dependency_tol) {
                // jn cannot be pivoted into the basis. If we do not have an
                // unbounded primal ray yet, then test if column jn yields one.
                // Compute the change in the primal objective that is caused by
                // a unit increase of x[jn] with corresponding adjustment of
                // free basic variables.
                if (!info->cols_inconsistent) {
                    const Vector& c = model_.c();
                    double delta_obj = c[jn];
                    auto update_delta_obj = [&](Int p, double f) {
                        Int j = basis_[p];
                        if (std::isinf(colweights[j]))
                            delta_obj -= c[j] * f;
                    };
                    for_each_nonzero(ftran, update_delta_obj);
                    if (std::abs(delta_obj) > dependency_tol) {
                        control_.Debug()
                            << Textline(
                                "Unbounded primal ray with objective change:")
                            << sci2(delta_obj) << '\n';
                        info->cols_inconsistent = true;
                    }
                }
                info->dependent_cols++;
                remaining.pop_back();
            } else {
                Int jb = basis_[pmax_nonfree];
                bool exchanged;
                info->errflag = ExchangeIfStable(jb, jn, ftran[pmax_nonfree],
                                                 -1, &exchanged);
                if (info->errflag)
                    return;
                if (!exchanged)     // factorization was unstable, try again
                    continue;
                remaining.pop_back();
                info->updates_start++;
            }
        }
        control_.IntervalLog()
            << " " << remaining.size() << " free variables remaining\n";
    }
    control_.Debug()
        << Textline("Number of free variables swapped for stability:")
        << stability_pivots << '\n';
}

void Basis::PivotFixedVariablesOutOfBasis(const double* colweights, Info* info){
    const Int m = model_.rows();
    const Int n = model_.cols();
    IndexedVector btran(m), row(n+m);
    const double dependency_tol = std::max(0.0, control_.dependency_tol());
    info->errflag = 0;
    info->dependent_rows = 0;
    Int stability_pivots = 0;

    // Maintain stack of fixed basic variables.
    std::vector<Int> remaining;
    for (Int j = 0; j < n; j++) {
        // Structural variables with zero weight must not be in the crash basis.
        if (colweights[j] == 0.0)
            assert(map2basis_[j] < 0);
    }
    for (Int j = n; j < n+m; j++) {
        if (colweights[j] == 0.0 && map2basis_[j] >= 0)
            remaining.push_back(j);
    }
    control_.Debug() << Textline("Number of fixed variables basic:")
                  << remaining.size() << '\n';

    control_.ResetPrintInterval();
    while (!remaining.empty()) {
        Int jb = remaining.back();
        assert(colweights[jb] == 0.0);
        assert(map2basis_[jb] >= 0);
        if ((info->errflag = control_.InterruptCheck()) != 0)
            return;

        TableauRow(jb, btran, row);
        Int jmax = -1;
        Int jmax_nonfixed = -1;
        double rmax = 0.0;
        double rmax_nonfixed = 0.0;
        auto update_max = [&](Int j, double r) {
            // Ignore structural variables with zero weight.
            if (j >= n || colweights[j] != 0.0) {
                r = std::abs(r);
                if (r > rmax) {
                    rmax = r;
                    jmax = j;
                }
                if (colweights[j] != 0.0) {
                    if (r > rmax_nonfixed) {
                        rmax_nonfixed = r;
                        jmax_nonfixed = j;
                    }
                }
            }
        };
        for_each_nonzero(row, update_max);

        if (rmax > 4.0 && rmax_nonfixed < 1.0) {
            assert(colweights[jmax] == 0.0);
            bool exchanged;
            info->errflag = ExchangeIfStable(jb, jmax, row[jmax], 1,
                                             &exchanged);
            if (info->errflag)
                return;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            remaining.pop_back();
            remaining.push_back(jmax);
            info->updates_start++;
            stability_pivots++;
        } else {
            if (rmax_nonfixed <= dependency_tol) {
                // jb cannot be pivoted out of the basis. If we do not have an
                // unbounded dual ray yet, then test if row jb-n yields one.
                // Compute the change in the dual objective that is caused by a
                // unit increase of y[jb-n] with corresponding adjustment of the
                // remaining y[i].
                if (!info->rows_inconsistent) {
                    double delta_obj = Dot(btran, model_.b());
                    if (std::abs(delta_obj) > dependency_tol) {
                        control_.Debug()
                            << Textline(
                                "Unbounded dual ray with objective change:")
                            << sci2(delta_obj) << '\n';
                        info->rows_inconsistent = true;
                    }
                }
                info->dependent_rows++;
                remaining.pop_back();
            } else {
                // jb can be exchanged by a non-fixed variable.
                // Among all numerically stable pivots, choose the one that
                // maximizes volume of the basis matrix.
                Int jmax_scaled = -1;
                double rmax_scaled = 0.0;
                auto update_max = [&](Int j, double r) {
                    r = std::abs(r);
                    double rscaled = r * colweights[j];
                    if (r >= 0.1*rmax_nonfixed && rscaled > rmax_scaled) {
                        rmax_scaled = rscaled;
                        jmax_scaled = j;
                    }
                };
                for_each_nonzero(row, update_max);
                assert(jmax_scaled >= 0);
                bool exchanged;
                double pivot = row[jmax_scaled];
                info->errflag = ExchangeIfStable(jb, jmax_scaled, pivot, 1,
                                                 &exchanged);
                if (info->errflag)
                    return;
                if (!exchanged)     // factorization was unstable, try again
                    continue;
                remaining.pop_back();
                info->updates_start++;
            }
        }
        control_.IntervalLog()
            << " " << remaining.size() << " fixed variables remaining\n";
    }
    control_.Debug()
        << Textline("Number of fixed variables swapped for stability:")
        << stability_pivots << '\n';
}

Vector CopyBasic(const Vector& x, const Basis& basis) {
    const Int m = basis.model().rows();
    Vector xbasic(m);
    for (Int p = 0; p < m; p++)
        xbasic[p] = x[basis[p]];
    return xbasic;
}

}  // namespace ipx
