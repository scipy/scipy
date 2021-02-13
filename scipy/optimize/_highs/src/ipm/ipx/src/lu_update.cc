// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "lu_update.h"

namespace ipx {

Int LuUpdate::Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                        const double* Bx, bool strict_abs_pivottol) {
    updates_ = 0;
    return _Factorize(Bbegin, Bend, Bi, Bx, strict_abs_pivottol);
}

void LuUpdate::GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                          Int* colperm, std::vector<Int>* dependent_cols) {
    _GetFactors(L, U, rowperm, colperm, dependent_cols);
}

void LuUpdate::SolveDense(const Vector& rhs, Vector& lhs, char trans) {
    _SolveDense(rhs, lhs, trans);
}

void LuUpdate::FtranForUpdate(Int nz, const Int* bi, const double* bx) {
    _FtranForUpdate(nz, bi, bx);
}

void LuUpdate::FtranForUpdate(Int nz, const Int* bi, const double* bx,
                              IndexedVector& lhs) {
    _FtranForUpdate(nz, bi, bx, lhs);
}

void LuUpdate::BtranForUpdate(Int p) {
    _BtranForUpdate(p);
}

void LuUpdate::BtranForUpdate(Int p, IndexedVector& lhs) {
    _BtranForUpdate(p, lhs);
}

Int LuUpdate::Update(double pivot) {
    updates_++;
    return _Update(pivot);
}

bool LuUpdate::NeedFreshFactorization() {
    return _NeedFreshFactorization();
}

double LuUpdate::fill_factor() const {
    return _fill_factor();
}

double LuUpdate::pivottol() const {
    return _pivottol();
}

void LuUpdate::pivottol(double new_pivottol) {
    _pivottol(new_pivottol);
}

Int LuUpdate::updates() const {
    return updates_;
}

}  // namespace ipx
