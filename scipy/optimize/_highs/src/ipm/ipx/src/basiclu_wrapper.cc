// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "basiclu_wrapper.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "basiclu.h"

namespace ipx {

BasicLu::BasicLu(const Control& control, Int dim) : control_(control) {
    static_assert(sizeof(Int) == sizeof(lu_int),
                  "IPX integer type does not match BASICLU integer type");
    istore_.resize(BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * dim);
    xstore_.resize(BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * dim);

    Int status = basiclu_initialize(dim, istore_.data(), xstore_.data());
    if (status != BASICLU_OK)
        throw std::logic_error("basiclu_initialize failed");

    // Set initial size of BASICLU work arrays to 1 element, so that data()
    // does not return NULL (not sure if it would if size is 0).
    Li_.resize(1);
    Lx_.resize(1);
    Ui_.resize(1);
    Ux_.resize(1);
    Wi_.resize(1);
    Wx_.resize(1);
    xstore_[BASICLU_MEMORYL] = 1;
    xstore_[BASICLU_MEMORYU] = 1;
    xstore_[BASICLU_MEMORYW] = 1;
}

Int BasicLu::_Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                        const double* Bx, bool strict_abs_pivottol) {
    Int status;
    if (strict_abs_pivottol) {
        xstore_[BASICLU_REMOVE_COLUMNS] = 1;
        xstore_[BASICLU_ABS_PIVOT_TOLERANCE] = kLuDependencyTol;
    } else {
        xstore_[BASICLU_REMOVE_COLUMNS] = 0;
        xstore_[BASICLU_ABS_PIVOT_TOLERANCE] = 1e-14; // BASICLU default
    }
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_factorize(istore_.data(), xstore_.data(),
                                   Li_.data(), Lx_.data(),
                                   Ui_.data(), Ux_.data(),
                                   Wi_.data(), Wx_.data(),
                                   Bbegin, Bend, Bi, Bx, ncall);
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK && status != BASICLU_WARNING_singular_matrix)
        throw std::logic_error("basiclu_factorize failed");

    Int matrix_nz = xstore_[BASICLU_MATRIX_NZ];
    Int lnz = xstore_[BASICLU_LNZ];
    Int unz = xstore_[BASICLU_UNZ];
    Int dim = xstore_[BASICLU_DIM];
    fill_factor_ = 1.0 * (lnz+unz+dim) / matrix_nz;

    double normLinv = xstore_[BASICLU_NORMEST_LINV];
    double normUinv = xstore_[BASICLU_NORMEST_UINV];
    double stability = xstore_[BASICLU_RESIDUAL_TEST];
    control_.Debug(3)
        << " normLinv = " << sci2(normLinv) << ','
        << " normUinv = " << sci2(normUinv) << ','
        << " stability = " << sci2(stability) << '\n';

    Int ret = 0;
    if (stability > kLuStabilityThreshold)
        ret |= 1;
    if (status == BASICLU_WARNING_singular_matrix)
        ret |= 2;
    return ret;
}

void BasicLu::_GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                          Int* colperm, std::vector<Int>* dependent_cols) {
    Int *Lbegin = nullptr, *Ubegin = nullptr;
    Int *Lindex = nullptr, *Uindex = nullptr;
    double *Lvalue = nullptr, *Uvalue = nullptr;
    Int dim = xstore_[BASICLU_DIM];
   
    if (L) {
        Int lnz = xstore_[BASICLU_LNZ];
        L->resize(dim, dim, dim+lnz);
        Lbegin = L->colptr();
        Lindex = L->rowidx();
        Lvalue = L->values();
    }
    if (U) {
        Int unz = xstore_[BASICLU_UNZ];
        U->resize(dim, dim, dim+unz);
        Ubegin = U->colptr();
        Uindex = U->rowidx();
        Uvalue = U->values();
    }
    Int status = basiclu_get_factors(istore_.data(), xstore_.data(),
                                     Li_.data(), Lx_.data(),
                                     Ui_.data(), Ux_.data(),
                                     Wi_.data(), Wx_.data(),
                                     rowperm, colperm,
                                     Lbegin, Lindex, Lvalue,
                                     Ubegin, Uindex, Uvalue);
    if (status != BASICLU_OK)
        throw std::logic_error("basiclu_get_factors failed");

    if (L) {
        // Remove unit diagonal from L.
        Int num_dropped = RemoveDiagonal(*L, nullptr);
        assert(num_dropped == dim);
        (void)(num_dropped);
    }
    if (dependent_cols) {
        // Dependent columns are at the end of the BASICLU pivot sequence.
        Int rank = xstore_[BASICLU_RANK];
        dependent_cols->clear();
        for (Int k = rank; k < dim; k++)
            dependent_cols->push_back(k);
    }
}

void BasicLu::_SolveDense(const Vector& rhs, Vector& lhs, char trans) {
    Int status = basiclu_solve_dense(istore_.data(), xstore_.data(),
                                     Li_.data(), Lx_.data(),
                                     Ui_.data(), Ux_.data(),
                                     Wi_.data(), Wx_.data(),
                                     &rhs[0], &lhs[0], trans);
    if (status != BASICLU_OK)
        throw std::logic_error("basiclu_solve_dense failed");
}

void BasicLu::_FtranForUpdate(Int nzrhs, const Int* bi, const double* bx) {
    Int status;
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_solve_for_update(istore_.data(), xstore_.data(),
                                          Li_.data(), Lx_.data(),
                                          Ui_.data(), Ux_.data(),
                                          Wi_.data(), Wx_.data(),
                                          nzrhs, bi, bx,
                                          nullptr, nullptr, nullptr, 'N');
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK)
        throw std::logic_error(
            "basiclu_solve_for_update (ftran without lhs) failed");
}

void BasicLu::_FtranForUpdate(Int nzrhs, const Int* bi, const double* bx,
                              IndexedVector& lhs) {
    Int status;
    Int nzlhs = 0;
    lhs.set_to_zero();
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_solve_for_update(istore_.data(), xstore_.data(),
                                          Li_.data(), Lx_.data(),
                                          Ui_.data(), Ux_.data(),
                                          Wi_.data(), Wx_.data(),
                                          nzrhs, bi, bx,
                                          &nzlhs, lhs.pattern(), lhs.elements(),
                                          'N');
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK)
        throw std::logic_error(
            "basiclu_solve_for_update (ftran with lhs) failed");
    lhs.set_nnz(nzlhs); 
}

void BasicLu::_BtranForUpdate(Int j) {
    Int status;
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_solve_for_update(istore_.data(), xstore_.data(),
                                          Li_.data(), Lx_.data(),
                                          Ui_.data(), Ux_.data(),
                                          Wi_.data(), Wx_.data(),
                                          0, &j, nullptr,
                                          nullptr, nullptr, nullptr, 'T');
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK)
        throw std::logic_error(
            "basiclu_solve_for_update (btran without lhs) failed");
}

void BasicLu::_BtranForUpdate(Int j, IndexedVector& lhs) {
    Int status;
    Int nzlhs = 0;
    lhs.set_to_zero();
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_solve_for_update(istore_.data(), xstore_.data(),
                                          Li_.data(), Lx_.data(),
                                          Ui_.data(), Ux_.data(),
                                          Wi_.data(), Wx_.data(),
                                          0, &j, nullptr,
                                          &nzlhs, lhs.pattern(), lhs.elements(),
                                          'T');
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK)
        throw std::logic_error(
            "basiclu_solve_for_update (btran with lhs) failed");
    lhs.set_nnz(nzlhs);
}

Int BasicLu::_Update(double pivot) {
    double max_eta_old = xstore_[BASICLU_MAX_ETA];
    Int status;
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_update(istore_.data(), xstore_.data(),
                                Li_.data(), Lx_.data(),
                                Ui_.data(), Ux_.data(),
                                Wi_.data(), Wx_.data(), pivot);
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }
    if (status != BASICLU_OK && status !=  BASICLU_ERROR_singular_update)
        throw std::logic_error("basiclu_update failed");
    if (status == BASICLU_ERROR_singular_update)
        return -1;

    // Print a debugging message if a new eta entry is large.
    double max_eta = xstore_[BASICLU_MAX_ETA];
    if (max_eta > 1e10 && max_eta > max_eta_old)
        control_.Debug(3) << " max eta = " << sci2(max_eta) << '\n';

    // stability check
    double pivot_error = xstore_[BASICLU_PIVOT_ERROR];
    if (pivot_error > kFtDiagErrorTol) {
        control_.Debug(3)
            << " relative error in new diagonal entry of U = "
            << sci2(pivot_error) << '\n';
        return 1;
    }
    return 0;
}

bool BasicLu::_NeedFreshFactorization() {
    Int dim = xstore_[BASICLU_DIM];
    Int nforrest = xstore_[BASICLU_NFORREST];
    double update_cost = xstore_[BASICLU_UPDATE_COST];

    return nforrest == dim || update_cost > 1.0;
}

double BasicLu::_fill_factor() const {
    return fill_factor_;
}

double BasicLu::_pivottol() const {
    return xstore_[BASICLU_REL_PIVOT_TOLERANCE];
}

void BasicLu::_pivottol(double new_pivottol) {
    xstore_[BASICLU_REL_PIVOT_TOLERANCE] = new_pivottol;
}

void BasicLu::Reallocate() {
    assert(Li_.size() == xstore_[BASICLU_MEMORYL]);
    assert(Lx_.size() == xstore_[BASICLU_MEMORYL]);
    assert(Ui_.size() == xstore_[BASICLU_MEMORYU]);
    assert(Ux_.size() == xstore_[BASICLU_MEMORYU]);
    assert(Wi_.size() == xstore_[BASICLU_MEMORYW]);
    assert(Wx_.size() == xstore_[BASICLU_MEMORYW]);

    if (xstore_[BASICLU_ADD_MEMORYL] > 0) {
        Int new_size = xstore_[BASICLU_MEMORYL] + xstore_[BASICLU_ADD_MEMORYL];
        new_size *= kReallocFactor;
        Li_.resize(new_size);
        Lx_.resize(new_size);
        xstore_[BASICLU_MEMORYL] = new_size;
    }
    if (xstore_[BASICLU_ADD_MEMORYU] > 0) {
        Int new_size = xstore_[BASICLU_MEMORYU] + xstore_[BASICLU_ADD_MEMORYU];
        new_size *= kReallocFactor;
        Ui_.resize(new_size);
        Ux_.resize(new_size);
        xstore_[BASICLU_MEMORYU] = new_size;
    }
    if (xstore_[BASICLU_ADD_MEMORYW] > 0) {
        Int new_size = xstore_[BASICLU_MEMORYW] + xstore_[BASICLU_ADD_MEMORYW];
        new_size *= kReallocFactor;
        Wi_.resize(new_size);
        Wx_.resize(new_size);
        xstore_[BASICLU_MEMORYW] = new_size;
    }
}

}  // namespace ipx
