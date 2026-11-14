#include <cassert>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "ipm/ipx/basiclu_wrapper.h"
#include "ipm/basiclu/basiclu.h"

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
    fill_factor_ = 0.0;
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
    bool has_logged = false;
    double last_log = basiclu_wallclock();
    for (Int ncall = 0; ; ncall++) {
        status = basiclu_factorize(istore_.data(), xstore_.data(),
                                   Li_.data(), Lx_.data(),
                                   Ui_.data(), Ux_.data(),
                                   Wi_.data(), Wx_.data(),
                                   Bbegin, Bend, Bi, Bx, ncall);
	// Possibly log factorization progress 
	logFactorization(status, has_logged, last_log);
        if (status != BASICLU_REALLOCATE)
            break;
        Reallocate();
    }

    if (status == BASICLU_WARNING_timeout) return IPX_ERROR_time_interrupt;

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
    for (;;) {
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
    for (;;) {
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
    for (;;) {
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
    for (;;) {
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
    for (;;) {
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

void BasicLu::_timeStart(double new_time_start) {
    xstore_[BASICLU_TIME_START] = new_time_start;
}

void BasicLu::_timeLimit(double new_time_limit) {
    xstore_[BASICLU_TIME_LIMIT] = new_time_limit;
}

void BasicLu::Reallocate() {
    assert(Li_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYL]));
    assert(Lx_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYL]));
    assert(Ui_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYU]));
    assert(Ux_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYU]));
    assert(Wi_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYW]));
    assert(Wx_.size() == static_cast<size_t>(xstore_[BASICLU_MEMORYW]));

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

void BasicLu::logFactorization(const Int status, bool& has_logged, double& last_log) {

  auto forceLogging = [&] () {
    if ((has_logged && status == BASICLU_OK) ||
	status == BASICLU_WARNING_singular_matrix ||
	status == BASICLU_WARNING_timeout ||
	status == BASICLU_ERROR_invalid_store ||
	status == BASICLU_ERROR_invalid_call ||
	status == BASICLU_ERROR_argument_missing ||
	status == BASICLU_ERROR_invalid_argument ||
	status == BASICLU_ERROR_maximum_updates ||
	status == BASICLU_ERROR_singular_update ||
	status == BASICLU_ERROR_invalid_object ||
	status == BASICLU_ERROR_out_of_memory) return true;
    return false;
  };

  auto statusToString = [&] () {
    if (status == BASICLU_OK) return "OK";
    if (status == BASICLU_REALLOCATE) return "Reallocate";
    if (status == BASICLU_WARNING_singular_matrix) return "Warning: singular matrix";
    if (status == BASICLU_WARNING_timeout) return "Warning: timeout";
    if (status == BASICLU_ERROR_invalid_store) return "Error: invalid store";
    if (status == BASICLU_ERROR_invalid_call) return "Error: invalid call";
    if (status == BASICLU_ERROR_argument_missing) return "Error: argument missing";
    if (status == BASICLU_ERROR_invalid_argument) return "Error: invalid argument";
    if (status == BASICLU_ERROR_maximum_updates) return "Error: maximum updates";
    if (status == BASICLU_ERROR_singular_update) return "Error: singular update";
    if (status == BASICLU_ERROR_invalid_object) return "Error: invalid object";
    if (status == BASICLU_ERROR_out_of_memory) return "Error: out of memory";
    return "Unknown status";
  };

  const double min_log_interval = 5.0;
  double wallclock = basiclu_wallclock();
  if (wallclock > last_log + min_log_interval || forceLogging()) {
    Int num_pivot = xstore_[BASICLU_RANK] + xstore_[BASICLU_RANKDEF];
    Int dim = xstore_[BASICLU_DIM];
    std::stringstream factorization_logging_stream;
    factorization_logging_stream.str(std::string());
    factorization_logging_stream << "   Factorization: pivoted on";
    factorization_logging_stream 
      << " "  << Format(num_pivot, 7) << " of " << Format(dim, 7) << " rows";
    
    if (!control_.timelessLog())
      factorization_logging_stream << "  " << time(control_.Elapsed());
    factorization_logging_stream << " " << statusToString() << "\n";
    control_.hLog(factorization_logging_stream);
     
    has_logged = true;
    last_log = wallclock;
  }
}

double basiclu_wallclock() {
  using namespace std::chrono;
  return duration_cast<duration<double> >(high_resolution_clock::now().time_since_epoch()).count();
}
  

}  // namespace ipx
