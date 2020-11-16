// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_BASICLU_WRAPPER_H_
#define IPX_BASICLU_WRAPPER_H_

#include "control.h"
#include "lu_update.h"

namespace ipx {

class BasicLu : public LuUpdate {
public:
    BasicLu(const Control& control, Int dim);
    ~BasicLu() = default;

private:
    Int _Factorize(const Int* Bbegin, const Int* Bend, const Int* Bi,
                   const double* Bx, bool strict_abs_pivottol) override;
    void _GetFactors(SparseMatrix* L, SparseMatrix* U, Int* rowperm,
                     Int* colperm, std::vector<Int>* dependent_cols) override;
    void _SolveDense(const Vector& rhs, Vector& lhs, char trans) override;
    void _FtranForUpdate(Int nz, const Int* bi, const double* bx) override;
    void _FtranForUpdate(Int nz, const Int* bi, const double* bx,
                         IndexedVector& lhs) override;
    void _BtranForUpdate(Int j) override;
    void _BtranForUpdate(Int j, IndexedVector& lhs) override;
    Int _Update(double pivot) override;
    bool _NeedFreshFactorization() override;
    double _fill_factor() const override;
    double _pivottol() const override;
    void _pivottol(double new_pivottol) override;

    // Reallocates (Li,Lx), (Ui,Ux) and/or (Wi,Wx) as requested by BASICLU.
    void Reallocate();

    // When memory is reallocated, allocate for kReallocFactor*required amount.
    static constexpr double kReallocFactor = 1.5;

    const Control& control_;
    std::vector<Int> istore_;
    std::vector<double> xstore_;
    std::vector<Int> Li_, Ui_, Wi_;
    std::vector<double> Lx_, Ux_, Wx_;
    double fill_factor_;
};

}  // namespace ipx

#endif  // IPX_BASICLU_WRAPPER_H_
