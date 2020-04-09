// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_BASICLU_KERNEL_H_
#define IPX_BASICLU_KERNEL_H_

#include "lu_factorization.h"

namespace ipx {

class BasicLuKernel : public LuFactorization {
private:
    void _Factorize(Int dim, const Int* Bbegin, const Int* Bend,
                    const Int* Bi, const double* Bx, double pivottol,
                    bool strict_abs_pivottol,
                    SparseMatrix* L, SparseMatrix* U,
                    std::vector<Int>* rowperm, std::vector<Int>* colperm,
                    std::vector<Int>* dependent_cols) override;
};

}  // namespace ipx

#endif  // IPX_BASICLU_KERNEL_H_
