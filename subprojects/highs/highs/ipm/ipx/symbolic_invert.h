#ifndef IPX_SYMBOLIC_INVERT_H_
#define IPX_SYMBOLIC_INVERT_H_

#include <vector>
#include "ipm/ipx/ipx_internal.h"
#include "ipm/ipx/model.h"

namespace ipx {

// Computes the # structural nonzeros per row and column of inverse(B), where
// B = AI[:,basis] is the m-by-m matrix defined by the model and basis.
//
// @basis must have size m and be such that B is structurally nonsingular.
//        (Otherwise an assertion will fail.)
// @rowcounts must either be NULL or an integer array of size m.
//            If not NULL, then on return rowcounts[p], 0 <= p < m, holds the
//            number of nonzeros in row p of inverse(B). Notice that row p of
//            inverse(B) corresponds to column p of B.
// @colcounts must either be NULL or an integer array of size m.
//            If not NULL, then on return colcounts[i], 0 <= i < m, holds the
//            number of nonzeros in column i of inverse(B). Notice that column
//            i of inverse(B) corresponds to row p of B.
//
void SymbolicInvert(const Model& model, const std::vector<Int>& basis,
                    Int* rowcounts, Int* colcounts);

}  // namespace ipx

#endif  // IPX_SYMBOLIC_INVERT_H_
