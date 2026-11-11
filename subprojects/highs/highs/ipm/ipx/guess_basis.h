#ifndef IPX_GUESS_BASIS_H_
#define IPX_GUESS_BASIS_H_

#include "ipm/ipx/control.h"
#include "ipm/ipx/model.h"

namespace ipx {

// Returns a vector of m column indices (m = model.rows()) as a guess for a
// basis. Columns with a larger weight are chosen preferably. The basis can be
// singular (and usually is).
//
// @colweights size n+m array (n = model.cols()) with entries in the closed
//             interval [0,INFINITY].
//
std::vector<Int> GuessBasis(const Control& control, const Model& model,
                            const double* colweights, bool& interrupt);

}  // namespace ipx

#endif // IPX_GUESS_BASIS_H_
