// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_MAXVOLUME_H_
#define IPX_MAXVOLUME_H_

#include "basis.h"
#include "control.h"

namespace ipx {

class Maxvolume {
public:
    explicit Maxvolume(const Control& control);
    Int RunSequential(const double* colscale, Basis& basis);
    Int RunHeuristic(const double* colscale, Basis& basis);

    Int updates() const;        // # basis updates
    Int skipped() const;        // # columns computed but basis not updated
    Int passes() const;         // # passes over tblmatrix (-1 for heuristic)
    Int slices() const;         // # slices of tblmatrix formed (heuristic only)
    double volinc() const;      // log2 of volume increase
    double time() const;        // runtime

    // The remaining statistics are only set by RunSequential() and count the
    // final pass over the tableau matrix. They are meaningful only if no
    // updates occured in the final pass; this is guaranteed if parameter
    // max_passes was larger than passes().

    Int tblnnz() const;              // nnz of tableau matrix
    double tblmax() const;           // max entry of scaled tableau matrix
    double frobnorm_squared() const; // Frob.norm squared of scaled tableau

private:
    static constexpr double kPivotZeroTol = 1e-7;

    struct Slice;
    void Reset();
    Int Driver(Basis& basis, Slice& slice);
    Int ScaleFtran(double colscale_jn, const Vector& invscale_basic,
                   IndexedVector& ftran);

    const Control& control_;
    Int updates_{0};
    Int skipped_{0};
    Int passes_{0};
    Int slices_{0};
    double volinc_{0.0};
    double time_{0.0};
    Int tblnnz_{0};
    double tblmax_{0.0};
    double frobnorm_squared_{0.0};
};

}  // namespace ipx

#endif  // IPX_MAXVOLUME_H_
