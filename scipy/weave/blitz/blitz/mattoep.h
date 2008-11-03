/***************************************************************************
 * blitz/mattoep.h      Declarations for Toeplitz matrices
 *
 * $Id: mattoep.h 1414 2005-11-01 22:04:59Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

#ifndef BZ_MATTOEP_H
#define BZ_MATTOEP_H

#ifndef BZ_MSTRUCT_H
 #error <blitz/mattoep.h> must be included via <blitz/mstruct.h>
#endif

BZ_NAMESPACE(blitz)

// Toeplitz matrix
// [ 0 1 2 3 ]
// [ 1 2 3 4 ]
// [ 2 3 4 5 ]
// [ 3 4 5 6 ]

class ToeplitzIterator {
public:
    ToeplitzIterator(unsigned rows, unsigned cols)
    {
        rows_ = rows;
        cols_ = cols;
        i_ = 0;
        j_ = 0;
        good_ = true;
        offset_ = 0;
    }

    operator bool() const { return good_; }

    void operator++()
    {
        ++offset_;
        if (i_ < rows_ - 1)
            ++i_;
        else if (j_ < cols_ - 1)
            ++j_;
        else
            good_ = false;
    }

    unsigned row() const
    { return i_; }

    unsigned col() const
    { return j_; }

    unsigned offset() const
    { return offset_; }

protected:
    unsigned offset_;
    unsigned i_, j_;
    unsigned rows_, cols_;
    bool     good_;
};

class Toeplitz : public GeneralMatrix {

public:
    typedef ToeplitzIterator T_iterator;

    Toeplitz()
        : rows_(0), cols_(0)
    { }

    Toeplitz(unsigned rows, unsigned cols)
        : rows_(rows), cols_(cols)
    { }

    unsigned columns() const
    { return cols_; }

    unsigned coordToOffset(unsigned i, unsigned j) const
    {
        BZPRECONDITION(inRange(i,j));
        return i + j;
    }

    unsigned firstInRow(unsigned i) const
    { return 0; }

    template<typename T_numtype>
    T_numtype get(const T_numtype * restrict data,
        unsigned i, unsigned j) const
    {
        BZPRECONDITION(inRange(i,j));
        return data[coordToOffset(i,j)];
    }

    template<typename T_numtype>
    T_numtype& get(T_numtype * restrict data, unsigned i, unsigned j)
    {
        BZPRECONDITION(inRange(i,j));
        return data[coordToOffset(i,j)];
    }

    unsigned lastInRow(const unsigned)  const { return cols_ - 1; }
    unsigned firstInCol(const unsigned) const { return 0; }
    unsigned lastInCol(const unsigned)  const { return rows_ - 1; }

    bool inRange(const unsigned i,const unsigned j) const { return (i<rows_) && (j<cols_); }

    unsigned numElements() const { return rows_ + cols_ - 1; }

    unsigned rows() const { return rows_; }

    void resize(const unsigned rows,const unsigned cols) {
        rows_ = rows;
        cols_ = cols;
    }

private:
    unsigned rows_, cols_;
};

BZ_NAMESPACE_END

#endif // BZ_MATSYMM_H

