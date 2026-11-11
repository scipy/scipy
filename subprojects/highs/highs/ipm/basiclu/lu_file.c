/* 
 * lu_file.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Data file implementation
 *
 * A data file stores lines of (index,value) pairs. Entries of each line are
 * contiguous in memory. Lines can be in any order in memory and there can be
 * gaps between consecutive lines.
 *
 * The data file implementation uses arrays
 *
 *  index[0:fmem-1], value[0:fmem-1],
 *  begin[0:nlines], end[0:nlines],
 *  next[0:nlines], prev[0:nlines]
 *
 *  index, value    storing (index,value) pairs
 *  begin[k]        pointer to first element in line 0 <= k < nlines,
 *  end[k]          pointer to one past the last element in line k.
 *  begin[nlines]   pointer to the first element of unused space
 *  end[nlines]     holds fmem
 *
 *  next, prev hold the line numbers (0..nlines-1) in a double linked list in
 *  the order in which they appear in memory. That is, for line 0 <= k < nlines,
 *  next[k] and prev[k] are the line which comes before respectively after line
 *  k in memory. next[nlines] and prev[nlines] are the first respectively last
 *  line in memory order.
 *
 * Methods:
 *
 *  lu_file_empty
 *  lu_file_reappend
 *  lu_file_compress
 *  lu_file_diff
 *
 */

#include "ipm/basiclu/lu_def.h"
#include "ipm/basiclu/lu_list.h"
#include "ipm/basiclu/lu_file.h"

/* ==========================================================================
 * lu_file_empty
 *
 * Initialize empty file with @fmem memory space.
 * ========================================================================== */

void lu_file_empty(
    lu_int nlines, lu_int *begin, lu_int *end, lu_int *next, lu_int *prev,
    lu_int fmem)
{
    lu_int i;

    begin[nlines] = 0;
    end[nlines] = fmem;
    for (i = 0; i < nlines; i++)
        begin[i] = end[i] = 0;
    for (i = 0; i < nlines; i++)
    {
        next[i] = i+1;
        prev[i+1] = i;
    }
    next[nlines] = 0;
    prev[0] = nlines;
}


/* ==========================================================================
 * lu_file_reappend
 *
 * Reappend line to file end and add @extra_space elements room. The file must
 * have at least length(line) + @extra_space elements free space.
 * ========================================================================== */

void lu_file_reappend(
    lu_int line, lu_int nlines, lu_int *begin, lu_int *end, lu_int *next,
    lu_int *prev, lu_int *index, double *value, lu_int extra_space)
{
    lu_int fmem, used, room, ibeg, iend, pos;

    fmem = end[nlines];
    used = begin[nlines];
    room = fmem-used;
    ibeg = begin[line];         /* old beginning of line */
    iend = end[line];
    begin[line] = used;         /* new beginning of line */
    assert(iend-ibeg <= room);
    for (pos = ibeg; pos < iend; pos++)
    {
        index[used] = index[pos];
        value[used++] = value[pos];
    }
    end[line] = used;
    room = fmem-used;
    assert(room >= extra_space);
    used += extra_space;
    begin[nlines] = used;       /* beginning of unused space */
    lu_list_move(line, 0, next, prev, nlines, NULL);
}


/* ==========================================================================
 * lu_file_compress
 *
 * Compress file to reuse memory gaps. The ordering of lines in the file is
 * unchanged. To each line with nz entries add @stretch*nz+@pad elements extra
 * space. Chop extra space if it would overlap the following line in memory.
 *
 * Return: number of entries in file
 * ========================================================================== */

lu_int lu_file_compress(
    lu_int nlines, lu_int *begin, lu_int *end, const lu_int *next,
    lu_int *index, double *value, double stretch, lu_int pad)
{
    lu_int i, pos, ibeg, iend, used, extra_space, nz = 0;

    used = 0; extra_space = 0;
    for (i = next[nlines]; i < nlines; i = next[i]) /* move line i */
    {
        ibeg = begin[i];
        iend = end[i];
        assert(ibeg >= used);
        used += extra_space;
        if (used > ibeg)
            used = ibeg;        /* chop extra space added before */
        begin[i] = used;
        for (pos = ibeg; pos < iend; pos++)
        {
            index[used] = index[pos];
            value[used++] = value[pos];
        }
        end[i] = used;
        extra_space = stretch*(iend-ibeg) + pad;
        nz += iend-ibeg;
    }
    assert(used <= begin[nlines]);
    used += extra_space;
    if (used > begin[nlines])
        used = begin[nlines];   /* never use more space than before */
    begin[nlines] = used;
    return nz;
}


/* ==========================================================================
 * lu_file_diff (for debugging)
 *
 * @begin_row, @end_row, @begin_col, @end_col are pointer into @index, @value,
 * defining lines of the "row file" and the "column file".
 *
 * Task:
 *
 *  val == NULL: count row file entries that are missing in column file.
 *  val != NULL: count row file entries that are missing in column file
 *               or which values are different.
 *
 * The method does a column file search for each row file entry. To check
 * consistency of rowwise and columnwise storage, the method must be called
 * twice with row pointers and column pointers swapped.
 * ========================================================================== */

lu_int lu_file_diff(
    lu_int nrow, const lu_int *begin_row, const lu_int *end_row,
    const lu_int *begin_col, const lu_int *end_col, const lu_int *index,
    const double *value)
{
    lu_int i, j, pos, where, ndiff = 0;

    for (i = 0; i < nrow; i++)
    {
        for (pos = begin_row[i]; pos < end_row[i]; pos++)
        {
            j = index[pos];
            for (where = begin_col[j]; where < end_col[j] &&
                     index[where] != i; where++)
                ;
            if (where == end_col[j] || (value && value[pos] != value[where]))
                ndiff++;
        }
    }
    return ndiff;
}
