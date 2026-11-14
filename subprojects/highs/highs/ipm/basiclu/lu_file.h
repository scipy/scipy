#ifndef _LU_FILE_H
#define _LU_FILE_H

void lu_file_empty(
    lu_int nlines, lu_int *begin, lu_int *end, lu_int *next, lu_int *prev,
    lu_int fmem);

void lu_file_reappend(
    lu_int line, lu_int nlines, lu_int *begin, lu_int *end, lu_int *next,
    lu_int *prev, lu_int *index, double *value, lu_int extra_space);

lu_int lu_file_compress(
    lu_int nlines, lu_int *begin, lu_int *end, const lu_int *next,
    lu_int *index, double *value, double stretch, lu_int pad);

lu_int lu_file_diff(
    lu_int nrow, const lu_int *begin_row, const lu_int *end_row,
    const lu_int *begin_col, const lu_int *end_col, const lu_int *index,
    const double *value);

#endif
