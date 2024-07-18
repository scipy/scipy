/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file sdiagonal.c
 * \brief Auxiliary routines to work with diagonal elements
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory
 * June 30, 2009
 * </pre>
 */

#include "slu_sdefs.h"

int sfill_diag(int n, NCformat *Astore)
/* fill explicit zeros on the diagonal entries, so that the matrix is not
   structurally singular. */
{
    float *nzval = (float *)Astore->nzval;
    int_t *rowind = Astore->rowind;
    int_t *colptr = Astore->colptr;
    int_t nnz = colptr[n];
    int fill = 0;
    float *nzval_new;
    float zero = 0.0;
    int_t *rowind_new;
    int i, j, diag;

    for (i = 0; i < n; i++)
    {
	diag = -1;
	for (j = colptr[i]; j < colptr[i + 1]; j++)
	    if (rowind[j] == i) diag = j;
	if (diag < 0) fill++;
    }
    if (fill)
    {
	nzval_new = floatMalloc(nnz + fill);
	rowind_new = intMalloc(nnz + fill);
	fill = 0;
	for (i = 0; i < n; i++)
	{
	    diag = -1;
	    for (j = colptr[i] - fill; j < colptr[i + 1]; j++)
	    {
		if ((rowind_new[j + fill] = rowind[j]) == i) diag = j;
		nzval_new[j + fill] = nzval[j];
	    }
	    if (diag < 0)
	    {
		rowind_new[colptr[i + 1] + fill] = i;
		nzval_new[colptr[i + 1] + fill] = zero;
		fill++;
	    }
	    colptr[i + 1] += fill;
	}
	Astore->nzval = nzval_new;
	Astore->rowind = rowind_new;
	SUPERLU_FREE(nzval);
	SUPERLU_FREE(rowind);
    }
    Astore->nnz += fill;
    return fill;
}

int sdominate(int n, NCformat *Astore)
/* make the matrix diagonally dominant */
{
    float *nzval = (float *)Astore->nzval;
    int_t *rowind = Astore->rowind;
    int_t *colptr = Astore->colptr;
    int_t nnz = colptr[n];
    int fill = 0;
    float *nzval_new;
    int_t *rowind_new;
    int_t i, j, diag;
    double s;

    for (i = 0; i < n; i++)
    {
	diag = -1;
	for (j = colptr[i]; j < colptr[i + 1]; j++)
	    if (rowind[j] == i) diag = j;
	if (diag < 0) fill++;
    }
    if (fill)
    {
	nzval_new = floatMalloc(nnz + fill);
	rowind_new = intMalloc(nnz+ fill);
	fill = 0;
	for (i = 0; i < n; i++)
	{
	    s = 1e-6;
	    diag = -1;
	    for (j = colptr[i] - fill; j < colptr[i + 1]; j++)
	    {
		if ((rowind_new[j + fill] = rowind[j]) == i) diag = j;
		s += fabs(nzval_new[j + fill] = nzval[j]);
	    }
	    if (diag >= 0) {
		nzval_new[diag+fill] = s * 3.0;
	    } else {
		rowind_new[colptr[i + 1] + fill] = i;
		nzval_new[colptr[i + 1] + fill] = s * 3.0;
		fill++;
	    }
	    colptr[i + 1] += fill;
	}
	Astore->nzval = nzval_new;
	Astore->rowind = rowind_new;
	SUPERLU_FREE(nzval);
	SUPERLU_FREE(rowind);
    }
    else
    {
	for (i = 0; i < n; i++)
	{
	    s = 1e-6;
	    diag = -1;
	    for (j = colptr[i]; j < colptr[i + 1]; j++)
	    {
		if (rowind[j] == i) diag = j;
		s += fabs(nzval[j]);
	    }
	    nzval[diag] = s * 3.0;
	}
    }
    Astore->nnz += fill;
    return fill;
}
